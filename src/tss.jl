function coverage(bam_file::String; norm=10000000, is_reversed=false, unique_mappings_only=true)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    ref_lengths = bam_chromosome_lengths(reader)
    ref_names = bam_chromosome_names(reader)
    coverage_f = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    coverage_r = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    count = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        ref = BAM.refname(record)
        pos = BAM.position(record)
        len = BAM.alignlength(record)
        (unique_mappings_only && (get_XA_tag(BAM.auxdata(record).data) != "-")) && continue
        count += 1
        (BAM.ispositivestrand(record) == !is_reversed)  ? (coverage_f[ref][pos:pos+len-1] .+= 1) : (coverage_r[ref][pos:pos+len-1] .-= 1)
    end
    close(reader)
    #println(bam_file, ": ", count)
    norm_factor = norm/count
    return Dict(ref => cov .* norm_factor for (ref,cov) in coverage_f), Dict(ref => cov .* norm_factor for (ref,cov) in coverage_r)
end

function diff(coverage::Vector{Float64})
    d = zeros(Float64,length(coverage))
    d[1] = coverage[1]
    d[2:end] = coverage[2:end] - coverage[1:end-1]
    return d
end

function tss(notex_fs::Vector{String}, notex_rs::Vector{String}, tex_fs::Vector{String}, tex_rs::Vector{String}; 
    min_step=10, min_ratio=1.3)
    chrs = get_chr_from_wig(notex_fs[1])
    result::Dict{String, DataFrame} = Dict(chr=>DataFrame(pos=Int[], val=Float64[]) for chr in chrs)
    for i in 1:length(notex_fs) # für 0.1, dann 2.0
        forward = read_wig(notex_fs[i]) # gibt vector mit replikaten zurück
        reverse = read_wig(notex_rs[i])
        forward_tex = read_wig(tex_fs[i])
        reverse_tex = read_wig(tex_rs[i])
        #write(open("/home/malte/Workspace/data/vibrio/reversetex_$(i).txt", "w"), join(["$j $d" for (j,d) in enumerate(reverse_tex[1]["NC_002505"])], '\n'))
        for chr in chrs
            (tex_f, notex_f) = join_replicates(forward, forward_tex, chr)
            (tex_r, notex_r) = join_replicates(reverse, reverse_tex, chr)
            d_forward = diff(tex_f)
            d_reverse = diff(tex_r)
            #write(open("/home/malte/Workspace/data/vibrio/dreverse_$(i)_$(chr).txt", "w"), join(["$j $d" for (j,d) in enumerate(d_reverse)], '\n'))
            #write(open("/home/malte/Workspace/data/vibrio/dreverse_$(i)_$(chr).txt", "w"), join(["$j $d" for (j,d) in enumerate(d_reverse)], '\n'))
            check_forward = circshift(((tex_f ./ notex_f) .>= min_ratio), 1) .& (d_forward .>= min_step)
            check_reverse = circshift(((tex_r ./ notex_r) .>= min_ratio), 1) .& (d_reverse .>= min_step)
            append!(result[chr], DataFrame(pos=findall(!iszero, check_forward), val=abs.(d_forward[check_forward])))
            append!(result[chr], DataFrame(pos=findall(!iszero, check_reverse) .* -1, val=abs.(d_reverse[check_reverse])))
            #CSV.write("/home/malte/Workspace/data/vibrio/test_$(i)_$(chr).csv", result[chr])
        end
    end
    for (chr, ts) in result
        sort!(ts, :pos)
    end
    #CSV.write("/home/malte/Workspace/data/vibrio/test.csv", result["NC_002505"])
    return result
end

function terms(coverage_fs::Vector{String}, coverage_rs::Vector{String}; min_step=10)
    chrs = get_chr_from_wig(coverage_fs[1])
    result::Dict{String, DataFrame} = Dict(chr=>DataFrame(pos=Int[], val=Float64[]) for chr in chrs)
    for i in 1:length(coverage_fs)
        forward = read_wig(coverage_fs[i])
        reverse = read_wig(coverage_rs[i])
        for chr in chrs
            f = join_replicates(forward, chr)
            r = join_replicates(reverse, chr)
            d_forward = diff(f)
            d_reverse = diff(r)
            check_forward = d_forward .<= -min_step
            check_reverse = d_reverse .<= -min_step
            append!(result[chr], DataFrame(pos=findall(!iszero, check_forward), val=abs.(d_forward[check_forward])))
            append!(result[chr], DataFrame(pos=findall(!iszero, check_reverse) .* -1, val=abs.(d_reverse[check_reverse])))
        end
    end
    for (chr, ts) in result
        sort!(ts, :pos)
    end
    #CSV.write("/home/malte/Workspace/data/vibrio/test.csv", result["NC_002505"])
    return result
end

function annotate_utrs!(annotations::Dict{String, DataFrame}, tss::Dict{String, DataFrame}, terms::Dict{String, DataFrame}; 
    max_distance=150, add_distance=150, guess_distance=150)
    for (chr, annotation) in annotations
        annotation[!, :fiveUTR] = annotation[!, :start] .- guess_distance
        annotation[!, :fiveType] = fill("guess", nrow(annotation))
        annotation[!, :threeUTR] = annotation[!, :stop] .+ guess_distance
        annotation[!, :threeType] = fill("guess", nrow(annotation))
        for row in eachrow(annotation)
            five_hits = tss[chr][row[:start]-max_distance .<= tss[chr][!,:pos] .<= row[:start], :]
            three_hits = terms[chr][row[:stop] .<= terms[chr][!,:pos] .<= row[:stop]+max_distance, :]
            isempty(five_hits) && (five_hits = tss[chr][row[:start]-max_distance-add_distance .<= tss[chr][!,:pos] .<= row[:start]-max_distance, :])
            isempty(three_hits) && (three_hits = terms[chr][row[:stop]+max_distance .<= terms[chr][!,:pos] .<= row[:stop]+max_distance+add_distance, :])
            isempty(five_hits) || (row[:fiveUTR]=five_hits[argmax(five_hits[!, :val]), :pos]; row[:fiveType]="max")
            isempty(three_hits) || (row[:threeUTR]=three_hits[argmax(three_hits[!, :val]), :pos]; row[:threeType]="max")
            #(row[:stop]==-372) && (println(three_hits))
        end
    end
end