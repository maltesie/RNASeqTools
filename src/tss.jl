using XAM
include("utils.jl")
include("io.jl")

function coverage(bam_file::String; norm=1000000, unique_mappings_only=true)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    ref_lengths = bam_chromosome_lengths(reader)
    ref_names = bam_chromosome_names(reader)
    coverage_f = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    coverage_r = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    count = 0
    while !eof(reader)
        count += 1
        read!(reader, record)
        !BAM.ismapped(record) && continue
        ref = BAM.refname(record)
        pos = BAM.position(record)
        len = BAM.alignlength(record)
        if unique_mappings_only
            aux = get_XA_tag(BAM.auxdata(record).data)
            (aux == "-") || continue
        end
        (strandint(record) == 1) ? (coverage_f[ref][pos:pos+len-1] .+= 1) : (coverage_r[ref][pos:pos+len-1] .-= 1)
    end
    close(reader)
    norm_factor = norm/count
    return Dict(ref => cov .* norm_factor for (ref,cov) in coverage_f), Dict(ref => cov .* norm_factor for (ref,cov) in coverage_r)
end

function diff(coverage::Vector{Float64})
    d = similar(coverage)
    d[1] = coverage[1]
    d[2] = coverage[2] - coverage[1]
    d[3] = maximum(coverage[3] .- coverage[1:2])
    for (i,val) in enumerate(coverage[4:end])
        d[i+3] = maximum(val .- coverage[i:i+2])
    end
    return d
end

function tss(notex_fs::Vector{String}, notex_rs::Vector{String}, tex_fs::Vector{String}, tex_rs::Vector{String}; 
    min_step=10, min_ratio=1.5)
    result::Dict{String, DataFrame} = Dict()
    for i in 1:length(notex_fs)
        forward = read_wig(notex_fs[i])[1]
        reverse = read_wig(notex_rs[i])[1]
        forward_tex = read_wig(tex_fs[i])[1]
        reverse_tex = read_wig(tex_rs[i])[1]
        for chr in keys(forward)
            (notex_f, tex_f) = make_same_length(forward[chr], forward_tex[chr])
            (notex_r, tex_r) = make_same_length(reverse[chr], reverse_tex[chr])
            result[chr] = DataFrame(pos=Int[], val=Float64[])
            d_forward = diff(tex_f)
            d_reverse = diff(tex_r)
            check_forward = ((tex_f ./ notex_f) .>= min_ratio) .& (d_forward .>= min_step)
            check_reverse = ((tex_r ./ notex_r) .>= min_ratio) .& (d_reverse .<= -min_step)
            append!(result[chr], DataFrame(pos=findall(!iszero, check_forward), val=d_forward[check_forward]))
            append!(result[chr], DataFrame(pos=findall(!iszero, check_reverse) .* -1, val=d_reverse[check_reverse]))
        end
    end
    for (chr, ts) in result
        sort!(ts, :pos)
    end
    return result
end

function terms(coverage_fs::Vector{String}, coverage_rs::Vector{String}; min_step=10)
    result::Dict{String, DataFrame} = Dict()
    for i in 1:length(coverage_fs)
        forward = read_wig(coverage_fs[i])[1]
        reverse = read_wig(coverage_rs[i])[1]
        for chr in keys(forward)
            result[chr] = DataFrame(pos=Int[], val=Float64[])
            d_forward = diff(forward[chr])
            d_reverse = diff(reverse[chr])
            check_forward = d_forward .<= -min_step
            check_reverse = d_reverse .>= min_step
            append!(result[chr], DataFrame(pos=findall(!iszero, check_forward), val=d_forward[check_forward]))
            append!(result[chr], DataFrame(pos=findall(!iszero, check_reverse) .* -1, val=d_reverse[check_reverse]))
        end
    end
    for (chr, ts) in result
        sort!(ts, :pos)
    end
    return result
end

function annotate_utrs!(annotations::Dict{String, DataFrame}, tss::Dict{String, DataFrame}, terms::Dict{String, DataFrame}; 
    max_distance=300, guess_distance=150)
    for (chr, annotation) in annotations
        annotation[!, :fiveUTR] = annotation[!, :start] -= guess_distance
        annotation[!, :fiveType] = fill("guess", nrows(annotation))
        annotation[!, :threeUTR] = annotation[!, :stop] += guess_distance
        annotation[!, :threeType] = fill("guess", nrows(annotation))
        for row in eachrow(annotation)
            five_hits = @view tss[chr][row[:start]-max_distance .<= tss[chr][!,:pos] .<= row[:start]]
            three_hits = @view terms[chr][row[:stop] .<= terms[chr] .<= row[:stop]+max_distance]
            isempty(five_hits) || (row[:fiveUTR]=maximum(five_hits[!, :pos]); row[:fiveType]="max")
            isempty(three_hits) || (row[:threeUTR]=maximum(three_hits[!, :pos]); row[:threeType]="max")
        end
    end
end