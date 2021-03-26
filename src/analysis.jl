function similarity(read1::LongDNASeq, read2::LongDNASeq; score_model=nothing)
    isnothing(score_model) && (score_model = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1))
    (length(read1) > length(read2)) ? ((short_seq, long_seq) = (read2, read1)) : ((short_seq, long_seq) = (read1, read2))
    aln = local_alignment(long_seq, short_seq, score_model)
    hasalignment(aln) || (return 0.0) 
    return count_matches(alignment(aln))/length(short_seq)
end

function similarity(reads::PairedReads; window_size=10, step_size=2)
    similarities = Dict{UInt, Float64}()
    score_model = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1)
    for (read1, read2) in reads
        push!(similarities, key=>similarity(read1, read2; score_model=score_model))
    end
    return similarities
end

function nucleotide_count(reads::Reads; normalize=true)
    max_length = maximum([length(read) for read in reads])
    count = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    nb_reads = length(reads)
    for read in reads
        (align==:left) ? 
        (index = 1:length(read)) : 
        (index = (max_length - length(read) + 1):max_length)
        for (i, n) in zip(index, read)
            count[n][i] += 1
        end
    end
    if normalize
        for (key, c) in count
            c /= length(reads)
        end
    end
    return count
end

function nucleotide_count(reads::PairedReads; normalize=true)
    max_length = maximum(vcat([[length(read1) length(read2)] for (read1, read2) in reads]...))
    count1 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    count2 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    nb_reads = length(reads)
    for (read1, read2) in reads
        (align==:left) ? 
        (index1 = 1:length(read1); index2 = 1:length(read2)) : 
        (index1 = (max_length - length(read1) + 1):max_length; index2 = (max_length - length(read2) + 1):max_length)
        for ((i1, n1),(i2, n2)) in zip(zip(index1, read1), zip(index2, read2))
            count1[n1][i1] += 1
            count2[n2][i2] += 1
        end
    end
    if normalize
        for ((key1, c1), (key2, c2)) in zip(count1, count2)
            c1 /= length(reads)
            c2 /= length(reads)
        end
    end
    return count1, count2
end

function diff(coverage::Vector{Float64})
    d = zeros(Float64,length(coverage))
    d[1] = coverage[1]
    d[2:end] = coverage[2:end] - coverage[1:end-1]
    return d
end

function tss(notex_f::Coverage, notex_r::Coverage, tex_f::Coverage, tex_r::Coverage; min_step=10, min_ratio=1.3)
    @assert notex_f.chroms == notex_r.chroms == tex_f.chroms == tex_r.chroms
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