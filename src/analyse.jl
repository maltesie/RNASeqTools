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