function similarity(read1::LongDNASeq, read2::LongDNASeq; window_size=10, step_size=2)
    (length(read1) > length(read2)) ? ((short_seq, long_seq) = (read2, read1)) : ((short_seq, long_seq) = (read1, read2))
    nb_windows = round(Int, (length(short_seq) - window_size) / 2, RoundUp)
    @assert nb_windows > 0
    window_counter = 0
    for window_start in 1:step_size:length(short_seq)-window_size
        approxoccursin(long_seq, short_seq[window_start:window_start+window_size]) && (window_counter += 1)
    end
    return window_counter / nb_windows
end

function similarity(reads::PairedReads; window_size=10, step_size=2)
    similarities = Dict{UInt, Float64}()
    for (key, (read1, read2)) in reads.dict
        push!(similarities, key=>similarity(read1, read2; window_size=window_size, step_size=step_size))
    end
    return similarities
end