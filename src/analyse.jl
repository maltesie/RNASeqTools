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
    for (key, (read1, read2)) in reads.dict
        push!(similarities, key=>similarity(read1, read2; score_model=score_model))
    end
    return similarities
end