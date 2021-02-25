function Base.filter!(f, reads::Reads)
    for (key, read) in reads.dict
        f(read) || delete!(reads, key)
    end
end

function Base.filter!(f, reads::PairedReads; both=false)
    for (key, (read1, read2)) in reads.dict
        both ? (f(read1) && f(read2)) || delete!(reads, key) : (f(read1) || f(read2)) || delete!(reads, key)
    end
end