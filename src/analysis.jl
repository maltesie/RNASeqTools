using XAM

function isproperpair(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_PROPER_PAIR) != 0
end

function areconcordant(pos1::Int, pos2::Int, chr1::String, chr2::String, aux1::String, aux2::String; 
    read1_reversed=false, read2_reversed=false)::Bool

    read1_reversed ? read1_factor = -1 : read1_factor = 1
    poss1::Array{Int,1} = [pos1]
    chrs1::Array{String,1} = [chr1]
    for alignment in split(aux1, ";")
        (isempty(alignment) | (alignment == "-")) && continue
        push!(poss1, read1_factor * parse(Int, split(alignment, ",")[2]))
        push!(chrs1, String(split(alignment, ",")[1]))
    end

    read2_reversed ? read2_factor = -1 : read2_factor = 1
    poss2::Array{Int,1} = [pos2]
    chrs2::Array{String,1} = [chr2]
    for alignment in split(aux2, ";")
        (isempty(alignment) | (alignment == "-")) && continue
        push!(poss2, read2_factor * parse(Int, split(alignment, ",")[2]))
        push!(chrs2, split(alignment, ",")[1])
    end

    for (p1::Int, c1::String) in zip(poss1, chrs1)
        for (p2::Int, c2::String) in zip(poss2, chrs2)
            c1 != c2 && continue
            (abs(p1-p2) < 1000) && (return true)
        end
    end
    return false
end

function single_fragment_set(bam_file::String; nb_reads::Int=-1)
    reader = open(BAM.Reader, bam_file)
    record::BAM.Record = BAM.Record()
    single_fragments::Array{String, 1} = []
    c::Int = 0 
    while !eof(reader)
        read!(reader, record)
        (isproperpair(record) & BAM.ismapped(record) & BAM.isnextmapped(record)) && push!(single_fragments, BAM.tempname(record))
        c += 1
        ((nb_reads > 0) & (c >= nb_reads)) && break 
    end
    close(reader)
    single_fragment_set::Set{String} = Set(single_fragments) 
    return single_fragment_set
end