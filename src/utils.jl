using FASTX
using CodecZlib
using BioSequences

function combine_gffs(gff_files::Array{String, 1}; out_file="out.gff3")
    writer = open(out_file, "w")
    for gff_in_file in gff_files
        reader = open(gff_in_file)
        write(writer, read(reader, String))
        close(reader)
    end
    close(writer)
end

function reverse_reads(in_file::String, out_file::String; complement=true)

    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    ide::String = ""
    des::String = ""
    seq::String = ""
    qua::String = ""

    if endswith(in_file, ".fastq.gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(in_file, "r")))
    elseif endswith(in_file, ".fastq")
        reader = FASTQ.Reader(open(in_file, "r"))
    else
        throw(ErrorException("Only .fastq and .fastq.gz files are supported."))
    end
    writer = FASTQ.Writer(GzipCompressorStream(open(out_file, "w")))
    
    while !eof(reader)
        read!(reader, record1)
        complement ? seq = convert(String, reverse_complement(LongDNASeq(sequence(record1)))) : seq = reverse(sequence(record1))
        record2 = FASTQ.Record(identifier(record1), description(record1), seq, reverse(quality(record1)))
        write(writer, record2)
    end
    close(reader)
    close(writer)
end

function complement_reads(in_file::String, out_file::String)

    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    ide::String = ""
    des::String = ""
    seq::String = ""
    qua::String = ""

    if endswith(in_file, ".fastq.gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(in_file, "r")))
    elseif endswith(in_file, ".fastq")
        reader = FASTQ.Reader(open(in_file, "r"))
    else
        throw(ErrorException("Only .fastq and .fastq.gz files are supported."))
    end
    writer = FASTQ.Writer(GzipCompressorStream(open(out_file, "w")))
    
    while !eof(reader)
        read!(reader, record1)
        seq = convert(String, complement(LongDNASeq(sequence(record1))))
        record2 = FASTQ.Record(identifier(record1), description(record1), seq, reverse(quality(record1)))
        write(writer, record2)
    end
    close(reader)
    close(writer)
end

function get_position(pos::Int, chr::String, aux::String; reversed=false)::Tuple{Int, String}
    poss::Array{Int,1} = [pos]
    chrs::Array{String,1} = [chr]
    reversed ? factor = -1 : factor = 1
    for alignment in split(aux, ";")
        (alignment == "-") && (return pos, chr)
        isempty(alignment) && continue
        c, p, q, e = split(alignment, ",")
        push!(poss, factor * parse(Int, p))
        push!(chrs, String(c))
    end
    ind = sortperm(abs.(poss))[1]
    spos = poss[ind]
    schr = chrs[ind]
    return spos, schr
end