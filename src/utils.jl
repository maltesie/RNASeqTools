using FASTX

function combine_gffs(gff_files::Array{String, 1}; out_file="out.gff3")
    writer = open(out_file, "w")
    for gff_in_file in gff_files
        reader = open(gff_in_file)
        write(writer, read(reader, String))
        close(reader)
    end
    close(writer)
end

function flip_reads(in_file::String, out_file::String)

    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()

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
        record2 = FASTQ.Record(identifier(record1), reverse(sequence(record1)), reverse(quality(record1)))
        write(writer, record2)
    end
    close(reader)
    close(writer)
end