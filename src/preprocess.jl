function split_libs(infile1::String, prefixfile1::String, infile2::String, barcodes::Array{String,1}, libnames::Array{String, 1}, 
    output_folder::String; stop_early=-1, report_file="demultiplexing_report.txt")

    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)
    output_files = [[abspath(output_folder, "$(name)_1.fastq.gz"),
                    abspath(output_folder, "$(name)_2.fastq.gz")] for name in libnames]
    stats::Array{Int,1} = zeros(Int, length(libnames))
    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    recordp::FASTQ.Record = FASTQ.Record()
    endswith(infile1, ".gz") ? reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r"))) : reader1 = FASTQ.Reader(open(infile1, "r"))
    endswith(infile2, ".gz") ? reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r"))) : reader2 = FASTQ.Reader(open(infile2, "r"))
    endswith(prefixfile1, ".gz") ? readerp = FASTQ.Reader(GzipDecompressorStream(open(prefixfile1, "r"))) : readerp = FASTQ.Reader(open(prefixfile1, "r"))
    writers = [[FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2)),
                FASTQ.Writer(GzipCompressorStream(open(outfile2, "w"), level=2))] 
                for (outfile1, outfile2) in output_files] 

    sleep(0.01)

    c = 0
    while !eof(reader1)
        read!(reader1, record1)
        read!(reader2, record2)
        read!(readerp, recordp)
        ((c >= stop_early) & (stop_early > 0)) && break
        c += 1
        prefix = LongDNASeq(FASTQ.sequence(recordp))
        (library_id, nb_errors) = demultiplex(dplxr, prefix)
        if (nb_errors == -1)
            read = LongDNASeq(FASTQ.sequence(record1))
            (library_id, nb_errors) = demultiplex(dplxr, read)
            (nb_errors == -1) && continue
        end
        stats[library_id] += 1
        write(writers[library_id][1], record1)
        write(writers[library_id][2], record2)
    end

    for (writer1, writer2) in writers
        close(writer1)
        close(writer2)
    end

    close(reader1)
    close(reader2)

    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libnames, stats)])
    count_string *= "not identifyable - $(c-sum(stats))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write_file(abspath(output_folder, report_file), count_string)

end

function split_libs(infile1::String, infile2::String, barcodes::Array{String,1}, libnames::Array{String, 1}, 
    output_folder::String; stop_early=-1, report_file="demultiplexing_report.txt")

    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)
    output_files = [[abspath(output_folder, "$(name)_1.fastq.gz"),
                    abspath(output_folder, "$(name)_2.fastq.gz")] for name in libnames]
    stats::Array{Int,1} = zeros(Int, length(libnames))
    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    endswith(infile1, ".gz") ? reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r"))) : reader1 = FASTQ.Reader(open(infile1, "r"))
    endswith(infile2, ".gz") ? reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r"))) : reader2 = FASTQ.Reader(open(infile2, "r"))
    writers = [[FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2)),
                FASTQ.Writer(GzipCompressorStream(open(outfile2, "w"), level=2))] 
                for (outfile1, outfile2) in output_files] 

    sleep(0.01)

    c = 0
    while !eof(reader1)
        read!(reader1, record1)
        read!(reader2, record2)
        ((c >= stop_early) & (stop_early > 0)) && break
        c += 1
        read = LongDNASeq(FASTQ.sequence(record1))
        (library_id, nb_errors) = demultiplex(dplxr, read)
        nb_errors == -1 && continue
        stats[library_id] += 1
        write(writers[library_id][1], record1)
        write(writers[library_id][2], record2)
    end

    for (writer1, writer2) in writers
        close(writer1)
        close(writer2)
    end

    close(reader1)
    close(reader2)

    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libnames, stats)])
    count_string *= "not identifyable - $(c-sum(stats))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write_file(abspath(output_folder, report_file), count_string)

end

function split_libs(infile::String, barcodes::Array{String,1}, libnames::Array{String, 1}, 
    output_folder::String; stop_early=-1, report_file="demultiplexing_report.txt")

    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)
    output_files = [abspath(output_folder, "$(name).fastq.gz") for name in libnames]
    stats::Array{Int,1} = zeros(Int, length(libnames))
    record::FASTQ.Record = FASTQ.Record()
    endswith(infile, ".gz") ? reader = FASTQ.Reader(GzipDecompressorStream(open(infile, "r"))) : reader = FASTQ.Reader(open(infile, "r"))
    writers = [FASTQ.Writer(GzipCompressorStream(open(outfile, "w"), level=2)) for outfile in output_files] 

    sleep(0.01)

    c = 0
    while !eof(reader)
        read!(reader, record)
        ((c >= stop_early) & (stop_early > 0)) && break
        c += 1
        read = LongDNASeq(FASTQ.sequence(record))
        (library_id, nb_errors) = demultiplex(dplxr, read)
        nb_errors == -1 && continue
        stats[library_id] += 1
        write(writers[library_id], record)
    end

    for writer in writers
        close(writer)
    end
    close(reader)

    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libnames, stats)])
    count_string *= "not identifyable - $(c-sum(stats))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write_file(abspath(output_folder, report_file), count_string)
end

function trim_fastp(input_files::Vector{String}, output_folder::String; 
    suffix="_trimmed", adapters=String[], fastp_bin="fastp",
    min_length=25, write_report=true)

    record = FASTQ.Record()
    report = ""
    output_files = [split(in_file, ".fastq")[1]*"$(suffix).fastq.gz" for in_file in input_files]
    isempty(adapters) && (adapters = fill("", length(input_files)))
    c = 0
    for (in_file, out_file, adapter) in zip(input_files, output_files, adapters)
        
        isnothing(adapter) ? 
        cmd = `$fastp_bin -i $in_file -o $out_file -M 25 -l $min_length --trim_poly_x 10 --cut_front --cut_tail` :
        cmd = `$fastp_bin -a $adapter -i $in_file -o $out_file -M 25 -l $min_length --trim_poly_x 10 --cut_front --cut_tail`
        run(cmd)
        cc = 0
        reader = FASTQ.Reader(GzipDecompressorStream(open(out_file, "r")))
        while !eof(reader)
            read!(reader, record)
            cc += 1
        end
        report *= "$in_file\t$(cc)\n"
        c += cc
    end 
    report *= "$(c) reads in total\n"
    write_report && write_file(abspath(output_folder, "trimming_report.txt"), report)
end

function trim_fastp(input_files::Vector{Vector{String}}, output_folder::String; 
    suffix="_trimmed", adapters=String[], fastp_bin="fastp",
    min_length=25, write_report=true)

    record = FASTQ.Record()
    report = ""
    output_files = [[split(in_file1, ".fastq")[1]*"$(suffix).fastq.gz",
            split(in_file2, ".fastq")[1]*"$(suffix).fastq.gz"] for (in_file1, in_file2) in input_files]
    isempty(adapters) && (adapters = fill("", length(input_files)))
    c = 0
    for ((in_file1, in_file2), (out_file1, out_file2), adapter) in zip(input_files, output_files, adapters)
        
        isnothing(adapter) ? 
        cmd = `$fastp_bin -i $in_file1 -o $out_file1 -I $in_file2 -O $out_file2 -M 25 -l $min_length --trim_poly_x 10 --cut_front --cut_tail` :
        cmd = `$fastp_bin -a $adapter $in_file1 -o $out_file1 -I $in_file2 -O $out_file2 -M 25 -l $min_length --trim_poly_x 10 --cut_front --cut_tail`
        run(cmd)
        cc = 0
        reader = FASTQ.Reader(GzipDecompressorStream(open(out_file, "r")))
        while !eof(reader)
            read!(reader, record)
            cc += 1
        end
        report *= "$in_file\t$(cc)\n"
        c += cc
    end 
    report *= "$(c) reads in total\n"
    write_report && write_file(abspath(output_folder, "trimming_report.txt"), report)
end