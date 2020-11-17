import BioSequences: Demultiplexer, LongDNASeq, demultiplex, sequence
import FASTX: FASTQ.Reader, FASTQ.Record, FASTQ.Writer, sequence
import CodecZlib: GzipCompressorStream, GzipDecompressorStream

function write_file(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

function split_libs(infile1::String, infile2::String, barcodes::Array{String,1}, libnames::Array{String, 1}, 
    output_folder::String; bc_in_file1=true, barcode_length=-1, stop_early=-1, write_report=false)

    barcode_length < 0 && (barcode_length = length(barcodes[1]))
    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)
    output_files = [[abspath(output_folder, "$(name)_1.fastq.gz"),
                    abspath(output_folder, "$(name)_2.fastq.gz")] for (in_file1, in_file2) in input_files]
    stats::Array{Int,1} = zeros(Int, length(input_files))
    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r")))
    reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r")))
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
        bc_in_file1 ? read = LongDNASeq(sequence(record1)) : read = LongDNASeq(sequence(record2))
        (library_id, nb_errors) = demultiplex(dplxr, read)
        nb_errors == -1 && continue
        stats[library_id] += 1
        write(writers[library_id][1], record_1)
        write(writers[library_id][2], record_2)
    end

    for (writer_1, writer_2) in writers
        close(writer_1)
        close(writer_2)
    end

    close(reader_1)
    close(reader_2)

    count_string = join(["$(names[i]) - $(stats_dict[i])\n" for i in 1:length(libnames)])
    count_string *= "not identifyable - $(c-sum(stats_dict))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write_report && write_file(abspath(output_folder, "demultiplexing_report.txt"), count_string)

end

function split_libs(infile::String, barcodes::Array{String,1}, libnames::Array{String, 1}, 
    output_folder::String; barcode_length=-1, stop_early=-1, write_report=false)

    barcode_length < 0 && (barcode_length = length(barcodes[1]))
    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)

    output_files = [abspath(output_folder, "$(name).fastq.gz") for in_file in input_files]

    stats::Array{Int,1} = zeros(Int, length(input_files))
    record::FASTQ.Record = FASTQ.Record()
    reader = FASTQ.Reader(GzipDecompressorStream(open(infile, "r")))
    writers = [FASTQ.Writer(GzipCompressorStream(open(outfile, "w"), level=2)) for outfile in output_files] 

    sleep(0.01)

    c = 0
    while !eof(reader)
        read!(reader, record)
        ((c >= stop_early) & (stop_early > 0)) && break
        c += 1
        read = LongDNASeq(sequence(record))
        (library_id, nb_errors) = demultiplex(dplxr, read)
        nb_errors == -1 && continue
        stats[library_id] += 1
        write(writers[library_id], record)
    end

    for writer in writers
        close(writer)
    end
    close(reader)

    count_string = join(["$(names[i]) - $(stats_dict[i])\n" for i in 1:length(libnames)])
    count_string *= "not identifyable - $(c-sum(stats_dict))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write_report && write_file(abspath(output_folder, "demultiplexing_report.txt"), count_string)
end

function trim_fastp(input_files::Array{String, 1}, output_folder::String; 
    suffix="_trimmed", adapters::Array{Union{String,Nothing}, 1}=[], fastp_bin="fastp",
    min_length=25)

    record = FASTQ.Record()
    report = ""
    output_files = [split(in_file, ".fastq")[1]*"$(suffix).fastq.gz" for in_file in input_files]
    isempty(adapters) && (adapters = Array{String, 1}(nothing, length(input_files)))
    c = 0
    for (in_file, out_file, adapter) in zip(input_files, output_files, adapters)
        
        isnothing(adapter) ? 
        cmd = `$fastp_bin -i $in_file -o $out_file -M 25 -l $min_length --trim_poly_g 10 --cut_front --cut_tail` :
        cmd = `$fastp_bin -a $adapter -i $in_file -o $out_file -M 25 -l $min_length --trim_poly_g 10 --cut_front --cut_tail`
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

function trim_fastp(input_files::Array{Array{String, 1}, 1}, output_folder::String; 
    suffix="_trimmed", adapters::Array{Union{String,Nothing}, 1}=[], fastp_bin="fastp",
    min_length=25)

    record = FASTQ.Record()
    report = ""
    output_files = [[split(in_file1, ".fastq")[1]*"$(suffix).fastq.gz",
            split(in_file2, ".fastq")[1]*"$(suffix).fastq.gz"] for (in_file1, in_file2) in input_files]
    isempty(adapters) && (adapters = Array{String, 1}(nothing, length(input_files)))
    c = 0
    for ((in_file1, in_file2), (out_file1, out_file2), adapter) in zip(input_files, output_files, adapters)
        
        isnothing(adapter) ? 
        cmd = `$fastp_bin -i $in_file1 -o $out_file1 -I $in_file2 -O $out_file2 -M 25 -l $min_length --trim_poly_g 10 --cut_front --cut_tail` :
        cmd = `$fastp_bin -a $adapter $in_file1 -o $out_file1 -I $in_file2 -O $out_file2 -M 25 -l $min_length --trim_poly_g 10 --cut_front --cut_tail`
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