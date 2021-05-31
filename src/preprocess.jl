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

function trim_fastp(input_files::Vector{Tuple{String, Union{String, Nothing}}}; 
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, umi_loc=:read1, min_length=nothing, 
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=nothing, filter_complexity=nothing,
    average_window_quality=nothing, skip_quality_filtering=false, overwrite_existing=false)

    for (in_file1, in_file2) in input_files
        @assert endswith(in_file1, ".fasta") || endswith(in_file1, ".fasta.gz") || endswith(in_file1, ".fastq") || endswith(in_file1, ".fastq.gz")
        isnothing(in_file2) || (@assert endswith(in_file1, ".fasta") || endswith(in_file1, ".fasta.gz") || endswith(in_file1, ".fastq") || endswith(in_file1, ".fastq.gz"))
    end
    @assert umi_loc in [:read1, :read2]

    params = []
    skip_quality_filtering && push!(params, "--disable_quality_filtering")
    !skip_quality_filtering && !isnothing(average_window_quality) && push!(params, "--cut_mean_quality=$average_window_quality")
    !skip_quality_filtering && !isnothing(filter_complexity) && append!(params, ["-y" "--complexity_threshold=$filter_complexity"])
    !skip_quality_filtering && !isnothing(trim_poly_x) && append!(params, ["-x", "--poly_x_min_len=$trim_poly_x"])
    !skip_quality_filtering && !isnothing(trim_poly_g) ? append!(params, ["-g", "--poly_g_min_len=$trim_poly_g"]) : push!(params,"-G")
    !skip_quality_filtering && cut_tail && push!(params, "--cut_tail")
    !skip_quality_filtering && cut_front && push!(params, "--cut_front")
    !isnothing(min_length) && push!(params, "--length_required=$min_length")
    !isnothing(umi) && append!(params, ["-U", "--umi_loc=$(String(umi_loc))", "--umi_len=$umi"])
    !isnothing(adapter) && push!(params, "--adapter_sequence=$adapter")

    for (in_file1, in_file2) in input_files
        startswith(basename(in_file1), prefix) && continue
        html_file = joinpath(dirname(in_file1), prefix * (endswith(in_file1, ".gz") ? basename(in_file1)[1:end-9] : basename(in_file1)[1:end-6]) * ".html")
        json_file = joinpath(dirname(in_file1), prefix * (endswith(in_file1, ".gz") ? basename(in_file1)[1:end-9] : basename(in_file1)[1:end-6]) * ".json")
        out_file1 = joinpath(dirname(in_file1), prefix * basename(in_file1))
        (isfile(out_file1) && !overwrite_existing) && continue
        out_file2 = !isnothing(in_file2) ? joinpath(dirname(in_file2), prefix * basename(in_file2)) : nothing
        push!(out_files, (out_file1, out_file2))
        file_params = ["--in1=$in_file1", "--out1=$out_file1", "--html=$(html_file)", "--json=$(json_file)"]
        !isnothing(in_file2) && append!(file_params, ["--in2=$in_file2", "--out2=$out_file2"])
        cmd = `$fastp_bin $file_params $params`
        run(cmd)
    end 

end

function trim_fastp(input_files::SingleTypeFiles; 
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, min_length=25, 
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=25, skip_quality_filtering=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}([(file, nothing) for file in input_files])
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, umi=umi, umi_loc=:read1, min_length=min_length,
        cut_front=cut_front, cut_tail=cut_tail, trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity, average_window_quality=average_window_quality, skip_quality_filtering=skip_quality_filtering)
    return SingleTypeFiles([joinpath(dirname(file),prefix*basename(file)) for file in input_files if !startswith(basename(file), prefix)], input_files.type)
end

function trim_fastp(input_files::PairedSingleTypeFiles; 
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, umi_loc=:read1, min_length=25, 
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=25, skip_quality_filtering=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}(input_files.list)
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, umi=umi, umi_loc=umi_loc, min_length=min_length, cut_front=cut_front,
        cut_tail=cut_tail, trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity, average_window_quality=average_window_quality, skip_quality_filtering=skip_quality_filtering)
    return PairedSingleTypeFiles([(joinpath(dirname(file1),prefix*basename(file1)), joinpath(dirname(file2),prefix*basename(file2))) for (file1, file2) in input_files if !startswith(basename(file1), prefix) | !startswith(basename(file2), prefix)], input_files.type, input_files.suffix1, input_files.suffix2)
end