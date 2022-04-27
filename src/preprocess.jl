function download_sra(input_file::String; output_path=dirname(input_file), fastqdump_bin="fastq-dump", prefetch_bin="prefetch", keep_sra=false, overwrite_existing=false)
    open(input_file) do f
        for accession_number in eachline(f)
            sra_file = joinpath(output_path, "$accession_number.sra")
            fastqgz_file = sra_file[1:end-3] * "fastq.gz"
            (isfile(sra_file) || isfile(fastqgz_file)) && !overwrite_existing && continue
            isfile(sra_file) && rm(sra_file)
            isfile(fastqgz_file) && rm(fastqgz_file)
            run(`$prefetch_bin --output-directory $output_path $accession_number`)
            run(`$fastqdump_bin --gzip --outdir $output_path $sra_file`)
            keep_sra || rm(sra_file)
        end
    end
end
function download_sra(sra_ids::Vector{String}, output_path::String; fastqdump_bin="fastq-dump", prefetch_bin="prefetch", keep_sra=false, overwrite_existing=false)
    f = tempname()
    write(f, join(sra_ids, "\n"))
    download_sra(f; output_path=output_path, fastqdump_bin=fastqdump_bin, prefetch_bin=prefetch_bin, keep_sra=keep_sra, overwrite_existing=overwrite_existing)
    rm(f)
end

function split_libs(infile1::String, prefixfile1::Union{String,Nothing}, infile2::Union{String,Nothing}, libname_to_barcode::Dict{String,LongDNASeq}, output_path::String)

    dplxr = Demultiplexer(collect(values(libname_to_barcode)), n_max_errors=1, distance=:hamming)
    bc_len = [length(v) for v in values(libname_to_barcode)][1]
    output_files = isnothing(infile2) ?
    [joinpath(output_path, "$(name).fastq.gz") for name in keys(libname_to_barcode)] :
    [[joinpath(output_path, "$(name)_1.fastq.gz"), joinpath(output_path, "$(name)_2.fastq.gz")] for name in keys(libname_to_barcode)]
    for file in output_files
        isnothing(infile2) ?
        (isfile(file) && return) :
        ((isfile(file[1]) || isfile(file[2])) && return)
    end
    stats::Array{Int,1} = zeros(Int, length(libname_to_barcode))
    record1::FASTQ.Record = FASTQ.Record()
    isnothing(infile2) || (record2::FASTQ.Record = FASTQ.Record())
    isnothing(prefixfile1) || (recordp::FASTQ.Record = FASTQ.Record())
    endswith(infile1, ".gz") ? reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r"))) : reader1 = FASTQ.Reader(open(infile1, "r"))
    isnothing(infile2) || (endswith(infile2, ".gz") ?
                            reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r"))) :
                            reader2 = FASTQ.Reader(open(infile2, "r")))
    isnothing(prefixfile1) || (endswith(prefixfile1, ".gz") ?
                                readerp = FASTQ.Reader(GzipDecompressorStream(open(prefixfile1, "r"))) :
                                readerp = FASTQ.Reader(open(prefixfile1, "r")))
    writers = isnothing(infile2) ?
                [FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2))
                for outfile1 in output_files] :
                [[FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2)),
                FASTQ.Writer(GzipCompressorStream(open(outfile2, "w"), level=2))]
                for (outfile1, outfile2) in output_files]
    unidentified_writer = FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified.fastq.gz"), "w"), level=2))
    c = 0
    while !eof(reader1)
        read!(reader1, record1)
        isnothing(infile2) || read!(reader2, record2)
        isnothing(prefixfile1) || read!(readerp, recordp)
        c += 1
        nb_errors = -1
        isnothing(prefixfile1) || (prefix = LongDNASeq(FASTQ.sequence(recordp)); (library_id, nb_errors) = demultiplex(dplxr, prefix))
        if (nb_errors == -1)
            read = LongDNASeq(FASTQ.sequence(record1))
            (library_id, nb_errors) = demultiplex(dplxr, @view(read[1:bc_len]))
            (nb_errors == -1) && (write(unidentified_writer, record1); isnothing(infile2) || write(unidentified_writer, record2); continue)
        end
        stats[library_id] += 1
        isnothing(infile2) ?
        write(writers[library_id], record1) :
        (write(writers[library_id][1], record1); write(writers[library_id][2], record2))
    end
    close(reader1)
    isnothing(infile2) || close(reader2)
    close(unidentified_writer)
    for w in writers
        isnothing(infile2) ?
        close(w) :
        (close(w[1]);close(w[2]))
    end
    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libname_to_barcode, stats)])
    count_string *= "\nnot identifyable - $(c-sum(stats))\n"
    count_string = "Counted $c entries in total\n\n$count_string\n"
    write(infile1 * ".log", count_string)
end

function split_libs(infile1::String, infile2::String, libname_to_barcode::Dict{String,LongDNASeq}, output_path::String)
    split_libs(infile1, nothing, infile2, libname_to_barcode, output_path)
end

function split_libs(infile::String, libname_to_barcode::Dict{String,LongDNASeq}, output_path::String)
    split_libs(infile, nothing, nothing, libname_to_barcode, output_path)
end

function trim_fastp(input_files::Vector{Tuple{String, Union{String, Nothing}}};
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, umi_loc=:read1, min_length=nothing,
    max_length=nothing, cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=nothing, filter_complexity=nothing,
    average_window_quality=nothing, skip_quality_filtering=false, overwrite_existing=false)

    for (in_file1, in_file2) in input_files
        @assert endswith(in_file1, ".fasta") || endswith(in_file1, ".fasta.gz") || endswith(in_file1, ".fastq") || endswith(in_file1, ".fastq.gz")
        isnothing(in_file2) || (@assert endswith(in_file1, ".fasta") || endswith(in_file1, ".fasta.gz") || endswith(in_file1, ".fastq") || endswith(in_file1, ".fastq.gz"))
    end
    @assert umi_loc in [:read1, :read2]

    params = []
    skip_quality_filtering && push!(params, "--disable_quality_filtering")
    !skip_quality_filtering && !isnothing(average_window_quality) && push!(params, "--cut_mean_quality=$average_window_quality")
    !skip_quality_filtering && cut_tail && push!(params, "--cut_tail")
    !skip_quality_filtering && cut_front && push!(params, "--cut_front")

    !isnothing(filter_complexity) && append!(params, ["-y" "--complexity_threshold=$filter_complexity"])
    !isnothing(trim_poly_x) && append!(params, ["-x", "--poly_x_min_len=$trim_poly_x"])
    !isnothing(trim_poly_g) ? append!(params, ["-g", "--poly_g_min_len=$trim_poly_g"]) : push!(params,"-G")
    !isnothing(min_length) && push!(params, "--length_required=$min_length")
    !isnothing(max_length) && push!(params, "--max_len1=$max_length")
    !isnothing(umi) && append!(params, ["-U", "--umi_loc=$(String(umi_loc))", "--umi_len=$umi"])
    !isnothing(adapter) && push!(params, "--adapter_sequence=$adapter")

    for (in_file1, in_file2) in input_files
        startswith(basename(in_file1), prefix) && continue
        html_file = joinpath(dirname(in_file1), prefix * (endswith(in_file1, ".gz") ? basename(in_file1)[1:end-9] : basename(in_file1)[1:end-6]) * ".html")
        json_file = joinpath(dirname(in_file1), prefix * (endswith(in_file1, ".gz") ? basename(in_file1)[1:end-9] : basename(in_file1)[1:end-6]) * ".json")
        out_file1 = joinpath(dirname(in_file1), prefix * basename(in_file1))
        (isfile(out_file1) && !overwrite_existing) && continue
        out_file2 = !isnothing(in_file2) ? joinpath(dirname(in_file2), prefix * basename(in_file2)) : nothing
        file_params = ["--in1=$in_file1", "--out1=$out_file1", "--html=$(html_file)", "--json=$(json_file)"]
        !isnothing(in_file2) && append!(file_params, ["--in2=$in_file2", "--out2=$out_file2"])
        cmd = `$fastp_bin $file_params $params`
        run(cmd)
    end
end

function trim_fastp(input_files::SingleTypeFiles;
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, min_length=25, max_length=nothing,
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=25, skip_quality_filtering=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}([(file, nothing) for file in input_files])
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, umi=umi, umi_loc=:read1,
                min_length=min_length, max_length=max_length, cut_front=cut_front, cut_tail=cut_tail,
                trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity,
                average_window_quality=average_window_quality, skip_quality_filtering=skip_quality_filtering)
    return SingleTypeFiles([joinpath(dirname(file),prefix*basename(file)) for file in input_files if !startswith(basename(file), prefix)], input_files.type)
end

function trim_fastp(input_files::PairedSingleTypeFiles;
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, umi=nothing, umi_loc=:read1, min_length=25, max_length=nothing,
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=25, skip_quality_filtering=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}(input_files.list)
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, umi=umi, umi_loc=umi_loc,
                min_length=min_length, max_length=max_length, cut_front=cut_front, cut_tail=cut_tail,
                trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity,
                average_window_quality=average_window_quality, skip_quality_filtering=skip_quality_filtering)
    return PairedSingleTypeFiles([(joinpath(dirname(file1),prefix*basename(file1)), joinpath(dirname(file2),prefix*basename(file2))) for (file1, file2) in input_files if !startswith(basename(file1), prefix) | !startswith(basename(file2), prefix)], input_files.type, input_files.suffix1, input_files.suffix2)
end

function split_reads(file::String, split_at::Int; overwrite_existing=false)
    @assert any([endswith(file, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])

    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")

    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()

    out_file_base = file[1:end-(is_zipped ? 9 : 6)]
    out_file1 = out_file_base * "_1" * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
    out_file2 = out_file_base * "_2" * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
    (isfile(out_file1) && isfile(out_file2) && !overwrite_existing) && (return out_file1, out_file2)
    out_f1 = is_zipped ? GzipCompressorStream(open(out_file1, "w"); level=2) : open(out_file1, "w")
    out_f2 = is_zipped ? GzipCompressorStream(open(out_file2, "w"); level=2) : open(out_file2, "w")
    writer1 = is_fastq ? FASTQ.Writer(out_f1) : FASTA.Writer(out_f1)
    writer2 = is_fastq ? FASTQ.Writer(out_f2) : FASTA.Writer(out_f2)
    out_record1 = is_fastq ? FASTQ.Record() : FASTA.Record()
    out_record2 = is_fastq ? FASTQ.Record() : FASTA.Record()

    while !eof(reader)
        read!(reader, record)
        out_seq = is_fastq ? FASTQ.sequence(record) : FASTA.sequence(record)
        id = is_fastq ? FASTQ.identifier(record) : FASTA.identifier(record)
        q = is_fastq ? FASTQ.quality(record) : nothing
        out_record1 = is_fastq ? FASTQ.Record(id, out_seq[1:split_at], q[1:split_at]) : FASTA.Record(id, out_seq[1:split_at-1])
        out_record2 = is_fastq ? FASTQ.Record(id, out_seq[split_at+1:end], q[split_at+1:end]) : FASTA.Record(id, out_seq[split_at:end])
        write(writer1, out_record1)
        write(writer2, out_record2)
    end
    sleep(0.5)
    close(reader)
    close(writer1)
    close(writer2)
    return out_file1, out_file2
end

function transform(file::String; to_dna=false, reverse=false, complement=false, overwrite_existing=false, is_rna=false)
    @assert any([endswith(file, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])

    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")

    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()

    out_file = joinpath(dirname(file), "transformed_"*basename(file))
    (isfile(out_file) && !overwrite_existing) && (return out_file)
    out_f = is_zipped ? GzipCompressorStream(open(out_file, "w")) : open(out_file, "w")
    writer = is_fastq ? FASTQ.Writer(out_f) : FASTA.Writer(out_f)
    out_record = is_fastq ? FASTQ.Record() : FASTA.Record()

    in_seq = (to_dna || is_rna) ? LongRNASeq(0) : LongDNASeq(0)
    out_seq = (is_rna && !to_dna) ? LongRNASeq(0) : LongDNASeq(0)

    while !eof(reader)
        read!(reader, record)
        in_seq = is_fastq ? FASTQ.sequence((to_dna || is_rna) ? LongRNASeq : LongDNASeq, record) : FASTA.sequence((to_dna || is_rna) ? LongRNASeq : LongDNASeq, record)
        out_seq = to_dna ? LongDNASeq(in_seq) : in_seq
        reverse && reverse!(out_seq)
        complement && complement!(out_seq)
        out_record = is_fastq ? FASTQ.Record(FASTQ.identifier(record), out_seq, FASTQ.quality(record)) : FASTA.Record(FASTA.identifier(record), out_seq)
        write(writer, out_record)
    end
    sleep(0.5)
    close(reader)
    close(writer)
    return out_file
end

function transform(input_files::SingleTypeFiles; to_dna=false, reverse=false, complement=false, overwrite_existing=false, is_rna=false, stop_at=-1)
    transformed_files = String[]
    for file in input_files
        startswith(file, "trafo_") && continue
        push!(transformed_files, transform(file; to_dna=to_dna, reverse=reverse, complement=complement, overwrite_existing=overwrite_existing, is_rna=is_rna))
    end
    return SingleTypeFiles(transformed_files, input_files.type)
end