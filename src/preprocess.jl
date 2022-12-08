function download_prefetch(input_file::String; output_path=dirname(input_file), fastqdump_bin="fastq-dump", prefetch_bin="prefetch",
                            keep_sra=false, split_files=true, split_spot=false, overwrite_existing=false)
    sra_ids = String[]
    open(input_file) do f
        for sra_id in eachline(f)
            push!(sra_ids, sra_id)
        end
    end
    download_prefetch(sra_ids, output_path; fastqdump_bin=fastqdump_bin, prefetch_bin=prefetch_bin, keep_sra=keep_sra,
        split_files=split_files, split_spot=split_spot, overwrite_existing=overwrite_existing)
end
function download_prefetch(sra_runids::Vector{String}, output_path::String; fastqdump_bin="fastq-dump", prefetch_bin="prefetch",
                            keep_sra=false, split_files=true, split_spot=false, overwrite_existing=false)
    for accession_number in sra_runids
        sra_file = joinpath(output_path, "$accession_number.sra")
        fastqgz_file = sra_file[1:end-3] * "fastq.gz"
        isfile(sra_file) && overwrite_existing && rm(sra_file)
        isfile(fastqgz_file) && overwrite_existing && rm(fastqgz_file)
        if !isfile(sra_file)
            run(`$prefetch_bin --output-directory $output_path $accession_number`)
            if isdir(joinpath(output_path, accession_number))
                mv(joinpath(output_path, accession_number, "$accession_number.sra"), joinpath(output_path, "$accession_number.sra"))
                rm(joinpath(output_path, accession_number))
            end
        end
        params = String["--gzip"]
        split_files && push!(params, "--split-files")
        split_spot && push!(params, "--split-spot")
        !isfile(fastqgz_file) && run(`$fastqdump_bin $params --outdir $output_path $sra_file`)
        keep_sra || rm(sra_file)
    end
end

function demultiplex(testseq::LongDNA{4}, barcode_queries::Vector{ApproximateSearchQuery{typeof(isequal), LongDNA{4}}}; k=1)
    for i in 1:length(barcode_queries)
        isnothing(findfirst(barcode_queries[i], k, testseq)) || (return i)
    end
    return -1
end
function split_libs(infile1::String, prefixfile::Union{String,Nothing}, infile2::Union{String,Nothing}, libname_to_barcode::Dict{String,String}, output_path::String=dirname(infile1);
                        bc_len=length(first(values(libname_to_barcode))), check_range=1:bc_len, allowed_barcode_distance=1, overwrite_existing=false)

    lib_to_barcode = Dict(k=>LongDNA{4}(v) for (k,v) in libname_to_barcode)
    barcode_queries = [ApproximateSearchQuery(v) for v in values(lib_to_barcode)]
    for bc in values(lib_to_barcode)
        sum(!isnothing(findfirst(barcode_queries[i], 2*allowed_barcode_distance, bc)) for i in 1:length(lib_to_barcode)) > length(lib_to_barcode) &&
            throw(AssertionError("The supplied barcodes do not support the allowed_barcode_distance!"))
    end

    output_files = isnothing(infile2) ?
    [joinpath(output_path, "$(name).fastq.gz") for name in keys(lib_to_barcode)] :
    [(joinpath(output_path, "$(name)_1.fastq.gz"), joinpath(output_path, "$(name)_2.fastq.gz")) for name in keys(lib_to_barcode)]
    for file in output_files
        if isnothing(infile2)
            if isfile(file)
                overwrite_existing || throw(AssertionError("File $file already exists."))
                rm(file)
            end
        else
            if isfile(file[1])
                overwrite_existing || throw(AssertionError("File $(file[1]) already exists."))
                rm(file[1])
            end
            if isfile(file[2])
                overwrite_existing || throw(AssertionError("File $(file[2]) already exists."))
                rm(file[2])
            end
        end
    end
    nb_stats = length(lib_to_barcode)+1
    stats::Vector{Int} = zeros(Int, nb_stats)
    record1::FASTQ.Record = FASTQ.Record()
    isnothing(infile2) || (record2::FASTQ.Record = FASTQ.Record())
    isnothing(prefixfile) || (recordp::FASTQ.Record = FASTQ.Record())
    endswith(infile1, ".gz") ? reader1 = FASTQ.Reader(GzipDecompressorStream(open(infile1, "r"))) : reader1 = FASTQ.Reader(open(infile1, "r"))
    isnothing(infile2) || (endswith(infile2, ".gz") ?
                            reader2 = FASTQ.Reader(GzipDecompressorStream(open(infile2, "r"))) :
                            reader2 = FASTQ.Reader(open(infile2, "r")))
    isnothing(prefixfile) || (endswith(prefixfile, ".gz") ?
                                readerp = FASTQ.Reader(GzipDecompressorStream(open(prefixfile, "r"))) :
                                readerp = FASTQ.Reader(open(prefixfile, "r")))
    writers = isnothing(infile2) ?
                [FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2))
                for outfile1 in output_files] :
                [[FASTQ.Writer(GzipCompressorStream(open(outfile1, "w"), level=2)),
                FASTQ.Writer(GzipCompressorStream(open(outfile2, "w"), level=2))]
                for (outfile1, outfile2) in output_files]
    isnothing(infile2) ? push!(writers, FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified.fastq.gz"), "w"), level=2))) :
        push!(writers, [FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified_1.fastq.gz"), "w"), level=2)),
                        FASTQ.Writer(GzipCompressorStream(open(joinpath(output_path, "unidentified_2.fastq.gz"), "w"), level=2))])
    c = 0
    while !eof(reader1) && !eof(reader2)
        read!(reader1, record1)
        isnothing(infile2) || read!(reader2, record2)
        isnothing(prefixfile) || read!(readerp, recordp)
        c += 1
        checkseq = isnothing(prefixfile) ? LongDNA{4}(view(record1.data,record1.sequence))[check_range] : LongDNA{4}(view(recordp.data,recordp.sequence))
        library_id = demultiplex(checkseq, barcode_queries; k=allowed_barcode_distance)
        library_id == -1 && (library_id = nb_stats)
        stats[library_id] += 1
        isnothing(infile2) ?
        write(writers[library_id], record1) :
        (write(writers[library_id][1], record1); write(writers[library_id][2], record2))
    end
    close(reader1)
    isnothing(infile2) || close(reader2)
    for w in writers
        isnothing(infile2) ?
        close(w) :
        (close(w[1]);close(w[2]))
    end
    libnames = String.(keys(lib_to_barcode))
    sorted_index = sortperm(libnames)
    count_string = join(["$(name) - $(stat)\n" for (name, stat) in zip(libnames[sorted_index], stats[sorted_index])])
    count_string *= "\nnot identifyable - $(stats[end])\n"
    count_string = "Counted $c entries in total:\n\n$count_string\n"
    write(infile1 * ".log", count_string)
    return isnothing(infile2) ? FastqgzFiles(output_files) : PairedSingleTypeFiles(output_files, ".fastq.gz", "_1", "_2")
end

function split_libs(infile1::String, infile2::String, libname_to_barcode::Dict{String,String}, output_path::String=dirname(infile1); bc_len=8, check_range=1:bc_len, overwrite_existing=false)
    split_libs(infile1, nothing, infile2, libname_to_barcode, output_path; bc_len=bc_len, check_range=check_range, overwrite_existing=overwrite_existing)
end

function split_libs(infile::String, libname_to_barcode::Dict{String,String}, output_path::String=dirname(infile); bc_len=8, check_range=1:bc_len, overwrite_existing=false)
    split_libs(infile, nothing, nothing, libname_to_barcode, output_path; bc_len=bc_len, check_range=check_range, overwrite_existing=overwrite_existing)
end

function trim_fastp(input_files::Vector{Tuple{String, Union{String, Nothing}}};
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, trim=nothing, trim_loc=:read1, min_length=nothing,
    max_length=nothing, cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=nothing, filter_complexity=nothing,
    average_window_quality=nothing, deduplicate=false, skip_quality_filtering=false, overwrite_existing=false)

    @assert trim_loc in [:read1, :read2]

    params = []
    !isnothing(trim) && push!(params, trim_loc === :read1 ? "--trim_front1=$trim" : "--trim_front2=$trim")
    skip_quality_filtering && push!(params, "--disable_quality_filtering")
    !skip_quality_filtering && !isnothing(average_window_quality) && push!(params, "--cut_mean_quality=$average_window_quality")
    !skip_quality_filtering && cut_tail && push!(params, "--cut_tail")
    !skip_quality_filtering && cut_front && push!(params, "--cut_front")

    deduplicate && push!(params, "--dedup")
    !isnothing(filter_complexity) && append!(params, ["-y" "--complexity_threshold=$filter_complexity"])
    !isnothing(trim_poly_x) && append!(params, ["-x", "--poly_x_min_len=$trim_poly_x"])
    !isnothing(trim_poly_g) ? append!(params, ["-g", "--poly_g_min_len=$trim_poly_g"]) : push!(params,"-G")
    !isnothing(min_length) && push!(params, "--length_required=$min_length")
    !isnothing(max_length) && push!(params, "--max_len1=$max_length")
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
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, trim=nothing, trim_loc=nothing, min_length=20,
    max_length=nothing, cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=20, deduplicate=false, skip_quality_filtering=false, overwrite_existing=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}([(file, nothing) for file in input_files])
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, trim=trim, trim_loc=:read1,
                min_length=min_length, max_length=max_length, cut_front=cut_front, cut_tail=cut_tail,
                trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity,
                average_window_quality=average_window_quality, deduplicate=deduplicate,
                skip_quality_filtering=skip_quality_filtering, overwrite_existing=overwrite_existing)
    return SingleTypeFiles([joinpath(dirname(file),prefix*basename(file)) for file in input_files if !startswith(basename(file), prefix)], input_files.type)
end

function trim_fastp(input_files::PairedSingleTypeFiles;
    fastp_bin="fastp", prefix="trimmed_", adapter=nothing, trim=nothing, trim_loc=:read1, min_length=20,
    max_length=nothing, cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=20, deduplicate=false, skip_quality_filtering=false, overwrite_existing=false)

    files = Vector{Tuple{String, Union{String, Nothing}}}(input_files.list)
    trim_fastp(files; fastp_bin=fastp_bin, prefix=prefix, adapter=adapter, trim=trim, trim_loc=trim_loc,
                min_length=min_length, max_length=max_length, cut_front=cut_front, cut_tail=cut_tail,
                trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x, filter_complexity=filter_complexity,
                average_window_quality=average_window_quality, deduplicate=deduplicate,
                skip_quality_filtering=skip_quality_filtering, overwrite_existing=overwrite_existing)
    return PairedSingleTypeFiles([(joinpath(dirname(file1),prefix*basename(file1)), joinpath(dirname(file2),prefix*basename(file2)))
        for (file1, file2) in input_files if !startswith(basename(file1), prefix) | !startswith(basename(file2), prefix)], input_files.type, input_files.suffix1, input_files.suffix2)
end

function split_each_read(file::String, split_after::Int; suffix1="_1", suffix2="_2", overwrite_existing=false)
    is_fastq = any([endswith(file, ending) for ending in FASTQ_TYPES])
    is_zipped = endswith(file, ".gz")

    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()

    out_file_base = file[1:end-(is_zipped ? 9 : 6)]
    out_file1 = out_file_base * suffix1 * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
    out_file2 = out_file_base * suffix2 * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
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
        out_record1 = is_fastq ? FASTQ.Record(id, out_seq[1:split_after], q[1:split_after]) : FASTA.Record(id, out_seq[1:split_after-1])
        out_record2 = is_fastq ? FASTQ.Record(id, out_seq[split_after+1:end], q[split_after+1:end]) : FASTA.Record(id, out_seq[split_after:end])
        write(writer1, out_record1)
        write(writer2, out_record2)
    end
    sleep(0.5)
    close(reader)
    close(writer1)
    close(writer2)
    return out_file1, out_file2
end
function split_each_read(files::SingleTypeFiles, split_after::Int; suffix1="_1", suffix2="_2", overwrite_existing=false)
    processed_files = Vector{Tuple{String,Sting}}()
    for file in files
        push!(processed_files, split_each_read(file, split_after; suffix1=suffix1, suffix2=suffix2, overwrite_existing=overwrite_existing))
    end
    return PairedSingleTypeFiles(processed_files)
end

function split_interleaved(file::String; suffix1="_1", suffix2="_2", force_same_name=true, overwrite_existing=false)
    is_fastq = any([endswith(file, ending) for ending in FASTQ_TYPES])
    is_zipped = endswith(file, ".gz")

    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()

    out_file_base = file[1:end-(is_zipped ? 9 : 6)]
    out_file1 = out_file_base * suffix1 * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
    out_file2 = out_file_base * suffix2 * (is_fastq ? ".fastq" : ".fasta") * (is_zipped ? ".gz" : "")
    (isfile(out_file1) && isfile(out_file2) && !overwrite_existing) && (return out_file1, out_file2)
    out_f1 = is_zipped ? GzipCompressorStream(open(out_file1, "w"); level=2) : open(out_file1, "w")
    out_f2 = is_zipped ? GzipCompressorStream(open(out_file2, "w"); level=2) : open(out_file2, "w")
    writer1 = is_fastq ? FASTQ.Writer(out_f1) : FASTA.Writer(out_f1)
    writer2 = is_fastq ? FASTQ.Writer(out_f2) : FASTA.Writer(out_f2)
    namebuffer = UInt8[]
    while !eof(reader)
        read!(reader, record)
        write(writer1, record)
        force_same_name && copy!(namebuffer, view(record.data, record.identifier))
        read!(reader, record)
        force_same_name && copy!(view(record.data, record.identifier), namebuffer)
        write(writer2, record)
    end
    sleep(0.5)
    close(reader)
    close(writer1)
    close(writer2)
    return out_file1, out_file2
end
function split_interleaved(files::SingleTypeFiles; suffix1="_1", suffix2="_2", force_same_name=true, overwrite_existing=false)
    processed_files = Vector{Tuple{String,String}}()
    for file in files
        push!(processed_files, split_interleaved(file; suffix1=suffix1, suffix2=suffix2, force_same_name=force_same_name, overwrite_existing=overwrite_existing))
    end
    return PairedSingleTypeFiles(processed_files)
end

function transform(file::String; prefix="trafo_", to_dna=false, reverse=false, complement=false, overwrite_existing=false, is_rna=false)
    is_fastq = any([endswith(file, ending) for ending in FASTQ_TYPES])
    is_zipped = endswith(file, ".gz")

    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()

    out_file = joinpath(dirname(file), prefix*basename(file))
    (isfile(out_file) && !overwrite_existing) && (return out_file)
    out_f = is_zipped ? GzipCompressorStream(open(out_file, "w")) : open(out_file, "w")
    writer = is_fastq ? FASTQ.Writer(out_f) : FASTA.Writer(out_f)
    out_record = is_fastq ? FASTQ.Record() : FASTA.Record()

    in_seq = (to_dna || is_rna) ? LongRNASeq(0) : LongDNA{4}(0)
    out_seq = (is_rna && !to_dna) ? LongRNASeq(0) : LongDNA{4}(0)
    while !eof(reader)
        read!(reader, record)
        in_seq = is_fastq ? FASTQ.sequence((to_dna || is_rna) ? LongRNASeq : LongDNA, record) : FASTA.sequence((to_dna || is_rna) ? LongRNASeq : LongDNA, record)
        out_seq = to_dna ? LongDNA{4}(in_seq) : in_seq
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

function transform(input_files::SingleTypeFiles; prefix="trafo_", to_dna=false, reverse=false, complement=false, overwrite_existing=false, is_rna=false)
    transformed_files = String[]
    for file in input_files
        startswith(file, prefix) && continue
        push!(transformed_files, transform(file; to_dna=to_dna, reverse=reverse, complement=complement, overwrite_existing=overwrite_existing, is_rna=is_rna))
    end
    return SingleTypeFiles(transformed_files, input_files.type)
end

"""
    Wrapper function for kraken2 taxonomic sequence classifier that assigns taxonomic labels to DNA sequences.

    # Output
    *.results.txt is a five fields tab-delimited text file.
    *.report.txt is in kraken2's modified report format (--report-minimizer-data).
    Report file is used for KronaTools plot.
"""
function align_kraken2(sequence_file::String, db_path::String;
    kraken_bin = "kraken2", threads = 6, report = true, results = false, quick = false, min_hit_groups = false, )

    output_file = split(sequence_file, ".")[1] * ".kraken2_results.txt"
    report_file = split(sequence_file, ".")[1] * ".report.txt"

    params = ["--db", db_path,
              "--threads", threads,
              sequence_file,
             ]

    # append additional options
    report && append!(params, ["--report", report_file, "--report-minimizer-data"])
    quick && push!(params, "--quick")
    min_hit_groups && push!(params, "--minimum-hit-groups")
    # if results file is unwanted send stdout to /dev/null
    !(results) && (output_file = devnull)

    run(pipeline(`$kraken_bin $params`, stdout=output_file))
end

function fnamefromseqfile(seqfilename::String, newtype::String=".bam")
    for ending in Iterators.flatten((FASTA_TYPES, FASTQ_TYPES))
        endswith(seqfilename, ending) && return seqfilename[1:end-length(ending)] * newtype
    end
    return seqfilename * newtype
end
function fnamefromseqfile(seqfilename1::String, seqfilename2::Union{Nothing,String}, newtype::String=".bam")
    isnothing(seqfilename2) && return fnamefromseqfile(seqfilename1, newtype)
    for ending in Iterators.flatten((FASTA_TYPES, FASTQ_TYPES))
        if endswith(seqfilename1, ending) && endswith(seqfilename2, ending)
            for (ii, (i,j)) in enumerate(zip(basename(seqfilename1), basename(seqfilename2)))
                ii == 1 && return seqfilename1[1:end-length(ending)] * newtype
                i != j && return joinpath(dirname(seqfilename1), basename(seqfilename1)[1:ii-1] * newtype)
            end
        end
    end
    return seqfilename1 * newtype
end

"""
    Core dispatch of align_minimap, which is a wrapper for minimap2. Takes ´in_file´ and ´genome_file´ and builds a shell
    command that runs minimap2 with the specified additional parameters and writes the resulting alignments in bam format into ´out_file´
    using samtools.

    # Arguments
    - ´min_score::Int´: defines the minimum score for bwa-mem2 to output alignments for
    - ´match::Int´, ´mismatch::Int´, ´gap_open::Int´, ´gamp_extend::Int´: define the affine gap model that is used for scoring
    - ´clipping_penalty::Int´: scoring penalty for clipping at the beginning or the end of the alignment.
    - ´unpair_penalty::Int´: scoring penalty only used in paired end mapping. Penalty for mappings that do not set the paired flag
    - ´unpair_rescue::Bool´: perform Smith-Waterman to try to rescue read pair if pair is lost.
    - ´min_seed_len::Int´: Minimum seed length. Matches shorter than INT will be missed.
    - ´reseeding_factor::Float64´: Trigger re-seeding for a MEM longer than minSeedLenFLOAT. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy.
"""
function align_minimap(in_file::String, genome_file::String, out_file::String=fnamefromseqfile(in_file);
    preset=nothing, match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_extend1=2, gap_extend2=1, minimizer_len=15,
    threads=6, all_secondary=false, skip_stats=false, sort_bam=true, minimap_bin="minimap2", sam_bin="samtools", overwrite_existing=false)

    !overwrite_existing && isfile(out_file) && throw(AssertinError("File $out_file already exists!"))
    !isnothing(preset) && !(preset in (:mapont, :asm5, :asm10, :asm20, :sr)) &&
        throw(AssertionError("Allowed preset values are: :map-ont, :asm5, :asm10, :asm20 and :sr"))

    params = isnothing(preset) ?
        Any["-a", "-A", match, "-B", mismatch, "-O", "$gap_open1,$gap_open2", "-E", "$gap_extend1,$gap_extend2", "-k", minimizer_len] :
        Any["-ax", preset]

    append!(params, ["-t", threads, "--eqx"])
    all_secondary && push!(params, "-P")
    fileparams = [genome_file, in_file]

    stats_file = out_file * ".log"
    run(pipeline(
        `$minimap_bin $params $fileparams`,
        stdout = pipeline(
            `$sam_bin view -u`,
            stdout = sort_bam ? pipeline(`$sam_bin sort -o $out_file`) : out_file)))
    skip_stats ?
    run(`$sam_bin index $out_file`) :
    run(pipeline(
        `$sam_bin index $out_file`,
        `$sam_bin stats $out_file`,
        stats_file))
end

"""
    Helper dispatch of align_minimap. Runs align_minimap on `read_files::SingleTypeFiles` against `genome::Genome`.
    This enables easy handling of whole folders containing sequence files and the manipulation of a genome using Genome.
"""
function align_minimap(sequence_files::SingleTypeFiles, genome::Genome;
    preset=nothing, match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_extend1=2, gap_extend2=1,
    minimizer_len=15, threads=6, all_secondary=false, sort_bam=true, minimap_bin="minimap2", sam_bin="samtools", overwrite_existing=false)
    sequence_files.type in Iterators.flatten((FASTA_TYPES, FASTQ_TYPES)) || throw(AssertionError("Unknown file type."))
    tmp_genome = tempname()
    write(tmp_genome, genome)
    outfiles = Vector{String}()
    for file in sequence_files
        out_file = isempty(sequence_files.type) ? fnamefromseqfile(file) : file[1:end-length(sequence_files.type)] * ".bam"
        push!(outfiles, out_file)
        (isfile(out_file) && !overwrite_existing) && continue
        align_minimap(file, out_file, tmp_genome;
            preset=preset, match=match, mismatch=mismatch, gap_open1=gap_open1, gap_open2=gap_open2,
            gap_extend1=gap_extend1, gap_extend2=gap_extend2, minimizer_len=minimizer_len, threads=threads, all_secondary=all_secondary,
            sort_bam=sort_bam, minimap_bin=minimap_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    return SingleTypeFiles(outfiles)
end

"""
    Core dispatch of align_mem, which is a wrapper for bwa-mem2. Takes ´in_file1´ and ´in_file2´ and ´genome_file´ and builds a shell
    command that runs bwa-mem2 with the specified additional parameters and writes the resulting alignments in bam format into ´out_file´
    using samtools.

    # Arguments
    - ´min_score::Int´: defines the minimum score for bwa-mem2 to output alignments for
    - ´match::Int´, ´mismatch::Int´, ´gap_open::Int´, ´gamp_extend::Int´: define the affine gap model that is used for scoring
    - ´clipping_penalty::Int´: scoring penalty for clipping at the beginning or the end of the alignment.
    - ´unpair_penalty::Int´: scoring penalty only used in paired end mapping. Penalty for mappings that do not set the paired flag
    - ´unpair_rescue::Bool´: perform Smith-Waterman to try to rescue read pair if pair is lost.
    - ´min_seed_len::Int´: Minimum seed length. Matches shorter than INT will be missed.
    - ´reseeding_factor::Float64´: Trigger re-seeding for a MEM longer than minSeedLenFLOAT. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy.
"""
function align_mem(in_file1::String, in_file2::Union{String,Nothing}, genome_file::String, out_file::String=fnamefromseqfile(in_file1, in_file2);
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false,
    min_seed_len=18, reseeding_factor=1.4, is_ont=false, threads=6, sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    params = ["-A", match, "-B", mismatch, "-O", gap_open, "-E", gap_extend, "-T", min_score, "-L", clipping_penalty, "-r", reseeding_factor, "-k", min_seed_len, "-t", threads]
    isnothing(in_file2) || append!(params, ["-U", unpair_penalty])
    is_ont && append!(params, ["-x", "ont2d"])
    unpair_rescue && push!(params, "-P")
    fileparams = isnothing(in_file2) ? [genome_file, in_file1] : [genome_file, in_file1, in_file2]
    stats_file = out_file * ".log"
    run(pipeline(
        `$bwa_bin mem -v 1 $params $fileparams`,
        stdout = pipeline(
            `$sam_bin view -u`,
            stdout = sort_bam ? pipeline(`$sam_bin sort -o $out_file`) : out_file)))
    run(pipeline(
        `$sam_bin index $out_file`,
        `$sam_bin stats $out_file`,
        stats_file))
end
align_mem(in_file::String, genome_file::String, out_file::String=fnamefromseqfile(in_file);
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, min_seed_len=18,
    reseeding_factor=1.4, is_ont=false, threads=6, sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools") =
    align_mem(in_file, nothing, genome_file, out_file;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, clipping_penalty=clipping_penalty, min_seed_len=min_seed_len,
        reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads, sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin)

"""
    Helper dispatch of align_mem, which is a wrapper for bwa-mem2. Runs align_mem on `read_files::T`
    where T is a FileCollection and aligns it against `genome::Genome`. This enables easy handling of
    whole folders containing sequence files and the manipulation of a genome using Genome.
"""
function align_mem(read_files::T, genome::Genome; min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1,
    clipping_penalty=5, unpair_penalty=9, reseeding_factor=1.4, min_seed_len=18, unpair_rescue=false, is_ont=false,
    threads=6, sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:FileCollection}

    tmp_genome = tempname()
    write(tmp_genome, genome)
    outfiles = Vector{String}()
    for file in read_files
        out_file = isa(read_files, SingleTypeFiles) ? file[1:end-length(read_files.type)] * ".bam" :
                    first(file)[1:end-length(read_files.type)-length(read_files.suffix1)] * ".bam"
        push!(outfiles, out_file)
        (isfile(out_file) && !overwrite_existing) && continue
        isa(read_files, SingleTypeFiles) ?
        align_mem(file, tmp_genome, out_file;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open,
                gap_extend=gap_extend, clipping_penalty=clipping_penalty,  min_seed_len=min_seed_len,
                reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads, sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(first(file), last(file), tmp_genome, out_file;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
                clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
                min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
                sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        isfile(tmp_genome * ending) && rm(tmp_genome * ending)
    end
    return SingleTypeFiles(outfiles)
end

"""
    Helper dispatch of align_mem, which is a wrapper for bwa-mem2. Runs align_mem on `reads::T` where T
    is a SequenceContainer and aligns it against `genomes::Vector{Genome}`. This enables easy handling
    of Sequences and their alignment against multiple genomes.
"""
function align_mem(reads::Sequences, genomes::Vector{Genome}, out_file::String="alignments.bam"; min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1,
    clipping_penalty=5, unpair_penalty=9, unpair_rescue=false, min_seed_len=18, reseeding_factor=1.4, is_ont=false, threads=6,
    sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false)

    (isfile(out_file) && !overwrite_existing) && return
    tmp_reads = tempname()
    tmp_reads2 = tempname()
    tmp_genome = tempname()
    ispaired = reads.tempnames[1:2:end] == reads.tempnames[2:2:end]
    if ispaired
        write(tmp_reads, tmp_reads2, reads)
    else
        write(tmp_reads, reads)
    end
    for (i,genome) in enumerate(genomes)
        write(tmp_genome, genome)
        this_out_file = out_file
        length(genomes) > 1 && (this_out_file = joinpath(dirname(out_file), "$(i)_" * basename(out_file)))

        !ispaired ?
        align_mem(tmp_reads, tmp_genome, this_out_file;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor,
            is_ont=is_ont, threads=threads, sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(tmp_reads, tmp_reads2, tmp_genome, this_out_file;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
            min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
            sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin)

        rm(tmp_genome)
    end
    rm(tmp_reads)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
    if length(genomes) > 1
        bam_files = []
        bai_files = []
        for i in 1:length(genomes)
            push!(bam_files, joinpath(dirname(out_file), "$(i)_" * basename(out_file)))
            push!(bai_files, joinpath(dirname(out_file), "$(i)_" * basename(out_file) * ".bai"))
        end
        run(`$sam_bin merge -X $out_file $bam_files $bai_files`)
        run(`$sam_bin index $out_file`)
        for (bam_file, bai_file) in zip(bam_files, bai_files)
            rm(bam_file)
            rm(bai_file)
        end
    end
end
align_mem(reads::T, genome::Genome, out_file::String;
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, unpair_penalty=9, min_seed_len=18, reseeding_factor=1.4,
    unpair_rescue=false, is_ont=false, threads=6, sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:Sequences} =
    align_mem(reads, [genome], out_file;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, unpair_penalty=unpair_penalty,
        unpair_rescue=unpair_rescue, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, threads=threads,
        sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin, overwrite_existing=overwrite_existing)

function preprocess_data(files::Union{SingleTypeFiles, PairedSingleTypeFiles}, genome::Genome;
    fastp_bin="fastp", trimmed_prefix="trimmed_", adapter=nothing, trim=nothing, trim_loc=:read1, min_length=20, max_length=nothing,
    cut_front=true, cut_tail=true, trim_poly_g=nothing, trim_poly_x=10, filter_complexity=nothing,
    average_window_quality=20, deduplicate=false, skip_quality_filtering=false,
    min_score=20, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false,
    min_seed_len=18, reseeding_factor=1.4, is_ont=false, threads=6, sort_bam=true, bwa_bin="bwa-mem2", sam_bin="samtools",
    norm=1000000, include_secondary_alignments=true, suffix_forward="_forward", suffix_reverse="_reverse",
    is_reverse_complement=false, overwrite_existing=false)

    println("Preprocessing files:")
    show(files)
    trimmed = trim_fastp(files; overwrite_existing=overwrite_existing, fastp_bin=fastp_bin, prefix=trimmed_prefix, adapter=adapter, trim=trim, trim_loc=trim_loc,
                        min_length=min_length, max_length=max_length, cut_front=cut_front, cut_tail=cut_tail, trim_poly_g=trim_poly_g, trim_poly_x=trim_poly_x,
                        filter_complexity=filter_complexity, average_window_quality=average_window_quality, deduplicate=deduplicate, skip_quality_filtering=skip_quality_filtering)
    println("Aligning files...")
    bams = align_mem(trimmed, genome; overwrite_existing=overwrite_existing, min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
                        clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor,
                        is_ont=is_ont, threads=threads, sort_bam=sort_bam, bwa_bin=bwa_bin, sam_bin=sam_bin)
    println("Computing coverage...")
    compute_coverage(bams; overwrite_existing=overwrite_existing, norm=norm, is_reverse_complement=is_reverse_complement, include_secondary_alignments=include_secondary_alignments,
                        suffix_forward=suffix_forward, suffix_reverse=suffix_reverse);
    println("Done.")
end