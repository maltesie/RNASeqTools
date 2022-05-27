#!/usr/bin/env julia
"""
    Wrapper function for kraken2 taxonomic sequence classifier that assigns taxonomic labels to DNA sequences.

    # Output
    *.results.txt is a five fields tab-delimited text file.
    *.report.txt is in kraken2's modified report format (--report-minimizer-data).
    Report file is used for KronaTools plot.
"""
function align_kraken2(
        db_location::String, sequence_file::String;
        kraken_bin = "kraken2",
        threads = 6, report = true, results = false,
        quick = false, min_hit_groups = false, )

    output_file = split(sequence_file, ".")[1] * ".kraken2_results.txt"
    report_file = split(sequence_file, ".")[1] * ".report.txt"

    params = ["--db", db_location,
              sequence_file,
              "--threads", threads,
             ]

    # append additional options
    report && append!(params, ["--report", report_file, "--report-minimizer-data"])
    quick && push!(params, "--quick")
    min_hit_groups && push!(params, "--minimum-hit-groups")
    # if results file is unwanted send stdout to /dev/null
    !(results) && (output_file = devnull)

    cmd = pipeline(`$kraken_bin $params`, stdout=output_file)
    run(cmd)
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
function align_minimap(in_file::String, out_file::String, genome_file::String;
    preset=nothing, match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_extend1=2, gap_extend2=1,
    minimizer_len=15, threads=6, all_secondary=false, skip_stats=false, minimap_bin="minimap2", sam_bin="samtools")

    !isnothing(preset) && !(preset in ("map-ont", "asm5", "asm10", "asm20", "sr")) &&
        throw(AssertionError("Allowed preset values are: map-ont, asm5, asm10, asm20"))

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
            stdout = pipeline(
                `$sam_bin sort -o $out_file`))))
    skip_stats ?
    run(`$sam_bin index $out_file`) :
    run(pipeline(
        `$sam_bin index $out_file`,
        `$sam_bin stats $out_file`,
        stats_file))
end
align_minimap(in_file::String, genome_file::String;
    preset=nothing, match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_extend1=2, gap_extend2=1,
    minimizer_len=15, threads=6, all_secondary=false, skip_stats=false, minimap_bin="minimap2", sam_bin="samtools") =
    align_minimap(in_file, joinpath(in_file[1:findlast('.', in_file)] * "bam"), genome_file;
    preset=preset, match=match, mismatch=mismatch, gap_open1=gap_open1, gap_open2=gap_open2,
    gap_extend1=gap_extend1, gap_extend2=gap_extend2, minimizer_len=minimizer_len, threads=threads, all_secondary=all_secondary,
    skip_stats=skip_stats, minimap_bin=minimap_bin, sam_bin=sam_bin)

"""
    Helper dispatch of align_minimap. Runs align_minimap on `read_files::SingleTypeFiles` against `genome::Genome`.
    This enables easy handling of whole folders containing sequence files and the manipulation of a genome using Genome.
"""
function align_minimap(sequence_files::SingleTypeFiles, genome::Genome;
    preset=nothing, match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_extend1=2, gap_extend2=1,
    minimizer_len=15, threads=6, all_secondary=false, minimap_bin="minimap2", sam_bin="samtools", overwrite_existing=false)
    sequence_files.type in (".fastq.gz", ".fastq", ".fasta.gz", ".fasta", ".fa", ".fna") || throw(AssertionError("Unknown file type."))
    tmp_genome = tempname()
    write(tmp_genome, genome)
    outfiles = Vector{String}()
    for file in sequence_files
        out_file = isa(sequence_files, SingleTypeFiles) ?
            file[1:end-length(sequence_files.type)] * ".bam" :
            first(file)[1:end-length(sequence_files.type)-length(sequence_files.suffix1)] * ".bam"
        push!(outfiles, out_file)
        (isfile(out_file) && !overwrite_existing) && continue
        align_minimap(file, out_file, tmp_genome;
            preset=preset, match=match, mismatch=mismatch, gap_open1=gap_open1, gap_open2=gap_open2,
            gap_extend1=gap_extend1, gap_extend2=gap_extend2, minimizer_len=minimizer_len, threads=threads, all_secondary=all_secondary,
            minimap_bin=minimap_bin, sam_bin=sam_bin)
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
function align_mem(in_file1::String, in_file2::Union{String,Nothing}, out_file::String, genome_file::String;
    min_score=25, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false,
    min_seed_len=19, reseeding_factor=1.5, is_ont=false, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    params = ["-A", match, "-B", mismatch, "-O", gap_open, "-E", gap_extend, "-T", min_score, "-L", clipping_penalty, "-r", reseeding_factor, "-k", min_seed_len]
    isnothing(in_file2) || append!(params, ["-U", unpair_penalty])
    is_ont && append!(params, ["-x", "ont2d"])
    unpair_rescue && push!(params, "-P")
    fileparams = isnothing(in_file2) ? [genome_file, in_file1] : [genome_file, in_file1, in_file2]
    stats_file = out_file * ".log"
    run(pipeline(
        `$bwa_bin mem -v 1 -t 6 $params $fileparams`,
        stdout = pipeline(
            `$sam_bin view -u`,
            stdout = pipeline(
                `$sam_bin sort -o $out_file`))))
    run(pipeline(
        `$sam_bin index $out_file`,
        `$sam_bin stats $out_file`,
        stats_file))
end
align_mem(in_file::String, out_file::String, genome_file::String;
    min_score=25, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, min_seed_len=19,
    reseeding_factor=1.5, is_ont=false, bwa_bin="bwa-mem2", sam_bin="samtools") =
    align_mem(in_file, nothing, out_file::String, genome_file::String;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, clipping_penalty=clipping_penalty, min_seed_len=min_seed_len,
        reseeding_factor=reseeding_factor, is_ont=is_ont, bwa_bin=bwa_bin, sam_bin=sam_bin)

"""
    Helper dispatch of align_mem, which is a wrapper for bwa-mem2. Runs align_mem on `read_files::T`
    where T is a FileCollection and aligns it against `genome::Genome`. This enables easy handling of
    whole folders containing sequence files and the manipulation of a genome using Genome.
"""
function align_mem(read_files::T, genome::Genome; min_score=25, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, reseeding_factor=1.5,
                    min_seed_len=19, unpair_rescue=false, is_ont=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:FileCollection}
    tmp_genome = tempname()
    write(tmp_genome, genome)
    outfiles = Vector{String}()
    for file in read_files
        out_file = isa(read_files, SingleTypeFiles) ? file[1:end-length(read_files.type)] * ".bam" : first(file)[1:end-length(read_files.type)-length(read_files.suffix1)] * ".bam"
        push!(outfiles, out_file)
        (isfile(out_file) && !overwrite_existing) && continue
        isa(read_files, SingleTypeFiles) ?
        align_mem(file, out_file, tmp_genome;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open,
                gap_extend=gap_extend, clipping_penalty=clipping_penalty,  min_seed_len=min_seed_len,
                reseeding_factor=reseeding_factor, is_ont=is_ont, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(first(file), last(file), out_file, tmp_genome;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
                clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
                min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, bwa_bin=bwa_bin, sam_bin=sam_bin)
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
function align_mem(reads::Sequences, genomes::Vector{Genome}, out_file::String; min_score=25, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5,
                unpair_penalty=9, unpair_rescue=false, min_seed_len=19, reseeding_factor=1.5, is_ont=false, bwa_bin="bwa-mem2", sam_bin="samtools",
                overwrite_existing=false)

    (isfile(out_file) && !overwrite_existing) && return
    tmp_reads = tempname()
    tmp_reads2 = tempname()
    tmp_genome = tempname()
    ispaired = reads.seqnames[1:2:end] == reads.seqnames[2:2:end]
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
        align_mem(tmp_reads, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor,
            is_ont=is_ont, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(tmp_reads, tmp_reads2, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue,
            min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont, bwa_bin=bwa_bin, sam_bin=sam_bin)

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
    min_score=25, match=1, mismatch=4, gap_open=6, gap_extend=1, unpair_penalty=9, min_seed_len=19, reseeding_factor=1.5,
    unpair_rescue=false, is_ont=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:Sequences} =
        align_mem(reads, [genome], out_file;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, unpair_penalty=unpair_penalty,
            unpair_rescue=unpair_rescue, min_seed_len=min_seed_len, reseeding_factor=reseeding_factor, is_ont=is_ont,
            bwa_bin=bwa_bin, sam_bin=sam_bin, overwrite_existing=overwrite_existing)

function samevalueintervals(d::Vector{T}) where T
    index::Int = 1
    rindex::Int = 1
    lindex::Int = 0
    n_unique = length(unique(d))
    ranges = Vector{UnitRange{Int}}(undef, n_unique)
    while index < length(d)
        if d[index] != d[index+1]
            ranges[rindex] = lindex+1:index
            lindex = index
            rindex += 1
        end
        index += 1
    end
    length(d) == 1 && (ranges[1] = 1:1)
    n_unique>1 && (ranges[end] = last(ranges[end-1])+1:length(d))
    return ranges
end

function partsindex(ranges::Vector{UnitRange{Int}}, sorted_index::Vector{Int}, read_leftpos::Vector{Int}, reads::Vector{Symbol})
    pindex = Vector{Int}(undef, length(read_leftpos))
    for r in ranges
        for i in view(sorted_index, r)
            c = 0
            for ii in view(sorted_index, r)
                i === ii && continue
                if reads[i] === :read1
                    reads[ii] === :read2 && continue
                    (read_leftpos[ii] > read_leftpos[i]) && continue
                else
                    (reads[ii] === :read2) && (read_leftpos[ii] > read_leftpos[i]) && continue
                end
                (reads[i] === reads[ii]) && (read_leftpos[i] === read_leftpos[ii]) && (i>ii) && continue
                c += 1
            end
            pindex[first(r) + c] = i
        end
    end
    return pindex
end

"""
    Helper function to extract start and stop of the alignment on the read from
    a cigar string. Also compute the stop position on the reference and the complete length of the read.
"""
function readpositions(cigar::AbstractString)
    seqstart = 1
    seqstop = 0
    pending_seqstop = 0
    relrefstop = 0
    inseq = false
    seqlen = 0
    n = 0
    for c in cigar
        if isdigit(c)
            n = n * 10 + convert(Int, c - '0')
        else
            seqlen += n
            op = BioAlignments.Operation(c)
            if BioAlignments.isinsertop(op)
                inseq || (seqstart += n)
                pending_seqstop += n
            elseif BioAlignments.isdeleteop(op)
                relrefstop += n
                seqlen -= n
            elseif BioAlignments.ismatchop(op)
                inseq = true
                seqstop += n + pending_seqstop
                relrefstop +=n
                pending_seqstop = 0
            end
            n = 0
        end
    end
    return seqstart, seqstop, relrefstop, seqlen
end

"""
    Helper function to extract start and stop of the alignment on the read from
    a XAM.BAM record. Also compute the stop position on the reference and the
    complete length of the read.
"""
function readpositions(record::BAM.Record)
    offset, nops = BAM.cigar_position(record)
    seqstart = 1
    seqstop = 0
    pending_seqstop = 0
    relrefstop = 0
    inseq = false
    seqlen = 0
    for i in offset:4:offset + (nops - 1) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        op = BioAlignments.Operation(x & 0x0F)
        n = x >> 4
        seqlen += n
        if BioAlignments.isinsertop(op)
            inseq || (seqstart += n)
            pending_seqstop += n
        elseif BioAlignments.isdeleteop(op)
            relrefstop += n
            seqlen -= n
        elseif BioAlignments.ismatchop(op)
            inseq = true
            seqstop += n + pending_seqstop
            relrefstop +=n
            pending_seqstop = 0
        end
    end
    return seqstart, seqstop, relrefstop, seqlen
end

isprimary(record::BAM.Record)::Bool = BAM.flag(record) & 0x900 == 0
ispaired(record::BAM.Record)::Bool = BAM.flag(record) & 0x001 != 0
isread2(record::BAM.Record)::Bool =  BAM.flag(record) & 0x080 != 0
ispositivestrand(record::BAM.Record)::Bool = BAM.flag(record) & 0x010 == 0
mateispositivestrand(record::BAM.Record)::Bool = BAM.flag(record) & 0x020 == 0
paironsamestrand(record::BAM.Record, invert::Symbol)::Bool =
    (ispositivestrand(record) != (invert in (:both, :read1))) == (mateispositivestrand(record) != (invert in (:both, :read2)))

@inline function translateddata(data::SubArray{UInt8,1})
    for i in 1:length(data)
        (data[i] == 0x00) && (return data[1:i-1])
    end
end

@inline function xatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        (0x00 == data[i]) & (0x58 == data[i+1]) & (0x41 == data[i+2]) &&
        (return translateddata(@view(data[i+4:end])))
    end
    return nothing
end

function hasxatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        (0x00 == data[i]) & (0x58 == data[i+1]) & (0x41 == data[i+2]) && (return true)
    end
    return false
end

function hasxatag(record::BAM.Record)
    return hasxatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
end

function xatag(record::BAM.Record)
    xa = xatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
    return isnothing(xa) ? nothing : String(xa)
end

function astag(record::BAM.Record)
    data = @view(record.data[BAM.auxdata_position(record):BAM.data_size(record)])
    for i in 1:length(data)-4
        (0x00 == data[i]) & (0x41 == data[i+1]) & (0x53 == data[i+2]) && (return data[i+4])
    end
end

function nmtag(record::BAM.Record)
    record.data[BAM.auxdata_position(record)] == UInt8('N') && record.data[BAM.auxdata_position(record)+1] == UInt8('M') || throw("auxdata does not start with NM tag.")
    t = record.data[BAM.auxdata_position(record)+2] == UInt8('C') ? UInt8 :
        record.data[BAM.auxdata_position(record)+2] == UInt8('S') ? UInt16 :
        record.data[BAM.auxdata_position(record)+2] == UInt8('I') ? UInt32 :
        throw("NM tag type not supported: $(Char(record.data[BAM.auxdata_position(record)+2]))")
    return unsafe_load(Ptr{t}(pointer(record.data, BAM.auxdata_position(record)+3)))
end

struct Alignments{T<:Union{String, UInt}}
    tempnames::Vector{T}
    leftpos::Vector{Int}
    rightpos::Vector{Int}
    read_leftpos::Vector{Int}
    read_rightpos::Vector{Int}
    reads::Vector{Symbol}
    nms::Vector{UInt32}
    refnames::Vector{String}
    strands::Vector{Strand}
    annames::Vector{String}
    antypes::Vector{String}
    anols::Vector{UInt8}
    anleftrel::Vector{UInt8}
    anrightrel::Vector{UInt8}
    ranges::Vector{UnitRange{Int}}
end

function Alignments(bam_file::String; include_secondary_alignments=true, include_alternative_alignments=false, is_reverse_complement=false, hash_id=true)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    ns = Vector{hash_id ? UInt : String}(undef, 10000)
    ls = Vector{Int}(undef, 10000)
    rs = Vector{Int}(undef, 10000)
    is = Vector{String}(undef, 10000)
    ss = Vector{Strand}(undef, 10000)
    nms = Vector{UInt32}(undef, 10000)
    rls = Vector{Int}(undef, 10000)
    rrs = Vector{Int}(undef, 10000)
    rds = Vector{Symbol}(undef, 10000)
    index::Int = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        !isprimary(record) && !include_secondary_alignments && continue
        hasxatag(record) && !include_alternative_alignments && continue
        current_read = (isread2(record) != is_reverse_complement) ? :read2 : :read1
        index += 1
        if index > length(ns)
            for z in (ns, ls, rs, is, ss, rls, rrs, rds, nms)
                resize!(z, length(z)+10000)
            end
        end
        n::(hash_id ? UInt : String) = hash_id ? hash(@view(record.data[1:max(BAM.seqname_length(record) - 1, 0)])) : BAM.tempname(record)
        (l,r) = (BAM.leftposition(record), BAM.rightposition(record))
        (ref,s) = (BAM.refname(record), BAM.ispositivestrand(record) != (current_read === :read2) ? STRAND_POS : STRAND_NEG)
        nm = nmtag(record)
        if hasxatag(record)
            xa = xatag(record)
            xastrings = string.(split(xa, ";")[1:end-1])
            for xastring in enumerate(xastrings)
                ap = AlignedPart(xastring)
                (ns[index], ls[index], rs[index], is[index], ss[index], rls[index], rrs[index], rds[index], nms[index]) =
                (name(ap), leftposition(ap), rightposition(ap), refname(ap), strand(ap), first(readrange(ap)), last(readrange(ap)), current_read, nms(ap))
            end
        end
        readstart, readstop, _, readlen = readpositions(record)
        (rl, rr) = (BAM.ispositivestrand(record) != (current_read === :read2)) ? (readstart,readstop) : (readlen-readstop+1,readlen-readstart+1)
        (ns[index], ls[index], rs[index], is[index], ss[index], rls[index], rrs[index], rds[index], nms[index]) = (n, l, r, ref, s, rl, rr, current_read, nm)
    end
    for z in (ns, ls, rs, is, ss, rls, rrs, rds, nms)
        resize!(z, index)
    end
    close(reader)
    nindex = sortperm(ns)
    ns = ns[nindex]
    ranges = samevalueintervals(ns)
    pindex = partsindex(ranges, nindex, rls, rds)
    return Alignments(ns, ls[pindex], rs[pindex], rls[pindex], rrs[pindex], rds[pindex], nms[pindex], is[pindex], ss[pindex],
                        Vector{String}(undef,length(ns)), Vector{String}(undef,length(ns)), zeros(UInt8,length(ns)),
                        ones(UInt8,length(ns)) .* 0xff, ones(UInt8,length(ns)) .* 0xff, ranges)
end

Base.length(alns::Alignments) = length(alns.tempnames)
readscount(alns::Alignments) = length(alns.ranges)

function Base.iterate(alns::Alignments)
    return isempty(alns.ranges) ? nothing : (AlignedRead(alns.ranges[1], alns), 2)
end
function Base.iterate(alns::Alignments, state::Int)
    return state > length(alns.ranges) ? nothing : (AlignedRead(alns.ranges[state], alns), state+1)
end

function Base.filter!(seqs::Sequences{T}, alns::Alignments{T}) where {T<:Union{String, UInt}}
    filter!(seqs, Set(alns.tempnames))
end

struct AlignedPart
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    nms::UInt32
    read::Symbol
end

AlignedPart(alns::Alignments, i::Int) =
AlignedPart(
    Interval(alns.refnames[i], alns.leftpos[i], alns.rightpos[i], alns.strands[i],
        AlignmentAnnotation(
            isassigned(alns.antypes, i) ? alns.antypes[i] : "",
            isassigned(alns.annames, i) ? alns.annames[i] : "",
            alns.anols[i]
        )
    ),
    alns.read_leftpos[i]:alns.read_rightpos[i],
    alns.nms[i],
    alns.reads[i]
)

AlignedPart(alnpart::AlignedPart; new_name::Union{Nothing, String}=nothing, new_type::Union{Nothing, String}=nothing) =
AlignedPart(
    Interval(refname(alnpart), leftposition(alnpart), rightposition(alnpart), strand(alnpart),
        AlignmentAnnotation(
            isnothing(new_type) ? type(alnpart) : new_type,
            isnothing(new_name) ? name(alnpart) : new_name,
            overlap(alnpart)
        )
    ),
    alnpart.seq,
    alnpart.nms,
    alnpart.read
)

Base.getindex(alns::Alignments, i::Int) = AlignedPart(alns, i)
Base.getindex(alns::Alignments, r::UnitRange{Int}) = [alns[i] for i::Int in r]

"""
    Constructor for the AlignedPart struct. Builds AlignedPart from a XA string, which is created by bwa-mem2
    for alternative mappings. Inverts strand, if `check_invert::Bool` is true.
"""
function AlignedPart(xapart::String; read=:read1)
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = ((refstart > 0) != read===:read2) ? STRAND_POS : STRAND_NEG
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = readpositions(cigar)
    seq_interval = ((refstart > 0)  != (current_read === :read2)) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    return AlignedPart(ref_interval, seq_interval, parse(UInt32, nm), read)
end

"""

    Base.summarize(part::AlignedPart)::IO

Generates string with information on the AlignedPart.
"""
function summarize(part::AlignedPart)
    s = "[$(first(part.seq)), $(last(part.seq))] on $(part.read) - "
    s *= "[$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) ($(part.ref.strand)) "
    s *= "with edit distance $(part.nms) - "
    s *= isempty(annotation(part)) ? "not annotated." : "$(type(part)):$(name(part)) ($(overlap(part))% in annotation)"
    return s
end
"""
    Base.show(part::AlignedPart)::IO

Function for structured printing of the content of AlignedPart
"""
function Base.show(part::AlignedPart)
    println(summarize(part))
end

"""
    read(aln::AlignedPart)::Symbol

Returns the read from which the alignment comes. If is_reverse_complement is true, the reads will be
inverted (:read1 will be :read2 and vice versa)....
"""
sameread(aln1::AlignedPart, aln2::AlignedPart) = aln1.read === aln2.read

"""
    annotation(aln::AlignedPart)::AlignmentAnnotation

Returns the AlignmentAnnotation struct that is used as the metadata of the Interval{AlignmentAnnotation} and contains
information on the reference sequence part of the alignment.
"""
annotation(aln::AlignedPart) = aln.ref.metadata

"""
    name(aln::AlignedPart)::String

Returns the annotation's name as found in the .gff annotation file by the name_key argument of the Features constructor.
Availble only after using annotate! on Alignments.
"""
name(aln::AlignedPart) = name(annotation(aln))

"""
    type(aln::AlignedPart)::String

Returns the annotation's type as found in the .gff annotation file as the feature type.
Availble only after using annotate! on Alignments.
"""
type(aln::AlignedPart) = type(annotation(aln))

"""
    overlap(aln::AlignedPart)::UInt8

Returns the percentage of the sequence on the read contained within the annotation assigned to it.
Availble only after using annotate! on Alignments.
"""
overlap(aln::AlignedPart) = overlap(annotation(aln))

"""
    editdistance(aln::AlignedPart)::UInt8

Returns the number of edit operations between the aligned sequence and the reference (NM tag in SAM/BAM).
"""
editdistance(aln::AlignedPart) = aln.nms

"""
    refinterval(aln::AlignedPart)::Interval{AlignmentAnnotation}

Returns the Interval{AlignmentAnnotation} that contains information on the reference sequence part of the alignment.
Availble only after using annotate! on Alignments.
"""
refinterval(aln::AlignedPart) = aln.ref

"""
    refname(aln::AlignedPart)::String

Returns the id of the reference sequence the annotation comes from as found in the .gff annotation file as the chromosome name.
Availble only after using annotate! on Alignments.
"""
refname(aln::AlignedPart) = aln.ref.seqname

"""
    refrange(aln::AlignedPart)::UnitRange

Returns the Interval of the alignment on the reference sequence as a UnitRange.
"""
refrange(aln::AlignedPart) = leftposition(aln):rightposition(aln)

"""
    refsequence(aln::AlignedPart, genome::Genome)::LongDNASeq

Returns the part of the reference sequence that the aln::AlignedPart belongs to.
"""
refsequence(aln::AlignedPart, genome::Genome)::LongDNASeq = genome[refname(aln)][refrange(aln)]

"""
    leftposition(aln::AlignedPart)::Int

Returns the leftmost position of the alignment on the reference sequence.
Availble only after using annotate! on Alignments.
"""
BioGenerics.leftposition(aln::AlignedPart) = aln.ref.first

"""
    rightposition(aln::AlignedPart)::Int

Returns the rightmost position of the alignment on the reference sequence.
Availble only after using annotate! on Alignments.
"""
BioGenerics.rightposition(aln::AlignedPart) = aln.ref.last

"""
    strand(aln::AlignedPart)::Strand

Returns the strand of the annotation.
Availble only after using annotate! on Alignments.
"""
GenomicFeatures.strand(aln::AlignedPart) = aln.ref.strand

"""
    strand(aln::AlignedPart)::Strand

Returns the strand of the annotation.
Availble only after using annotate! on Alignments.
"""
ispositivestrand(aln::AlignedPart) = aln.ref.strand === STRAND_POS

"""
    readrange(aln::AlignedPart)::UnitRange

Returns the Interval of the alignment on the read sequence as a UnitRange.
"""
readrange(aln::AlignedPart) = aln.seq

"""
    nms(aln::AlignedPart)::UInt32

Returns the edit distance of the AlignedPart.
"""
nms(aln::AlignedPart) = aln.nms

function readsequence(aln::AlignedPart, seqs::Sequences)::LongDNASeq
    r = searchsorted(seqs.seqnames, name(aln))
    if length(r) === 2
        s = aln.read === :read1 ? seqs.seq[seqs.ranges[first(r)]] : seqs.seq[seqs.ranges[last(r)]]
        return s[readrange(aln)]
    elseif length(r) === 1
        return seqs.seq[seqs.ranges[first(r)]][readrange(aln)]
    else
        throw(KeyError)
    end
end

"""
    leftposition(rang::UnitRange)::Int

convenience function for consistency between Interval and UnitRange.
"""
BioGenerics.leftposition(rang::UnitRange) = first(rang)

"""
    rightposition(rang::UnitRange)::Int

convenience function for consistency between Interval and UnitRange.
"""
BioGenerics.rightposition(rang::UnitRange) = last(rang)

"""
    hasannotation(rang::UnitRange)::Bool

Checks if an annotation was assigned to the AlignedPart.
"""
hasannotation(aln::AlignedPart) = !isempty(annotation(aln))

"""
    sameannotation(aln1::AlignedPart, aln2::AlignedPart)::Bool

Checks if two AlignedParts share the same name and type in their annotation.
"""
sameannotation(aln1::AlignedPart, aln2::AlignedPart) = samename(aln1, aln2) && sametype(aln1, aln2)

"""
    samename(aln1::AlignedPart, aln2::AlignedPart)::Bool

Checks if two AlignedParts share the same name and type in their annotation.
"""
samename(aln1::AlignedPart, aln2::AlignedPart) = name(aln1) == name(aln2)

"""
    sametype(aln1::AlignedPart, aln2::AlignedPart)::Bool

Checks if two AlignedParts share the same name and type in their annotation.
"""
sametype(aln1::AlignedPart, aln2::AlignedPart) = type(aln1) == type(aln2)

"""
    distanceonread(aln1::AlignedPart, aln2::AlignedPart)::Int

Computes the distance between two alignments on the same read. Returns negative numbers
if the two parts overlap.
"""
distanceonread(aln1::AlignedPart, aln2::AlignedPart) = max(first(aln1.seq), first(aln2.seq)) - min(last(aln1.seq), last(aln2.seq))

isfirstread(part::AlignedPart) = part.read === :read1

struct AlignedRead
    range::UnitRange{Int}
    alns::Alignments
end

function Base.iterate(alnread::AlignedRead)
    return isnothing(alnread.range) ? nothing : (AlignedPart(alnread.alns, first(alnread.range)), first(alnread.range)+1)
end
function Base.iterate(alnread::AlignedRead, state::Int)
    return state > last(alnread.range) ? nothing : (AlignedPart(alnread.alns, state), state+1)
end
Base.getindex(alnread::AlignedRead, i::Int) = alnread.alns[alnread.range[i]]
Base.getindex(alnread::AlignedRead, r::UnitRange{Int}) = AlignedRead(alnread.range[r], alnread.alns)
Base.getindex(alnread::AlignedRead, b::Vector{Bool}) = [alnread.alns[index] for (i::Int,index::Int) in enumerate(alnread.range) if b[i]]
Base.isempty(alnread::AlignedRead) = isempty(alnread.range)
Base.length(alnread::AlignedRead) = length(alnread.range)

readid(alnread::AlignedRead) = alnread.alns.tempnames[first(alnread.range)]
parts(alnread::AlignedRead) = [AlignedPart(alnread.alns, i) for i in alnread.range]
typein(alnread::AlignedRead, types::Vector{String}) = any(alnread.alns.antypes[i] in types for i::Int in alnread.range if isassigned(alnread.alns.antypes, i))
hastype(alnread::AlignedRead, t::String) = any(alnread.alns.antypes[i] === t for i::Int in alnread.range if isassigned(alnread.alns.antypes, i))
namein(alnread::AlignedRead, names::Vector{String}) = any(alnread.alns.annames[i] in names for i::Int in alnread.range if isassigned(alnread.alns.annames, i))
hasname(alnread::AlignedRead, n::String) = any(alnread.alns.annames[i] === n for i::Int in alnread.range if isassigned(alnread.alns.annames, i))
hasannotation(alnread::AlignedRead) = any(isassigned(alnread.alns.annames, i) && isassigned(alnread.alns.antypes, i) for i in alnread.range)
annotatedcount(alnread::AlignedRead) = sum(isassigned(alnread.alns.annames, i) && isassigned(alnread.alns.antypes, i) for i in alnread.range)
annotationcount(alnread::AlignedRead) = length(Set(name(part) for part in alnread))
isfullyannotated(alnread::AlignedRead) = all(isassigned(alnread.alns.annames, i) && isassigned(alnread.alns.antypes, i) for i in alnread.range)
function BioGenerics.leftposition(alnread::AlignedRead)
    check_refname = alnread.alns.refnames[first(alnread.range)]
    all(v .== check_refname for v in view(alnread.alns.refnames, alnread.range)) || throw(AssertionError("AlignmentParts are not on the same reference sequence."))
    minimum(view(alnread.alns.leftpos, alnread.range))
end
function BioGenerics.rightposition(alnread::AlignedRead)
    check_refname = alnread.alns.refnames[first(alnread.range)]
    all(v .== check_refname for v in view(alnread.alns.refnames, alnread.range)) || throw(AssertionError("AlignmentParts are not on the same reference sequence."))
    maximum(view(alnread.alns.rightpos, alnread.range))
end

function GenomicFeatures.strand(alnread::AlignedRead)
    length(alnread) > 0 || (return STRAND_NA)
    #println(view(alnread.alns.strands, alnread.range), "\n",alnread.alns.strands[first(alnread.range)])
    check_strand = alnread.alns.strands[first(alnread.range)]
    return all(s === check_strand for s in view(alnread.alns.strands, alnread.range)) ? check_strand : STRAND_BOTH
end

function ispositivestrand(alnread::AlignedRead)
    s = strand(alnread)
    s === STRAND_NA && throw(AssertionError("Empty Alignment does not have a strand."))
    s === STRAND_BOTH && throw(AssertionError("Alignments are not on the same strand."))
    return s === STRAND_POS
end

summarize(alnread::AlignedRead) = "Alignment with $(length(alnread)) part(s):\n   " * join([summarize(part) for part in alnread], "\n   ")
function Base.show(alnread::AlignedRead)
    println(summarize(alnread))
end

refname(i::Interval) = i.seqname
Base.length(i::Interval{T}) where {T<:AnnotationStyle} = rightposition(i) - leftposition(i) + 1

"""
    overlapdistance(i1::Interval, i2::Interval)::Float64

Returns the negative distance between two AlignedParts on the reference sequence or the positive overlap if it exists.
Returns -Inf if the alignments do not share the same reference id or lie on different strands.
"""
function overlapdistance(i1::Interval{T}, i2::Interval{I})::Float64 where {T,I}
    (!(strand(i1) === strand(i2)) || !(refname(i1) === refname(i2))) && return -Inf
    return min(
        min(rightposition(i1) - leftposition(i1) + 1, rightposition(i2) - leftposition(i2) + 1),
        min(rightposition(i1) - leftposition(i2) + 1, rightposition(i2) - leftposition(i1) + 1)
    )
end

"""
    distance(i1::Interval, i2::Interval)::Float64

Returns the distance between two AlignedParts on the reference sequence. Returns Inf if the alignments do not share
the same reference id or lie on different strands
"""
distance(i1::Interval, i2::Interval)::Float64 = -min(-0.0, overlapdistance(i1,i2))

function distance(l1::Int, r1::Int, l2::Int, r2::Int)::Int
    l2>r1 && (return l2-r1)
    l1>r2 && (return l1-r2)
    return 0
end

function countchimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true)
    length(alnread) > 1 || (return 0)
    c = 0
    for (i1, i2) in combinations(alnread.range, 2)
        (check_annotation && isassigned(alnread.alns.annames, i1) &&  isassigned(alnread.alns.annames, i2) &&
            (alnread.alns.annames[i1] === alnread.alns.annames[i2])) && continue
        ((alnread.alns.refnames[i1] === alnread.alns.refnames[i2]) && (alnread.alns.strands[i1] === alnread.alns.strands[i2]) &&
            distance(alnread.alns.leftpos[i1], alnread.alns.rightpos[i1], alnread.alns.leftpos[i2], alnread.alns.rightpos[i2]) < min_distance) && continue
        c += 1
    end
    return c
end

function ischimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true)
    length(alnread) > 1 || (return false)
    for (i1, i2) in combinations(alnread.range, 2)
        (check_annotation && isassigned(alnread.alns.annames, i1) &&  isassigned(alnread.alns.annames, i2) &&
            (alnread.alns.annames[i1] === alnread.alns.annames[i2])) && continue
        ((alnread.alns.refnames[i1] === alnread.alns.refnames[i2]) && (alnread.alns.strands[i1] === alnread.alns.strands[i2]) &&
            distance(alnread.alns.leftpos[i1], alnread.alns.rightpos[i1], alnread.alns.leftpos[i2], alnread.alns.rightpos[i2]) < min_distance) && continue
        return true
    end
    return false
end

function ischimeric(part1::AlignedPart, part2::AlignedPart; min_distance=1000, check_annotation=true)
    check_annotation && hasannotation(part1) && hasannotation(part2) && (name(part1) == name(part2)) && (return false)
    return distance(refinterval(part1), refinterval(part2)) > min_distance
end

function ismulti(alnread::AlignedRead; method=:distance)
    return method === :annotation ? (annotationcount(alnread) >= 3) :
        method === :distance ? (countchimeric(alnread; check_annotation=false) >= 3) :
        method === :both ? (countchimeric(alnread; check_annotation=true) >= 3) :
        throw(AssertionError("method must be :distance, :annotation or :both"))
end

function Base.empty!(alns::Alignments)
    empty!(alns.tempnames)
    empty!(alns.leftpos)
    empty!(alns.rightpos)
    empty!(alns.refnames)
    empty!(alns.strands)
    empty!(alns.nms)
    empty!(alns.annames)
    empty!(alns.antypes)
    empty!(alns.anols)
    empty!(alns.ranges)
    empty!(alns.read_leftpos)
    empty!(alns.read_rightpos)
    empty!(alns.anleftrel)
    empty!(alns.anrightrel)
    empty!(alns.reads)
end

function occurences(test_sequence::LongDNASeq, bam_file::String, similarity_cut::Float64; score_model=nothing, ignore_mapped=true)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    c = 0
    while !eof(reader)
        read!(reader, record)
        (ignore_mapped && BAM.ismapped(record)) && continue
        c += similarity(test_sequence, BAM.sequence(record), score_model=score_model) > similarity_cut
    end
    return c
end

function annotate!(features::Features, files::SingleTypeFiles; only_unique_alignments=true, is_reverse_complement=false, count_key="Count", abs_lfc_cut=1, pvalue_cut=0.05)
    if files.type == ".bam"
        for (i,bam_file) in enumerate(files)
            reader = BAM.Reader(open(bam_file); index=bam_file*".bai")
            for feature in features
                c = 0
                for record in eachoverlap(reader, feature)
                    BAM.ismapped(record) || continue
                    hasxatag(record) && only_unique_alignments && continue
                    invert = is_reverse_complement != isread2(record)
                    same_strand = ispositivestrand(record) == (strand(feature) === STRAND_POS)
                    (same_strand == invert) && continue
                    c += 1
                end
                params(feature)["$count_key$i"] = "$c"
            end
            close(reader)
        end

    elseif files.type == ".csv"

        for file in files
            expname = basename(file)[1:end-4]
            table = CSV.read(file, DataFrame; missingstring="NA")
            replace!(table.padj, missing => 2.0)
            replace!(table.log2FoldChange, missing => 0.0)
            filter!(row->((row.padj <= pvalue_cut) && (abs(row.log2FoldChange) >= abs_lfc_cut)), table)
            for feature in features
                if name(feature) in table.Gene
                    row = table[findfirst(table.Gene .== name(feature)), :]
                    params(feature)["$(expname)_lfc"] = "$(row.log2FoldChange)"
                    params(feature)["$(expname)_fdr"] = "$(row.padj)"
                end
            end
        end
    end
end

annotate_filter(a::Interval{Annotation}, b::Interval{Nothing}) = strand(a) === strand(b)
function annotate!(alns::Alignments, features::Features{Annotation}; prioritize_type=nothing, overwrite_type=nothing)
    myiterators = Dict(refn=>GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(annotate_filter), Annotation}(
                                    annotate_filter,
                                    GenomicFeatures.ICTreeIntersection{Annotation}(),
                                    features.list.trees[refn],
                                    Interval("", 1, 1, STRAND_NA, Annotation()))
                            for refn in refnames(features))

    for i::Int in 1:length(alns)
        alns.refnames[i] in keys(myiterators) || continue
        myiterator = myiterators[alns.refnames[i]]
        myiterator.query = Interval(alns.refnames[i], alns.leftpos[i], alns.rightpos[i], alns.strands[i])
        for feature_interval in myiterator
            olp = round(UInt8, (min(min(feature_interval.last - feature_interval.first + 1, alns.rightpos[i] - alns.leftpos[i] + 1),
                                    min(feature_interval.last - alns.leftpos[i] + 1, alns.rightpos[i] - feature_interval.first + 1)) /
                                (alns.rightpos[i] - alns.leftpos[i] + 1)) * 100)

            priority = !isnothing(prioritize_type) && (type(feature_interval) === prioritize_type) && (olp > 80)
            overwrite = !isnothing(overwrite_type) && isassigned(alns.antypes, i) && (alns.antypes[i] === overwrite_type)

            if  priority || overwrite || alns.anols[i]<olp
                alns.antypes[i] = type(feature_interval)
                alns.annames[i] = name(feature_interval)
                feature_length = rightposition(feature_interval) - leftposition(feature_interval) + 1
                alns.anleftrel[i] = alns.leftpos[i] > leftposition(feature_interval) ?
                    round(UInt8, (alns.leftpos[i] - leftposition(feature_interval) + 1) / feature_length * 100) : 0x00
                alns.anrightrel[i] = alns.rightpos[i] < rightposition(feature_interval) ?
                    round(UInt8, (alns.rightpos[i] - leftposition(feature_interval) + 1) / feature_length * 100) : 0x65
                strand(feature_interval) === STRAND_NEG && ((alns.anleftrel[i], alns.anrightrel[i]) = (0x65 - alns.anrightrel[i], 0x65 - alns.anleftrel[i]))
                alns.anols[i] = olp
                priority && break
            end
        end
    end
end

function annotate!(features::Features, feature_alignments::Alignments; key_gen=typenamekey)
    for feature in features
        key = key_gen(feature)
        push!(params(feature), "Count"=>"0")
        push!(params(feature), "Key"=>key)
        if key in keys(feature_alignments)
            refs = Set(refname(aln) for aln in feature_alignments[key])
            merge!(params(feature), Dict(ref=>join(("$(readrange(aln))" for aln in feature_alignments[key] if refname(aln)==ref), ",") for ref in refs))
            push!(params(feature), "Count"=>"$(length(refs))")
        end
    end
end

Base.getindex(genome::Genome, ap::AlignedPart) = refsequence(ap, genome)
Base.getindex(sequence::LongDNASeq, ap::AlignedPart) = sequence[readrange(ap)]

struct GenomeComparison
    alns::Vector{PairwiseAlignment}
    fromto::Vector{Tuple{Interval,Interval}}
end

function GenomeComparison(refgenome_file::String, compgenome_file::String, bam_file::String)
    pairwise_alignments = PairwiseAlignment[]
    fromto = Tuple{Interval,Interval}[]
    refgenome = Genome(refgenome_file)
    compgenome = Genome(compgenome_file)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    while !eof(reader)
        read!(reader, record)
        seqstart, seqstop, _, _ = readpositions(record)
        refseq = refgenome[BAM.refname(record)]#[BAM.leftposition(record):BAM.rightposition(record)]
        compseq = compgenome[BAM.seqname(record)]#[seqstart:seqstop]
        BAM.ispositivestrand(record) || reverse_complement!(compseq)
        aln = BAM.alignment(record)
        #filter!(x -> x.op === OP_HARD_CLIP || x.op === OP_SOFT_CLIP, aln.anchors)
        alnseq = AlignedSequence(refseq, aln)
        pwa = PairwiseAlignment(alnseq, compseq)
        push!(pairwise_alignments, pwa)
        push!(fromto, (Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), STRAND_NA),
                        Interval(BAM.tempname(record), seqstart, seqstop, BAM.ispositivestrand(record) ? STRAND_POS : STRAND_NEG)))
    end
    return GenomeComparison(pairwise_alignments, fromto)
end

Base.length(genomecomp::GenomeComparison) = length(genomecomp.alns)
Base.iterate(genomecomp::GenomeComparison) = ((genomecomp.alns[1], genomecomp.fromto[1]...), 2)
Base.iterate(genomecomp::GenomeComparison, state::Int) = state > length(genomecomp) ? nothing :
                                                            ((genomecomp.alns[state], genomecomp.fromto[state]...), state+1)

cutfill(str::String, len::Int) = length(str) > len ? str[1:len] : str * repeat(" ", len - length(str))
function summarize(genomecomp::GenomeComparison)
    s = "refname\t\tfrom\tto\tstrand\tcompname\tfrom\tto\tmatches\tsnps\tinserts\tdels\n"
    for (i, (pwa, from, to)) in enumerate(genomecomp)
        lref, rref = cutfill(string(leftposition(from)), 8), cutfill(string(rightposition(from)), 8)
        lcomp, rcomp = cutfill(string(leftposition(to)), 8), cutfill(string(rightposition(to)), 8)
        refn, compn = cutfill(refname(from), 16), cutfill(refname(to), 16)
        matches = cutfill(string(count_matches(pwa)), 8)
        snp = cutfill(string(count_mismatches(pwa)), 8)
        inserts = cutfill(string(count_insertions(pwa)), 8)
        dels = cutfill(string(count_deletions(pwa)), 8)
        st = cutfill(strand(to) === STRAND_POS ? "+" : "-", 8)
        s *= "$refn$lref$rref$st$compn$lcomp$rcomp$matches$snp$inserts$dels\n"
    end
    return s
end

Base.show(genomecomp::GenomeComparison) = println(summarize(genomecomp))