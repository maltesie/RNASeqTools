#!/usr/bin/env julia
"""
Wrapper function for kraken2 taxonomic sequence classifier that assigns taxonomic labels to DNA sequences.

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

function align_mem(in_file1::String, in_file2::Union{String,Nothing}, out_file::String, genome_file::String;
    min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false, reseeding_factor=1.5,
    bwa_bin="bwa-mem2", sam_bin="samtools")

    tmp_bwa = tempname()
    tmp_view = tempname()

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    params = ["-A", match, "-B", mismatch, "-O", gap_open, "-E", gap_extend, "-T", min_score, "-L", clipping_penalty, "-r", reseeding_factor]
    isnothing(in_file2) || append!(params, ["-U", unpair_penalty])
    unpair_rescue && push!(params, "-P")
    fileparams = isnothing(in_file2) ? [genome_file, in_file1] : [genome_file, in_file1, in_file2]
    cmd = pipeline(`$bwa_bin mem -v 1 -t 6 $params $fileparams`, stdout=tmp_bwa)
    run(cmd)
    cmd = pipeline(`$sam_bin view -u $tmp_bwa`, stdout=tmp_view)
    run(cmd)
    cmd = pipeline(`$sam_bin sort $tmp_view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)
    stats_file = out_file * ".log"
    cmd = pipeline(`$sam_bin stats $out_file`, stdout=stats_file)
    run(cmd)

    rm(tmp_bwa)
    rm(tmp_view)
end
align_mem(in_file::String, out_file::String, genome_file::String;
    min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, reseeding_factor=1.5, bwa_bin="bwa-mem2", sam_bin="samtools") =
    align_mem(in_file, nothing, out_file::String, genome_file::String;
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, clipping_penalty=clipping_penalty, reseeding_factor=reseeding_factor,
        bwa_bin=bwa_bin, sam_bin=sam_bin)

function align_mem(read_files::T, genome::Genome; min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, reseeding_factor=1.5,
                unpair_rescue=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where {T<:FileCollection}
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
                gap_extend=gap_extend, clipping_penalty=clipping_penalty, reseeding_factor=reseeding_factor, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(first(file), last(file), out_file, tmp_genome;
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
                clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue, reseeding_factor=reseeding_factor, bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        isfile(tmp_genome * ending) && rm(tmp_genome * ending)
    end
    return SingleTypeFiles(outfiles)
end

function align_mem(reads::T, genomes::Vector{Genome}, out_file::String; min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5,
                unpair_penalty=9, unpair_rescue=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where T<:SequenceContainer
    (isfile(out_file) && !overwrite_existing) && return
    tmp_reads = tempname()
    tmp_reads2 = tempname()
    tmp_genome = tempname()
    write(tmp_reads, reads)
    for (i,genome) in enumerate(genomes)
        write(tmp_genome, genome)
        this_out_file = out_file
        length(genomes) > 1 && (this_out_file = joinpath(dirname(out_file), "$(i)_" * basename(out_file)))

        isa(reads, Sequences) ?
        align_mem(tmp_reads, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(tmp_reads, tmp_reads2, this_out_file, tmp_genome;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend,
            clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue, bwa_bin=bwa_bin, sam_bin=sam_bin)

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
    min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, unpair_penalty=9,
    unpair_rescue=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false) where T<:SequenceContainer =
        align_mem(reads, [genome], out_file;
            min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, unpair_penalty=unpair_penalty,
            unpair_rescue=unpair_rescue, bwa_bin=bwa_bin, sam_bin=sam_bin, overwrite_existing=overwrite_existing)

struct AlignedPart
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    nms::UInt8
    read::Symbol
end

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
        elseif BioAlignments.ismatchop(op)
            inseq = true
            seqstop += n + pending_seqstop
            relrefstop +=n
            pending_seqstop = 0
        end
    end
    return seqstart, seqstop, relrefstop, seqlen
end

function AlignedPart(xapart::Union{String, SubString{String}}; invert=false, read=:read1)
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = (refstart > 0) == !invert ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = readpositions(cigar)
    seq_interval = (refstart > 0) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    return AlignedPart(ref_interval, seq_interval, parse(Int8, nm), read)
end

annotation(aln::AlignedPart) = aln.ref.metadata
name(aln::AlignedPart) = annotation(aln).name
type(aln::AlignedPart) = annotation(aln).type
overlap(aln::AlignedPart) = annotation(aln).overlap
refinterval(aln::AlignedPart) = aln.ref
refname(aln::AlignedPart) = aln.ref.seqname
BioGenerics.leftposition(aln::AlignedPart) = aln.ref.first
BioGenerics.rightposition(aln::AlignedPart) = aln.ref.last
GenomicFeatures.strand(aln::AlignedPart) = aln.ref.strand
readinterval(aln::AlignedPart) = aln.seq
BioGenerics.leftposition(rang::UnitRange) = first(rang)
BioGenerics.rightposition(rang::UnitRange) = last(rang)
hasannotation(aln::AlignedPart) = !isempty(annotation(aln))
sameannotation(aln1::AlignedPart, aln2::AlignedPart) = samename(aln1, aln2) && sametype(aln1, aln2)
samename(aln1::AlignedPart, aln2::AlignedPart) = name(aln1) == name(aln2)
sametype(aln1::AlignedPart, aln2::AlignedPart) = type(aln1) == type(aln2)
#core(aln::AlignedPart) = aln.score
#isprimary(aln::AlignedPart) = aln.isprimary

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
    return record.data[BAM.auxdata_position(record)+3]
end

struct AlignedRead
    alns::Vector{AlignedPart}
end

function AlignedRead()
    return AlignedRead([])
end

alignments(alnread::AlignedRead) = alnread.alns
count(alnread::AlignedRead) = length(alnread.alns)
typein(alnread::AlignedRead, types::Vector{String}) = any(type(alnpart) in types for alnpart in alnread)
hastype(alnread::AlignedRead, t::String) = typein(alnread, [t])
namein(alnread::AlignedRead, names::Vector{String}) = any(name(alnpart) in names for alnpart in alnread)
hasname(alnread::AlignedRead, n::String) = namein(alnread, [n])
hasannotation(alnread::AlignedRead) = any(hasannotation(alnpart) for alnpart in alnread)
annotatedcount(alnread::AlignedRead) = sum(hasannotation(alnpart) for alnpart in alnread)
annotationcount(alnread::AlignedRead) = length(Set(name(part) for part in alnread))
isfullyannotated(alnread::AlignedRead) = all(hasannotation(alnpart) for alnpart in alnread)

function merge!(alnread1::AlignedRead, alnread2::AlignedRead)
    for part in alnread2
        push!(alnread1, part)
    end
end

function name(alnread::AlignedRead)
    names = Set(name(part) for part in alnread)
    @assert length(names) == 1
    for n in names
        return n
    end
end
function name(alnread1::AlignedRead, alnread2::AlignedRead)
    name1 = name(alnread1)
    name2 = name(alnread2)
    @assert name1 === name2
    return name1
end

Base.length(alnread::AlignedRead) = length(alnread.alns)
Base.isempty(alnread::AlignedRead) = isempty(alnread.alns)
Base.iterate(alnread::AlignedRead) = iterate(alnread.alns)
Base.iterate(alnread::AlignedRead, state::Int) = iterate(alnread.alns, state)
Base.getindex(alnread::AlignedRead, i::Int64) = alnread.alns[i]

function Base.show(alnread::AlignedRead)
    println("\nAlignment with $(length(alnread.alns)) part(s):")
    for part in alnread.alns
        print("   $(part.ref.strand): [$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) - ($(first(part.seq)), $(last(part.seq))) on $(part.read) with $(part.nms) $(part.nms == 1 ? "missmatch" : "missmatches") - ")
        isempty(annotation(part)) ? println("not annotated.") : println("annotation: $(type(part)):$(name(part)) [$(overlap(part))]")
    end
end

function hasannotation(alnread::AlignedRead, annotation_name::String)
    for part in alnread
        name(part) == annotation_name && (return true)
    end
    return false
end

refname(i::Interval) = i.seqname
Base.length(i::Interval{T}) where {T<:AnnotationStyle} = rightposition(i) - leftposition(i) + 1
function overlapdistance(i1::Interval, i2::Interval)::Float64
    strand(i1) != strand(i2) && return -Inf
    return min(min(length(i1), length(i2)), min(rightposition(i1) - leftposition(i2), rightposition(i2) - leftposition(i1)))
end
distance(i1::Interval, i2::Interval)::Float64 = -min(-0.0, overlapdistance(i1,i2))

function ischimeric(part1::AlignedPart, part2::AlignedPart; min_distance=1000, check_annotation=true)
    check_annotation && hasannotation(part1) && hasannotation(part2) && (name(part1) == name(part2)) && (return false)
    return distance(refinterval(part1), refinterval(part2)) > min_distance
end

function countchimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true)
    length(alnread) > 1 ?
    sum(ischimeric(part, otherpart; min_distance=min_distance, check_annotation=check_annotation) for (part, otherpart) in combinations(alnread.alns, 2)) :
    0
end

function ischimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true)
    countchimeric(alnread; min_distance=min_distance, check_annotation=check_annotation) > 0
end

function ismulti(alnread::AlignedRead; min_distance=1000, check_annotation=true)
    countchimeric(alnread; min_distance=min_distance, check_annotation=check_annotation) >= length(alnread)
end

struct Alignments{T} <: AlignmentContainer
    dict::Dict{T, AlignedRead}
end

Base.empty!(alignments::Alignments) = empty!(alignments.dict)
Base.getindex(alignments::Alignments, key::Union{String,UInt}) = alignments.dict[key]
Base.length(alignments::Alignments) = length(alignments.dict)
Base.keys(alignments::Alignments) = keys(alignments.dict)
Base.values(alignments::Alignments) = values(alignments.dict)
function Base.iterate(alignments::Alignments)
    dictiteration = iterate(alignments.dict)
    isnothing(dictiteration) && (return nothing)
    ((_, aln), state) = dictiteration
    return (aln, state)
end
function Base.iterate(alignments::Alignments, state::Int)
    dictiteration = iterate(alignments.dict, state)
    isnothing(dictiteration) && (return nothing)
    ((_, aln), state) = dictiteration
    return (aln, state)
end

function Alignments(bam_file::String; min_templength=nothing, only_unique=true, invert_strand=:none, reverse_order=false, hash_id=true)
    alignments= read_bam(bam_file; min_templength=min_templength, only_unique=only_unique, invert_strand=invert_strand, reverse_order=reverse_order, hash_id=hash_id)
    Alignments(alignments)
end

function Alignments(bam_file1::String, bam_file2::String; min_templength=nothing, only_unique=true, invert_strand=:none, reverse_order=false, hash_id=true)
    alignments = read_bam(bam_file1; min_templength=min_templength, only_unique=only_unique, invert_strand=invert_strand, reverse_order=reverse_order, hash_id=hash_id)
    read_bam!(alignments, bam_file2; min_templength=min_templength, only_unique=only_unique, invert_strand=invert_strand in (:read2, :both) ? :read1 : :none, reverse_order=reverse_order)
    Alignments(alignments)
end

function ispaired(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x001) != 0
end

function isproperpair(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x002) != 0
end

function isprimary(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x100) == 0
end

function isread1(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x040) != 0
end

function isread2(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x080) != 0
end

function ispositivestrand(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x010) == 0
end

function mateispositivestrand(record::BAM.Record)::Bool
    BAM.flag(record) & UInt(0x020) == 0
end

function paironsamestrand(record::BAM.Record, invert_strand::Symbol)::Bool
    (ispositivestrand(record) != (invert_strand in (:both, :read1))) == (mateispositivestrand(record) != (invert_strand in (:both, :read2)))
end

function is_bitstring_bam(file::String)
    reader = BAM.Reader(open(file))
    tempname = ""
    for record in reader
        tempname = BAM.tempname(record)
        break
    end
    close(reader)
    ((length(tempname) == 64) && all([c in ['0', '1'] for c in tempname])) && (return true)
    return false
end

function Base.push!(l::Vector{AlignedPart}, item::AlignedPart; reverse_order=false)
    for (i, part) in enumerate(l)
        if item.read === part.read
            first(item.seq) < first(part.seq) && (return insert!(l, i, item))
        elseif (item.read === :read1) != reverse_order
            return insert!(l, i, item)
        end
    end
    insert!(l, length(l)+1, item)
end

function read_bam!(reads::Dict{T, AlignedRead}, bam_file::String; min_templength=nothing, only_unique=true, invert_strand=:none, reverse_order=false) where T<:Union{UInt, String}
    @assert invert_strand in [:read1, :read2, :both, :none]
    hash_id = T <: UInt
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    is_bitstring = is_bitstring_bam(bam_file)
    check_templen = isnothing(min_templength) ? false : 0 < min_templength
    invert2 = invert_strand in (:read2, :both)
    invert1 = invert_strand in (:read1, :both)
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        is_unique = !hasxatag(record)
        !is_unique && only_unique && continue
        (check_templen && (0 < abs(BAM.templength(record)) < min_templength) && paironsamestrand(record, invert_strand) && isproperpair(record)) && continue
        id = hash_id ? (is_bitstring ? parse(UInt, BAM.tempname(record); base=2) : hash(@view(record.data[1:BAM.seqname_length(record)]))) : BAM.tempname(record)
        current_read = isread2(record) ? :read2 : :read1
        invert = (current_read === :read2 && invert2) || (current_read === :read1 && invert1)
        foundit = false
        nms = nmtag(record)
        if !is_unique && !only_unique
            xa = xatag(record)
            leftest = 0
            leftestpos = BAM.leftposition(record)
            xastrings = split(xa, ";")[1:end-1]
            for (i,xapart) in enumerate(xastrings)
                _, pos, _, nm = split(xapart, ",")
                pnm = parse(UInt8, nm)
                pnm > nms && continue
                intpos = abs(parse(Int, pos))
                intpos < leftestpos && (leftest = i; leftestpos=intpos; nms=pnm)
            end
            if leftest != 0
                alnpart = AlignedPart(xastrings[leftest]; invert=invert, read=current_read)
                foundit = true
            end
        end
        if !foundit
            readstart, readstop, _, readlen = readpositions(record)
            seq_interval = BAM.ispositivestrand(record) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
            st = BAM.ispositivestrand(record) == !invert ? Strand('+') : Strand('-')
            ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), st, AlignmentAnnotation())
            alnpart = AlignedPart(ref_interval, seq_interval, nms, current_read)
        end
        id in keys(reads) ? push!(reads[id].alns, alnpart; reverse_order=reverse_order) : push!(reads, id=>AlignedRead([alnpart]))
    end
    close(reader)
    return reads
end
function read_bam(bam_file::String; min_templength=nothing, only_unique=true, invert_strand=:none, reverse_order=false, hash_id=true)
    reads = Dict{hash_id ? UInt : String, AlignedRead}()
    read_bam!(reads, bam_file; min_templength=min_templength, only_unique=only_unique, invert_strand=invert_strand, reverse_order=reverse_order)
    return reads
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

function annotate!(features::Features, files::SingleTypeFiles; only_unique=true, invert_strand=:none, count_key="Count", abs_lfc_cut=1, pvalue_cut=0.05)
    if files.type == ".bam"

        @assert invert_strand in [:read1, :read2, :both, :none]
        for (i,bam_file) in enumerate(files)
            reader = BAM.Reader(open(bam_file); index=bam_file*".bai")
            for feature in features
                c = 0
                for record in eachoverlap(reader, feature)
                    BAM.ismapped(record) || continue
                    hasxatag(record) && only_unique && continue
                    invert = isread2(record) && ispaired(record) ? invert_strand in (:read2, :both) : invert_strand in (:read1, :both)
                    (ispositivestrand(record) == (strand(feature) === STRAND_POS) == invert) && continue
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

function annotate!(aln::AlignedPart, feature_interval::Interval{Annotation}; prioritize_type=nothing, overwrite_type=nothing)
    olp = round(UInt8, (overlapdistance(feature_interval, refinterval(aln)) / length(refinterval(aln))) * 100)
    foundpriority = (!isnothing(prioritize_type) && (type(feature_interval) == prioritize_type))
    overwrite = (!isnothing(overwrite_type) && (type(annotation(aln)) == overwrite_type))
    if foundpriority || overwrite || overlap(aln)<olp
        annotation(aln).type = feature_interval.metadata.type
        annotation(aln).name = feature_interval.metadata.name
        annotation(aln).overlap = olp
    end
    return foundpriority
end

function annotate!(alns::Alignments, features::Features; prioritize_type=nothing, overwrite_type=nothing, remove_empty=false)
    myiterators = Dict(refn=>GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter,
        GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refn], Interval("", 1, 1, STRAND_NA, Annotation("", "", Dict{String,String}()))) for refn in refnames(features))
    for alignedread in alns
        for aln in alignedread
            myiterator = myiterators[refname(aln)]
            myiterator.query = refinterval(aln)
            for feature_interval in myiterator
                annotate!(aln, feature_interval; prioritize_type=prioritize_type, overwrite_type=overwrite_type) && break
            end
        end
        remove_empty && deleteat!(alignedread.alns, findall(x->!hasannotation(x), alignedread.alns))
    end
end

function annotate!(features::Features, feature_alignments::Alignments{String}; key_gen=typenamekey)
    for feature in features
        key = key_gen(feature)
        push!(params(feature), "Count"=>"0")
        push!(params(feature), "Key"=>key)
        if key in keys(feature_alignments)
            refs = Set(refname(aln) for aln in feature_alignments[key])
            merge!(params(feature), Dict(ref=>join(("$(readinterval(aln))" for aln in feature_alignments[key] if refname(aln)==ref), ",") for ref in refs))
            push!(params(feature), "Count"=>"$(length(refs))")
        end
    end
end
