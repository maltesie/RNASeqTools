function align_mem(in_file1::String, in_file2::Union{String,Nothing}, out_file::String, genome_file::String; 
    min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, unpair_rescue=false, 
    bwa_bin="bwa-mem2", sam_bin="samtools")

    tmp_bwa = tempname()
    tmp_view = tempname()

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    params = ["-A", match, "-B", mismatch, "-O", gap_open, "-E", gap_extend, "-T", min_score, "-L", clipping_penalty]
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
    min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, bwa_bin="bwa-mem2", sam_bin="samtools") = 
    align_mem(in_file, nothing, out_file::String, genome_file::String; 
        min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, clipping_penalty=clipping_penalty, bwa_bin=bwa_bin, sam_bin=sam_bin)

function align_mem(read_files::T, genome::Genome; min_score=30, match=1, mismatch=4, gap_open=6, gap_extend=1, clipping_penalty=5, unpair_penalty=9, 
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
                gap_extend=gap_extend, clipping_penalty=clipping_penalty, bwa_bin=bwa_bin, sam_bin=sam_bin) :
        align_mem(first(file), last(file), out_file, tmp_genome; 
                min_score=min_score, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, 
                clipping_penalty=clipping_penalty, unpair_penalty=unpair_penalty, unpair_rescue=unpair_rescue, bwa_bin=bwa_bin, sam_bin=sam_bin)
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

function local_alignment(reference_sequence::LongDNASeq, query_sequence::LongDNASeq, scoremodel::AffineGapScoreModel)
    res = pairalign(LocalAlignment(), query_sequence, reference_sequence, scoremodel)
    return res
end

struct AlignedPart
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    isprimary::Bool
    score::UInt8
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

function AlignedPart(xapart::Union{String, SubString{String}}; invertstrand=false)
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = (refstart > 0) == !invertstrand ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = readpositions(cigar)
    seq_interval = (refstart > 0) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    return AlignedPart(ref_interval, seq_interval, false, 0)
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
score(aln::AlignedPart) = aln.score
isprimary(aln::AlignedPart) = aln.isprimary

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
merge!(alnread1::AlignedRead, alnread2::AlignedRead) = append!(alnread1.alns, alnread2.alns)
hasannotation(alnread::AlignedRead) = any(hasannotation(alnpart) for alnpart in alnread)
isfullyannotated(alnread::AlignedRead) = all(hasannotation(alnpart) for alnpart in alnread)

Base.length(alnread::AlignedRead) = length(alnread.alns)
Base.isempty(alnread::AlignedRead) = isempty(alnread.alns)
Base.iterate(alnread::AlignedRead) = iterate(alnread.alns)
Base.iterate(alnread::AlignedRead, state::Int) = iterate(alnread.alns, state)
Base.getindex(alnread::AlignedRead, i::Int64) = alnread.alns[i]

function Base.show(alnread::AlignedRead)
    println("Alignment with $(length(alnread.alns)) part(s):")
    for part in alnread.alns
        print("   $(part.ref.strand): [$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) - ($(first(part.seq)), $(last(part.seq))) on read - ")
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

function countconcordant(alnread1::AlignedRead, alnread2::AlignedRead; min_distance=100, check_annotation=true)
    c = 0
    for part in alnread1, otherpart in alnread2
        if check_annotation 
            if hasannotation(part) && hasannotation(otherpart) 
                if (name(part) == name(otherpart))
                    c+=1
                    continue
                end
            end
        end
        distance(refinterval(part), refinterval(otherpart)) < min_distance && (c+=1)
    end
    return c
end

function ischimeric(alnread::AlignedRead)
    return count(alnread) > 1 ? true : false
end

function ischimeric(alnread1::AlignedRead, alnread2::AlignedRead; min_distance=100, check_annotation=true, per_read=false)
    per_read && (ischimeric(alnread1) || ischimeric(alnread2)) && (return true)
    return (count(alnread1) + count(alnread2) - countconcordant(alnread1, alnread2; min_distance=min_distance, check_annotation=check_annotation)) >= 2
end

function istriplet(alnread1::AlignedRead, alnread2::AlignedRead; min_distance=100, check_annotation=true)
    2 < (count(alnread1) + count(alnread2)) < 5 || (return false)
    return (count(alnread1) + count(alnread2) - countconcordant(alnread1, alnread2; min_distance=min_distance, check_annotation=check_annotation)) == 3
end

struct Alignments{T} <: AlignmentContainer
    dict::Dict{T, AlignedRead}
end

Base.getindex(alignments::Alignments, key::Union{String,UInt}) = alignments.dict[key]
Base.length(alignments::Alignments) = length(alignments.dict)
Base.keys(alignments::Alignments) = keys(alignments.dict)
Base.values(alignments::Alignments) = values(alignments.dict)
function Base.iterate(alignments::Alignments) 
    dictiteration = iterate(alignments.dict)
    isnothing(dictiteration) && (return nothing)
    ((key, aln), state) = dictiteration
    return (aln, state)
end
function Base.iterate(alignments::Alignments, state::Int) 
    dictiteration = iterate(alignments.dict, state)
    isnothing(dictiteration) && (return nothing)
    ((key, aln), state) = dictiteration
    return (aln, state)
end

function Alignments(bam_file::String; min_templength=nothing, only_unique=true, rev_comp=:none, hash_id=true)
    alignments1, alignments2 = read_bam(bam_file; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp, hash_id=hash_id)
    @assert isempty(alignments2)
    Alignments(alignments1)
end

struct PairedAlignments{T} <: AlignmentContainer
    read1::Dict{T, AlignedRead}
    read2::Dict{T, AlignedRead}
    keys::Set{T}
end

Base.getindex(alignments::PairedAlignments, key::Union{String,UInt}) = (alignments.read1[key], alignments.read2[key])
Base.keys(alignments::PairedAlignments) = alignments.keys
Base.length(alignments::PairedAlignments) = length(alignments.keys)
function Base.iterate(alignments::PairedAlignments) 
    dictiteration = iterate(alignments.keys)
    isnothing(dictiteration) && (return nothing)
    key, state = dictiteration
    return ((alignments.read1[key], alignments.read2[key]), state)
end
function Base.iterate(alignments::PairedAlignments, state::Int) 
    dictiteration = iterate(alignments.keys, state)
    isnothing(dictiteration) && (return nothing)
    key, state = dictiteration
    return ((alignments.read1[key], alignments.read2[key]), state)
end

function PairedAlignments(bam_file1::String, bam_file2::String; min_templength=nothing, only_unique=true, rev_comp=:none, hash_id=true)
    alignments1, alignments_e1 = read_bam(bam_file1; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp, hash_id=hash_id)
    alignments2, alignments_e2 = read_bam(bam_file2; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp, hash_id=hash_id)
    @assert isempty(alignments_e1) && isempty(alignments_e2)
    PairedAlignments(alignments1, alignments2, intersect(Set(keys(alignments1)), Set(keys(alignments2))))
end

function PairedAlignments(pebam_file::String; min_templength=nothing, only_unique=true, rev_comp=:none, hash_id=true)
    alignments1, alignments2 = read_bam(pebam_file; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp, hash_id=hash_id)
    PairedAlignments(alignments1, alignments2, intersect(Set(keys(alignments1)), Set(keys(alignments2))))
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

function paironsamestrand(record::BAM.Record, invertstrand::Symbol)::Bool
    (ispositivestrand(record) != (invertstrand in (:both, :read1))) == (mateispositivestrand(record) != (invertstrand in (:both, :read2)))
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

function read_bam(bam_file::String; min_templength=nothing, only_unique=true, invertstrand=:none, hash_id=true)
    @assert invertstrand in [:read1, :read2, :both, :none]
    reads1 = Dict{hash_id ? UInt : String, AlignedRead}()
    reads2 = Dict{hash_id ? UInt : String, AlignedRead}()
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    is_bitstring = is_bitstring_bam(bam_file)
    check_templen = isnothing(min_templength) ? false : 0 < min_templength
    c = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        (check_templen && (0 < abs(BAM.templength(record)) < min_templength) && paironsamestrand(record, invertstrand) && isproperpair(record)) && continue
        c += 1
        id = hash_id ? (is_bitstring ? parse(UInt, BAM.tempname(record); base=2) : hash(@view(record.data[1:BAM.seqname_length(record)]))) : BAM.tempname(record)
        current_read_dict = isread2(record) && ispaired(record) ? reads2 : reads1
        invert = (current_read_dict === reads2 && invertstrand in (:read2, :both)) || (current_read_dict === reads1 && invertstrand in (:read1, :both))
        foundit = false
        if hasxatag(record) && !only_unique
            xa = xatag(record)
            leftest = 0
            leftestpos = BAM.leftposition(record)
            nms = nmtag(record)
            xastrings = split(xa, ";")[1:end-1]
            for (i,xapart) in enumerate(xastrings)
                _, pos, _, nm = split(xapart, ",")
                parse(UInt8, nm) > nms && continue
                intpos = abs(parse(Int, pos))
                intpos < leftestpos && (leftest = i; leftestpos=intpos; nms=nm)
            end
            if leftest != 0 
                alnpart = AlignedPart(xastrings[leftest]; invertstrand=invert)
                foundit = true
            end
        end
        if !foundit
            readstart, readstop, _, readlen = readpositions(record)
            seq_interval = BAM.ispositivestrand(record) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
            st = BAM.ispositivestrand(record) == !invert ? Strand('+') : Strand('-')
            ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), st, AlignmentAnnotation())
            alnpart = AlignedPart(ref_interval, seq_interval, BAM.isprimary(record), astag(record)) #record["AS"]::UInt8)
        end
        id in keys(current_read_dict) ? push!(current_read_dict[id].alns, alnpart) : push!(current_read_dict, id=>AlignedRead([alnpart]))
    end
    close(reader)
    return reads1, reads2
end

function annotate!(features::Features, bam_file::String; only_unique=true, invert_strand=:none, count_key="Count")
    @assert invert_strand in [:read1, :read2, :both, :none]
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
        params(feature)[count_key] = "$c"
    end
    close(reader)
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

function annotate!(alns::Alignments, features::Features; prioritize_type=nothing, overwrite_type=nothing)
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
    end
end

function annotate!(alns::PairedAlignments, features::Features; prioritize_type=nothing, overwrite_type=nothing)
    myiterators = Dict(refn=>GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter, 
        GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refn], Interval("", 1, 1, STRAND_NA, Annotation("", "", Dict{String, String}()))) for refn in refnames(features)) 
    for (alignment1, alignment2) in alns
        for alignedread in (alignment1, alignment2)
            for aln in enumerate(alignedread)
                myiterator = myiterators[refname(aln)]
                myiterator.query = refinterval(aln)
                for feature_interval in myiterator
                    annotate!(aln, feature_interval; prioritize_type=prioritize_type, overwrite_type=overwrite_type) && break
                end
            end
        end
    end
end

function conservedfeatures(features::Features, feature_alignments::Alignments{String}; key_gen=typenamekey)
    my_features = Vector{Interval{Annotation}}()
    for feature in features
        key = key_gen(feature)
        if key in keys(feature_alignments)
            refs = Set(refname(aln) for aln in feature_alignments[key])
            para = Dict(ref=>join(("$(readinterval(aln))" for aln in feature_alignments[key] if refname(aln)==ref), ",") for ref in refs)
            push!(para, "Count"=>"$(length(refs))")
            push!(para, "Name"=>key)
            annot = Annotation("ALN", key, para)
            push!(my_features, Interval(refname(feature), leftposition(feature), rightposition(feature), strand(feature), annot))
        end
    end
    return Features(my_features)
end