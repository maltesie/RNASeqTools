#!/usr/bin/env julia

function bam_chromosome_lengths(reader::BAM.Reader)
    chr_lengths = Int[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_lengths, parse(Int, meta["LN"]))
    end
    return chr_lengths
end

function bam_chromosome_names(reader::BAM.Reader)
    chr_names = String[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_names, meta["SN"])
    end
    return chr_names
end

function samevalueintervals(d::Vector{UInt}, index::Vector{Int})
    rindex::Int = 1
    lindex::Int = 1
    n_unique = length(Set(d))
    ranges = Vector{UnitRange{Int}}(undef, n_unique)
    subd = view(d, index)
    for i::Int in 1:(length(d)-1)
        if subd[i] != subd[i+1]
            ranges[rindex] = lindex:i
            lindex = i+1
            rindex += 1
        end
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
            if BioAlignments.isinsertop(op) || (op == BioAlignments.OP_HARD_CLIP)
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
        if BioAlignments.isinsertop(op) || (op == BioAlignments.OP_HARD_CLIP)
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
        (data[i] == 0x00) && (return @view data[1:i-1])
    end
end

@inline function xatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        ((0x00 == data[i]) & (0x58 == data[i+1]) & (0x41 == data[i+2])) &&
        (return translateddata(@view(data[i+4:end])))
    end
    return nothing
end

function hasxatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        ((0x00 == data[i]) & (0x58 == data[i+1]) & (0x41 == data[i+2])) && (return true)
    end
    return false
end

function hasxatag(record::BAM.Record)
    return hasxatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
end

xatag(record::BAM.Record) = xatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))

function astag(record::BAM.Record)
    data = @view(record.data[BAM.auxdata_position(record):BAM.data_size(record)])
    for i in 1:length(data)-4
        (0x00 == data[i]) & (0x41 == data[i+1]) & (0x53 == data[i+2]) && (return data[i+4])
    end
end

const nmdict = Dict{UInt8,Int}(UInt8('C') => 24, UInt8('S') => 16, UInt8('I') => 0)
function nmtag(record::BAM.Record)
    ((record.data[BAM.auxdata_position(record)] == UInt8('N')) && (record.data[BAM.auxdata_position(record)+1] == UInt8('M'))) || throw("auxdata does not start with NM tag.")
    t = nmdict[record.data[BAM.auxdata_position(record)+2]]
    return unsafe_load(Ptr{UInt32}(pointer(record.data, BAM.auxdata_position(record)+3))) << t >> t
end

"""
AlignedReads(bam_file::String)

    Constructor for AlignedReads. Reads alignment records and stores alignment positions on reference and read and additional information
    such as reference name, strand and edit distance of the alignment sorted according to the read they originate from.
"""
function AlignedReads(bam_file::String; include_secondary_alignments=true, include_alternative_alignments=false,
        is_reverse_complement=false, min_mapping_quality=0)
    record = BAM.Record()
    ns = Vector{UInt}(undef, 1000000)
    ls = Vector{Int}(undef, 1000000)
    rs = Vector{Int}(undef, 1000000)
    is = Vector{String}(undef, 1000000)
    ss = Vector{Strand}(undef, 1000000)
    nms = Vector{UInt32}(undef, 1000000)
    rls = Vector{Int}(undef, 1000000)
    rrs = Vector{Int}(undef, 1000000)
    rds = Vector{Symbol}(undef, 1000000)
    mmq = UInt8(min_mapping_quality)
    chromosome_list = Vector{Tuple{String,Int}}()
    BAM.Reader(open(bam_file)) do reader
        append!(chromosome_list, [n for n in zip(bam_chromosome_names(reader), bam_chromosome_lengths(reader))])
        index::Int = 0
        while !eof(reader)
            read!(reader, record)
            BAM.mappingquality(record) >= mmq || continue
            has_alternatives = hasxatag(record)
            (!BAM.ismapped(record) || (!isprimary(record) && !include_secondary_alignments) || (has_alternatives && !include_alternative_alignments)) && continue
            current_read::Symbol = (isread2(record) != is_reverse_complement) ? :read2 : :read1
            index += 1
            xas = has_alternatives ? StringView(xatag(record)) : nothing
            ((index + (has_alternatives ? count(';', xas) : 0)) > length(ns)) && resize!.((ns, ls, rs, is, ss, rls, rrs, rds, nms), length(ns)+1000000)
            n = hash(@view record.data[1:Int(BAM.seqname_length(record))-1])
            (l::Int, r::Int) = (BAM.leftposition(record), BAM.rightposition(record))
            (ref::String, s::Strand) = (BAM.refname(record), BAM.ispositivestrand(record) != (current_read === :read2) ? STRAND_POS : STRAND_NEG)
            nm::UInt32 = nmtag(record)
            if has_alternatives
                v::Int = 0
                for _ in 1:count(';', xas)
                    nextv = findnext(';', xas, v+1)
                    ap = AlignedInterval(xas[(v+1):(nextv-1)]; read=current_read)
                    ns[index], ls[index], rs[index], is[index], ss[index], rls[index], rrs[index], rds[index], nms[index] =
                    n, leftposition(ap), rightposition(ap), refname(ap), strand(ap), first(readrange(ap)), last(readrange(ap)), current_read, editdistance(ap)
                    index += 1
                    v = nextv
                end
            end
            readstart, readstop, _, readlen = readpositions(record)
            (rl, rr) = (BAM.ispositivestrand(record) != (current_read === :read2)) ? (readstart,readstop) : (readlen-readstop+1,readlen-readstart+1)
            ns[index], ls[index], rs[index], is[index], ss[index], rls[index], rrs[index], rds[index], nms[index] = n, l, r, ref, s, rl, rr, current_read, nm
        end
        resize!.((ns, ls, rs, is, ss, rls, rrs, rds, nms), index)
    end
    nindex = sortperm(ns)
    ranges = samevalueintervals(ns, nindex)
    pindex = partsindex(ranges, nindex, rls, rds)
    return AlignedReads(chromosome_list, ns, ls, rs, rls, rrs, rds, nms, is, ss,
                        Vector{String}(undef,length(ns)), Vector{String}(undef,length(ns)), zeros(UInt8,length(ns)),
                        fill(0xff,length(ns)), fill(0xff,length(ns)), pindex, zeros(Bool, length(ns)), ranges)
end

Base.length(alns::AlignedReads) = length(alns.ranges)
ninterval(alns::AlignedReads) = length(alns.tempnames)
nread(alns::AlignedReads) = length(alns)

@inline function Base.iterate(alns::AlignedReads)
    return isempty(alns.ranges) ? nothing : (alns[1], 2)
end
@inline function Base.iterate(alns::AlignedReads, state::Int)
    return state > length(alns.ranges) ? nothing : (alns[state], state+1)
end

function Base.filter!(alns::AlignedReads, tempnames::Set{UInt})
    bitindex = [alns.tempnames[alns.pindex[first(r)]] in tempnames for r in alns.ranges]
    n_seqs = sum(bitindex)
    alns.ranges[1:n_seqs] = alns.ranges[bitindex]
    resize!(alns.ranges, n_seqs)
end

function sync!(seqs::Sequences, alns::AlignedReads)
    filter!(seqs, Set(alns.tempnames))
    sortbyname!(seqs)
    filter!(alns, Set(seqs.tempnames))
end
sync!(alns::AlignedReads, seqs::Sequences) = sync!(seqs, alns)

function print_reads(alns::AlignedReads; n=-1, only_chimeric=false, filter_name=nothing, filter_type=nothing)
    c = 0
    for aln in alns
        c == n && break
        isnothing(filter_name) || hasname(aln, filter_name) || continue
        isnothing(filter_type) || hastype(aln, filter_type) || continue
        only_chimeric && !ischimeric(aln) && continue
        show(aln)
        println("---")
        c += 1
    end
end

function AlignedInterval(alns::AlignedReads, i::Int)
    i_p = alns.pindex[i]
    AlignedInterval(
        Interval(alns.refnames[i_p], alns.leftpos[i_p], alns.rightpos[i_p], alns.strands[i_p],
            alns.annotated[i_p] ? AlignmentAnnotation(alns.antypes[i_p], alns.annames[i_p], alns.anols[i_p]) : AlignmentAnnotation()
        ),
        alns.read_leftpos[i_p]:alns.read_rightpos[i_p],
        alns.nms[i_p],
        alns.reads[i_p]
    )
end

AlignedInterval(alnpart::AlignedInterval; new_name::Union{Nothing, String}=nothing, new_type::Union{Nothing, String}=nothing) =
AlignedInterval(
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

Base.getindex(alns::AlignedReads, i::Int) = AlignedRead(alns.ranges[i], alns)
Base.getindex(alns::AlignedReads, r::UnitRange{Int}) = AlignedReads(alns.chroms, alns.tempnames, alns.leftpos, alns.rightpos, alns.read_leftpos, alns.read_rightpos,
                                                                alns.reads, alns.nms, alns.refnames, alns.strands, alns.annames, alns.antypes,
                                                                alns.anols, alns.anleftrel, alns.anrightrel, alns.pindex, alns.annotated, alns.ranges[r])

"""
AlignedInterval(xapart::StringView)

    Constructor for the AlignedInterval struct. Builds AlignedInterval from a XA string containing alternative mappings.
    Inverts strand, if `check_invert::Bool` is true.
"""
function AlignedInterval(xapart::StringView; read=:read1)
    chr, pos, cigar = (i for (i,c) in enumerate(xapart) if c === ',')
    refstart = parse(Int, xapart[(chr+1):(pos-1)])
    is_forward = ((refstart > 0) != (read===:read2))
    strand = is_forward ? STRAND_POS : STRAND_NEG
    readstart, readstop, relrefstop, readlen = readpositions(xapart[(pos+1):(cigar-1)])
    seq_interval = is_forward ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
    refstart *= sign(refstart)
    ref_interval = Interval(xapart[chr+1:pos-1], refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    return AlignedInterval(ref_interval, seq_interval, parse(UInt32, xapart[(cigar+1):end]), read)
end

"""
    Base.summarize(part::AlignedInterval, readseq::LongDNA)::IO

Generates string with information on the AlignedInterval combined with the read sequence.
"""
function summarize(part::AlignedInterval, readseqpart::Union{LongSubSeq, LongSequence})
    s = "[$(first(part.seq)), $(last(part.seq))] on $(part.read) - "
    s *= "[$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) ($(part.ref.strand)) "
    s *= "with edit distance $(part.nms) - "
    s *= isempty(annotation(part)) ? "not annotated." : "$(type(part)):$(name(part)) ($(overlap(part))% in annotation)\n"
    s *= "\t$(first(part.seq)) - " * join(["$(readseqpart[l:r-1])" for (l,r) in partition(1:100:length(readseqpart), 2, 1)], "\n\t") *
        (isempty([true for (l,r) in partition(1:100:length(readseqpart), 2, 1)]) ? "" : "\n\t") *
        "$(readseqpart[last(1:100:length(readseqpart)):length(readseqpart)]) - $(last(part.seq))"
    return s
end

"""
    Base.summarize(part::AlignedInterval)::IO

Generates string with information on the AlignedInterval.
"""
function summarize(part::AlignedInterval)
    s = "[$(first(part.seq)), $(last(part.seq))] on $(part.read) - "
    s *= "[$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) ($(part.ref.strand)) "
    s *= "with edit distance $(part.nms) - "
    s *= isempty(annotation(part)) ? "not annotated." : "$(type(part)):$(name(part)) ($(overlap(part))% in annotation)"
    return s
end

"""
    Base.show(part::AlignedInterval)::IO

Function for structured printing of the content of AlignedInterval
"""
function Base.show(io::IO, part::AlignedInterval, readseqpart::Union{LongSubSeq, LongSequence})
    print(io, summarize(part, readseqpart))
end

"""
    Base.show(part::AlignedInterval)::IO

Function for structured printing of the content of AlignedInterval
"""
function Base.show(io::IO, part::AlignedInterval)
    print(io, summarize(part))
end

"""
    read(aln::AlignedInterval)::Symbol

Returns the read from which the alignment comes. If is_reverse_complement is true, the reads will be
inverted (:read1 will be :read2 and vice versa)....
"""
sameread(aln1::AlignedInterval, aln2::AlignedInterval) = aln1.read === aln2.read

"""
    annotation(aln::AlignedInterval)::AlignmentAnnotation

Returns the AlignmentAnnotation struct that is used as the metadata of the Interval{AlignmentAnnotation} and contains
information on the reference sequence part of the alignment.
"""
annotation(aln::AlignedInterval) = aln.ref.metadata

"""
    name(aln::AlignedInterval)::String

Returns the annotation's name as found in the .gff annotation file by the name_key argument of the Features constructor.
Availble only after using annotate! on AlignedReads.
"""
name(aln::AlignedInterval) = name(annotation(aln))

"""
    type(aln::AlignedInterval)::String

Returns the annotation's type as found in the .gff annotation file as the feature type.
Availble only after using annotate! on AlignedReads.
"""
type(aln::AlignedInterval) = type(annotation(aln))

"""
    overlap(aln::AlignedInterval)::UInt8

Returns the percentage of the sequence on the read contained within the annotation assigned to it.
Availble only after using annotate! on AlignedReads.
"""
overlap(aln::AlignedInterval) = overlap(annotation(aln))

"""
    editdistance(aln::AlignedInterval)::UInt8

Returns the number of edit operations between the aligned sequence and the reference (NM tag in SAM/BAM).
"""
editdistance(aln::AlignedInterval) = aln.nms

"""
    refinterval(aln::AlignedInterval)::Interval{AlignmentAnnotation}

Returns the Interval{AlignmentAnnotation} that contains information on the reference sequence part of the alignment.
Availble only after using annotate! on AlignedReads.
"""
refinterval(aln::AlignedInterval) = aln.ref

"""
    refname(aln::AlignedInterval)::String

Returns the id of the reference sequence the annotation comes from as found in the .gff annotation file as the chromosome name.
Availble only after using annotate! on AlignedReads.
"""
refname(aln::AlignedInterval) = aln.ref.seqname

"""
    refrange(aln::AlignedInterval)::UnitRange

Returns the Interval of the alignment on the reference sequence as a UnitRange.
"""
refrange(aln::AlignedInterval) = leftposition(aln):rightposition(aln)

"""
    refsequence(aln::AlignedInterval, genome::Genome)::LongDNA

Returns the part of the reference sequence that the aln::AlignedInterval belongs to.
"""
refsequence(aln::AlignedInterval, genome::Genome)::LongDNA{4} = genome[refname(aln)][refrange(aln)]

"""
    leftposition(aln::AlignedInterval)::Int

Returns the leftmost position of the alignment on the reference sequence.
Availble only after using annotate! on AlignedReads.
"""
BioGenerics.leftposition(aln::AlignedInterval) = aln.ref.first

"""
    rightposition(aln::AlignedInterval)::Int

Returns the rightmost position of the alignment on the reference sequence.
Availble only after using annotate! on AlignedReads.
"""
BioGenerics.rightposition(aln::AlignedInterval) = aln.ref.last

"""
    strand(aln::AlignedInterval)::Strand

Returns the strand of the annotation.
Availble only after using annotate! on AlignedReads.
"""
GenomicFeatures.strand(aln::AlignedInterval) = aln.ref.strand

"""
    strand(aln::AlignedInterval)::Strand

Returns the strand of the annotation.
Availble only after using annotate! on AlignedReads.
"""
ispositivestrand(aln::AlignedInterval) = aln.ref.strand === STRAND_POS

"""
    readrange(aln::AlignedInterval)::UnitRange

Returns the Interval of the alignment on the read sequence as a UnitRange.
"""
readrange(aln::AlignedInterval) = aln.seq

"""
    nms(aln::AlignedInterval)::UInt32

Returns the edit distance of the AlignedInterval.
"""
nms(aln::AlignedInterval) = aln.nms

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

Checks if an annotation was assigned to the AlignedInterval.
"""
hasannotation(aln::AlignedInterval) = !isempty(annotation(aln))

"""
    sameannotation(aln1::AlignedInterval, aln2::AlignedInterval)::Bool

Checks if two AlignedReads share the same name and type in their annotation.
"""
sameannotation(aln1::AlignedInterval, aln2::AlignedInterval) = samename(aln1, aln2) && sametype(aln1, aln2)

"""
    samename(aln1::AlignedInterval, aln2::AlignedInterval)::Bool

Checks if two AlignedReads share the same name and type in their annotation.
"""
samename(aln1::AlignedInterval, aln2::AlignedInterval) = name(aln1) == name(aln2)

"""
    sametype(aln1::AlignedInterval, aln2::AlignedInterval)::Bool

Checks if two AlignedReads share the same name and type in their annotation.
"""
sametype(aln1::AlignedInterval, aln2::AlignedInterval) = type(aln1) == type(aln2)

"""
    distanceonread(aln1::AlignedInterval, aln2::AlignedInterval)::Int

Computes the distance between two alignments on the same read. Returns negative numbers
if the two parts overlap.
"""
distanceonread(aln1::AlignedInterval, aln2::AlignedInterval) = max(first(aln1.seq), first(aln2.seq)) - min(last(aln1.seq), last(aln2.seq))

isoverlapping(aln1::AlignedInterval, aln2::AlignedInterval) = strand(aln1) === strand(aln2) ? isoverlapping(aln1.ref, aln2.ref) : false

isfirstread(part::AlignedInterval) = part.read === :read1

@inline function Base.iterate(alnread::AlignedRead)
    return isnothing(alnread.range) ? nothing : (alnread[1], 2)
end
@inline function Base.iterate(alnread::AlignedRead, state::Int)
    return state > length(alnread.range) ? nothing : (alnread[state], state+1)
end

Base.lastindex(alnread::AlignedRead) = length(alnread.range)
function Base.getindex(alnread::AlignedRead, i::Int)
    f::Int = i+first(alnread.range)-1
    return f in alnread.range ? AlignedInterval(alnread.alns,  f) : throw(BoundsError(alnread, i))
end
Base.getindex(alnread::AlignedRead, r::UnitRange{Int}) = AlignedRead(alnread.range[r], alnread.alns)
Base.getindex(alnread::AlignedRead, b::Vector{Bool}) = [alnread.alns[alnread.prindex[index]] for (i::Int,index::Int) in enumerate(alnread.range) if b[i]]
Base.isempty(alnread::AlignedRead) = isempty(alnread.range)
Base.length(alnread::AlignedRead) = length(alnread.range)

readid(alnread::AlignedRead) = alnread.alns.tempnames[alnread.alns.pindex[first(alnread.range)]]
parts(alnread::AlignedRead) = [AlignedInterval(alnread.alns, i) for i in alnread.range]
annotatedparts(alnread::AlignedRead) = [AlignedInterval(alnread.alns, i) for i in alnread.range if alnread.alns.annotated[alnread.alns.pindex[i]]]
typein(alnread::AlignedRead, types::Vector{String}) = any(alnread.alns.antypes[i] in types for i::Int in view(alnread.alns.pindex, alnread.range) if alnread.alns.annotated[i])
hastype(alnread::AlignedRead, t::String) = any(alnread.alns.antypes[i] === t for i::Int in view(alnread.alns.pindex, alnread.range) if alnread.alns.annotated[i])
namein(alnread::AlignedRead, names::Vector{String}) = any(alnread.alns.annames[i] in names for i::Int in view(alnread.alns.pindex, alnread.range) if alnread.alns.annotated[i])
hasname(alnread::AlignedRead, n::String) = any(alnread.alns.annames[i] === n for i::Int in view(alnread.alns.pindex, alnread.range) if alnread.alns.annotated[i])
hasannotation(alnread::AlignedRead) = any(alnread.alns.annotated[i] for i in view(alnread.alns.pindex, alnread.range))
annotatedcount(alnread::AlignedRead) = sum(alnread.alns.annotated[i] for i in view(alnread.alns.pindex, alnread.range))
nnames(alnread::AlignedRead) = length(Set(name(part) for part in alnread))
ntypes(alnread::AlignedRead) = length(Set(type(part) for part in alnread))

function BioGenerics.leftposition(alnread::AlignedRead)
    check_refname = alnread.alns.refnames[alnread.alns.pindex[first(alnread.range)]]
    all(v .== check_refname for v in view(alnread.alns.refnames, view(alnread.alns.pindex, alnread.range))) || throw(AssertionError("AlignmentParts are not on the same reference sequence."))
    minimum(view(alnread.alns.leftpos, view(alnread.alns.pindex, alnread.range)))
end
function BioGenerics.rightposition(alnread::AlignedRead)
    check_refname = alnread.alns.refnames[alnread.alns.pindex[first(alnread.range)]]
    all(v .== check_refname for v in view(alnread.alns.refnames, view(alnread.alns.pindex, alnread.range))) || throw(AssertionError("AlignmentParts are not on the same reference sequence."))
    maximum(view(alnread.alns.rightpos, view(alnread.alns.pindex, alnread.range)))
end

function GenomicFeatures.strand(alnread::AlignedRead)
    length(alnread) > 0 || (return STRAND_NA)
    check_strand = alnread.alns.strands[alnread.alns.pindex[first(alnread.range)]]
    return all(s === check_strand for s in view(alnread.alns.strands, view(alnread.alns.pindex, alnread.range))) ? check_strand : STRAND_BOTH
end

function ispositivestrand(alnread::AlignedRead)
    s = strand(alnread)
    s === STRAND_NA && throw(AssertionError("Empty Alignment does not have a strand."))
    s === STRAND_BOTH && throw(AssertionError("AlignedRead parts are not on the same strand."))
    return s === STRAND_POS
end

summarize(alnread::AlignedRead) = (ischimeric(alnread) ? "Chimeric" : "Single") *
                                    " Alignment with $(length(alnread)) part(s):\n   " *
                                    join([summarize(part) for part in alnread], "\n   ") * "\n"

function summarize(alnread::AlignedRead, readseq::Union{LongSubSeq, LongSequence})
    length(Set(alnread.alns.reads[i] for i in view(alnread.alns.pindex, alnread.range))) === 1 || throw(AssertionError("Need two sequences for paired end reads!"))
    (typeof(alnread.alns) == AlignedReads{String} ? "[$(alnread.alns.tempnames[alnread.alns.pindex[first(alnread.range)]])] " : "") *
    (ischimeric(alnread) ? "Chimeric" : "Single") *
    " Alignment with $(length(alnread)) part(s) on $(length(readseq)) nt read:\n    1:\t" *
    join([summarize(part, readseq[readrange(part)])*(i < length(alnread) ? "\n    $(i+1):\t" : "") for (i,part) in enumerate(alnread)]) * "\n"
end

function summarize(alnread::AlignedRead, read1seq::Union{LongSubSeq, LongSequence}, read2seq::Union{LongSubSeq, LongSequence})
    (typeof(alnread.alns) == AlignedReads{String} ? "[$(alnread.alns.tempnames[alnread.alns.pindex[first(alnread.range)]])] " : "") *
    (ischimeric(alnread) ?  "Chimeric" : "Single") *
    " Alignment with $(length(alnread)) part(s) on $(length(read1seq)) nt read1 and $(length(read2seq)) nt read2:\n    1:\t" *
    join([summarize(part, isfirstread(part) ? read1seq[readrange(part)] : read2seq[readrange(part)])*(i < length(alnread) ? "\n    $(i+1):\t" : "") for (i,part) in enumerate(alnread)]) * "\n"
end

function Base.show(io::IO, alnread::AlignedRead)
    print(io, summarize(alnread))
end
function Base.show(io::IO, alnread::AlignedRead, readseq::Union{LongSubSeq, LongSequence})
    print(io, summarize(alnread, readseq))
end
function Base.show(io::IO, alnread::AlignedRead, read1seq::Union{LongSubSeq, LongSequence}, read2seq::Union{LongSubSeq, LongSequence})
    print(io, summarize(alnread, read1seq, read2seq))
end

"""
    overlapdistance(i1::Interval, i2::Interval)::Float64

Returns the negative distance between two AlignedReads on the reference sequence or the positive overlap if it exists.
Returns -Inf if the alignments do not share the same reference id or lie on different strands.
"""
function overlapdistance(i1::Interval{T}, i2::Interval{I})::Float64 where {T,I}
    (!(strand(i1) === strand(i2)) || !(refname(i1) === refname(i2))) && return -Inf
    return min(
        min(rightposition(i1) - leftposition(i1) + 1, rightposition(i2) - leftposition(i2) + 1),
        min(rightposition(i1) - leftposition(i2) + 1, rightposition(i2) - leftposition(i1) + 1)
    )
end

function distance(l1::Int, r1::Int, l2::Int, r2::Int; check_order=false)::Float64
    check_order && l1 > l2 && return Inf
    l2>r1 && return l2-r1+1
    l1>r2 && return l1-r2+1
    return 0
end

"""
    distance(i1::Interval, i2::Interval)::Float64

Returns the distance between two AlignedReads on the reference sequence. Returns Inf if the alignments do not share
the same reference id or lie on different strands
"""
function distance(i1::Interval, i2::Interval; check_order=false)::Float64
    ((refname(i1) != refname(i2)) || (strand(i1) != strand(i2))) && return Inf
    d::Float64 = distance(leftposition(i1), rightposition(i1), leftposition(i2), rightposition(i2))
    return check_order && d>0 && first(i1) > first(i2) ? Inf : d
end

function nchimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true, check_order=false)
    length(alnread) > 1 || return 0
    c = 0
    for i in alnread.range, j in i+1:last(alnread.range)
        c += ischimeric(AlignedInterval(alnread.alns, i), AlignedInterval(alnread.alns, j); min_distance=min_distance, check_annotation=check_annotation, check_order=check_order)
    end
    return c
end

function ischimeric(alnread::AlignedRead; min_distance=1000, check_annotation=true, check_order=false)
    length(alnread) > 1 || return false
    for i in alnread.range, j in i+1:last(alnread.range)
        ischimeric(AlignedInterval(alnread.alns, i), AlignedInterval(alnread.alns, j); min_distance=min_distance, check_annotation=check_annotation, check_order=check_order) && return true
    end
    return false
end

function ischimeric(part1::AlignedInterval, part2::AlignedInterval; min_distance=1000, check_annotation=true, check_order=false)
    check_annotation && hasannotation(part1) && hasannotation(part2) && (name(part1) == name(part2)) && (return false)
    return distance(refinterval(part1), refinterval(part2); check_order=check_order) > min_distance
end

function ismulti(alnread::AlignedRead; method=:distance)
    return method === :annotation ? (nnames(alnread) >= 3) :
        method === :distance ? (nchimeric(alnread; check_annotation=false) >= 3) :
        method === :both ? (nchimeric(alnread; check_annotation=true) >= 3) :
        throw(AssertionError("method must be :distance, :annotation or :both"))
end

summarize(alns::AlignedReads) = "$(typeof(alns)) with $(ninterval(alns)) aligned parts from $(nread(alns)) reads."

Base.show(io::IO, alns::AlignedReads) = println(io, summarize(alns))

function Base.empty!(alns::AlignedReads)
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
    empty!(alns.pindex)
    empty!(alns.annotated)
    empty!(alns.reads)
end

nannotated(alns::AlignedReads) = sum(alns.annotated)

annotate_filter(a::Interval{Annotation}, b::Interval{Annotation}) = strand(a) === strand(b)
function annotate!(alns::AlignedReads, features::Features{Annotation}; prioritize_type=nothing, min_prioritize_overlap=80, overwrite_type=nothing)
    ann = Annotation()
    intersection = IntervalTrees.Intersection{Int, Interval{Annotation}, 64}()
    for i::Int in 1:ninterval(alns)
        alns.refnames[i] in keys(features.list.trees) || continue
        query = Interval(alns.refnames[i], alns.leftpos[i], alns.rightpos[i], alns.strands[i], ann)
        IntervalTrees.firstintersection!(features.list.trees[alns.refnames[i]].root, query, nothing, intersection, annotate_filter)
        while intersection.index > 0
            feature_interval = intersection.node.entries[intersection.index]
            olp = round(UInt8, (min(feature_interval.last, alns.rightpos[i]) - max(feature_interval.first, alns.leftpos[i]) + 1) /
                                (alns.rightpos[i] - alns.leftpos[i] + 1) * 100)

            priority = !isnothing(prioritize_type) && (type(feature_interval) === prioritize_type) && (olp > min_prioritize_overlap)
            overwrite = !isnothing(overwrite_type) && alns.annotated[i] && (alns.antypes[i] === overwrite_type)

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
                alns.annotated[i] = true
                priority && break
            end
            IntervalTrees.nextintersection!(intersection.node, intersection.index, query, intersection, annotate_filter)
        end
    end
end

Base.getindex(genome::Genome, ap::AlignedInterval) = refsequence(ap, genome)
Base.getindex(sequence::LongDNA{4}, ap::AlignedInterval) = view(sequence, readrange(ap))

mutable struct BAMIterator
    reader::BAM.Reader
    record::BAM.Record
    stop::Union{Nothing,Int}
end

@inline Base.iterate(it::BAMIterator, state=1) = (eof(it.reader) || (!isnothing(it.stop) && state>it.stop)) ? (close(it.reader); nothing) : (read!(it.reader, it.record), state+1)
eachbamrecord(fname::String; stopat=nothing) = BAMIterator(BAM.Reader(open(fname)), BAM.Record(), stopat)