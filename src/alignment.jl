function align_mem(in_file::String, out_file::String, genome_file::String; 
    z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")

    tmp_bwa = joinpath(dirname(in_file), "tmp.bwa")
    tmp_view = joinpath(dirname(in_file), "tmp.view")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -v 1 -t 6 $genome_file $in_file`, stdout=tmp_bwa)
    run(cmd)
    cmd = pipeline(`$sam_bin view -u $tmp_bwa`, stdout=tmp_view)
    run(cmd)
    cmd = pipeline(`$sam_bin sort $tmp_view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm(tmp_bwa)
    rm(tmp_view)
end

function align_mem(in_file1::String, in_file2::String, out_file::String, genome_file::String; 
    z_score=100, unpair_panelty=9, unpair_rescue=false, bwa_bin="bwa-mem2", sam_bin="samtools")

    tmp_bwa = joinpath(dirname(in_file1), "tmp.bwa")
    tmp_view = joinpath(dirname(in_file1), "tmp.view")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = unpair_rescue ? 
        pipeline(`$bwa_bin mem -d $z_score -v 1 -U $unpair_panelty -P -t 6 $genome_file $in_file1 $in_file2`, stdout=tmp_bwa) : 
        pipeline(`$bwa_bin mem -d $z_score -v 1 -U $unpair_panelty -t 6 $genome_file $in_file1 $in_file2`, stdout=tmp_bwa)
    run(cmd)
    cmd = pipeline(`$sam_bin view -u $tmp_bwa`, stdout=tmp_view)
    run(cmd)
    cmd = pipeline(`$sam_bin sort $tmp_view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm(tmp_bwa)
    rm(tmp_view)
end

function align_mem(read_files::SingleTypeFiles, genome::Genome; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false)
    tmp_genome = joinpath(dirname(read_files.list[1]), "tmp_genome.fa")
    write(tmp_genome, genome)
    for file in read_files
        out_file = file[1:end-length(read_files.type)] * ".bam"
        (isfile(out_file) && !overwrite_existing) && continue
        align_mem(file, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        isfile(tmp_genome * ending) && rm(tmp_genome * ending)
    end
end

function align_mem(read_files::PairedSingleTypeFiles, genome::Genome; z_score=100, unpair_panelty=9, unpair_rescue=false, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false)
    tmp_genome = joinpath(dirname(read_files.list[1][1]), "tmp_genome.fa")
    write(tmp_genome, genome)
    for (file1, file2) in read_files
        out_file = file1[1:end-length(read_files.type)] * ".bam"
        (isfile(out_file) && !overwrite_existing) && continue
        align_mem(file1, file2, out_file, tmp_genome; z_score=z_score, unpair_panelty=unpair_panelty, unpair_rescue=unpair_rescue, bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        isfile(tmp_genome * ending) && rm(tmp_genome * ending)
    end
end

function align_mem(reads::Reads, genome::Genome, out_file::String; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")
    tmp_reads = joinpath(dirname(out_file), "tmp_reads.fasta")
    tmp_genome = joinpath(dirname(out_file), "tmp_genome.fa")
    write(tmp_genome, genome)
    write(tmp_reads, reads)
    align_mem(tmp_reads, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_reads)
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
end

function align_mem(reads::PairedReads, genome::Genome, out_file::String; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")
    tmp_reads1 = joinpath(dirname(out_file), "temp1.fasta")
    tmp_reads2 = joinpath(dirname(out_file), "temp2.fasta")
    tmp_genome = joinpath(dirname(out_file), "tmp_genome.fa")
    write(tmp_reads1, tmp_reads2, reads)
    write(tmp_genome, genome)
    align_mem(tmp_reads1, tmp_reads2, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_reads1)
    rm(tmp_reads2)
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
end

function local_alignment(reference_sequence::LongDNASeq, query_sequence::LongDNASeq, scoremodel::AffineGapScoreModel)
    res = pairalign(LocalAlignment(), query_sequence, reference_sequence, scoremodel)
    return res
end

struct MyAlignment
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    isprimary::Bool
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

function MyAlignment(xapart::Union{String, SubString{String}}; invertstrand=false)
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = (refstart > 0) == !invertstrand ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = readpositions(cigar)
    seq_interval = (refstart > 0) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    return MyAlignment(ref_interval, seq_interval, false)
end

annotation(aln::MyAlignment) = aln.ref.metadata
annotationname(aln::MyAlignment) = annotation(aln).name
annotationtype(aln::MyAlignment) = annotation(aln).type
annotationoverlap(aln::MyAlignment) = annotation(aln).overlap
refinterval(aln::MyAlignment) = aln.ref
refname(aln::MyAlignment) = aln.ref.seqname
BioGenerics.leftposition(aln::MyAlignment) = aln.ref.first
BioGenerics.rightposition(aln::MyAlignment) = aln.ref.last
GenomicFeatures.strand(aln::MyAlignment) = aln.ref.strand
readinterval(aln::MyAlignment) = aln.seq
hasannotation(aln::MyAlignment) = !isempty(annotation(aln))

@inline function translateddata(data::SubArray{UInt8,1})
    for i in 1:length(data)
        (data[i] == 0x00) && (return data[1:i-1])
    end
end

@inline function nmtag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
      (0x4d == data[i]) & (0x43 == data[i+1]) && (return Int(data[i+2]))
    end
    return nothing
end

function hasxatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        (0x00 == data[i]) & (UInt8('X') == data[i+1]) & (UInt8('A') == data[i+2]) && (return true)
    end
    return false
end

@inline function xatag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        (0x00 == data[i]) & (UInt8('X') == data[i+1]) & (UInt8('A') == data[i+2]) && 
        (return translateddata(@view(data[i+4:end])))
    end
    return nothing
end

function hasxatag(record::BAM.Record)
    return hasxatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
end

function xatag(record::BAM.Record)
    xa = xatag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
    return isnothing(xa) ? nothing : String(xa)
end

function nmtag(record::BAM.Record)
    return nmtag(@view(record.data[BAM.auxdata_position(record):BAM.data_size(record)]))
end

struct AlignedRead
    alns::Vector{MyAlignment}
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
        isempty(annotation(part)) ? println("not annotated.") : println("annotation: $(annotationtype(part)):$(annotationname(part)) [$(annotationoverlap(part))]")
    end
end

function hasannotation(alnread::AlignedRead, annotation_name::String)
    for part in alnread
        annotationname(part) == annotation_name && (return true)
    end
    return false
end

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
                if (annotationname(part) == annotationname(otherpart))
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

struct Alignments <: AlignmentContainer
    dict::Dict{UInt, AlignedRead}
end

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

function Alignments(bam_file::String; min_templength=nothing, only_unique=true, rev_comp=:none)
    alignments1, alignments2 = read_bam(bam_file; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp)
    @assert isempty(alignments2)
    Alignments(alignments1)
end

struct PairedAlignments <: AlignmentContainer
    read1::Dict{UInt, AlignedRead}
    read2::Dict{UInt, AlignedRead}
    keys::Set{UInt}
end

Base.getindex(alignments::PairedAlignments, key::UInt) = (alignments.read1[key], alignments.read2[key])
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

function PairedAlignments(bam_file1::String, bam_file2::String; min_templength=nothing, only_unique=true, rev_comp=:none)
    alignments1, alignments_e1 = read_bam(bam_file1; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp)
    alignments2, alignments_e2 = read_bam(bam_file2; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp)
    @assert isempty(alignments_e1) && isempty(alignments_e2)
    PairedAlignments(alignments1, alignments2, intersect(Set(keys(alignments1)), Set(keys(alignments2))))
end

function PairedAlignments(pebam_file::String; min_templength=nothing, only_unique=true, rev_comp=:none)
    alignments1, alignments2 = read_bam(pebam_file; min_templength=min_templength, only_unique=only_unique, invertstrand=rev_comp)
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

function read_bam(bam_file::String; min_templength=nothing, only_unique=true, invertstrand=:none)
    @assert invertstrand in [:read1, :read2, :both, :none]
    reads1 = Dict{UInt, AlignedRead}()
    reads2 = Dict{UInt, AlignedRead}()
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
        id = is_bitstring ? parse(UInt, BAM.tempname(record); base=2) : hash(@view(record.data[1:BAM.seqname_length(record)]))
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
                chr, pos, cigar, nm = split(xapart, ",")
                parse(Int, nm) > nms && continue
                intpos = abs(parse(Int, pos))
                intpos < leftestpos && (leftest = i; leftestpos=intpos; nms=nm)
            end
            if leftest != 0 
                alnpart = MyAlignment(xastrings[leftest]; invertstrand=invert)
                foundit = true
            end
        end
        if !foundit
            readstart, readstop, relrefstop, readlen = readpositions(record)
            seq_interval = BAM.ispositivestrand(record) ? (readstart:readstop) : (readlen-readstop+1:readlen-readstart+1)
            strand = BAM.ispositivestrand(record) == !invert ? Strand('+') : Strand('-')
            ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), strand, AlignmentAnnotation())
            alnpart = MyAlignment(ref_interval, seq_interval, BAM.isprimary(record))
        end
        id in keys(current_read_dict) ? push!(current_read_dict[id].alns, alnpart) : push!(current_read_dict, id=>AlignedRead([alnpart]))
    end
    close(reader)
    return reads1, reads2
end

function annotate!(alns::Alignments, features::Features; prioritize_type=nothing)
    myiterators = Dict(refn=>GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter, 
        GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refn], Interval("", 1, 1, STRAND_NA, Annotation("", ""))) for refn in refnames(features)) 
    for alignedread in alns
        for (i,alignment) in enumerate(alignedread)
            foundpriority = false
            myiterator = myiterators[refname(alignment)]
            myiterator.query = refinterval(alignment)
            for feature_interval in myiterator
                olp = round(UInt8, (overlapdistance(feature_interval, refinterval(alignment)) / length(refinterval(alignment))) * 100)
                foundpriority = (!isnothing(prioritize_type) && (annotationtype(feature_interval) == prioritize_type))
                if foundpriority || annotationoverlap(alignment) < olp
                    annotation(alignedread[i]).type = feature_interval.metadata.type
                    annotation(alignedread[i]).name = feature_interval.metadata.name
                    annotation(alignedread[i]).overlap = olp
                    foundpriority && break
                end
            end
        end
    end
end

function annotate!(alns::PairedAlignments, features::Features; prioritize_type=nothing, overwrite_type=nothing)
    myiterators = Dict(refn=>GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter, 
        GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refn], Interval("", 1, 1, STRAND_NA, Annotation("", ""))) for refn in refnames(features)) 
    for (alignment1, alignment2) in alns
        for alignedread in (alignment1, alignment2)
            for (i,alignment) in enumerate(alignedread)
                foundpriority = false
                myiterator = myiterators[refname(alignment)]
                myiterator.query = refinterval(alignment)
                for feature_interval in myiterator
                    olp = round(UInt8, (overlapdistance(feature_interval, refinterval(alignment)) / length(refinterval(alignment))) * 100)
                    foundpriority = (!isnothing(prioritize_type) && (annotationtype(feature_interval) == prioritize_type))
                    overwrite = (!isnothing(overwrite_type) && (annotation(alignedread[i]).type == overwrite_type))
                    if foundpriority || annotationoverlap(alignment) < olp || overwrite
                        annotation(alignedread[i]).type = feature_interval.metadata.type
                        annotation(alignedread[i]).name = feature_interval.metadata.name
                        annotation(alignedread[i]).overlap = olp
                        foundpriority && break
                    end
                end
            end
        end
    end
end