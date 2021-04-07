function align_backtrack(in_file::String, out_file::String, genome_file::String; 
    max_miss=2, bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin aln -n $max_miss -t 6 -R 500 $genome_file $in_file`, stdout="tmp.sai")
    run(cmd)
    cmd = pipeline(`$bwa_bin samse $genome_file tmp.sai $in_file`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.sai")
    rm("tmp.bwa")
    rm("tmp.view")
end

function align_backtrack(in_file1::String, in_file2::String, out_file::String, genome_file::String; 
    max_miss=2, bwa_bin="bwa", sam_bin="samtools")
        
    cmd = pipeline(`./bin/bwa index -a is $genome_file`, stdout=nothing)
    run(cmd)
    cmd = pipeline(`./bin/bwa aln -n $max_miss -t 6 -R 500 $genome_file $in_file1`, stdout="tmp1.sai")
    run(cmd)
    cmd = pipeline(`./bin/bwa aln -n $max_miss -t 6 -R 500 $genome_file $in_file2`, stdout="tmp2.sai")
    run(cmd)
    cmd = pipeline(`./bin/bwa sampe -a 1500 -P $genome_file tmp1.sai tmp2.sai $in_file1 $in_file2`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp1.sai")
    rm("tmp2.sai")
    rm("tmp.bwa")
    rm("tmp.view")
end

function align_backtrack(reads::Reads, out_file::String, genome_file::String; max_miss=2, bwa_bin="bwa", sam_bin="samtools")
    tmp_file = joinpath(dirname(out_file), "temp.fasta")
    write(tmp_file, reads)
    align_backtrack(tmp_file, out_file, genome_file; max_miss=max_miss, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_file)
end

function align_mem(in_file::String, out_file::String, genome_file::String; 
    z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -v 1 -t 6 $genome_file $in_file`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.bwa")
    rm("tmp.view")
end

function align_mem(in_file1::String, in_file2::String, out_file::String, genome_file::String; 
    z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -v 1 -t 6 $genome_file $in_file1 $in_file2`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.bwa")
    rm("tmp.view")
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

function align_mem(read_files::PairedSingleTypeFiles, genome::Genome; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools", overwrite_existing=false)
    tmp_genome = joinpath(dirname(read_files.list[1][1]), "tmp_genome.fa")
    write(tmp_genome, genome)
    for (file1, file2) in read_files
        out_file = file1[1:end-length(read_files.type)] * ".bam"
        (isfile(out_file) && !overwrite_existing) && continue
        align_mem(file1, file2, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
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

struct Alignment
    ref::Interval{AlignmentAnnotation}
    seq::Interval{Nothing}
    isprimary::Bool
end

function positions(cigar::AbstractString)
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

function Alignment(xapart::Union{String, SubString{String}})
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = refstart > 0 ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = positions(cigar)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, AlignmentAnnotation())
    seq_interval = strand === STRAND_POS ? Interval("read", readstart, readstop) : Interval("read", readlen-readstop+1, readlen-readstart+1)
    return Alignment(ref_interval, seq_interval, false)
end

annotation(aln::Alignment) = aln.ref.metadata
annotationname(aln::Alignment) = annotation(aln).name
annotationtype(aln::Alignment) = annotation(aln).type
annotationoverlap(aln::Alignment) = annotation(aln).overlap
refinterval(aln::Alignment) = aln.ref
readinterval(aln::Alignment) = aln.seq
hasannotation(aln::Alignment) = !isempty(annotation(aln))

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
    alns::Vector{Alignment}
end

function AlignedRead(record::BAM.Record; invertstrand=false)
    BAM.ismapped(record) || return AlignedRead([])
    xa = xatag(record)
    if !isnothing(xa)
        leftest = 0
        leftestpos = BAM.leftposition(record)
        nms = nmtag(record)
        xastrings = split(xa, ";")[1:end-1]
        for (i,xapart) in enumerate(xastrings)
            chr, pos, cigar, nm = split(xapart, ",")
            parse(Int, nm) > nms && continue
            intpos = abs(parse(Int, pos))
            intpos < leftestpos && (leftest = i; leftestpos=intpos)
        end
        leftest != 0 && (return AlignedRead([Alignment(xastrings[leftest])]))
    end
    readstart, readstop, relrefstop, readlen = positions(BAM.cigar(record))
    strand = BAM.ispositivestrand(record) ? Strand('+') : Strand('-')
    ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), strand, AlignmentAnnotation())
    seq_interval = strand === STRAND_POS ? Interval("read", readstart, readstop) : Interval("read", readlen-readstop+1, readlen-readstart+1)
    aln_part = Alignment(ref_interval, seq_interval, BAM.isprimary(record))
    return AlignedRead([aln_part])
end

alignments(alnread::AlignedRead) = alnread.alns
count(alnread::AlignedRead) = length(alnread.alns)
merge!(alnread1::AlignedRead, alnread2::AlignedRead) = append!(alnread1.alns, alnread2.alns)
hasannotation(alnread::AlignedRead) = any([!isempty(annotation(alnpart)) for alnpart in alnread])
isfullyannotated(alnread::AlignedRead) = all([!isempty(annotation(alnpart)) for alnpart in alnread])

Base.isempty(alnread::AlignedRead) = isempty(alnread.alns)
Base.iterate(alnread::AlignedRead) = iterate(alnread.alns)
Base.iterate(alnread::AlignedRead, state::Int) = iterate(alnread.alns, state)
Base.getindex(alnread::AlignedRead, i::Int64) = alnread.alns[i]

function Base.show(alnread::AlignedRead)
    println("Alignment with $(length(alnread.alns)) part(s):")
    for part in alnread.alns
        print("   $(part.ref.strand): [$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) - ($(part.seq.first), $(part.seq.last)) on read - ")
        isempty(annotation(part)) ? println("not annotated.") : println("annotation: $(annotationtype(part)):$(annotationname(part)) [$(annotationoverlap(part))]")
    end
end

function hasannotation(alnread::AlignedRead, annotation_name::String)
    for part in alnread
        annotationname(part) == annotation_name && (return true)
    end
    return false
end

Base.length(i::Interval) = rightposition(i) - leftposition(i) + 1
function overlapdistance(i1::Interval, i2::Interval)::Float64
    strand(i1) != strand(i2) && return -Inf
    return min(min(length(i1), length(i2)), min(rightposition(i1) - leftposition(i2), rightposition(i2) - leftposition(i1)))
end
distance(i1::Interval, i2::Interval)::Float64 = -min(-0.0, overlapdistance(i1,i2))

function countconcordant(alnread1::AlignedRead, alnread2::AlignedRead; max_distance=100)
    c = 0
    for part in alnread1, otherpart in alnread2
        if hasannotation(part) && hasannotation(otherpart)
            annotationname(part) == annotationname(otherpart) && (c+=1)
        else
            distance(refinterval(part), refinterval(otherpart)) < max_distance && (c+=1)
        end
    end
    return c
end

function ischimeric(alnread::AlignedRead)
    return count(alnread) > 1 ? true : false
end

function ischimeric(alnread1::AlignedRead, alnread2::AlignedRead; max_distance=100)
    (ischimeric(alnread1) || ischimeric(alnread2)) && (return true)
    return (count(alnread1) + count(alnread2) - countconcordant(alnread1, alnread2; max_distance=max_distance)) >= 2
end

function istriplet(alnread1::AlignedRead, alnread2::AlignedRead; max_distance=100)
    2 < (count(alnread1) + count(alnread2)) < 5 || (return false)
    return (count(alnread1) + count(alnread2) - countconcordant(alnread1, alnread2; max_distance=max_distance)) == 3
end

struct Alignments <: AlignmentContainer
    dict::Dict{UInt, AlignedRead}
    name::Union{String,Nothing}
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

function Alignments(bam_file::String; stop_at=nothing, name=nothing)
    alignments1, alignments2 = read_bam(bam_file; stop_at=stop_at)
    @assert isempty(alignments2)
    Alignments(alignments1, name)
end

struct PairedAlignments <: AlignmentContainer
    dict::Dict{UInt, Tuple{AlignedRead, AlignedRead}}
    name::Union{String, Nothing}
end

Base.length(alignments::PairedAlignments) = length(alignments.dict)
Base.keys(alignments::PairedAlignments) = keys(alignments.dict)
Base.values(alignments::PairedAlignments) = values(alignments.dict)
function Base.iterate(alignments::PairedAlignments) 
    dictiteration = iterate(alignments.dict)
    isnothing(dictiteration) && (return nothing)
    ((key, (aln1, aln2)), state) = dictiteration
    return ((aln1, aln2), state)
end
function Base.iterate(alignments::PairedAlignments, state::Int) 
    dictiteration = iterate(alignments.dict, state)
    isnothing(dictiteration) && (return nothing)
    ((key, (aln1, aln2)), state) = dictiteration
    return ((aln1, aln2), state)
end

function PairedAlignments(bam_file1::String, bam_file2::String; stop_at=nothing, name=nothing)
    alignments1, alignments_e1 = read_bam(bam_file1; stop_at=stop_at)
    alignments2, alignments_e2 = read_bam(bam_file2; stop_at=stop_at)
    @assert isempty(alignments_e1) && isempty(alignments_e2)
    alignments = Dict(key=>(alignments1[key], alignments2[key]) for key in intersect(Set(keys(alignments1)), Set(keys(alignments2))))
    PairedAlignments(alignments, name)
end

function PairedAlignments(pebam_file::String; stop_at=nothing, name=nothing)
    alignments1, alignments2 = read_bam(pebam_file; stop_at=stop_at)
    alignments = Dict(key=>(alignments1[key], alignments2[key]) for key in intersect(Set(keys(alignments1)), Set(keys(alignments2))))
    PairedAlignments(alignments, name)
end

function ispaired(record::BAM.Record)::Bool
    return BAM.isfilled(record) && (BAM.flag(record) & UInt16(0x001) != 0)
end

function isread1(record::BAM.Record)::Bool
    return BAM.isfilled(record) && (BAM.flag(record) & UInt16(0x040) != 0)
end

function isread2(record::BAM.Record)::Bool
    return BAM.isfilled(record) && (BAM.flag(record) & UInt16(0x080) != 0)
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

function read_bam(bam_file::String; stop_at=nothing)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    reads1 = Dict{UInt, AlignedRead}()
    reads2 = Dict{UInt, AlignedRead}()
    is_bitstring = is_bitstring_bam(bam_file)
    c = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        id = is_bitstring ? parse(UInt, BAM.tempname(record); base=2) : hash(record.data[1:BAM.seqname_length(record)])
        current_read_dict = isread2(record) && ispaired(record) ? reads2 : reads1
        id in keys(current_read_dict) ? merge!(current_read_dict[id], AlignedRead(record)) : push!(current_read_dict, id=>AlignedRead(record))
        c += 1
        isnothing(stop_at) || ((c >= stop_at) && break) 
    end
    close(reader)
    return reads1, reads2
end

strand_filter(a::Interval, b::Interval) = strand(a) == strand(b)

function annotate!(alns::Alignments, features::Features)
    for alignedread in alns
        for alignment in alignedread
            for feature_interval in eachoverlap(features.list, refinterval(alignment), filter=strand_filter)
                olp = round(UInt8, (overlapdistance(feature_interval, refinterval(alignment)) / length(refinterval(alignment))) * 100)
                if annotationoverlap(alignment) < olp
                    alignment.ref.metadata.name = feature_interval.metadata.name
                    alignment.ref.metadata.type = feature_interval.metadata.type
                    alignment.ref.metadata.overlap = olp
                end
            end
        end
    end
end

function annotate!(alns::PairedAlignments, features::Features)
    for (alignment1, alignment2) in alns
        for alignedread in [alignment1, alignment2]
            for alignment in alignedread
                for feature_interval in eachoverlap(features.list, refinterval(alignment), filter=strand_filter)
                    olp = round(UInt8, (overlapdistance(feature_interval, refinterval(alignment)) / length(refinterval(alignment))) * 100)
                    if annotationoverlap(alignment) < olp
                        alignment.ref.metadata.name = feature_interval.metadata.name
                        alignment.ref.metadata.type = feature_interval.metadata.type
                        alignment.ref.metadata.overlap = olp
                    end
                end
            end
        end
    end
end