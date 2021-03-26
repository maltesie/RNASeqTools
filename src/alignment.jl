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

struct AlignmentPart
    ref::Interval{Vector{Annotation}}
    seq::Interval{Nothing}
    isprimary::Bool
end

function AlignmentPart(xapart::Union{String, SubString{String}})
    chr, pos, cigar, nm = split(xapart, ",")
    refstart = parse(Int, pos)
    strand = refstart > 0 ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop, readlen = positions(cigar)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, Vector{Annotation}())
    seq_interval = strand === STRAND_POS ? Interval("read", readstart, readstop) : Interval("read", readlen-readstop+1, readlen-readstart+1)
    return AlignmentPart(ref_interval, seq_interval, false)
end

annotations(alnpart::AlignmentPart) = alnpart.ref.metadata
annotationnames(alnpart::AlignmentPart) = Set(annotation.name for annotation in annotations(alnpart))
refinterval(alnpart::AlignmentPart) = alnpart.ref
readinterval(alnpart::AlignmentPart) = alnpart.seq

hasannotation(alnpart::AlignmentPart) = !isempty(annotations(alnpart))
function hasannotation(alnpart::AlignmentPart, annotation_name::String)
    for annot in annotations(alnpart)
        annot.name == annotation_name && (return true)
    end
    return false
end

function findfirst(alnpart::AlignmentPart, type::String)
    for annotation in annotations(alnpart)
        annotation.type == type && return annotation
    end
    return nothing
end

struct ReadAlignment
    parts::Vector{AlignmentPart}
end

function ReadAlignment(record::BAM.Record; invertstrand=false)
    BAM.ismapped(record) || return ReadAlignment([])
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
        leftest != 0 && (return ReadAlignment([AlignmentPart(xastrings[leftest])]))
    end
    readstart, readstop, relrefstop, readlen = positions(BAM.cigar(record))
    strand = BAM.ispositivestrand(record) ? Strand('+') : Strand('-')
    ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), strand, Vector{Annotation}())
    seq_interval = strand === STRAND_POS ? Interval("read", readstart, readstop) : Interval("read", readlen-readstop+1, readlen-readstart+1)
    aln_part = AlignmentPart(ref_interval, seq_interval, BAM.isprimary(record))
    return ReadAlignment([aln_part])
end

parts(readalignment::ReadAlignment) = readalignment.parts
count(readalignment::ReadAlignment) = length(readalignment.parts)
merge!(readaln1::ReadAlignment, readaln2::ReadAlignment) = append!(readaln1.parts, readaln2.parts)
Base.iterate(readalignment::ReadAlignment) = iterate(readalignment.parts)
Base.iterate(readalignment::ReadAlignment, state::Int) = iterate(readalignment.parts, state)
Base.getindex(readaln::ReadAlignment, i::Int64) = readaln.parts[i]

function Base.show(readaln::ReadAlignment)
    println("Alignment with $(length(readaln.parts)) part(s):")
    for part in readaln.parts
        print("   $(part.ref.strand): [$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) - ($(part.seq.first), $(part.seq.last)) on read - ")
        annotation_string = join([anno.type*":"*anno.name for anno in part.ref.metadata], ";")
        isempty(annotation_string) ? println("not annotated.") : println("annotations: $annotation_string")
    end
end

hasannotation(readaln::ReadAlignment) = !all([isempty(annotations(alnpart)) for alnpart in readaln])
isfullyannotated(readaln::ReadAlignment) = !any([isempty(annotations(alnpart)) for alnpart in readaln])
function hasannotation(readaln::ReadAlignment, annotation_name::String)
    for part in readaln
        hasannotation(part, annotation_name) && (return true)
    end
    return false
end

function overlapdistance(i1::Interval, i2::Interval)::Float64
    strand(i1) != strand(i2) && return -Inf
    return min(rightposition(i1) - leftposition(i2), rightposition(i2) - leftposition(i1))
end
distance(i1::Interval, i2::Interval)::Float64 = -min(0, overlapdistance(i1,i2))

function annotationoverlap(part1::AlignmentPart, part2::AlignmentPart)
    c = 0
    for a1 in annotations(part1), a2 in annotations(part2)
        a1.name == a2.name && (c+=1)
    end
    return c
end

function countconcordant(readaln1::ReadAlignment, readaln2::ReadAlignment; max_distance=100)
    c = 0
    for part in readaln1, otherpart in readaln2
        if hasannotation(part) && hasannotation(otherpart)
            annotationoverlap(part, otherpart) > 0 && (c+=1)
        else
            distance(refinterval(part), refinterval(otherpart)) > max_distance && (c+=1)
        end
    end
    return c
end

function ischimeric(readaln::ReadAlignment)
    return count(readaln) > 1 ? true : false
end

function ischimeric(readaln1::ReadAlignment, readaln2::ReadAlignment; max_distance=100)
    (ischimeric(readaln1) || ischimeric(readaln2)) && (return true)
    return (count(readaln1) + count(readaln2) - countconcordant(readaln1, readaln2; max_distance=max_distance)) >= 2
end

function istriplet(readaln1::ReadAlignment, readaln2::ReadAlignment; max_distance=100)
    2 < (count(readaln1) + count(readaln2)) < 5 || (return false)
    return (count(readaln1) + count(readaln2) - countconcordant(readaln1, readaln2; max_distance=max_distance)) == 3
end

Base.isempty(readalignment::ReadAlignment) = isempty(readalignment.parts)

struct Alignments <: AlignmentContainer
    dict::Dict{UInt, ReadAlignment}
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
    dict::Dict{UInt, Tuple{ReadAlignment, ReadAlignment}}
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

function read_bam(bam_file::String; stop_at=nothing)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    reads1 = Dict{UInt, ReadAlignment}()
    reads2 = Dict{UInt, ReadAlignment}()
    is_bitstring = is_bitstring_bam(bam_file)
    c = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        id = is_bitstring ? parse(UInt, BAM.tempname(record); base=2) : hash(record.data[1:BAM.seqname_length(record)])
        #println(record.data)
        #println(record.data[1:BAM.templength(record)])
        current_read_dict = isread2(record) && ispaired(record) ? reads2 : reads1
        id in keys(current_read_dict) ? merge!(current_read_dict[id], ReadAlignment(record)) : push!(current_read_dict, id=>ReadAlignment(record))
        c += 1
        isnothing(stop_at) || ((c >= stop_at) && break) 
    end
    close(reader)
    return reads1, reads2
end

strand_filter(a::Interval, b::Interval) = strand(a) == strand(b)

function annotate!(alns::Alignments, features::Features)
    for alignment in alns
        for part in alignment
            for feature_interval in eachoverlap(features.list, part.ref, filter=strand_filter)
                push!(part.ref.metadata, feature_interval.metadata)
            end
        end
    end
end

function annotate!(alns::PairedAlignments, features::Features)
    for (alignment1, alignment2) in alns
        for part in alignment1
            for feature in eachoverlap(features.list, part.ref; filter=strand_filter)
                push!(part.ref.metadata, feature.metadata)
            end
        end
        for part in alignment2
            for feature in eachoverlap(features.list, part.ref; filter=strand_filter)
                push!(part.ref.metadata, feature.metadata)
            end
        end
    end
end