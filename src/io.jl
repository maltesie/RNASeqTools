function write_file(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

struct Annotation
    type::String
    name::String
end

Base.isempty(annotation::Annotation) = isempty(annotation.type) && isempty(annotation.name)
Base.isnothing(annotation::Annotation) = isnothing(annotation.type) && isnothing(annotation.name)
name(interval::Interval{Annotation}) = interval.metadata.name
type(interval::Interval{Annotation}) = interval.metadata.type

struct Features <: AnnotationContainer
    list::IntervalCollection
    description::Union{String, Nothing}
end

function Features(gff_file::String, type::Vector{String}, name_key::String; description=nothing)
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    for feature in features
        GFF3.featuretype(feature) in type || continue
        seqname = GFF3.seqid(feature)
        annotation = Annotation(GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","))
        push!(intervals, Interval(seqname, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annotation))
    end
    return Features(IntervalCollection(intervals, true), description)
end

function Features(gff_file::String, type::String, name_key::String; description=nothing)
    return Features(gff_file, [type], name_key; description=description)
end

Base.push!(features::Features, interval::Interval) = push!(features.list, interval)
function append!(features1::Features, features2::Features) 
    for feature in features2
        push!(features1, feature)
    end
end

Base.iterate(features::Features) = iterate(features.list)
Base.iterate(features::Features, state) = iterate(features.list, state)
Base.length(features::Features) = length(features.list)

struct AlignmentPart
    ref::Interval{Vector{Annotation}}
    seq::Interval{Nothing}
    isprimary::Bool
end

function AlignmentPart(xa_part::String)
    chr, pos, cigar, nm = split(xa_part, ",")
    refstart = parse(Int, pos)
    strand = refstart > 0 ? Strand('+') : Strand('-')
    refstart *= sign(refstart)
    readstart, readstop, relrefstop = positions(cigar)
    ref_interval = Interval(chr, refstart, refstart+relrefstop, strand, Vector{Annotation}())
    seq_interval = Interval(readstart, readstop)
    return AlignmentPart(ref_interval, seq_interval, false)
end

struct ReadAlignment
    parts::Vector{AlignmentPart}
    alt::Vector{Vector{AlignmentPart}}
end

function ReadAlignment(record::BAM.Record; invertstrand=false)
    BAM.ismapped(record) || return ReadAlignment([], [])
    readstart, readstop, relrefstop = positions(BAM.cigar(record))
    strand = BAM.ispositivestrand(record) ? Strand('+') : Strand('-')
    ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), strand, Vector{Annotation}())
    seq_interval = Interval("read", readstart, readstop)
    aln_part = AlignmentPart(ref_interval, seq_interval, BAM.isprimary(record))
    alts = hasxatag(record) ? [AlignmentPart(xa) for xa in split(String(xatag(record)))[1:end-1]] : AlignmentPart[]
    return ReadAlignment([aln_part], alts)
end 

Base.iterate(readalignment::ReadAlignment) = iterate(readalignment.parts)
Base.iterate(readalignment::ReadAlignment, state::Int) = iterate(readalignment.parts, state)
function Base.show(readaln::ReadAlignment)
    println("Alignment with $(length(readaln.parts)) part(s):")
    for part in readaln.parts
        print("   $(part.ref.strand): [$(part.ref.first), $(part.ref.last)] on $(part.ref.seqname) - ($(part.seq.first), $(part.seq.last)) on read - ")
        annotation_string = join([anno.type*":"*anno.name for anno in part.ref.metadata], ";")
        isempty(annotation_string) ? println("not annotated.") : println("annotations: $annotation_string")
    end
end

function Base.push!(readalignment::ReadAlignment, record::BAM.Record)
    BAM.ismapped(record) || return nothing
    readstart, readstop, relrefstop = positions(BAM.cigar(record))
    strand = BAM.ispositivestrand(record) ? Strand('+') : Strand('-')
    ref_interval = Interval(BAM.refname(record), BAM.leftposition(record), BAM.rightposition(record), strand, Vector{Annotation}())
    seq_interval = Interval("read", readstart, readstop)
    aln_part = AlignmentPart(ref_interval, seq_interval, BAM.isprimary(record))
    push!(readalignment.parts, aln_part)
    hasxatag(record) && push!(readalignment.alt, [AlignmentPart(xa) for xa in split(String(xatag(record)))[1:end-1]])
end

function primaryalignmentpart(readalignment::ReadAlignment)
    for aln_part in readalignment
        aln_part.isprimary && (return aln_part)
    end
end

function hasannotation(readaln::ReadAlignment, annotation_name::String)
    for part in readaln
        for annot in part.ref.metadata
            annot.name == annotation_name && (return true)
        end
    end
    return false
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
        id in keys(current_read_dict) ? push!(current_read_dict[id], record) : push!(current_read_dict, id=>ReadAlignment(record))
        c += 1
        isnothing(stop_at) || ((c >= stop_at) && break) 
    end
    close(reader)
    return reads1, reads2
end

struct Genome <: SequenceContainer
    seq::LongDNASeq
    chrs::Dict{String, UnitRange{Int}}
    name::Union{String, Nothing}
end

function Genome(sequences::Vector{LongDNASeq}, names::Vector{String}; description=nothing)
    seq = LongDNASeq()
    chrs = Dict{String, UnitRange{Int}}()
    sequence_position = 1
    for (sequence, name) in zip(sequences, names)
        seq *= sequence
        push!(chrs, name=>sequence_position:sequence_position+length(sequence)-1)
        sequence_position += length(sequence)
    end
    return Genome(seq, chrs, description)
end

function Genome(sequence::LongDNASeq, name::String)
    return Genome([sequence], [name])
end

function Genome(genome_fasta::String)
    (name, sequences) = read_genomic_fasta(genome_fasta)
    chrs::Dict{String,UnitRange{Int}} = Dict()
    total_seq = ""
    temp_start = 1
    for (chr, chr_seq) in sequences
        chrs[chr] = temp_start:(temp_start+length(chr_seq)-1)
        temp_start += length(chr_seq)
        total_seq *= chr_seq
    end
    Genome(LongDNASeq(total_seq), chrs, name)
end

Base.length(genome::Genome) = length(genome.seq)

function chomosomecount(genome::Genome)
    return length(genome.chrs)
end

function Base.iterate(genome::Genome)
    (chr, slice) = first(genome.chrs)
    ((chr, genome.seq[slice]), 1)
end

function Base.iterate(genome::Genome, state::Int)
    state += 1
    state > genome.chrs.count && (return nothing)
    for (i, (chr, slice)) in enumerate(genome.chrs)
        (i == state) && (return ((chr, genome.seq[slice]), state))
    end
end

function Base.:*(genome1::Genome, genome2::Genome)
    new_description = isnothing(genome1.name) ? genome2 : (isnothing(genome2.name) ? genome1.name : genome1.name * "+" * genome2.name)
    return Genome(genome1.seq*genome2.seq, merge(genome1.chrs, Dict(key=>(range .+ length(genome1)) for (key, range) in genome2.chrs)), new_description)
end

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>String(seq) for (chr, seq) in genome), file; name=genome.name)
end

function read_genomic_fasta(fasta_file::String)
    genome::Dict{String, String} = Dict()
    chrs = String[]
    start_ids = Int[]
    name = ""
    open(fasta_file, "r") do file
        lines = readlines(file)
        startswith(lines[1], ">") && (name = join(split(lines[1])[2:end]))
        for (i,line) in enumerate(lines)
            startswith(line, ">") &&  (push!(chrs, split(line," ")[1][2:end]); push!(start_ids, i))
        end
        push!(start_ids, length(lines)+1)
        for (chr, (from,to)) in zip(chrs, [@view(start_ids[i:i+1]) for i in 1:length(start_ids)-1])
            genome[chr] = join(lines[from+1:to-1])
        end
    end
    return name, genome
end

function write_genomic_fasta(genome::Dict{String, String}, fasta_file::String; name=nothing, chars_per_row=80)
    open(fasta_file, "w") do file
        for (i, (chr, seq)) in enumerate(genome)
            s = String(seq)
            l = length(s)
            !isnothing(name) ? println(file, ">$chr") : println(file, ">$chr $name")
            for i in 0:length(seq)÷chars_per_row
                ((i+1)*chars_per_row > l) ? println(file, s[i*chars_per_row+1:end]) : println(file, s[i*chars_per_row+1:(i+1)*chars_per_row])
            end
        end
    end
end

struct PairedReads <: SequenceContainer
    dict::Dict{UInt, LongDNASeqPair}
    name::Union{String, Nothing}
end

Base.length(reads::PairedReads) = length(reads.dict)
Base.keys(reads::PairedReads) = keys(reads.dict)
Base.values(reads::PairedReads) = values(reads.dict)
function Base.iterate(reads::PairedReads) 
    it = iterate(reads.dict)
    isnothing(it) ? (return nothing) : ((key, (read1, read2)), state) = it
    return ((read1, read2), state)
end
function Base.iterate(reads::PairedReads, state::Int) 
    it = iterate(reads.dict, state)
    isnothing(it) ? (return nothing) : ((key, (read1, read2)), state) = it
    return ((read1, read2), state)
end

function PairedReads(file1::String, file2::String; description=nothing, stop_at=nothing)
    reads1 = read_reads(file1; nb_reads=stop_at)
    reads2 = read_reads(file2; nb_reads=stop_at)
    @assert length(reads1) == length(reads2)
    @assert all([haskey(reads2, key) for key in keys(reads1)])
    PairedReads(Dict(key=>(reads1[key], reads2[key]) for key in keys(reads1)), description)
end

function Base.write(fasta_file1::String, fasta_file2::String, reads::PairedReads)
    f1 = endswith(fasta_file1, ".gz") ? GzipCompressorStream(open(fasta_file1, "w")) : open(fasta_file1, "w")
    f2 = endswith(fasta_file2, ".gz") ? GzipCompressorStream(open(fasta_file2, "w")) : open(fasta_file2, "w")
    for (key, (read1, read2)) in reads.dict
        str_key = bitstring(key)
        write(f1, ">$str_key\n$(String(read1))\n")
        write(f2, ">$str_key\n$(String(read2))\n")
    end
    close(f1)
    close(f2)
end

struct Reads <: SequenceContainer
    dict::Dict{UInt, LongDNASeq}
    name::Union{String, Nothing}
end

Base.length(reads::Reads) = length(reads.dict)
Base.keys(reads::Reads) = keys(reads.dict)
Base.values(reads::Reads) = values(reads.dict)
function Base.iterate(reads::Reads) 
    it = iterate(reads.dict)
    isnothing(it) ? (return nothing) : ((key, read), state) = it
    return (read, state)
end
function Base.iterate(reads::Reads, state::Int) 
    it = iterate(reads.dict, state)
    isnothing(it) ? (return nothing) : ((key, read), state) = it
    return (read, state)
end

function Base.write(fasta_file::String, reads::Reads)
    f = endswith(fasta_file, ".gz") ? GzipCompressorStream(open(fasta_file, "w")) : open(fasta_file, "w")
    for (key, read) in reads.dict
        write(f, ">$(bitstring(key))\n$(String(read))\n")
    end
    close(f)
end

function Reads(file::String; description="", stop_at=nothing)
    reads = read_reads(file, nb_reads=stop_at)
    Reads(reads, description)
end

function Reads(f, paired_reads::PairedReads; use_when_tied=:none)
    @assert use_when_tied in [:none, :read1, :read2]
    reads = Dict{UInt, LongDNASeq}()
    for (key, (read1, read2)) in paired_reads
        if use_when_tied == :read1 
            f(read1) ? push!(reads, key=>copy(read1)) : (f(read2) && push!(reads, key=>copy(read2)))
        elseif use_when_tied == :read2
            f(read2) ? push!(reads, key=>copy(read2)) : (f(read1) && push!(reads, key=>copy(read1)))
        elseif use_when_tied == :none
            check1, check2 = f(read1), f(read2)
            check1 && check2 && continue
            check1 && push!(reads, key=>copy(read1))
            check2 && push!(reads, key=>copy(read2))
        end
    end
    Reads(reads, paired_reads.name)
end

function read_reads(file::String; nb_reads=nothing)
    @assert any([endswith(file, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])
    reads::Dict{UInt, LongDNASeq} = Dict()
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    is_bitstring = is_bitstring_fasta(file)
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()
    sequencer = is_fastq ? FASTQ.sequence : FASTA.sequence
    read_counter = 0
    while !eof(reader)
        read!(reader, record)
        id = is_bitstring ? parse(UInt, identifier(record); base=2) : hash(record.data[record.identifier])
        push!(reads, id => LongDNASeq(record.data[record.sequence]))
        read_counter += 1
        isnothing(nb_reads) || (read_counter >= nb_reads && break)
    end
    close(reader)
    return reads
end

struct SingleTypeFiles <: FileCollection
    list::Vector{String}
    type::String
end

function SingleTypeFiles(files::Vector{String})
    endings = [fname[findlast('.', fname):end] for fname in files]
    @assert length(unique(endings)) == 1
    SingleTypeFiles(files, endings[1])
end

function SingleTypeFiles(folder::String, type::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if endswith(fname, type)], type)
end

function SingleTypeFiles(folder::String, type::String, prefix::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if (endswith(fname, type) && startswith(fname, prefix))], type)
end

Base.length(files::SingleTypeFiles) = length(files.list)

function Base.iterate(files::SingleTypeFiles)
    isempty(files.list) && (return nothing)
    return (files.list[1], 1)
end

function Base.iterate(files::SingleTypeFiles, state::Int)
    state + 1 > length(files.list) && (return nothing)
    return (files.list[state+1], state + 1)
end

function hassingledir(files::SingleTypeFiles)
    return length(unique(dirname(file) for file in files)) == 1
end

function Base.dirname(files::SingleTypeFiles)
    @assert hassingledir(files)
    return dirname(files.list[1])
end

struct PairedSingleTypeFiles <: FileCollection
    list::Vector{Tuple{String,String}}
    type::String
    suffix1::Union{String,Nothing}
    suffix2::Union{String,Nothing}
end

function PairedSingleTypeFiles(files1::Vector{String}, files2::Vector{String})
    endingsa = [fname[findlast(fname, "."):end] for fname in files1]
    endingsb = [fname[findlast(fname, "."):end] for fname in files2]
    @assert (length(unique(endingsa)) == 1) && (unique(endingsa) == unique(endingsb))
    PairedSingleTypeFiles(collect(zip(files1, files2)), endingsa[1], nothing, nothing)
end

function PairedSingleTypeFiles(folder::String, type::String; suffix1="_1", suffix2="_2", prefix=nothing)
    type_files = [joinpath(folder, fname) for fname in readdir(folder) if isnothing(prefix) ? endswith(fname, type) : endswith(fname, type) && startswith(fname, prefix)]
    names1 = [f[1:end-(length(type)+length(suffix1))] for f in type_files if f[end-(length(type)+length(suffix1)-1):end-length(type)] == suffix1]
    names2 = [f[1:end-(length(type)+length(suffix2))] for f in type_files if f[end-(length(type)+length(suffix2)-1):end-length(type)] == suffix2]
    @assert Set(names1) == Set(names2)
    PairedSingleTypeFiles([(joinpath(folder, name * suffix1 * type), joinpath(folder, name * suffix2 * type)) for name in names1], type, suffix1, suffix2)
end

Base.length(files::PairedSingleTypeFiles) = length(files.list)

function Base.iterate(files::PairedSingleTypeFiles)
    isempty(files.list) && (return nothing)
    return (files.list[1], 1)
end

function Base.iterate(files::PairedSingleTypeFiles, state::Int)
    state + 1 > length(files.list) && (return nothing)
    return (files.list[state+1], state + 1)
end

function hassingledir(files::PairedSingleTypeFiles)
    dirs1 = unique([dirname(file[1]) for file in files])
    dirs2 = unique([dirname(file[2]) for file in files])
    return length(dirs1) == 1 && Set(dirs1) == Set(dirs2)
end

function Base.dirname(files::PairedSingleTypeFiles)
    @assert hassingledir(files)
    return dirname(files.list[1][1])
end