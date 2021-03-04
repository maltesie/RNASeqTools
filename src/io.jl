function write_file(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

struct Alignments <: AlignmentContainer
    dict::Dict{UInt, BAM.Record}
end

function Alignments(bam_file::String; uniques_only=false, stop_at=nothing, use_if_pe=:read1)
    @assert use_if_pe in [:read1, :read2]
    alignments = read_bam(bam_file; uniques_only=uniques_only, stop_at=stop_at, use_if_pe=use_if_pe)
    Alignments(alignments, length(alignments))
end

struct PairedAlignments <: AlignmentContainer
    dict::Dict{UInt, Tuple{BAM.Record, BAM.Record}}
    count::Int
end

function PairedAlignments(bam_file1::String, bam_file2::String; uniques_only=false, stop_at=nothing)
    alignments1 = read_bam(bam_file1; uniques_only=uniques_only, stop_at=stop_at)
    alignments2 = read_bam(bam_file2; uniques_only=uniques_only, stop_at=stop_at)
    alignments = Dict(key=>(alignments1[key], alignments2[key]) for key in intersect(Set(keys(alignments1)), Set(keys(alignments2))))
    PairedAlignments(alignments, length(alignments))
end

function PairedAlignments(pebam_file::String; uniques_only=false, stop_at=nothing)
    alignments = read_bam(pebam_file; uniques_only=uniques_only, stop_at=stop_at, use_if_pe=:both)
    PairedAlignments(alignments, length(alignments))
end

function strandint(record::BAM.Record; is_rev=false)
    BAM.ispositivestrand(record) ? strand = 1 : strand = -1
    is_rev ? (return strand * -1) : (return strand)
end

@inline function translated_data(data::SubArray{UInt8,1})
    for i in 1:length(data)
        (data[i] == 0x00) && (return data[1:i-1])
    end
end

@inline function get_NM_tag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
      (0x4d == data[i]) & (0x43 == data[i+1]) && (return Int(data[i+2]))
    end
    return nothing
end

@inline function get_XA_tag(data::SubArray{UInt8,1})
    for i in 1:length(data)-2
        (0x00 == data[i]) & (UInt8('X') == data[i+1]) & (UInt8('A') == data[i+2]) && 
        (return translated_data(@view(data[i+4:end])))
    end
    return nothing
end

function read_bam(bam_file::String; uniques_only=false, stop_at=nothing, use_if_pe=:both)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    reads1 = Dict{UInt, BAM.Record}()
    reads2 = Dict{UInt, BAM.Record}()
    c = 0
    while !eof(reader)
        read!(reader, record)
        !BAM.ismapped(record) && continue
        start, stop = BAM.position(record)*strandint(record), BAM.rightposition(record)*strandint(record)
        start > stop && ((start, stop) = (stop, start))
        slice = BAM.auxdata_position(record):BAM.data_size(record)
        xa = get_XA_tag(@view(record.data[slice]))
        (uniques_only && !isnothing(xa)) && continue
        xa_string = isnothing(xa) ? nothing : String(xa)
        nms = get_NM_tag(@view(record.data[slice]))
        new_alignment = Alignment(start, stop, BAM.refname(record), nms, xa_string, "")
        push!(aligned_reads, hash(record.data[1:BAM.seqname_length(record)])=>new_alignment)
        c += 1
        isnothing(stop_at) || ((c >= stop_at) && break) 
    end
    close(reader)
    return aligned_reads
end

function read_wig(wig_file::String)
    coverages::Vector{Dict{String,Vector{Float64}}} = []
    temp_collect = Dict()
    track = ""
    span = 0
    chr = ""
    text = ""
    open(wig_file, "r") do file
        lines = readlines(file)
        for line in lines
            if startswith(line, "track")
                track = split(split(line, " ")[2],"=")[2]
                temp_collect[track] = Dict()
            elseif startswith(line, "variableStep")
                span = parse(Int, split(split(line, " ")[3],"=")[2]) - 1
                chr = split(split(line, " ")[2],"=")[2]
                temp_collect[track][chr] = Tuple{Int, Float64}[]
            else
                str_index, str_value = split(line, " ")
                index = parse(Int, str_index)
                value = parse(Float64, str_value)
                for i in index:index+span
                    push!(temp_collect[track][chr], (i, value))
                end
            end
        end
        for (track, collection) in temp_collect
            coverage = Dict()
            for (chr, points) in collection
                coverage[chr] = zeros(Float64, points[end][1])
                for (index, value) in points
                    coverage[chr][index] = value
                end
            end
            push!(coverages, coverage)
        end
    end
    return coverages
end

struct Genome <: SequenceContainer
    seq::LongDNASeq
    chrs::Dict{String, UnitRange{Int}}
    name::String
    file::Union{String, Nothing}
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
    Genome(LongDNASeq(total_seq), chrs, name, genome_fasta)
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

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>String(genome.seq[s]) for (chr, s) in genome.chrs), file; name=genome.name)
end

function Base.write(genome::Genome)
    @assert !isnothing(genome.file)
    write_genomic_fasta(Dict(chr=>String(genome.seq[s]) for (chr, s) in genome.chrs), genome.file; name=genome.name)
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

function write_genomic_fasta(genome::Dict{String, String}, fasta_file::String; name="", chars_per_row=80)
    open(fasta_file, "w") do file
        for (i, (chr, seq)) in enumerate(genome)
            s = String(seq)
            l = length(s)
            !isempty(name) ? println(file, ">$chr chromosome $i") : println(file, ">$chr $name chromosome $i")
            for i in 0:Int(length(seq)/chars_per_row)
                ((i+1)*chars_per_row > l) ? println(s[i*chars_per_row+1:end]) : println(s[i*chars_per_row+1:(i+1)*chars_per_row])
            end
        end
    end
end

struct PairedReads <: SequenceContainer
    dict::Dict{UInt64, LongDNASeqPair}
    name::Union{String, Nothing}
end

Base.length(reads::PairedReads) = Base.length(reads.dict)
Base.keys(reads::PairedReads) = Base.keys(reads.dict)
Base.values(reads::PairedReads) = Base.values(reads.dict)
Base.iterate(reads::PairedReads) = Base.iterate(reads.dict)[2]
Base.iterate(reads::PairedReads, state::Int) = Base.iterate(reads.dict, state)[2]

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
        bs = bitstring(key)
        write(f1, ">$bs\n$(String(read1))\n")
        write(f2, ">$bs\n$(String(read2))\n")
    end
    close(f1)
    close(f2)
end

struct Reads <: SequenceContainer
    dict::Dict{UInt64, LongDNASeq}
    name::Union{String, Nothing}
end

Base.length(reads::Reads) = Base.length(reads.dict)
Base.keys(reads::Reads) = Base.keys(reads.dict)
Base.values(reads::Reads) = Base.values(reads.dict)
Base.iterate(reads::Reads) = Base.iterate(reads.dict)[2]
Base.iterate(reads::Reads, state::Int) = Base.iterate(reads.dict, state)[2]

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
    reads::Dict{UInt64, LongDNASeq} = Dict()
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    is_bitstring = is_bitstring_fasta(file)
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    read_counter = 0
    while !eof(f)
        id = is_bitstring ? parse(UInt, readline(f)[2:end]; base=2) : hash(split(readline(f)[2:end])[1])
        push!(reads, id => LongDNASeq(readline(f)))
        is_fastq && (skiplines(f, 2))
        read_counter += 1
        isnothing(nb_reads) || (read_counter >= nb_reads && break)
    end
    close(f)
    return reads
end

struct SingleTypeFiles <: FileCollection
    list::Vector{String}
    type::String
end

function SingleTypeFiles(files::Vector{String})
    endings = [fname[findlast(fname, "."):end] for fname in files]
    @assert length(unique(endings)) == 1
    SingleTypeFiles(files, endings[1])
end

function SingleTypeFiles(folder::String, type::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if endswith(fname, type)], type)
end

function SingleTypeFiles(folder::String, type::String, prefix::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if (endswith(fname, type) && startswith(fname, prefix))], type)
end

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