function write_file(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

struct Alignment <: Annotation
    start::Int
    stop::Int
    data::Vector{UInt8}
end

struct Alignments <: SequenceContainer
    dict::Dict{UInt, Alignment}
    chrs::Dict{Int, String}
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

function read_bam(bam_file::String; uniques_only=false, stop_at=nothing)
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    chrs = Dict(i=>name for (i,name) in enumerate(bam_chromosome_names(reader)))
    aligned_reads = Dict()
    c = 0
    while !eof(reader)
        read!(reader, record)
        !BAM.ismapped(record) && continue
        start, stop = BAM.position(record)*strandint(record), BAM.rightposition(record)*strandint(record)
        start > stop && ((start, stop) = (stop, start))
        slice = BAM.auxdata_position(record):BAM.data_size(record)
        xa = get_XA_tag(@view(record.data[slice]))
        (uniques_only && !isnothing(xa)) && continue
        nms = get_NM_tag(@view(record.data[slice]))
        new_row = DataFrame(name=BAM.tempname(record), start=start, stop=stop, chr=BAM.refname(record), 
                            nm=nms, xa=xa_string, cigar=BAM.cigar(record))
        append!(aligned_reads, new_row)
        ((nb_reads > 0) & (c >= nb_reads)) && break 
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
    spec::String
end

function Genome(genome_fasta::String)
    (name, sequences) = read_genomic_fasta(genome_fasta)
    chrs::Dict{String,UnitRange{Int}} = Dict()
    seq = ""
    temp_start = 1
    for (chr, sequence) in sequences
        chrs[chr] = temp_start:(temp_start+length(seq)-1)
        temp_start += length(seq)
        seq *= sequence
    end
    return Genome(LongDNASeq(seq), chrs, name)
end

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>String(genome.seq[s]) for (chr, s) in genome.chrs), file; name=genome.spec)
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
    name::String
    count::Int
end

function PairedReads(file1::String, file2::String; description="", stop_at=nothing)
    reads1 = read_reads(file1; nb_reads=stop_at)
    reads2 = read_reads(file2; nb_reads=stop_at)
    @assert length(reads1) == length(reads2)
    @assert all([haskey(reads2, key) for key in keys(reads1)])
    PairedReads(Dict(key=>(reads1[key], reads2[key]) for key in keys(reads1)), description, length(reads1))
end

function Base.write(file1::String, file2::String, reads::PairedReads)
    writer1 = FASTA.Writer(GzipCompressorStream(open(file1, "w")))
    writer2 = FASTA.Writer(GzipCompressorStream(open(file2, "w")))
    for (key, (read1, read2)) in reads.dict
        record1 = FASTA.Record(key, read1)
        record2 = FASTA.Record(key, read2)
        write(writer1, record1)
        write(writer2, record2)
    end
    close(writer1)
    close(writer2)
end

struct Reads <: SequenceContainer
    dict::Dict{UInt64, LongDNASeq}
    name::String
    count::Int
end

function Base.write(file::String, reads::Reads)
    writer = FASTA.Writer(GzipCompressorStream(open(file, "w")))
    for (key, (read1, read2)) in reads.dict
        record = FASTA.Record(key, read1)
        write(writer, record)
    end
    close(writer)
end

function Reads(file::String; description="", stop_at=nothing)
    reads = read_reads(file, nb_reads=stop_at)
    Reads(reads, description, length(reads))
end

function read_reads(file::String; nb_reads=nothing)
    @assert any([endswith(file, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])
    reads::Dict{UInt64, LongDNASeq} = Dict()
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    (is_zipped && is_fastq) && (reader = FASTQ.Reader(GzipDecompressorStream(open(file, "r"))))
    (is_zipped && !is_fastq) && (reader = FASTA.Reader(GzipDecompressorStream(open(file, "r"))))
    (!is_zipped && is_fastq) && (reader = FASTQ.Reader(open(file, "r")))
    (!is_zipped && !is_fastq) && (reader = FASTA.Reader(open(file, "r"))) 
    is_fastq ? record = FASTQ.Record() : record = FASTA.Record()
    my_seq = LongDNASeq()
    c = 0
    while !eof(reader)
        c += 1
        read!(reader, record)
        my_seq = LongDNASeq(record.data[record.sequence])
        push!(reads, hash(record.data[record.identifier])=>my_seq)
        isnothing(nb_reads) || ((c >= nb_reads) && break)
    end
    close(reader)
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

