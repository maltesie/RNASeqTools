struct Genome
    seq::LongDNASeq
    chroms::Dict{String, UnitRange{Int}}
end

function Genome(sequences::Vector{LongDNASeq}, names::Vector{String})
    seq = LongDNASeq(0)
    chrs = Dict{String, UnitRange}()
    sequence_position = 1
    for (sequence, name) in zip(sequences, names)
        seq *= sequence
        push!(chrs, name=>sequence_position:sequence_position+length(sequence)-1)
        sequence_position += length(sequence)
    end
    return Genome(seq, chrs)
end

function Genome(sequence::LongDNASeq, name::String)
    return Genome([sequence], [name])
end

function Genome(genome_fasta::String)
    sequences = read_genomic_fasta(genome_fasta)
    chrs::Dict{String, UnitRange} = Dict()
    total_seq = ""
    temp_start = 1
    for (chr, chr_seq) in sequences
        chrs[chr] = temp_start:(temp_start+length(chr_seq)-1)
        temp_start += length(chr_seq)
        total_seq *= chr_seq
    end
    Genome(LongDNASeq(total_seq), chrs)
end

Base.length(genome::Genome) = length(genome.seq)
Base.getindex(genome::Genome, key::String) = genome.seq[genome.chroms[key]]

function chomosomecount(genome::Genome)
    return length(genome.chroms)
end

function Base.iterate(genome::Genome)
    (chr, slice) = first(genome.chroms)
    ((chr, genome.seq[slice]), 1)
end

function Base.iterate(genome::Genome, state::Int)
    state += 1
    state > genome.chroms.count && (return nothing)
    for (i, (chr, slice)) in enumerate(genome.chroms)
        (i == state) && (return ((chr, genome.seq[slice]), state))
    end
end

function merge(genomes::Genome ...)
    length(genomes) == 1 && return genomes[1]
    merged = genomes[1]
    for genome in genomes[2:end]
        merged *= genome
    end
    return merged
end

function Base.:*(genome1::Genome, genome2::Genome)
    return Genome(genome1.seq*genome2.seq, merge(genome1.chroms, Dict(key=>(range .+ length(genome1)) for (key, range) in genome2.chroms)))
end

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>String(seq) for (chr, seq) in genome), file)
end

function read_genomic_fasta(fasta_file::String)
    genome::Dict{String, String} = Dict()
    chrs = String[]
    start_ids = Int[]
    open(fasta_file, "r") do file
        lines = readlines(file)
        for (i,line) in enumerate(lines)
            startswith(line, ">") &&  (push!(chrs, split(line," ")[1][2:end]); push!(start_ids, i))
        end
        push!(start_ids, length(lines)+1)
        for (chr, (from,to)) in zip(chrs, [@view(start_ids[i:i+1]) for i in 1:length(start_ids)-1])
            genome[chr] = join(lines[from+1:to-1])
        end
    end
    return genome
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

struct Sequences{T}
    seq::LongDNASeq
    seqnames::Vector{T}
    ranges::Vector{UnitRange{Int}}
end

function Sequences()
    Sequences(LongDNASeq(""), UInt64[], UnitRange{Int}[])
end

function Sequences(seqs::Vector{LongDNASeq}; seqnames=Vector{T}[]) where {T <: Union{String, UInt}}
    if isempty(seqnames)
        seqnames = UInt.(1:length(seqs))
    else
        (length(seqnames) == length(seqs) == length(unique(seqnames))) || throw(AssertionError("number of unique names must match the number of sequences!"))
    end
    ranges = Vector{UnitRange{Int}}(undef, length(seqs))
    seq = LongDNASeq("")
    resize!(seq, sum(length(s) for s in seqs))
    current_range = 1:0
    for (i,s) in enumerate(seqs)
        current_range = last(current_range)+1:last(current_range)+length(s)
        seq[current_range] = s
        ranges[i] = current_range
    end
    sortindex = sortperm(T isa String ? hash.(seqnames) : seqnames)
    return Sequences(seq, seqnames[sortindex], ranges[sortindex])
end

function Sequences(file::String; is_reverse_complement=false)
    read_reads(file; is_reverse_complement=is_reverse_complement)
end

function Sequences(file1::String, file2::String; is_reverse_complement=false)
    read_reads(file1, file2; is_reverse_complement=is_reverse_complement)
end

function Base.getindex(seqs::Sequences, index::UInt) 
    r = searchsorted(seqs.seqnames, index)
    if length(r) === 2
        return (seqs.seq[seqs.ranges[first(r)]], seqs.seq[seqs.ranges[last(r)]]) 
    elseif length(r) === 1
        return seqs.seq[seqs.ranges[first(r)]]
    else
        throw(KeyError)
    end
end

function Base.getindex(seqs::Sequences, index::String) 
    seqs[hash(UInt8.(index))]
end
Base.length(seqs::Sequences) = length(seqs.seqnames)
Base.iterate(seqs::Sequences) = (seqs.seq[seqs.ranges[1]], 2)
Base.iterate(seqs::Sequences, state::Int) = state === length(seqs.ranges) ? nothing : (seqs.seq[seqs.ranges[state]], state+1)
eachpair(seqs::Sequences) = partition(seqs, 2)
Base.empty!(seqs::Sequences) = (empty!(seqs.seq); empty!(seqs.seqnames); empty!(seqs.ranges))

function Base.filter!(seqs::Sequences, ids::Vector{UInt})
    bitindex = (in).(seqs.seqnames, Ref(ids))
    n_seqs = sum(bitindex)
    seqs.ranges[1:n_seqs] = seqs.ranges[bitindex]
    seqs.seqnames[1:n_seqs] = seqs.seqnames[bitindex]
    resize!(seqs.ranges, n_seqs)
    resize!(seqs.seqnames, n_seqs)
end

function Base.write(fasta_file::String, seqs::Sequences)
    f = endswith(fasta_file, ".gz") ? GzipCompressorStream(open(fasta_file, "w")) : open(fasta_file, "w")
    for (i, read) in enumerate(seqs)
        write(f, ">$(seqs.seqnames[i])\n$(String(read))\n")
    end
    close(f)
end

function Base.write(fasta_file1::String, fasta_file2::String, seqs::Sequences)
    f1 = endswith(fasta_file1, ".gz") ? GzipCompressorStream(open(fasta_file1, "w")) : open(fasta_file1, "w")
    f2 = endswith(fasta_file2, ".gz") ? GzipCompressorStream(open(fasta_file2, "w")) : open(fasta_file2, "w")
    for (i, (read1, read2)) in enumerate(eachpair(seqs))
        h = seqs.seqnames[2*i]
        write(f1, ">$h\n$(String(read1))\n")
        write(f2, ">$h\n$(String(read2))\n")
    end
    close(f1)
    close(f2)
end

function read_reads(file::String; is_reverse_complement=false)::Sequences
    seqs::Sequences = Sequences()
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()
    i::Int64 = 0
    current_range::UnitRange{Int64} = 1:0 
    while !eof(reader)
        read!(reader, record)
        id = hash(@view(record.data[record.identifier]))
        s = LongDNASeq(@view(record.data[record.sequence]))
        is_reverse_complement && reverse_complement!(s)
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(record.sequence)
        length(seqs.seq) < last(current_range) && resize!(seqs.seq, length(seqs.seq)+max(1000000, last(current_range)-length(seqs.seq)))
        length(seqs.seqnames) < i && (resize!(seqs.seqnames, length(seqs.seqnames)+10000);resize!(seqs.ranges, length(seqs.ranges)+10000))
        seqs.seqnames[i] = id
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s
    end
    resize!(seqs.seq, last(current_range))
    resize!(seqs.seqnames, i)
    resize!(seqs.ranges, i)
    sort_index = sortperm(seqs.seqnames)
    seqs.seqnames[1:end] = seqs.seqnames[sort_index]
    seqs.ranges[1:end] = seqs.ranges[sort_index]
    close(reader)
    return seqs
end

function read_reads(file1::String, file2::String; is_reverse_complement=false)::Sequences
    any([endswith(file1, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]]) && 
    any([endswith(file2, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]]) ||
    throw(AssertionError("Accepted filetypes are: .fastq, .fastq.gz, .fasta and .fasta.gz"))
    seqs::Sequences = Sequences()
    
    is_fastq1 = any([endswith(file1, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped1 = endswith(file1, ".gz")
    f1 = is_zipped1 ? GzipDecompressorStream(open(file1, "r")) : open(file1, "r")
    reader1 = is_fastq1 ? FASTQ.Reader(f1) : FASTA.Reader(f1)
    record1 = is_fastq1 ? FASTQ.Record() : FASTA.Record()
    
    is_fastq2 = any([endswith(file2, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped2 = endswith(file2, ".gz")
    f2 = is_zipped2 ? GzipDecompressorStream(open(file2, "r")) : open(file2, "r")
    reader2 = is_fastq2 ? FASTQ.Reader(f2) : FASTA.Reader(f2)
    record2 = is_fastq2 ? FASTQ.Record() : FASTA.Record()
    
    i::Int64 = 0
    current_range::UnitRange{Int64} = 1:0 
    while !eof(reader1) && !eof(reader2)
        read!(reader1, record1)
        read!(reader2, record2)
        id1 = hash(@view record1.data[record1.identifier])
        id2 = hash(@view record2.data[record2.identifier])
        id1 === id2 || throw(
            AssertionError("tempnames of read1 and read2 do not match: $(String(record1.data[record1.identifier])) != $(String(record2.data[record2.identifier]))")
        )
        s1 = LongDNASeq(@view record1.data[record1.sequence])
        s2 = LongDNASeq(@view record2.data[record2.sequence])
        is_reverse_complement && ((s1,s2) = (s2,s1))
        reverse_complement!(s2)
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(s1)
        paired_end = last(current_range) + length(s2)
        length(seqs.seq) < paired_end && resize!(seqs.seq, length(seqs.seq)+max(1000000, paired_end-length(seqs.seq)))
        length(seqs.seqnames) < i+1 && (resize!(seqs.seqnames, length(seqs.seqnames)+10000);resize!(seqs.ranges, length(seqs.ranges)+10000))
        seqs.seqnames[i] = id1
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s1
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(s2)
        seqs.seqnames[i] = id2
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s2
    end
    resize!(seqs.seq, last(current_range))
    resize!(seqs.seqnames, i)
    resize!(seqs.ranges, i)
    si = sortperm(seqs.seqnames[1:2:end])
    sort_index = zeros(Int, 2*length(si))
    sort_index[1:2:end] = si .* 2 .- 1
    sort_index[2:2:end] = si .* 2
    seqs.seqnames[1:end] = seqs.seqnames[sort_index]
    seqs.ranges[1:end] = seqs.ranges[sort_index]
    close(reader1)
    close(reader2)
    return seqs
end

function cut!(read::LongDNASeq, pos::Int; keep=:left, from=:left)
    0 <= pos <= length(read) || resize!(read, 0)

    if (from == :left) && (keep == :left)
        resize!(read, pos)

    elseif (from == :left) && (keep == :right)
        reverse!(resize!(reverse!(read), length(read)-pos, true))

    elseif (from == :right) && (keep == :left)
        resize!(read, length(read)-pos)

    elseif (from == :right) && (keep == :right)
        reverse!(resize!(reverse!(read), pos, true))
    end
end

function cut!(read::LongDNASeq, int::Tuple{Int, Int})
    (0 <= first(int) < last(int) <= length(read)) || resize!(read, 0)
    reverse!(resize!(reverse!(resize!(read, last(int), true)), length(read)-first(int)+1, true))
end

function approxoccursin(s1::LongDNASeq, s2::LongDNASeq; k=1)
    return approxsearch(s2, s1, k) != 0:-1
end

occurences(test_sequence::LongDNASeq, seqs::Sequences, similarity_cut::Float64; score_model=AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1)) =
    sum(similarity(test_sequence, seq; score_model=score_model) > similarity_cut for seq in seqs)

function similarity(read1::LongDNASeq, read2::LongDNASeq; score_model=nothing)
    isnothing(score_model) && (score_model = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1))
    (length(read1) > length(read2)) ? ((short_seq, long_seq) = (read2, read1)) : ((short_seq, long_seq) = (read1, read2))
    aln = pairalign(OverlapAlignment(), long_seq, short_seq, score_model; score_only=false)
    return max(BioAlignments.score(aln), 0.0)/length(short_seq)
end

function nucleotidecount(seqs::Sequences; normalize=true)
    max_length = maximum([length(read) for read in seqs])
    count = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    for read in seqs
        (align==:left) ?
        (index = 1:length(read)) :
        (index = (max_length - length(read) + 1):max_length)
        for (i, n) in zip(index, read)
            count[n][i] += 1
        end
    end
    if normalize
        for (_, c) in count
            c /= length(seqs)
        end
    end
    return count
end