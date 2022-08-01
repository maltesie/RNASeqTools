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

function Base.merge(genomes::Vector{Genome})
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

function write_genomic_fasta(genome::Dict{String, String}, fasta_file::String; name=nothing, chars_per_row=60)
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

function Sequences()
    Sequences(LongDNASeq(""), UInt[], UnitRange{Int}[])
end

function Sequences(::Type{T}) where {T<:Union{String, UInt}}
    Sequences(LongDNASeq(""), T[], UnitRange{Int}[])
end

function Sequences(seqs::Vector{LongDNASeq}, seqnames::Vector{UInt})
    (length(seqnames) == length(seqs)) || throw(AssertionError("number of names must match the number of sequences!"))
    (length(seqnames) == length(unique(seqnames))) || throw(AssertionError("names must be unique!"))
    ranges = Vector{UnitRange{Int}}(undef, length(seqs))
    seq = LongDNASeq("")
    resize!(seq, sum(length(s) for s in seqs))
    current_range = 1:0
    for (i,s) in enumerate(seqs)
        current_range = last(current_range)+1:last(current_range)+length(s)
        seq[current_range] = s
        ranges[i] = current_range
    end
    sortindex = sortperm(seqnames)
    return Sequences(seq, seqnames[sortindex], ranges[sortindex])
end
Sequences(seqs::Vector{LongDNASeq}) = Sequences(seqs, Vector{UInt}(1:length(seqs)))

function Sequences(seqs::Vector{LongDNASeq}, seqnames::Vector{String})
    (length(seqnames) == length(seqs)) || throw(AssertionError("number of names must match the number of sequences!"))
    (length(seqnames) == length(unique(seqnames))) || throw(AssertionError("names must be unique!"))
    ranges = Vector{UnitRange{Int}}(undef, length(seqs))
    seq = LongDNASeq("")
    resize!(seq, sum(length(s) for s in seqs))
    current_range = 1:0
    for (i,s) in enumerate(seqs)
        current_range = last(current_range)+1:last(current_range)+length(s)
        seq[current_range] = s
        ranges[i] = current_range
    end
    sortindex = sortperm(hash.(seqnames))
    return Sequences(seq, seqnames[sortindex], ranges[sortindex])
end

function Sequences(file::String; is_reverse_complement=false, hash_id=true)
    read_reads(file; is_reverse_complement=is_reverse_complement, hash_id=hash_id)
end

function Sequences(file1::String, file2::String; is_reverse_complement=false, hash_id=true)
    read_reads(file1, file2; is_reverse_complement=is_reverse_complement, hash_id=hash_id)
end

Sequences(genomes::Vector{Genome}) =
    Sequences([genome.seq for genome in genomes], UInt.(i for i in 1:length(genomes)), [i:i for i in 1:length(genomes)])

Base.getindex(seqs::Sequences, index::Int) = seqs.seq[seqs.ranges[index]]
Base.getindex(seqs::Sequences, range::Union{StepRange{Int, Int}, UnitRange{Int}}) = Sequences(seqs.seq, seqs.seqnames[range], seqs.ranges[range])

function Base.getindex(seqs::Sequences, index::Union{UInt,String})
    r = searchsorted(seqs.seqnames, index)
    if length(r) === 2
        return (seqs.seq[seqs.ranges[first(r)]], seqs.seq[seqs.ranges[last(r)]])
    elseif length(r) === 1
        return seqs.seq[seqs.ranges[first(r)]]
    else
        throw(KeyError)
    end
end

Base.length(seqs::Sequences) = length(seqs.seqnames)
Base.iterate(seqs::Sequences) = (seqs.seq[seqs.ranges[1]], 2)
Base.iterate(seqs::Sequences, state::Int) = state > length(seqs.ranges) ? nothing : (seqs.seq[seqs.ranges[state]], state+1)
eachpair(seqs::Sequences) = partition(seqs, 2)
Base.empty!(seqs::Sequences) = (empty!(seqs.seq); empty!(seqs.seqnames); empty!(seqs.ranges))

function Base.filter!(seqs::Sequences{T}, ids::Set{T}) where {T<:Union{String, UInt}}
    bitindex = (in).(seqs.seqnames, Ref(ids))
    n_seqs = sum(bitindex)
    seqs.ranges[1:n_seqs] = seqs.ranges[bitindex]
    seqs.seqnames[1:n_seqs] = seqs.seqnames[bitindex]
    resize!(seqs.ranges, n_seqs)
    resize!(seqs.seqnames, n_seqs)
end

function read_reads(file::String; is_reverse_complement=false, hash_id=true)::Sequences
    seqs::Sequences{hash_id ? UInt : String} = hash_id ? Sequences(UInt) : Sequences(String)
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()
    i::Int64 = 0
    current_range::UnitRange{Int64} = 1:0
    while !eof(reader)
        read!(reader, record)
        id = hash_id ? hash(@view(record.data[record.identifier])) : FASTX.identifier(record)
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

function read_reads(file1::String, file2::String; is_reverse_complement=false, hash_id=true)::Sequences
    any([endswith(file1, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]]) &&
    any([endswith(file2, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]]) ||
    throw(AssertionError("Accepted filetypes are: .fastq, .fastq.gz, .fasta and .fasta.gz"))
    seqs::Sequences{hash_id ? UInt : String} = hash_id ? Sequences(UInt) : Sequences(String)

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
        id1 = hash_id ? hash(@view record1.data[record1.identifier]) : FASTX.identifier(record1)
        id2 = hash_id ? hash(@view record2.data[record2.identifier]) : FASTX.identifier(record2)
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

function translation_dict(from_sequence::LongDNASeq, to_sequence::LongDNASeq)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
    res = alignment(pairalign(GlobalAlignment(), from_sequence, to_sequence, scoremodel))
    trans_dict::Dict{Int, Union{Nothing, Int}} = Dict()
    from_pos = 1
    to_pos = 1
    for (b1, b2) in collect(res)
        if (b1 == '-')
            to_pos += 1
        elseif (b2 == '-')
            transdict[from_pos] = nothing
            from_pos += 1
        else
            trans_dict[from_pos] = to_pos
            to_pos += 1
            from_pos += 1
        end
    end
    return trans_dict
end

function cut!(read::LongDNASeq, pos::Int; keep=:left, from=:left)
    0 <= pos <= length(read) || resize!(read, 0)
    if (from == :left) && (keep == :left)
        resize!(read, pos)
    elseif (from == :left) && (keep == :right)
        reverse!(resize!(reverse!(read), length(read)-pos))
    elseif (from == :right) && (keep == :left)
        resize!(read, length(read)-pos)
    elseif (from == :right) && (keep == :right)
        reverse!(resize!(reverse!(read), pos))
    end
end

function cut!(read::LongDNASeq, int::Tuple{Int, Int})
    (0 <= first(int) < last(int) <= length(read)) || resize!(read, 0)
    reverse!(resize!(reverse!(resize!(read, last(int))), length(read)-first(int)+1))
end

function approxoccursin(s1::LongDNASeq, s2::LongDNASeq; k=1)
    return approxsearch(s2, s1, k) != 0:-1
end

approxcount(test_sequence::LongDNASeq, genomes::Vector{Genome}; k=1) = sum(approxoccursin(test_sequence, genome.seq; k=k) for genome in genomes)

approxcount(test_sequences::Sequences, genome::Genome; k=1) = sum(approxoccursin(s, genome.seq; k=k) for s in test_sequences)

similarcount(test_sequence::LongDNASeq, seqs::Sequences, min_score::Float64; score_model=AffineGapScoreModel(match=1, mismatch=-4, gap_open=-5, gap_extend=-1)) =
    sum(normalized_alignment_score(test_sequence, seq; score_model=score_model) >= min_score for seq in seqs)

hassimilar(genome::Genome, test_sequence::LongDNASeq, min_score::Float64; score_model=AffineGapScoreModel(match=1, mismatch=-4, gap_open=-5, gap_extend=-1)) =
    normalized_alignment_score(test_sequence, genome.seq; score_model=score_model) >= min_score

function normalized_alignment_score(read1::LongDNASeq, read2::LongDNASeq;
    score_model=AffineGapScoreModel(match=1, mismatch=-4, gap_open=-5, gap_extend=-1), alignment_type=OverlapAlignment())

    (length(read1) > length(read2)) ? ((short_seq, long_seq) = (read2, read1)) : ((short_seq, long_seq) = (read1, read2))
    aln = pairalign(alignment_type, long_seq, short_seq, score_model; score_only=false)
    return max(BioAlignments.score(aln), 0.0)/length(short_seq)
end

function nucleotidecount(seqs::Sequences; normalize=true, align=:left)
    align in (:left, :right) || throw(AssertionError("align has to be :left or :right"))
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

