function Genome(sequence_dict::Dict{String, LongDNA{4}})
    seq = LongDNA{4}(undef, sum(length(s) for s in values(sequence_dict)))
    chrs = Dict{String, UnitRange}()
    si = 1
    for (name, sequence) in sequence_dict
        srange = si:si+length(sequence)-1
        seq[srange] = sequence
        push!(chrs, name=>srange)
        si = last(srange) + 1
    end
    return Genome(seq, chrs)
end
Genome(sequences::Vector{LongDNA{4}}, names::Vector{String}) = Genome(Dict(n=>s for (n,s) in zip(sequences, names)))
Genome(sequence::LongDNA{4}, name::String) = Genome([sequence], [name])
Genome(genome_file::String) = Genome(read_genomic_fasta(genome_file))

Base.length(genome::Genome) = length(genome.seq)
Base.getindex(genome::Genome, key::String) = genome.seq[genome.chroms[key]]
Base.getindex(genome::Genome, key::Pair{String, UnitRange}) = genome[first(key)][last(key)]

function nchromosome(genome::Genome)
    return length(genome.chroms)
end

function Base.iterate(genome::Genome)
    (chr, slice) = first(genome.chroms)
    ((chr, genome.seq[slice]), 1)
end

function Base.iterate(genome::Genome, state::Int)
    state += 1
    state > length(genome.chroms) && (return nothing)
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
    write_genomic_fasta(Dict(chr=>seq for (chr, seq) in genome), file)
end

function read_genomic_fasta(fasta_file::String)
    genome::Dict{String, LongDNA{4}} = Dict()
    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            genome[FASTA.identifier(record)] = FASTA.sequence(LongDNA{4}, record)
        end
    end
    return genome
end

function write_genomic_fasta(genome::Dict{String, T}, fasta_file::String) where T <: BioSequence
    open(FASTA.Writer, fasta_file) do writer
        for (chr, seq) in genome
            write(writer, FASTA.Record(chr, seq))
        end
    end
end

function Sequences()
    Sequences(LongDNA{4}(""), UInt[], UnitRange{Int}[])
end

function Sequences(::Type{T}) where {T<:Union{String, UInt}}
    Sequences(LongDNA{4}(""), T[], UnitRange{Int}[])
end

function Sequences(seqs::Vector{LongDNA{4}}, seqnames::Vector{UInt})
    (length(seqnames) == length(seqs)) || throw(AssertionError("number of names must match the number of sequences!"))
    (length(seqnames) == length(unique(seqnames))) || throw(AssertionError("names must be unique!"))
    ranges = Vector{UnitRange{Int}}(undef, length(seqs))
    seq = LongDNA{4}("")
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
Sequences(seqs::Vector{LongDNA{4}}) = Sequences(seqs, Vector{UInt}(1:length(seqs)))

function Sequences(seqs::Vector{LongDNA{4}}, seqnames::Vector{String})
    (length(seqnames) == length(seqs)) || throw(AssertionError("number of names must match the number of sequences!"))
    (length(seqnames) == length(unique(seqnames))) || throw(AssertionError("names must be unique!"))
    ranges = Vector{UnitRange{Int}}(undef, length(seqs))
    seq = LongDNA{4}("")
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
Base.getindex(seqs::Sequences, range::Union{StepRange{Int, Int}, UnitRange{Int}}) = Sequences(seqs.seq, seqs.tempnames[range], seqs.ranges[range])

function Base.getindex(seqs::Sequences, index::Union{UInt,String})
    r = searchsorted(seqs.tempnames, index)
    if length(r) === 2
        return (seqs.seq[seqs.ranges[first(r)]], seqs.seq[seqs.ranges[last(r)]])
    elseif length(r) === 1
        return seqs.seq[seqs.ranges[first(r)]]
    else
        throw(KeyError)
    end
end

Base.length(seqs::Sequences) = length(seqs.tempnames)
nread(seqs::Sequences) = length(unique(seqs.tempnames))
Base.iterate(seqs::Sequences) = (seqs.seq[seqs.ranges[1]], 2)
Base.iterate(seqs::Sequences, state::Int) = state > length(seqs.ranges) ? nothing : (seqs.seq[seqs.ranges[state]], state+1)
eachpair(seqs::Sequences) = partition(seqs, 2)
Base.empty!(seqs::Sequences) = (empty!(seqs.seq); empty!(seqs.tempnames); empty!(seqs.ranges))

function Base.filter!(seqs::Sequences{T}, tempnames::Set{T}) where {T<:Union{String, UInt}}
    bitindex = (in).(seqs.tempnames, Ref(tempnames))
    n_seqs = sum(bitindex)
    seqs.ranges[1:n_seqs] = seqs.ranges[bitindex]
    seqs.tempnames[1:n_seqs] = seqs.tempnames[bitindex]
    resize!(seqs.ranges, n_seqs)
    resize!(seqs.tempnames, n_seqs)
end

function read_reads(file::String; is_reverse_complement=false, hash_id=true)::Sequences
    seqs::Sequences{hash_id ? UInt : String} = hash_id ? Sequences(UInt) : Sequences(String)
    is_fastq = any([endswith(file, ending) for ending in FASTQ_TYPES])
    is_zipped = endswith(file, ".gz")
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()
    i::Int64 = 0
    current_range::UnitRange{Int64} = 1:0
    while !eof(reader)
        read!(reader, record)
        id = hash_id ? hash(@view(record.data[record.identifier])) : FASTX.identifier(record)
        s = LongDNA{4}(@view(record.data[record.sequence]))
        is_reverse_complement && reverse_complement!(s)
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(record.sequence)
        length(seqs.seq) < last(current_range) && resize!(seqs.seq, length(seqs.seq)+max(1000000, last(current_range)-length(seqs.seq)))
        length(seqs.tempnames) < i && (resize!(seqs.tempnames, length(seqs.tempnames)+10000);resize!(seqs.ranges, length(seqs.ranges)+10000))
        seqs.tempnames[i] = id
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s
    end
    resize!(seqs.seq, last(current_range))
    resize!(seqs.tempnames, i)
    resize!(seqs.ranges, i)
    sort_index = sortperm(seqs.tempnames)
    seqs.tempnames[1:end] = seqs.tempnames[sort_index]
    seqs.ranges[1:end] = seqs.ranges[sort_index]
    close(reader)
    return seqs
end

function read_reads(file1::String, file2::String; is_reverse_complement=false, hash_id=true)::Sequences
    seqs::Sequences{hash_id ? UInt : String} = hash_id ? Sequences(UInt) : Sequences(String)

    is_fastq1 = any([endswith(file1, ending) for ending in FASTQ_TYPES])
    is_zipped1 = endswith(file1, ".gz")
    f1 = is_zipped1 ? GzipDecompressorStream(open(file1, "r")) : open(file1, "r")
    reader1 = is_fastq1 ? FASTQ.Reader(f1) : FASTA.Reader(f1)
    record1 = is_fastq1 ? FASTQ.Record() : FASTA.Record()

    is_fastq2 = any([endswith(file2, ending) for ending in FASTQ_TYPES])
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
        s1 = LongDNA{4}(@view record1.data[record1.sequence])
        s2 = LongDNA{4}(@view record2.data[record2.sequence])
        is_reverse_complement && ((s1,s2) = (s2,s1))
        reverse_complement!(s2)
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(s1)
        paired_end = last(current_range) + length(s2)
        length(seqs.seq) < paired_end && resize!(seqs.seq, length(seqs.seq)+max(1000000, paired_end-length(seqs.seq)))
        length(seqs.tempnames) < i+1 && (resize!(seqs.tempnames, length(seqs.tempnames)+10000);resize!(seqs.ranges, length(seqs.ranges)+10000))
        seqs.tempnames[i] = id1
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s1
        i += 1
        current_range = last(current_range)+1:last(current_range)+length(s2)
        seqs.tempnames[i] = id2
        seqs.ranges[i] = current_range
        seqs.seq[current_range] = s2
    end
    resize!(seqs.seq, last(current_range))
    resize!(seqs.tempnames, i)
    resize!(seqs.ranges, i)
    si = sortperm(seqs.tempnames[1:2:end])
    sort_index = zeros(Int, 2*length(si))
    sort_index[1:2:end] = si .* 2 .- 1
    sort_index[2:2:end] = si .* 2
    seqs.tempnames[1:end] = seqs.tempnames[sort_index]
    seqs.ranges[1:end] = seqs.ranges[sort_index]
    close(reader1)
    close(reader2)
    return seqs
end

function translation_dict(from_sequence::LongDNA{4}, to_sequence::LongDNA{4})
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

function cut!(read::LongDNA{4}, pos::Int; keep=:left, from=:left)
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

function cut!(read::LongDNA{4}, int::Tuple{Int, Int})
    (0 <= first(int) < last(int) <= length(read)) || resize!(read, 0)
    reverse!(resize!(reverse!(resize!(read, last(int))), length(read)-first(int)+1))
end

struct MatchIterator
    sq::Union{ExactSearchQuery,ApproximateSearchQuery}
    seq::LongDNA{4}
    max_dist::Int
end
MatchIterator(test_sequence::LongDNA{4}, seq::LongDNA{4}; k=1) =
    MatchIterator(k>0 ? ApproximateSearchQuery(test_sequence) : ExactSearchQuery(test_sequence), seq, k)
function Base.iterate(m::MatchIterator, state=1)
    current_match = m.max_dist > 0 ? findnext(m.sq, m.max_dist, m.seq, state) : findnext(m.sq, m.seq, state)
    isnothing(current_match) && return nothing
    current_match, last(current_match)
end
Base.eachmatch(test_sequence::LongDNA{4}, seq::LongDNA{4}; k=1) = MatchIterator(test_sequence, seq; k=k)

struct GenomeMatchIterator
    mi_f::MatchIterator
    mi_r::MatchIterator
    chroms::Dict{String, UnitRange}
end
GenomeMatchIterator(test_sequence::LongDNA{4}, genome::Genome; k=1) =
    GenomeMatchIterator(MatchIterator(test_sequence, genome.seq; k=k), MatchIterator(reverse_complement(test_sequence), genome.seq; k=k), genome.chroms)
function Base.iterate(gmi::GenomeMatchIterator, (state,rev)=(1,false))
    i = iterate(rev ? gmi.mi_r : gmi.mi_f, state)
    isnothing(i) && return rev ? nothing : iterate(gmi, (1,true))
    for (chr, r) in gmi.chroms
        (first(r) <= first(first(i)) <= last(r)) && return Interval(chr, first(i) .- first(r), rev ? STRAND_NEG : STRAND_POS), (last(i), rev)
    end
end
Base.eachmatch(test_sequence::LongDNA{4}, genome::Genome; k=1) = GenomeMatchIterator(test_sequence, genome; k=k)

function Base.collect(m::T) where T<:Union{MatchIterator,GenomeMatchIterator}
    re = Interval{Nothing}[]
    for match in m
        push!(re, match)
    end
    return re
end

function nucleotidedistribution(seqs::Sequences; normalize=true, align=:left)
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
        for c in values(count)
            c ./= length(seqs)
        end
    end
    return count
end

information_content(m::Matrix{Float64}; n=Inf) = m .* (2 .- [-1 * sum(x > 0 ? x * log2(x) : 0 for x in r) for r in eachcol(m)] .- (3/(2*log(2)*n)))'

function logo(seqs::Sequences)
    nc = nucleotidedistribution(seqs)
    m = hcat(nc[DNA_A], nc[DNA_T], nc[DNA_G], nc[DNA_C])
    information_content(m; n=length(seqs))
end