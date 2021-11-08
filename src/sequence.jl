struct Genome
    seq::LongDNASeq
    chrs::Dict{String, UnitRange}
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
    (name, sequences) = read_genomic_fasta(genome_fasta)
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
Base.getindex(genome::Genome, key::String) = genome.seq[genome.chrs[key]]

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
    return Genome(genome1.seq*genome2.seq, merge(genome1.chrs, Dict(key=>(range .+ length(genome1)) for (key, range) in genome2.chrs)))
end

function Base.write(file::String, genome::Genome)
    write_genomic_fasta(Dict(chr=>String(seq) for (chr, seq) in genome), file)
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

function is_bitstring_fasta(file::String)
    (endswith(file, ".fasta") || endswith(file, ".fasta.gz")) || (return false)
    f = endswith(file, ".fasta.gz") ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    first_line = readline(f)
    close(f)
    ((length(first_line) == 65) && all([c in ['0', '1'] for c in first_line[2:end]])) && (return true)
    return false
end

struct PairedSequences{T} <: SequenceContainer
    dict::Dict{T, LongDNASeqPair}
end

PairedSequences(::T) where T  = PairedSequences(Dict{T, LongDNASeqPair}())

Base.getindex(reads::PairedSequences, index::UInt) = reads.dict[index]
Base.getindex(reads::PairedSequences, index::String) = reads.dict[index]
Base.length(reads::PairedSequences) = length(reads.dict)
Base.keys(reads::PairedSequences) = keys(reads.dict)
Base.values(reads::PairedSequences) = values(reads.dict)
Base.iterate(reads::PairedSequences) = iterate(reads.dict)
Base.iterate(reads::PairedSequences, state::Int) = iterate(reads.dict, state)
Base.empty!(seqs::PairedSequences) = empty!(seqs.dict)

function PairedSequences(file1::String, file2::String; stop_at=nothing, is_reverse_complement=false, hash_id=true)
    read_reads(file1, file2; nb_reads=stop_at, is_reverse_complement=is_reverse_complement, hash_id=hash_id)
end

function Base.write(fasta_file1::String, fasta_file2::String, reads::PairedSequences)
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

struct Sequences{T} <: SequenceContainer
    dict::Dict{T, LongDNASeq}
end

Base.getindex(reads::Sequences, index::UInt) = reads.dict[index]
Base.getindex(reads::Sequences, index::String) = reads.dict[index]
Base.length(reads::Sequences) = length(reads.dict)
Base.keys(reads::Sequences) = keys(reads.dict)
Base.values(reads::Sequences) = values(reads.dict)
Base.iterate(reads::Sequences) = iterate(reads.dict)
Base.iterate(reads::Sequences, state::Int) = iterate(reads.dict, state)
Base.empty!(seqs::Sequences) = empty!(seqs.dict)

function Base.write(fasta_file::String, reads::Sequences{T}) where T
    f = endswith(fasta_file, ".gz") ? GzipCompressorStream(open(fasta_file, "w")) : open(fasta_file, "w")
    for (key, read) in reads.dict
        T === String ? write(f, ">$key\n$(String(read))\n") : write(f, ">$(bitstring(key))\n$(String(read))\n")
    end
    close(f)
end

Sequences(seqs::Vector{LongDNASeq}) = Sequences(Dict(i=>seq for (i::UInt,seq) in enumerate(seqs)))
Sequences(::T) where T  = Sequences(Dict{T, LongDNASeq}())

function Sequences(file::String; stop_at=nothing, is_reverse_complement=false, hash_id=true)
    reads = read_reads(file, nb_reads=stop_at, is_reverse_complement=is_reverse_complement, hash_id=hash_id)
    Sequences(reads)
end

function Sequences(f, paired_reads::PairedSequences{T}; use_when_tied=:none) where T
    @assert use_when_tied in [:none, :read1, :read2]
    reads = Dict{T, LongDNASeq}()
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
    Sequences(reads)
end

function read_reads(file::String; nb_reads=nothing, is_reverse_complement=false, hash_id=true)
    @assert any([endswith(file, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])
    reads::Dict{hash_id ? UInt : String, LongDNASeq} = Dict()
    is_fastq = any([endswith(file, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped = endswith(file, ".gz")
    is_bitstring = is_bitstring_fasta(file)
    f = is_zipped ? GzipDecompressorStream(open(file, "r")) : open(file, "r")
    reader = is_fastq ? FASTQ.Reader(f) : FASTA.Reader(f)
    record = is_fastq ? FASTQ.Record() : FASTA.Record()
    read_counter = 0
    while !eof(reader)
        read!(reader, record)
        id = !hash_id ? identifier(record) : 
                (is_bitstring ? parse(UInt, identifier(record); base=2) : hash(@view(record.data[record.identifier])))
        push!(reads, id => is_reverse_complement ? reverse_complement(LongDNASeq(record.data[record.sequence])) : LongDNASeq(record.data[record.sequence]))
        read_counter += 1
        isnothing(nb_reads) || (read_counter >= nb_reads && break)
    end
    close(reader)
    return reads
end

function read_reads(file1::String, file2::String; nb_reads=nothing, is_reverse_complement=false, hash_id=true)
    @assert any([endswith(file1, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])
    reads::Dict{hash_id ? UInt : String, LongDNASeqPair} = Dict()
    
    is_fastq1 = any([endswith(file1, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped1 = endswith(file1, ".gz")
    is_bitstring1 = is_bitstring_fasta(file1)
    f1 = is_zipped1 ? GzipDecompressorStream(open(file1, "r")) : open(file1, "r")
    reader1 = is_fastq1 ? FASTQ.Reader(f1) : FASTA.Reader(f1)
    record1 = is_fastq1 ? FASTQ.Record() : FASTA.Record()
    
    @assert any([endswith(file2, ending) for ending in [".fastq", ".fastq.gz", ".fasta", ".fasta.gz"]])
    is_fastq2 = any([endswith(file2, ending) for ending in [".fastq", ".fastq.gz"]])
    is_zipped2 = endswith(file2, ".gz")
    is_bitstring2 = is_bitstring_fasta(file2)
    f2 = is_zipped2 ? GzipDecompressorStream(open(file2, "r")) : open(file2, "r")
    reader2 = is_fastq2 ? FASTQ.Reader(f2) : FASTA.Reader(f2)
    record2 = is_fastq2 ? FASTQ.Record() : FASTA.Record()
    
    read_counter = 1
    while !eof(reader1)
        read!(reader1, record1)
        id1 = !hash_id ? identifier(record1) : 
                (is_bitstring1 ? parse(UInt, identifier(record1); base=2) : hash(@view(record1.data[record1.identifier])))
        read!(reader2, record2)
        id2 = !hash_id ? identifier(record2) : 
                (is_bitstring2 ? parse(UInt, identifier(record2); base=2) : hash(@view(record2.data[record2.identifier])))
        id1 == id2 || throw(AssertionError("entry identifiers do not match for entry $read_counter."))
        is_reverse_complement ?
        push!(reads, id1 => (LongDNASeq(record2.data[record2.sequence]), reverse_complement(LongDNASeq(record1.data[record1.sequence])))) : 
        push!(reads, id1 => (LongDNASeq(record1.data[record1.sequence]), reverse_complement(LongDNASeq(record2.data[record2.sequence]))))
        read_counter += 1
        isnothing(nb_reads) || (read_counter > nb_reads && break)
    end
    close(reader1)
    close(reader2)
    return reads
end

function rev_comp!(reads::Sequences)
    for read in reads
        BioSequences.reverse_complement!(read)
    end
end

function rev_comp!(reads::PairedSequences; treat=:both)
    @assert treat in [:both, :read1, :read2]
    for (read1, read2) in reads
        treat in [:both, :read1] && BioSequences.reverse_complement!(read1)
        treat in [:both, :read2] && BioSequences.reverse_complement!(read2)
    end
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

function cut!(reads::Sequences, pos::Int; keep=:left, from=:left)
    for read in reads
        cut!(read, pos; keep=keep, from=from)
    end
end

function cut!(reads::PairedSequences, pos::Int; keep=:left, from=:left)
    for (read1, read2) in reads
        ((pos > length(read1)) && (pos > length(read2))) && continue
        cut!(read1, pos; keep=keep, from=from)
        cut!(read2, pos; keep=keep, from=from)
    end
end

function cut!(reads::Sequences, seq::LongDNASeq; keep=:left_of_query, from=:left)
    @assert keep in [:left_of_query, :right_of_query, :left_and_query, :right_and_query]
    @assert from in [:left, :right]
    for read in reads
        s = from == :left ? findfirst(seq, read) : findlast(seq, read)
        isnothing(s) && continue
        (start, stop) = s
        if keep == :right_of_query
            cut!(read, stop; keep=:right)
        elseif keep == :left_of_query
            cut!(read, start-1; keep=:left)
        elseif keep == :right_and_query
            cut!(read, start-1; keep=:right)
        elseif keep == :left_and_query
            cut!(read, stop; keep=:left)
        end
    end
end

function cut!(reads::PairedSequences, seq::LongDNASeq; keep=:left_of_query, from=:left, treat=:both)
    @assert keep in [:left_of_query, :right_of_query, :left_and_query, :right_and_query]
    @assert treat in [:read1, :read2, :both]
    @assert from in [:left, :right]
    for (read1, read2) in reads
        s1 = treat in [:both, :read1] ? (from == :left ? findfirst(seq, read1) : findlast(seq, read1)) : nothing
        s2 = treat in [:both, :read2] ? (from == :left ? findfirst(seq, read2) : findlast(seq, read2)) : nothing
        if !isnothing(s1)
            start1, stop1 = s1
            if keep == :right_of_query
                cut!(read1, stop1; keep=:right)
            elseif keep == :left_of_query
                cut!(read1, start1-1; keep=:left)
            elseif keep == :right_and_query
                cut!(read1, start1-1; keep=:right)
            elseif keep == :left_and_query
                cut!(read1, stop1; keep=:left)
            end
        end
        if !isnothing(s2)
            start2, stop2 = s2
            if keep == :right_of_query
                cut!(read2, stop2; keep=:right)
            elseif keep == :left_of_query
                cut!(read2, start2-1; keep=:left)
            elseif keep == :right_and_query
                cut!(read2, start2-1; keep=:right)
            elseif keep == :left_and_query
                cut!(read2, stop2; keep=:left)
            end
        end
    end
end

function approxoccursin(s1::LongDNASeq, s2::LongDNASeq; k=1)
    return approxsearch(s2, s1, k) != 0:-1
end

function Base.filter!(f, reads::Sequences)
    for (key, read) in reads.dict
        f(read) || delete!(reads.dict, key)
    end
end

function Base.filter!(f, reads::PairedSequences; logic=:or)
    @assert logic in [:or, :xor, :and]
    for (key, (read1, read2)) in reads.dict
        if logic == :and
            f(read1) && f(read2) || delete!(reads.dict, key)
        elseif logic == :or
            f(read1) || f(read2) || delete!(reads.dict, key)
        elseif logic == :xor
            check1, check2 = f(read1), f(read2)
            ((check1 && !check2) || (!check1 && check2)) || delete!(reads.dict, key)
        end
    end
end

occurences(test_sequence::LongDNASeq, seqs::Sequences{T}, similarity_cut::Float64; score_model=AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1)) where T =
    sum(similarity(test_sequence, seq; score_model=score_model) > similarity_cut for seq in seqs)

function similarity(read1::LongDNASeq, read2::LongDNASeq; score_model=nothing)
    isnothing(score_model) && (score_model = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1))
    (length(read1) > length(read2)) ? ((short_seq, long_seq) = (read2, read1)) : ((short_seq, long_seq) = (read1, read2))
    aln = pairalign(OverlapAlignment(), long_seq, short_seq, score_model; score_only=false)
    return max(BioAlignments.score(aln), 0.0)/length(short_seq)
end

function similarity(reads::PairedSequences{T}) where T
    similarities = Dict{T, Float64}()
    score_model = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-1, gap_extend=-1)
    for (read1, read2) in reads
        push!(similarities, key=>similarity(read1, read2; score_model=score_model))
    end
    return similarities
end

function nucleotidecount(reads::Sequences; normalize=true)
    max_length = maximum([length(read) for read in reads])
    count = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    for read in reads
        (align==:left) ?
        (index = 1:length(read)) :
        (index = (max_length - length(read) + 1):max_length)
        for (i, n) in zip(index, read)
            count[n][i] += 1
        end
    end
    if normalize
        for (_, c) in count
            c /= length(reads)
        end
    end
    return count
end

function nucleotidecount(reads::PairedSequences; normalize=true)
    max_length = maximum(vcat([[length(read1) length(read2)] for (read1, read2) in reads]...))
    count1 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    count2 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    nb_reads = length(reads)
    for (read1, read2) in reads
        (align==:left) ?
        (index1 = 1:length(read1); index2 = 1:length(read2)) :
        (index1 = (max_length - length(read1) + 1):max_length; index2 = (max_length - length(read2) + 1):max_length)
        for ((i1, n1),(i2, n2)) in zip(zip(index1, read1), zip(index2, read2))
            count1[n1][i1] += 1
            count2[n2][i2] += 1
        end
    end
    if normalize
        for ((key1, c1), (key2, c2)) in zip(count1, count2)
            c1 /= length(reads)
            c2 /= length(reads)
        end
    end
    return count1, count2
end
