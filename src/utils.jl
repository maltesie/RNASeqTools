function reverse_complement!(reads::Reads)
    for (key, read) in reads.dict
        BioSequences.reverse_complement!(read)
    end
end

function reverse_complement!(reads::PairedReads; use_read1=true)
    for (key, (read1, read2)) in reads.dict
        use_read1 ? BioSequences.reverse_complement!(read1) : BioSequences.reverse_complement!(read2)
    end
end

function cut!(read::LongDNASeq, pos::Int; keep_left=true, from_left=true)
    from_left ? my_pos = pos : my_pos = length(read) - pos 
    println(read)
    keep_left ? read = read[1:my_pos-1] : read = read[my_pos+1:end]
    println(read)
end

function cut!(reads::Reads, pos::Int; keep_left=true, from_left=true)
    for (key, read) in reads.dict
        (pos > length(read)) || cut!(read, pos; keep_left=keep_left, from_left=from_left)
    end
end

function cut!(reads::PairedReads, pos::Int; keep_left=true, from_left=true)
    for (key, (read1, read2)) in reads.dict
        (pos > length(read1)) || cut!(read1, pos; keep_left=keep_left, from_left=from_left)
        (pos > length(read2)) || cut!(read2, pos; keep_left=keep_left, from_left=from_left)
        println(read1, " ", read2)
    end
end

function cut!(reads::Reads, seq::LongDNASeq; keep_seq=false, keep_left=true)
    for (key, read) in reads.dict
        slice = approxsearch(read, seq, 1)
        slice == 0:-1 && continue
        ((keep_seq && keep_left) || (!keep_seq && !keep_left)) ? pos = last(slice) : pos = first(slice)
        (pos > length(read)) || cut!(read, pos; keep_left=keep_left)
    end
end

function cut!(reads::PairedReads, seq::LongDNASeq; keep_seq=false, keep_left=true)
    for (key, (read1, read2)) in reads.dict
        slice1 = approxsearch(read1, seq, 1)
        slice2 = approxsearch(read2, seq, 1)
        if slice1 != 0:-1
            ((keep_seq && keep_left) || (!keep_seq && !keep_left)) ? pos = last(slice1) : pos = first(slice1)
            (pos > length(read1))|| cut!(read1, pos; keep_left=keep_left)
        end
        if slice2 != 0:-1
            ((keep_seq && keep_left) || (!keep_seq && !keep_left)) ? pos = last(slice2) : pos = first(slice2)
            (pos > length(read2)) || cut!(read2, pos; keep_left=keep_left)
        end
    end
end

function Base.filter!(f, reads::Reads)
    for (key, read) in reads.dict
        f(read) || delete!(reads.dict, key)
    end
end

function Base.filter!(f, reads::PairedReads; both=false)
    for (key, (read1, read2)) in reads.dict
        both ? (f(read1) && f(read2)) || delete!(reads.dict, key) : (f(read1) || f(read2)) || delete!(reads.dict, key)
    end
end

function join_replicates(coverage_notex::Vector{Dict{String,Vector{Float64}}}, coverage_tex::Vector{Dict{String,Vector{Float64}}}, chr::String)
    max_tex = maximum([length(c[chr]) for c in coverage_tex])
    max_notex = maximum([length(c[chr]) for c in coverage_notex])
    l = max(max_tex, max_notex)
    nb_reps = length(coverage_tex)
    tex = zeros(Float64, (l,nb_reps))
    no_tex = zeros(Float64, (l,nb_reps))
    for i in 1:nb_reps
        tex[1:length(coverage_tex[i][chr]),i] = coverage_tex[i][chr]
        no_tex[1:length(coverage_notex[i][chr]),i] = coverage_notex[i][chr]
    end
    return vec(mean(tex, dims=2)), vec(mean(no_tex, dims=2))
end

function join_replicates(coverage_term::Vector{Dict{String,Vector{Float64}}}, chr::String)
    nb_reps = length(coverage_term)
    l = maximum([length(c[chr]) for c in coverage_term])
    term = zeros(Float64, (l,nb_reps))
    for i in 1:nb_reps
        term[1:length(coverage_term[i][chr]),i] = coverage_term[i][chr]
    end
    return vec(mean(term, dims=2))
end

function get_chr_from_wig(wig_file::String)
    chrs::Set{String} = Set()
    open(wig_file, "r") do file
        for line in eachline(file)
            startswith(line, "variableStep") && push!(chrs, split(split(line, " ")[2],"=")[2])
        end
    end
    return chrs
end

function approxoccursin(a::LongDNASeq, b::LongDNASeq; k=1)
    return approxsearch(a, b, k) != 0:-1
end


