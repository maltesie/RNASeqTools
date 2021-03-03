function rev_comp!(reads::Reads)
    for (key, read) in reads.dict
        BioSequences.reverse_complement!(read)
    end
end

function rev_comp!(reads::PairedReads; treat=:both)
    @assert treat in [:both, :read1, :read2]
    for (key, (read1, read2)) in reads.dict
        treat in [:both, :read1] && BioSequences.reverse_complement!(read1)
        treat in [:both, :read2] && BioSequences.reverse_complement!(read2)
    end
end

function cut!(read::LongDNASeq, pos::Int; keep=:left, from=:left)
    0 <= pos <= length(read) || (return nothing)
    (pos == 0 || pos == length(read)) && (copy!(read, LongDNASeq("")); return nothing)
    if (from == :left) && (keep == :left)
        copy!(read, @view(read[1:pos]))
    
    elseif (from == :left) && (keep == :right)
        copy!(read, @view(read[pos+1:end]))
    
    elseif (from == :right) && (keep == :left)
        copy!(read, @view(read[1:length(read)-pos]))
    
    elseif (from == :right) && (keep == :right)
        copy!(read, @view(read[length(read)-pos+1:end]))
    end
end

function cut!(reads::Reads, pos::Int; keep=:left, from=:left)
    for (key, read) in reads.dict
        cut!(read, pos; keep=keep, from=from)
    end
end

function cut!(reads::PairedReads, pos::Int; keep=:left, from=:left)
    for (key, (read1, read2)) in reads.dict
        ((pos > length(read1)) && (pos > length(read2))) && continue
        cut!(read1, pos; keep=keep, from=from)
        cut!(read2, pos; keep=keep, from=from)
    end
end

function cut!(reads::Reads, seq::LongDNASeq; keep=:left_of_query, from=:left)
    @assert keep in [:left_of_query, :right_of_query, :left_and_query, :right_and_query]
    @assert from in [:left, :right]
    for (key, read) in reads.dict
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

function cut!(reads::PairedReads, seq::LongDNASeq; keep=:left_of_query, treat=:both)
    @assert keep in [:left_of_query, :right_of_query, :left_and_query, :right_and_query]
    @assert treat in [:read1, :read2, :both]
    for (key, (read1, read2)) in reads.dict
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
            start2, stop2 = s1
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

function approxoccursin(s1::LongDNASeq, s2::LongDNASeq; k=1, check_indels=false)
    length(s2) < length(s1) && (return false)
    if !check_indels
        for i in 1:length(s2)-length(s1)+1
            (sum(@view(s2[i:i+length(s1)-1]) .!== s1) <= k) && (return true) 
        end
        return false
    else
        return approxsearch(s2, s1, k) != 0:-1
    end
end

function Base.filter!(f, reads::Reads)
    for (key, read) in reads.dict
        f(read) || delete!(reads.dict, key)
    end
end

function Base.filter!(f, reads::PairedReads; logic=:or)
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
