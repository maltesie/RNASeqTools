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

function cut!(read::LongDNASeq, pos::Int; keep=:left, from=:left)
    0 < pos < length(read) || (return nothing)
    
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
        #reads.dict[key] = (read1, read2)
        #println(read1, " ", read2)
    end
end

function cut!(reads::Reads, seq::LongDNASeq; keep=:left_of_query)
    @assert keep in [:left_of_query, :right_of_query, :left_and_query, :right_and_query]
    for (key, read) in reads.dict
        (start, stop) = approxsearch(read, seq, 1)
        start == 0 && continue
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
        (treat in [:both, :read1]) ? (start1, stop1) = approxsearch(read1, seq, 1) : (start1, stop1) = (0,-1)
        (treat in [:both, :read2]) ? (start2, stop2)  = approxsearch(read2, seq, 1) : (start2, stop2) = (0,-1)
        if (start1 != 0)
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
        if (start2 != 0)
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
        #(start1 != 0 || start2 != 0) && (reads.dict[key] = (read1, read2))
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
        return approxsearch(a, b, k) != 0:-1
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
