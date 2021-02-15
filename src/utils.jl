function combine_gffs(gff_files::Array{String, 1}; out_file="out.gff3")
    writer = open(out_file, "w")
    for gff_in_file in gff_files
        reader = open(gff_in_file)
        write(writer, read(reader, String))
        close(reader)
    end
    close(writer)
end

function reverse_complement_reads(in_file::String, out_file::String)

    record1::FASTQ.Record = FASTQ.Record()
    record2::FASTQ.Record = FASTQ.Record()
    ide::String = ""
    des::String = ""
    seq::String = ""
    qua::String = ""

    if endswith(in_file, ".fastq.gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(in_file, "r")))
    elseif endswith(in_file, ".fastq")
        reader = FASTQ.Reader(open(in_file, "r"))
    else
        throw(ErrorException("Only .fastq and .fastq.gz files are supported."))
    end
    writer = FASTQ.Writer(GzipCompressorStream(open(out_file, "w")))
    
    while !eof(reader)
        read!(reader, record1)
        seq = convert(String, reverse_complement(LongDNASeq(FASTQ.sequence(record1))))
        record2 = FASTQ.Record(identifier(record1), description(record1), seq, reverse(quality(record1)))
        write(writer, record2)
    end
    close(reader)
    close(writer)
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

function get_position(pos::Int, chr::String, aux::String; reversed=false)::Tuple{Int, String}
    poss::Array{Int,1} = [pos]
    chrs::Array{String,1} = [chr]
    reversed ? factor = -1 : factor = 1
    for alignment in split(aux, ";")
        (alignment == "-") && (return pos, chr)
        isempty(alignment) && continue
        c, p, q, e = split(alignment, ",")
        push!(poss, factor * parse(Int, p))
        push!(chrs, String(c))
    end
    ind = sortperm(abs.(poss))[1]
    spos = poss[ind]
    schr = chrs[ind]
    return spos, schr
end

function bam_chromosome_lengths(reader::BAM.Reader)
    chr_lengths = Int[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_lengths, parse(Int, meta["LN"]))
    end
    return chr_lengths
end

function bam_chromosome_names(reader::BAM.Reader)
    chr_names = String[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_names, meta["SN"])
    end
    return chr_names
end

function translate_positions!(pos::Vector{Int}, trans_dict::Dict{Int, Union{Int, Nothing}})
    for p in pos
        p = trans_dict[p]
    end
end

@inline function isproperpair(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_PROPER_PAIR) != 0
end

@inline function areconcordant(pos1::Int, pos2::Int, chr1::String, chr2::String, aux1::String, aux2::String)::Bool
    
    poss1::Array{Int,1} = [pos1]
    chrs1::Array{String,1} = [chr1]
    for alignment in split(aux1, ";")
        (isempty(alignment) | (alignment == "-")) && continue
        push!(poss1, parse(Int, split(alignment, ",")[2]))
        push!(chrs1, String(split(alignment, ",")[1]))
    end
    
    poss2::Array{Int,1} = [pos2]
    chrs2::Array{String,1} = [chr2]
    for alignment in split(aux2, ";")
        (isempty(alignment) | (alignment == "-")) && continue
        push!(poss2, parse(Int, split(alignment, ",")[2]))
        push!(chrs2, split(alignment, ",")[1])
    end
    
    for (p1::Int, c1::String) in zip(poss1, chrs1)
        for (p2::Int, c2::String) in zip(poss2, chrs2)
            c1 != c2 && continue
            (abs(p1-p2) <1000) && (return true)
        end
    end
    return false
end

function get_single_fragment_set(bam_file::String; nb_reads::Int=-1)::Set{String}
    reader = open(BAM.Reader, bam_file)
    record::BAM.Record = BAM.Record()
    single_fragments::Array{String, 1} = []
    c::Int = 0 
    while !eof(reader)
        read!(reader, record)
        isproperpair(record) && push!(single_fragments, BAM.tempname(record))
        c += 1
        ((nb_reads > 0) & (c >= nb_reads)) && break 
    end
    close(reader)
    single_fragment_set::Set{String} = Set(single_fragments) 
    return single_fragment_set
end