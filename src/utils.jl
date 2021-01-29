using FASTX
using CodecZlib
using BioSequences

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
    #open("/home/malte/Workspace/data/log$chr.txt", "a") do file
    #    println(file, join(["$i\t$a1\t$a2\t$m\t$((a2+a1)/2)" for (i, (a1,a2, m)) in enumerate(zip(tex[1:2050, 1], tex[1:2050, 2], vec(mean(tex[1:2050, :], dims=2))))], '\n'))
    #end
    #open("/home/malte/Workspace/data/log2$chr.txt", "a") do file
    #    println(file, join(["$i\t$a1\t$a2\t$((a2+a1)/2)" for (i, (a1,a2)) in enumerate(zip(coverage_tex[1][chr][1:2050], coverage_tex[2][chr][1:2050]))], '\n'))
    #end
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