function combine_gffs(gff_files::Array{String, 1}; out_file="out.gff3")
    writer = open(out_file, "w")
    for gff_in_file in gff_files
        reader = open(gff_in_file)
        write(writer, read(reader, String))
        close(reader)
    end
    close(writer)
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

function datatable_columns(df::SubDataFrame)
    [(id=column, name=column) for column in names(df)[1:end-3]]
end

function datatabele_data(df::SubDataFrame)
    [Dict(name=>row[Symbol(name)] for name in names(df)[1:end-3]) for row in eachrow(df)]
end

function translate_positions!(pos::Vector{Int}, trans_dict::Dict{Int, Union{Int, Nothing}})
    for p in pos
        p = trans_dict[p]
    end
end

@inline function isproperpair(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_PROPER_PAIR) != 0
end

@inline function areconcordant(pos1::Int, pos2::Int, chr1::String, chr2::String, aux1::String, aux2::String; distance=1000)::Bool
    
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
            (abs(p1-p2) < distance) && (return true)
        end
    end
    return false
end

function translate_position_mg1655_to_bw25113(pos::Int)::Union{Nothing, Int}
    
    if abs(pos) < 66533
        return pos
    
    elseif (66533 <= abs(pos) <= 70072) 
        return nothing
    elseif (70072 < abs(pos) < 257900)
        return sign(pos) * (abs(pos) - (3540-27)) 
    
    elseif (257900 <= abs(pos) <= 258675)
        return nothing
    elseif (258675 < abs(pos) < 364356)
        return sign(pos) * (abs(pos) - ((3540-27)+776)) 
    
    elseif (364356 <= abs(pos) <= 366418)
        return nothing
    elseif (366418 < abs(pos) <= 371758)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585))) 
    elseif (371758 < abs(pos) < 1299495)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223)) 
    
    elseif (1299495 <= abs(pos) <= 1300693)
        return nothing
    elseif (1300693 < abs(pos) < 1978495)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199))

    elseif (1978495 <= abs(pos) <= 1979270)
        return nothing
    elseif (1979270 < abs(pos) < 2521007)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776))

    elseif (2173363 <= abs(pos) <= 2173364)
        return nothing
    elseif (2173364 < abs(pos) < 2521007)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2))

    elseif (2521007 <= abs(pos) <= 2521126)
        return nothing
    elseif (2521126 < abs(pos) <= 3560456)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120))
    elseif (3560456 < abs(pos) < 4093991)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1))

    elseif (4093991 <= abs(pos) <= 4097451)
        return nothing
    elseif (4097451 < abs(pos) < 4296282)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1+(3461-29)))

    elseif (4296282 <= abs(pos) <= 4296392)
        return nothing
    elseif (4296392 < abs(pos))
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1+(3461-29)+111))
    end
end

function is_bitstring_fasta(file::String)
    (endswith(file, ".fasta") || endswith(file, ".fasta.gz")) || (return false)
    f = endswith(file, ".fasta.gz") ? GzipCompressorStream(open(file, "r")) : open(file, "r")
    first_line = readline(f)
    close(f)
    ((length(first_line) == 65) && all([c in ['0', '1'] for c in first_line[2:end]])) && (return true)
    return false
end

function translation_dict(from_sequence::String, to_sequence::String)
    s1 = LongDNASeq(from_sequence)
    s2 = LongDNASeq(to_sequence)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
    res = alignment(pairalign(GlobalAlignment(), s1, s2, scoremodel))
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
