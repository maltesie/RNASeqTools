const ANNOTATION_VCH = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff"
const GENOME_VCH = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"

function embl_to_gff(embl_file::String, gff_file::String, chrs::Vector{String})

    chroms = copy(reverse(chrs))
    writer = GFF3.Writer(open(gff_file, "w"))
    write(writer, GFF3.Record("##gff-version 3.2.1"))
    lines = readlines(embl_file)
    current_refname = pop!(chroms)
    for i in 1:length(lines)-1
        line = lines[i]
        next_line = lines[i+1]
        if startswith(line, "FT")

            if occursin(line, "source")
                _, _, r = split(line)
                _, l = split(r, "..")
                write(writer, GFF3.Record("##sequence-region $current_refname 1 $l"))
                write(writer, GFF3.Record("$current_refname\t.\tsource\t1\t$l\t.\t.\t.\tName=$current_refname"))
            end
            if occursin(line, "gene")
                _, _, r = split(line)
                strand = "+"
                startswith(r, "complement") && (r = r[12:end-1]; strand="-")
                l, r = split(r, "..")
                name = split(next_line, "locus_tag=")[end]
                write(writer, GFF3.Record("$current_refname\t.\tCDS\t$l\t$r\t.\t$strand\t.\tName=$name"))
            end
        end

    end
    close(writer)
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