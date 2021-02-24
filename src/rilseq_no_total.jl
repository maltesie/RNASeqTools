using BioSequences
using FASTX
using CodecZlib
using XAM
using Combinatorics
using HypothesisTests
using DataFrames
using CSV
using XLSX
using JSON

function trim_fastp(folder::String, names::Array{String,1}; fastp_bin="fastp")

    record = FASTQ.Record()
    in_files = [[abspath(folder, "demultiplexed", "$(name)_1.fastq.gz"), 
                abspath(folder, "demultiplexed", "$(name)_2.fastq.gz")] 
                for name in names]
    out_files = [[abspath(folder, "trimmed", "$(name)_1.fastq.gz"), 
                    abspath(folder, "trimmed", "$(name)_2.fastq.gz")] 
                    for name in names]
    
    reporter = open(abspath(folder, "reports", "quality.txt"), "w")
    
    c = 0
    for ((file1, file2), (out1, out2), name) in zip(in_files, out_files, names)
        
        cmd = `$fastp_bin --adapter_sequence AGATCGGAAG --adapter_sequence_r2 AGATCGGAAG -i $file1 -I $file2 -o $out1 -O $out2 -M 25 -l 25 --trim_poly_g 10 --cut_front --cut_tail --umi --umi_loc read1 --umi_len 9`# --trim_poly_g 10 -r --cut_right_window_size=6`
        run(cmd)
        cc = 0
        reader = FASTQ.Reader(GzipDecompressorStream(open(out1, "r")))
        while !eof(reader)
            read!(reader, record)
            c += 1
            cc += 1
        end
        report = "$name\t$(cc) * 2\n"
        write(reporter, report)
    end 
    report = "$(c) * 2 reads in total after quality trimming in $folder\n"
    write(reporter, report)
    close(reporter)

end

function isread1(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_READ1) != 0
end

function isreverse(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_REVERSE) != 0
end

function cut_fastq_reads(pe_bam::String, out_fastq1::String, out_fastq2::String; 
        rev_complement=false, cut_len::Int=25)

    fastq::FASTQ.Record = FASTQ.Record()
    fastq2::FASTQ.Record = FASTQ.Record()
    bam::BAM.Record = BAM.Record() 
    bam_reader = open(BAM.Reader, pe_bam)
    fastq_writer1::FASTQ.Writer = FASTQ.Writer(open(out_fastq1, "w"))
    fastq_writer2::FASTQ.Writer = FASTQ.Writer(open(out_fastq2, "w"))
    writer = [fastq_writer1, fastq_writer2]
    des::String, seq::String, qua::Array{UInt8, 1} = "", "", []
    
    while !eof(bam_reader)
        read!(bam_reader, bam)
        des, seq, qua = BAM.tempname(bam), BAM.sequence(bam), BAM.quality(bam)
        writeto = 1
        if rev_complement == isread1(bam)
            if !isreverse(bam)
                seq = convert(String, reverse_complement(LongDNASeq(seq[end-cut_len+1:end])))
                qua = reverse(qua[end-cut_len+1:end])
            end
            writeto = 2
        else
            if isreverse(bam)
                seq = convert(String, reverse_complement(LongDNASeq(seq[1:cut_len-1])))
                qua = reverse(qua[1:cut_len-1])
            end
        end
        fastq = FASTQ.Record(des, seq, qua)
        write(writer[writeto], fastq)
    end
    close(bam_reader)
    close(fastq_writer1)
    close(fastq_writer2)
end

function align(folder::String, fasta_genome::String, 
                names::Array{String,1}; rev_complement=false, index_genome=true, 
                se_miss=3, pe_miss=3, se_length=25, bwa_bin="bwa", sam_bin="samtools")

    in_files = [[abspath(folder, "trimmed", "$(name)_1.fastq.gz"), 
                    abspath(folder, "trimmed", "$(name)_2.fastq.gz")] 
                    for name in names]
    
    pe_files = [abspath(folder, "pe_bams", "$(name).bam") for name in names]
    
    se_files = [[abspath(folder, "se_bams", "$(name)_1.bam"), 
                    abspath(folder, "se_bams", "$(name)_2.bam")] 
                    for name in names]
    
    reporter = open(abspath(folder, "reports", "pe_alignment.txt"), "w")
    record1::BAM.Record = BAM.Record()
    
    for ((file1, file2), pe_bam, (se_bam1, se_bam2), name) in zip(in_files, pe_files, se_files, names)
        
        cmd = pipeline(`$bwa_bin index -a is $fasta_genome`, stdout=nothing)
        index_genome && run(cmd)
        cmd = pipeline(`$bwa_bin aln -n $pe_miss -t 8 -R 200 $fasta_genome $file1`, stdout="tmp1.sai")
        run(cmd)
        cmd = pipeline(`$bwa_bin aln -n $pe_miss -t 8 -R 200 $fasta_genome $file2`, stdout="tmp2.sai")
        run(cmd)
        cmd = pipeline(`$bwa_bin sampe -a 1500 -P $fasta_genome tmp1.sai tmp2.sai $file1 $file2`, 
            stdout="tmp.bwa")
        run(cmd)
        cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
        run(cmd)
        cmd = pipeline(`$sam_bin sort tmp.view -o $pe_bam`, stdout=nothing)
        run(cmd)
        cmd = pipeline(`$sam_bin index $pe_bam`, stdout=nothing)
        run(cmd)
        
        cut_fastq_reads(pe_bam, "tmp1.fastq", "tmp2.fastq";
            rev_complement=rev_complement, cut_len=se_length)
        
        cmd = pipeline(`$bwa_bin aln -n $se_miss -t 6 -R 200 $fasta_genome tmp1.fastq`, stdout="tmp1.sai")
        run(cmd)
        cmd = pipeline(`$bwa_bin samse $fasta_genome tmp1.sai tmp1.fastq`, stdout="tmp.bwa")
        run(cmd)
        cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
        run(cmd)
        cmd = pipeline(`$sam_bin sort tmp.view -o $se_bam1`)
        run(cmd)
        cmd = pipeline(`$sam_bin index $se_bam1`)
        run(cmd)

        cmd = pipeline(`$bwa_bin aln -n $se_miss -t 6 -R 200 $fasta_genome tmp2.fastq`, stdout="tmp2.sai")
        run(cmd)
        cmd = pipeline(`$bwa_bin samse $fasta_genome tmp2.sai tmp2.fastq`, stdout="tmp.bwa")
        run(cmd)
        cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
        run(cmd)
        cmd = pipeline(`$sam_bin sort tmp.view -o $se_bam2`)
        run(cmd)
        cmd = pipeline(`$sam_bin index $se_bam2`)
        run(cmd)
        
        rm("tmp1.sai")
        rm("tmp2.sai")
        rm("tmp.bwa")
        rm("tmp.view")
        rm("tmp1.fastq")
        rm("tmp2.fastq")
        
    end 
end

@inline function get_position(pos::Int, chr::String, aux::String, nm::Int)::Tuple{Int, String}
    (aux == "-") && (return pos, chr)
    spos::Int = pos
    schr::String = chr
    temp_pos::Int = pos
    temp_chr::String = chr
    temp_nm::Int = 10000
    ct = 1
    i_t = 1
    for (i::Int, c::Char) in enumerate(aux)
        if c == ','
            if (ct == 1)
                temp_chr = aux[i_t:i-1]
            elseif (ct == 2) 
                temp_pos = parse(Int, aux[i_t:i-1])
            elseif (ct==4)
                temp_nm = parse(Int, aux[i_t:i-1])
            end
            i_t = i + 1
            ct += 1
        elseif c == ';'
            if (abs(temp_pos) < abs(spos)) && (temp_nm <= nm)
                spos = temp_pos
                schr = temp_chr
            end
            ct = 1
            i_t = i + 1
        end
    end
    #startswith(schr, "0;") && println(spos, " ", schr, " ", aux)
    return spos, schr
end

function write_chimeric_fragments(read1_names::Array{String,1}, read1_poss::Array{Int,1}, 
        read1_chrs::Array{String,1}, read1_auxs::Array{String,1}, read1_nms::Array{Int,1}, 
        read2_names::Array{String,1}, read2_poss::Array{Int,1}, read2_chrs::Array{String,1}, 
        read2_auxs::Array{String,1}, read2_nms::Array{Int,1}, single_fragment_set::Set{String}, 
        out_file::String)
    
    nb_single1::Int, nb_single2::Int, nb_chimeric::Int = 0, 0, 0
    read2_trans::Dict{String, Int} = Dict(n => i for (i,n) in enumerate(read2_names))
    writer = GzipCompressorStream(open(out_file, "w"), level=2)
    for (read1_name, read1_pos, read1_chr, read1_aux, read1_nm) in 
        zip(read1_names, read1_poss, read1_chrs, read1_auxs, read1_nms)
        read1_name in keys(read2_trans) || continue
        i2 = read2_trans[read1_name]
        read2_pos::Int, read2_chr::String, read2_aux::String, read2_nm::Int = read2_poss[i2], read2_chrs[i2], read2_auxs[i2], read2_nms[i2]
        sp1, sc1 = get_position(read1_pos, read1_chr, read1_aux, read1_nm)
        sp2, sc2 = get_position(read2_pos, read2_chr, read2_aux, read2_nm)
        if (read1_name in single_fragment_set)
            write(writer, "$read1_name\tsingle\t$sp1\t$sp2\t$sc1\t$sc2\n")
            nb_single1 += 1
        elseif areconcordant(read1_pos, read2_pos, read1_chr, read2_chr, read1_aux, read2_aux)
            write(writer, "$read1_name\tsingle\t$sp1\t$sp2\t$sc1\t$sc2\n")
            nb_single2 += 1
        else
            write(writer, "$read1_name\tchimeric\t$sp1\t$sp2\t$sc1\t$sc2\n")
            nb_chimeric += 1
        end
    end
    close(writer)
    return nb_single1, nb_single2, nb_chimeric
end

@inline function isproperpair(record::BAM.Record)::Bool
    return (BAM.flag(record) & SAM.FLAG_PROPER_PAIR) != 0
end

function write_file(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

function strandint(record::BAM.Record; is_rev=false)
    BAM.ispositivestrand(record) ? strand = 1 : strand = -1
    is_rev ? (return strand * -1) : (return strand)
end

@inline function translated_data(data::Array{UInt8,1})::String
    for i in 1:length(data)
        (data[i] == 0x00) && (return String(data[1:i-1]))
    end
end

@inline function get_NM_tag(data::Array{UInt8,1})::Int
    for i in 1:length(data)-2
      (0x4d == data[i]) & (0x43 == data[i+1]) && (return Int(data[i+2]))
    end
    return -1
end

@inline function get_XA_tag(data::Array{UInt8,1})::String
    t = UInt8['X', 'A']
    for i in 1:length(data)-2
        (0x00 == data[i]) & (t[1] == data[i+1]) & (t[2] == data[i+2]) && 
        (return translated_data(data[i+4:end]))
    end
    return "-"
end

function merge_bam_row!(a::DataFrameRow, b::DataFrameRow)
    (start1, stop1, start2, stop2) = sort([a[:start], a[:stop], b[:start], b[:stop]])
    a[:start] = start1
    a[:stop] = stop2
    a[:nm] = a[:nm] + b[:nm] + start2 - stop1
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

function read_bam(bam_file::String; nb_reads::Int = -1)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    aligned_reads = DataFrame(name=String[], start=Int[], stop=Int[], chr=String[], nm=Int[], aux=String[], cigar=String[])
    c::Int = 0
    while !eof(reader)
        read!(reader, record)
        !BAM.ismapped(record) && continue
        (start, stop) = sort([BAM.position(record)*strandint(record), BAM.rightposition(record)*strandint(record)])
        aux_data = BAM.auxdata(record).data
        new_row = DataFrame(name=BAM.tempname(record), start=start, stop=stop, chr=BAM.refname(record), 
                            nm=get_NM_tag(aux_data), aux=get_XA_tag(aux_data), cigar=BAM.cigar(record))
        append!(aligned_reads, new_row)
        ((nb_reads > 0) & (c >= nb_reads)) && break 
    end
    close(reader)
    return aligned_reads[:,:name], aligned_reads[:,:start], aligned_reads[:,:chr], aligned_reads[:,:aux], aligned_reads[:,:nm]
end

function all_interactions(folder::String, fasta_genome::String, names::Array{String}; 
                        generate_bams=true, cut_len=25, read1_reversed=false, read2_reversed=false)
    
    read1_names::Array{String,1}, read1_chrs::Array{String,1}, read1_auxs::Array{String,1} = [], [], []
    read1_poss::Array{Int,1}, read1_nms::Array{Int,1} = [], []
    read2_names::Array{String,1}, read2_chrs::Array{String,1}, read2_auxs::Array{String,1} = [], [], []
    read2_poss::Array{Int,1}, read2_nms::Array{Int,1} = [], []
    
    single_fragment_set::Set{String} = Set()

    fastq_files = [[abspath(folder, "trimmed", "$(name)_1.fastq.gz"), 
                    abspath(folder, "trimmed", "$(name)_2.fastq.gz")] 
                    for name in names]
    
    pe_files = [abspath(folder, "pe_bams", "$(name).bam") for name in names]
    
    se_files = [[abspath(folder, "se_bams", "$(name)_1.bam"), 
                    abspath(folder, "se_bams", "$(name)_2.bam")] 
                    for name in names]
    
    chimeric_files = [abspath(folder, "results", "$(name)_raw.txt.gz") for name in names]
    
    reporter = open(abspath(folder, "reports", "all_pairs.txt"), "w")
    total_count::Int = 0
    
    for (pe_file, (se_file1, se_file2), (fastq_file1, fastq_file2), chimeric_file, name) in 
        zip(pe_files, se_files, fastq_files, chimeric_files, names)
        
        single_fragment_set = Set()#get_single_fragment_set(pe_file)

        read1_names, read1_poss, read1_chrs, read1_auxs, read1_nms = read_bam(se_file1)

        read2_names, read2_poss, read2_chrs, read2_auxs, read2_nms = read_bam(se_file2)

        single_count1, single_count2, chimeric_count = write_chimeric_fragments(read1_names, 
            read1_poss, read1_chrs, read1_auxs, read1_nms, read2_names, read2_poss, read2_chrs, 
            read2_auxs, read2_nms, single_fragment_set, chimeric_file)

        total_count += single_count1 + single_count2 + chimeric_count
        report = "$name\t$single_count1\t$single_count2\t$chimeric_count\n"
        write(reporter, report)
    end
    report = "Found $total_count transcripts in total for $folder.\n"
    write(reporter, report)
    close(reporter)
end

Coord = Pair{Int, String}
Interact = Pair{Coord, Coord}
MyMatrix = Array{Array{Int,1}, 1}

function min_pv_region(interaction::Interact, sum_pairs::Dict{Interact, Int}, count_reads1::Dict{Coord,Int},
    count_reads2::Dict{Coord,Int}, tot_count::Int; seglen=100, minints=5, maxsegs=5, minodds=1.0)

    (coord1, coord2) = interaction
    nb_coords = 2*maxsegs+1
    coordset1 = [Coord(coord1[1]+(c*seglen), coord1[2]) for c in -maxsegs+1:maxsegs-1] 
    coordset2 = [Coord(coord2[1]+(c*seglen), coord2[2]) for c in -maxsegs+1:maxsegs-1]
    best_parameter::Array{Float64,1} = [2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    f = FisherExactTest(1,1,1,1)
    for (length1, length2) in with_replacement_combinations(1:maxsegs, 2)
        for (start1, start2) in with_replacement_combinations(1:maxsegs-1, 2)
            count1, count2, count_int = 1, 1, 1
            for sliding_coord1 in coordset1[start1:start1+length1]
                sliding_coord1 in keys(count_reads1) && (count1 += count_reads1[sliding_coord1])
            end
            count1 == 0 && continue
            for sliding_coord2 in coordset2[start2:start2+length2]
                sliding_coord2 in keys(count_reads2) && (count2 += count_reads2[sliding_coord2])
                for sliding_coord1 in coordset1[start1:start1+length1]
                    sliding_interaction = Interact(sliding_coord1, sliding_coord2)
                    sliding_interaction in keys(sum_pairs) && (count_int += sum_pairs[sliding_interaction])
                end
            end
            count_int == 0 && continue
            a = count_int 
            b = count1-count_int
            c = count2-count_int
            d = tot_count-count1-count2+count_int
            #println("$a,$b,$c,$d")
            b == 0 | c == 0 ? odds = Inf : odds = (a*d)/(b*c)
            odds < minodds && continue
            try
                f = FisherExactTest(a,b,c,d)
            catch e
                continue
            end
            pv = pvalue(f)
            pv > best_parameter[1] && continue
            (pv == best_parameter[1]) & (best_parameter[3] > count_int) && continue
            best_parameter = [pv, odds, count_int, b, c, d, 
                coord1[1]-(maxsegs+start1)*seglen, length1*seglen, 
                coord2[1]-(maxsegs+start2)*seglen, length2*seglen]
        end
    end
    return best_parameter 
end

function get_annotations(gff_file::String; utr_length=150)::Dict{String, DataFrame}
    reader = open(gff_file)
    gff = read(reader, String)
    results::Dict{String, DataFrame} = Dict()

    for line in split(gff, "\n")
        startswith(line, "#") | isempty(line) && continue
        chr, _, typ, start, stop, _, strand, _, aux = String.(split(line, "\t"))
        (typ in ["Gene", "gene"]) || continue
        if strand == "-"
            start_int = parse(Int, stop) * -1
            stop_int = parse(Int, start) * -1
        else
            start_int = parse(Int, start)
            stop_int = parse(Int, stop)
        end
        if (typ == "gene") 
            name = String(split(aux, ";")[2][6:end])
            if startswith(name, "CCNA_") 
                (name[6] == 'R') && (typ = "srna")
            end
        end
        (typ == "Gene") && (name = String.(split(aux, ";")[1][6:end]); typ="gene")
        (occursin("locus_tag", aux) || ("VC" == aux[6:7])) || (typ = "srna")
        row = DataFrame(name=name, start=start_int, stop=stop_int, typ=typ)
        chr in keys(results) ? append!(results[chr], row) : results[chr] = row
    end

    for (chr, data) in results
        sort!(data, :start)
        tmp_data = data[data.typ .!= "srna", :]
        utr5 = DataFrame(name="U5." .* tmp_data.name, start=tmp_data.start .- utr_length, stop=tmp_data.start, typ=tmp_data.typ)
        utr3 = DataFrame(name="U3." .* tmp_data.name, start=tmp_data.stop, stop=tmp_data.stop .+ utr_length, typ=tmp_data.typ)
        igr_type = Array{String, 1}(undef, length(tmp_data.name)-1) .= "igr"
        igr = DataFrame(name="IG." .* tmp_data.name[1:end-1] .* "|" .* tmp_data.name[2:end], 
            start=tmp_data.stop[1:end-1] .+ utr_length, stop=tmp_data.start[2:end] .- utr_length, typ=igr_type)
        igr= igr[igr.start .< igr.stop, :]
        #igr = DataFrame(name=igr.name[index], start=igr.start[index], stop=igr.stop[index])
        append!(data, utr5)
        append!(data, utr3)
        append!(data, igr)
        sort!(data, :start)
    end
    
    return results
end

@inline function get_name(chr::String, pos::Int, annotations::Dict{String, DataFrame})::String
    srna = annotations[chr].name[(annotations[chr].start .<= pos .<= annotations[chr].stop) .& (annotations[chr].typ .=="srna")]
    isempty(srna) || return srna[1]
    names = annotations[chr].name[annotations[chr].start .<= pos .<= annotations[chr].stop]
    isempty(names) ? (return "not found") : (return names[1])
end

@inline function get_gene(chr::String, pos_start::Int, pos_end::Int,
    annotations::Dict{String, DataFrame}; utr_length=150)::Tuple{String, String, Float64}
    gene = "not found"
    type = "not found"
    srna = "none"
    max_overlap = 0.
    best_pos = -1.
    for annotation in eachrow(annotations[chr])
        (startswith(annotation[:name], "U3.") || startswith(annotation[:name], "U5.")) && continue
        if startswith(annotation[:name], "IG.")
            ((pos_end <= annotation[:start]) || (pos_start >= annotation[:stop])) && continue
            overlap = min(pos_end, annotation[:stop]) - max(pos_start, annotation[:start]) + 1
            position = min(max((round(Int, (pos_start + pos_end) / 2.0) - annotation[:start]) / (annotation[:stop] - annotation[:start]), 0.0), 1.0)
        else
            ((pos_end <= annotation[:start]-utr_length) || (pos_start >= annotation[:stop]+utr_length)) && continue
            overlap = min(pos_end, annotation[:stop]+utr_length) - max(pos_start, annotation[:start]-utr_length) + 1
            position = min(max((round(Int, (pos_start + pos_end) / 2.0) - annotation[:start] - utr_length) / (annotation[:stop] - annotation[:start] + 2*utr_length), 0.0), 1.0)
        end
        if (overlap > max_overlap)
            max_overlap = overlap
            gene = annotation[:name] 
            type = annotation[:typ]
            best_pos = position
        end
        ((annotation[:typ] == "srna") && (overlap>0)) && (srna = annotation[:name])
    end
    (srna != "none") && return (srna, "srna", -1.0)
    return (gene, type, best_pos)
end

function get_full_name(chr::String, pos_start::Int, pos_end::Int,
        annotations::Dict{String, DataFrame}; utr_length=150)::Tuple{String, String, String, Float64}
    
    gene::String, type::String, position::Float64 = get_gene(chr, pos_start, pos_end, annotations; utr_length=utr_length)
    
    name1::String = get_name(chr, pos_start, annotations)
    name2::String = get_name(chr, pos_end, annotations)
    (name1 == name2) ? (return name1, gene, type, position) : (return "$name1 -> $name2", gene, type, position)
end

function get_rrna(gff_file::String)::Tuple{Array{Int, 2}, Dict{String,Int}}

    reader = open(gff_file)
    gff = read(reader, String)
    chr_trans::Dict{String,Int} = Dict()
    strand_trans::Dict{Bool, Int} = Dict(false=>-1, true=>1)
    results = Array{Int, 2}(undef, 0, 3)

    for line in split(gff, "\n")
        startswith(line, "#") | isempty(line) && continue
        chr, _, typ, start, stop, _, strand, _, _ = split(line, "\t")
        chr in keys(chr_trans) || (chr_trans[chr] = length(chr_trans)+1)
        (typ == "rRNA") | (typ == "tRNA") || continue
        strand=="+" ? results = [results; [chr_trans[chr], parse(Int, start), parse(Int, stop)]'] :
                        results = [results; [chr_trans[chr], parse(Int, stop)*-1, parse(Int, start)*-1]']
    end
    close(reader)
    #print(results)
    return results, chr_trans
end

function check_rrna(pos1::Int, pos2::Int, chr1::String, chr2::String, 
        rrna::Array{Int, 2}, chr_trans::Dict{String, Int})::Bool
    
    if chr1 in keys(chr_trans)
        any((chr_trans[chr1] .== rrna[:,1]) .& (pos1 + 2 .>= rrna[:,2]) .& (pos1 .<= rrna[:,3])) && (return true)
    end
    if chr2 in keys(chr_trans)
        any((chr_trans[chr2] .== rrna[:,1]) .& (pos2 + 2 .>= rrna[:,2]) .& (pos2 .<= rrna[:,3])) && (return true)
    end
    return false
end

function process_interactions(interaction_file::String, rrna::Array{Int, 2}, chr_trans::Dict{String, Int};
    seglen=100, minints=5)

    reader = GzipDecompressorStream(open(interaction_file, "r"))
    region_interactions::Dict{Interact, Array{Pair{Int, Int},1}} = Dict()
    count_reads1::Dict{Coord, Int} = Dict()
    count_reads2::Dict{Coord, Int} = Dict()
    count_total = 0
    count_rrna = 0
    #for (key, value) in region_interactions
    #    ((1293000 < abs(key[1][1]) < 1294000) || (1293000 < abs(key[2][1]) < 1294000)) && println("inside1 $key $(length(value))")
    #end
    while !eof(reader)
        line = String.(split(readline(reader)))
        name = line[1] 
        typ = line[2]
        (typ == "chimeric") || continue
        pos1 = parse(Int, line[3]) 
        pos2  = parse(Int, line[4]) 
        chr1 = line[5] 
        chr2 = line[6] 
        #((1293900 < abs(pos1) < 1294000) || (1293900 < abs(pos2) < 1294000)) && println("before check pos 1 $pos1 pos 2 $pos2") 
        check_rrna(pos1, pos2, chr1, chr2, rrna, chr_trans) && (continue; count_rrna+=1)
        #((1293900 < abs(pos1) < 1294000) || (1293900 < abs(pos2) < 1294000)) && println("after check pos 1 $pos1 pos 2 $pos2") 
        count_total += 1
        seg1::Int = (pos1÷seglen)*seglen
        seg1<0 && (seg1-=seglen)
        seg2::Int = (pos2÷seglen)*seglen
        seg2<0 && (seg2-=seglen)
        coord1 = Coord(seg1, chr1)
        coord2 = Coord(seg2, chr2)
        interaction = Interact(coord1, coord2)
        #((interaction[1][1] == 1293900) || (interaction[2][1] == 1293900)) && println(interaction, " ", pos1, " ", pos2)
        
        (interaction in keys(region_interactions)) ? push!(region_interactions[interaction], pos1=>pos2) : (region_interactions[interaction] = [pos1=>pos2])

        coord1 in keys(count_reads1) ? 
        (count_reads1[coord1] += 1) : 
        push!(count_reads1, coord1=>1)
        coord2 in keys(count_reads2) ?
        (count_reads2[coord2] += 1) : 
        push!(count_reads2, coord2=>1)

        #(((361800 => "NC_002505") => (1293900 => "NC_002505")) in keys(region_interactions)) && println("hey")
    end
    close(reader)
    println(count_total, " ", count_rrna)
    sum_pairs::Dict{Interact, Int} = Dict()
    for interaction in keys(region_interactions)
        nb_interactions = length(region_interactions[interaction])
        (nb_interactions >= minints) && (sum_pairs[interaction] = nb_interactions)
    end

    return region_interactions, sum_pairs, count_reads1, count_reads2, count_total
end


function significant_chimeras(folder::String, gff_genome::String, names::Array{String, 1}; 
            seglen=100, minints=5, minodds=1.0, utr_length=150)
    
    report = ""
    rrna::Array{Int, 2}, chr_trans::Dict{String,Int} = get_rrna(gff_genome)

    annotations = get_annotations(gff_genome, utr_length=utr_length)

    header = ["Name RNA1", "pos1", "Name RNA2", "pos2", "# of chimeric fragments", "library", "norm. odds ration", "p-value",
            "IP interactions RNA1", "IP interactions RNA2", "RNA1 chromosome", "RNA1 first start", "RNA1 last start",
            "RNA2 chromosome", "RNA2 first start", "RNA2 last start", "b", "c", "d"]
    
    interaction_files =[abspath(folder, "results", "$(name)_raw.txt.gz") for name in names]
    
    out_files =[abspath(folder, "results", "$(name).csv") for name in names]
    c=0

    for (ip_interaction_file, out_file) in zip(interaction_files, out_files)
        
        ip_interactions, ip_interaction_count, ip_reads1_count, ip_reads2_count, ip_total_count = 
        process_interactions(ip_interaction_file, rrna, chr_trans)
        
        #for (key, value) in ip_interactions
        #    ((1293000 < abs(key[1][1]) < 1294000) || (1293000 < abs(key[2][1]) < 1294000)) && println("outside $key $(length(value))")
        #end

        results = DataFrame(name1=String[], pos1s=Int[], name2=String[], pos2s=Int[], nb_fragments=Int[], libs=String[], 
            norm_odds=Float64[], p_value=Float64[], ip_int1=Float64[], ip_int2=Float64[], chromosome1=String[], first_start1=Int[], last_start1=Int[], chromosome2=String[], 
            first_start2=Int[], last_start2=Int[], b=Int[], c=Int[], d=Int[])
        
        for (interaction, interaction_count) in sort(collect(ip_interaction_count), by=x->x[2], rev=true)

            (coord1, coord2) = interaction
            (pos1, chr1), (pos2, chr2) = coord1, coord2
            
            a = interaction_count
            b = ip_reads1_count[coord1] - interaction_count
            c = ip_reads2_count[coord2] - interaction_count
            d = ip_total_count - ip_reads1_count[coord1] - ip_reads2_count[coord2] + interaction_count
            
            b == 0 | c == 0 ? odds = Inf : odds = (a*d)/(b*c)
            odds < minodds && continue
            
            f = FisherExactTest(a,b,c,d)
            norm_odds = round(f.ω, digits=2)
            pv = round(pvalue(f), digits=10)
            #(pv > 0.05) && continue
            
            name1 = get_name(chr1, pos1, annotations)
            
            name2 = get_name(chr2, pos2, annotations)
            
            lib = "tbd"
            
            ip_int1 = ip_reads1_count[coord1]/ip_total_count
            ip_int2 = ip_reads2_count[coord2]/ip_total_count

            starts1 = sort(ip_interactions[interaction], by=x->x[1])
            first_start1, last_start1 = starts1[1][1], starts1[end][1]
            starts2 = sort(ip_interactions[interaction], by=x->x[2])
            first_start2, last_start2 = starts2[1][2], starts2[end][2]

            row = DataFrame(name1=name1, pos1s=pos1, name2=name2, pos2s=pos2, nb_fragments=interaction_count, 
                libs=lib, norm_odds=norm_odds, p_value=pv, ip_int1=ip_int1, ip_int2=ip_int2, 
                chromosome1=chr1, first_start1=first_start1, last_start1=last_start1, 
                chromosome2=chr2, first_start2=first_start2, last_start2=last_start2, 
                b=b, c=c, d=d)

            append!(results, row)
        end
        CSV.write(out_file, results, header=header)
    end
end

function exaustive_close_indices(pos1::Int, pos2::Int, poss::MyMatrix; expansion_step=200)
    diff1::Array{Int, 1} = [-expansion_step, expansion_step]
    diff2::Array{Int, 1} = [-expansion_step, expansion_step]
    indices::Array{Int, 1} = []
    done = false
    while !done
        indices_check = length(indices)
        for (i,(a, b)) in enumerate(poss)
            if (diff1[1] <= pos1-a <= diff1[2]) & (diff2[1] <= pos2-b <= diff2[2])
                !(i in indices) && push!(indices, i)
            end
        end
        indices_check == length(indices) && (done = true)
        #println(poss[indices])
        d1::Array{Int,1} = [p[1]-pos1 for p in poss[indices]]
        diff1 = [minimum(d1)-expansion_step,  maximum(d1)+expansion_step]
        d2::Array{Int,1} = [p[2]-pos2 for p in poss[indices]]
        diff2 = [minimum(d2)-expansion_step,  maximum(d2)+expansion_step]
        #println(diff1, " ", diff2)
    end
    return indices
end

function merged_table(table::CSV.File, annotations::Dict{String,DataFrame}; utr_length=150)::DataFrame
    poss::MyMatrix = [[row[:pos1], row[:pos2]] for row in table]
    already_found::Set{Int} = Set() 
    merged_data = DataFrame(name1=String[], gene1=String[], name2=String[], gene2=String[], 
        type1=String[], type2=String[], relpos1=Float64[], relpos2=Float64[], nb_fragments=Int[], 
        libs=String[], norm_odds=Float64[], p_value=Float64[], chromosome1=String[], 
        first_start1=Int[], last_start1=Int[], chromosome2=String[], first_start2=Int[], last_start2=Int[])
    for (i, row) in enumerate(table)
        (i in already_found) && continue
        pos1, pos2 = row[:pos1], row[:pos2]
        indices = exaustive_close_indices(pos1, pos2, poss)
        for ind in indices
            push!(already_found, ind)
        end
        
        nb = sum([r[Symbol("# of chimeric fragments")] for r in table[indices]])
        lib = "tbd"
        chr1 = [r[Symbol("RNA1 chromosome")] for r in table[indices]][1]
        chr2 = [r[Symbol("RNA2 chromosome")] for r in table[indices]][1]
        first1 = minimum([r[Symbol("RNA1 first start")] for r in table[indices]])
        last1 = maximum([r[Symbol("RNA1 last start")] for r in table[indices]])
        first2 = minimum([r[Symbol("RNA2 first start")] for r in table[indices]])
        last2 = maximum([r[Symbol("RNA2 last start")] for r in table[indices]])
        ips1 = sum([r[Symbol("IP interactions RNA1")] for r in table[indices]])
        ips2 = sum([r[Symbol("IP interactions RNA2")] for r in table[indices]])
            
        n1, gene1, type1, rel_pos1 = get_full_name(chr1, first1, last1, annotations; utr_length=utr_length)
        n2, gene2, type2, rel_pos2 = get_full_name(chr2, first2, last2, annotations; utr_length=utr_length)
    
        b = sum([r[:b] for r in table[indices]])
        c = sum([r[:c] for r in table[indices]])
        d = ([r[:d] for r in table[indices]][1] + 
            [r[:c] for r in table[indices]][1] +
            [r[:b] for r in table[indices]][1] +
            [r[Symbol("# of chimeric fragments")] for r in table[indices]][1]) - (nb + b + c)
        
        f = FisherExactTest(nb,b,c,d)
        norm_odds = round(f.ω, digits=2)
        pv = round(pvalue(f), digits=10)
        #(pv > 0.05) && continue
        
        row = DataFrame(name1=n1, name2=n2, gene1=gene1, gene2=gene2, type1=type1,
        type2=type2, relpos1=rel_pos1, relpos2=rel_pos2, nb_fragments=nb, libs=lib, 
        norm_odds=norm_odds, p_value=pv, chromosome1=chr1, first_start1=first1, 
        last_start1=last1, chromosome2=chr2, first_start2=first2, last_start2=last2)
        append!(merged_data, row)
    end
    return merged_data
end

function match_id(first1::Int, last1::Int, first2::Int, last2::Int, ranges::MyMatrix)
    for (i, range) in enumerate(ranges)
        ((last1 <= range[1]) || (first1 >= range[2])) && continue
        ((last2 <= range[3]) || (first2 >= range[4])) && continue
        overlap1 = min(last1, range[2]) - max(first1, range[1]) + 1
        overlap2 = min(last2, range[4]) - max(first2, range[3]) + 1
        combined_length1 = max(last1, range[2]) - min(first1, range[1]) + 1
        combined_length2 = max(last2, range[4]) - min(first2, range[3]) + 1
        (((overlap1/combined_length1) > 0.3) && ((overlap2/combined_length2) > 0.3)) && (return i)
    end
    return -1
end

function unified_table(table1::CSV.File, table2::CSV.File, annotations::Dict{String,DataFrame}; utr_length=150)::DataFrame
    
    merged_data1 = merged_table(table1, annotations; utr_length=utr_length)
    merged_data2 = merged_table(table2, annotations; utr_length=utr_length)
    range1 = [[merged_data1[i,:first_start1], merged_data1[i,:last_start1],
            merged_data1[i,:first_start2], merged_data1[i,:last_start2]] for i in 1:length(merged_data1[!,:name1])]
    range2 = [[merged_data2[i,:first_start1], merged_data2[i,:last_start1],
            merged_data2[i,:first_start2], merged_data2[i,:last_start2]] for i in 1:length(merged_data2[!,:name1])]
    
    for row in eachrow(merged_data1)
        !isempty(range2) ? (id2 = match_id(row[:first_start1], row[:last_start1], row[:first_start2], row[:last_start2], range2)) : (id2=-1)
        if (id2 >= 0)
            row[:libs] = "1,2"
            row[:nb_fragments] += merged_data2[id2,:nb_fragments]
            row[:nb_fragments] = round(Int, row[:nb_fragments] / 2.0)
            row[:first_start1] = min(row[:first_start1], merged_data2[id2,:first_start1])
            row[:last_start1] = max(row[:last_start1], merged_data2[id2,:last_start1])
            row[:first_start2] = min(row[:first_start2], merged_data2[id2,:first_start2])
            row[:last_start2] = max(row[:last_start2], merged_data2[id2,:last_start2])
            row[:p_value] = max(row[:p_value], merged_data2[id2,:p_value])
            row[:norm_odds] = min(row[:norm_odds], merged_data2[id2,:norm_odds])
        else
            row[:libs] = "1"
        end
    end
    for row in eachrow(merged_data2)
        !isempty(range1) ? (id1 = match_id(row[:first_start1], row[:last_start1], row[:first_start2], row[:last_start2], range1)) : (id1=-1)
        (id1 == -1) && (row[:libs] = "2"; push!(merged_data1, row))
    end
    sort!(merged_data1, :nb_fragments, rev=true)
    return merged_data1
    #CSV.write(out_file, merged_data1)
end

function postprocess(folder::String, gff_genome::String, names::Array{String, 1}, fname::String; utr_length=150, export_single_csv=true)
     
    table_files =[[abspath(folder, "results", "$(names[i]).csv"),
                    abspath(folder, "results", "$(names[i+1]).csv")] for i in 1:2:length(names)-1]
    
    unified_sheet_names = ["$(names[i])" for i in 1:2:length(names)-1]
    
    annotations = get_annotations(gff_genome, utr_length=utr_length)
    XLSX.openxlsx(abspath(fname), mode="w") do xf
        for (i, ((file_rep1, file_rep2), sheet_name)) in enumerate(zip(table_files, unified_sheet_names))
            table1, table2 = CSV.File(file_rep1), CSV.File(file_rep2)
            unified_tab = unified_table(table1, table2, annotations; utr_length=utr_length)
            export_single_csv && CSV.write(abspath(folder, sheet_name * ".csv"), unified_tab)
            (XLSX.sheetcount(xf) < i) &&  XLSX.addsheet!(xf, "new")
            sheet = xf[i]
            XLSX.rename!(sheet, sheet_name)
            XLSX.writetable!(sheet, [c for c in eachcol(unified_tab)], DataFrames.names(unified_tab))
        end
    end
end

function rilseq_analysis(lib_names::Array{String,1}, rilseq_folder::String, fasta_genome::String, gff_genome::String; 
    utr_length=150, stop_early=-1, skip_trimming=false, skip_aligning=false, skip_interactions=false, export_single_csv=true, 
    skip_significant=false, skip_postprocess=false, bwa_bin="bwa", sam_bin="samtools", fastp_bin="fastp")

    isdir(joinpath(rilseq_folder, "demultiplexed")) || mkpath(joinpath(rilseq_folder, "demultiplexed"))
    isdir(joinpath(rilseq_folder, "trimmed")) || mkpath(joinpath(rilseq_folder, "trimmed"))
    isdir(joinpath(rilseq_folder, "pe_bams")) || mkpath(joinpath(rilseq_folder, "pe_bams"))
    isdir(joinpath(rilseq_folder, "se_bams")) || mkpath(joinpath(rilseq_folder, "se_bams"))
    isdir(joinpath(rilseq_folder, "results")) || mkpath(joinpath(rilseq_folder, "results"))
    isdir(joinpath(rilseq_folder, "reports")) || mkpath(joinpath(rilseq_folder, "reports"))
    
    skip_trimming || trim_fastp(rilseq_folder, lib_names, fastp_bin=fastp_bin)
    skip_aligning || align(rilseq_folder, fasta_genome, lib_names; rev_complement=true, se_miss=3, pe_miss=3, bwa_bin=bwa_bin, sam_bin=sam_bin)
    skip_interactions || all_interactions(rilseq_folder, fasta_genome, lib_names)
    skip_significant || significant_chimeras(rilseq_folder, gff_genome, lib_names; utr_length=utr_length)
    skip_postprocess || postprocess(rilseq_folder, gff_genome, lib_names, abspath(rilseq_folder, "combined_results.xlsx"); export_single_csv=export_single_csv, utr_length=150)
end

function run_rilseq_analysis()
    barcodes_cc = [
        "AATAATGT",
        "CAACACTT",
        "ATAATTCT",
        "GTCCATAT"
    ]
    barcodes_vc = [
        "CAAGTGAT",
        "CGACTTGG",
        "GCGAGTTG",
        "AAGACGGG"
    ]
    libnames_cc = [
        "CC1",
        "CC2",
        "CC3",
        "CC4"
    ]
    libnames_vc = [
        "VC1",
        "VC2",
        "VC3",
        "VC4"
    ]
    folder_vc = "/home/abc/Data/vibrio/rilseq/library_rilseq/"
    folder_cc = "/home/abc/Data/caulo/rilseq/"
    gff_cc = "/home/abc/Data/caulo/annotation/NC_011916.gff"
    gff_vc = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"
    fasta_cc = "/home/abc/Data/caulo/genome/GCF_000022005.1_ASM2200v1_genomic.fna"
    fasta_vc = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"

    #rilseq_analysis(libnames_vc, folder_vc, fasta_vc, gff_vc; skip_trimming=true, skip_aligning=true, skip_interactions=false, skip_significant=false, skip_postprocess=false)
    rilseq_analysis(libnames_cc, folder_cc, fasta_cc, gff_cc; skip_trimming=true, skip_aligning=true, skip_interactions=true, skip_significant=false, skip_postprocess=false)
end

run_rilseq_analysis() 

#a = get_annotations("/home/abc/Data/caulo/annotation/NC_011916.gff")
#println(a)
#println([sum(a[chr][!, :typ] .== "srna") for chr in keys(a)])