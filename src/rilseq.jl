using BioSequences
using FASTX
using CodecZlib
using XAM
using Combinatorics
using HypothesisTests
using DataFrames
using CSV
using XLSX

function preprocess(input_files::Array{Array{String,1},1}, project_folders::Array{String,1}, 
    barcodes::Array{String,1}, names::Array{String,1}; stop_early::Int=-1)

    barcode_length::UInt8 = length(barcodes[1])+1
    dplxr = Demultiplexer(LongDNASeq.(barcodes), n_max_errors=1, distance=:hamming)

    stats_dict::Dict{String,Array{Int,1}} = Dict(folder=>zeros(Int,12) for folder in project_folders)
    record_1::FASTQ.Record = FASTQ.Record()
    record_2::FASTQ.Record = FASTQ.Record()

    for ((reads_fastq_file_1, reads_fastq_file_2), folder) in zip(input_files, project_folders)
        reader_1 = FASTQ.Reader(GzipDecompressorStream(open(reads_fastq_file_1, "r")))
        reader_2 = FASTQ.Reader(GzipDecompressorStream(open(reads_fastq_file_2, "r")))
        
        out_files = [[abspath(folder, "demultiplexed", "$(name)_1.fastq.gz"), 
                    abspath(folder, "demultiplexed", "$(name)_2.fastq.gz")] 
                    for name in names]
        
        writers = [[FASTQ.Writer(GzipCompressorStream(open(out1, "w"), level=2)),
                    FASTQ.Writer(GzipCompressorStream(open(out2, "w"), level=2))] 
                    for (out1, out2) in out_files] 
        sleep(0.01)
        
        reporter = open(abspath(folder, "reports", "demultiplexing.txt"), "w")
        c = 0
        while !eof(reader_1)
            read!(reader_1, record_1)
            read!(reader_2, record_2)
            ((c >= stop_early) & (stop_early > 0)) && break
            c += 1

            read_1 = LongDNASeq(FASTQ.sequence(record_1))
            (library_id, nb_errors) = demultiplex(dplxr, read_1)
            nb_errors == -1 && continue
            stats_dict[folder][library_id] += 1
            #println(sequence(record_1), " $library_id")
            record_1.sequence = (record_1.sequence.start + barcode_length):record_1.sequence.stop
            record_1.quality = (record_1.quality.start + barcode_length):record_1.quality.stop
            #record_2.sequence = record_2.sequence.start:(record_2.sequence.stop-barcode_length)
            #record_2.quality = record_2.quality.start:(record_2.quality.stop-barcode_length)
            #println(sequence(record_1), "\n")
            write(writers[library_id][1], record_1)
            write(writers[library_id][2], record_2)
        end
        
        for (writer_1, writer_2) in writers
            close(writer_1)
            close(writer_2)
        end
        
        close(reader_1)
        close(reader_2)

        count_string = join(["$(names[i]) - $(stats_dict[folder][i])\n" for i in 1:12])
        count_string *= "not identifyable - $(c-sum(stats_dict[folder]))\n"
        write(reporter, "Counted $c entries in $reads_fastq_file_1\n$count_string\n")
        close(reporter)
    end
end

function trim_fastp(project_folders::Array{String,1}, names::Array{String,1})

    record = FASTQ.Record()
    
    for folder in project_folders
    
        in_files = [[abspath(folder, "demultiplexed", "$(name)_1.fastq.gz"), 
                    abspath(folder, "demultiplexed", "$(name)_2.fastq.gz")] 
                    for name in names]
        out_files = [[abspath(folder, "trimmed", "$(name)_1.fastq.gz"), 
                      abspath(folder, "trimmed", "$(name)_2.fastq.gz")] 
                      for name in names]
        
        reporter = open(abspath(folder, "reports", "quality.txt"), "w")
        
        c = 0
        for ((file1, file2), (out1, out2), name) in zip(in_files, out_files, names)
            
            cmd = `./bin/fastp -a AGATCGGAAGAGC -i $file1 -I $file2 -o $out1 -O $out2 -M 25 -l 25 --trim_poly_g 10 --cut_front --cut_tail`# --trim_poly_g 10 -r --cut_right_window_size=6`
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

function align(project_folders::Array{String,1}, fasta_genome::String, 
                names::Array{String,1}; rev_complement=false, index_genome=true, 
                se_miss=2, pe_miss=0, se_length=25)

    for folder in project_folders
        
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
            
            cmd = pipeline(`./bin/bwa index -a is $fasta_genome`, stdout=nothing)
            index_genome && run(cmd)
            cmd = pipeline(`./bin/bwa aln -n $pe_miss -t 8 -R 200 $fasta_genome $file1`, stdout="tmp1.sai")
            run(cmd)
            cmd = pipeline(`./bin/bwa aln -n $pe_miss -t 8 -R 200 $fasta_genome $file2`, stdout="tmp2.sai")
            run(cmd)
            cmd = pipeline(`./bin/bwa sampe -a 1500 -P $fasta_genome tmp1.sai tmp2.sai $file1 $file2`, 
                stdout="tmp.bwa")
            run(cmd)
            cmd = pipeline(`./bin/samtools view -u tmp.bwa`, stdout="tmp.view")
            run(cmd)
            cmd = pipeline(`./bin/samtools sort tmp.view -o $pe_bam`, stdout=nothing)
            run(cmd)
            cmd = pipeline(`./bin/samtools index $pe_bam`, stdout=nothing)
            run(cmd)
            
            cut_fastq_reads(pe_bam, "tmp1.fastq", "tmp2.fastq";
                rev_complement=rev_complement, cut_len=se_length)
            
            cmd = pipeline(`./bin/bwa aln -n $se_miss -t 6 -R 200 $fasta_genome tmp1.fastq`, stdout="tmp1.sai")
            run(cmd)
            cmd = pipeline(`./bin/bwa samse $fasta_genome tmp1.sai tmp1.fastq`, stdout="tmp.bwa")
            run(cmd)
            cmd = pipeline(`./bin/samtools view -u tmp.bwa`, stdout="tmp.view")
            run(cmd)
            cmd = pipeline(`./bin/samtools sort tmp.view -o $se_bam1`)
            run(cmd)
            cmd = pipeline(`./bin/samtools index $se_bam1`)
            run(cmd)

            cmd = pipeline(`./bin/bwa aln -n $se_miss -t 6 -R 200 $fasta_genome tmp2.fastq`, stdout="tmp2.sai")
            run(cmd)
            cmd = pipeline(`./bin/bwa samse $fasta_genome tmp2.sai tmp2.fastq`, stdout="tmp.bwa")
            run(cmd)
            cmd = pipeline(`./bin/samtools view -u tmp.bwa`, stdout="tmp.view")
            run(cmd)
            cmd = pipeline(`./bin/samtools sort tmp.view -o $se_bam2`)
            run(cmd)
            cmd = pipeline(`./bin/samtools index $se_bam2`)
            run(cmd)
            
            rm("tmp1.sai")
            rm("tmp2.sai")
            rm("tmp.bwa")
            rm("tmp.view")
            rm("tmp1.fastq")
            rm("tmp2.fastq")
            
        end 
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

function all_interactions(project_folders::Array{String}, fasta_genome::String, names::Array{String}; 
                        generate_bams=true, cut_len=25, read1_reversed=false, read2_reversed=false)
    
    read1_names::Array{String,1}, read1_chrs::Array{String,1}, read1_auxs::Array{String,1} = [], [], []
    read1_poss::Array{Int,1}, read1_nms::Array{Int,1} = [], []
    read2_names::Array{String,1}, read2_chrs::Array{String,1}, read2_auxs::Array{String,1} = [], [], []
    read2_poss::Array{Int,1}, read2_nms::Array{Int,1} = [], []
    
    single_fragment_set::Set{String} = Set()
    
    for folder in project_folders
        
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
            
            single_fragment_set = get_single_fragment_set(pe_file)
    
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

function get_annotations(gff_file::String)::Dict{String, DataFrame}
    reader = open(gff_file)
    gff = read(reader, String)
    results::Dict{String, DataFrame} = Dict()

    for line in split(gff, "\n")
        startswith(line, "#") | isempty(line) && continue
        chr, _, typ, start, stop, _, strand, _, aux = split(line, "\t")
        typ == "CDS" || continue
        if strand == "-"
            start_int = parse(Int, stop) * -1
            stop_int = parse(Int, start) * -1
        else
            start_int = parse(Int, start)
            stop_int = parse(Int, stop)
        end
        name = split(aux, ";")[1][6:end]
        row = DataFrame(name=name, start=start_int, stop=stop_int)
        chr in keys(results) ? append!(results[chr], copy(row)) : results[chr] = copy(row)
    end

    for (chr, data) in results
        sort!(data, :start)
        utr5 = DataFrame(name="U5." .* data.name, start=data.start .- 100, stop=copy(data.start))
        utr3 = DataFrame(name="U3." .* data.name, start=copy(data.stop), stop=data.stop .+ 100)
        igr = DataFrame(name="IG<" .* data.name[1:end-1] .* "|" .* data.name[2:end] .* ">", 
            start=data.stop[1:end-1] .+ 100, stop=data.start[2:end] .- 100)
        igr= igr[igr.start .< igr.stop, :]
        #igr = DataFrame(name=igr.name[index], start=igr.start[index], stop=igr.stop[index])
        append!(data, utr5)
        append!(data, utr3)
        append!(data, igr)
        sort!(data, :start)
    end
    
    return results
end

function get_name(chr::String, pos::Int, annotations::Dict{String, DataFrame})::String
    names = annotations[chr].name[annotations[chr].start .< pos .< annotations[chr].stop]
    isempty(names) ? (return "not found") : (return names[1])
end

function get_full_name(chr::String, pos_start::Int, pos_end::Int,
        annotations::Dict{String, DataFrame})::String
    name1::String = get_name(chr, pos_start, annotations)
    name2::String = get_name(chr, pos_end, annotations)
    (name1 == name2) ? (return name1) : (return "$name1 -> $name2")
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
    
    return results, chr_trans
end

function check_rrna(pos1::Int, pos2::Int, chr1::String, chr2::String, 
        rrna::Array{Int, 2}, chr_trans::Dict{String, Int})::Bool
    
    if chr1 in keys(chr_trans)
        any((chr_trans[chr1] .== rrna[:,1]) .& (pos1 .>= rrna[:,2]) .& (pos1 .<= rrna[:,3])) && (return true)
    end
    if chr2 in keys(chr_trans)
        any((chr_trans[chr2] .== rrna[:,1]) .& (pos2 .>= rrna[:,2]) .& (pos2 .<= rrna[:,3])) && (return true)
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

    while !eof(reader)
        line::Array{SubString{String},1} = split(readline(reader))
        #println(line)
        name::String = line[1] 
        typ::String = line[2]
        pos1::Int = parse(Int, line[3]) 
        pos2::Int = parse(Int, line[4]) 
        chr1::String = line[5] 
        chr2::String = line[6] 
        check_rrna(pos1, pos2, chr1, chr2, rrna, chr_trans) && (count_rrna += 1; continue)
        count_total += 1
        seg1::Int = (pos1÷seglen)*seglen
        seg1<0 && (seg1-=seglen)
        seg2::Int = (pos2÷seglen)*seglen
        seg2<0 && (seg2-=seglen)
        coord1 = Coord(seg1, chr1)
        coord2 = Coord(seg2, chr2)
        interaction = Interact(coord1, coord2)

        if (typ == "chimeric") 
            interaction in keys(region_interactions) ? 
            push!(region_interactions[interaction], pos1=>pos2) :
            region_interactions[interaction] = [pos1=>pos2]
        end

        coord1 in keys(count_reads1) ? 
        (count_reads1[coord1] += 1) : 
        push!(count_reads1, coord1=>1)
        coord2 in keys(count_reads2) ?
        (count_reads2[coord2] += 1) : 
        push!(count_reads2, coord2=>1)

    end
    close(reader)
    #println(count_rrna, " rRNA vs. ", count_total, " rest")
    sum_pairs::Dict{Interact, Int} = Dict()
    for interaction in keys(region_interactions)
        nb_interactions = length(region_interactions[interaction])
        (nb_interactions > minints) && (sum_pairs[interaction] = nb_interactions)
    end
    
    return region_interactions, sum_pairs, count_reads1, count_reads2, count_total
end


function significant_chimeras(project_folders::Array{String, 1}, gff_genome::String, names::Array{String, 1}; 
            seglen=100, minints=5, minodds=1.0)
    
    report = ""
    rrna::Array{Int, 2}, chr_trans::Dict{String,Int} = get_rrna(gff_genome)

    annotations = get_annotations(gff_genome)

    header = ["Name RNA1", "pos1", "Name RNA2", "pos2", "# of chimeric fragments", "library", "norm. odds ration", "p-value",
            "IP interactions RNA1", "IP interactions RNA2", "Total interactions RNA1", "Total interactions RNA2",
            "RNA1 chromosome", "RNA1 first start", "RNA1 last start",
            "RNA2 chromosome", "RNA2 first start", "RNA2 last start", "b", 
            "c", "d"]
    
    interaction_files =[[abspath(project_folders[1], "results", "$(name)_raw.txt.gz"),
                        abspath(project_folders[2], "results", "$(name)_raw.txt.gz")] for name in names]
    
    out_files =[abspath(project_folders[2], "results", "$(name).csv") for name in names]
    c=0

    for ((tot_interaction_file, ip_interaction_file), out_file) in zip(interaction_files, out_files)
        
        ip_interactions, ip_interaction_count, ip_reads1_count, ip_reads2_count, ip_total_count = 
        process_interactions(ip_interaction_file, rrna, chr_trans)
        tot_interactions, tot_interaction_count, tot_reads1_count, tot_reads2_count, tot_total_count = 
        process_interactions(tot_interaction_file, rrna, chr_trans)
        
        results = DataFrame(name1=String[], pos1s=Int[], name2=String[], pos2s=Int[], nb_fragments=Int[], libs=String[], 
            norm_odds=Float64[], p_value=Float64[], ip_int1=Float64[], ip_int2=Float64[],
            total_int1=Float64[], total_int2=Float64[],
            chromosome1=String[], first_start1=Int[], last_start1=Int[], chromosome2=String[], 
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
            (pv > 0.05) && continue
            
            name1 = get_name(chr1, pos1, annotations)
            
            name2 = get_name(chr2, pos2, annotations)
            
            lib = "tbd"
            
            ip_int1 = ip_reads1_count[coord1]/ip_total_count
            coord1 in keys(tot_reads1_count) ? 
                    (total_int1 = tot_reads1_count[coord1]/tot_total_count) :
                    (total_int1 = 0)
            ip_int2 = ip_reads2_count[coord2]/ip_total_count
            coord2 in keys(tot_reads2_count) ? 
                    (total_int2 = tot_reads2_count[coord2]/tot_total_count) :
                    (total_int2 = 0)
            starts1 = sort(ip_interactions[interaction], by=x->x[1])
            first_start1, last_start1 = starts1[1][1], starts1[end][1]
            starts2 = sort(ip_interactions[interaction], by=x->x[2])
            first_start2, last_start2 = starts2[1][2], starts2[end][2]

            row = DataFrame(name1=name1, pos1s=pos1, name2=name2, pos2s=pos2, nb_fragments=interaction_count, 
                libs=lib, norm_odds=norm_odds, p_value=pv, ip_int1=ip_int1, total_int1=total_int1, 
                ip_int2=ip_int2, total_int2=total_int2, 
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

function merged_table(table::CSV.File, annotations::Dict{String,DataFrame})::DataFrame
    poss::MyMatrix = [[row[:pos1], row[:pos2]] for row in table]
    already_found::Set{Int} = Set() 
    merged_data = DataFrame(name1=String[], name2=String[], nb_fragments=Int[], 
        libs=String[], norm_odds=Float64[], p_value=Float64[], ip_total_ratio1=Float64[], 
        ip_total_ratio2=Float64[], chromosome1=String[], first_start1=Int[], last_start1=Int[], 
        chromosome2=String[], first_start2=Int[], last_start2=Int[])
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
        tot1 = sum([r[Symbol("Total interactions RNA1")] for r in table[indices]])
        tot2 = sum([r[Symbol("Total interactions RNA2")] for r in table[indices]])
        rat1, rat2 = round(ips1/tot1, digits=2), round(ips2/tot2, digits=2)
            
        n1 = get_full_name(chr1, first1, last1, annotations)
        n2 = get_full_name(chr2, first2, last2, annotations)
    
        b = sum([r[:b] for r in table[indices]])
        c = sum([r[:c] for r in table[indices]])
        d = ([r[:d] for r in table[indices]][1] + 
            [r[:c] for r in table[indices]][1] +
            [r[:b] for r in table[indices]][1] +
            [r[Symbol("# of chimeric fragments")] for r in table[indices]][1]) - (nb + b + c)
        
        f = FisherExactTest(nb,b,c,d)
        norm_odds = round(f.ω, digits=2)
        pv = round(pvalue(f), digits=10)
        
        row = DataFrame(name1=n1, name2=n2, nb_fragments=nb, 
        libs=lib, norm_odds=norm_odds, p_value=pv, ip_total_ratio1=rat1, 
        ip_total_ratio2=rat2, chromosome1=chr1, first_start1=first1, last_start1=last1, 
        chromosome2=chr2, first_start2=first2, last_start2=last2)
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

function unified_table(table1::CSV.File, table2::CSV.File, annotations::Dict{String,DataFrame})::DataFrame
    
    merged_data1 = merged_table(table1, annotations)
    merged_data2 = merged_table(table2, annotations)
    range1 = [[merged_data1[!,:first_start1][i], merged_data1[!,:last_start1][i],
            merged_data1[!,:first_start2][i], merged_data1[!,:last_start2][i]] for i in 1:length(merged_data1[!,:name1])]
    range2 = [[merged_data2[!,:first_start1][i], merged_data2[!,:last_start1][i],
            merged_data2[!,:first_start2][i], merged_data2[!,:last_start2][i]] for i in 1:length(merged_data2[!,:name1])]
    
    for row in eachrow(merged_data1)
        !isempty(range2) ? (id2 = match_id(row[:first_start1], row[:last_start1], row[:first_start2], row[:last_start2], range2)) : (id2=-1)
        if (id2 >= 0)
            row[:libs] = "1,2"
            row[:nb_fragments] += merged_data2[!,:nb_fragments][id2]
            row[:first_start1] = min(row[:first_start1], merged_data2[!,:first_start1][id2])
            row[:last_start1] = max(row[:last_start1], merged_data2[!,:last_start1][id2])
            row[:first_start2] = min(row[:first_start2], merged_data2[!,:first_start2][id2])
            row[:last_start2] = max(row[:last_start2], merged_data2[!,:last_start2][id2])
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

function cytoscape_json(data::DataFrame, fname::String)
    json = []
    nodes::Set{String} = Set()
    for row in eachrow(data)
        (row.name1 in nodes) || (push!(nodes, row.name1); push!(json, Dict(
            "data"=>Dict(
                "id" => row.name1, 
                "idInt" => length(nodes),
                "name" => row.name1,
                "score" => 1.0
                ),
            "position" => Dict(),
            "group" => "nodes",
            "removed" => false,
            "selected" => false,
            "selectable" => true,
            "locked" => false,
            "grabbed" => false,
            "grabbable" => true
            )))
        (row.name2 in nodes) || (push!(nodes, row.name2); push!(json, Dict(
            "data"=>Dict(
                "id" => row.name2, 
                "idInt" => length(nodes),
                "name" => row.name2,
                "score" => 1.0
                ),
            "position" => Dict(),
            "group" => "nodes",
            "removed" => false,
            "selected" => false,
            "selectable" => true,
            "locked" => false,
            "grabbed" => false,
            "grabbable" => true
            )))
        push!(json, Dict(
            "data"=>Dict(
                "source" => row.name1, 
                "target" => row.name2,
                "id" => row.name1 * "+" * row.name2,
                "weight" => row.nb_fragments
                ),
            "position" => Dict(),
            "group" => "edges",
            "removed" => false,
            "selected" => false,
            "selectable" => true,
            "locked" => false,
            "grabbed" => false,
            "grabbable" => true
            ))
    end

    open(fname,"w") do f
        JSON.print(f, json)
    end
end

function postprocess(folder::String, gff_genome::String, names::Array{String, 1}; export_json=true)
     
    table_files =[[abspath(folder, "results", "$(names[i]).csv"),
                    abspath(folder, "results", "$(names[i+1]).csv")] for i in 1:2:length(names)-1]
    
    unified_sheet_names = ["$(join(split(names[i], "_")[[1,3]], "_"))" for i in 1:2:length(names)-1]
    
    annotations = get_annotations(gff_genome)
    XLSX.openxlsx(abspath("combined_results.xlsx"), mode="w") do xf
        for (i, ((file_rep1, file_rep2), sheet_name)) in enumerate(zip(table_files, unified_sheet_names))
            table1, table2 = CSV.File(file_rep1), CSV.File(file_rep2)
            unified_tab = unified_table(table1, table2, annotations)
            export_json && cytoscape_json(unified_tab, sheet_name * ".json")
            (XLSX.sheetcount(xf) < i) &&  XLSX.addsheet!(xf, "new")
            sheet = xf[i]
            XLSX.rename!(sheet, sheet_name)
            XLSX.writetable!(sheet, [c for c in eachcol(unified_tab)], DataFrames.names(unified_tab))
        end
    end
end

function rilseq_analysis(lib_names::Array{String,1}, barcodes::Array{String,1}, rilseq_reads1::String, rilseq_reads2::String, 
    total_rna_reads1::String, total_rna_reads2::String, rilseq_folder::String, total_rna_folder::String, fasta_genome::String,
    gff_genome::String; stop_early=-1, skip_preprocessing=false, skip_trimming=false, skip_aligning=false, skip_interactions=false,
    export_json=true)

    input_files = [[total_rna_reads1, total_rna_reads2], [rilseq_reads1, rilseq_reads2]]
    project_folders = [total_rna_folder, rilseq_folder]

    for folder in project_folders
        isdir("$(folder)demultiplexed") || mkpath("$(folder)demultiplexed")
        isdir("$(folder)trimmed") || mkpath("$(folder)trimmed")
        isdir("$(folder)pe_bams") || mkpath("$(folder)pe_bams")
        isdir("$(folder)se_bams") || mkpath("$(folder)se_bams")
        isdir("$(folder)results") || mkpath("$(folder)results")
        isdir("$(folder)reports") || mkpath("$(folder)reports")
    end

    skip_preprocessing || preprocess(input_files, project_folders, barcodes, lib_names, stop_early=stop_early)
    skip_trimming || trim_fastp(project_folders, lib_names)
    skip_aligning || align(project_folders, fasta_genome, lib_names; rev_complement=true, se_miss=1, pe_miss=3)
    skip_interactions || all_interactions(project_folders, fasta_genome, lib_names)
    skip_interactions || significant_chimeras(project_folders, gff_genome, lib_names)
    postprocess(rilseq_folder, gff_genome, lib_names; export_json=export_json)
end