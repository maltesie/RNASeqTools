using RNASeqTools

function write_utrs_fasta(annotations::Dict{String,DataFrame}, genome::Dict{String,String}, threeUTR_fasta::String, fiveUTR_fasta::String)
    record = FASTA.Record()
    five_writer = FASTA.Writer(GzipCompressorStream(open(fiveUTR_fasta, "w")))
    three_writer = FASTA.Writer(GzipCompressorStream(open(threeUTR_fasta, "w")))
    for (chr, sequence) in genome
        for (i,row) in enumerate(eachrow(annotations[chr]))
            ((row[:threeType] == "max") || (row[:fiveType] == "max")) || continue
            name = row["name"]
            if (row[:threeType] == "max") 
                (row["start"] < 0) ? 
                threeUTR = reverse_complement(LongDNASeq(sequence[-row["threeUTR"]:-row["stop"]])) : 
                threeUTR = sequence[row["stop"]:row["threeUTR"]]
                if length(threeUTR) > 20
                    record = FASTA.Record("$(name)_threeUTR", threeUTR)
                    write(three_writer, record)
                end
            end
            if (row["fiveType"] == "max")
                (row["start"] < 0) ?
                fiveUTR = reverse_complement(LongDNASeq(sequence[-row["start"]:-row["fiveUTR"]])) :
                fiveUTR = sequence[row["fiveUTR"]:row["start"]]
                if length(fiveUTR) > 20
                    record = FASTA.Record("$(name)_fiveUTR", fiveUTR)
                    write(five_writer, record)
                end
            end
        end
    end
    close(five_writer)
    close(three_writer)
end

function align_mem2(sequence_fasta::String, genome_files::Vector{String}, out_folder::String; bwa_bin="bwa", sam_bin="samtools")
    for genome in genome_files
        out_file = joinpath(out_folder, join(split(basename(genome), ".")[1:end-1],".") * ".bam")
        align_mem(sequence_fasta, out_file, genome; bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
end

function align_local(sequence_fasta::String, genome_files::Vector{String}, out_file::String)
    occursin(".fasta", sequence_fasta) ? reads = FastaReads(sequence_fasta) : reads = FastqReads(sequence_fasta)
    scoremodel = AffineGapScoreModel(match=5, mismatch=-1, gap_open=-3, gap_extend=-3)
    alignment_table = DataFrame(name=String[], length=Int[], sequence=String[])
    index_dict::Dict{String, Int} = Dict()
    for genome_file in genome_files
        genome = Genome(genome_file)
        alignment_table[!, Symbol(genome.spec * "_seqaln")] = fill("", nrow(alignment_table))
        alignment_table[!, Symbol(genome.spec * "_refaln")] = fill("", nrow(alignment_table))
        alignment_table[!, Symbol(genome.spec * "_score")] = Vector{Union{Missing,Float64}}(fill(missing, nrow(alignment_table)))
        for (i, (name, seq)) in enumerate(reads.seqs)
            pairwise_result = local_alignment(genome.seq, seq, scoremodel)
            if hasalignment(pairwise_result)
                a = alignment(pairwise_result)
                (name in keys(index_dict)) ? index = index_dict[name] : (index = nrow(alignment_table)+1; push!(index_dict, name=>index))
                (index > nrow(alignment_table)) && append!(alignment_table, 
                                                            DataFrame(merge(Dict("name"=>name, "length"=>length(seq), "sequence"=>String(seq)), 
                                                            Dict((endswith(column_name, "score") ? column_name=>0.0 : column_name=>"") for column_name in names(alignment_table) 
                                                            if !(column_name in ["name", "length", "sequence"])))))
                a.a.aln.anchors[1] = AlignmentAnchor(0, a.a.aln.anchors[1].refpos-a.a.aln.anchors[1].seqpos, '0')
                skip_end = length(seq) - a.a.aln.anchors[end].seqpos
                a.a.aln.anchors[end] = AlignmentAnchor(a.a.aln.anchors[end].seqpos + skip_end, a.a.aln.anchors[end].refpos + skip_end, a.a.aln.anchors[end].op)
                pairs = collect(a)
                qa = join([p[1] for p in pairs])
                ra = join([p[2] for p in pairs])
                alignment_table[index, Symbol(genome.spec * "_seqaln")] = qa
                alignment_table[index, Symbol(genome.spec * "_refaln")] = ra
                alignment_table[index, Symbol(genome.spec * "_score")] = score(pairwise_result)/length(seq)
            end
        end
    end
    CSV.write(out_file, alignment_table)
end

