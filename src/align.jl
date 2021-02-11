function align_backtrack(in_file::String, out_file::String, genome_file::String; 
    max_miss=2, bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin aln -n $max_miss -t 6 -R 500 $genome_file $in_file`, stdout="tmp.sai")
    run(cmd)
    cmd = pipeline(`$bwa_bin samse $genome_file tmp.sai $in_file`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.sai")
    rm("tmp.bwa")
    rm("tmp.view")
end

function align_backtrack(in_file1::String, in_file2::String, out_file::String, genome_file::String; 
    max_miss=2, bwa_bin="bwa", sam_bin="samtools")
        
    cmd = pipeline(`./bin/bwa index -a is $genome_file`, stdout=nothing)
    run(cmd)
    cmd = pipeline(`./bin/bwa aln -n $max_miss -t 6 -R 500 $genome_file $in_file1`, stdout="tmp1.sai")
    run(cmd)
    cmd = pipeline(`./bin/bwa aln -n $max_miss -t 6 -R 500 $genome_file $in_file2`, stdout="tmp2.sai")
    run(cmd)
    cmd = pipeline(`./bin/bwa sampe -a 1500 -P $genome_file tmp1.sai tmp2.sai $in_file1 $in_file2`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp1.sai")
    rm("tmp2.sai")
    rm("tmp.bwa")
    rm("tmp.view")
end

function align_mem(in_file::String, out_file::String, genome_file::String; 
    bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -t 6 $genome_file $in_file`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.bwa")
    rm("tmp.view")
end

function align_mem(in_file1::String, in_file2::String, out_file::String, genome_file::String; 
    bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -t 6 $genome_file $in_file1 $in_file2`, stdout="tmp.bwa")
    run(cmd)
    cmd = pipeline(`$sam_bin view -u tmp.bwa`, stdout="tmp.view")
    run(cmd)
    cmd = pipeline(`$sam_bin sort tmp.view -o $out_file`)
    run(cmd)
    cmd = pipeline(`$sam_bin index $out_file`)
    run(cmd)

    rm("tmp.bwa")
    rm("tmp.view")
end

function align_mem(sequence_fasta::String, genome_files::Vector{String}, out_folder::String; bwa_bin="bwa", sam_bin="samtools")
    for genome in genome_files
        out_file = joinpath(out_folder, join(split(basename(genome), ".")[1:end-1],".") * ".bam")
        align_mem(sequence_fasta, out_file, genome; bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
end

function local_alignment(reference_sequence::LongDNASeq, query_sequence::LongDNASeq, scoremodel::AffineGapScoreModel)
    res = pairalign(LocalAlignment(), query_sequence, reference_sequence, scoremodel)
    return res
end

function align_local(sequence_fasta::String, genome_files::Vector{String}, out_file::String)
    occursin(".fasta", sequence_fasta) ? reads = FastaReads(sequence_fasta) : reads = FastqReads(sequence_fasta)
    scoremodel = AffineGapScoreModel(match=5, mismatch=-1, gap_open=-3, gap_extend=-3)
    alignment_table = DataFrame(name=String[], start=Int[], stop=Int[], length=Int[], sequence=String[])
    for genome_file in genome_files
        genome = Genome(genome_file)
        alignment_table[Symbol(genome.spec * "_sequence")] = fill("", nrow(alignment_table))
        alignment_table[Symbol(genome.spec * "_score")]::Union{Missing, Int} = fill(missing, nrow(alignment_table))
        for (i, (name, seq)) in enumerate(reads.seqs)
            pairwise_result = local_alignment(genome.seq, seq, scoremodel)
            if hasalignment(pairwise_result)
                a = alignment(pairwise_result)
                pairs = collect(a)
                qa = join([p[1] for p in pairs])
                ra = join([p[1] for p in pairs])
            end
        end
    end
end
