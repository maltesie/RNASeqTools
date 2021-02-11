using BioSequences

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

function align_mem(sequence_fasta::String, out_folder::String, genome_files::Vector{String}; bwa_bin="bwa", sam_bin="samtools")
    for genome in genome_files
        out_file = joinpath(out_folder, join(split(basename(genome), ".")[1:end-1],".") * ".bam")
        align_mem(sequence_fasta, out_file, genome; bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
end

function local_alignment(reference_sequence::LongDNASeq, query_sequence::LongDNASeq)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
    res = pairalign(LocalAlignment(), reference_sequence, query_sequence, scoremodel)
    return res
end

function align_local(sequence_fasta::String, out_folder::String, genome_files::Vector{String})
    occursin(sequence_fasta, ".fasta") ? reads = FastaReads(sequence_fasta) : reads = FastqReads(sequence_fasta)
    for genome_file in genome_files
        genome = Genome(genome_file)
        already_found = Int[]
        for (chr, ref_seq) in genome.seqs
            for (i, (name, seq)) in enumerate(reads.seqs)
                (i in already_found) && continue
                @time pairwise_result = local_alignment(ref_seq, seq)
                hasalignment(pairwise_result) && print(alignment(pairwise_result))
            end
        end
    end
end
