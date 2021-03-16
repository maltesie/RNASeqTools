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

function align_backtrack(reads::Reads, out_file::String, genome_file::String; max_miss=2, bwa_bin="bwa", sam_bin="samtools")
    tmp_file = joinpath(dirname(out_file), "temp.fasta")
    write(tmp_file, reads)
    align_backtrack(tmp_file, out_file, genome_file; max_miss=max_miss, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_file)
end

function align_mem(in_file::String, out_file::String, genome_file::String; 
    z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -v 1 -t 6 $genome_file $in_file`, stdout="tmp.bwa")
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
    z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -t 6 -v 1 $genome_file $in_file1 $in_file2`, stdout="tmp.bwa")
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

function align_mem(read_files::SingleTypeFiles, genome::Genome; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")
    tmp_genome = joinpath(dirname(files.list[1]), "tmp_genome.fa")
    write(tmp_genome, genome)
    for file in read_files
        out_file = file[1:end-length(read_files.type)] * ".bam"
        align_mem(file, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    end
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
end

function align_mem(reads::Reads, genome::Genome, out_file::String; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")
    tmp_reads = joinpath(dirname(out_file), "tmp_reads.fasta")
    tmp_genome = joinpath(dirname(out_file), "tmp_genome.fa")
    write(tmp_genome, genome)
    write(tmp_reads, reads)
    align_mem(tmp_reads, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_reads)
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
end

function align_mem(reads::PairedReads, genome::Genome, out_file::String; z_score=100, bwa_bin="bwa-mem2", sam_bin="samtools")
    tmp_reads1 = joinpath(dirname(out_file), "temp1.fasta")
    tmp_reads2 = joinpath(dirname(out_file), "temp2.fasta")
    tmp_genome = joinpath(dirname(out_file), "tmp_genome.fa")
    write(tmp_reads1, tmp_reads2, reads)
    write(tmp_genome, genome)
    align_mem(tmp_reads1, tmp_reads2, out_file, tmp_genome; z_score=z_score, bwa_bin=bwa_bin, sam_bin=sam_bin)
    rm(tmp_reads1)
    rm(tmp_reads2)
    rm(tmp_genome)
    for ending in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
        rm(tmp_genome * ending)
    end
end

function local_alignment(reference_sequence::LongDNASeq, query_sequence::LongDNASeq, scoremodel::AffineGapScoreModel)
    res = pairalign(LocalAlignment(), query_sequence, reference_sequence, scoremodel)
    return res
end