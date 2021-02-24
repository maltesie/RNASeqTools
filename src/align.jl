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
    z_score=100, bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -t 6 $genome_file $in_file`, stdout="tmp.bwa")
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
    z_score=100, bwa_bin="bwa", sam_bin="samtools")

    cmd = pipeline(`$bwa_bin index -a is $genome_file`)
    run(cmd)
    cmd = pipeline(`$bwa_bin mem -d $z_score -t 6 $genome_file $in_file1 $in_file2`, stdout="tmp.bwa")
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


