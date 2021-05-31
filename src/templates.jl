function prepare_data(data_path::String, genome_file::String; files=FastqgzFiles)
    genome = Genome(genome_file)
    files = files(data_path)
    trimmed = trim_fastp(files)
    align_mem(trimmed, genome)
    bams = BamFiles(data_path)
    compute_coverage(bams)
end

function analyze_deg(data_path::String, annotation_file::String)
end