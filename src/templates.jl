function prepare_data(data_path::String, genome::Genome; files=FastqgzFiles)
    files = files(data_path)
    trimmed = trim_fastp(files)
    bams = align_mem(trimmed, genome)
    compute_coverage(bams)
end

function analyze_deg(coverage_files::PairedSingleTypeFiles, features::Features, conditions::Dict{String, UnitRange}, results_path::String; add_keys=["BaseValueFrom", "BaseValueTo", "LogFoldChange", "PValue", "AdjustedPValue"])
    coverages = [Coverage(file1, file2) for (file1, file2) in coverage_files]
    for ((name1, range1), (name2, range2)) in combinations(conditions, 2)
        annotate!(features, coverages[range1], coverages[range2])
        write(joinpath(results_path, name1 * "_vs_" * name2 * ".csv"), asdataframe(features, add_keys=add_keys))
    end
end

function raw_counts(bam_files::SingleTypeFiles, features::Features, conditions::Dict{String, UnitRange}, results_file::String; invert_strand=:none)
    for (name1, range1) in conditions
        for (i, j) in enumerate(range1)
            annotate!(features, bam_files[j]; invert_strand=invert_strand, count_key="$name1$i")
        end
    end
    write(results_file, asdataframe(features))
end

function feature_ratio(coverage_files::PairedSingleTypeFiles, features::Features, results_file::String)
    write(results_file, join(["$(basename(file1)[1:end-11])\t$(covratio(features, Coverage(file1, file2)))" for (file1, file2) in coverage_files], "\n"))
end