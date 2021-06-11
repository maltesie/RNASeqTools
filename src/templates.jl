function prepare_data(data_path::String, genome::Genome; files=FastqgzFiles)
    files = files(data_path)
    trimmed = trim_fastp(files)
    bams = align_mem(trimmed, genome)
    compute_coverage(bams)
end

function de_genes(features::Features, coverage_files::PairedSingleTypeFiles, conditions::Dict{String, UnitRange{Int}}, results_path::String; add_keys=["BaseValueFrom", "BaseValueTo", "LogFoldChange", "PValue", "AdjustedPValue"])
    coverages = [Coverage(file1, file2) for (file1, file2) in coverage_files]
    for ((name1, range1), (name2, range2)) in combinations(conditions, 2)
        annotate!(features, coverages[range1], coverages[range2])
        write(joinpath(results_path, name1 * "_vs_" * name2 * ".csv"), asdataframe(features, add_keys=add_keys))
    end
end

function raw_counts(features::Features, bam_files::SingleTypeFiles, conditions::Dict{String, UnitRange{Int}}, results_file::String; invert_strand=:none)
    expnames = String[]
    for (name1, range1) in conditions
        for (i, j) in enumerate(range1)
            annotate!(features, bam_files[j]; invert_strand=invert_strand, count_key="$name1$i")
            push!(conditions, "$name1$i")
        end
    end
    write(results_file, asdataframe(features; add_keys=expnames))
end

function feature_ratio(features::Features, coverage_files::PairedSingleTypeFiles, results_file::String)
    write(results_file, join(["$(basename(file1)[1:end-11])\t$(covratio(features, Coverage(file1, file2)))" for (file1, file2) in coverage_files], "\n"))
end

function annotated_utrs(features::Features, results_file::String; tex=nothing, notex=nothing, term=nothing)
end

function conserved_features(features::Features, genome::Genome, targets::SingleTypeFiles, results_file::String)
end

function interaction_graph(features::Features, bams::SingleTypeFiles, results_path::String)
end