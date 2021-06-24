function prepare_data(data_path::String, genome::Genome; files=FastqgzFiles)
    files = files(data_path)
    trimmed = trim_fastp(files)
    bams = align_mem(trimmed, genome;)
    compute_coverage(bams)
end

function de_genes(features::Features, coverages::Vector{Coverage}, conditions::Dict{String, UnitRange{Int}}, results_path::String; between_conditions=nothing, add_keys=["BaseValueFrom", "BaseValueTo", "LogFoldChange", "PValue", "AdjustedPValue"])
    between_conditions = isnothing(between_conditions) ? combinations(collect(conditions), 2) : [((a,conditions[a]), (b,conditions[b])) for (a,b) in between_conditions]
    for ((name1, range1), (name2, range2)) in between_conditions
        annotate!(features, coverages[range1], coverages[range2])
        write(joinpath(results_path, name1 * "_vs_" * name2 * ".csv"), asdataframe(features, add_keys=add_keys))
    end
end

function raw_counts(features::Features, coverages::Vector{Coverage}, conditions::Dict{String, UnitRange{Int}}, results_file::String)
    expnames = String[]
    for (name, range) in conditions
        annotate!(features, coverages[range]; count_key="$name")
        append!(expnames, ["$name$i" for i in 1:length(range)])
    end
    write(results_file, asdataframe(features; add_keys=expnames))
end

function raw_counts(features::Features, bams::SingleTypeFiles, conditions::Dict{String, UnitRange{Int}}, results_file::String)
    expnames = String[]
    for (name, range) in conditions
        annotate!(features, bams[range]; count_key="$name")
        append!(expnames, ["$name$i" for i in 1:length(range)])
    end
    write(results_file, asdataframe(features; add_keys=expnames))
end

function feature_ratio(features::Features, coverage_files::PairedSingleTypeFiles, results_file::String)
    result_string = "filename\t" * join([t for t in features.types], "\t") * "\n"
    split_features = split(features)
    for (file1,file2) in coverage_files
        coverage = Coverage(file1,file2)
        result_string *= basename(file1)[1:end-11] * "\t" * join([covratio(f, coverage) for f in split_features], "\t") * "\n"
    end
    write(results_file, result_string)
end

function annotated_utrs(features::Features, results_file::String; tex=nothing, notex=nothing, term=nothing)
end

function conserved_features(features::Features, genome::Genome, targets::SingleTypeFiles, results_file::String)
end

function interaction_graph(features::Features, bams::SingleTypeFiles, results_path::String)
end