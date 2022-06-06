struct Counts
    conditions::Dict{String, UnitRange{Int}}
    values::Matrix{Float64}
    features::Features{Annotation}
end

function Counts(features::Features, samples::Vector{Coverage}, conditions::Dict{String, UnitRange{Int}}; aggregation=maximum, normalization_method=:none)
    normalization_method in (:none, :tpm, :tpkm, :tmm) || raise(AssertionError("No method implemented for $normalization_method"))
    c = coveragecount(features, samples; aggregation=aggregation)
    normalization_method in (:tmm, :tpm) && normalize!(c; normalization_method=normalization_method)
    normalization_method === :tpkm && normalize!(c, features)
    return Counts(conditions, c, features)
end

function Counts(features::Features, samples::PairedSingleTypeFiles, conditions::Dict{String, UnitRange{Int}}; aggregation=:maximum, normalization_method=:none)
    samples.type === ".bw" || throw(AssertionError("File type has to be .bw"))
    coverages = [Coverage(file_forward, file_reverse) for (file_forward, file_reverse) in samples]
    Counts(features, coverages, conditions; aggregation=aggregation, normalization_method=normalization_method)
end

function Counts(features::Features, samples::SingleTypeFiles, conditions::Dict{String, UnitRange{Int}};
        normalization_method=:none, include_secondary_alignments=true, include_alternative_alignments=false)
    normalization_method in (:none, :tpm, :tpkm, :tmm) || raise(AssertionError("No method implemented for $normalization_method"))
    samples.type === ".bam" || throw(AssertionError("File type has to be .bam"))
    counts = zeros(Float64, length(features), sum(length(v) for v in values(conditions)))
    feature_trans = Dict{String, Int}(name(feature)*type(feature)=>i for (i, feature) in enumerate(features))
    for range in values(conditions)
        mybams = samples[range]
        for (i, bam_file) in enumerate(mybams)
            alignments = Alignments(bam_file; include_secondary_alignments=include_secondary_alignments, include_alternative_alignments=include_alternative_alignments)
            annotate!(alignments, features)
            for alignment in alignments
                for part in alignment
                    hasannotation(part) || continue
                    counts[feature_trans[name(part)*type(part)], range[i]] += 1.0
                end
            end
        end
    end
    normalization_method in (:tmm, :tpm) && normalize!(counts; normalization_method=normalization_method)
    normalization_method === :tpkm && normalize!(counts, features)
    Counts(conditions, counts, features)
end

function getindex(counts::Counts, condition::String)
    return counts.values[!, counts.conditions[condition]]
end

conditionrange(counts::Counts, condition::String) = counts.conditions[condition]

function coveragecount(features::Features, samples::Vector{Coverage}; aggregation=maximum)
    vals = [values(coverage) for coverage in samples]
    averages = zeros(Float64, length(features), length(samples))
    for (j,rep) in enumerate(vals)
        for (i,feature) in enumerate(features)
            strand_vals = strand(feature) === STRAND_NEG ? last(rep[refname(feature)]) : first(rep[refname(feature)])
            averages[i,j] = aggregation(strand_vals[leftposition(feature):rightposition(feature)])
        end
    end
    return averages
end

function normalize!(m::Matrix{Float64}; normalization_method=:tmm)
    if normalization_method === :tmm
        avg_sample::Vector{Float64} = [geomean(m[i, :]) for i in 1:first(size(m))]
        norm_factors = 2 .^ [median(log2.(m[:, i]) .- log2.(avg_sample)) for i in 1:last(size(m))]
        m ./= norm_factors'
    elseif normalization_method === :tpm
        m ./= (1000000 ./ sum(m; dims=1))'
    else
        raise(AssertionError("No method implemented for $normalization_method"))
    end
    return m
end

function normalize!(m::Matrix{Float64}, features::Features)
    normalize!(m; normalization_method=:tpm)
    for feature in features
        m[i, !] ./= (length(feature) / 1000)
    end
end

function dgetable(counts::Counts, control_condition::String, experiment_condition::String; method=:ttest)
    if method === :ttest
        avg_control = mean(@view(counts[control_condition]), dims=2)
        avg_experiment = mean(@view(counts[experiment_condition]), dims=2)
        ps = pvalue.(UnequalVarianceTTest.(@view(counts[control_condition]), @view(counts[experiment_condition])))
        fc = log2(avg_experiment) - log2(avg_control)
    elseif method === :glm
        throw(NotImplementedError)
    else
        throw(AssertionError("method must be one of :ttest, :glm"))
    end
    adjps = adjust(PValues(ps), BenjaminiHochberg())
    return DataFrame(
        feature=[name(feature) for feature in count.features],
        type=[type(feature) for feature in count.features],
        strand=[strand(feature) for feature in count.features],
        left=[leftposition(feature) for feature in count.features],
        right=[rightposition(feature) for feature in count.features],
        foldChange=round.(2 .^ fc; digits=3),
        log2FoldChange=round.(fc; digits=3),
        pValue=ps,
        fdr=adjps,
        averageControl=avg_control,
        averageExperiment=avg_experiment
    )
end

function correlation(counts::Counts)
    return cor(counts.values, dims=1)
end

function mincorrelation(counts::Counts)
    corr = correlation(counts)
    min_corr = 1.0
    for i in first(size(corr)), j in 1:i
        min_corr = min(matrix[i,j], min_corr)
    end
    return min_corr
end

function feature_count(features::Features, coverages::Vector{Coverage}, conditions::Dict{String, UnitRange{Int}}, results_path::String; between_conditions=nothing)
    expnames = Dict{String,Vector{String}}()
    for (name, range) in conditions
        annotate!(features, coverages[range]; count_key="$name")
        expnames[name] = ["$name$i" for i in 1:length(range)]
    end
    write(joinpath(results_path, "all_counts.csv"), asdataframe(features; add_keys=vcat([val for val in values(expnames)]...)))
    if !isnothing(between_conditions)
        for (cond1, cond2) in between_conditions
            exps = [expnames[cond1]...,expnames[cond2]...]
            write(joinpath(results_path, "$(cond1)_vs_$cond2.csv"), asdataframe(features; add_keys=exps))
        end
    end
end

function feature_count(features::Features, bams::SingleTypeFiles, conditions::Dict{String, UnitRange{Int}}, results_path::String; between_conditions=nothing, is_reverse_complement=false, only_unique_alignments=true)
    expnames = Dict{String,Vector{String}}()
    mybams = copy(bams)
    for (name, range) in conditions
        mybams = bams[range]
        annotate!(features, mybams; count_key="$name", is_reverse_complement=is_reverse_complement, only_unique_alignments=only_unique_alignments)
        expnames[name] = ["$name$i" for i in 1:length(range)]
        println("Finished counting in $name")
    end
    write(joinpath(results_path, "all_counts.csv"), asdataframe(features; add_keys=vcat([val for val in values(expnames)]...)))
    if !isnothing(between_conditions)
        for (cond1, cond2) in between_conditions
            exps = [expnames[cond1]...,expnames[cond2]...]
            write(joinpath(results_path, "$(cond1)_vs_$cond2.csv"), asdataframe(features; add_keys=exps))
        end
    end
end