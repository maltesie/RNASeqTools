struct Counts{T<:Union{Float64, Int}}
    conditions::Dict{String, UnitRange{Int}}
    values::Matrix{T}
    features::Features{Annotation}
end

function Counts(features::Features, samples::Vector{Coverage}, conditions::Dict{String, UnitRange{Int}}; normalize_counts=true)
    c = rawcount(features, samples)
    normalize_counts && normalize!(c)
    return Counts(conditions, c, features)
end

function getindex(counts::Counts, condition::String)
    return counts.values[!, counts.conditions[condition]]
end

conditionrange(counts::Counts, condition::String) = counts.conditions[condition]

function rawcount(features::Features, samples::Vector{Coverage}; aggregation=maximum)
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

function normalize!(m::Matrix)
    avg_sample::Vector{Float64} = [geomean(m[i, :]) for i in 1:first(size(m))]
    norm_factors = 2 .^ [median(log2.(m[:, i]) .- log2.(avg_sample)) for i in 1:last(size(m))]
    m ./= norm_factors'
    return m
end

function diffexptable(counts::Counts, control_condition::String, experiment_condition::String; method=:ttest)

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