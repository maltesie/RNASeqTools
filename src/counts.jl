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

function FeatureCounts(features::Features, samples::Vector{Coverage}; conditions=Dict("sample"=>collect(1:length(samples))), aggregation=maximum, normalization_method=:none)
    normalization_method in (:none, :tpm, :tpkm, :tmm) || raise(AssertionError("No method implemented for $normalization_method"))
    c = coveragecount(features, samples; aggregation=aggregation)
    normalization_method in (:tmm, :tpm) && normalize!(c; normalization_method=normalization_method)
    normalization_method === :tpkm && normalize!(c, features)
    return FeatureCounts(conditions, c, features)
end

function FeatureCounts(features::Features, samples::PairedSingleTypeFiles; conditions=groups(samples), aggregation=:maximum, normalization_method=:none)
    samples.type === ".bw" || throw(AssertionError("File type has to be .bw"))
    coverages = [Coverage(file_forward, file_reverse) for (file_forward, file_reverse) in samples]
    FeatureCounts(features, coverages; conditions=conditions, aggregation=aggregation, normalization_method=normalization_method)
end

function FeatureCounts(features::Features, samples::SingleTypeFiles; conditions=groups(samples),
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
    FeatureCounts(conditions, counts, features)
end

function Base.getindex(counts::T, condition::String) where {T<:CountContainer}
    return counts.values[:, counts.conditions[condition]]
end
Base.length(counts::T) where {T<:CountContainer} = size(counts.values)[1]

conditionrange(counts::T, condition::String) where {T<:CountContainer} = counts.conditions[condition]

function gmean(a::Array{Float64})
    s = 0.0
    n = length(a)
    for i = 1 : n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end

function normalize!(m::Matrix{Float64}; normalization_method=:rle, q=0.5)
    if normalization_method === :rle
        avg_sample::Vector{Float64} = [gmean(m[i, :]) for i in 1:first(size(m))]
        logdiffs = [log2.(m[:, i]) .- log2.(avg_sample) for i in 1:last(size(m))]
        norm_factors = 2 .^ [quantile(logdiff[(!).(isnan.(logdiff))], q) for logdiff in logdiffs]
        m ./= norm_factors'
    elseif normalization_method === :tpm
        m ./= (sum(m; dims=1) ./ 1000000)
    else
        raise(AssertionError("No method implemented for $normalization_method"))
    end
    return m
end
normalize!(counts::T; normalization_method=:rle) where {T<:CountContainer} =
    normalize!(counts.values; normalization_method=normalization_method)

function normalize!(m::Matrix{Float64}, features::Features; normalization_method=:tpkm)
    normalization_method != :tpkm && raise(AssertionError("No method implemented for $normalization_method"))
    normalize!(m; normalization_method=:tpm)
    for feature in features
        m[i, !] ./= (length(feature) / 1000)
    end
end
normalize!(counts::T, features::Features; normalization_method=:rpkm) where {T<:CountContainer} =
    normalize!(counts.values, features; normalization_method=normalization_method)

function cqn_offsets(m::Matrix{Float64})

end

function normalize!(m::Matrix{Float64}, features::Features, genome::Genome; normalization_method=:cqn)
    normalization_method != :cqn && raise(AssertionError("No method implemented for $normalization_method"))
end
normalize!(counts::T, features::Features, genome::Genome; normalization_method=:cqn) where {T<:CountContainer} =
    normalize!(counts.values, features, genome; normalization_method=normalization_method)

function difference_ttest(counts::T, control_condition::String, experiment_condition::String; within_sample=false, tail=:both) where {T<:CountContainer}
    avg_control = vec(mean(counts[control_condition], dims=2))
    avg_experiment = vec(mean(counts[experiment_condition], dims=2))
    m_control = counts[control_condition]
    m_experiment = counts[experiment_condition]
    m_control = within_sample ? [log.(m_control[i,:]) for i in 1:length(counts)] : [m_control[i,:] for i in 1:length(counts)]
    m_experiment = within_sample ? [log.(m_experiment[i,:]) for i in 1:length(counts)] : [m_experiment[i,:] for i in 1:length(counts)]
    ps = pvalue.(within_sample ? OneSampleTTest.(m_control, m_experiment) : UnequalVarianceTTest.(m_control, m_experiment) ; tail=tail)
    fc = log2.(avg_experiment) .- log2.(avg_control)
    padj = adjust(PValues(ps), BenjaminiHochberg())
    return (avg_control, fc, padj)
end

function difference_glm(counts::T, control_condition::String, experiment_condition::String; d=NegativeBinomial(), tail=:both) where {T<:CountContainer}
    control_range = counts.conditions[control_condition]
    exp_range = counts.conditions[experiment_condition]
    design_matrix = ones(Int, (length(control_range)+length(exp_range), 2))
    design_matrix[1:length(control_range), 2] .= 0
    data_matrix = Int.(floor.(hcat(counts[control_condition], counts[experiment_condition])))
    bases = Vector{Float64}(undef, length(counts))
    fcs = Vector{Float64}(undef, length(counts))
    ps = Vector{Float64}(undef, length(counts))
    compute_pvalue = tail === :both ? x -> 2 * ccdf(Normal(), abs(x)) :
                        tail === :right ? x -> ccdf(Normal(), x) :
                        tail === :left ? x -> cdf(Normal(), x) :
                        raise(AssertionError("tail can be :right, :left or :both"))
    l = LogLink()
    for (i, y) in enumerate(eachrow(data_matrix))
        try
            g = glm(design_matrix, y, d, l)
            base, fc = coef(g)
            _, fc_se = stderror(g)
            z = fc / fc_se
            p = compute_pvalue(z)
            isnan(p) && (p=1.0)
            bases[i] = base
            fcs[i] = fc / log(2)
            ps[i] = p
        catch e
            bases[i] = mean(y[1:length(control_range)])
            fcs[i] = log2(mean(y[end-length(exp_range):end]) / mean(y[1:length(control_range)]))
            ps[i] = 1.0
        end
    end
    println(sum(isnan(v) for v in ps))
    padj = adjust(PValues(ps), BenjaminiHochberg())
    return (exp.(bases), fcs, padj)
end

function differential_counts_table(counts::Counts, control_condition::String, experiment_condition::String; method=:glm)
    method in (:glm, :ttest) || throw(AssertionError("Method must be eather :ttest or :glm"))

    dge_function = method === :glm ? dge_glm : dge_ttest
    (base, fc, padj) = dge_function(counts, control_condition, experiment_condition)

    return DataFrame(
        feature=[name(feature) for feature in count.features],
        type=[type(feature) for feature in count.features],
        strand=[strand(feature) for feature in count.features],
        left=[leftposition(feature) for feature in count.features],
        right=[rightposition(feature) for feature in count.features],
        base=base,
        foldChange=round.(2 .^ fc; digits=3),
        log2FoldChange=round.(fc; digits=3),
        fdr=padj,
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