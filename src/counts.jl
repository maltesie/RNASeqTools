function coveragecount(features::Features{Annotation}, samples::Vector{Coverage}; aggregation=:max)
    aggregation in (:max, :mean) || throw(AssertionError("aggretation must be :max or :mean"))
    vals = [values(coverage) for coverage in samples]
    averages = zeros(Float64, length(features), length(samples))
    aggfn = aggregation == :max ? maximum : mean
    for (j,rep) in enumerate(vals)
        for (i,feature) in enumerate(features)
            strand_vals = strand(feature) === STRAND_NEG ? last(rep[refname(feature)]) : first(rep[refname(feature)])
            averages[i,j] = aggfn(abs.(strand_vals[leftposition(feature):rightposition(feature)]))
        end
    end
    return averages
end

function FeatureCounts(features::Features{Annotation}, samples::Vector{Coverage}; aggregation=:max, normalization_method=:none)
    normalization_method in (:none, :tpm, :tpkm, :tmm) || throw(AssertionError("No method implemented for $normalization_method"))
    c = coveragecount(features, samples; aggregation=aggregation)
    normalization_method in (:tmm, :tpm) && normalize!(c; normalization_method=normalization_method)
    normalization_method === :tpkm && normalize!(c, features)
    return FeatureCounts(features, c)
end

function FeatureCounts(features::Features{Annotation}, samples::PairedSingleTypeFiles; aggregation=:max, normalization_method=:none)
    samples.type === ".bw" || throw(AssertionError("Only .bw files are supported."))
    coverages = [Coverage(file_forward, file_reverse) for (file_forward, file_reverse) in samples]
    FeatureCounts(features, coverages; aggregation=aggregation, normalization_method=normalization_method)
end

function FeatureCounts(features::Features{Annotation}, samples::SingleTypeFiles; normalization_method=:none, include_secondary_alignments=true,
        include_alternative_alignments=false, is_reverse_complement=false)
    normalization_method in (:none, :tpm, :tpkm, :rle) || raise(AssertionError("No method implemented for $normalization_method"))
    samples.type === ".bam" || throw(AssertionError("Only .bam files are supported"))
    counts = zeros(Float64, length(features), length(samples))
    feature_trans = Dict{UInt, Int}(hash(name(feature), hash(type(feature)))=>i for (i, feature) in enumerate(features))
    for (i, bam_file) in enumerate(samples)
        alignments = AlignedReads(bam_file; include_secondary_alignments=include_secondary_alignments,
                                include_alternative_alignments=include_alternative_alignments, is_reverse_complement=is_reverse_complement)
        annotate!(alignments, features)
        for alignment in alignments
            for part in alignment
                hasannotation(part) || continue
                counts[feature_trans[hash(name(feature), hash(type(feature)))], range[i]] += 1.0
            end
        end
    end
    normalization_method in (:rle, :tpm) && normalize!(counts; normalization_method=normalization_method)
    normalization_method === :tpkm && normalize!(counts, features)
    FeatureCounts(features, counts)
end

Base.length(counts::T) where {T<:CountContainer} = size(counts.values)[1]
nsamples(counts::T) where {T<:CountContainer} = size(counts.values)[2]

summarize(counts::FeatureCounts) = "FeatureCounts with $(length(counts)) features and $(nsamples(counts)) samples."
Base.show(io::IO, counts::FeatureCounts) = print(io, summarize(counts))

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
    elseif normalization_method === :rpm
        m ./= (sum(m; dims=1) ./ 1000000)
    else
        raise(AssertionError("No method implemented for $normalization_method"))
    end
    return m
end
normalize!(counts::T; normalization_method=:rle) where {T<:CountContainer} =
    normalize!(counts.values; normalization_method=normalization_method)

function normalize!(m::Matrix{Float64}, features::Features{Annotation}; normalization_method=:rpkm)
    normalization_method != :tpkm && raise(AssertionError("No method implemented for $normalization_method"))
    normalize!(m; normalization_method=:rpm)
    for feature in features
        m[i, !] ./= (length(feature) / 1000)
    end
    return m
end
normalize!(counts::T, features::Features{Annotation}; normalization_method=:rpkm) where {T<:CountContainer} =
    normalize!(counts.values, features; normalization_method=normalization_method)

function normalize!(m::Matrix{Float64}, features::Features{Annotation}, genome::Genome; normalization_method=:cqn)
    normalization_method != :cqn && raise(AssertionError("No method implemented for $normalization_method"))
end
normalize!(counts::T, features::Features{Annotation}, genome::Genome; normalization_method=:cqn) where {T<:CountContainer} =
    normalize!(counts.values, features, genome; normalization_method=normalization_method)

function difference_glm(counts::T, ctrl_index::Vector{Int}, exp_index::Vector{Int}; d=NegativeBinomial(), tail=:both) where {T<:CountContainer}
    design_matrix = ones(Int, (length(ctrl_index)+length(exp_index), 2))
    design_matrix[1:length(ctrl_index), 2] .= 0
    data_matrix = Int.(round.(Int, hcat(counts.values[:, ctrl_index], counts.values[:, exp_index])))
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
            bases[i] = mean(y[1:length(ctrl_index)])
            fcs[i] = log2(mean(y[end-length(exp_index):end]) / mean(y[1:length(ctrl_index)]))
            ps[i] = 1.0
        end
    end
    padj = adjust(PValues(ps), BenjaminiHochberg())
    return (exp.(bases), fcs, padj)
end
difference_glm(counts::T, ctrl_index::UnitRange{Int}, exp_index::UnitRange{Int}; d=NegativeBinomial(), tail=:both) where {T<:CountContainer} =
    difference_glm(counts, collect(ctrl_index), collect(exp_index); d=d, tail=tail)

function correlation(counts::T) where T <: CountContainer
    return cor(counts.values, dims=1)
end