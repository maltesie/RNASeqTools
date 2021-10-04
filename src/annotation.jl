mutable struct Annotation <: AnnotationStyle
    type::String
    name::String
    params::Dict{String, String}
end

function Annotation()
    Annotation("", "", Dict{String, String}())
end

function Annotation(type::String, name::String)
    Annotation(type, name, Dict{String, String}())
end

Base.isempty(annotation::Annotation) = isempty(annotation.type) && isempty(annotation.name)

mutable struct AlignmentAnnotation <: AnnotationStyle
    type::String
    name::String
    overlap::UInt8
end

function AlignmentAnnotation()
    AlignmentAnnotation("", "", 0)
end

Base.isempty(annotation::AlignmentAnnotation) = isempty(annotation.type) && isempty(annotation.name) && (annotation.overlap==0)

name(annot::T) where {T<:AnnotationStyle} = annot.name
type(annot::T) where {T<:AnnotationStyle} = annot.type

struct Features{T} <: AnnotationContainer
    list::IntervalCollection{T}
end

function Features()
    return Features([])
end

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true))
end

function Features(gff_file::String, type::Vector{String}; name_key="Name", fallback_key=nothing, same_name_rule=:all)
    @assert same_name_rule in (:first, :all, :none)
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    names = Dict{Tuple{String,String},Int}()
    for feature in features
        (GFF3.featuretype(feature) in type || isempty(type)) || continue
        seqn = GFF3.seqid(feature)
        name = ("NA", "NA")
        try
            name = (GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","))
        catch ex
            if ex isa KeyError
                try
                    isnothing(fallback_key) || (name = (GFF3.featuretype(feature), join(GFF3.attributes(feature, fallback_key), ",")))
                catch _
                    continue
                end
            else
              rethrow(ex)
            end
        end
        name in keys(names) ? (names[name]+=1) : (names[name]=1)
        (same_name_rule === :first && names[name] > 1) && continue
        names[name] > 1 ? (n = name[2] * "$(names[name])") : (n = name[2])
        annot = Annotation(GFF3.featuretype(feature), n, Dict(pair[1] => join(pair[2], ",") for pair in GFF3.attributes(feature)))
        push!(intervals, Interval(seqn, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annot))
    end
    same_name_rule === :none && (return Features([i for i in intervals if ((type(i), name(i)) in names && names[(type(i), name(i))] == 1)]))
    return Features(intervals)
end

function Features(gff_file::String, type::String; name_key="Name", fallback_key=nothing)
    return Features(gff_file, [type], name_key=name_key, fallback_key=fallback_key)
end

function Features(gff_file::String; name_key="Name", fallback_key=nothing)
    return Features(gff_file, String[], name_key=name_key, fallback_key=fallback_key)
end

function Features(coverage::Coverage, type::String)
    return Features([Interval(refname(i), leftposition(i), rightposition(i), strand(i), Annotation(type, "", Dict{String,String}("Value"=>"$(value(i))"))) for i in coverage])
end

type(features::Features) = Set(type(f) for f in features)
refnames(features::Features) = collect(keys(features.list.trees))
function overlaps(alignmentinterval::Interval{AlignmentAnnotation}, feature::Interval{Annotation})
    seqname(feature) == seqname(alignmentinterval) || return false
    strand(feature) == strand(alignmentinterval) || return false
    return leftposition(feature) <= rightposition(alignmentinterval) && leftposition(alignmentinterval) <= rightposition(feature)
end

Base.push!(features::Features, interval::Interval) = push!(features.list, interval)
function Base.merge!(features1::Features, features2::Features)
    for feature in features2
        push!(features1, feature)
    end
end
function Base.merge(features1::Features, features2::Features)
    re = Features()
    for feature in features1
        push!(re, feature)
    end
    for feature in features2
        push!(re, feature)
    end
    return re
end
Base.:*(featuresa::Features, featuresb::Features) = merge(featuresa, featuresb)

Base.iterate(features::Features) = iterate(features.list)
Base.iterate(features::Features, state::Tuple{Int64,GenomicFeatures.ICTree{Annotation},GenomicFeatures.ICTreeIteratorState{Annotation}}) = iterate(features.list, state)
Base.length(features::Features) = length(features.list)
Base.split(features::Features) = [Features([feature for feature in features if type(feature)==t], [t]) for t in types(features)]

Base.convert(::Type{Interval{Float64}}, i::Interval{T}) where T<:AnnotationStyle = Interval(refname(i), leftposition(i), rightposition(i), strand(i), 0.0)
strand_filter(a::Interval, b::Interval)::Bool = strand(a) == strand(b)
function GenomicFeatures.eachoverlap(features::I, feature::Interval{T}) where {I<:AnnotationContainer,T}
    t = T
    if I === Coverage && T === Annotation
        feature = Interval{Float64}(feature)
        t = Float64
    end
    haskey(features.list.trees, refname(feature)) ?
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), t}(strand_filter,
                GenomicFeatures.ICTreeIntersection{t}(), features.list.trees[refname(feature)], feature)) :
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), t}(strand_filter,
                GenomicFeatures.ICTreeIntersection{t}(), GenomicFeatures.ICTree{t}(), feature))
end

function hasoverlap(features::Features, feature::Interval)
    for _ in eachoverlap(features, feature)
        return true
    end
    return false
end
function firstoverlap(features::Features, feature::Interval)
    for int in eachoverlap(features, feature)
        return int
    end
    return nothing
end

function paramstring(params::Dict{String,String};priority=("Name", "Count"))
    ps = join(("$key=$(params[key])" for key in priority if key in keys(params)), ";")
    os = join(("$key=$value" for (key,value) in params if !(key in priority)), ";")
    return ps * ((isempty(ps) || isempty(os)) ? "" : ";") * os
end

function Base.write(file::String, features::Features)
    writer = GFF3.Writer(open(file, "w"))
    for feature in features
        line = "$(feature.seqname)\t.\t$(type(feature))\t$(feature.first)\t$(feature.last)\t.\t$(feature.strand)\t.\t$(paramstring(params(feature)))"
        record = GFF3.Record(line)
        write(writer, record)
    end
    close(writer)
end

function maxsignalposition(leftpos::Int, rightpos::Int, chr::String, st::Strand, coverage_dict::Dict{String,Coverage}, modify::Symbol)
    maxsignal::Float64 = 0.0
    pos::Int = -1
    l = "guess"
    for (lib, coverage) in coverage_dict
        for signal in eachoverlap(coverage, Interval(chr, leftpos, rightpos, st, Annotation()))
            value(signal) > maxsignal && (maxsignal=value(signal); pos=rightposition(signal); l=lib)
        end

    end
    if maxsignal != 0.0
        return modify === :left ? (pos, rightpos, l) : (leftpos, pos, l)
    else
        return (leftpos, rightpos, l)
    end
end

function addutrs!(features::Features; tss_positions::Union{Dict{String,Coverage},Nothing}=nothing, term_positions::Union{Dict{String,Coverage},Nothing}=nothing,
                    cds_type="CDS", five_type="5UTR", three_type="3UTR", max_utr_length=250, min_utr_length=25, guess_missing=true)

    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (type(feature)==cds_type) && (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (type(feature)==cds_type) && (feature.strand == STRAND_NEG)]
    first_feature, last_feature = base_features_pos[1], base_features_pos[end]
    stop, start = leftposition(first_feature), rightposition(last_feature)
    push!(new_features, Interval(refname(last_feature), start+1, start+max_utr_length, STRAND_POS, Annotation(three_type, name(last_feature), params(last_feature))))
    push!(new_features, Interval(refname(first_feature), max(1, stop-max_utr_length), stop-1, STRAND_POS, Annotation(five_type, name(first_feature), params(first_feature))))
    first_feature, last_feature = base_features_neg[1], base_features_neg[end]
    stop, start = leftposition(first_feature), rightposition(last_feature)
    push!(new_features, Interval(refname(last_feature), start+1, start+max_utr_length, STRAND_NEG, Annotation(five_type, name(last_feature), params(last_feature))))
    push!(new_features, Interval(refname(first_feature), max(1, stop-max_utr_length), stop-1, STRAND_NEG, Annotation(three_type, name(first_feature), params(first_feature))))

    for base_features in (base_features_pos, base_features_neg)
        nb_features = length(base_features)
        stran = base_features === base_features_pos ? STRAND_POS : STRAND_NEG
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            threeref, fiveref, threename, fivename = refname(feature), refname(next_feature), name(feature), name(next_feature)
            stop, start = leftposition(next_feature), rightposition(feature)
            threestart::Int, threestop::Int, fivestart::Int, fivestop::Int = 0, 0, 0, 0
            if refname(feature) != refname(next_feature)
                threestart, threestop = start+1, start+max_utr_length
                fivestart, fivestop = max(1, stop-max_utr_length), stop-1
            elseif  stop - start > 2 * max_utr_length + 1
                threestart, threestop = start+1, start+max_utr_length
                fivestart, fivestop = stop-max_utr_length, stop-1
            elseif stop - start > 2 * min_utr_length
                new_utr_length = floor(Int, (stop-start)/2)
                threestart, threestop = start+1, start+new_utr_length
                fivestart, fivestop = stop-new_utr_length, stop-1
            else
                continue
            end

            stran === STRAND_NEG && ((threestart, fivestart, threestop, fivestop, threename, threeref) = (fivestart, threestart, fivestop, threestop, fivename, fiveref))

            isnothing(tss_positions) ? (found_five_in="guess") :
            ((fivestart, fivestop, found_five_in) = maxsignalposition(fivestart, fivestop, fiveref, stran, tss_positions, stran === STRAND_POS ? :left : :right))
            isnothing(term_positions) ? (found_three_in="guess") :
            ((threestart, threestop, found_three_in) = maxsignalposition(threestart, threestop, threeref, stran, term_positions, stran === STRAND_POS ? :right : :left))
            (guess_missing || fixed_five) && push!(new_features, Interval(fiveref, fivestart, fivestop, stran,
                Annotation(five_type, fivename, merge(Dict("source"=>found_five_in), copy(stran === STRAND_NEG ? params(feature) : params(next_feature))))))
            (guess_missing || fixed_three) && push!(new_features, Interval(threeref, threestart, threestop, stran,
                Annotation(three_type, threename, merge(Dict("source"=>found_three_in), copy(stran === STRAND_NEG ? params(next_feature) : params(feature))))))
        end
    end
    for feature in new_features
        push!(features, feature)
    end
end

function addigrs!(features::Features; igr_type="IGR", min_igr_length=50)
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (feature.strand == STRAND_NEG)]
    nb_features_pos = length(base_features_pos)
    nb_features_neg = length(base_features_neg)
    for base_features in (base_features_pos, base_features_neg)
        nb_features = base_features === base_features_pos ? nb_features_pos : nb_features_neg
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            refname(feature) == refname(next_feature) || continue
            stop, start = leftposition(next_feature), rightposition(feature)
            (stop-1) - (start + 1) > min_igr_length || continue
            igr = Interval(refname(feature), start+1, stop-1, base_features === base_features_pos ? STRAND_POS : STRAND_NEG,
                                        Annotation(igr_type, name(feature)*":"*name(next_feature),
                                        Dict(key=>param(feature, key)*":"*param(next_feature, key) for key in keys(params(feature)) if key in keys(params(next_feature)))))
            push!(new_features, igr)
        end
    end

    for feature in new_features
        push!(features, feature)
    end
end

type(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.type
name(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.name
refname(feature::Interval{T}) where T<:AnnotationStyle = feature.seqname
params(feature::Interval{Annotation}) = feature.metadata.params
param(feature::Interval{Annotation}, key::String) = feature.metadata.params[key]
param(feature::Interval{Annotation}, key::String, ::Type{I}) where {I} = parse(I, feature.metadata.params[key])
hasannotationkey(feature::Interval{Annotation}, key::String) = key in keys(params(feature))
annotation(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata

typenamekey(feature::Interval{Annotation}) = type(feature) * ":" * name(feature)
function featureseqs(features::Features, genome::Genome; key_gen=typenamekey)
    seqs = Dict{String, LongDNASeq}()
    for feature in features
        seq = genome[refname(feature)][leftposition(feature):rightposition(feature)]
        strand(feature) == STRAND_NEG && reverse_complement!(seq)
        push!(seqs, key_gen(feature) => seq)
    end
    return Sequences{String}(seqs)
end

function annotate!(features::Features, from_reps::Vector{Coverage}, to_reps::Vector{Coverage})
    averages = zeros(Float64, length(features), length(from_reps)+length(to_reps))
    ps = zeros(Float64, length(features))
    fc = zeros(Float64, length(features))
    stop_from = length(from_reps)
    start_to = stop_from + 1
    stop_to = stop_from + length(to_reps)
    averages[:,1:stop_from] = normalizedcount(features, from_reps)
    averages[:,start_to:stop_to] = normalizedcount(features, to_reps)
    noise = rand!(similar(averages)) .* 0.0000001
    averages += noise
    avf = mean(@view(averages[:,1:stop_from]), dims=2)
    avt = mean(@view(averages[:,start_to:stop_to]), dims=2)
    for i in 1:length(features)
        t = UnequalVarianceTTest(@view(averages[i,1:stop_from]), @view(averages[i,start_to:stop_to]))
        ps[i] = pvalue(t)
        isnan(ps[i]) && (ps[i] = 1.0)
        fc[i] = log2(mean(@view(averages[i,start_to:stop_to]))/mean(@view(averages[i,1:stop_from])))
        isnan(fc[i]) && (fc[i] = 0.0)
    end
    adjps = adjust(PValues(ps), BenjaminiHochberg())
    for (feature, p, adj_p, fold_change, average_from, average_to) in zip(features, ps, adjps, fc, avf, avt)
        params(feature)["LogFoldChange"] = "$(round(fold_change; digits=3))"
        params(feature)["FoldChange"] = "$(round(2^fold_change; digits=3))"
        params(feature)["PValue"] = "$p"
        params(feature)["AdjustedPValue"] = "$adj_p"
        params(feature)["BaseValueFrom"] = "$average_from"
        params(feature)["BaseValueTo"] = "$average_to"
    end
end

function annotate!(features::Features, samples::Vector{Coverage}; count_key="Count")
    vals = [values(coverage) for coverage in samples]
    averages = zeros(Float64, length(features), length(samples))
    for (i,feature) in enumerate(features)
        for (j,rep) in enumerate(vals)
            strand_vals = strand(feature) === STRAND_NEG ? last(rep[refname(feature)]) : first(rep[refname(feature)])
            averages[i,j] = sum(strand_vals[leftposition(feature):rightposition(feature)])
        end
    end
    for (feature, vals) in zip(features, eachrow(averages))
        for (i,val) in enumerate(vals)
            params(feature)["$count_key$i"] = "$val"
        end
    end
end

function covratio(features::Features, coverage::Coverage)
    vals = values(coverage)
    total = 0.0
    for val in values(vals)
        total += sum(first(val)) + sum(last(val))
    end
    s = 0.0
    for feature in features
        picker = strand(feature) === STRAND_NEG ? last : first
        rightposition(feature) > length(picker(vals[refname(feature)])) && continue
        s += sum(picker(vals[refname(feature)])[leftposition(feature):rightposition(feature)])
    end
    return s/total
end

function normalizedcount(features::Features, samples::Vector{Coverage})
    vals = [values(coverage) for coverage in samples]
    averages = zeros(Float64, length(features), length(samples))
    for (i,feature) in enumerate(features)
        for (j,rep) in enumerate(vals)
            strand_vals = strand(feature) === STRAND_NEG ? last(rep[refname(feature)]) : first(rep[refname(feature)])
            averages[i,j] = mean(strand_vals[leftposition(feature):rightposition(feature)])
        end
    end
    avg_sample::Vector{Float64} = [geomean(averages[i, :]) for i in 1:length(features)] .+ 0.000000001
    norm_factors = [median(averages[:, i] ./ avg_sample) for i in 1:length(samples)]
    averages ./= norm_factors'
    return averages
end

function correlation(features::Features, coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    averages = normalizedcount(features, collect(coverages))
    correlations = cor(averages, dims=1)
    return correlations
end
correlation(features::Features, coverages::Vector{Coverage}) = correlation(features, coverages...)

function mincorrelation(features::Features, coverages::Coverage ...)
    corr = correlation(features, coverages)
    min_corr = 1.0
    for i in first(size(corr)), j in 1:i
        min_corr = min(matrix[i,j], min_corr)
    end
    return min_corr
end
mincorrelation(features::Features, coverages::Vector{Coverage}) = mincorrelation(features, coverages...)

function asdataframe(features::Features; add_keys=:all)
    add_keys === :none && (add_keys = Set())
    add_keys === :all && (add_keys = Set(key for feature in features for key in keys(params(feature))))
    df = DataFrame(name=repeat([""], length(features)), refname=repeat([""], length(features)), type=repeat([""], length(features)),
                    left=repeat([-1], length(features)), right=repeat([-1], length(features)), strand=repeat(['*'], length(features)))
    for key in add_keys
        df[!, Symbol(key)] = repeat([""], length(features))
    end
    for (i,feature) in enumerate(features)
        df[i, :name] = name(feature)
        df[i, :refname] = refname(feature)
        df[i, :type] = type(feature)
        df[i, :left] = leftposition(feature)
        df[i, :right] = rightposition(feature)
        df[i, :strand] = strand(feature) === STRAND_NEG ? '-' : '+'
        pa = params(feature)
        for key in keys(pa)
            (key in add_keys || add_keys === :all) && (df[i, Symbol(key)] = pa[key])
        end
    end
    return df
end

Base.write(fname::String, df::DataFrame) = CSV.write(fname, df)

function Base.show(features::Features)
    chrs = refnames(features)
    ts = type(features)
    stats = Dict(chr=>Dict(t=>0 for t in ts) for chr in chrs)
    nb_pos = 0
    nb_neg = 0
    for f in features
        stats[refname(f)][type(f)] += 1
        nb_pos += strand(f) === STRAND_POS
        nb_neg += strand(f) === STRAND_NEG
    end
    dt = Int(floor(maximum(length(t) for t in ts)/8))+1
    printstring = "\n$(nb_pos+nb_neg) features in total with $nb_pos on + strand and $nb_neg on - strand:\n\n"
    printstring *= "type$(repeat("\t",dt))$(join([chr[1:min(8,length(chr))]*repeat(" ", 9-min(8,length(chr))) for chr in chrs], "\t"))\n\n"
    for t in ts
        printstring *= "$(t)$(repeat("\t",dt-Int(floor(length(t)/8))))$(join([stats[chr][t] for chr in chrs], "\t\t"))\n"
    end
    println(printstring)
end
