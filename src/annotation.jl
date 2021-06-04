struct Annotation <: AnnotationStyle
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

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true))
end

function Features(gff_file::String, type::Vector{String}; name_key="Name")
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    for feature in features
        (GFF3.featuretype(feature) in type || isempty(type)) || continue
        seqn = GFF3.seqid(feature)
        annot = Annotation(GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","), Dict(pair[1] => join(pair[2], ",") for pair in GFF3.attributes(feature)))
        push!(intervals, Interval(seqn, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annot))
    end
    return Features(intervals)
end

function Features(gff_file::String, type::String; name_key="Name")
    return Features(gff_file, [type], name_key=name_key)
end

function Features(gff_file::String; name_key="Name")
    return Features(gff_file, String[], name_key=name_key)
end

function Features(coverage::Coverage, type::String)
    return Features([Interval(refname(i), leftposition(i), rightposition(i), strand(i), Annotation(type, "", Dict{String,String}("Value"=>"$(value(i))"))) for i in coverage])
end

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

Base.iterate(features::Features) = iterate(features.list)
Base.iterate(features::Features, state::Tuple{Int64,GenomicFeatures.ICTree{Annotation},GenomicFeatures.ICTreeIteratorState{Annotation}}) = iterate(features.list, state)
Base.length(features::Features) = length(features.list)

strand_filter(a::Interval, b::Interval)::Bool = strand(a) == strand(b)
function GenomicFeatures.eachoverlap(features::I, feature::Interval{T}) where {I<:AnnotationContainer, T<:AnnotationStyle}
    haskey(features.list.trees, refname(feature)) ?
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter, 
                GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refname(feature)], feature)) :
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter, 
                GenomicFeatures.ICTreeIntersection{Annotation}(), GenomicFeatures.ICTree{Annotation}(), feature))
end

function hasoverlap(features::Features, feature::Interval)
    for int in eachoverlap(features, feature)
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
    return ps * (isempty(ps) ? "" : ";") * join(("$key=$value" for (key,value) in params if !(key in priority)), ";")
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

function maxsignalposition(coverage::Coverage, from::Int, to::Int, strand::Strand)
    maxsignal::Float32 = 0.0
    pos::Int = -1
    for term in eachoverlap(coverage, Interval(refname(feature), from, to, strand, Annotation()))
        value(term) > maxsignal && (maxsignal=value(term); pos=rightposition(term))
    end
    return maxsignal, pos
end

function addutrs!(features::Features, tss_coverage::Union{Coverage,Nothing}, term_coverage::Union{Coverage,Nothing}; cds_typ="CDS", max_utr_length=150, min_utr_length=25, guess_missing=true)
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (type(feature)==cds_typ) && (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (type(feature)==cds_typ) && (feature.strand == STRAND_NEG)]
    first_feature, last_feature = base_features_pos[1], base_features_pos[end]
    stop, start = leftposition(first_feature), rightposition(last_feature)
    push!(new_features, Interval(refname(last_feature), start+1, start+utr_length, STRAND_POS, Annotation("3UTR", name(last_feature), Dict{String,String}())))
    push!(new_features, Interval(refname(first_feature), max(1, stop-utr_length), stop-1, STRAND_POS, Annotation("5UTR", name(first_feature), Dict{String,String}())))
    first_feature, last_feature = base_features_neg[1], base_features_neg[end]
    stop, start = leftposition(first_feature), rightposition(last_feature)
    push!(new_features, Interval(refname(last_feature), start+1, start+utr_length, STRAND_NEG, Annotation("5UTR", name(last_feature), Dict{String,String}())))
    push!(new_features, Interval(refname(first_feature), max(1, stop-utr_length), stop-1, STRAND_NEG, Annotation("3UTR", name(first_feature), Dict{String,String}())))
    
    for base_features in (base_features_pos, base_features_neg)
        nb_features = length(base_features)
        stran = base_features === base_features_pos ? STRAND_POS : STRAND_NEG
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            threeref, fiveref, threename, fivename = refname(feature), refname(next_feature), name(feature), name(next_feature)
            stop, start = leftposition(next_feature), rightposition(feature)
            threestart::Int, threestop::Int, fivestart::Int, fivestop::Int, maxsignal::Float32  = 0, 0, 0, 0, 0.0
            if refname(feature) != refname(next_feature)
                threestart, threestop = start+1, start+utr_length
                fivestart, fivestop = max(1, stop-utr_length), stop-1
            elseif  stop - start > 2 * utr_length + 1 
                threestart, threestop = start+1, start+utr_length
                fivestart, fivestop = stop-utr_length, stop-1
            elseif stop - start > 2 * min_utr_length
                new_utr_length = floor(Int, (stop-start)/2)
                threestart, threestop = start+1, start+new_utr_length
                fivestart, fivestop = stop-new_utr_length, stop-1
            else
                continue
            end
            base_features === base_features_neg && ((threestart, fivestart, threestop, fivestop) = (fivestart, threestart, fivestop, threestop))
            isnothing(tss_coverage) || (maxsignal, fivestop = maxsignalposition(tss_coverage, fivestart, fivestop, STRAND_POS))
            isnothing(term_coverage) || (maxsignal, threestop = maxsignalposition(term_coverage, threestart, threestop, STRAND_POS))
            if guess_missing || max_signal != 0.0
                push!(new_features, Interval(threeref, threestart, threestop, stran, Annotation("3UTR", threename, Dict{String,String}())))
                push!(new_features, Interval(fiveref, fivestart, fivestop, stran, Annotation("5UTR", fivename, Dict{String,String}())))
            end
        end
    end
    for feature in new_features
        push!(features, feature)
    end
end

function addigrs!(features::Features; fiveutr_type="5UTR", threeutr_type="3UTR")
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (type(feature) in [fiveutr_type, threeutr_type]) && (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (type(feature) in [fiveutr_type, threeutr_type]) && (feature.strand == STRAND_NEG)]
    nb_features_pos = length(base_features_pos)
    nb_features_neg = length(base_features_neg)
    for base_features in (base_features_pos, base_features_neg)
        nb_features = base_features === base_features_pos ? nb_features_pos : nb_features_neg
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            
            base_features === base_features_pos ? 
            type(feature) == threeutr_type && type(next_feature) == fiveutr_type && refname(feature) == refname(next_feature) || continue :
            type(feature) == fiveutr_type && type(next_feature) == threeutr_type && refname(feature) == refname(next_feature) || continue
            
            stop, start = leftposition(next_feature), rightposition(feature)
            (start + 1) < stop || continue
            igr = Interval(refname(feature), start+1, stop-1, base_features === base_features_pos ? STRAND_POS : STRAND_NEG, Annotation("IGR", name(feature)*":"*name(next_feature), Dict{String,String}()))
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
    averages = zeros(Float32, length(features), length(from_reps)+length(to_reps))
    ps = zeros(Float32, length(features))
    fc = zeros(Float32, length(features))
    stop_from = length(from_reps)
    start_to = stop_from + 1
    stop_to = stop_from + length(to_reps)
    averages[:,1:stop_from] = normalizedcount(features, from_reps)
    averages[:,start_to:stop_to] = normalizedcount(features, to_reps)
    avf = mean(@view(averages[:,1:stop_from]), dims=2)
    avt = mean(@view(averages[:,start_to:stop_to]), dims=2)
    for i in 1:length(features)
        t = OneSampleTTest(@view(averages[i,1:stop_from]), @view(averages[i,start_to:stop_to]))
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
    averages = zeros(Float32, length(features), length(samples))
    for (i,feature) in enumerate(features)
        for (j,rep) in enumerate(vals)
            strand_vals = strand(feature) === STRAND_NEG ? last(rep[refname(feature)]) : first(rep[refname(feature)])
            averages[i,j] = mean(strand_vals[leftposition(feature):rightposition(feature)])
        end
    end
    avg_sample::Vector{Float32} = [geomean(averages[i, :]) for i in 1:length(features)] .+ 0.000001
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
    df = DataFrame(name=String[], refname=String[], type=String[], left=Int[], right=Int[], strand=Char[])
    add_row =  DataFrame(name="", refname="", type="", left=0, right=0, strand='+')
    for key in add_keys
        df[!, Symbol(key)] = String[]
        add_row[!, Symbol(key)] = [""]
    end
    for feature in features
        add_row[1, :name] = name(feature)
        add_row[1, :refname] = refname(feature)
        add_row[1, :type] = type(feature)
        add_row[1, :left] = leftposition(feature)
        add_row[1, :right] = rightposition(feature)
        add_row[1, :strand] = strand(feature) === STRAND_NEG ? '-' : '+'
        pa = params(feature)
        for key in keys(pa)
            (key in add_keys || add_keys === :all) && (add_row[1, Symbol(key)] = pa[key])
        end
        append!(df, add_row)
    end
    return df
end

Base.write(fname::String, df::DataFrame) = CSV.write(fname, df)