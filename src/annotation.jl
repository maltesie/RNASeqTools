struct Annotation <: AnnotationStyle
    type::String
    name::String
    params::Dict{String, String}
end

function Annotation()
    Annotation("", "", Dict{String, String}())
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

struct Features{T} <: AnnotationContainer
    list::IntervalCollection{T}
end

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true))
end

function Features(gff_file::String, type::Vector{String}, name_key::String)
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    for feature in features
        GFF3.featuretype(feature) in type || continue
        seqname = GFF3.seqid(feature)
        annotation = Annotation(GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","), Dict(pair[1] => join(pair[2], ",") for pair in GFF3.attributes(feature)))
        push!(intervals, Interval(seqname, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annotation))
    end
    return Features(intervals)
end

function Features(gff_file::String, type::String, name_key::String)
    return Features(gff_file, [type], name_key)
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
        nb_features_pos = length(base_features)
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
            base_features === base_features_neg && (threestart, fivestart, threestop, fivestop = fivestart, threestart, fivestop, threestop)
            isnothing(tss_coverage) || maxsignal, fivestop = maxsignalposition(tss_coverage, fivestart, fivestop, STRAND_POS)
            isnothing(term_coverage) || maxsignal, threestop = maxsignalposition(term_coverage, threestart, threestop, STRAND_POS)
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
        stan = base_features === base_features_pos ? STRAND_POS : STRAND_NEG
        for i in 1:nb_features_-1
            feature, next_feature = base_features[i], base_features[i+1]
            
            base_features === base_features_pos ? 
            type(feature) == threeutr_type && type(next_feature) == fiveutr_type && refname(feature) == refname(next_feature) || continue :
            type(feature) == fiveutr_type && type(next_feature) == threeutr_type && refname(feature) == refname(next_feature) || continue
            
            stop, start = leftposition(next_feature), rightposition(feature)
            (start + 1) < stop || continue
            igr = Interval(refname(feature), start+1, stop-1, stran, Annotation("IGR", name(feature)*":"*name(next_feature), Dict{String,String}()))
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

typenamekey(feature::Interval{Annotation}) = type(feature) * ":" * name(feature)
function featureseqs(features::Features, genome::Genome; key_gen=typenamekey)
    seqs = Dict{String, LongDNASeq}()
    for feature in features
        sequence = genome[refname(feature)][leftposition(feature):rightposition(feature)]
        strand(feature) == STRAND_NEG && reverse_complement!(sequence)
        push!(seqs, key_gen(feature) => sequence)
    end
    return Sequences{String}(seqs)
end