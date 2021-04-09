struct Annotation <: AnnotationStyle
    type::String
    name::String
end

function Annotation()
    Annotation("", "")
end

Base.isempty(annotation::Annotation) = isempty(annotation.type) && isempty(annotation.name)

struct AlignmentAnnotation <: AnnotationStyle
    type::String
    name::String
    overlap::UInt8
end

function AlignmentAnnotation()
    AlignmentAnnotation("", "", 0)
end

Base.isempty(annotation::AlignmentAnnotation) = isempty(annotation.type) && isempty(annotation.name) && (annotation.overlap==0)

struct Features <: AnnotationContainer
    list::IntervalCollection{Annotation}
end

function Features(gff_file::String, type::Vector{String}, name_key::String)
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    for feature in features
        GFF3.featuretype(feature) in type || continue
        seqname = GFF3.seqid(feature)
        annotation = Annotation(GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","))
        push!(intervals, Interval(seqname, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annotation))
    end
    return Features(IntervalCollection(intervals, true))
end

function Features(gff_file::String, type::String, name_key::String)
    return Features(gff_file, [type], name_key)
end

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true))
end

Base.push!(features::Features, interval::Interval) = push!(features.list, interval)
function merge!(features1::Features, features2::Features) 
    for feature in features2
        push!(features1, feature)
    end
end

Base.iterate(features::Features) = iterate(features.list)
Base.iterate(features::Features, state::Tuple{Int64,IntervalBTree{Int64,Interval{Annotation},64},IntervalBTreeIteratorState{Int64,Interval{Annotation},64}}) = iterate(features.list, state)
Base.length(features::Features) = length(features.list)

function Base.write(file::String, features::Features)
    writer = GFF3.Writer(open(file, "w"))
    for feature in features
        line = "$(feature.seqname)\t.\t$(annotationtype(feature))\t$(feature.first)\t$(feature.last)\t.\t$(feature.strand)\t.\tName=$(annotationname(feature))"
        record = GFF3.Record(line)
        write(writer, record)
    end
    close(writer)
end

function addutrs!(features::Features, typ::String; utr_length=150, min_utr_length=25)
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (annotationtype(feature)==typ) && (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (annotationtype(feature)==typ) && (feature.strand == STRAND_NEG)]
    nb_features_pos = length(base_features_pos)
    nb_features_neg = length(base_features_neg)
    for i in 1:nb_features_pos-1
        feature, next_feature = base_features_pos[i], base_features_pos[i+1]
        stop, start = leftposition(next_feature), rightposition(feature)
        if refname(feature) != refname(next_feature)
            push!(new_features, Interval(refname(feature), start+1, start+utr_length, STRAND_POS, Annotation("3UTR", annotationname(feature))))
            push!(new_features, Interval(refname(next_feature), max(1, stop-utr_length), stop-1, STRAND_POS, Annotation("5UTR", annotationname(next_feature))))
        elseif  stop - start > 2 * utr_length + 1 
            push!(new_features, Interval(refname(next_feature), stop-utr_length, stop-1, STRAND_POS, Annotation("5UTR", annotationname(next_feature))))
            push!(new_features, Interval(refname(feature), start+1, start+utr_length, STRAND_POS, Annotation("3UTR", annotationname(feature))))
            push!(new_features, Interval(refname(feature), start+utr_length+1, stop-utr_length-1, STRAND_POS, Annotation("IGR", annotationname(feature)*":"*annotationname(next_feature))))
        elseif stop - start > 2 * min_utr_length
            new_utr_length = floor(Int, (stop-start)/2)
            push!(new_features, Interval(refname(next_feature), stop-new_utr_length, stop-1, STRAND_POS, Annotation("5UTR", annotationname(feature))))
            push!(new_features, Interval(refname(feature), start+1, start+new_utr_length, STRAND_POS, Annotation("3UTR", annotationname(feature))))
        end
    end
    for i in 1:nb_features_neg-1
        feature, next_feature = base_features_neg[i], base_features_neg[i+1]
        stop, start = leftposition(next_feature), rightposition(feature)
        if refname(feature) != refname(next_feature) 
            push!(new_features, Interval(refname(feature), start+1, start+utr_length, STRAND_NEG, Annotation("5UTR", annotationname(feature))))
            push!(new_features, Interval(refname(next_feature), max(1, stop-utr_length), stop-1, STRAND_NEG, Annotation("3UTR", annotationname(next_feature))))
        elseif  stop - start > 2 * utr_length + 1
            push!(new_features, Interval(refname(next_feature), stop-utr_length, stop-1, STRAND_NEG, Annotation("3UTR", annotationname(next_feature))))
            push!(new_features, Interval(refname(feature), start+1, start+utr_length, STRAND_NEG, Annotation("3UTR", annotationname(feature))))
            push!(new_features, Interval(refname(feature), start+1, start+utr_length, STRAND_NEG, Annotation("IGR", annotationname(feature)*":"*annotationname(next_feature))))
        elseif stop - start > 2 * min_utr_length
            new_utr_length = floor(Int, (stop-start)/2)
            push!(new_features, Interval(refname(next_feature), stop-new_utr_length, stop-1, STRAND_NEG, Annotation("5UTR", annotationname(next_feature))))
            push!(new_features, Interval(refname(feature), start+1, start+new_utr_length, STRAND_NEG, Annotation("3UTR", annotationname(feature))))
        end
    end
    for feature in new_features
        push!(features, feature)
    end
end

annotationtype(feature::Interval{T}) where {T<:AnnotationStyle} = feature.metadata.type
annotationname(feature::Interval{T}) where {T<:AnnotationStyle} = feature.metadata.name
refname(feature::Interval{T}) where {T<:AnnotationStyle} = feature.seqname