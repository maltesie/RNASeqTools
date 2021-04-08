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

type(feature::Interval{T}) where T <: AnnotationStyle = feature.metadata.type
name(feature::Interval{T}) where T <: AnnotationStyle = feature.metadata.name
refname(feature::Interval{T}) where T <: AnnotationStyle = feature.seqname
Base.isempty(annotation::AlignmentAnnotation) = isempty(annotation.type) && isempty(annotation.name) && (annotation.overlap==0)

struct Features <: AnnotationContainer
    list::IntervalCollection{Annotation}
    description::Union{String, Nothing}
end

function Features(gff_file::String, type::Vector{String}, name_key::String; description=nothing)
    features = open(collect, GFF3.Reader, gff_file)
    intervals = Vector{Interval{Annotation}}()
    for feature in features
        GFF3.featuretype(feature) in type || continue
        seqname = GFF3.seqid(feature)
        annotation = Annotation(GFF3.featuretype(feature), join(GFF3.attributes(feature, name_key), ","))
        push!(intervals, Interval(seqname, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annotation))
    end
    return Features(IntervalCollection(intervals, true), description)
end

function Features(gff_file::String, type::String, name_key::String; description=nothing)
    return Features(gff_file, [type], name_key; description=description)
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
        line = "$(feature.seqname)\t.\t$(type(feature))\t$(feature.first)\t$(feature.last)\t.\t$(feature.strand)\t.\tName=$(name(feature))"
        record = GFF3.Record(line)
        write(writer, record)
    end
    close(writer)
end

function addutrs!(features::Features, typ::String)
    new_features = Vector{Interval{Annotation}}()
    for feature in features
        type(feature) != typ && continue
        start, stop = leftposition(feature), rightposition(feature)
        if strand(feature) == STRAND_POS
            global fiveutr = Interval(refname(feature), start-150, start-1, STRAND_POS, Annotation("5UTR", name(feature)))
            global threeutr = Interval(refname(feature), stop+1, stop+150, STRAND_POS, Annotation("3UTR", name(feature)))
        else
            global fiveutr = Interval(refname(feature), stop+1, stop+150, STRAND_NEG, Annotation("5UTR", name(feature)))
            global threeutr = Interval(refname(feature), start-150, start-1, STRAND_NEG, Annotation("3UTR", name(feature)))
        end
        push!(new_features, fiveutr)
        push!(new_features, threeutr)
    end
    for feature in new_features
        push!(features, feature)
    end
end