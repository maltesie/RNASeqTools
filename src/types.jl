abstract type SequenceContainer end

abstract type AnnotationStyle end

abstract type AnnotationContainer end

abstract type AlignmentContainer end

abstract type FileCollection end

abstract type CountContainer end

struct Sequences
    seq::LongDNA{4}
    tempnames::Vector{UInt}
    ranges::Vector{UnitRange{Int}}
end

struct Logo
    nseqs::Int
    weights::Matrix{Float64}
    alphabet::Vector{DNA}
end

struct Genome
    seq::LongDNA{4}
    chroms::Dict{String, UnitRange{Int}}
end

struct Annotation <: AnnotationStyle
    type::String
    name::String
    params::Dict{String, String}
end

struct AlignmentAnnotation <: AnnotationStyle
    type::String
    name::String
    overlap::UInt8
end

struct AlignedReads
    chroms::Vector{Tuple{String,Int}}
    tempnames::Vector{UInt}
    leftpos::Vector{Int}
    rightpos::Vector{Int}
    read_leftpos::Vector{Int}
    read_rightpos::Vector{Int}
    reads::Vector{Symbol}
    nms::Vector{UInt32}
    refnames::Vector{String}
    strands::Vector{Strand}
    annames::Vector{String}
    antypes::Vector{String}
    anols::Vector{UInt8}
    anleftrel::Vector{UInt8}
    anrightrel::Vector{UInt8}
    pindex::Vector{Int}
    annotated::Vector{Bool}
    ranges::Vector{UnitRange{Int}}
end

struct AlignedRead
    range::UnitRange{Int}
    alns::AlignedReads
end

struct AlignedInterval
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    nms::UInt32
    read::Symbol
end

struct Features{T} <: AnnotationContainer
    list::IntervalCollection{T}
    chroms::Dict{String, Int}
end

struct FeatureCounts <: CountContainer
    features::Features{Annotation}
    values::Matrix{Float64}
end

struct Coverage <: AnnotationContainer
    list::IntervalCollection{Float64}
    chroms::Vector{Tuple{String, Int}}
end

mutable struct SingleTypeFiles <: FileCollection
    list::Vector{String}
    type::String
end

mutable struct PairedSingleTypeFiles <: FileCollection
    list::Vector{Tuple{String,String}}
    type::String
    suffix1::String
    suffix2::String
end