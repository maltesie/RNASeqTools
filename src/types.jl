abstract type SequenceContainer end

abstract type AnnotationStyle end

abstract type AnnotationContainer end

abstract type AlignmentContainer end

abstract type FileCollection end

abstract type InteractionContainer end

abstract type CountContainer end

LongDNAPair = Tuple{LongSequence, LongSequence}
CoverageValues = Dict{String, Vector{Float64}}

struct Sequences{T<:Union{String, UInt}}
    seq::LongSequence
    seqnames::Vector{T}
    ranges::Vector{UnitRange{Int}}
end

struct SequenceLogo
    alphabet::Alphabet
    bits::Matrix{Float64}
end

struct Genome
    seq::LongSequence
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

struct Alignments{T<:Union{String, UInt}}
    chroms::Vector{Tuple{String,Int}}
    tempnames::Vector{T}
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
    ranges::Vector{UnitRange{Int}}
end

struct AlignedPart
    ref::Interval{AlignmentAnnotation}
    seq::UnitRange{Int}
    nms::UInt32
    read::Symbol
end

struct AlignedRead
    range::UnitRange{Int}
    alns::Alignments
end

struct GenomeComparison
    alns::Vector{PairwiseAlignment}
    fromto::Vector{Tuple{Interval,Interval}}
end

struct Features{T} <: AnnotationContainer
    list::IntervalCollection{T}
    chroms::Dict{String, Int}
end

struct Counts <: CountContainer
    conditions::Dict{String, Vector{Int}}
    values::Matrix{Float64}
end

struct FeatureCounts <: CountContainer
    conditions::Dict{String, Vector{Int}}
    values::Matrix{Float64}
    features::Features{Annotation}
end

struct Coverage <: AnnotationContainer
    list::IntervalCollection{Float64}
    chroms::Vector{Tuple{String, Int}}
end

struct BaseCoverage
    genome::Genome
    fcount::Dict{String, Dict{Symbol, Vector{Int}}}
    rcount::Dict{String, Dict{Symbol, Vector{Int}}}
end

struct BaseAnnotation <: AnnotationStyle
    type::String
    name::String
    ref::Vector{Int}
    a::Vector{Int}
    t::Vector{Int}
    g::Vector{Int}
    c::Vector{Int}
    gap::Vector{Int}
    ins::Vector{Int}
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