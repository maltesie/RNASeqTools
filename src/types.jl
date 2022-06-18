abstract type SequenceContainer end

abstract type AnnotationStyle end

abstract type AnnotationContainer end

abstract type AlignmentContainer end

abstract type FileCollection end

abstract type InteractionContainer end

abstract type CountContainer end

LongDNASeqPair = Tuple{LongDNASeq, LongDNASeq}
CoverageValues = Dict{String, Vector{Float64}}