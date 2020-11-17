module RNASeqTools

export align_mem, align_backtrack 
export trim_fastp, split_libs 
export read_bam, areconcordant, single_fragment_set

include("align.jl")
include("analysis.jl")
include("preprocess.jl")
include("io.jl")
include("utils.jl")

end # module
