module RNASeqTools

using XAM
using FASTX
using BioAlignments
using CSV
using XLSX
using DataFrames
using BioSequences
using CodecZlib

export align_mem, align_backtrack, align_local
export trim_fastp, split_libs 
export read_bam, areconcordant, single_fragment_set
export rilseq_analysis

include("types.jl")
include("align.jl")
include("analysis.jl")
include("preprocess.jl")
include("io.jl")
include("utils.jl")
include("rilseq.jl")

end # module
