module RNASeqTools

import XAM, BioSequences, FASTX, DataFrames

export align_mem, align_backtrack, trim_fastp, demultiplex, read_bam, areconcordant, single_fragment_set

include("align.jl")
include("analysis.jl")
include("preprocess.jl")
include("io.jl")

end # module
