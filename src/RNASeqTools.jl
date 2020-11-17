module RNASeqTools

import XAM, BioAlignments

export align_mem, align_backtrack, trim_fastp, demultiplex, read_bam, areconcordant, single_fragment_set

inlude("align.jl")
include("analysis.jl")
include("preprocess.jl")
include("io.jl")

end # module
