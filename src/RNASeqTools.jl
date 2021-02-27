module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib
using BioAlignments, BioSequences, GenomicFeatures
using DataFrames
using Plots, Measures
using Dash, DashHtmlComponents, DashCoreComponents, DashTable

export align_mem, align_backtrack, align_local
export trim_fastp, split_libs 
export Genome, FastqReads, FastaReads, Alignments, PairedAlignments
export cut!
export tss, terms, annotate_utrs!
export similarity
export dashboard, hist_length_distribution, hist_similarity, line_nucleotide_distribution
export Reads, PairedReads

include("types.jl")
include("misc.jl")
include("preprocess.jl")
include("align.jl")

include("io.jl")
include("utils.jl")
include("annotate.jl")
include("analyse.jl")
include("visualize.jl")

end
