module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3
using BioAlignments, BioSequences, GenomicFeatures
using DataFrames
using Plots, Measures
using Dash, DashHtmlComponents, DashCoreComponents, DashTable

export align_mem, align_backtrack, align_local
export trim_fastp, split_libs 
export Genome, Reads, PairedReads, Alignments, PairedAlignments, SingleTypeFiles, PairedSingleTypeFiles, Features
export cut!, rev_comp!, approxoccursin, annotate!
export similarity
export dashboard, hist_length_distribution, hist_similarity, line_nucleotide_distribution

include("types.jl")
include("misc.jl")

include("io.jl")
include("preprocess.jl")
include("align.jl")
include("utils.jl")
include("annotate.jl")
include("analyse.jl")
include("visualize.jl")

end
