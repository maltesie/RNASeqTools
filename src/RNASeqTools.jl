module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
using Plots, Measures, Statistics
using Dash, DashHtmlComponents, DashCoreComponents, DashTable
import IntervalTrees.IntervalBTree, IntervalTrees.IntervalBTreeIteratorState

export align_mem, local_alignment
export trim_fastp, split_libs 
export Genome, Reads, PairedReads, Alignments, PairedAlignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage, Annotation, AlignmentAnnotation
export cut!, rev_comp!, rev_comp, approxoccursin, annotate!
export hasannotation, ischimeric, istriplet, refinterval, readinterval, annotation, annotationname, annotationtype, annotationoverlap, count, alignments, refname
export values, addutrs!, addigrs!, hasoverlap, firstoverlap
export similarity
export dashboard, hist_length_distribution, hist_similarity, line_nucleotide_distribution

include("types.jl")
include("misc.jl")

include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("annotation.jl")
include("alignment.jl")
include("coverage.jl")
include("analysis.jl")
include("visualization.jl")

end
