module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
using Plots, Measures, Statistics
using Dash, DashHtmlComponents, DashCoreComponents, DashTable
import IntervalTrees.IntervalBTree, IntervalTrees.IntervalBTreeIteratorState

export align_mem, local_alignment
export trim_fastp, split_libs 
export Genome, Sequences, PairedSequences, Alignments, PairedAlignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage, Annotation, AlignmentAnnotation
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles
export cut!, rev_comp!, rev_comp, approxoccursin, annotate!, featureseqs, conservedfeatures
export hasannotation, ischimeric, istriplet, refinterval, readinterval, annotation, name, type, overlap, count, alignments, refname
export values, addutrs!, addigrs!, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation
export similarity, tsss, terms, nucleotidecount
export dashboard, lengthhist, similarityhist, nucleotidedist

include("types.jl")
include("misc.jl")

include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("annotation.jl")
include("alignment.jl")
include("coverage.jl")
include("visualization.jl")

end
