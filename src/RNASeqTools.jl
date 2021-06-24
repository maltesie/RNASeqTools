module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
using Plots, Measures, Statistics, HypothesisTests, MultipleTesting, DataFrames, StatsBase, MultivariateStats, Combinatorics, MetaGraphs, Random
using Dash, DashHtmlComponents, DashCoreComponents, DashTable
import IntervalTrees.IntervalBTree, IntervalTrees.IntervalBTreeIteratorState

export align_mem, local_alignment
export trim_fastp, split_libs 
export Genome, Sequences, PairedSequences, Alignments, PairedAlignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage, Annotation, AlignmentAnnotation
export FastaFiles, PairedFastaFiles, FastagzFiles, PairedFastagzFiles, FastqFiles, PaireFastqFiles, FastqgzFiles, PairedFastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles
export cut!, rev_comp!, rev_comp, approxoccursin, annotate!, featureseqs, conservedfeatures, asdataframe
export hasannotation, ischimeric, istriplet, refinterval, readinterval, annotation, hasannotation, name, type, overlap, count, alignments, refname, params, param
export values, addutrs!, addigrs!, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, mincorrelation, normalizedcount, covratio
export similarity, tsss, terms, nucleotidecount
export dashboard, lengthhist, similarityhist, nucleotidedist, expressionpca
export feature_ratio, raw_counts, de_genes, prepare_data

include("types.jl")
include("misc.jl")

include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("coverage.jl")
include("annotation.jl")
include("alignment.jl")
include("visualization.jl")
include("templates.jl")

end
