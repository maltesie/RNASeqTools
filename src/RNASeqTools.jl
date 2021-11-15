module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig, DelimitedFiles, JLD2, BGZFStreams
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
import Plots: plot, histogram, scatter
import DataFrames: DataFrame, sort, nrow, names, eachcol
using Measures, Statistics, HypothesisTests, MultipleTesting, StatsBase, MultivariateStats, Combinatorics, Random
using LightGraphs
using Dash, IterTools

export align_mem, local_alignment
export trim_fastp, split_libs
export Genome, Sequences, PairedSequences, Alignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage, Interactions, Annotation, AlignmentAnnotation, ErrorCoverage, ErrorFeatures
export FastaFiles, PairedFastaFiles, FastagzFiles, PairedFastagzFiles, FastqFiles, PaireFastqFiles, FastqgzFiles, PairedFastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles
export cut!, rev_comp!, rev_comp, approxoccursin, annotate!, featureseqs, conservedfeatures, asdataframe, transform, ispositivestrand, refpercentage, totalvalues, refvalues
export hasannotation, ischimeric, istriplet, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread, name, type, overlap, count, parts, refname, params, param, setparam, hastype, hasname, typein, namein, distanceonread
export values, addutrs!, addigrs!, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, mincorrelation, normalizedcount, covratio
export similarity, tsss, terms, nucleotidecount, occurences, hasannotationkey
export dashboard, lengthhist, similarityhist, nucleotidedist, expressionpca, similarity
export feature_ratio, feature_count, de_genes, prepare_data, chimeric_alignments, remove_features, unmapped_reads, transcriptional_startsites, full_annotation, prepare_data
export ANNOTATION_VCH, GENOME_VCH


include("types.jl")
include("misc.jl")

include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("coverage.jl")
include("annotation.jl")
include("alignment.jl")
include("graph.jl")
include("visualization.jl")
include("templates.jl")

end
