module RNASeqTools

using XAM, FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig, DelimitedFiles, JLD2, BGZFStreams
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
import Plots: plot, histogram, scatter
import DataFrames: DataFrame, sort, nrow, names, innerjoin
using Measures, Statistics, HypothesisTests, MultipleTesting, StatsBase, MultivariateStats, Combinatorics, Random
using LightGraphs, ElasticArrays
using Dash, IterTools

export align_mem, local_alignment, align_minimap, align_kraken2
export trim_fastp, split_libs, download_sra, split_reads
export Genome, Sequences, PairedSequences, Alignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage, Interactions, Annotation, AlignmentAnnotation, BaseAnnotation, BaseCoverage, Counts, GenomeComparison
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles, GraphFiles
export PairedFastaFiles, PairedFastagzFiles, PaireFastqFiles, PairedFastqgzFiles
export cut!, rev_comp!, rev_comp, approxoccursin, nucleotidecount, similarcount, approxcount, hassimilar, annotate!, featureseqs, conservedfeatures, asdataframe, transform, ispositivestrand, refpercentage, totalvalues, refvalues
export hasannotation, ischimeric, istriplet, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread, name, type, overlap, count, parts, refname, params, param, setparam, hastype, hasname, typein, namein, distanceonread
export values, addutrs!, addigrs!, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, mincorrelation, covratio
export similarity, tsss, terms, hasannotationkey, readid, summarize, missmatchcount, eachpair, isfirstread, testsignificance!, addrelativepositions!
export dashboard, lengthhist, similarityhist, nucleotidedist, expressionpca, similarity, diffexptable, tss_annotation
export checkinteractions, uniqueinteractions, mismatchfractions, ismulti, mismatchpositions
export feature_ratio, feature_count, de_genes, prepare_data, chimeric_alignments, remove_features, unmapped_reads, transcriptional_startsites, full_annotation, prepare_data, deseq2_R, direct_rna_pipeline
export ANNOTATION_VCH, GENOME_VCH

include("types.jl")
include("misc.jl")

include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("coverage.jl")
include("annotation.jl")
include("counts.jl")
include("alignment.jl")
include("chimeric.jl")
include("visualization.jl")
include("templates.jl")

end
