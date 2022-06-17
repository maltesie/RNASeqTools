module RNASeqTools

import XAM: BAM
using FASTX, CSV, XLSX, CodecZlib, GFF3, BigWig, DelimitedFiles, JLD2, BGZFStreams
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
import DataFrames: DataFrame, sort, nrow, names, innerjoin
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions
using LightGraphs, ElasticArrays, IterTools

export align_mem, align_minimap, align_kraken2
export trim_fastp, split_libs, download_sra
export Genome, Sequences, PairedSequences, Alignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage
export Interactions, Annotation, AlignmentAnnotation, BaseAnnotation, BaseCoverage, Counts, GenomeComparison
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles, GraphFiles
export PairedFastaFiles, PairedFastagzFiles, PaireFastqFiles, PairedFastqgzFiles
export cut!, approxoccursin, nucleotidecount, similarcount, approxcount, hassimilar, annotate!, featureseqs, asdataframe, transform, ispositivestrand
export hasannotation, ischimeric, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread
export name, type, overlap, count, parts, refname, params, param, setparam, hastype, hasname, typein, namein, distanceonread
export values, addutrs!, addigrs!, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, mincorrelation, covratio
export similarity, tsss, terms, hasannotationkey, readid, summarize, missmatchcount, eachpair, isfirstread, testsignificance!, addrelativepositions!
export checkinteractions, uniqueinteractions, mismatchfractions, ismulti, mismatchpositions, deletionpositions
export feature_ratio, feature_count, preprocess_data, chimeric_alignments, remove_features, unmapped_reads, full_annotation
export preprocess_data, deseq2_R, direct_rna_pipeline, tss_annotation
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
include("templates.jl")

end
