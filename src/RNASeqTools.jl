module RNASeqTools
import XAM: BAM
using FASTX, CodecZlib, GFF3, BigWig, DelimitedFiles, BGZFStreams
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
import DataFrames: DataFrame, sort, nrow, names, innerjoin
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, GLM, StatsBase
using ElasticArrays, IterTools

export align_mem, align_minimap, align_kraken2
export trim_fastp, split_libs, download_prefetch, download_fasterq
export Genome, Sequences, PairedSequences, Alignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage
export Interactions, Annotation, AlignmentAnnotation, BaseAnnotation, BaseCoverage, Counts, FeatureCounts, GenomeComparison
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles, GraphFiles
export PairedFastaFiles, PairedFastagzFiles, PaireFastqFiles, PairedFastqgzFiles
export cut!, approxoccursin, nucleotidecount, similarcount, approxcount, hassimilar, annotate!, featureseqs, asdataframe, transform, ispositivestrand
export hasannotation, annotatedcount, annotationcount, alignednucleotidescount, ischimeric, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread
export name, type, overlap, count, parts, refname, featureparams, featureparam, setfeatureparam, hastype, hasname, typein, namein, distanceonread
export values, hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, mincorrelation, covratio, conditionsdict, normalize!
export similarity, transcriptionalstartsites, terminationsites, hasannotationkey, readid, summarize, missmatchcount, eachpair, isfirstread
export preprocess_data, groupfiles, difference_table
export ANNOTATION_VCH, GENOME_VCH

include("types.jl")
include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("coverage.jl")
include("annotation.jl")
include("counts.jl")
include("alignment.jl")

const ANNOTATION_VCH = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff"
const GENOME_VCH = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"

end
