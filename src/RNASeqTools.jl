module RNASeqTools
import XAM: BAM
using FASTX, CodecZlib, GFF3, BigWig, DelimitedFiles, BGZFStreams
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
import DataFrames: DataFrame, sort, nrow, names, innerjoin
using Statistics, HypothesisTests, MultipleTesting, Combinatorics, Random, Distributions, GLM, StatsBase
using ElasticArrays, IterTools

export align_mem, align_minimap, align_kraken2
export preprocess_data, trim_fastp, split_libs, download_prefetch, download_fasterq, split_paired_reads_file, transform
export Genome, Sequences, PairedSequences, Alignments, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage
export Annotation, AlignmentAnnotation, BaseAnnotation, BaseCoverage, Counts, FeatureCounts, GenomeComparison
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles
export PairedFastaFiles, PairedFastagzFiles, PaireFastqFiles, PairedFastqgzFiles
export cut!, nucleotidecount, annotate!, featureseqs, asdataframe, ispositivestrand,hasannotation, annotatedcount
export annotationcount, alignednucleotidescount, ischimeric, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread
export name, type, overlap, parts, refname, featureparams, featureparam, setfeatureparam, hastype, hasname, typein, namein, distanceonread
export hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, covratio, normalize!, hasannotationkey, readid, summarize
export mismatchcontexthist, eachpair, isfirstread
export groupfiles, difference_table
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
const FASTQ_TYPES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")

end
