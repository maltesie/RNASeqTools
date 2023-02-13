module RNASeqTools

import XAM: BAM
using FASTX, CodecZlib, GFF3, BigWig
using BioAlignments, BioSequences, GenomicFeatures, BioGenerics
using IterTools, StringViews, Distributions, GLM, MultipleTesting

export align_mem, align_minimap, align_kraken2
export preprocess_data, trim_fastp, split_libs, download_prefetch, transform, split_interleaved, split_each_read
export Genome, Sequences, PairedSequences, AlignedReads, AlignedInterval, AlignedRead, SingleTypeFiles, PairedSingleTypeFiles, Features, Coverage
export Annotation, AlignmentAnnotation, Logo, FeatureCounts
export FastaFiles, FastagzFiles, FastqFiles, FastqgzFiles, BamFiles, GenomeFiles, GffFiles, CoverageFiles, CsvFiles
export PairedFastaFiles, PairedFastagzFiles, PaireFastqFiles, PairedFastqgzFiles
export nucleotidedistribution, annotate!, featureseqs, ispositivestrand, hasannotation, nannotated, editdistance
export ischimeric, ismulti, refinterval, readrange, refrange, annotation, hasannotation, ispositivestrand, sameread, nread
export name, type, overlap, parts, refname, params, param, setparam!, hasparam, hastype, hasname, typein, namein, distanceonread
export hasoverlap, firstoverlap, compute_coverage, merge!, merge, correlation, covratio, normalize!, readid, summarize
export add5utrs!, add3utrs!, addutrs!, addigrs!, maxdiffpositions, consensusseq, consensusbits, ninterval, eachbamrecord, eachfastqrecord, eachfastarecord
export eachpair, isfirstread, sync!
export groupfiles, filesexist, difference_glm

include("types.jl")
include("files.jl")
include("preprocess.jl")
include("sequence.jl")
include("coverage.jl")
include("annotation.jl")
include("alignment.jl")
include("counts.jl")

const FASTQ_TYPES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
const FASTA_TYPES = (".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz")

end
