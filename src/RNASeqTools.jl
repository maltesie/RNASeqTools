module RNASeqTools

using XAM
using FASTX
using BioAlignments
using CSV
using XLSX
using DataFrames
using BioSequences
using CodecZlib
using GenomicFeatures

export align_mem, align_backtrack, align_local
export trim_fastp, split_libs 
export read_bam, read_wig, read_gff, Genome, FastqReads, FastaReads
export tss, terms, annotate_utrs!
export visualize_conserved_utrs
export Reads, PairedReads

include("types.jl")
include("align.jl")
include("preprocess.jl")
include("io.jl")
include("utils.jl")
include("annotate.jl")

end
