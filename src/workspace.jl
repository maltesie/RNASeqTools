using RNASeqTools
using BioSequences
using CodecZlib
using GenomicFeatures
using GFF3
using XAM

VC_GENOME_FASTA = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"
VC_GENOME_GFF = "/home/abc/Data/vibrio/annotation/NC_002505_6_rnaseq.gff3"

function my_read_bam(bam_file::String, features::Features; only_unique=true, stop_at=nothing, invertstrand=:none)
    @assert invertstrand in [:read1, :read2, :both, :none]
    reads1 = Dict{UInt, AlignedRead}()
    reads2 = Dict{UInt, AlignedRead}()
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    is_bitstring = is_bitstring_bam(bam_file)
    c = 0
    for feature in features
        for record in eachoverlap(reader, feature)
        end
    end
    close(reader)
    return reads1, reads2
end

function test_bam()
    features = Features(VC_GENOME_GFF)
    @time my_read_bam("/home/abc/Data/caulo/rilseq/trimmed_CC3_1.bam", features)
end
