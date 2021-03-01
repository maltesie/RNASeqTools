using RNASeqTools
using BioSequences
using CodecZlib

VC_GENOME_FASTA = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"
VC_GENOME_GFF = "/home/abc/Data/vibrio/annotion/NC_002505_6.gff3"

function test_mem_aligner()
    analysis_folder = "/home/abc/Workspace/RILSeq/library_rilseq/"
    genome = Genome(VC_GENOME_FASTA)
    sequences = Dict{UInt, LongDNASeq}()
    for (i, (chr, refseq)) in enumerate(genome)
        merge!(sequences, Dict(UInt(i*100+j)=>refseq[index1:index1+40]*refseq[index2:index2+60] 
            for (j,(index1,index2)) in enumerate([(rand(1:length(refseq)-40), rand(1:length(refseq)-60)) for k in 0:99])))
    end
    reads1 = Reads(sequences, "test", length(sequences))
    bam_file1 = "/home/abc/Workspace/RILSeq/library_rilseq/test.bam"
    align_mem(reads1, bam_file1, VC_GENOME_FASTA)
    alignments1 = Alignments(bam_file1)

    sequences2 = Dict{UInt, LongDNASeq}()
    for (i, (chr, refseq)) in enumerate(genome)
        merge!(sequences, Dict(UInt(i*100+j)=>refseq[index1:index1+100]
            for (j,index1) in enumerate([rand(1:length(refseq)-100) for k in 0:99])))
    end
    reads2 = Reads(sequences, "test", length(sequences))
    bam_file2 = "/home/abc/Workspace/RILSeq/library_rilseq/test2.bam"
    align_mem(reads2, bam_file2, VC_GENOME_FASTA)
end

#reads = Reads("/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC1_1.fastq.gz"; stop_at=100)
test_mem_aligner()
#reads = Reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz")

#genome = Genome(VC_GENOME_FASTA)

#reads = read_reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz"; nb_reads=100)

#align_mem(reads, bam_file, VC_GENOME_FASTA)

