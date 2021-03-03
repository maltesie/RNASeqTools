using RNASeqTools
using BioSequences
using CodecZlib
using GenomicFeatures
using GFF3
using XAM

VC_GENOME_FASTA = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"
VC_GENOME_GFF = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"



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
#test_mem_aligner()
#reads = Reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz")

#genome = Genome(VC_GENOME_FASTA)

#reads = read_reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz"; nb_reads=100)

#align_mem(reads, bam_file, VC_GENOME_FASTA)
function run_reads_filter()
    mypairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_1.fastq.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_2.fastq.gz")
    rev_comp!(mypairedreads; treat=:read1)
    query = dna"TTTCTTTGATGTCCC"
    filter!(s->occursin(query, s), mypairedreads; logic=:or)
    write("/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_1_rybb.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_2_rybb.fasta.gz", mypairedreads)
end

#run_reads_filter()

function run_reads_split()
    mypairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_1_rybb.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC3_2_rybb.fasta.gz")
    query = dna"TTTCTTTGATGTCCCCATTTTGTG"
    #println(mypairedreads.count)
    partnerreads = Reads(s->!occursin(query, s), mypairedreads)
    #rev_comp!(partnerreads)
    #cut!(partnerreads, 60; from=:left, keep=:left)
    rybbreads = Reads(s->occursin(query, s), mypairedreads)
    cut!(rybbreads, query, keep=:left_of_query)
    cut!(rybbreads, 9, from=:right, keep=:right)
    println(length(rybbreads.dict))
    bam_file = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam"
    bam_index = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam.bai"
    align_mem(partnerreads, bam_file, VC_GENOME_FASTA)
    features = open(collect, GFF3.Reader, VC_GENOME_GFF)
    filter!(x -> GFF3.featuretype(x) == "Gene", features)
    reader = open(BAM.Reader, bam_file, index=bam_index)
    counts = Dict{String, Int}()
    nucleotides = Dict{String, Dict{LongDNASeq, Int}}()
    found_ids = []
    for feature in features
        for record in eachoverlap(reader, feature)
            id = parse(UInt, BAM.tempname(record); base=2)
            println
            push!(found_ids, id)
            seq = haskey(rybbreads.dict, id) ? rybbreads.dict[id] : LongDNASeq("N")
            name = GFF3.attributes(feature, "Name")[1]
            if haskey(nucleotides, name)
                haskey(nucleotides[name], seq) ? nucleotides[name][seq] += 1 : push!(nucleotides[name], seq=>1)
            else
                push!(nucleotides, name=>Dict(seq=>1))
            end
            haskey(counts, name) ? counts[name]+=1 : push!(counts, name=>1)
        end
    end
    all_count = 0
    out_string = ""
    for (feature, count) in counts
        nc_string = join(["$seq:$c" for (seq, c) in nucleotides[feature]], ";")
        out_string *= "$feature,$count,$nc_string\n"
        all_count += count
    end
    for (key, partnerread) in partnerreads.dict
        !(key in found_ids) && println(partnerread)
    end
    out_string = "all,$all_count,\n" * out_string
    RNASeqTools.write_file("/home/abc/Workspace/RILSeq/library_rilseq/library_interactions.csv", out_string)
    #align_mem(mysinglereads, bam_file, VC_GENOME_FASTA)
    alignments = Alignments(bam_file)
    println(alignments.count)
    #for alignment in alignments.dict
    #    println(alignment)
    #end

end

#run_reads_split()

function run_find_rybb()
    genome = Genome(VC_GENOME_FASTA)
    @time println(approxsearch(genome.seq, reverse_complement(dna"TTTCTTTGATGTCCCCA"), 1))
end

files = PairedSingleTypeFiles("/home/abc/Data/vibrio/rilseq/library_rilseq/reads", ".fastq.gz")
trim_fastp(files)
