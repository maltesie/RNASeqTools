using RNASeqTools
using BioSequences
using CodecZlib
using GenomicFeatures
using GFF3
using XAM

VC_GENOME_FASTA = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"
VC_GENOME_GFF = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"


#reads = Reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz")

#genome = Genome(VC_GENOME_FASTA)

#reads = read_reads("/home/abc/Workspace/RILSeq/library_rilseq/synthetic_reads.fasta.gz"; nb_reads=100)

#align_mem(reads, bam_file, VC_GENOME_FASTA)
function run_reads_filter()
    mypairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/trimmed_VC3_1.fastq.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/trimmed_VC3_2.fastq.gz")
    rev_comp!(mypairedreads; treat=:read1)
    query = dna"TTTCTTTGATGTCCC"
    filter!(s->occursin(query, s), mypairedreads; logic=:or)
    write("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_1.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_2.fasta.gz", mypairedreads)
end

#run_reads_filter()

function run_reads_split()
    mypairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_1.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_2.fasta.gz")
    query = dna"TTTCTTTGATGTCCC"
    #println(mypairedreads.count)
    #partnerreads = Reads(s->!occursin(query, s), mypairedreads)
    #rev_comp!(partnerreads)
    #cut!(partnerreads, 60; from=:left, keep=:left)
    rybbreads = Reads(s->occursin(query, s), mypairedreads; use_when_tied=:read1)
    #bindingsitereads = Reads(s->occursin(query, s), mypairedreads; use_when_tied=:read2)
    #cut!(bindingsitereads, query; keep=:left_of_query)
    #cut!(bindingsitereads, 9; from=:right, keep=:left)
    #cut!(bindingsitereads, 14; from=:right, keep=:right)
    #cut!(rybbreads, query; keep=:left_of_query)
    #cut!(rybbreads, 9; from=:right, keep=:right)
    #cut!(mypairedreads, query; keep=:left_of_query)
    #cut!(mypairedreads, 9; from=:right, keep=:left)
    bam_file = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam"
    bam_index = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam.bai"
    genome = Genome(VC_GENOME_FASTA)
    align_mem(mypairedreads, genome, bam_file)
    features = open(collect, GFF3.Reader, VC_GENOME_GFF)
    filter!(x -> GFF3.featuretype(x) == "Gene", features)
    reader = open(BAM.Reader, bam_file, index=bam_index)
    counts = Dict{String, Int}()
    nucleotides = Dict{String, Dict{LongDNASeq, Int}}()
    bindingnucleotides = Dict{String, Dict{LongDNASeq, Int}}()
    found_ids = []
    for feature in features
        for record in eachoverlap(reader, feature)
            id = BAM.tempname(record)
            push!(found_ids, id)
            seq = haskey(rybbreads.dict, id) ? rybbreads.dict[id] : LongDNASeq("N")
            bseq = haskey(bindingsitereads.dict, id) ? bindingsitereads.dict[id] : LongDNASeq("N")
            name = GFF3.attributes(feature, "Name")[1]
            if haskey(nucleotides, name)
                haskey(nucleotides[name], seq) ? nucleotides[name][seq] += 1 : push!(nucleotides[name], seq=>1)
            else
                push!(nucleotides, name=>Dict(seq=>1))
            end
            if haskey(bindingnucleotides, name)
                haskey(bindingnucleotides[name], bseq) ? bindingnucleotides[name][bseq] += 1 : push!(bindingnucleotides[name], bseq=>1)
            else
                push!(bindingnucleotides, name=>Dict(bseq=>1))
            end
            haskey(counts, name) ? counts[name]+=1 : push!(counts, name=>1)
        end
    end
    all_count = 0
    out_string1 = ""
    out_string2 = ""
    for (feature, count) in counts
        nc_string = join(["$seq:$c" for (seq, c) in nucleotides[feature]], ";")
        bnc_string = join(["$seq:$c" for (seq, c) in bindingnucleotides[feature]], ";")
        out_string1 *= "$feature,$count,$nc_string\n"
        out_string2 *= "$feature,$count,$bnc_string\n"
        all_count += count
    end
    #for (key, partnerread) in partnerreads.dict
    #    !(key in found_ids) && println(partnerread)
    #end
    out_string1 = "all,$all_count,\n" * out_string1
    RNASeqTools.write_file("/home/abc/Workspace/RILSeq/library_rilseq/library_interactions_rybb.csv", out_string1)
    out_string2 = "all,$all_count,\n" * out_string2
    RNASeqTools.write_file("/home/abc/Workspace/RILSeq/library_rilseq/library_interactions_other.csv", out_string2)
    #align_mem(mysinglereads, bam_file, VC_GENOME_FASTA)
    #alignments = Alignments(bam_file)
    #println(alignments.count)
    #for alignment in alignments.dict
    #    println(alignment)
    #end

end

#run_reads_split()

function run_find_rybb()
    genome = Genome(VC_GENOME_FASTA)
    @time println(approxsearch(genome.seq, reverse_complement(dna"TTTCTTTGATGTCCCCA"), 1))
end

#files = PairedSingleTypeFiles("/home/abc/Data/vibrio/rilseq/library_rilseq/reads", ".fastq.gz")
#trim_fastp(files; umi=9)
#genome = Genome(VC_GENOME_FASTA)
#pairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/trimmed_VC3_1.fastq.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/trimmed_VC3_2.fastq.gz")
#bam_file = "/home/abc/Workspace/RILSeq/library_rilseq/test.bam"
#align_mem(pairedreads, genome, bam_file)
#@time alns = PairedAlignments(bam_file)
#for (i, (key, (a1, a2))) in enumerate(alns.dict)
#    #println(BAM.alignlength(record) == BAM.rightposition(record) - BAM.position(record) + 1)
#    println(RNASeqTools.isread1(a1) || RNASeqTools.isread2(a2))
#    println(LongDNASeq(BAM.sequence(a1)) == pairedreads.dict[key][1])
#    i == 5 && break
#end
#println(BAM.header(reader))
#pairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_1.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_2.fasta.gz")
#reads = Reads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_1.fasta.gz")
#for (i, (read1, read2)) in enumerate(pairedreads)
#    println(read1)
#    i == 5 && break
#end

function te()
    a = collect(1:20)
    b = Int[]
    c = Int[]
    for i in a
        current_vec = isodd(i) ? b : c
        push!(current_vec, i)
    end
    println(b, "\n", c)
end
#te()

function test_mem_aligner()
    analysis_folder = "/home/abc/Workspace/RILSeq/library_rilseq/"
    genome = Genome(VC_GENOME_FASTA)
    #write("/home/abc/Data/test/genome.fa", genome)
    sequences = Dict{String, LongDNASeq}()
    for (i, (chr, refseq)) in enumerate(genome)
        merge!(sequences, Dict("$(i*100+j)"=>refseq[index1:index1+30]*refseq[index2:index2+40]*refseq[index3:index3+50] 
            for (j,(index1,index2,index3)) in enumerate([(rand(1:length(refseq)-30), rand(1:length(refseq)-40), rand(1:length(refseq)-50)) for k in 0:99])))
    end
    reads1 = Reads(sequences, "test")
    #println(length(reads1), "\n", reads1)
    reads1_file = "/home/abc/Data/test/chimeric_reads.fasta.gz"
    write(reads1_file, reads1)
    reads1 = Reads(reads1_file)
    #println(length(reads1), "\n", reads1)
    bam_file1 = "/home/abc/Workspace/RILSeq/library_rilseq/test.bam"
    align_mem(reads1, genome, bam_file1)

    sequences2 = Dict{String, LongDNASeq}()
    for (i, (chr, refseq)) in enumerate(genome)
        merge!(sequences, Dict("$(i*100+j)"=>refseq[index1:index1+100]
            for (j,index1) in enumerate([rand(1:length(refseq)-100) for k in 0:99])))
    end
    reads2 = Reads(sequences, "test")
    bam_file2 = "/home/abc/Workspace/RILSeq/library_rilseq/test2.bam"
    #align_mem(reads2, genome, bam_file2)
    #bam_file = "/home/abc/Data/vibrio/rilseq/micha_rilseq/pe_bams/hfq_1_0.2.bam"
    #reader = BAM.Reader(open(bam_file1))
    #c1 = 0
    #c2 = 0
    #for (i, record) in enumerate(reader)
    #    BAM.isprimary(record) && (println(BAM.cigar(record)))
    #    if haskey(record, "XA") 
    #        println("XA ",record["XA"]::String)
    #        println("SA ",record["SA"]::String)
    #    end
    #end#
    #println(c1, " ", c2)
    @time alignments = Alignments(bam_file1)
    for entry in alignments
        println(entry)
    end
    #println(length(collect(keys(alignments))))
end

#reads = Reads("/home/abc/Data/vibrio/rilseq/library_rilseq/trimmed/VC1_1.fastq.gz"; stop_at=100)
#test_mem_aligner()

#RNASeqTools.seqpositions("72S51M")

function test_abc()
    bam = "/home/abc/Data/vibrio/phage_inside/trimmed_MS-10_S37_R1_001.bam"
    phage_gff = "/home/abc/Data/vp882/annotation/GCF_000868005.1_ViralProj18851_genomic.gff"
    genome_gff = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"
    vc_features = Features(genome_gff; type="Gene")
    phage_features = Features(phage_gff; type="gene")
    #@time alns = Alignments(bam, stop_at=5000000)
    #@time annotate!(alns, phage_features)
    #@time annotate!(alns, vc_features)
    
    for (i,al) in enumerate(alns)
        i % 200000 == 0 && println(al)
    end
end

#test_abc()

function test_annotation()
    features1 = Features("/home/abc/Data/vibrio/annotation/NC_002505_6.gff3")
    features2 = Features("/home/abc/Data/vp882/annotation/GCF_000868005.1_ViralProj18851_genomic.gff")
    #@time reads = Reads("/home/abc/Data/vibrio/phage_inside/trimmed_MS-10_S37_R1_001.fastq"; stop_at=1000000)
    @time alignments = Alignments("/home/abc/Data/vibrio/phage_inside/trimmed_MS-10_S37_R1_001.bam"; stop_at=5000000)
    @time annotate!(alignments, features1)
    @time annotate!(alignments, features2)
    #print(length(alignments))
    #features = open(collect, GFF3.Reader, "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3")
    for (i,alignment) in enumerate(alignments)
        println(alignment)
        i == 2 && break
    end

    # Keep mRNA features.
    #filter!(x -> GFF3.featuretype(x) == "Gene", features)

    # Open a BAM file and iterate over records overlapping mRNA transcripts.
    #reader = open(BAM.Reader, "/home/abc/Data/vibrio/phage_inside/trimmed_MS-10_S37_R1_001.bam", index = "/home/abc/Data/vibrio/phage_inside/trimmed_MS-10_S37_R1_001.bam.bai")
    #@time for feature in features
    #    for record in eachoverlap(reader, feature)
    #        # `record` overlaps `feature`.
    #        # ...
    #    end
    #end
    #close(reader)
    #ints = IntervalCollection([Interval("test", 1, 10, '-'), Interval("test", 5, 15, '-'), Interval("test", 15, 20, '-')])
    #cov = coverage(ints)
    #println(cov)
    #writer = BigWig.Writer("/home/abc/Workspace/test.bigwig")
end

#test_annotation()

function run_reads_split2()
    mypairedreads = PairedReads("/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_1.fasta.gz", "/home/abc/Data/vibrio/rilseq/library_rilseq/reads/rybb_VC3_2.fasta.gz")
    query = dna"TTTCTTTGATGTCCC"
    #rybbreads = Reads(s->occursin(query, s), mypairedreads; use_when_tied=:read1)
    bam_file = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam"
    bam_index = "/home/abc/Workspace/RILSeq/library_rilseq/chimeras.bam.bai"
    genome = Genome(VC_GENOME_FASTA)
    align_mem(mypairedreads, genome, bam_file)
    @time alignments = PairedAlignments(bam_file)
    features = open(collect, GFF3.Reader, VC_GENOME_GFF)
    println(typeof(features))
    filter!(x -> GFF3.featuretype(x) == "Gene", features)
    for feature in features
        println(GFF3.attributes(feature, "Name"))
        break
    end
    count1 = 0
    count2 = 0
    betweens = Dict{LongDNASeq, Int}()
    for (key, (alignment1, alignment2)) in alignments.dict
        if !isempty(alignment1)
            read = mypairedreads.dict[key][1] 
            if occursin(query, read)
                count1+=1
                slice = findfirst(query, read)
                #println(RNASeqTools.primaryalignmentpart(alignment1).readstart, " ", RNASeqTools.primaryalignmentpart(alignment1).readstop, " ", slice[1], "\n", read)
                between = read[RNASeqTools.primaryalignmentpart(alignment1).readstop:slice[1]]
                between in keys(betweens) ? betweens[between] += 1 : push!(betweens, between=>1)
            end
        end
        if !isempty(alignment2) 
            read = mypairedreads.dict[key][2] 
            if occursin(query, read)
                count2+=1
                slice = findfirst(query, read)
                between = read[RNASeqTools.primaryalignmentpart(alignment2).readstop:slice[1]]
                between in keys(betweens) ? betweens[between] += 1 : push!(betweens, between=>1)
            end
        end
    end
    println(betweens)
end

#run_reads_split2()

function invert_reads1()
    trimmed_files = PairedSingleTypeFiles("/home/abc/Data/vibrio/library_rilseq/", ".fastq.gz"; prefix="trimmed")
    for (file1, file2) in trimmed_files
        reads = Reads(file2)
        #rev_comp!(reads)
        out = file2[1:end-4] * "a.gz"
        write(out, reads)
    end
end

#invert_reads1()

function align_library_rilseq()
    trimmed_files = PairedSingleTypeFiles("/home/abc/Data/vibrio/library_rilseq/", ".fasta.gz"; prefix="trimmed")
    genome = Genome("/home/abc/Data/vibrio/genome/NC_002505_6.fa")
    rybb_chr = Genome(LongDNASeq("TTTCTTTGATGTCCCCATTTTGTGGAGCCCATCAACCCCGCCATTTCGGTTCAAGGTTGATGGGTTTTTT"), "rybb")
    combined_genome = genome * rybb_chr
    align_mem(trimmed_files, combined_genome)
end

#align_library_rilseq()

function library_rilseq()
    features = Features("/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"; type="Gene")
    push!(features, Interval("rybb", 1, 70, Strand('+'), RNASeqTools.Annotation("rybb1", "rybb2"))) 
    push!(features, Interval("rybb", 1, 70, Strand('-'), RNASeqTools.Annotation("rybb1", "rybb2")))
    @time alignments = PairedAlignments("/home/abc/Data/vibrio/library_rilseq/trimmed_VC3_1.bam"; stop_at=100000)
    @time annotate!(alignments, features)    
    for (i,(alignment1, alignment2)) in enumerate(alignments)
        if RNASeqTools.hasannotation(alignment1, "rybb2") || RNASeqTools.hasannotation(alignment2, "rybb2")
            show(alignment1)
            show(alignment2)
            println("")
        end
    end
end

#align_library_rilseq()
#library_rilseq()

function test_reads()
    reads = Reads("/home/abc/Data/vibrio/library_rilseq/rybb_VC3_1.fasta.gz")
    @time for read in reads
    end
end

#test_reads()

using GFF3, CSV, DataFrames, GenomicFeatures, BioGenerics

function csv_to_feature(file::String)
    table = CSV.read(file, DataFrame)
    intervals = Interval{RNASeqTools.Annotation}[]
    for row in eachrow(table)
        strand = row.start > 0 ? Strand('+') : Strand('-')
        (start, stop) = row.start > 0 ? (row.start, row.stop) : (-row.stop, -row.start)
        push!(intervals, Interval(row.chr, row.start, row.stop, strand, RNASeqTools.Annotation("srna", row.name)))
    end
    return Features(IntervalCollection(intervals, true), nothing)
end

function combine_annotations()
    records = GFF3.Record[]
    genome_gff = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"
    all = Features("/home/abc/Data/vibrio/annotation/NC_002505_6.gff3", "Gene", "Name")
    all_set = Set(f.metadata.name for f in all)
    srnas = Set(f.metadata.name for f in csv_to_feature("/home/abc/Data/vibrio/annotation/NC_002505_6_srna.csv"))
    trnas = Set(f.metadata.name[2:end-1] for f in Features(genome_gff, "tRNA", "gene"))
    rrnas = Set(f.metadata.name[2:end-1] for f in Features(genome_gff, "rRNA", "gene"))
    mrnas = setdiff(all_set, union(srnas, trnas, rrnas))
    #println(length(srnas))
    #println(rrnas)
    #println(trnas)
    #println(mrnas)
    utr5 = CSV.read("/home/abc/Workspace/ConservedUTRs/vibrio_first_round/utrNC_002505.csv", DataFrame)
    utr6 = CSV.read("/home/abc/Workspace/ConservedUTRs/vibrio_first_round/utrNC_002506.csv", DataFrame)
    for utr in [utr5, utr6]
        chr = utr === utr5 ? "NC_002505" : "NC_002506"
        for row in eachrow(utr)
            row.name in mrnas || continue
            is_negative_strand = row.start < 0
            five_start, five_stop, three_start, three_stop = sort(abs.([row.start, row.stop, row.threeUTR, row.fiveUTR]))
            is_negative_strand && ((five_start, five_stop, three_start, three_stop) = (three_start, three_stop, five_start, five_stop))
            five_stop -= 1
            three_start += 1
            strand = is_negative_strand ? '-' : '+'
            five_type = row.fiveType == "guess" ? "guess" : "rnaseq"
            three_type = row.threeType == "guess" ? "guess" : "rnaseq"
            line5 = "$chr\t.\t5UTR\t$five_start\t$five_stop\t.\t$strand\t.\tName=$(row.name);From=$five_type"
            line3 = "$chr\t.\t3UTR\t$three_start\t$three_stop\t.\t$strand\t.\tName=$(row.name);From=$three_type"
            push!(records, GFF3.Record(line5))
            push!(records, GFF3.Record(line3))
        end
    end
    for feature in all
        if RNASeqTools.name(feature) in srnas
            type = "sRNA"
        elseif RNASeqTools.name(feature) in rrnas
            type = "rRNA"
        elseif RNASeqTools.name(feature) in trnas
            type = "tRNA"
        elseif RNASeqTools.name(feature) in mrnas
            type = "mRNA"
        end
        start = leftposition(feature)
        stop = rightposition(feature)
        strand = GenomicFeatures.strand(feature)
        chr = seqname(feature)
        name = RNASeqTools.name(feature)
        line = "$chr\t.\t$type\t$start\t$stop\t.\t$strand\t.\tName=$name"
        push!(records, GFF3.Record(line))
    end
    records5 = [record for record in records if GFF3.seqid(record) == "NC_002505"]
    records6 = [record for record in records if GFF3.seqid(record) == "NC_002506"]
    sort!(records5, by=x->GFF3.seqstart(x))
    sort!(records6, by=x->GFF3.seqstart(x))
    writer = GFF3.Writer(open("/home/abc/Data/vibrio/annotation/NC_002505_6_utrs.gff3", "w"))
    for record in records5
        write(writer, record)
    end
    for record in records6
        write(writer, record)
    end
    close(writer)
end
combine_annotations()