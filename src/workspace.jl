using RNASeqTools
using BioSequences
using CodecZlib
using GenomicFeatures
using GFF3
using XAM

VC_GENOME_FASTA = "/home/abc/Data/vibrio/genome/NC_002505_6.fa"
VC_GENOME_GFF = "/home/abc/Data/vibrio/annotation/NC_002505_6.gff3"

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
#combine_annotations()

function check_rilseq()
    features = Features("/home/abc/Data/vibrio/annotation/NC_002505_6_utrs.gff3", ["mRNA", "5UTR", "3UTR", "tRNA", "rRNA"], "Name")
    push!(features, Interval("rybb", 1, 70, Strand('+'), RNASeqTools.Annotation("rybb1", "rybb2"))) 
    push!(features, Interval("rybb", 1, 70, Strand('-'), RNASeqTools.Annotation("rybb1", "rybb2")))
    @time alignments = PairedAlignments("/home/abc/Data/vibrio/library_rilseq/trimmed_VC3_1.bam")
    #reads = PairedReads("/home/abc/Data/vibrio/library_rilseq/trimmed_VC3_1.fasta.gz", "/home/abc/Data/vibrio/library_rilseq/trimmed_VC3_2.fasta.gz")
    annotate!(alignments, features)  
    c = 0  
    c2 = 0
    results = Dict{String, Set{LongDNASeq}}()
    counts = Dict{LongDNASeq, Int}()
    @time for (key,(alignment1, alignment2)) in alignments.dict
        !ischimeric(alignment1, alignment2) && continue
        c2 += 1
        #hasannotation(alignment2, "23Sc") && println("hey")
        if istriplet(alignment1, alignment2)
            #read1, read2 = reads.dict[key]
            #println(read1)
            #show(alignment1)
            #show(alignment2)
            #println(read2)
            #println("")
            c+=1
        end
    end
    println(c, " von ", c2)
end

#check_rilseq()