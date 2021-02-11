using DataFrames
using CSV
using XLSX
using Statistics
using XAM
using BioAlignments
using CodecZlib

function convert_annotation(xls_annot::String, csv_annot::String)
    df = DataFrame(XLSX.readtable(xls_annot, 1)...)
    chrs = [val==1 ? "NC_002505" : "NC_002506" for val in df[!, :Chromosome]]
    starts = [startswith(cont, "compl") ? -1*parse(Int, split(split(cont,"..")[2],")")[1]) : parse(Int, split(cont,"..")[1]) for cont in df[!, :Region]]
    stops =  [startswith(cont, "compl") ? -1*parse(Int, split(split(cont,"..")[1],"(")[2]) : parse(Int, split(cont,"..")[2]) for cont in df[!, :Region]]
    out_df = DataFrame(chr=chrs, name=df[!, :Name], start=starts, stop=stops)
    sort!(out_df, :start)
    CSV.write(csv_annot, out_df);
end

#convert_annotation("/home/malte/Workspace/data/V.ch._sRNAs.xlsx", "/home/malte/Workspace/data/vibrio_srna.csv")

function convert_annotation_to_positive!(annotation::DataFrame)
    annotation[!, :strand] = fill("+", nrow(annotation))
    for row in eachrow(annotation)
        if (row[:start] < 0) 
            row[:start], row[:stop] = abs(row[:stop]), abs(row[:start])
            row[:strand] = "-"
        end
    end
    sort!(annotation, [:chr,:start])
end

function convert_hmm_to_table(hmm_file_chr1::String, hmm_file_chr2::String, srna_annotations::String, out_table::String)
    out_df = CSV.read(srna_annotations, DataFrame)
    convert_annotation_to_positive!(out_df)
    df = Dict("NC_002505"=>CSV.read(hmm_file_chr1, DataFrame), "NC_002506"=>CSV.read(hmm_file_chr2, DataFrame))
    temp_pos = Int[]
    temp_class = Int[]
    temp_reads = Int[]
    XLSX.openxlsx(abspath(out_table), mode="w") do xf
        sheet = xf[1]
        out_df[!,:essential] = fill("-", nrow(out_df))
        out_df[!,:growth] = fill("-", nrow(out_df))
        out_df[!, :average_reads] = fill("-", nrow(out_df))
        for row in eachrow(out_df)
            index = row[:start] .<= df[row[:chr]][!,Symbol("TA start")] .<= row[:stop]
            temp_pos = df[row[:chr]][index, Symbol("TA start")]
            temp_reads = df[row[:chr]][index, :Reads]
            row[:chr] == "NC_002505" ? temp_class = df[row[:chr]][index, Symbol("Iteration 6")] : temp_class = df[row[:chr]][index, Symbol("Iternation 2")] 
            (1 in temp_class) && (row[:essential] = "$(sum(temp_class .== 1))/$(length(temp_class))")
            (3 in temp_class) && (row[:growth] = "$(sum(temp_class .== 3))/$(length(temp_class))")
            (sum(temp_reads) > 0) && (sum(temp_class .== 1) == length(temp_class)) && (row[:essential] *= "*")
            row[:average_reads] = "$(round(mean(temp_reads), digits=1)) pm $(round(std(temp_reads), digits=1))"
        end
        XLSX.writetable!(sheet, [c for c in eachcol(out_df)], DataFrames.names(out_df)) 
    end 
end

#hmm1 = "/home/malte/Workspace/data/HMMposterior_probabilities_chr1.txt"
#hmm2 = "/home/malte/Workspace/data/HMMposterior_probabilities_chr2.txt"
#ann = "/home/malte/Workspace/data/vibrio_srna.csv"
#convert_hmm_to_table(hmm1, hmm2, ann, "../data/hmm_output.xlsx")

function convert_edges(edges_file::String, nodes_file::String, out_file::String)
    out_df = CSV.read(edges_file, DataFrame)
    nodes_df = CSV.read(nodes_file, DataFrame)
    nodes = Dict(row[:altName]=>row[:color] for row in eachrow(nodes_df))
    
    out_df[!, :fromColor] = fill("black", nrow(out_df))
    out_df[!, :toColor] = fill("black", nrow(out_df))

    for row in eachrow(out_df)
        row[:fromColor] = nodes[row[:fromAltName]] 
        row[:toColor] = nodes[row[:toAltName]] 
    end

    CSV.write(out_file, out_df)
end

#edges = "/home/malte/Workspace/Dash/CoexpressionGraph/assets/edges.txt"
#nodes = "/home/malte/Workspace/Dash/CoexpressionGraph/assets/nodes.txt"
#outedges = "/home/malte/Workspace/Dash/CoexpressionGraph/assets/edges_new.txt"
#convert_edges(edges, nodes, outedges)

function motif_search(motif_file::String, out_fasta::String, out_csv::String)
    df = CSV.read(motif_file, DataFrame, header=38, datarow=39, delim="\t")
    df = df[df[!,:evidence_level] .== "Strong", :]
    sequences = ">id0\n" * join([uppercase(s[1:end-21])*"\n>id$i" for (i,s) in enumerate(df[!, :sequence])], "\n")
    CSV.write(out_csv, df)
    open(out_fasta,"w") do io
        println(io, sequences[1:end-7])
    end
end

function motif_search(motif_file::String, out_fasta::Array{String, 1})
    df = CSV.read(motif_file, DataFrame)
    sequences = ""
    sequences_1 = ""
    sequences_2 = ""
    for (i, row) in enumerate(eachrow(df))
        println(row)
        sequences *= ">id$i\n" * uppercase(row[:sequence][1:end-21]) * "\n"
        (row[:class] == 1) && (sequences_1 *= ">id$i\n" * uppercase(row[:sequence][1:end-21]) * "\n")
        (row[:class] == 2) && (sequences_2 *= ">id$i\n" * uppercase(row[:sequence][1:end-21]) * "\n")
    end
    print(sequences)
    print("\n\n")
    print(sequences_2)
    print("\n\n")
    print(sequences_1)
    open(out_fasta[1],"w") do io
        println(io, sequences_1)
    end
    open(out_fasta[2],"w") do io
        println(io, sequences_2)
    end
    open(out_fasta[3],"w") do io
        println(io, sequences)
    end
end

#motif_search("/home/malte/Downloads/Promoter_class.csv", ["/home/malte/Downloads/Promoter1.fasta", "/home/malte/Downloads/Promoter2.fasta", "/home/malte/Downloads/Promoter3.fasta"])

function get_reference_info(bam_file::String)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")

    for meta in findall(BAM.header(reader), "SQ")
        println(meta["SN"])
    end

    while !eof(reader)
        read!(reader, record)
        break
    end
    close(reader)
end

#bam_file = "/home/malte/Workspace/data/test.bam"

#get_reference_info(bam_file)

#function align(sequence_file::String, reference_file::String)

#end

function read_coverage(wig_file::String)
    coverage = Dict()
    temp_collect = Dict()
    span = 0
    chr = ""
    open(wig_file, "r") do file
        for line in split(read(file, String), "\n")[1:end-1]
            startswith(line, "track") && continue
            if startswith(line, "variableStep")
                span = parse(Int, split(split(line, " ")[3],"=")[2])
                chr = split(split(line, " ")[2],"=")[2]
                temp_collect[chr] = Tuple{Int, Float64}[]
            else
                str_index, str_value = split(line, " ")
                index = parse(Int, str_index)
                value = parse(Float64, str_value)

                for i in index:index+span
                    push!(temp_collect[chr], (i, value))
                end
            end
        end
    end
    for (chr, points) in temp_collect
        coverage[chr] = zeros(Float64, points[end][1])
        for (index, value) in points
            coverage[chr][index] = value
        end
    end
    return coverage
end

#wig = "/home/malte/Workspace/dRNASeq/data/reademption_campbelli/output/coverage/coverage-tnoar_min_normalized/Vcamp-luxR-Rep1_S34_R1_001_trimmed_div_by_12872450.0_multi_by_11533765.0_forward.wig"

#read_coverage(wig)

#function test_utr_annotation(coverage_files::Array{String, 2})

#end

function run_preprocess_drnaseq()
    read_files = ["/home/malte/Workspace/data/vibrio/drnaseq/tex_01_1.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_01_2.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_1.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_2.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_1.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_2.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_1.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_2.fastq.gz"]

    adapter = fill("AATGATACGGCGACCACCG", 8)

    output_folder = "/home/malte/Workspace/data/vibrio/drnaseq/"

    fastp_bin_path = "/home/malte/Tools/fastp"

    trim_fastp(read_files, output_folder; adapters=adapter, fastp_bin=fastp_bin_path)
end

#run_preprocess_drnaseq()

function run_align_drnaseq()
    read_files = ["/home/malte/Workspace/data/vibrio/drnaseq/tex_01_1_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_01_2_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_1_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_2_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_1_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_2_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_1_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_2_trimmed.fastq.gz"]

    out_files = ["/home/malte/Workspace/data/vibrio/drnaseq/tex_01_1.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_01_2.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_1.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_2.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_1.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_2.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_1.bam",
    "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_2.bam"]

    reference = "/home/malte/Workspace/data/vibrio/annotation/NC_002505_6.fa"

    bwa = "/home/malte/Tools/bwa/bwa"

    sam = "/home/malte/Tools/samtools/bin/samtools"

    for (in_file, out_file) in zip(read_files, out_files)
        align_backtrack(in_file, out_file, reference; max_miss=3, bwa_bin=bwa, sam_bin=sam)
    end
end

#run_align_drnaseq()

function run_coverage_drnaseq()

    folder = "/home/malte/Workspace/data/vibrio/drnaseq/"

    bam_files = [(joinpath(folder, "tex_01_1.bam"), joinpath(folder, "tex_01_2.bam")),
    (joinpath(folder, "tex_20_1.bam"), joinpath(folder, "tex_20_2.bam")),
    (joinpath(folder, "notex_01_1.bam"), joinpath(folder, "notex_01_2.bam")),
    (joinpath(folder, "notex_20_1.bam"), joinpath(folder, "notex_20_2.bam"))]

    wig_files = [(joinpath(folder, "tex_01_f.wig"), joinpath(folder, "tex_01_r.wig")),
    (joinpath(folder, "tex_20_f.wig"), joinpath(folder, "tex_20_r.wig")),
    (joinpath(folder, "notex_01_f.wig"), joinpath(folder, "notex_01_r.wig")),
    (joinpath(folder, "notex_20_f.wig"), joinpath(folder, "notex_20_r.wig"))]

    for ((bam1, bam2),(wigf, wigr)) in zip(bam_files, wig_files)
        (coverage_f1, coverage_r1) = coverage(bam1)
        (coverage_f2, coverage_r2) = coverage(bam2)
        write_wig([coverage_f1, coverage_f2], wigf)
        write_wig([coverage_r1, coverage_r2], wigr)
    end
end

#run_coverage_drnaseq()

function run_preprocess_termseq()
    read_files = ["/home/malte/Workspace/data/vibrio/termseq/term_1.fastq.gz",
    "/home/malte/Workspace/data/vibrio/termseq/term_2.fastq.gz",
    "/home/malte/Workspace/data/vibrio/termseq/term_3.fastq.gz",]

    adapter = fill("AGATCGGAAGAG", 3)

    output_folder = "/home/malte/Workspace/data/vibrio/termseq/"

    fastp_bin_path = "/home/malte/Tools/fastp"

    trim_fastp(read_files, output_folder; adapters=adapter, fastp_bin=fastp_bin_path)
end

#run_preprocess_termseq()

function run_align_termseq()
    read_files = ["/home/malte/Workspace/data/vibrio/termseq/term_1_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/termseq/term_2_trimmed.fastq.gz",
    "/home/malte/Workspace/data/vibrio/termseq/term_3_trimmed.fastq.gz",]

    out_files = ["/home/malte/Workspace/data/vibrio/termseq/term_1.bam",
    "/home/malte/Workspace/data/vibrio/termseq/term_2.bam",
    "/home/malte/Workspace/data/vibrio/termseq/term_3.bam"]

    reference = "/home/malte/Workspace/data/vibrio/annotation/NC_002505_6.fa"

    bwa = "/home/malte/Tools/bwa/bwa"

    sam = "/home/malte/Tools/samtools/bin/samtools"

    for (in_file, out_file) in zip(read_files, out_files)
        align_backtrack(in_file, out_file, reference; max_miss=3, bwa_bin=bwa, sam_bin=sam)
    end
end

#run_align_termseq()

function run_coverage_termseq()

    folder = "/home/malte/Workspace/data/vibrio/termseq/"

    bam_files = [(joinpath(folder, "term_1.bam"), joinpath(folder, "term_2.bam"), joinpath(folder, "term_3.bam"))]

    wig_files = [(joinpath(folder, "term_f.wig"), joinpath(folder, "term_r.wig"))]

    for ((bam1, bam2, bam3),(wigf, wigr)) in zip(bam_files, wig_files)
        (coverage_f1, coverage_r1) = coverage(bam1, is_reversed=true)
        (coverage_f2, coverage_r2) = coverage(bam2, is_reversed=true)
        (coverage_f3, coverage_r3) = coverage(bam3, is_reversed=true)
        write_wig([coverage_f1, coverage_f2, coverage_f3], wigf)
        write_wig([coverage_r1, coverage_r2, coverage_r3], wigr)
    end
end

#run_coverage_termseq()

function run_utr_annotation()
    gff = "/home/malte/Workspace/data/vibrio/annotation/NC_002505_6.gff3"
    drna_folder = "/home/malte/Workspace/data/vibrio/drnaseq/"
    term_folder = "/home/malte/Workspace/data/vibrio/termseq/"

    drna_notex_fs = [joinpath(drna_folder, "notex_01_f.wig"), joinpath(drna_folder, "notex_20_f.wig")]
    drna_tex_fs = [joinpath(drna_folder, "tex_01_f.wig"), joinpath(drna_folder, "tex_20_f.wig")]
    drna_notex_rs = [joinpath(drna_folder, "notex_01_r.wig"), joinpath(drna_folder, "notex_20_r.wig")]
    drna_tex_rs = [joinpath(drna_folder, "tex_01_r.wig"), joinpath(drna_folder, "tex_20_r.wig")]

    term_fs = [joinpath(term_folder, "term_f.wig")]
    term_rs = [joinpath(term_folder, "term_r.wig")]

    tsss = tss(drna_notex_fs, drna_notex_rs, drna_tex_fs, drna_tex_rs)
    termss = terms(term_fs, term_rs)
    
    annotation = read_annotations(gff)

    annotate_utrs!(annotation, tsss, termss)
    for (chr, ann) in annotation
        CSV.write("/home/malte/Workspace/data/vibrio/utr$chr.csv", ann)
    end
end

#run_utr_annotation()
#drna_folder = "/home/malte/Workspace/data/vibrio/drnaseq/"
#get_chr_from_wig(joinpath(drna_folder, "tex_01_r.wig"))

#run_utr_annotation()

function translate_position_mg1655_to_bw25113(pos::Int)::Union{Nothing, Int}
    
    if abs(pos) < 66533
        return pos
    
    elseif (66533 <= abs(pos) <= 70072) 
        return nothing
    elseif (70072 < abs(pos) < 257900)
        return sign(pos) * (abs(pos) - (3540-27)) 
    
    elseif (257900 <= abs(pos) <= 258675)
        return nothing
    elseif (258675 < abs(pos) < 364356)
        return sign(pos) * (abs(pos) - ((3540-27)+776)) 
    
    elseif (364356 <= abs(pos) <= 366418)
        return nothing
    elseif (366418 < abs(pos) <= 371758)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585))) 
    elseif (371758 < abs(pos) < 1299495)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223)) 
    
    elseif (1299495 <= abs(pos) <= 1300693)
        return nothing
    elseif (1300693 < abs(pos) < 1978495)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199))

    elseif (1978495 <= abs(pos) <= 1979270)
        return nothing
    elseif (1979270 < abs(pos) < 2521007)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776))

    elseif (2173363 <= abs(pos) <= 2173364)
        return nothing
    elseif (2173364 < abs(pos) < 2521007)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2))

    elseif (2521007 <= abs(pos) <= 2521126)
        return nothing
    elseif (2521126 < abs(pos) <= 3560456)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120))
    elseif (3560456 < abs(pos) < 4093991)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1))

    elseif (4093991 <= abs(pos) <= 4097451)
        return nothing
    elseif (4097451 < abs(pos) < 4296282)
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1+(3461-29)))

    elseif (4296282 <= abs(pos) <= 4296392)
        return nothing
    elseif (4296392 < abs(pos))
        return sign(pos) * (abs(pos) - ((3540-27)+776+(2063-585)-1223+1199+776+2+120-1+(3461-29)+111))
    end
end

#println(translate_position_mg1655_to_bw25113.([4602509,4233758,1652331,1652331,1652331]))

function run_translation()
    in_file = "/home/malte/Workspace/Dash/FindFitness/assets/ecoli/srna_mg1655.csv"
    out_file = "/home/malte/Workspace/Dash/FindFitness/assets/ecoli/srna_bw25113.csv"
    annotations = CSV.read(in_file, DataFrame)

    annotations[!,:Left] = translate_position_mg1655_to_bw25113.(annotations[!,:Left])
    annotations[!,:Right] = translate_position_mg1655_to_bw25113.(annotations[!,:Right])

    CSV.write(out_file, annotations)
end

function translation_dict(from_sequence::String, to_sequence::String)
    s1 = LongDNASeq(from_sequence)
    s2 = LongDNASeq(to_sequence)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
    res = alignment(pairalign(GlobalAlignment(), s1, s2, scoremodel))
    trans_dict::Dict{Int, Union{Nothing, Int}} = Dict()
    from_pos = 1
    to_pos = 1
    for (b1, b2) in collect(res)
        if (b1 == '-')
            to_pos += 1
        elseif (b2 == '-')
            transdict[from_pos] = nothing
            from_pos += 1
        else
            trans_dict[from_pos] = to_pos
            to_pos += 1
            from_pos += 1
        end
    end
    return trans_dict
end

function local_alignment(long_sequence::String, short_sequence::String)
    s1 = LongDNASeq(long_sequence)
    s2 = LongDNASeq(short_sequence)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
    res = alignment(pairalign(LocalAlignment(), s1, s2, scoremodel))
    return collect(res)
end

function run_local_align()
    check_sequence = "ATTTCTCTGAGATGTTCGCAAGCGGGCCAGTCCCCTGAGCCGATATTTCATACCACAAGAATGTGGCGCTCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGAGCATGCTCGGTCCGTCCGAGAAGCCTAAAAACTGCGACGACACATTCACCTTGAACCAAGGGTTCAAGGGTTAAAAAACAGCCTGCGGCGGCATCTCGGAGATTC"
    @time genome = read_genomic_fasta("/home/malte/Workspace/data/ecoli/annotations/mg1655.fna")

    for (chr1, seq1) in genome
        @time test = align_local(seq1, check_sequence)
        for (a, b) in test
            println("$a $b")
        end
    end
end

#run_local_align()

function run_translation2()
    in_file = "/home/malte/Workspace/Dash/FindFitness/assets/ecoli/sRNAs_BW25113_coordinates.csv"
    out_file = "/home/malte/Workspace/Dash/FindFitness/assets/ecoli/srna_mg1655_new.csv"

    genome1 = read_genomic_fasta("/home/malte/Workspace/data/ecoli/annotations/bw25113.fna")
    genome2 = read_genomic_fasta("/home/malte/Workspace/data/ecoli/annotations/mg1655.fna")

    for (chr1, seq1) in genome1, (chr2, seq2) in genome2
        @time td = translation_dict(seq1[1:5000], seq2[50:5050])
        #annotations = CSV.read(in_file, DataFrame)

        #translate_positions!(annotations[!,:start])
        #translate_positions!(annotations[!,:stop])

        #CSV.write(out_file, annotations)
    end
end

function align_levenshtein(in_file::String, genome_file::String)
    genome = read_genomic_fasta(genome_file)
    reads = read_reads_fastq(in_file)

    for (chr, genome_sequence) in genome, (name,read) in reads
        @time test = pairalign(LevenshteinDistance(), genome_sequence, read; distance_only=true)
    end
end

function run_write_utrs_fasta(annotations::Dict{String,DataFrame}, genome::Dict{String,String})
    genome_file = "/home/malte/Workspace/data/vibrio/annotation/NC_002505_6.fa"
    genome = read_genomic_fasta(genome_file)
    annotation_5 = "/home/malte/Workspace/data/vibrio/utrNC_002505.csv"
    annotation_6 = "/home/malte/Workspace/data/vibrio/utrNC_002506.csv"
    annotations = Dict("NC_002505"=>CSV.read(annotation_5, DataFrame), "NC_002506"=>CSV.read(annotation_6, DataFrame))
    record = FASTA.Record()
    five_file = "/home/malte/Workspace/data/vibrio/fiveUTR.fasta.gz"
    three_file = "/home/malte/Workspace/data/vibrio/threeUTR.fasta.gz"
    write_utrs_fasta()
end

#save_utrs()

function run_multi_genome_align()
    genome_folder = "/home/abc/Data/vibrio/genomes/"
    genomes = [joinpath(genome_folder, file) for file in readdir(genome_folder) if endswith(file, ".fna")]

    fiveutrs = "/home/abc/Workspace/ConservedUTRs/fiveUTR.fasta.gz"
    threeutrs = "/home/abc/Workspace/ConservedUTRs/threeUTR.fasta.gz"

    fivefolder = "/home/abc/Workspace/ConservedUTRs/fiveUTRs/"
    threefolder = "/home/abc/Workspace/ConservedUTRs/threeUTRs/"

    align_mem(fiveutrs, fivefolder, genomes)
    align_mem(threeutrs, threefolder, genomes)
end

#run_multi_genome_align()

function conservation_table(bam_files::Vector{String})
    name_trans = Dict(
     "GCF_000011805.1_ASM1180v1_genomic" => "Aliivibrio fischeri ES114",
     "GCF_000024825.1_ASM2482v1_genomic" => "Vibrio antiquarius",
     "GCF_000039765.1_ASM3976v1_genomic" => "Vibrio vulnificus CMCP6",
     "GCF_000184325.1_ASM18432v1_genomic" => "Vibrio furnissii NCTC 11218",
     "GCF_000196095.1_ASM19609v1_genomic" => "Vibrio parahaemolyticus RIMD 2210633",
     "GCF_000196495.1_ASM19649v1_genomic" => "Aliivibrio salmonicida LFI1238",
     "GCF_000217675.1_ASM21767v1_genomic" => "Vibrio anguillarum 775",
     "GCF_000241385.1_ASM24138v1_genomic" => "Vibrio sp. EJY3",
     "GCF_001558435.2_ASM155843v2_genomic" => "Vibrio harveyi FDAARGOS_107",
     "GCF_011212705.1_ASM1121270v1_genomic" => "Vibrio coralliilyticus OCN008",
     "GCF_900205735.1_N16961_v2_genomic" => "Vibrio cholerae N16961",
    )
    names = [name_trans[name] for name in [join(split(basename(file),".")[1:end-1], ".") for file in bam_files]]
    table = DataFrame(merge(Dict("name"=>String[], "score"=>Float64[]), Dict(name=>Union{Missing,Int}[] for name in names)))
    for (file, name) in zip(bam_files, names)
        alignments = read_bam(file)
        for row in eachrow(alignments)
            (row[:name] in table[!,:name]) || append!(table, DataFrame(merge(Dict("name"=>row[:name], "score"=>0.0), Dict(name=>missing for name in names))))
            table[table.name .== row[:name], Symbol(name)] .= row[:nm]
            table[table.name .== row[:name], :score] .+= (1.0 - Float64(row[:nm])/50.0)
            #(table[table.name .== row[:name], :score][1] > 10) &&  println((1.0 - Float64(row[:nm])/50.0), " $name $(row[:name])")
            #((1.0 - Float64(row[:nm])/50.0) > 1) && println((1.0 - Float64(row[:nm])/50.0), " $name $(row[:name])")
            #(sum(table.name .== row[:name]) > 1) && println("$name $nb_hits")
        end
    end
    return table
end

function run_constervation_table()
    five_folder = "/home/abc/Workspace/ConservedUTRs/fiveUTRs/"
    five_files = [joinpath(five_folder, file) for file in readdir(five_folder) if endswith(file, ".bam")]
    three_folder = "/home/abc/Workspace/ConservedUTRs/threeUTRs/"
    three_files = [joinpath(three_folder, file) for file in readdir(three_folder) if endswith(file, ".bam")]

    five_table = conservation_table(five_files)
    three_table = conservation_table(three_files)
    CSV.write("/home/abc/Workspace/ConservedUTRs/fiveUTR_table.csv", five_table)
    CSV.write("/home/abc/Workspace/ConservedUTRs/threeUTR_table.csv", three_table)
end

#run_constervation_table()

function run_local_alignment()
    fasta = "/home/abc/Workspace/ConservedUTRs/threeUTR.fasta.gz"
    outfolder = "/home/abc/Workspace/ConservedUTRs/alignments"
    genome = "/home/abc/Data/vibrio/genome/NC_016602_3.fna"
    align_local(fasta, outfolder, [genome])
end

run_local_alignment()