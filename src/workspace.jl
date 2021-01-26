using DataFrames
using CSV
using XLSX
using Statistics
using XAM
using BioAlignments

include("tss.jl")
include("preprocess.jl")

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

function align(sequence_file::String, reference_file::String)

end

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

function test_utr_annotation(coverage_files::Array{String, 2})

end

function run_preprocess()
    read_files = ["/home/malte/Workspace/data/vibrio/drnaseq/tex_01_1.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/tex_01_2.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_1.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/tex_20_2.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_1.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/notex_01_2.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_1.fastq.gz",
            "/home/malte/Workspace/data/vibrio/drnaseq/notex_20_2.fastq.gz"]

    output_folder = "/home/malte/Workspace/data/vibrio/drnaseq/"

    fastp_bin_path = "/home/malte/Tools/fastp"

    trim_fastp(read_files, output_folder; fastp_bin=fastp_bin_path)
end

run_preprocess()