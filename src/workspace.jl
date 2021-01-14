using DataFrames
using CSV
using XLSX
using Statistics

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

hmm1 = "/home/malte/Workspace/data/HMMposterior_probabilities_chr1.txt"
hmm2 = "/home/malte/Workspace/data/HMMposterior_probabilities_chr2.txt"
ann = "/home/malte/Workspace/data/vibrio_srna.csv"
convert_hmm_to_table(hmm1, hmm2, ann, "../data/hmm_output.xlsx")