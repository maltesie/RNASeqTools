using XAM

function strandint(record::BAM.Record; is_rev=false)
    BAM.ispositivestrand(record) ? strand = 1 : strand = -1
    is_rev ? (return strand * -1) : (return strand)
end

function read_bam(bam_file::String; is_rev=false, nb_reads::Int = -1)
    record::BAM.Record = BAM.Record()
    reader = open(BAM.Reader, bam_file)
    read_names::Array{String, 1} = []
    read_poss::Array{Int, 1} = [] 
    read_chrs::Array{String, 1} = []
    read_auxs::Array{String, 1} = []
    is_rev ? rev_factor = -1 : rev_factor = 1
    aux::String = ""
    c::Int = 0
    cc::Int = 0
    while !eof(reader)
        read!(reader, record)
        !BAM.ismapped(record) && continue
        push!(read_poss, BAM.position(record)*strandint(record, is_rev=is_rev))
        push!(read_chrs, BAM.refname(record))
        push!(read_names, BAM.tempname(record))
        "XA" in Set{String}(l for (l::String,) in BAM.auxdata(record)) ? 
            (aux = BAM.auxdata(record)["XA"]; cc+=1) : 
            aux = "-"
        push!(read_auxs, aux)
        c += 1
        ((nb_reads > 0) & (c >= nb_reads)) && break 
    end
    close(reader)
    return read_names, read_poss, read_chrs, read_auxs
end

function read_annotations(gff_file::String, annotation_types = ["CDS"])::Dict{String, DataFrame}
    
    open(gff_file) do reader
        gff = read(reader, String)
    end
    results::Dict{String, DataFrame} = Dict()

    for line in split(gff, "\n")
        startswith(line, "#") | isempty(line) && continue
        chr, _, typ, start, stop, _, strand, _, aux = split(line, "\t")
        typ in annotation_types || continue
        if strand == "-"
            start_int = parse(Int, stop) * -1
            stop_int = parse(Int, start) * -1
        else
            start_int = parse(Int, start)
            stop_int = parse(Int, stop)
        end
        name = split(aux, ";")[1][6:end]
        row = DataFrame(name=name, start=start_int, stop=stop_int)
        chr in keys(results) ? append!(results[chr], row) : results[chr] = row
    end

    for (chr, data) in results
        DataFrames.sort!(data, :start)
        utr5 = DataFrame(name="U5." .* data.name, start=data.start .- 100, stop=copy(data.start))
        utr3 = DataFrame(name="U3." .* data.name, start=copy(data.stop), stop=data.stop .+ 100)
        igr = DataFrame(name="IG<" .* data.name[1:end-1] .* "|" .* data.name[2:end] .* ">", 
            start=data.stop[1:end-1] .+ 100, stop=data.start[2:end] .- 100)
        igr= igr[igr.start .< igr.stop, :]
        DataFrames.append!(data, utr5)
        DataFrames.append!(data, utr3)
        DataFrames.append!(data, igr)
        DataFrames.sort!(data, :start)
    end
    
    return results
end

function read_rrna(gff_file::String; rrna_types=["rRNA", "tRNA"])::Dict{String, DataFrame}

    open(gff_file) do reader
        gff = read(reader, String)
    end
    results::Dict{String, DataFrame} = Dict()
    factor::Int = 1

    for line in split(gff, "\n")
        startswith(line, "#") | isempty(line) && continue
        chr, _, typ, start, stop, _, strand, _, _ = split(line, "\t")
        chr in keys(chr_trans) || (chr_trans[chr] = length(chr_trans)+1)
        typ in rrna_types || continue
        if strand == "-"
            start_int = parse(Int, stop) * -1
            stop_int = parse(Int, start) * -1
        else
            start_int = parse(Int, start)
            stop_int = parse(Int, stop)
        end
        name = split(aux, ";")[1][6:end]
        row = DataFrame(name=name, start=start_int, stop=stop_int)
        chr in keys(results) ? append!(results[chr], row) : results[chr] = row
    end
    
    return results
end