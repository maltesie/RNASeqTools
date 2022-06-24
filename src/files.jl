function Base.write(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

function SingleTypeFiles(files::Vector{String})
    endings = [fname[findlast('.', fname):end] for fname in files]
    length(unique(endings)) <= 1
    SingleTypeFiles(files, isempty(endings) ? "" : endings[1])
end

function SingleTypeFiles(folder::String, type::String; prefix=nothing)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if isnothing(prefix) ? endswith(fname, type) : endswith(fname, type) && startswith(fname, prefix)], type)
end

function Base.:*(filesa::SingleTypeFiles, filesb::SingleTypeFiles)
    @assert type(filesa) == type(filesb)
    return SingleTypeFiles(vcat(filesa.list, filesb.list), type(filesa))
end

function hassingledir(files::SingleTypeFiles)
    return length(unique(dirname(file) for file in files)) == 1
end

function Base.dirname(files::SingleTypeFiles)
    @assert hassingledir(files)
    return dirname(files.list[1])
end

function PairedSingleTypeFiles(files1::Vector{String}, files2::Vector{String})
    endingsa = unique([fname[findlast('.', fname):end] for fname in files1])
    endingsb = unique([fname[findlast('.', fname):end] for fname in files2])
    @assert (length(endingsa) == 1) && (endingsa == endingsb)
    typ = endingsa[1]
    PairedSingleTypeFiles(collect(zip(files1, files2)), typ, "", "")
end

function PairedSingleTypeFiles(list::Vector{Tuple{String,String}})
    endingsa = unique([fname[1][findlast('.', fname[1]):end] for fname in list])
    endingsb = unique([fname[2][findlast('.', fname[2]):end] for fname in list])
    @assert (length(endingsa) == 1) && (endingsa == endingsb)
    typ = endingsa[1]
    PairedSingleTypeFiles(list, typ, "", "")
end

function PairedSingleTypeFiles(folder::String, type::String; suffix1="_1", suffix2="_2", prefix=nothing)
    folder = abspath(folder)
    type_files = [joinpath(folder, fname) for fname in readdir(folder) if isnothing(prefix) ? endswith(fname, type) : endswith(fname, type) && startswith(fname, prefix)]
    names1 = [f[1:end-(length(type)+length(suffix1))] for f in type_files if f[end-(length(type)+length(suffix1)-1):end-length(type)] == suffix1]
    names2 = [f[1:end-(length(type)+length(suffix2))] for f in type_files if f[end-(length(type)+length(suffix2)-1):end-length(type)] == suffix2]
    @assert Set(names1) == Set(names2)
    PairedSingleTypeFiles([(joinpath(folder, name * suffix1 * type), joinpath(folder, name * suffix2 * type)) for name in names1], type, suffix1, suffix2)
end

type(files::T) where {T<:FileCollection} = files.type
Base.length(files::T) where {T<:FileCollection} = length(files.list)
Base.iterate(files::T) where {T<:FileCollection} = iterate(files.list)
Base.iterate(files::T, state::Int) where {T<:FileCollection} = iterate(files.list, state)
Base.copy(files::SingleTypeFiles) = SingleTypeFiles(files.list, files.type)
Base.copy(files::PairedSingleTypeFiles) = PairedSingleTypeFiles(files.list, files.type, files.suffix1, files.suffix2)
Base.getindex(files::T, i::Int) where {T<:FileCollection} = files.list[i]
Base.getindex(files::SingleTypeFiles, u::UnitRange{Int}) = SingleTypeFiles(files.list[u], files.type)
Base.getindex(files::SingleTypeFiles, index::Vector{Int}) = SingleTypeFiles(files.list[index], files.type)
Base.getindex(files::PairedSingleTypeFiles, u::UnitRange{Int}) = PairedSingleTypeFiles(files.list[u], files.type, files.suffix1, files.suffix2)
Base.getindex(files::PairedSingleTypeFiles, index::Vector{Int}) = PairedSingleTypeFiles(files.list[index], files.type, files.suffix1, files.suffix2)
function Base.:*(filesa::PairedSingleTypeFiles, filesb::PairedSingleTypeFiles)
    @assert type(filesa) == type(filesb)
    suffix1 = filesa.suffix1 == filesb.suffix1 ? filesa.suffix1 : ""
    suffix2 = filesa.suffix2 == filesb.suffix2 ? filesa.suffix2 : ""
    return PairedSingleTypeFiles(vcat(filesa.list, filesb.list), type(filesa), suffix1, suffix2)
end

function hassingledir(files::PairedSingleTypeFiles)
    dirs1 = unique([dirname(file[1]) for file in files])
    dirs2 = unique([dirname(file[2]) for file in files])
    return length(dirs1) == 1 && Set(dirs1) == Set(dirs2)
end

function Base.dirname(files::PairedSingleTypeFiles)
    hassingledir(files) || throw(AssertionError("Files have to be in a single folder."))
    return dirname(files.list[1][1])
end

function Base.write(fname::String, files::SingleTypeFiles)
    files.type in (".csv",) || throw(AssertionError("File type has to be .csv"))
    if files.type == ".csv"
        tables = Vector{Tuple{String,Vector{Any},Vector{String}}}()
        for file in files
            sheetname = basename(file)[1:end-length(files.type)]
            dataframe = DataFrame(CSV.File(file; stringtype=String))
            push!(tables,(sheetname,collect(eachcol(dataframe)), names(dataframe)))
        end
        XLSX.writetable(fname, tables; overwrite=true)
    end
end

function Base.show(files::SingleTypeFiles)
    s = "SingleTypeFiles with file type $(files.type) and $(length(files)) entries:\n\n"
    for (i,file) in enumerate(files)
        s *= "$i:\t$file\n"
    end
    println(s)
end

function Base.show(files::PairedSingleTypeFiles)
    s = "PairedSingleTypeFiles with file type $(files.type) and $(length(files)) entries:\n\n"
    for (i,(file1, file2)) in enumerate(files)
        s *= "$(i):\t$file1\n\t$file2\n"
    end
    println(s)
end

CsvFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".csv"; prefix=prefix)
FastqFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fastq"; prefix=prefix)
FastqgzFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fastq.gz"; prefix=prefix)
FastaFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fasta"; prefix=prefix)
FastagzFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fasta.gz"; prefix=prefix)
GenomeFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fna"; prefix=prefix)
FaFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".fa"; prefix=prefix)
BamFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".bam"; prefix=prefix)
GffFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".gff"; prefix=prefix)
GraphFiles(folder::String; prefix=nothing) = SingleTypeFiles(folder, ".jld2"; prefix=prefix)
PairedFastqFiles(folder::String; suffix1="_1", suffix2="_2", prefix=nothing) = PairedSingleTypeFiles(folder, ".fastq"; suffix1=suffix1, suffix2=suffix2, prefix=prefix)
PairedFastqgzFiles(folder::String; suffix1="_1", suffix2="_2", prefix=nothing) = PairedSingleTypeFiles(folder, ".fastq.gz"; suffix1=suffix1, suffix2=suffix2, prefix=prefix)
PairedFastaFiles(folder::String; suffix1="_1", suffix2="_2", prefix=nothing) = PairedSingleTypeFiles(folder, ".fasta"; suffix1=suffix1, suffix2=suffix2, prefix=prefix)
PairedFastagzFiles(folder::String; suffix1="_1", suffix2="_2", prefix=nothing) = PairedSingleTypeFiles(folder, ".fasta.gz"; suffix1=suffix1, suffix2=suffix2, prefix=prefix)
CoverageFiles(folder::String; suffix1="_forward", suffix2="_reverse", prefix=nothing) = PairedSingleTypeFiles(folder, ".bw"; suffix1=suffix1, suffix2=suffix2, prefix=prefix)