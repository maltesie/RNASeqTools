function Base.write(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

struct SingleTypeFiles <: FileCollection
    list::Vector{String}
    type::String
end

function SingleTypeFiles(files::Vector{String})
    endings = [fname[findlast('.', fname):end] for fname in files]
    @assert length(unique(endings)) == 1
    SingleTypeFiles(files, endings[1])
end

function SingleTypeFiles(folder::String, type::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if endswith(fname, type)], type)
end

function SingleTypeFiles(folder::String, type::String, prefix::String)
    SingleTypeFiles([joinpath(folder, fname) for fname in readdir(folder) if (endswith(fname, type) && startswith(fname, prefix))], type)
end

Base.length(files::SingleTypeFiles) = length(files.list)

function Base.iterate(files::SingleTypeFiles)
    isempty(files.list) && (return nothing)
    return (files.list[1], 1)
end

function Base.iterate(files::SingleTypeFiles, state::Int)
    state + 1 > length(files.list) && (return nothing)
    return (files.list[state+1], state + 1)
end

function hassingledir(files::SingleTypeFiles)
    return length(unique(dirname(file) for file in files)) == 1
end

function Base.dirname(files::SingleTypeFiles)
    @assert hassingledir(files)
    return dirname(files.list[1])
end

struct PairedSingleTypeFiles <: FileCollection
    list::Vector{Tuple{String,String}}
    type::String
    suffix1::Union{String,Nothing}
    suffix2::Union{String,Nothing}
end

function PairedSingleTypeFiles(files1::Vector{String}, files2::Vector{String})
    endingsa = [fname[findlast(fname, "."):end] for fname in files1]
    endingsb = [fname[findlast(fname, "."):end] for fname in files2]
    @assert (length(unique(endingsa)) == 1) && (unique(endingsa) == unique(endingsb))
    PairedSingleTypeFiles(collect(zip(files1, files2)), endingsa[1], nothing, nothing)
end

function PairedSingleTypeFiles(folder::String, type::String; suffix1="_1", suffix2="_2", prefix=nothing)
    type_files = [joinpath(folder, fname) for fname in readdir(folder) if isnothing(prefix) ? endswith(fname, type) : endswith(fname, type) && startswith(fname, prefix)]
    names1 = [f[1:end-(length(type)+length(suffix1))] for f in type_files if f[end-(length(type)+length(suffix1)-1):end-length(type)] == suffix1]
    names2 = [f[1:end-(length(type)+length(suffix2))] for f in type_files if f[end-(length(type)+length(suffix2)-1):end-length(type)] == suffix2]
    @assert Set(names1) == Set(names2)
    PairedSingleTypeFiles([(joinpath(folder, name * suffix1 * type), joinpath(folder, name * suffix2 * type)) for name in names1], type, suffix1, suffix2)
end

Base.length(files::PairedSingleTypeFiles) = length(files.list)

function Base.iterate(files::PairedSingleTypeFiles)
    isempty(files.list) && (return nothing)
    return (files.list[1], 1)
end

function Base.iterate(files::PairedSingleTypeFiles, state::Int)
    state + 1 > length(files.list) && (return nothing)
    return (files.list[state+1], state + 1)
end

function hassingledir(files::PairedSingleTypeFiles)
    dirs1 = unique([dirname(file[1]) for file in files])
    dirs2 = unique([dirname(file[2]) for file in files])
    return length(dirs1) == 1 && Set(dirs1) == Set(dirs2)
end

function Base.dirname(files::PairedSingleTypeFiles)
    @assert hassingledir(files)
    return dirname(files.list[1][1])
end