function Base.write(filename::String, content::String)
    open(filename, "w") do f
        write(f, content)
    end
end

function commonstart(strs::Vector{String})
    for c in 1:minimum(length(s) for s in strs)
        all(s[c] == strs[1][c] for s in strs) || return strs[1][1:c-1]
    end
    return sortperm(strs)[1]
end

function commonprefix(files::Vector{String}; delimiter='_')
    common_start = commonstart(files)
    i = findfirst(delimiter, common_start)
    return isnothing(i) ? common_start : common_start[1:i]
end

function commonend(strs::Vector{String})
    for c in 0:minimum(length(s) for s in strs)-1
        all(s[end-c] == strs[1][end-c] for s in strs) || (return strs[1][end-c+1:end])
    end
    return strs[sortperm(strs)[1]]
end

function commonfiletype(files::Vector{String}; delimiter='.')
    common_end = commonend(files)
    i = findfirst(Char(delimiter), common_end)
    return isnothing(i) ? "" : common_end[i:end]
end

function commonsuffix(files::Vector{String}; delimiter_suffix='_', delimiter_filetype='.')
    common_filetype_len = length(commonfiletype(files; delimiter=delimiter_filetype))
    common_end = commonend([f[1:end-common_filetype_len] for f in files])
    i = findlast(delimiter_suffix, common_end)
    return isnothing(i) ? common_end : common_end[i:end]
end

function groups(files::Vector{String}; delimiter_prefix='_', delimiter_suffix='_')
    stretches = String[]
    offset = length(commonprefix(files; delimiter=delimiter_prefix)) + 1
    stripped_strs = [s[offset:end] for s in files]
    for stripped_str in stripped_strs
        check_strs = copy(stripped_strs)
        for (i, c) in enumerate(stripped_str)
            check_index = [((i <= length(check_str)) && (c === check_str[i])) for check_str in check_strs]
            if sum(check_index) == 1
                push!(stretches, stripped_str[1:i-(stripped_str[i-1] == delimiter_suffix ? 2 : 1)])
                break
            end
            check_strs = check_strs[check_index]
        end
    end
    common_end = commonend(stretches)
    i = findlast(delimiter_suffix, common_end)
    common_end_len = isnothing(i) ? length(common_end) : length(common_end)-i+1
    conds = [s[1:end-common_end_len] for s in stretches]
    return Dict(c => collect(findall(c .== conds)) for c in unique(conds))
end

function SingleTypeFiles(files::Vector{String})
    ending = commonfiletype(files)
    isempty(ending) && throw(AssertionError("Could not detect common file ending."))
    SingleTypeFiles(files, ending)
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
    ending = commonfiletype(vcat(files1, files2))
    isempty(ending) && throw(AssertionError("Could not detect common file ending."))
    suffix1 = commonsuffix(files1)
    suffix2 = commonsuffix(files2)
    PairedSingleTypeFiles(collect(zip(files1, files2)), ending, suffix1, suffix2)
end

function PairedSingleTypeFiles(list::Vector{Tuple{String,String}})
    PairedSingleTypeFiles([f[1] for f in list], [f[2] for f in list])
end

function PairedSingleTypeFiles(folder::String, type::String; suffix1="_1", suffix2="_2", prefix=nothing)
    folder = abspath(folder)
    type_files = [joinpath(folder, fname) for fname in readdir(folder) if isnothing(prefix) ? endswith(fname, type) : endswith(fname, type) && startswith(fname, prefix)]
    names1 = [f[1:end-(length(type)+length(suffix1))] for f in type_files if f[end-(length(type)+length(suffix1)-1):end-length(type)] == suffix1]
    names2 = [f[1:end-(length(type)+length(suffix2))] for f in type_files if f[end-(length(type)+length(suffix2)-1):end-length(type)] == suffix2]
    @assert Set(names1) == Set(names2)
    PairedSingleTypeFiles([(joinpath(folder, name * suffix1 * type), joinpath(folder, name * suffix2 * type)) for name in names1], type, suffix1, suffix2)
end

groupfiles(files::SingleTypeFiles) = groups([basename(f) for f in files.list])
function groupfiles(files::PairedSingleTypeFiles)
    g1 = groups([basename(f[1]) for f in files.list])
    g2 = groups([basename(f[2]) for f in files.list])
    g1 == g2 || throw(AssertionError("Groups made from first and second file in pairs do not match!"))
    return g1
end

function filesexist(files::SingleTypeFiles)
    for file in files
        isfile(file) || throw(AssertionError("$file does not exist!"))
    end
end

function filesexist(files::PairedSingleTypeFiles)
    for (file1, file2) in files
        isfile(file1) || throw(AssertionError("$file1 does not exist!"))
        isfile(file2) || throw(AssertionError("$file2 does not exist!"))
    end
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

function Base.show(files::SingleTypeFiles)
    s = "SingleTypeFiles with file type $(files.type) and $(length(files)) entries:\n\n"
    for (i,file) in enumerate(files)
        s *= "$i:\t$file\n"
    end
    print(s)
end

function Base.show(files::PairedSingleTypeFiles)
    s = "PairedSingleTypeFiles with file type $(files.type) and $(length(files)) entries:\n\n"
    for (i,(file1, file2)) in enumerate(files)
        s *= "$(i):\t$file1\n\t$file2\n"
    end
    print(s)
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
CoverageFiles(folder::String; suffix_forward="_forward", suffix_reverse="_reverse", prefix=nothing) =
    PairedSingleTypeFiles(folder, ".bw"; suffix1=suffix_forward, suffix2=suffix_reverse, prefix=prefix)