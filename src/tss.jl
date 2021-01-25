using XAM
include("utils.jl")
include("io.jl")

function coverage(bam_file::String; norm=1000000, max_misses=2, min_length=25, unique_mappings_only=true)::Dict
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    ref_lengths = bam_chromosome_lengths(reader)
    ref_names = bam_chromosome_names(reader)
    coverage = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    
    count = 0
    while !eof(reader)
        count += 1
        read!(reader, record)
        !BAM.ismapped(record) && continue
        ref = BAM.refname(record)
        pos = BAM.position(record) 
        len = BAM.alignlength(record)
        if unique_mappings_only
            aux = get_XA_tag(BAM.auxdata(record).data)
            (aux == "-") || continue
        end
        coverage[ref][pos:pos+len] .+= 1
    end
    close(reader)
    norm_factor = norm/count
    normalized_coverage = Dict(ref => cov .* norm_factor for (ref,cov) in coverage)
    return normalized_coverage
end

function diff(coverage::Array{Float64,1})::Array{Float64,1}
    d = zeros(length(coverage)-1)
    d[1] = coverage[2] - coverage[1]
    d[2] = maximum(coverage[3] .- coverage[1:2])
    for (i,val) in enumerate(coverage[4:end])
        d[i+2] = maximum(val .- coverage[i-3:i-1])
    end
    return d
end

function tss(coverage_fs::Vector{String}, coverage_rs::Vector{String}, tex_fs::Vector{String}, tex_rs::Vector{String}; 
    max_distance=300, min_step=10, min_ratio=1.1)
    tss = DataFrame()
    for i in 1:length(coverage_fs)
        d_f = {key=>diff(values) for (key, values) in coverage_fs[i]}
        d_r = {key=>diff(values) for (key, values) in coverage_rs[i]}
        tex_d_f = {key=>diff(values) for (key, values) in tex_fs[i]}
        tex_d_r = {key=>diff(values) for (key, values) in tex_rs[i]}
    end
end

function terms(coverage_fs::Vector{String}, coverage_rs::Vector{String}; max_distance=300, min_diff=10)
    d_f = {key=>diff(values) for (key, values) in coverage_f}
    d_r = {key=>diff(values) for (key, values) in coverage_r}
end

function annotate_utrs!(annotations::Dict{String, DataFrame}, tss::Dict{String, DataFrame}, terms::Dict{String, DataFrame}; 
    max_distance=300)::Dict{String, DataFrame}
    
end

#coverage("/home/malte/Workspace/data/test.bam")