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

function annotate_utrs!(annotations::DataFrame; max_distance=300)::DataFrame

end
coverage("/home/malte/Workspace/data/test.bam")