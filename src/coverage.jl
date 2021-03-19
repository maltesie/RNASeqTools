function coverage(bam_file::String; norm=10000000, is_reversed=false, unique_mappings_only=true)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    ref_lengths = bam_chromosome_lengths(reader)
    ref_names = bam_chromosome_names(reader)
    coverage_f = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    coverage_r = Dict(ref=>zeros(Int, len) for (ref, len) in zip(ref_names, ref_lengths))
    count = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        ref = BAM.refname(record)
        pos = BAM.position(record)
        len = BAM.alignlength(record)
        (unique_mappings_only && (get_XA_tag(BAM.auxdata(record).data) != "-")) && continue
        count += 1
        (BAM.ispositivestrand(record) == !is_reversed)  ? (coverage_f[ref][pos:pos+len-1] .+= 1) : (coverage_r[ref][pos:pos+len-1] .-= 1)
    end
    close(reader)
    norm_factor = norm/count
    return Dict(ref => cov .* norm_factor for (ref,cov) in coverage_f), Dict(ref => cov .* norm_factor for (ref,cov) in coverage_r)
end