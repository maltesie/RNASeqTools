function bam_chromosome_lengths(reader::BAM.Reader)
    chr_lengths = Int[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_lengths, parse(Int, meta["LN"]))
    end
    return chr_lengths
end

function bam_chromosome_names(reader::BAM.Reader)
    chr_names = String[]
    for meta in findall(BAM.header(reader), "SQ")
        push!(chr_names, meta["SN"])
    end
    return chr_names
end

function compute_coverage(bam_file::String; norm=1000000, unique_mappings_only=true)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    chromosome_list = [n for n in zip(bam_chromosome_names(reader), bam_chromosome_lengths(reader))]
    intervals_f = Vector{Interval{Nothing}}()
    intervals_r = Vector{Interval{Nothing}}()
    count = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        (unique_mappings_only && hasxatag(record)) && continue
        ref = BAM.refname(record)
        left = BAM.position(record)
        right = BAM.rightposition(record)
        count += 1
        BAM.ispositivestrand(record) ? push!(intervals_f, Interval(ref, left, right)) : push!(intervals_r, Interval(ref, left, right))
    end
    close(reader)
    coverage_f = coverage(IntervalCollection(intervals_f, true))
    coverage_r = coverage(IntervalCollection(intervals_r, true))
    norm_factor = norm/count
    filename_f = bam_file[1:end-4] * "_forward.bw"
    filename_r = bam_file[1:end-4] * "_reverse.bw"
    writer_f = BigWig.Writer(open(filename_f, "w"), chromosome_list)
    writer_r = BigWig.Writer(open(filename_r, "w"), chromosome_list)
    for interval in coverage_f
        write(writer_f, (interval.seqname, interval.first, interval.last, interval.metadata*norm_factor))
    end
    for interval in coverage_r
        write(writer_r, (interval.seqname, interval.first, interval.last, interval.metadata*norm_factor))
    end
    close(writer_f)
    close(writer_r)
end

function compute_coverage(files::SingleTypeFiles; norm=1000000, unique_mappings_only=true, skip_existing_files=true)
    @assert files.type == ".bam"
    for file in files
        filename_f = file[1:end-4] * "_forward.bw"
        filename_r = file[1:end-4] * "_reverse.bw"
        (skip_existing_files && isfile(filename_f) && isfile(filename_r)) && continue
        compute_coverage(file; norm=norm, unique_mappings_only=unique_mappings_only)
    end
end

struct Coverage <: AnnotationContainer
    list::IntervalCollection{Float32}
    chroms::Vector{Tuple{String, Int}}
end

function Coverage(bigwig_file::String, direction::Symbol)
    @assert direction in [:forward, :reverse]
    reader = open(BigWig.Reader, bigwig_file)
    stran = direction==:forward ? Strand('+') : Strand('-')
    chrlist = BigWig.chromlist(reader)
    intervals = Vector{Interval{Float32}}()
    for record in reader
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), stran, BigWig.value(record)))
    end
    close(reader)
    return Coverage(IntervalCollection(intervals, true), chrlist)
end

function Coverage(bigwig_forward_file::String, bigwig_reverse_file::String)
    reader_f = open(BigWig.Reader, bigwig_forward_file)
    reader_r = open(BigWig.Reader, bigwig_reverse_file)
    chrlist_f = BigWig.chromlist(reader_f)
    chrlist_r = BigWig.chromlist(reader_r)
    @assert chrlist_f == chrlist_r
    intervals = Vector{Interval{Float32}}()
    for record in reader_f
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_POS, BigWig.value(record)))
    end
    for record in reader_r
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_NEG, BigWig.value(record)))
    end
    close(reader_f)
    close(reader_r)
    return Coverage(IntervalCollection(intervals, true), chrlist_f)
end

function Coverage(values_f::CoverageValues, values_r::CoverageValues, chroms::Vector{Tuple{String, Int}})
    new_intervals = Vector{Interval{Float32}}()
    for vals in (values_f, values_r)
        current_pos = 1
        current_end = 1
        current_refid = 1
        current_ref = first(chroms[current_refid])
        current_len = last(chroms[current_refid])
        while true
            if current_end >= current_len
                current_refid += 1
                if current_refid > length(chroms)
                    break
                end
                current_ref = first(chroms[current_refid])
                current_len = last(chroms[current_refid])
                current_pos = 1
                current_end = 1
            end
            current_value = vals[current_ref][current_pos]
            while vals[current_ref][current_end+1] == current_value
                current_end += 1
                current_end == current_len && break
            end
            (current_pos == 1 && current_end == current_len) || push!(new_intervals, Interval(current_ref, current_pos, current_end, 
                                                                    vals===vals_f ? STRAND_POS : STRAND_NEG, current_value/nb_coverages))
            current_pos = current_end + 1
            current_end = current_pos
        end
    end
    return Coverage(IntervalCollection(new_intervals, true), chroms)
end

function Coverage(paired_files::PairedSingleTypeFiles)
    @assert paired_files.type == ".bw"
    coverages = Vector{Coverage}()
    for (file_forward, file_reverse) in paired_files
        push!(coverages, Coverage(file_forward, file_reverse))
    end
    return merge(coverages...)
end

Base.iterate(coverage::Coverage) = iterate(coverage.list)
Base.iterate(coverage::Coverage, state) =iterate(coverage.list, state)

function Base.write(filename_f::String, filename_r::String, coverage::Coverage)
    @assert endswith(filename_f, ".bw") && endswith(filename_f, ".bw")
    writer_f = BigWig.Writer(open(filename_f, "w"), chromosome_list)
    writer_r = BigWig.Writer(open(filename_r, "w"), chromosome_list)
    for interval in coverage
        strand(interval) == STRAND_POS ? 
        write(writer_f, (interval.seqname, interval.first, interval.last, interval.metadata)) :
        write(writer_r, (interval.seqname, interval.first, interval.last, interval.metadata))
    end
    close(writer_f)
    close(writer_r)
end
function Base.write(filename::String, coverage::Coverage, gff_type::String)
    writer = GFF3.Writer(open(filename, "w"))
    for interval in coverage
        line = "$(interval.seqname)\t.\t$gff_type\t$(interval.first)\t$(interval.last)\t.\t$(interval.strand)\t.\tValue=$(interval.metadata)"
        record = GFF3.Record(line)
        write(writer, record)
    end
    close(writer)
end

function Base.values(coverage::Coverage, interval::Interval)
    len = interval.last - interval.first + 1
    vals = zeros(Float32, len)
    offset = interval.first-1
    start = interval.first
    stop = interval.last
    for olint in eachoverlap(coverage.list, interval; filter=strand_filter)
        vals[max(start,olint.first)-offset:min(olint.last, stop)-offset] .+= olint.metadata
    end
    return vals
end

function Base.values(coverage::Coverage)
    vals = Dict{String,Tuple{Vector{Float32},Vector{Float32}}}()
    for (chr, len) in coverage.chroms
        interval_f = Interval(chr, 1, len, STRAND_POS)
        interval_r = Interval(chr, 1, len, STRAND_NEG)
        push!(vals, chr=>(values(coverage, interval_f), values(coverage, interval_r)))
    end
    return vals
end

function Base.merge(coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    nb_coverages = length(coverages)
    vals_f = Dict(chr=>zeros(Float32, len) for (chr,len) in coverages[1].chroms)
    vals_r = Dict(chr=>zeros(Float32, len) for (chr,len) in coverages[1].chroms)
    for coverage in coverages
        for interval in coverage
            strand(interval) == STRAND_POS ? 
            (vals_f[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata; has_pos=true) :
            (vals_r[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata; has_neg=true)
        end
    end
    return Coverage(vals_f, vals_r, coverage[1].chroms)
end
Base.merge(coverages::Vector{Coverage}) = merge(coverages...)

function differential(coverage1::Coverage, coverage2::Coverage)
    @assert coverage1.chroms == coverage2.chroms
    vals1 = values(coverage1)
    vals2 = values(coverage2)
    return Coverage(Dict(key=>vals1[key] ./ vals2[key] for key in keys(vals1)))
end

function correlation(coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    value_arrays = Dict{String,Matrix{Float32}}(chr => Matrix{Float32}(undef, length(coverages), len) for (chr, len) in coverages[1].chroms)
    for (i,coverage) in enumerate(coverages)  
        for (chr, (values_f, values_r)) in values(coverage)
            value_arrays[chr][i,1:end] = values_f .+ values_r
        end
    end

    correlations = Dict{String,Matrix{Float32}}(chr => Matrix{Float32}(undef, length(coverages), length(coverages)) for (chr, len) in coverages[1].chroms)
    for (chr, arr) in value_arrays
        correlations[chr] = cor(arr, dims=2)
    end
    return correlations
end
correlation(coverages::Vector{Coverage}) = correlation(coverages...)

function mincorrelation(coverages::Coverage ...)
    corr = correlation(coverages)
    min_corr = 1.0
    for (chr, matrix) in corr
        for i in first(size(matrix)), j in 1:i
            min_corr = min(matrix[i,j], min_corr)
        end
    end
    return min_corr
end
mincorrelation(coverages::Vector{Coverage}) = mincorrelation(coverages...)

function rolling_sum(a, n::Int)
    @assert 1<=n<=length(a)
    out = similar(a, length(a)-n+1)
    out[1] = sum(a[1:n])
    for i in eachindex(out)[2:end]
        out[i] = out[i-1]-a[i-1]+a[i+n-1]
    end
    return out
end

function diff(coverage::Vector{Float32}, invert::Bool, window_size::Int, min_background_increase::Float64)
    d = zeros(Float32,length(coverage))
    filtered_d = zeros(Float32,length(coverage))
    invert ? d[1:end-1] = @view(coverage[1:end-1]) .- @view(coverage[2:end])  : d[2:end] = @view(coverage[2:end]) .- @view(coverage[1:end-1])
    half_window_size = floor(Int, window_size/2)
    filtered_d[half_window_size:end-half_window_size] = rolling_sum(d,window_size)
    for i in 1:length(d)-window_size+1
        mv, mi = findmax(@view(d[i:i+window_size-1]))
        for j in 1:window_size
            mi != j && (d[i+j-1]=0.0)
        end
    end
    filtered_d[findall(iszero, d)] .= 0.0
    for i in findall(!iszero, filtered_d)
        check_range = invert ? (i:i+half_window_size) : (i-half_window_size+1:i)
        lim = filtered_d[i]/min_background_increase
        any(@view(coverage[check_range]) .< lim) || (filtered_d[i] = 0.0)
    end
    return filtered_d
end

function tsss(notex::Coverage, tex::Coverage; min_step=10, min_ratio=1.3, min_background_increase=0.2, window_size=10)
    @assert notex.chroms == tex.chroms
    chrs = [chr[1] for chr in notex.chroms]
    vals_notex = values(notex)
    vals_tex = values(tex)
    intervals = Vector{Interval{Float32}}()
    for chr in chrs
        notex_f, notex_r = vals_notex[chr]
        tex_f, tex_r = vals_tex[chr]
        d_forward = diff(tex_f, false, window_size, min_background_increase)
        d_reverse = diff(tex_r, true, window_size, min_background_increase)
        check_forward = ((tex_f ./ notex_f) .>= min_ratio) .& (d_forward .>= min_step)
        check_reverse = ((tex_r ./ notex_r) .>= min_ratio) .& (d_reverse .>= min_step)
        for (pos, val) in zip(findall(!iszero, check_forward), abs.(d_forward[check_forward]))
            push!(intervals, Interval(chr, pos, pos, STRAND_POS, val))
        end
        for (pos, val)  in zip(findall(!iszero, check_reverse), abs.(d_reverse[check_reverse]))
            push!(intervals, Interval(chr, pos, pos, STRAND_NEG, val))
        end
    end
    return Coverage(IntervalCollection(intervals, true), tex.chroms)
end

function terms(coverage::Coverage; min_step=10, window_size=10, min_background_increase=0.1)
    vals = values(coverage)
    intervals = Vector{Interval{Float32}}()
    chrs = [chr[1] for chr in coverage.chroms]
    for chr in chrs
        f, r = vals[chr]
        d_forward = diff(f, true, window_size, min_background_increase)
        d_reverse = diff(r, false, window_size, min_background_increase)
        check_forward = d_forward .>= min_step
        check_reverse = d_reverse .>= min_step
        for (pos, val) in zip(findall(!iszero, check_forward), abs.(d_forward[check_forward]))
            push!(intervals, Interval(chr, pos, pos, STRAND_POS, val))
        end
        for (pos, val)  in zip(findall(!iszero, check_reverse), abs.(d_reverse[check_reverse]))
            push!(intervals, Interval(chr, pos, pos, STRAND_NEG, val))
        end
    end
    return Coverage(IntervalCollection(intervals, true), coverage.chroms)
end

