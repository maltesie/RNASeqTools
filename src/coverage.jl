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

Base.push!(coverage::Coverage, interval::Interval) = push!(coverage.list, interval)
function Base.merge(coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    new_intervals = Vector{Interval{Float32}}()
    for coverage in coverages
        for interval in coverage
            push!(new_intervals, interval)
        end
    end
    return Coverage(IntervalCollection(new_intervals, true), coverages[1].chroms)
end
Base.merge(coverages::Vector{Coverage}) = merge(coverages...)

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

function diff(coverage::Vector{Float32})
    d = zeros(Float32,length(coverage))
    d[1] = coverage[1]
    d[2:end] = coverage[2:end] - coverage[1:end-1]
    return d
end

function tss(notex::Coverage, tex::Coverage; min_step=10, min_ratio=1.3)
    @assert notex.chroms == tex.chroms
    chrs = [chr[1] for chr in notex.chroms]
    vals_notex = values(notex)
    vals_tex = values(tex)
    intervals = Vector{Interval{Float32}}()
    for chr in chrs
        notex_f = first(vals_notex[chr])
        notex_r = last(vals_notex[chr])
        tex_f = first(vals_tex[chr])
        tex_r = last(vals_tex[chr])
        d_forward = diff(tex_f)
        d_reverse = diff(tex_r)
        check_forward = circshift(((tex_f ./ notex_f) .>= min_ratio), 1) .& (d_forward .>= min_step)
        check_reverse = circshift(((tex_r ./ notex_r) .>= min_ratio), 1) .& (d_reverse .>= min_step)
        for (pos, val) in zip(findall(!iszero, check_forward), abs.(d_forward[check_forward]))
            push!(intervals, Interval(chr, pos, pos, STRAND_POS, val))
        end
        for (pos, val)  in zip(findall(!iszero, check_reverse), abs.(d_reverse[check_reverse]))
            push!(intervals, Interval(chr, pos, pos, STRAND_NEG, val))
        end
    end
    return Coverage(IntervalCollection(intervals, true), tex.chroms)
end

function terms(coverage::Coverage; min_step=10)
    f, r = values(coverage)
    intervals = Vector{Interval{Float32}}()
    chrs = [chr[1] for chr in coverage.chroms]
    for chr in chrs
        d_forward = diff(f[chr])
        d_reverse = diff(r[chr])
        check_forward = d_forward .<= -min_step
        check_reverse = d_reverse .<= -min_step
        for (pos, val) in zip(findall(!iszero, check_forward), abs.(d_forward[check_forward]))
            push!(intervals, Interval(chr, pos, pos, STRAND_POS, val))
        end
        for (pos, val)  in zip(findall(!iszero, check_reverse), abs.(d_reverse[check_reverse]))
            push!(intervals, Interval(chr, pos, pos, STRAND_NEG, val))
        end
    end
    return Coverage(IntervalCollection(intervals, true), coverage.chroms)
end