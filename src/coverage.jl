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

function compute_coverage(bam_file::String; norm=1000000, unique_mappings_only=false)
    record::BAM.Record = BAM.Record()
    reader = BAM.Reader(open(bam_file), index=bam_file*".bai")
    chromosome_list = [n for n in zip(bam_chromosome_names(reader), bam_chromosome_lengths(reader))]
    intervals_f = Interval[]
    intervals_r = Interval[]
    count = 0
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        (unique_mappings_only && hasxatag(record)) && continue
        ref = BAM.refname(record)
        left = BAM.position(record)
        right = BAM.rightposition(record)
        count += 1
        BAM.ispositivestrand(record) ? push!(coverage_f, Interval(ref, left, right)) : push!(coverage_r, Interval(ref, left, right))
    end
    close(reader)
    coverage_f = coverage(IntervalCollection(intervals_f))
    coverage_r = coverage(IntervalCollection(intervals_r))
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

function compute_coverage(files::SingleTypeFiles; norm=1000000, unique_mappings_only=false, skip_existing_files=true)
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
    intervals = Vector{Interval{Float32}}()
    for record in reader_f
        push!(intervals, Interval(BigWig.chromid(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_POS, BigWig.value(record)))
    end
    for record in reader_r
        push!(intervals, Interval(BigWig.chromid(record), BigWig.chromstart(record), BigWig.chromend(record), STRNAD_NEG, BigWig.value(record)))
    end
    close(reader_f)
    close(reader_r)
    return Coverage(IntervalCollection(intervals, true), chrlist)
end

function Coverage(paired_files::PairedSingleTypeFiles)
    @assert paired_files.type = ".bw"
    coverages = Vector{Coverage}()
    for (file_forward, file_reverse) in paired_files
        push!(coverages, Coverage(file_forward, file_revers))
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
    for olint in eachoverlap(coverage.list, interval)
        vals[max(start,olint.first)-offset:min(olint.last, stop)-offset] .+= olint.metadata
    end
    return vals
end

function Base.values(coverage::Coverage)
    vals_f = Dict{String,Vector{Float32}}()
    vals_r = Dict{String,Vector{Float32}}()
    for (chr, len) in coverage.chroms
        interval_f = Interval(chr, 1, len, STRAND_POS)
        push!(vals_f, chr=>values(coverage, interval_f))
        interval_r = Interval(chr, 1, len, STRAND_NEG)
        push!(vals_r, chr=>values(coverage, interval_r))
    end
    return vals_f, vals_r
end

Base.push!(coverage::Coverage, interval::Interval) = push!(coverage.list, interval)
function merge!(coverage1::Coverage, coverage2::Coverage) 
    for interval in coverage2
        push!(coverage1, feature)
    end
end

function merge(coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    new_intervals = Vector{Interval{Float32}}()
    for coverage in coverages
        for interval in coverage
            push!(new_intervals, interval)
        end
    end
    return Coverage(IntervalCollection(new_intervals, true), coverages[1].chroms)
end
merge(coverages::Vector{Coverage}) = merge(coverages...)

function correlation(coverages::Coverage ...)
    @assert all(coverages[1].chroms == c.chroms for c in coverages[2:end])
    value_arrays = Dict{String,Matrix{Float32}}(chr => Matrix{Float32}(undef, length(coverages), len) for (chr, len) in coverages[1].chroms)
    for (i,coverage) in enumerate(coverages)  
        for (chr, values) in values(coverage)
            value_arrays[chr][i,!] = values
        end
    end
    correlations = Dict{String,Matrix{Float32}}(chr => Matrix{Float32}(undef, length(coverages), length(coverages)) for (chr, len) in coverages[1].chroms)
    for (chr, arr) in value_arrays
        correlations[chr] = cor(arr, dims=2)
    end
    return correlations
end
correlation(coverages::Vector{Coverage}) = correlation(coverages...)