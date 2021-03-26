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
    list::IntervalCollection{Float64}
    chroms::Vector{Tuple{String, Int}}
    description::Union{Nothing, String}
end

function Coverage(bigwig_file::String, direction::Symbol; description=nothing)
    @assert direction in [:forward, :reverse]
    reader = open(BigWig.Reader, bigwig_file)
    stran = direction==:forward ? Strand('+') : Strand('-')
    chrlist = BigWig.chromlist(reader)
    intervals = Interval[]
    for record in reader
        push!(intervals, Interval(BigWig.chromid(record), BigWig.chromstart(record), BigWig.chromend(record), stran, BigWig.value(record)))
    end
    close(reader)
    return Coverage(IntervalCollection(intervals, false), chrlist, description)
end

function Coverage(bigwig_forward_file::String, bigwig_reverse_file::String; description=nothing)
    @assert direction in [:forward, :reverse]
    reader_f = open(BigWig.Reader, bigwig_forward_file)
    reader_r = open(BigWig.Reader, bigwig_reverse_file)
    chrlist_f = BigWig.chromlist(reader_f)
    chrlist_r = BigWig.chromlist(reader_r)
    intervals = Interval[]
    for record in reader_f
        push!(intervals, Interval(BigWig.chromid(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_POS, BigWig.value(record)))
    end
    for record in reader_r
        push!(intervals, Interval(BigWig.chromid(record), BigWig.chromstart(record), BigWig.chromend(record), STRNAD_NEG, BigWig.value(record)))
    end
    close(reader_f)
    close(reader_r)
    return Coverage(IntervalCollection(intervals, true), chrlist, description)
end

Base.iterate(coverage::Coverage) = iterate(coverage.list)
Base.iterate(coverage::Coverage, state) =iterate(coverage.list, state)

function values(coverage::Coverage, interval::Interval)
    len = interval.last - interval.first + 1
    vals = zeros(Float64, len)
    offset = interval.first-1
    start = interval.first
    stop = interval.last
    for olint in eachoverlap(coverage.list, interval)
        vals[max(start,olint.first)-offset:min(olint.last, stop)-offset] .+= olint.metadata
    end
    return vals
end

function values(coverage::Coverage)
    vals_f = Dict{String,Vector{Float64}}()
    vals_r = Dict{String,Vector{Float64}}()
    for (chr, len) in coverage.chrlist
        interval_f = Interval(chr, 1, len, STRAND_POS)
        push!(vals_f, values(coverage, interval_f))
        interval_r = Interval(chr, 1, len, STRAND_NEG)
        push!(vals_r, values(coverage, interval_r))
    end
    return vals_f, vals_r
end

Base.push!(coverage::Coverage, interval::Interval) = push!(coverage.list, interval)
function merge!(coverage1::Coverage, coverage2::Coverage) 
    for interval in coverage2
        push!(coverage1, feature)
    end
end

function merge(coverages::Vector{Coverage})
    @assert all(coverages[1].chrlist == c.chrlist for c in coverages[2:end])
    new_intervals = Interval[]
    factor = 1/length(coverages)
    for coverage in coverages
        for interval in coverage
            interval.
            push!(new_intervals, interval)
        end
    end
    return Coverage(IntervalCollection(new_intervals, true), coverages[1].chrlist, join([c.description for c in coverages if !isnothing c.description], ";"))
end