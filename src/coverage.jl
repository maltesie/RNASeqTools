function compute_coverage(bam_file::String; norm=1000000, include_secondary_alignments=false,
                                            is_reverse_complement=false, max_temp_length=1000, include_chimeric=false,
                                            only_chimeric=false, overwrite_existing=false,
                                            suffix_forward="_forward", suffix_reverse="_reverse")

    filename_f = bam_file[1:end-4] * suffix_forward * ".bw"
    filename_r = bam_file[1:end-4] * suffix_reverse * ".bw"

    (isfile(filename_f) && isfile(filename_r)) && !overwrite_existing && return (filename_f, filename_r)

    alns = Alignments(bam_file; include_secondary_alignments=include_secondary_alignments, is_reverse_complement=is_reverse_complement)
    coverage = Coverage(alns; norm=norm, include_chimeric=include_chimeric, only_chimeric=only_chimeric, max_temp_length=max_temp_length)

    writer_f = BigWig.Writer(open(filename_f, "w"), alns.chroms)
    writer_r = BigWig.Writer(open(filename_r, "w"), alns.chroms)
    for interval in coverage
        writer = strand(interval) === STRAND_POS ? writer_f : writer_r
        write(writer, (interval.seqname, interval.first, interval.last, interval.metadata))
    end
    close(writer_f)
    close(writer_r)
    return (filename_f, filename_r)
end

function compute_coverage(files::SingleTypeFiles; norm=1000000, include_secondary_alignments=true, overwrite_existing=false,
                                                    is_reverse_complement=false, max_temp_length=1000, include_chimeric=false,
                                                    only_chimeric=false, suffix_forward="_forward", suffix_reverse="_reverse")
    files.type == ".bam" || throw(AssertionError("Only .bam files accepted for alignments."))
    bw_files = [compute_coverage(file; norm=norm, include_secondary_alignments=include_secondary_alignments, overwrite_existing=overwrite_existing,
                                        is_reverse_complement=is_reverse_complement, max_temp_length=max_temp_length,
                                        include_chimeric=include_chimeric, only_chimeric=only_chimeric,
                                        suffix_forward=suffix_forward, suffix_reverse=suffix_reverse)
                                        for file in files]
    return PairedSingleTypeFiles(bw_files, ".bw", "_forward", "_reverse")
end

function Coverage(alignments::Alignments; norm=1000000, max_temp_length=1000, include_chimeric=false, only_chimeric=false)
    vals_f = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in alignments.chroms)
    vals_r = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in alignments.chroms)
    for alignment in alignments
        if ischimeric(alignment; check_annotation = false, min_distance = max_temp_length)
            include_chimeric || only_chimeric || continue
            for part in alignment
                ref = refname(part)
                left = leftposition(part)
                right = rightposition(part)
                ispositivestrand(part) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .-= 1.0)
            end
        elseif !only_chimeric
            ref = refname(alignment[1])
            left = leftposition(alignment)
            right = rightposition(alignment)
            ispositivestrand(alignment) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .-= 1.0)
        end
    end
    Coverage(vals_f, vals_r, alignments.chroms) * (norm > 0 ? norm/length(alignments) : 1.0)
end

function Coverage(bigwig_file::String, direction::Symbol)
    direction in [:forward, :reverse] || throw(Assertion("Valid values for direction are :forward and :reverse"))
    reader = open(BigWig.Reader, bigwig_file)
    stran = direction==:forward ? STRAND_POS : STRAND_NEG
    chrlist = BigWig.chromlist(reader)
    intervals = Vector{Interval{Float64}}()
    for record in reader
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), stran, convert(Float64, BigWig.value(record))))
    end
    close(reader)
    return Coverage(IntervalCollection(intervals, true), chrlist)
end

function Coverage(bigwig_forward_file::String, bigwig_reverse_file::String)
    reader_f = open(BigWig.Reader, bigwig_forward_file)
    reader_r = open(BigWig.Reader, bigwig_reverse_file)
    chrlist_f = BigWig.chromlist(reader_f)
    chrlist_r = BigWig.chromlist(reader_r)
    (chrlist_f == chrlist_r) || throw(Assertion("Both strand .bw files have to be defined for the same reference sequences."))
    intervals = Vector{Interval{Float64}}()
    for record in reader_f
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_POS, convert(Float64, BigWig.value(record))))
    end
    for record in reader_r
        push!(intervals, Interval(BigWig.chrom(record), BigWig.chromstart(record), BigWig.chromend(record), STRAND_NEG, convert(Float64, BigWig.value(record))))
    end
    close(reader_f)
    close(reader_r)
    return Coverage(IntervalCollection(intervals, true), chrlist_f)
end

function Coverage(values_f::CoverageValues, values_r::CoverageValues, chroms::Vector{Tuple{String, Int}})
    new_intervals = Vector{Interval{Float64}}()
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
                                                                    vals===values_f ? STRAND_POS : STRAND_NEG, current_value))
            current_pos = current_end + 1
            current_end = current_pos
        end
    end
    return Coverage(IntervalCollection(new_intervals, true), chroms)
end

function Coverage(paired_files::PairedSingleTypeFiles)
    (paired_files.type == ".bw") || throw(Assertion("Only .bw files are accepted for coverage."))
    coverages = [Coverage(file_forward, file_reverse) for (file_forward, file_reverse) in paired_files]
    return merge(coverages...)
end

function Base.:*(coverage::Coverage, factor::Float64)
    vals = values(coverage)
    Coverage(Dict(chr=>(val[1] .* factor) for (chr,val) in vals), Dict(chr=>(val[2] .* factor) for (chr,val) in vals), coverage.chroms)
end

value(interval::Interval{Float64}) = interval.metadata
Base.length(coverage::Coverage) = length(coverage.list)
Base.iterate(coverage::Coverage) = iterate(coverage.list)
Base.iterate(coverage::Coverage, state) =iterate(coverage.list, state)

function Base.write(filename_f::String, filename_r::String, coverage::Coverage)
    (endswith(filename_f, ".bw") && endswith(filename_f, ".bw")) || throw(Assertion("Both files have to end in .bw\n1: $filename_f\n2: $filename_r"))
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

function Base.values(coverage::Coverage, interval::Interval; invert=false)
    len = interval.last - interval.first + 1
    vals = zeros(Float64, len)
    offset = interval.first-1
    start = interval.first
    stop = interval.last
    for olint in eachoverlap(coverage.list, interval; filter=strand_filter)
        vals[max(start,olint.first)-offset:min(olint.last, stop)-offset] .+= olint.metadata
    end
    return invert ? vals .* -1.0 : vals
end

function Base.values(coverage::Coverage; invert=false)
    vals = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    for (chr, len) in coverage.chroms
        interval_f = Interval(chr, 1, len, STRAND_POS)
        interval_r = Interval(chr, 1, len, STRAND_NEG)
        push!(vals, chr=>(values(coverage, interval_f; invert=invert), values(coverage, interval_r; invert=invert)))
    end
    return vals
end

function Base.merge(coverages::Vector{Coverage}; average=true)
    all(coverages[1].chroms == c.chroms for c in coverages[2:end]) || throw(Assertion("All Coverage objects have to be defined for the same reference seuqences."))
    vals_f = Dict(chr=>zeros(Float64, len) for (chr,len) in coverages[1].chroms)
    vals_r = Dict(chr=>zeros(Float64, len) for (chr,len) in coverages[1].chroms)
    for cv in coverages
        for interval in cv
            strand(interval) == STRAND_POS ?
            vals_f[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata :
            vals_r[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata
        end
    end
    merged = Coverage(vals_f, vals_r, coverages[1].chroms)
    average && (merged *= 1/length(coverages))
    return merged
end

function correlation(coverages::Coverage ...)
    vals = [values(c) for c in coverages]
    Dict(chr=>(cor(hcat([v[chr][1] for v in vals]...)), cor(hcat([v[chr][2] for v in vals]...))) for chr in keys(vals[1]))
end
correlation(coverages::Vector{Coverage}) = correlation(coverages...)

function localmaxdiffindex(coverage::Vector{Float64}; rev=false, min_diff=2, min_ratio=1.1, compute_within=5, circular=true)
    all(coverage .>= 0) || all(coverage .<= 0) || throw(AssertionError("Coverage values have to be positive for + strand and negative for - strand!"))
    rev && (coverage = reverse(coverage))
    circular && (coverage = vcat(coverage[end-compute_within+1:end], coverage, coverage[1:compute_within]))
    is_negative_strand = all(coverage .<= 0)
    d = zeros(Float64, length(coverage))
    d[is_negative_strand ? (1:length(d)-1) : (2:length(d))] = @view(coverage[2:end]) .- @view(coverage[1:end-1])
    peak_index = zeros(Bool, length(d))
    mima = zeros(Float64, length(d), 2)
    for i in compute_within+1:length(d)-compute_within
        ((d[i] >= min_diff) && (d[i] === maximum(view(d, i-compute_within:i+compute_within)))) || continue
        mi, ma = (minimum(view(coverage, i-compute_within:i-1)), maximum(view(coverage, i+1:i+compute_within)))
        peak_index[i] = is_negative_strand ? mi/ma >= min_ratio : ma/mi >= min_ratio
        mima[i, 1] = is_negative_strand ? -ma : mi
        mima[i, 2] = is_negative_strand ? -mi : ma
    end
    circular && (peak_index = peak_index[compute_within+1:end-compute_within])
    circular && (mima = mima[compute_within+1:end-compute_within, :])
    mima = mima[peak_index, :]
    rev && reverse!(peak_index)
    rev && reverse!(mima; dims=1)
    return peak_index, mima
end

function maxdiffpositions(coverage::Coverage; type="DIFF", rev=false, min_diff=2, min_ratio=1.1, compute_within=5, circular=true)
    coverage_values = values(coverage)
    features = Interval{Annotation}[]
    for chr in keys(coverage_values)
        findex, _ = localmaxdiffindex(coverage_values[chr][1]; rev=rev, min_diff=min_diff, min_ratio=min_ratio, 
                                                            compute_within=compute_within, circular=circular)
        rindex, _ = localmaxdiffindex(coverage_values[chr][2]; rev=rev, min_diff=min_diff, min_ratio=min_ratio, 
                                                            compute_within=compute_within, circular=circular)
        for (ii, i) in enumerate(findall(findex))
            push!(features, Interval(chr, i, i, STRAND_POS, Annotation(type, "forward_$ii")))
        end
        for (ii, i) in enumerate(findall(rindex))
            push!(features, Interval(chr, i, i, STRAND_NEG, Annotation(type, "reverse_$ii")))
        end
    end
    return Features(features)
end