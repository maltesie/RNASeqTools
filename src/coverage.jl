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

struct ErrorCoverage <: AnnotationContainer
    ref_seq::LongDNASeq
    chroms::Dict{String, UnitRange}
    fcount::Dict{DNA, Vector{Int}}
    rcount::Dict{DNA, Vector{Int}}
end

function ErrorCoverage(bam_file::String, genome::Genome; is_reverse_complement=false)
    function add_nucleotide(c::Matrix{Int}, base::DNA, i::Int)
        base === DNA_A && (c[1,i]+=1;return)
        base === DNA_T && (c[2,i]+=1;return)
        base === DNA_C && (c[3,i]+=1;return)
        base === DNA_G && (c[4,i]+=1;return)
        base === DNA_Gap && (c[5,i]+=1;return)
    end
    fcount = zeros(Int, 5, length(genome.seq))
    rcount = zeros(Int, 5, length(genome.seq))
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    seq = LongDNASeq(0)
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        hasxatag(record) && continue
        reverse = (isread2(record) != is_reverse_complement)
        ref_name::String = BAM.refname(record)
        starting_pos = first(genome.chroms[ref_name]) - 1
        seq = BAM.sequence(record)
        reverse && reverse_complement!(seq)
        ispositive = BAM.ispositivestrand(record) != reverse
        count = ispositive ? fcount : rcount
        offset, nops = BAM.cigar_position(record)
        current_ref::Int = ispositive ? BAM.leftposition(record) : BAM.rightposition(record)
        current_seq::Int = ispositive ? 1 : length(seq)
        for i in offset:4:offset + (nops - 1) * 4
            x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
            op = BioAlignments.Operation(x & 0x0F)
            n = x >> 4
            (starting_pos + current_ref + (ispositive ? n : 0)) > length(genome.seq) && continue
            r::StepRange{Int} = ispositive ? (0:1:n-1) : (0:-1:-n+1)
            if op == OP_INSERT || op == OP_SOFT_CLIP
                ispositive ? (current_seq += n) : (current_seq -= n)
            elseif BioAlignments.isdeleteop(op)
                for ii in r
                    ref_pos = starting_pos + current_ref + ii
                    add_nucleotide(count, DNA_Gap, ref_pos)
                end
                ispositive ? (current_ref += n) : (current_ref -= n)
            elseif BioAlignments.ismatchop(op)
                for ii in r
                    ref_pos = starting_pos + current_ref + ii
                    seq_pos = current_seq + ii
                    add_nucleotide(count, seq[seq_pos], ref_pos)
                end
                ispositive ? (current_ref += n) : (current_ref -= n)
                ispositive ? (current_seq += n) : (current_seq -= n)
            end
        end
    end
    fdict = Dict{DNA, Vector{Int}}(s=>fcount[i,:] for (i,s) in enumerate((DNA_A, DNA_T, DNA_C, DNA_G, DNA_Gap)))
    rdict = Dict{DNA, Vector{Int}}(s=>rcount[i,:] for (i,s) in enumerate((DNA_A, DNA_T, DNA_C, DNA_G, DNA_Gap)))
    return ErrorCoverage(genome.seq, genome.chroms, fdict, rdict)
end

struct ErrorAnnotation <: AnnotationStyle
    type::String
    name::String
    ref::Vector{Int}
    a::Vector{Int}
    t::Vector{Int}
    g::Vector{Int}
    c::Vector{Int}
    del::Vector{Int}
end

struct ErrorFeatures <: AnnotationContainer
    list::IntervalCollection{ErrorAnnotation}
    chroms::Dict{String, UnitRange{Int}}
end

function ErrorFeatures(gff_file::String, bam_file::String, genome::Genome)
    new_intervals = Interval{ErrorAnnotation}[]
    features = Features(gff_file)
    errorcov = ErrorCoverage(bam_file, genome)
    for feature in features
        left = leftposition(feature)
        right = rightposition(feature)
        seq = genome[refname(feature)][left:right]
        ispositivestrand(feature) || reverse_complement!(seq)
        count = ispositivestrand(feature) ? errorcov.fcount : errorcov.rcount
        r = ispositivestrand(feature) ? (left:right) : (right:-1:left)
        ref = Int[(seq[i] in (DNA_A, DNA_T, DNA_G, DNA_C)) ? count[seq[i]][ii] : 0 for (i, ii) in enumerate(r)]
        annotation = ErrorAnnotation(type(feature), name(feature), ref, count[DNA_A][r], count[DNA_T][r], count[DNA_G][r], count[DNA_C][r], count[DNA_Gap][r]) 
        push!(new_intervals, Interval(refname(feature), left, right, strand(feature), annotation))
    end
    return ErrorFeatures(IntervalCollection(new_intervals, true), genome.chroms)
end

function totalvalues(feature::Interval{ErrorAnnotation})
    return [feature.metadata.a[i]+feature.metadata.t[i]+feature.metadata.g[i]+feature.metadata.c[i]+feature.metadata.del[i] for i in 1:length(feature)]
end

refvalues(feature::Interval{ErrorAnnotation}) = feature.metadata.ref
refpercentage(feature::Interval{ErrorAnnotation}) = refvalues(feature) ./ totalvalues(feature)

function compute_coverage(bam_file::String; norm=1000000, only_unique_alignments=true, is_reverse_complement=false, max_temp_length=500)
    reader = BAM.Reader(open(bam_file), index = bam_file*".bai") # needed for names
    chromosome_list = [n for n in zip(
        bam_chromosome_names(reader), bam_chromosome_lengths(reader)
    )]
    vals_f = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in chromosome_list)
    vals_r = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in chromosome_list)
    count = 0
    for alignment in Alignments(bam_file; only_unique_alignments = only_unique_alignments, is_reverse_complement=is_reverse_complement)
        if ischimeric(alignment; check_annotation = false, min_distance = max_temp_length)
            for part in alignment
                ref = refname(part)
                left = leftposition(part)
                right = rightposition(part)
                ispositivestrand(part) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .+= 1.0)
                count += 1
            end
        else
            ref = refname(alignment[1])
            left = leftposition(alignment)
            right = rightposition(alignment)
            ispositivestrand(alignment) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .+= 1.0)
            count += 1
        end
    end
    norm_factor = norm > 0 ? norm/count : 1.0
    coverage = Coverage(vals_f, vals_r, chromosome_list)
    filename_f = bam_file[1:end-4] * "_forward.bw"
    filename_r = bam_file[1:end-4] * "_reverse.bw"
    writer_f = BigWig.Writer(open(filename_f, "w"), chromosome_list)
    writer_r = BigWig.Writer(open(filename_r, "w"), chromosome_list)
    for interval in coverage
        writer = strand(interval) === STRAND_POS ? writer_f : writer_r
        write(writer,
            (interval.seqname, interval.first,
                interval.last, interval.metadata * norm_factor))
    end
    close(writer_f)
    close(writer_r)
end

function compute_coverage(files::SingleTypeFiles; norm=0, only_unique_alignments=true, overwrite_existing=false, is_reverse_complement=false)
    files.type == ".bam" || throw(AssertionError("Only .bam files accepted for alignments."))
    bw_files = Vector{Tuple{String, String}}()
    for file in files
        filename_f = file[1:end-4] * "_forward.bw"
        filename_r = file[1:end-4] * "_reverse.bw"
        push!(bw_files, (filename_f, filename_r))
        (!overwrite_existing && isfile(filename_f) && isfile(filename_r)) && continue
        compute_coverage(file; norm=norm, only_unique_alignments=only_unique_alignments, is_reverse_complement=is_reverse_complement)
    end
    return PairedSingleTypeFiles(bw_files, ".bw", "_forward", "_reverse")
end

struct Coverage <: AnnotationContainer
    list::IntervalCollection{Float64}
    chroms::Vector{Tuple{String, Int}}
end

function Coverage(bigwig_file::String, direction::Symbol)
    direction in [:forward, :reverse] || throw(Assertion("Valid values for direction are :forward and :reverse"))
    reader = open(BigWig.Reader, bigwig_file)
    stran = direction==:forward ? STRAND_POS : STRAND_NEG
    chrlist = BigWig.chromlist(reader)
    intervals = Vector{Interval{Float64}}()
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

function Base.values(coverage::Coverage, interval::Interval)
    len = interval.last - interval.first + 1
    vals = zeros(Float64, len)
    offset = interval.first-1
    start = interval.first
    stop = interval.last
    for olint in eachoverlap(coverage.list, interval; filter=strand_filter)
        vals[max(start,olint.first)-offset:min(olint.last, stop)-offset] .+= olint.metadata
    end
    return vals
end

function Base.values(coverage::Coverage)
    vals = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    for (chr, len) in coverage.chroms
        interval_f = Interval(chr, 1, len, STRAND_POS)
        interval_r = Interval(chr, 1, len, STRAND_NEG)
        push!(vals, chr=>(values(coverage, interval_f), values(coverage, interval_r)))
    end
    return vals
end

function Base.merge(coverages::Coverage ...)
    all(coverages[1].chroms == c.chroms for c in coverages[2:end]) || throw(Assertion("All Coverage objects have to be defined for the same reference seuqences."))
    n = length(coverages)
    vals_f = Dict(chr=>zeros(Float64, len) for (chr,len) in coverages[1].chroms)
    vals_r = Dict(chr=>zeros(Float64, len) for (chr,len) in coverages[1].chroms)
    for cv in coverages
        for interval in cv
            strand(interval) == STRAND_POS ?
            vals_f[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata/n :
            vals_r[refname(interval)][leftposition(interval):rightposition(interval)] .+= interval.metadata/n
        end
    end
    return Coverage(vals_f, vals_r, coverages[1].chroms)
end
Base.merge(coverages::Vector{Coverage}) = merge(coverages...)

function rolling_sum(a, n::Int)
    (1<=n<=length(a)) || throw(Assertion("n has to be betwen 1 and $(length(a))."))
    out = similar(a, length(a)-n+1)
    out[1] = sum(a[1:n])
    for i in eachindex(out)[2:end]
        out[i] = out[i-1]-a[i-1]+a[i+n-1]
    end
    return out
end

function Base.diff(coverage::Vector{Float64}, invert::Bool, window_size::Int; min_background_ratio=1.2, filter_local_max=true, filter_background_ratio=true)
    d = zeros(Float64,length(coverage))
    min_background_increase = min_background_ratio - 1
    filtered_d = zeros(Float64,length(coverage))
    invert ? d[1:end-1] = @view(coverage[1:end-1]) .- @view(coverage[2:end])  : d[2:end] = @view(coverage[2:end]) .- @view(coverage[1:end-1])
    half_window_size = floor(Int, window_size/2)
    filtered_d[half_window_size:end-half_window_size] = rolling_sum(d,window_size)
    if filter_local_max
        for i in 1:length(d)-window_size+1
            _, mi = findmax(@view(d[i:i+window_size-1]))
            for j in 1:window_size
                mi != j && (d[i+j-1]=0.0)
            end
        end
    end
    filtered_d[findall(x -> x <= 0, d)] .= 0.0
    if filter_background_ratio
        for i in findall(!iszero, filtered_d)
            check_range = invert ? (i:i+half_window_size) : (i-half_window_size+1:i)
            lim = filtered_d[i]/min_background_increase
            any(@view(coverage[check_range]) .< lim) || (filtered_d[i] = 0.0)
        end
    end
    return filtered_d
end

function tsss(notex::Coverage, tex::Coverage; min_tex_ratio=1.3, min_step=10, window_size=10, min_background_ratio=1.2)
    (notex.chroms == tex.chroms) || throw(Assertion("All Coverage objects have to be defined for the same reference sequences."))
    chrs = [chr[1] for chr in notex.chroms]
    vals_notex = values(notex)
    vals_tex = values(tex)
    intervals = Vector{Interval{Float64}}()
    for chr in chrs
        notex_f, notex_r = vals_notex[chr]
        tex_f, tex_r = vals_tex[chr]
        d_forward = diff(tex_f, false, window_size; min_background_ratio=min_background_ratio)
        d_reverse = diff(tex_r, true, window_size; min_background_ratio=min_background_ratio)
        check_d_forward = diff(notex_f, false, window_size; min_background_ratio=min_background_ratio, filter_local_max=false, filter_background_ratio=false)
        check_d_reverse = diff(notex_r, true, window_size; min_background_ratio=min_background_ratio, filter_local_max=false, filter_background_ratio=false)
        check_forward = ((d_forward ./ check_d_forward) .>= min_tex_ratio) .& (d_forward .>= min_step)
        check_reverse = ((d_reverse ./ check_d_reverse) .>= min_tex_ratio) .& (d_reverse .>= min_step)
        for (pos, val) in zip(findall(!iszero, check_forward), abs.(d_forward[check_forward]))
            push!(intervals, Interval(chr, pos, pos, STRAND_POS, val))
        end
        for (pos, val)  in zip(findall(!iszero, check_reverse), abs.(d_reverse[check_reverse]))
            push!(intervals, Interval(chr, pos, pos, STRAND_NEG, val))
        end
    end
    return Coverage(IntervalCollection(intervals, true), tex.chroms)
end

function terms(coverage::Coverage; min_step=10, window_size=10, min_background_ratio=1.2)
    vals = values(coverage)
    intervals = Vector{Interval{Float64}}()
    chrs = [chr[1] for chr in coverage.chroms]
    for chr in chrs
        f, r = vals[chr]
        d_forward = diff(f, true, window_size; min_background_ratio=min_background_ratio)
        d_reverse = diff(r, false, window_size; min_background_ratio=min_background_ratio)
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
