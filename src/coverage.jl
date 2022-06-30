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

function BaseCoverage(bam_file::String, genome::Genome; include_secondary_alignments=true, is_reverse_complement=false,
                                                        only_positive_strand=false, quality_cut = 0x01)
    fcount = Dict(chr=>zeros(Int, 7, length(genome.chroms[chr])) for chr in keys(genome.chroms))
    rcount = Dict(chr=>zeros(Int, 7, length(genome.chroms[chr])) for chr in keys(genome.chroms))
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    seq = LongDNASeq(0)
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        hasxatag(record) && continue
        !isprimary(record) && !include_secondary_alignments && continue
        is_reverse = is_reverse_complement #(isread2(record) != is_reverse_complement)
        ref_name::String = BAM.refname(record)
        seq = BAM.sequence(record)
        qual = BAM.quality(record)
        length(seq) > 0 || continue
        is_positive = (BAM.ispositivestrand(record) != is_reverse) || only_positive_strand
        is_reverse && reverse_complement!(seq)
        is_positive || complement!(seq)
        count = is_positive ? fcount : rcount
        offset, nops = BAM.cigar_position(record)
        current_ref::Int = BAM.leftposition(record)
        current_seq::Int = 0
        for data_index in offset:4:offset + (nops - 1) * 4
            x = unsafe_load(Ptr{UInt32}(pointer(record.data, data_index)))
            op = BioAlignments.Operation(x & 0x0F)
            n = x >> 4
            r::UnitRange{Int} = current_ref:current_ref+n-1
            if op === OP_INSERT
                current_seq += n
                count[ref_name][7, current_ref] += 1
            elseif op === OP_SOFT_CLIP
                current_seq += n
            elseif BioAlignments.isdeleteop(op)
                for ref_pos in r
                    count[ref_name][5, ref_pos] += 1
                end
                current_ref += n
            elseif BioAlignments.ismatchop(op)
                for (ii, ref_pos) in enumerate(r)
                    if qual[current_seq + ii] >= quality_cut
                        base = seq[current_seq + ii]
                        base === DNA_A ? count[ref_name][1, ref_pos] += 1 :
                        base === DNA_T ? count[ref_name][2, ref_pos] += 1 :
                        base === DNA_G ? count[ref_name][3, ref_pos] += 1 :
                        base === DNA_C ? count[ref_name][4, ref_pos] += 1 :
                        count[ref_name][6, ref_pos] += 1
                    else
                        count[ref_name][6, ref_pos] += 1
                    end
                end
                current_ref += n
                current_seq += n
            end
        end
    end
    fdict = Dict(chr=>Dict(dna=>fcount[chr][i, :] for (i, dna) in enumerate((:A, :T, :G, :C, :Gap, :N, :Ins))) for chr in keys(genome.chroms))
    rdict = Dict(chr=>Dict(dna=>rcount[chr][i, :] for (i, dna) in enumerate((:A, :T, :G, :C, :Gap, :N, :Ins))) for chr in keys(genome.chroms))
    return BaseCoverage(genome, fdict, rdict)
end

function mismatchfractions(base_coverage::BaseCoverage, from::Symbol, to::Symbol; direction=:both)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    d = Dict(:A=>DNA_A, :T=>DNA_T, :G=>DNA_G, :C=>DNA_C, :Gap=>DNA_Gap, :N=>DNA_N)
    base_from = d[from]
    direction in (:forward, :reverse, :both) || throw(AssertionError("direction has to be :forward, :reverse or :both"))
    fstats = Dict(chr=>Float64[] for chr in keys(base_coverage.genome.chroms))
    rstats = Dict(chr=>Float64[] for chr in keys(base_coverage.genome.chroms))
    for chr in keys(base_coverage.genome.chroms)
        findex = base_coverage.genome[chr] .=== base_from
        rindex = BioSequences.complement(base_coverage.genome[chr]) .=== base_from
        f = base_coverage.fcount[chr][to][findex] ./ sum(base_coverage.fcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))[findex]
        r = base_coverage.rcount[chr][to][rindex] ./ sum(base_coverage.rcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))[rindex]
        append!(fstats[chr], filter!(x->(!isnan(x) && !iszero(x)), f))
        append!(rstats[chr], filter!(x->(!isnan(x) && !iszero(x)), r))
    end
    return fstats, rstats
end

function mismatchpositions(base_coverage::BaseCoverage, from::Symbol, to::Symbol; ratio_cut=0.9)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    d = Dict(:A=>DNA_A, :T=>DNA_T, :G=>DNA_G, :C=>DNA_C, :Gap=>DNA_Gap, :N=>DNA_N)
    base_from = d[from]
    pos = Interval{Tuple{Int, Float64}}[]
    for chr in keys(base_coverage.genome.chroms)
        index = base_coverage.genome[chr] .=== base_from
        ratios = base_coverage.fcount[chr][to] ./ sum(base_coverage.fcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))
        ratio_index = (ratios .> ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .=== true))
            mc = base_coverage.fcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_POS, (mc, round(r; digits=3))))
        end
        index = BioSequences.complement(base_coverage.genome[chr]) .=== base_from
        ratios = base_coverage.rcount[chr][to] ./ sum(base_coverage.rcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))
        ratio_index = (ratios .> ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .=== true))
            mc = base_coverage.rcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_NEG, (mc, round(r; digits=3))))
        end
    end
    return pos
end

function mismatchpositions(base_coverage::BaseCoverage; check=(:A, :T, :G, :C), ratio_cut=0.9)
    pos = Interval{Annotation}[]
    for from in check, to in check
        from === to && continue
        mpos = mismatchpositions(base_coverage, from, to; ratio_cut=ratio_cut)
        append!(pos, [Interval(refname(i), leftposition(i), leftposition(i), strand(i),
                        Annotation("SNP", "$from->$to", Dict("mutation"=>"$from->$to",
                                                            "count"=>"$(i.metadata[1])",
                                                            "frequency"=>"$(i.metadata[2])"))) for i in mpos])
    end
    return Features(pos)
end

function deletionpositions(base_coverage::BaseCoverage; check=(:A, :T, :G, :C), ratio_cut=0.9)
    pos = Interval{Annotation}[]
    for from in check
        mpos = mismatchpositions(base_coverage, from, :Gap; ratio_cut=ratio_cut)
        append!(pos, [Interval(refname(i), leftposition(i), leftposition(i), strand(i),
                        Annotation("DEL", "$from", Dict("base"=>"$from",
                                                        "count"=>"$(i.metadata[1])",
                                                        "frequency"=>"$(i.metadata[2])"))) for i in mpos])
    end
    merged_pos = Interval{Annotation}[]
    merger = Interval{Annotation}[]
    for check in (true, false)
        for pos_interval in pos
            ispositivestrand(pos_interval) == check || continue
            show(pos_interval)
            if !isempty(merger) && (leftposition(pos_interval) != (rightposition(merger[end])+1))
                push!(merged_pos, length(merger) == 1 ? merger[1] : Interval(refname(merger[end]),
                                                                                leftposition(merger[1]),
                                                                                rightposition(merger[end]),
                                                                                strand(merger[end]),
                                                                                Annotation("DEL", join(name(m) for m in merger),
                                                                                            Dict("base"=>join(name(m) for m in merger),
                                                                                                "count"=>"$(mean(param(m, "count", Int) for m in merger))",
                                                                                                "frequency"=>"$(mean(param(m, "frequency", Float64) for m in merger))"))))
                empty!(merger)
            end
            push!(merger, pos_interval)
        end
    end
    return Features(merged_pos)
end

all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
function mismatchcontexthist(base_coverage::BaseCoverage, from::Symbol, to::Symbol; pm=2, bins=20, ratio_cut=0.0)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    kmer_histograms = Dict(LongDNASeq(c)=>zeros(Int, bins) for c in all_perm([DNA_A, DNA_T, DNA_G, DNA_C], 2*pm+1))
    mpos = mismatchpositions(base_coverage, from, to; ratio_cut)
    for mismatch_interval in mpos
        p = leftposition(mismatch_interval)
        freq_index = Int(floor(mismatch_interval.metadata[2] * (bins-1))) + 1
        ref = refname(mismatch_interval)
        kmer = base_coverage.genome[ref][p-pm:p+pm]
        kmer_histograms[kmer][freq_index] += 1
    end
    return kmer_histograms
end

function coverage(feature::Interval{BaseAnnotation})
    return feature.metadata.a .+ feature.metadata.t .+ feature.metadata.g .+ feature.metadata.c .+ feature.metadata.gap
end

refcount(feature::Interval{BaseAnnotation}) = feature.metadata.ref

function compute_coverage(bam_file::String; norm=1000000, include_secondary_alignments=false,
                                            is_reverse_complement=false, max_temp_length=500,
                                            suffix_forward="_forward", suffix_reverse="_reverse")

    filename_f = bam_file[1:end-4] * suffix_forward * ".bw"
    filename_r = bam_file[1:end-4] * suffix_reverse * ".bw"
    reader = BAM.Reader(open(bam_file), index = bam_file*".bai")
    chromosome_list = [n for n in zip(
        bam_chromosome_names(reader), bam_chromosome_lengths(reader)
    )]
    vals_f = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in chromosome_list)
    vals_r = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in chromosome_list)
    count = 0
    for alignment in Alignments(bam_file; include_secondary_alignments=include_secondary_alignments,
                                            is_reverse_complement=is_reverse_complement)
        if ischimeric(alignment; check_annotation = false, min_distance = max_temp_length)
            for part in alignment
                ref = refname(part)
                left = leftposition(part)
                right = rightposition(part)
                ispositivestrand(part) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .-= 1.0)
                count += 1
            end
        else
            ref = refname(alignment[1])
            left = leftposition(alignment)
            right = rightposition(alignment)
            ispositivestrand(alignment) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .-= 1.0)
            count += 1
        end
    end
    norm_factor = norm > 0 ? norm/count : 1.0
    coverage = Coverage(vals_f, vals_r, chromosome_list)

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

function compute_coverage(files::SingleTypeFiles;
            norm=1000000, include_secondary_alignments=true, overwrite_existing=false, is_reverse_complement=false)
    files.type == ".bam" || throw(AssertionError("Only .bam files accepted for alignments."))
    bw_files = Vector{Tuple{String, String}}()
    for file in files
        filename_f = file[1:end-4] * "_forward.bw"
        filename_r = file[1:end-4] * "_reverse.bw"
        push!(bw_files, (filename_f, filename_r))
        (!overwrite_existing && isfile(filename_f) && isfile(filename_r)) && continue
        compute_coverage(file;
            norm=norm, include_secondary_alignments=include_secondary_alignments, is_reverse_complement=is_reverse_complement)
    end
    return PairedSingleTypeFiles(bw_files, ".bw", "_forward", "_reverse")
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

function Coverage(bc::BaseCoverage)
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

function Base.merge(coverages::Vector{Coverage})
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

function rolling_sum(a::Vector{Float64}, n::Int; circular=true)
    (1<=n<=length(a)) || throw(Assertion("n has to be betwen 1 and $(length(a))."))
    out = similar(a)
    out[1] = sum(a[1:n])
    for i in 2:length(a)-n+1
        out[i] = out[i-1]-a[i-1]+a[i+n-1]
    end
    for (ii,i) in enumerate(length(a)-n+2:length(a))
        out[i] = out[i-1]-a[i-1]+ (circular ? a[ii] : 0.0)
    end
    return circshift(out, Int((n-1)/2))
end

function maxdiffsignalpositions(coverage::Vector{Float64}; max_ppk=5, compute_max_within=5, circular=true)
    all(coverage .>= 0) || all(coverage .<= 0) || throw(AssertionError("Coverage values have to be positive for + strand and negative for - strand!"))
    is_negative_strand = all(coverage .<= 0)
    max_nb_peaks = length(coverage) / 1000 * max_ppk
    d = zeros(Float64,length(coverage))
    d[is_negative_strand ? (1:length(d)-1) : (2:length(d))] = @view(coverage[2:end]) .- @view(coverage[1:end-1])
    if circular
        d[is_negative_strand ? length(d) : 1] = coverage[1]-coverage[end]
        d = vcat(d[end-compute_max_within+1:end], d, d[1:compute_max_within])
    end
    peak_index = [i-compute_max_within for i in compute_max_within+1:length(d)-compute_max_within if all((d[i] .> d[i-compute_max_within:i-1]) .& (d[i] .> d[i+1:i+compute_max_within]))]
    d = d[compute_max_within+1:end-compute_max_within]
    for cut in 1:Int(floor(maximum(d)))
        length(peak_index) < max_nb_peaks && break
        peak_index = peak_index[d[peak_index] .>= cut]
    end
    return peak_index
end
function maxdiffsignalpositions(coverages::Vector{Coverage}; max_ppk=5, compute_max_within=3, circular=true, invert=false)
    merged_coverage_values = values(merge(coverages); invert=invert)
    Dict{String, Tuple{Vector{Int}, Vector{Int}}}(
        chr=>(
            maxdiffsignalpositions(merged_coverage_values[chr][1]; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular),
            maxdiffsignalpositions(merged_coverage_values[chr][2]; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular)
        )
        for chr in sort(collect(keys(merged_coverage_values)))
    )
end

function background_plateau_pairs(coverage::Vector{Float64}, peak_index::Vector{Int}; compute_step_within=3, circular=true)
    all(coverage .>= 0) || all(coverage .<= 0) || throw(AssertionError("Coverage values have to be positive for + strand and negative for - strand!"))
    is_negative_strand = all(coverage .<= 0)
    pairs = Matrix{Int}(undef, length(peak_index), 2)
    mycoverage = circular ? vcat(coverage[end-compute_step_within+1:end], coverage, coverage[1:compute_step_within]) : coverage
    top, bottom = is_negative_strand ? (minimum, maximum) : (maximum, minimum)
    for (ii, index) in enumerate(peak_index)
        circular && (index += compute_step_within)
        pairs[ii, 1] = bottom(@view(mycoverage[(index-compute_step_within):(index+compute_step_within)]))
        pairs[ii, 2] = top(@view(mycoverage[(index-compute_step_within):(index+compute_step_within)]))
    end
    is_negative_strand && (pairs .*= -1)
    return pairs
end
function background_plateau_pairs(s_ls::Vector{Coverage}, s_is::Vector{Coverage}, peaks::Dict{String, Tuple{Vector{Int}, Vector{Int}}};
                                    compute_step_within=1, circular=true)
    valss = vcat([values(s_l) for s_l in s_ls], [values(s_i) for s_i in s_is])
    hcat(
        [vcat(
            [vcat(
                background_plateau_pairs(vals[chr][1], peaks[chr][1]; compute_step_within=compute_step_within, circular=circular),
                background_plateau_pairs(vals[chr][2], peaks[chr][2]; compute_step_within=compute_step_within, circular=circular)
            )
            for chr in sort(collect(keys(peaks)))]...
        )
        for vals in valss]...
    )
end

function normalize_counts_between!(counts_between::Counts, nb_samples1::Int, nb_samples2::Int)
    avg_sample_notex::Vector{Float64} = [gmean(counts_between.values[i, 1:nb_samples1]) for i in 1:length(counts_between)]
    logdiffs_notex = [log.(counts_between.values[:, i]) .- log.(avg_sample_notex) for i in 1:nb_samples1]
    counts_between.values[:, 1:nb_samples1] ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_notex])'

    avg_sample_tex::Vector{Float64} = [gmean(counts_between.values[i, nb_samples1+1:nb_samples1+nb_samples2]) for i in 1:length(counts_between)]
    logdiffs_tex = [log.(counts_between.values[:, i]) .- log.(avg_sample_tex) for i in nb_samples1+1:nb_samples1+nb_samples2]
    counts_between.values[:, nb_samples1+1:nb_samples1+nb_samples2] ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_tex])'

    avg_log_diffs = log.(avg_sample_notex) .- log.(avg_sample_tex)
    h = fit(Histogram, avg_log_diffs; nbins=100)
    max_index = argmax(h.weights)
    max_value = (h.edges[1][max_index]+h.edges[1][max_index+1])/2.0
    counts_between.values[:, nb_samples1+1:nb_samples1+nb_samples2] ./= exp(max_value)
end

function normalized_notex_tex_increase(notex_background::Matrix{Float64}, notex_plateau::Matrix{Float64}, tex_background::Matrix{Float64}, tex_plateau::Matrix{Float64})
    avg_notex_plateau::Vector{Float64} = [gmean(notex_plateau[i,:]) for i in 1:size(notex_plateau)[1]]
    logdiffs_notex = [log.(notex_plateau[:, i]) .- log.(avg_notex_plateau) for i in 1:size(notex_plateau)[2]]
    notex_background ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_notex])'
    notex_plateau ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_notex])'

    avg_tex_plateau::Vector{Float64} = [gmean(tex_plateau[i,:]) for i in 1:size(tex_plateau)[1]]
    logdiffs_notex = [log.(tex_plateau[:, i]) .- log.(avg_tex_plateau) for i in 1:size(tex_plateau)[2]]
    tex_background ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_notex])'
    tex_plateau ./= exp.([median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs_notex])'

    avg_log_diffs = log.(avg_notex_plateau) .- log.(avg_tex_plateau)
    h = fit(Histogram, avg_log_diffs; nbins=100)
    max_index = argmax(h.weights)
    max_value = (h.edges[1][max_index]+h.edges[1][max_index+1])/2.0
    tex_plateau ./= exp(max_value)
    tex_background ./= exp(max_value)
    return Counts(Dict("notex_increase"=>1:size(notex_plateau)[2], "tex_increase"=>size(notex_plateau)[2]+1:size(notex_plateau)[2]+size(tex_plateau)[2]),
                    hcat(notex_plateau .- notex_background, tex_plateau .- tex_background) .+ 1)
end

function normalize_counts_within!(counts_within::Counts, nb_samples::Int)
    avg_sample::Vector{Float64} = [gmean(counts_within.values[i, nb_samples+1:end]) for i in 1:length(counts_within)]
    logdiffs = [log2.(counts_within.values[:, nb_samples+i]) .- log2.(avg_sample) for i in 1:nb_samples]
    for (i, nf) in enumerate(2 .^ [median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs])
        counts_within.values[:, [i,i+nb_samples]] ./= nf
    end
end

function transcriptionalstartsites(notexs::Vector{Coverage}, texs::Vector{Coverage};
                fdr_tex_increase=0.05, fdr_background_step=0.05, min_step_fc=1.2, min_step_reads=10, compute_step_within=1, compute_max_within=5, max_ppk=5, circular=true, source="drna_seq")
    chroms = vcat([s.chroms for s in notexs], [s.chroms for s in texs])
    all(chroms[1] == chrom for chrom in chroms) || throw(Assertion("All coverages have to be defined on the same reference sequences."))

    peaks = maxdiffsignalpositions(texs; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular)
    background_plateau_matrix = background_plateau_pairs(notexs, texs, peaks; compute_step_within=compute_step_within, circular=circular) .+ 1

    notex_background_index = collect(1:2:2*length(notexs))
    notex_plateau_index = notex_background_index .+ 1
    tex_background_index = collect(2*length(notexs)+1:2:2*(length(notexs)+length(texs)))
    tex_plateau_index = tex_background_index .+ 1

    counts_notex_tex_increase = normalized_notex_tex_increase(Float64.(background_plateau_matrix[:, notex_background_index]),
                                                                Float64.(background_plateau_matrix[:, notex_plateau_index]),
                                                                Float64.(background_plateau_matrix[:, tex_background_index]),
                                                                Float64.(background_plateau_matrix[:, tex_plateau_index]))

    counts_tex_step = Counts(Dict("tex_background"=>1:length(texs),
                                    "tex_plateau"=>length(texs)+1:(2*length(texs))),
                                    background_plateau_matrix[:,vcat(tex_background_index,tex_plateau_index)])
    normalize_counts_within!(counts_tex_step, length(texs))


    (_, _, adjp_between) = dge_glm(counts_notex_tex_increase, "notex_increase", "tex_increase"; tail=:right)
    (basevals, fc, adjp_within) = dge_glm(counts_tex_step, "tex_background", "tex_plateau"; tail=:right)

    offset = 0
    intervals = Vector{Interval{Annotation}}()
    for chr in sort(collect(keys(peaks)))
        for (cp, s) in zip(peaks[chr], (STRAND_POS, STRAND_NEG))
            append!(intervals, [Interval(chr, p, p, s,
                                    Annotation("TSS", "",
                                        Dict("fdr_background_step"=>"$(round(adjp_within[offset+i], digits=8))",
                                            "fdr_tex_increase"=>"$(round(adjp_between[offset+i], digits=8))",
                                            "avg_tex_plateau"=>"$(round(basevals[offset+i] * 2^fc[offset+i], digits=2))",
                                            "avg_tex_fc"=>"$(round(2^fc[offset+i], digits=2))",
                                            "avg_tex_increase"=>"$(round(basevals[offset+i] * 2^fc[offset+i] - basevals[offset+i], digits=2))",
                                            "source"=>source))) for (i,p) in enumerate(cp)
                                                if ((adjp_between[offset+i]<=fdr_tex_increase) &&
                                                    (adjp_within[offset+i]<=fdr_background_step) &&
                                                    (2^fc[offset+i] > min_step_fc) &&
                                                    ((basevals[offset+i] * 2^fc[offset+i] - basevals[offset+i]) > min_step_reads))])
            offset += length(cp)
        end
    end
    return Features(intervals)
end

all_fisher_test_combinations(s_ls::Vector{Tuple{Int,Int}}, s_hs::Vector{Tuple{Int,Int}}) =
    vec([pvalue(FisherExactTest(nt...,t...); tail=:right) for t in s_hs, nt in s_ls])
function all_fisher_test_combinations(s_ls::Matrix{Int}, s_hs::Matrix{Int})
    s_hs_rows, s_hs_samples = size(s_hs)
    s_ls_rows, s_ls_samples = size(s_ls)
    s_hs_rows == s_ls_rows || throw(AssertionError("s_hs and s_ls have to have the same number of rows!"))
    s_hs_samples % 2 == 0 && s_ls_samples % 2 == 0 || throw(AssertionError("s_hs and s_ls have to have even amount of columns!"))
    s_ls_pairs = [[(i,j) for (i,j) in partition(r, 2)] for r in eachrow(s_ls)]
    s_hs_pairs = [[(i,j) for (i,j) in partition(r, 2)] for r in eachrow(s_hs)]
    all_fisher_test_combinations.(s_ls_pairs, s_hs_pairs)
end

fisher_combined_pvalue(pvalues::Vector{Float64}) = 1-cdf(Chisq(2*length(pvalues)), -2*sum(log.(pvalues)))

test_multiple_fisher_combined(s_ls::Matrix{Int}, s_hs::Matrix{Int}) =
    adjust(PValues(fisher_combined_pvalue.(all_fisher_test_combinations(s_ls, s_hs))), BenjaminiHochberg())

function terminationsites(bcms::Vector{Coverage}, nobcms::Vector{Coverage};
    fdr_bcm_increase=0.05, fdr_background_step=0.05, min_step_fc=1.2, min_step_reads=10, compute_step_within=1, compute_max_within=5, max_ppk=5, circular=true, source="term_seq")

    chroms = vcat([s.chroms for s in bcms], [s.chroms for s in nobcms])
    all(chroms[1] == chrom for chrom in chroms) || throw(Assertion("All coverages have to be defined on the same reference sequences."))

    peaks = maxdiffsignalpositions(nobcms; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular, invert=true)
    background_plateau_matrix = background_plateau_pairs(bcms, nobcms, peaks; compute_step_within=compute_step_within, circular=circular) .+ 1

    nobcm_background_index = collect(2*length(bcms)+1:2:2*(length(bcms)+length(nobcms)))
    nobcm_plateau_index = nobcm_background_index .+ 1

    counts_within = Counts(Dict("nobcm_background"=>1:length(nobcms),
                                "nobcm_plateau"=>length(nobcms)+1:(2*length(nobcms))),
                            background_plateau_matrix[:,vcat(nobcm_background_index,nobcm_plateau_index)])
    normalize_counts_within!(counts_within, length(nobcms))

    (basevals, fc, adjp_background_step) = dge_glm(counts_within, "nobcm_background", "nobcm_plateau"; tail=:right)

    adjp_bcm_increase = ones(length(adjp_background_step))
    adjp_bcm_increase[adjp_background_step .<= fdr_background_step] =
        test_multiple_fisher_combined(background_plateau_matrix[:, 1:2*length(bcms)][adjp_background_step .<= fdr_background_step, :],
        background_plateau_matrix[:, 2*length(bcms)+1:2*(length(bcms)+length(nobcms))][adjp_background_step .<= fdr_background_step, :])

    offset = 0
    intervals = Vector{Interval{Annotation}}()
    for chr in sort(collect(keys(peaks)))
        for (cp, s) in zip(peaks[chr], (STRAND_POS, STRAND_NEG))
            append!(intervals, [Interval(chr, p, p, s,
                                    Annotation("TERM", "",
                                        Dict("fdr_background_step"=>"$(round(adjp_background_step[offset+i], digits=8))",
                                            "fdr_bcm_increase"=>"$(round(adjp_bcm_increase[offset+i], digits=8))",
                                            "rho_dependent"=> (adjp_bcm_increase[offset+i] <= fdr_bcm_increase) ? "true" : "false",
                                            "avg_nobcm_plateau"=>"$(round(basevals[offset+i] * 2^fc[offset+i], digits=2))",
                                            "avg_nobcm_fc"=>"$(round(2^fc[offset+i], digits=2))",
                                            "avg_nobcm_increase"=>"$(round(basevals[offset+i] * 2^fc[offset+i] - basevals[offset+i], digits=2))",
                                            "source"=>source))) for (i,p) in enumerate(cp)
                                                if ((adjp_background_step[offset+i]<=fdr_background_step) &&
                                                    (2^fc[offset+i] > min_step_fc) &&
                                                    ((basevals[offset+i] * 2^fc[offset+i] - basevals[offset+i]) > min_step_reads))])
            offset += length(cp)
        end
    end
    return Features(intervals)
end
