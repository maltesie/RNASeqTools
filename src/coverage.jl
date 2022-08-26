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
    seq = LongDNA{4}(0)
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
    kmer_histograms = Dict(LongDNA{4}(c)=>zeros(Int, bins) for c in all_perm([DNA_A, DNA_T, DNA_G, DNA_C], 2*pm+1))
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
                                            is_reverse_complement=false, max_temp_length=1000, include_chimeric=false,
                                            overwrite_existing=false, suffix_forward="_forward", suffix_reverse="_reverse")

    filename_f = bam_file[1:end-4] * suffix_forward * ".bw"
    filename_r = bam_file[1:end-4] * suffix_reverse * ".bw"

    (isfile(filename_f) && isfile(filename_r)) && !overwrite_existing && return (filename_f, filename_r)

    alns = Alignments(bam_file; include_secondary_alignments=include_secondary_alignments, is_reverse_complement=is_reverse_complement)
    coverage = Coverage(alns; norm=norm, include_chimeric=include_chimeric, max_temp_length=max_temp_length)

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
                                                    suffix_forward="_forward", suffix_reverse="_reverse")
    files.type == ".bam" || throw(AssertionError("Only .bam files accepted for alignments."))
    bw_files = [compute_coverage(file; norm=norm, include_secondary_alignments=include_secondary_alignments, overwrite_existing=overwrite_existing,
                                        is_reverse_complement=is_reverse_complement, max_temp_length=max_temp_length,
                                        include_chimeric=include_chimeric, suffix_forward=suffix_forward, suffix_reverse=suffix_reverse)
                                        for file in files]
    return PairedSingleTypeFiles(bw_files, ".bw", "_forward", "_reverse")
end

function Coverage(alignments::Alignments; norm=1000000, max_temp_length=1000, include_chimeric=false)
    vals_f = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in alignments.chroms)
    vals_r = CoverageValues(chr=>zeros(Float64, len) for (chr, len) in alignments.chroms)
    for alignment in alignments
        if ischimeric(alignment; check_annotation = false, min_distance = max_temp_length)
            include_chimeric || continue
            for part in alignment
                ref = refname(part)
                left = leftposition(part)
                right = rightposition(part)
                ispositivestrand(part) ? (vals_f[ref][left:right] .+= 1.0) : (vals_r[ref][left:right] .-= 1.0)
            end
        else
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

function maxdiffpositions(coverage::Vector{Float64}; max_ppk=5, compute_max_within=5, circular=true)
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
function maxdiffpositions(coverage::Coverage; max_ppk=5, compute_max_within=3, circular=true, invert=false)
    coverage_values = values(coverage; invert=invert)
    Dict{String, Tuple{Vector{Int}, Vector{Int}}}(
        chr=>(
            maxdiffpositions(coverage_values[chr][1]; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular),
            maxdiffpositions(coverage_values[chr][2]; max_ppk=max_ppk, compute_max_within=compute_max_within, circular=circular)
        )
        for chr in sort(collect(keys(coverage_values)))
    )
end