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

struct BaseCoverage
    genome::Genome
    fcount::Dict{String, Dict{Symbol, Vector{Int}}}
    rcount::Dict{String, Dict{Symbol, Vector{Int}}}
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

struct BaseAnnotation <: AnnotationStyle
    type::String
    name::String
    ref::Vector{Int}
    a::Vector{Int}
    t::Vector{Int}
    g::Vector{Int}
    c::Vector{Int}
    gap::Vector{Int}
    ins::Vector{Int}
end

function coverage(feature::Interval{BaseAnnotation})
    return feature.metadata.a .+ feature.metadata.t .+ feature.metadata.g .+ feature.metadata.c .+ feature.metadata.gap
end

refcount(feature::Interval{BaseAnnotation}) = feature.metadata.ref

function compute_coverage(bam_file::String; norm=1000000, include_secondary_alignments=false, is_reverse_complement=false, max_temp_length=500,
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
    for alignment in Alignments(bam_file; include_secondary_alignments=include_secondary_alignments, is_reverse_complement=is_reverse_complement)
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

function compute_coverage(files::SingleTypeFiles; norm=1000000, include_secondary_alignments=true, overwrite_existing=false, is_reverse_complement=false)
    files.type == ".bam" || throw(AssertionError("Only .bam files accepted for alignments."))
    bw_files = Vector{Tuple{String, String}}()
    for file in files
        filename_f = file[1:end-4] * "_forward.bw"
        filename_r = file[1:end-4] * "_reverse.bw"
        push!(bw_files, (filename_f, filename_r))
        (!overwrite_existing && isfile(filename_f) && isfile(filename_r)) && continue
        compute_coverage(file; norm=norm, include_secondary_alignments=include_secondary_alignments, is_reverse_complement=is_reverse_complement)
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

all_fisher_test_combinations(notex::Vector{Tuple{Int,Int}}, tex::Vector{Tuple{Int,Int}}) = vec([pvalue(FisherExactTest(nt...,t...); tail=:right) for t in tex, nt in notex])
function all_fisher_test_combinations(notex::Matrix{Int}, tex::Matrix{Int})
    tex_rows, tex_samples = size(tex)
    notex_rows, notex_samples = size(notex)
    tex_rows == notex_rows || throw(AssertionError("tex and notex have to have the same number of rows!"))
    tex_samples % 2 == 0 && notex_samples % 2 == 0 || throw(AssertionError("tex and notex have to have even amount of columns!"))
    notex_pairs = [[(i,j) for (i,j) in partition(r, 2)] for r in eachrow(notex)]
    tex_pairs = [[(i,j) for (i,j) in partition(r, 2)] for r in eachrow(tex)]
    all_fisher_test_combinations.(notex_pairs, tex_pairs)
end

fisher_combined_pvalue(pvalues::Vector{Float64}) = 1-cdf(Chisq(2*length(pvalues)), -2*sum(log.(pvalues)))

test_multiple_fisher_combined(notex::Matrix{Int}, tex::Matrix{Int}) = adjust(PValues(fisher_combined_pvalue.(all_fisher_test_combinations(notex, tex))), BenjaminiHochberg())

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

function maxdiffsignalpositions(coverage::Vector{Float64}; max_ppk=5, circular=true)
    all(coverage .>= 0) || all(coverage .<= 0) || throw(AssertionError("Coverage values have to be positive for + strand and negative for - strand!"))
    is_negative_strand = all(coverage .<= 0)
    max_nb_peaks = length(coverage) / 1000 * max_ppk
    d = zeros(Float64,length(coverage))
    d[is_negative_strand ? (1:length(d)-1) : (2:length(d))] = @view(coverage[2:end]) .- @view(coverage[1:end-1])
    circular && (d[is_negative_strand ? length(d) : 1] = coverage[1]-coverage[end])
    peak_index = [i for i in 2:length(d)-1 if d[i-1]<d[i]>d[i+1]]
    if circular
        if d[end]<d[1]>d[2]
            push!(peak_index, 1)
        elseif d[end-1]<d[end]>d[1]
            push!(peak_index, length(d))
        end
    else
        d[1]>d[2] && push!(peak_index, 1)
        d[end]>d[end-1] && push!(peak_index, length(d))
    end
    for cut in 1:Int(floor(maximum(d)))
        length(peak_index) < max_nb_peaks && break
        peak_index = peak_index[d[peak_index] .>= cut]
    end
    return peak_index
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

function tsss(notexs::Vector{Coverage}, texs::Vector{Coverage}; max_fdr=0.05, compute_step_within=1, max_ppk=5, circular=true, source="NA")
    chroms = vcat([s.chroms for s in notexs], [s.chroms for s in texs])
    all(chroms[1] == chrom for chrom in chroms) || throw(Assertion("All coverages have to be defined on the same reference sequences."))
    chr_names = collect(chr[1] for chr in chroms[1])
    merged_tex_values = values(merge(texs))
    peaks = Dict{String, Tuple{Vector{Int}, Vector{Int}}}(
        chr=>(
            maxdiffsignalpositions(merged_tex_values[chr][1]; max_ppk=max_ppk, circular=circular),
            maxdiffsignalpositions(merged_tex_values[chr][2]; max_ppk=max_ppk, circular=circular)
        )
        for chr in chr_names
    )
    #println(peaks["NC_002505"][1][1:10])
    valss = vcat([values(notex) for notex in notexs], [values(tex) for tex in texs])
    background_plateau_matrix =
    hcat(
        [vcat(
            [vcat(
                background_plateau_pairs(vals[chr][1], peaks[chr][1]; compute_step_within=compute_step_within),
                background_plateau_pairs(vals[chr][2], peaks[chr][2]; compute_step_within=compute_step_within)
            )
            for chr in chr_names]...
        )
        for vals in valss]...
    ) .+ 1
    #println(background_plateau_matrix[1:10, :])
    #adjp_fisher = test_multiple_fisher_combined(background_plateau_matrix[:, 1:(2*length(notexs))], background_plateau_matrix[:, (2*length(notexs))+1:end])
    counts_between = Counts(Dict("notex_plateau"=>1:length(notexs), 
                                 "tex_plateau"=>length(notexs)+1:(length(texs) + length(notexs))), 
                            background_plateau_matrix[:,2:2:end])
    tex_background_index = collect(2*length(notexs)+1:2:2*(length(notexs)+length(texs)))
    tex_plateau_index = tex_background_index .+ 1
    counts_within = Counts(Dict("tex_background"=>1:length(texs), 
                                "tex_plateau"=>length(texs)+1:(2*length(texs))), 
                            background_plateau_matrix[:,vcat(tex_background_index,tex_plateau_index)])

    normalize!(counts_between; normalization_method=:rle)
    #println(counts_between.values[1:10, :])
    (basevals, _, adjp_between) = dge_glm(counts_between, "notex_plateau", "tex_plateau"; tail=:right)

    avg_sample::Vector{Float64} = [gmean(counts_within.values[i, length(texs)+1:end]) for i in 1:length(counts_within)]
    logdiffs = [log2.(counts_within.values[:, length(texs)+i]) .- log2.(avg_sample) for i in 1:length(texs)]
    for (i, nf) in enumerate(2 .^ [median(logdiff[(!).(isnan.(logdiff))]) for logdiff in logdiffs])
        counts_within.values[:, [i,i+length(texs)]] ./= nf
    end
    (_, _, adjp_within) = dge_glm(counts_within, "tex_background", "tex_plateau"; tail=:right)
    #println(adjp_between[1:10])
    #println(adjp_within[1:10])
    offset = 0
    intervals = Vector{Interval{Annotation}}()
    for chr in chr_names
        for (cp, s) in zip(peaks[chr], (STRAND_POS, STRAND_NEG))
            append!(intervals, [Interval(chr, p, p, s,
                                    Annotation("TSS", "", 
                                        Dict("fdr_background_ratio"=>"$(round(adjp_within[offset+i], digits=8))",
                                            "fdr_tex_ratio"=>"$(round(adjp_between[offset+i], digits=8))",
                                            "height"=>"$(round(mean(basevals[offset+i]), digits=2))",
                                            "source"=>source))) for (i,p) in enumerate(cp) if adjp_between[offset+i]<=max_fdr && adjp_within[offset+i]<=max_fdr])
            #offset == 0 && println(intervals[1:10])
            offset += length(cp)
        end
    end
    return Features(intervals)
end

function terms(coverage::Coverage; min_step=10, half_window_size=5, min_background_ratio=1.2)
    vals = values(coverage)
    intervals = Vector{Interval{Float64}}()
    chrs = [chr[1] for chr in coverage.chroms]
    for chr in chrs
        f, r = vals[chr]
        d_forward = diff(f, true; half_window_size=half_window_size, min_background_ratio=min_background_ratio)
        d_reverse = diff(r, false; half_window_size=half_window_size, min_background_ratio=min_background_ratio)
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
