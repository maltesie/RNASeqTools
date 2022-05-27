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
        is_reverse = (isread2(record) != is_reverse_complement)
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
        f = base_coverage.fcount[chr][to][findex] ./ sum(values(base_coverage.fcount[chr]))[findex]
        r = base_coverage.rcount[chr][to][rindex] ./ sum(values(base_coverage.rcount[chr]))[rindex]
        append!(fstats[chr], filter!(x->!isnan(x), f))
        append!(rstats[chr], filter!(x->!isnan(x), r))
    end
    return fstats, rstats
end

function mismatchpositions(base_coverage::BaseCoverage, from::Symbol, to::Symbol; ratio_cut=0.9)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!")) 
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    d = Dict(:A=>DNA_A, :T=>DNA_T, :G=>DNA_G, :C=>DNA_C, :Gap=>DNA_Gap, :N=>DNA_N)
    base_from = d[from]
    pos = Interval{Annotation}[]
    for chr in keys(base_coverage.genome.chroms)
        index = base_coverage.genome[chr] .=== base_from
        ratios = base_coverage.fcount[chr][base_to] ./ sum(values(base_coverage.fcount[chr]))
        ratio_index = (ratios .>= ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .> 0))
            mc = base_coverage.fcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_POS, (mc, round(r; digits=3))))
        end
        index = BioSequences.complement(base_coverage.genome[chr]) .=== base_from
        ratios = base_coverage.rcount[chr][to] ./ sum(values(base_coverage.rcount[chr]))
        ratio_index = (ratios .>= ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .> 0))
            mc = base_coverage.rcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_NEG, (mc, round(r; digits=3))))
        end
    end
    return pos
end

function mismatchpositions(base_coverage::BaseCoverage; check=(:A, :T, :G, :C), ratio_cut=0.9)
    pos = Interval{Annotation}[]
    for from in check, to in check
        base_from === base_to && continue
        mpos = mismatchpositions(base_coverage, from, to; ratio_cut=ratio_cut)
        append!(pos, [Interval(rename(i), leftposition(i):leftpostion(i), strand(i), 
                        Annotation("SNP", "$from->$to", Dict("mutation"=>"$from->$to", 
                                                            "count"=>i.metatdata[1], 
                                                            "frequency"=>i.metadata[2]))) for i in mpos])
    end
    return Features(pos)
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
                                overwrite_existing=false, suffix_forward="_forward", suffix_reverse="_reverse")

    filename_f = bam_file[1:end-4] * suffix_forward * ".bw"
    filename_r = bam_file[1:end-4] * suffix_reverse * ".bw"
    !overwrite_existing && (isfile(filename_f) || isfile(filename_r)) && (return nothing)
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

function rolling_sum(a, n::Int)
    (1<=n<=length(a)) || throw(Assertion("n has to be betwen 1 and $(length(a))."))
    out = similar(a, length(a)-n+1)
    out[1] = sum(a[1:n])
    for i in eachindex(out)[2:end]
        out[i] = out[i-1]-a[i-1]+a[i+n-1]
    end
    return out
end

function Base.diff(coverage::Vector{Float64}, invert::Bool, window_size::Int; min_background_ratio=nothing, filter_local_max=true)
    d = zeros(Float64,length(coverage))
    isnothing(min_background_ratio) || (min_background_increase = min_background_ratio - 1)
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
    if !isnothing(min_background_ratio)
        for i in findall(!iszero, filtered_d)
            check_range = invert ? (i:i+half_window_size) : (i-half_window_size+1:i)
            lim = filtered_d[i]/min_background_increase
            any(@view(coverage[check_range]) .< lim) || (filtered_d[i] = 0.0)
        end
    end
    return filtered_d
end
#using CairoMakie
function tsss(notex::Coverage, tex::Coverage; min_tex_ratio=1.3, min_step=10, window_size=10, min_background_ratio=1.2)
    (notex.chroms == tex.chroms) || throw(Assertion("tex and notex coverage are defined on different reference sequences."))
    chrs = [chr[1] for chr in notex.chroms]
    vals_notex = values(notex)
    vals_tex = values(tex)
    intervals = Vector{Interval{Float64}}()
    for chr in chrs
        notex_f, notex_r = vals_notex[chr]
        tex_f, tex_r = vals_tex[chr]
        d_forward = diff(tex_f, false, window_size; min_background_ratio=min_background_ratio)
        d_reverse = diff(tex_r, true, window_size; min_background_ratio=min_background_ratio)
        check_d_forward = diff(notex_f, false, window_size; filter_local_max=false)
        check_d_reverse = diff(notex_r, true, window_size; filter_local_max=false)
        af = d_forward ./ check_d_forward
        ar = d_reverse ./ check_d_reverse
        fig = Figure(resolution=(1200,400))
		axes = [Axis(fig[1,i]) for i in 1:2]
		hist!(axes[1], abs.(af[abs.(af) .< 2]); bins=50)
		hist!(axes[2], abs.(ar[abs.(ar) .< 2]); bins=50)
        println("infs forward: ", sum(af .=== Inf), "\ninfs reverse: ", sum(ar .=== Inf))
        save("/home/abc/Workspace/TssTerms/$chr.pdf", fig)
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
