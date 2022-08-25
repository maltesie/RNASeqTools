type(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.type
name(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata.name
refname(feature::Interval{T}) where T<:AnnotationStyle = feature.seqname
featureparams(feature::Interval{Annotation}) = feature.metadata.params
featureparam(feature::Interval{Annotation}, key::String) = feature.metadata.params[key]
featureparam(feature::Interval{Annotation}, key::String, ::Type{T}) where {T} = parse(T, feature.metadata.params[key])
setfeatureparam(feature::Interval{Annotation}, key::String, value::String) = feature.metadata.params[key] = value
hasannotationkey(feature::Interval{Annotation}, key::String) = key in keys(featureparams(feature))
annotation(feature::Interval{T}) where T<:AnnotationStyle = feature.metadata
hasannotation(feature::Interval{T}) where T<:AnnotationStyle = !isempty(annotation(feature))
ispositivestrand(feature::Interval{T}) where T<:AnnotationStyle = strand(feature) === STRAND_POS

function summarize(feature::Interval{Annotation})
    s = "Feature: [$(leftposition(feature)), $(rightposition(feature))] on $(refname(feature)) ($(strand(feature))) - "
    s *= hasannotation(feature) ? "annotation: $(type(feature)):$(name(feature))" : "not annotated"
    s *= "\nparameters: $(paramstring(feature))\n"
end
function Base.show(feature::Interval{Annotation})
    println(summarize(feature))
end

function Annotation()
    Annotation("", "", Dict{String, String}())
end

function Annotation(type::String, name::String; args...)
    Annotation(type, name, Dict{String, String}("$key"=>value for (key,value) in args))
end

Base.isempty(annotation::Annotation) = isempty(annotation.type) && isempty(annotation.name)

function AlignmentAnnotation()
    AlignmentAnnotation("", "", 0)
end

Base.isempty(annotation::AlignmentAnnotation) = isempty(annotation.type) && isempty(annotation.name) && (annotation.overlap==0)

name(annot::T) where {T<:AnnotationStyle} = annot.name
type(annot::T) where {T<:AnnotationStyle} = annot.type
overlap(annot::AlignmentAnnotation) = annot.overlap

function BaseAnnotation(feature::Interval{Annotation}, base_coverage::BaseCoverage)
    left = leftposition(feature)
    right = rightposition(feature)
    seq = base_coverage.ref_seq[left:right]
    ispositivestrand(feature) || reverse_complement!(seq)
    count = ispositivestrand(feature) ? errorcov.fcount : errorcov.rcount
    r = ispositivestrand(feature) ? (left:right) : (right:-1:left)
    ref = Int[(seq[i] in (DNA_A, DNA_T, DNA_G, DNA_C)) ? count[seq[i]][ii] : 0 for (i, ii) in enumerate(r)]
    BaseAnnotation(type(feature), name(feature), ref, count[ref][:A][r], count[:T][r], count[:G][r], count[:C][r], count[:Gap][r], count[:Ins][r])
end

Features() = Features(Annotation)

function Features(::Type{T}) where T
    return Features(Vector{Interval{T}}())
end

function Features(feature_list::Vector{Interval{Annotation}})
    return Features(IntervalCollection(feature_list, true), Dict{String, Int}())
end

function Features(feature_list::Vector{Interval{Annotation}}, chroms::Dict{String, Int})
    return Features(IntervalCollection(feature_list, true), chroms)
end

function Features(gff_file::String, type::Vector{String}; name_keys=["Name"], same_name_rule=:all)
    same_name_rule in (:first, :all, :none) || throw(AssertionError("same_name_rule must be :first, :all, or :none"))
    features = GFF3.Reader(open(gff_file), skip_directives=false)
    intervals = Vector{Interval{Annotation}}()
    names = Dict{Tuple{String,String},Int}()
    chroms = Dict{String, Int}()
    for feature in features
        if GFF3.isdirective(feature)
            line = GFF3.content(feature)
            if startswith(line, "sequence-region")
                seqid, _, rp = split(line)[2:4]
                push!(chroms, seqid => parse(Int, rp))
            end
            continue
        end
        GFF3.featuretype(feature) == "Source" && continue
        (GFF3.featuretype(feature) in type || isempty(type)) || continue
        seqn = GFF3.seqid(feature)
        name = ("NA", "NA")
        as = GFF3.attributes(feature)
        for name_key in reverse(name_keys)
            for (key, values) in as
                name_key === key || continue
                name = (GFF3.featuretype(feature), join(values, ","))
            end
        end
        name in keys(names) ? (names[name]+=1) : (names[name]=1)
        (same_name_rule === :first && names[name] > 1) && continue
        n = names[name] > 1 ? name[2] * "$(names[name])" : name[2]
        annot = Annotation(GFF3.featuretype(feature), n, Dict(pair[1] => join(pair[2], ",") for pair in GFF3.attributes(feature)))
        push!(intervals, Interval(seqn, GFF3.seqstart(feature), GFF3.seqend(feature), GFF3.strand(feature), annot))
    end
    same_name_rule === :none && (return Features([i for i in intervals if ((type(i), name(i)) in names && names[(type(i), name(i))] == 1)]))
    return Features(intervals, chroms)
end

function Features(gff_file::String, type::String; name_keys=["Name"], same_name_rule=:all)
    return Features(gff_file, [type], name_keys=name_keys, same_name_rule=same_name_rule)
end

function Features(gff_file::String; name_keys=["Name"], same_name_rule=:all)
    return Features(gff_file, String[], name_keys=name_keys, same_name_rule=same_name_rule)
end

function Features(gff_file::String, bam_file::String, genome::Genome)
    new_intervals = Interval{BaseAnnotation}[]
    features = Features(gff_file)
    errorcov = BaseCoverage(bam_file, genome)
    for feature in features
        annotation = BaseAnnotation(feature, errorcov)
        push!(new_intervals, Interval(refname(feature), leftposition(feature), rightposition(feature), strand(feature), annotation))
    end
    return Features(IntervalCollection(new_intervals, true), genome.chroms)
end

types(features::Features) = Set(type(f) for f in features)
refnames(features::Features) = collect(keys(features.list.trees))
function overlaps(alignmentinterval::Interval{T}, feature::Interval{Annotation}) where T<:AnnotationStyle
    refname(feature) == refname(alignmentinterval) || return false
    strand(feature) == strand(alignmentinterval) || return false
    return leftposition(feature) <= rightposition(alignmentinterval) && leftposition(alignmentinterval) <= rightposition(feature)
end

Base.push!(features::Features, interval::Interval) = push!(features.list, interval)
function Base.merge!(features1::Features{T}, features2::Features{T}) where T
    for feature in features2
        push!(features1, feature)
    end
end
function Base.merge(features1::Features{T}, features2::Features{T}) where T
    re = Features(T)
    for feature in features1
        push!(re, feature)
    end
    for feature in features2
        push!(re, feature)
    end
    return re
end
Base.:*(features1::Features, features2::Features) = merge(features1, features2)

Base.iterate(features::T) where T<:AnnotationContainer = iterate(features.list)
Base.iterate(features::T, state::Tuple{Int64,GenomicFeatures.ICTree{I},GenomicFeatures.ICTreeIteratorState{I}}) where {T<:AnnotationContainer,I<:AnnotationStyle} = iterate(features.list, state)
Base.length(features::T) where T<:AnnotationContainer = length(features.list)
Base.split(features::T) where T<:AnnotationContainer = [Features([feature for feature in features if type(feature)==t], [t]) for t in types(features)]

function Base.getindex(features::T, i::Int) where T<:AnnotationContainer
    (i < 1 || i > length(features)) && throw(AssertionError("Trying to access $(length(features))-element $T at index $i"))
    for (c, f) in enumerate(features)
        c === i && return f
    end
end

Base.convert(::Type{Interval{Annotation}}, i::Interval{Nothing}) = Interval(refname(i), leftposition(i), rightposition(i), strand(i), Annotation())
strand_filter(a::Interval, b::Interval)::Bool = strand(a) === strand(b)
function GenomicFeatures.eachoverlap(feature::Interval{T}, features::Features{Annotation}) where T
    T !== Annotation && (feature = Interval{Annotation}(feature))
    haskey(features.list.trees, refname(feature)) ?
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter,
                GenomicFeatures.ICTreeIntersection{Annotation}(), features.list.trees[refname(feature)], feature)) :
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), Annotation}(strand_filter,
                GenomicFeatures.ICTreeIntersection{Annotation}(), GenomicFeatures.ICTree{Annotation}(), feature))
end

function hasoverlap(features::Features, feature::Interval)
    for _ in eachoverlap(feature, features)
        return true
    end
    return false
end
function firstoverlap(features::Features, feature::Interval)
    for int in eachoverlap(feature, features)
        return int
    end
    return nothing
end
function lastoverlap(features::Features, feature::Interval)
    olpinterval = copy(feature)
    for int in eachoverlap(feature, features)
        olpinterval = int
    end
    return olpinterval
end

function paramstring(params::Dict{String,String}; priority=("Name",))
    ps = join(("$key=$(params[key])" for key in priority if key in keys(params)), ";")
    os = join(("$key=$(params[key])" for key in sort(collect(keys(params))) if !(key in priority)), ";")
    return ps * ((isempty(ps) || isempty(os)) ? "none" : "; ") * os
end
paramstring(feature::Interval{Annotation}; priority=("Name",)) = paramstring(featureparams(feature); priority=priority)

function Base.write(file::String, features::Features; zip=false, tabix=false)
    chroms = copy(features.chroms)
    writer = GFF3.Writer(open(file, "w"))
    write(writer, GFF3.Record("##gff-version 3.2.1"))
    for feature in features
        #println(feature)
        if refname(feature) in keys(chroms)
            write(writer, GFF3.Record("##sequence-region $(refname(feature)) 1 $(chroms[refname(feature)])"))
            delete!(chroms, refname(feature))
        end
        write(writer, GFF3.Record("$(refname(feature))\t.\t$(type(feature))\t$(feature.first)\t$(feature.last)\t.\t$(feature.strand)\t.\t$(paramstring(featureparams(feature)))"))
    end
    b = close(writer)
    sleep(0.5)
    zip && run(`bgzip --force $file`)
    tabix && run(`tabix $(file).gz`)
    return b
end

typenamekey(feature::Interval{Annotation}) = type(feature) * ":" * name(feature)
function featureseqs(features::Features, genome::Genome; key_gen=typenamekey)
    seqs = Vector{LongDNA}()
    names = Vector{String}()
    for feature in features
        seq = genome[refname(feature)][leftposition(feature):rightposition(feature)]
        strand(feature) == STRAND_NEG && reverse_complement!(seq)
        push!(seqs, seq)
        push!(names, key_gen(feature))
    end
    return Sequences(seqs, names)
end

function covratio(features::Features, coverage::Coverage; round_digits=2)
    vals = values(coverage)
    total = 0.0
    for (values_forward, values_reverse) in values(vals)
        total += sum(values_forward) - sum(values_reverse)
    end
    ratios = Pair{String, Float64}[]
    for t in types(features)
        s = 0.0
        check_index = Dict(chr=>(zeros(Bool, length(v1)), zeros(Bool, length(v2))) for (chr, (v1,v2)) in vals)
        for feature in features
            type(feature) == t || continue
            is_negative_strand = strand(feature) === STRAND_NEG
            check_index[refname(feature)][is_negative_strand ? 2 : 1][leftposition(feature):rightposition(feature)] .= true
        end
        for (chr, (v1, v2)) in vals
            s += sum(v1[check_index[chr][1]])
            s -= sum(v2[check_index[chr][2]])
        end
        push!(ratios, t=>round(s/total; digits=round_digits))
    end
    return ratios
end
covratio(features::Features, coverages::PairedSingleTypeFiles; round_digits=2) =
    [covratio(features, Coverage(f1, f2); round_digits=round_digits) for (f1, f2) in coverages]

function asdataframe(features::Features; add_keys=:all)
    add_keys === :none && (add_keys = Set())
    add_keys === :all && (add_keys = Set(key for feature in features for key in keys(featureparams(feature))))
    df = DataFrame(name=repeat([""], length(features)), refname=repeat([""], length(features)), type=repeat([""], length(features)),
                    left=repeat([-1], length(features)), right=repeat([-1], length(features)), strand=repeat(['*'], length(features)))
    for key in add_keys
        df[!, Symbol(key)] = repeat([""], length(features))
    end
    for (i,feature) in enumerate(features)
        df[i, :name] = name(feature)
        df[i, :refname] = refname(feature)
        df[i, :type] = type(feature)
        df[i, :left] = leftposition(feature)
        df[i, :right] = rightposition(feature)
        df[i, :strand] = strand(feature) === STRAND_NEG ? '-' : '+'
        pa = featureparams(feature)
        for key in sort(collect(keys(pa)))
            (key in add_keys || add_keys === :all) && (df[i, Symbol(key)] = pa[key])
        end
    end
    return df
end

function Base.show(features::Features)
    chrs = refnames(features)
    ts = types(features)
    stats = Dict(chr=>Dict(t=>0 for t in ts) for chr in chrs)
    nb_pos = 0
    nb_neg = 0
    for f in features
        stats[refname(f)][type(f)] += 1
        nb_pos += strand(f) === STRAND_POS
        nb_neg += strand(f) === STRAND_NEG
    end
    dt = Int(floor(maximum(length(t) for t in ts)/8))+1
    printstring = "\n$(nb_pos+nb_neg) features in total with $nb_pos on + strand and $nb_neg on - strand:\n\n"
    printstring *= "type$(repeat("\t",dt))$(join([chr[1:min(8,length(chr))]*repeat(" ", 9-min(8,length(chr))) for chr in chrs], "\t"))\n\n"
    for t in ts
        printstring *= "$(t)$(repeat("\t",dt-Int(floor(length(t)/8))))$(join([stats[chr][t] for chr in chrs], "\t\t"))\n"
    end
    println(printstring)
end

function embl_to_gff(embl_file::String, gff_file::String, chrs::Vector{String})
    chroms = copy(reverse(chrs))
    writer = GFF3.Writer(open(gff_file, "w"))
    write(writer, GFF3.Record("##gff-version 3.2.1"))
    lines = readlines(embl_file)
    current_refname = pop!(chroms)
    for i in 1:length(lines)-1
        line = lines[i]
        next_line = lines[i+1]
        if startswith(line, "FT")

            if occursin(line, "source")
                _, _, r = split(line)
                _, l = split(r, "..")
                write(writer, GFF3.Record("##sequence-region $current_refname 1 $l"))
                write(writer, GFF3.Record("$current_refname\t.\tsource\t1\t$l\t.\t.\t.\tName=$current_refname"))
            end
            if occursin(line, "gene")
                _, _, r = split(line)
                strand = "+"
                startswith(r, "complement") && (r = r[12:end-1]; strand="-")
                l, r = split(r, "..")
                name = split(next_line, "locus_tag=")[end]
                write(writer, GFF3.Record("$current_refname\t.\tCDS\t$l\t$r\t.\t$strand\t.\tName=$name"))
            end
        end

    end
    close(writer)
end