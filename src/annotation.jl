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
    s = "Feature (Interval{Annotation}):\n\n[$(leftposition(feature)), $(rightposition(feature))]"
    s *= "on $(refname(feature)) ($(strand(feature))) with annotation $(type(feature)):$(name(feature))\n"
    s *= "parameters: $(paramstring(feature))"
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
        names[name] > 1 ? (n = name[2] * "$(names[name])") : (n = name[2])
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

function Features(coverage::Coverage; type="COV")
    chroms = Dict(seq=>val for (seq,val) in coverage.chroms)
    return Features([Interval(refname(i), leftposition(i), rightposition(i), strand(i), Annotation(type, "", Dict{String,String}("Coverage"=>"$(value(i))"))) for i in coverage], chroms)
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

type(features::Features) = Set(type(f) for f in features)
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

Base.convert(::Type{Interval{Float64}}, i::Interval{T}) where T<:AnnotationStyle = Interval(refname(i), leftposition(i), rightposition(i), strand(i), 0.0)
strand_filter(a::Interval, b::Interval)::Bool = strand(a) === strand(b)
function GenomicFeatures.eachoverlap(features::I, feature::Interval{T}) where {I<:AnnotationContainer,T}
    t = T
    if I === Coverage && T === Annotation
        feature = Interval{Float64}(feature)
        t = Float64
    end
    haskey(features.list.trees, refname(feature)) ?
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), t}(strand_filter,
                GenomicFeatures.ICTreeIntersection{t}(), features.list.trees[refname(feature)], feature)) :
    (return GenomicFeatures.ICTreeIntervalIntersectionIterator{typeof(strand_filter), t}(strand_filter,
                GenomicFeatures.ICTreeIntersection{t}(), GenomicFeatures.ICTree{t}(), feature))
end

function hasoverlap(features::Features, feature::Interval)
    for _ in eachoverlap(features, feature)
        return true
    end
    return false
end
function firstoverlap(features::Features, feature::Interval)
    for int in eachoverlap(features, feature)
        return int
    end
    return nothing
end
function lastoverlap(features::Features, feature::Interval)
    olpinterval = copy(feature)
    for int in eachoverlap(features, feature)
        olpinterval = int
    end
    return olpinterval
end

function paramstring(params::Dict{String,String};priority=("Name",))
    ps = join(("$key=$(params[key])" for key in priority if key in keys(params)), ";")
    os = join(("$key=$(params[key])" for key in sort(collect(keys(params))) if !(key in priority)), ";")
    return ps * ((isempty(ps) || isempty(os)) ? "" : ";") * os
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

nooverlapfromright(features::Features, f::Interval{Annotation}) = hasoverlap(features, f) ? min(rightposition(f), maximum(rightposition(i)+1 for i in eachoverlap(features, f))) : leftposition(f)
nooverlapfromleft(features::Features, f::Interval{Annotation}) = hasoverlap(features, f) ? max(leftposition(f), minimum(leftposition(i)-1 for i in eachoverlap(features, f))) : rightposition(f)

function findutr(f::Interval{Annotation}, features::Features{Annotation}, signals::Features{Annotation}, direction::Symbol;
                        utr_type="UTR", min_utr_length=25, max_utr_length=250,
                        guess_missing=true, expression_key=nothing, source_key=nothing, extra_key=nothing)

    direction in (:left, :right) || throw(AssertionError("direction has to be either :right or :left"))
    left, right =  direction === :right ? (rightposition(f)+1,rightposition(f)+max_utr_length) : (max(leftposition(f)-max_utr_length, 1), leftposition(f)-1)
    olpinterval =  Interval(refname(f), left, right, strand(f), Annotation())
    olps = collect(eachoverlap(signals, olpinterval))
    if isempty(olps)
        if guess_missing
            i = Interval(refname(f),
                            direction === :left ? nooverlapfromright(features, olpinterval) : left,
                            direction === :right ? nooverlapfromleft(features, olpinterval) : right,
                            strand(f), Annotation(utr_type, ""; Name=name(f), source="guess", expression="N/A"))
            return length(i) < min_utr_length ? nothing : i
        else
            return nothing
        end
    else
        maxint = isnothing(expression_key) ? (direction === :left ? olps[end] : olps[1]) : olps[argmax([featureparam(olp, expression_key, Float64) for olp in olps])]
        expression = isnothing(expression_key) ? "N/A" : featureparam(maxint, expression_key)
        source = isnothing(source_key) ? "N/A" : featureparam(maxint, source_key)
        resint = Interval(refname(f),
                        direction === :left ? leftposition(maxint) : left,
                        direction === :right ? rightposition(maxint) : right,
                        strand(f), Annotation(utr_type, ""; Name=name(f), source=source, expression=expression))
        isnothing(extra_key) || setfeatureparam(resint, extra_key, featureparam(maxint, extra_key))
        #utr_type == "3UTR" && println(resint)
        return resint
    end
end

function addutrs!(features::Features, signals::Features, t::Symbol; cds_type="CDS", utr_type="5UTR", min_utr_length=50, max_utr_length=150,
                                                                    guess_missing=true, expression_key=nothing, source_key=nothing, extra_key=nothing)
    t in (:five, :three) || throw(AssertionError("$t not supported!"))
    new_features = Vector{Interval{Annotation}}()
    for feature in features
        type(feature) != cds_type && continue
        direction = ((t === :five) != (strand(feature) === STRAND_POS)) ? :right : :left
        newutr = findutr(feature, features, signals, direction; utr_type=utr_type, min_utr_length=min_utr_length, max_utr_length=max_utr_length,
                                                                    guess_missing=guess_missing, expression_key=expression_key, source_key=source_key, extra_key=extra_key)
        !isnothing(newutr) && push!(new_features, newutr)
    end
    merge!(features, Features(new_features))
end

add5utrs!(features::Features, tss::Features; cds_type="CDS", utr_type="5UTR", min_utr_length=50, max_utr_length=150,
            guess_missing=true, expression_key="avg_tex_plateau", source_key="source") =
addutrs!(features, tss, :five; cds_type=cds_type, utr_type=utr_type, min_utr_length=min_utr_length, max_utr_length=max_utr_length,
                                guess_missing=guess_missing, expression_key=expression_key, source_key=source_key)

add3utrs!(features::Features, terms::Features; cds_type="CDS", utr_type="3UTR", min_utr_length=50, max_utr_length=150,
            guess_missing=true, expression_key="avg_nobcm_plateau", source_key="source", rho_dependent_key="rho_dependent") =
addutrs!(features, terms, :three; cds_type=cds_type, utr_type=utr_type, min_utr_length=min_utr_length, max_utr_length=max_utr_length,
                                    guess_missing=guess_missing, expression_key=expression_key, source_key=source_key, extra_key=rho_dependent_key)


function addigrs!(features::Features; igr_type="IGR", min_igr_length=20)
    new_features = Vector{Interval{Annotation}}()
    base_features_pos = [feature for feature in features if (feature.strand == STRAND_POS)]
    base_features_neg = [feature for feature in features if (feature.strand == STRAND_NEG)]
    for base_features in (base_features_pos, base_features_neg)
        nb_features = base_features === base_features_pos ? length(base_features_pos) : length(base_features_neg)
        for i in 1:nb_features-1
            feature, next_feature = base_features[i], base_features[i+1]
            refname(feature) == refname(next_feature) || continue
            stop, start = leftposition(next_feature), rightposition(feature)
            (stop-1) - (start + 1) > min_igr_length || continue
            igr = Interval(refname(feature), start+1, stop-1, base_features === base_features_pos ? STRAND_POS : STRAND_NEG,
                                        Annotation(igr_type, name(feature)*":"*name(next_feature),
                                        Dict(key=>param(feature, key)*":"*param(next_feature, key) for key in keys(featureparams(feature)) if key in keys(featureparams(next_feature)))))
            push!(new_features, igr)
        end
    end

    for feature in new_features
        push!(features, feature)
    end
end

typenamekey(feature::Interval{Annotation}) = type(feature) * ":" * name(feature)
function featureseqs(features::Features, genome::Genome; key_gen=typenamekey)
    seqs = Vector{LongDNASeq}()
    names = Vector{String}()
    for feature in features
        seq = genome[refname(feature)][leftposition(feature):rightposition(feature)]
        strand(feature) == STRAND_NEG && reverse_complement!(seq)
        push!(seqs, seq)
        push!(names, key_gen(feature))
    end
    return Sequences(seqs, names)
end

function covratio(features::Features, coverage::Coverage)
    vals = values(coverage)
    total = 0.0
    for val in values(vals)
        total += sum(first(val)) + sum(last(val))
    end
    s = 0.0
    for feature in features
        picker = strand(feature) === STRAND_NEG ? last : first
        rightposition(feature) > length(picker(vals[refname(feature)])) && continue
        s += sum(picker(vals[refname(feature)])[leftposition(feature):rightposition(feature)])
    end
    return s/total
end

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

Base.write(fname::String, df::DataFrame) = CSV.write(fname, df)

function Base.show(features::Features)
    chrs = refnames(features)
    ts = type(features)
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