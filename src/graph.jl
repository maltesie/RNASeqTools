struct Interactions <: InteractionContainer
    graph::SimpleDiGraph
    nodes::DataFrame
    edges::DataFrame
    edgestats::ElasticArray{Int}
    replicate_ids::Vector{Symbol}
end

"""
Method of write function which saves the Interactions struct in a jld2 file.
"""
function Base.write(filepath::String, interactions::Interactions)
    if !endswith(filepath, ".jld2")
        throw(ArgumentError("Append '.jld2' to filepath"))
    else
        save(filepath, "interactions", interactions)
    end
end

function leftmostposition(alnpart::AlignedPart, alnread::AlignedRead)
    minimum(alnread.alns.leftpos[i] for i::Int in alnread.range if (isassigned(alnread.alns.annames, i) && alnread.alns.annames[i] === name(alnpart)))
end

function rightmostposition(alnpart::AlignedPart, alnread::AlignedRead)
    maximum(alnread.alns.rightpos[i] for i::Int in alnread.range if (isassigned(alnread.alns.annames, i) && alnread.alns.annames[i] === name(alnpart)))
end

function leftmostrelposition(alnpart::AlignedPart, alnread::AlignedRead; five_type="5UTR", three_type="3UTR")
    minimum(alnread.alns.antypes[i] === five_type ? 0x00 : (alnread.alns.antypes[i] === three_type ? 0x65 : alnread.alns.anleftrel[i])
        for i::Int in alnread.range if (isassigned(alnread.alns.annames, i) && alnread.alns.annames[i] === name(alnpart)))
end

function rightmostrelposition(alnpart::AlignedPart, alnread::AlignedRead; five_type="5UTR", three_type="3UTR")
    maximum(alnread.alns.antypes[i] === five_type ? 0x00 : (alnread.alns.antypes[i] === three_type ? 0x65 : alnread.alns.anrightrel[i])
        for i::Int in alnread.range if (isassigned(alnread.alns.annames, i) && alnread.alns.annames[i] === name(alnpart)))
end

function myhash(part::AlignedPart; use_type=false)
    return use_type ? hash(name(part)*type(part)) : hash(name(part))
end

function Base.append!(interactions::Interactions, alignments::Alignments, replicate_id::Symbol;
                        min_distance=1000, filter_types=[], merge_annotation_types=true, merge_type="merge", cds_type="CDS", five_type="5UTR", three_type="3UTR")
    if !(String(replicate_id) in interactions.replicate_ids)
        interactions.edges[:, replicate_id] = repeat([false], nrow(interactions.edges))
        push!(interactions.replicate_ids, replicate_id)
    end
    trans = Dict{UInt, Int}(interactions.nodes[i, :hash]=>i for i in 1:nrow(interactions.nodes))
    trans_edges = Dict{Tuple{Int,Int},Int}((interactions.edges[i, :src],interactions.edges[i, :dst])=>i for i in 1:nrow(interactions.edges))

    for alignment in alignments
        !isempty(filter_types) && typein(alignment, filter_types) && continue
        is_chimeric = ischimeric(alignment; min_distance=min_distance)
        is_multi = is_chimeric ? ismulti(alignment) : false
        alnparts = parts(alignment)

        for (i,part) in enumerate(alnparts)
            hasannotation(part) || continue
            any(samename(part, formerpart) for formerpart in alnparts[1:i-1]) && continue

            h = myhash(part; use_type=!merge_annotation_types)
            if !(h in keys(trans))
                trans[h] = length(trans) + 1
                add_vertex!(interactions.graph)
                push!(interactions.nodes, (name(part), (merge_annotation_types && (type(part) in (cds_type, three_type, five_type))) ? merge_type : type(part), refname(part), 0, 0, strand(part), h))
            end
            is_chimeric || (interactions.nodes[trans[h], :nb_single] += 1)
        end

        for (part1, part2) in combinations(alnparts[collect(!any(samename(part, formerpart)
                                                                        for formerpart in alnparts[1:i-1])
                                                                            for (i, part) in enumerate(alnparts))], 2)
            (hasannotation(part1) && hasannotation(part2)) || continue
            ischimeric(part1, part2; min_distance=min_distance) || continue
            a, b = trans[myhash(part1; use_type=!merge_annotation_types)], trans[myhash(part2; use_type=!merge_annotation_types)]
            interactions.nodes[a, :nb_ints] += 1
            interactions.nodes[b, :nb_ints] += 1
            if !has_edge(interactions.graph, a, b)
                add_edge!(interactions.graph, a, b)
                trans_edges[(a,b)] = ne(interactions.graph)
            end
            iindex = trans_edges[(a, b)]
            left1, right1 = leftmostposition(part1, alignment), rightmostposition(part1, alignment)
            left2, right2 = leftmostposition(part2, alignment), rightmostposition(part2, alignment)
            statindex1, statindex2 = rightmostrelposition(part1, alignment) + 0x01, leftmostrelposition(part2, alignment) + 0x67
            nms1, nms2 = missmatchcount(part1), missmatchcount(part2)
            if iindex > nrow(interactions.edges)
                push!(interactions.edges, (a, b, 0, 0, left1, right1, left2, right2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                            (nms1>1 ? 1.0 : 0.0), (nms2>1 ? 1.0 : 0.0), (false for i in 1:length(interactions.replicate_ids))...))
                append!(interactions.edgestats, zeros(Int, 204))
            end
            interactions.edges[iindex, :nb_ints] += 1
            interactions.edgestats[statindex1, iindex] += 1
            interactions.edgestats[statindex2, iindex] += 1
            is_multi && (interactions.edges[iindex, :nb_multi] += 1)
            interactions.edges[iindex, :minleft1] = min(interactions.edges[iindex, :minleft1], left1)
            interactions.edges[iindex, :maxright1] = max(interactions.edges[iindex, :maxright1], right1)
            interactions.edges[iindex, :minleft2] = min(interactions.edges[iindex, :minleft2], left2)
            interactions.edges[iindex, :maxright2] = max(interactions.edges[iindex, :maxright2], right2)
            for (s,v) in zip((:meanleft1, :meanright1, :meanleft2, :meanright2, :meanlength1, :meanlength2, :meanmiss1, :meanmiss2, :nms1, :nms2),
                             (left1, right1, left2, right2, right1 - left1 + 1, right2 - left2 + 1, nms1, nms2, Int(nms1>0), Int(nms2>0)))
                interactions.edges[iindex, s] += (v - interactions.edges[iindex, s]) / interactions.edges[iindex, :nb_ints]
            end
            interactions.edges[iindex, replicate_id] = true
        end
    end
    return interactions
end

function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_ints=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :minleft1=>Int[], :maxright1=>Int[], :minleft2=>Int[],
                        :maxright2=>Int[], :meanlength1=>Float64[], :meanlength2=>Float64[], :meanleft1=>Float64[], :meanleft2=>Float64[],
                        :meanright1=>Float64[], :meanright2=>Float64[], :nms1=>Float64[], :nms2=>Float64[], :meanmiss1=>Float64[], :meanmiss2=>Float64[])
    return Interactions(SimpleDiGraph(), nodes, edges, ElasticArray{Int}(undef, 204, 0), Symbol[])
end

function Interactions(alignments::Alignments; replicate_id=:first, min_distance=1000, filter_types=[])
    append!(Interactions(), alignments, replicate_id, min_distance=min_distance, filter_types=filter_types)
end

"""
Load Interactions struct from jld2 file.
"""
Interactions(filepath::String) = load(filepath, "interactions")

function checkinteractions(ints::Interactions, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05, check_uppercase=true)
	verified_dict = merge(Dict(pair=>zeros(2) for pair in verified_pairs),
							Dict(reverse(pair)=>zeros(2) for pair in verified_pairs))
    check_uppercase && (verified_dict = Dict((uppercase(p[1]), uppercase(p[2]))=>v for (p,v) in verified_dict))
	df = asdataframe(ints; min_reads=min_reads, max_fdr=max_fdr)
    for row in eachrow(df)
        key = (row[:name1], row[:name2])
        check_uppercase && (key = (uppercase(key[1]), uppercase(key[2])))
        if key in keys(verified_dict)
            verified_dict[key][1] = row[:in_libs]
            verified_dict[key][2] = row[:nb_ints]
        end
    end

	sorted_keys = vcat([[pair, reverse(pair)] for pair in verified_pairs]...)
    check_uppercase && (sorted_keys = [(uppercase(p[1]), uppercase(p[2])) for p in sorted_keys])
	m = reduce(hcat, [verified_dict[key] for key in sorted_keys])'
	verified_stats = DataFrame(
		name1=String[n[1] for n in verified_pairs],
		name2=String[n[2] for n in verified_pairs],
		libs=max.(m[:,1][1:2:end], m[:,1][2:2:end]), count=m[:,2][1:2:end] .+ m[:,2][2:2:end]
		)
	return verified_stats
end
function checkinteractions(conditions::Vector{Interactions}, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05)
    verified_stats = DataFrame(name1=String[p[1] for p in verified_pairs], name2=String[p[2] for p in verified_pairs])
    for ints in conditions
        verified_stats = innerjoin(verified_stats, checkinteractions(ints, verified_pairs; min_reads=min_reads, max_fdr=max_fdr); on=[:name1, :name2], makeunique=true)
    end
    return verified_stats
end
function checkinteractions(graph_files::SingleTypeFiles, verified_pairs::Vector{Tuple{String,String}}; min_reads=5, max_fdr=0.05)
    graph_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(graph_file) for graph_file in graph_files]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr=max_fdr)
end
function checkinteractions(graph_files::SingleTypeFiles, verified_pairs_file::String; min_reads=5, max_fdr=0.05)
    graph_files.type === ".jld2" || throw(AssertionError("Can only read .jld2 files!"))
    conds = [Interactions(graph_file) for graph_file in graph_files]
    df = DataFrame(CSV.File(verified_pairs_file; stringtype=String))
    verified_pairs = [(a,b) for (a,b) in eachrow(df)]
    return checkinteractions(conds, verified_pairs; min_reads=min_reads, max_fdr=max_fdr)
end

function uniqueinteractions(ints::Interactions; min_reads=5, max_fdr=0.05)
    df = asdataframe(ints; min_reads=min_reads, max_fdr=max_fdr)
    Set(Set((a,b)) for (a,b) in zip(df[!, :name1], df[!, :name2]))
end

function addpvalues!(interactions::Interactions; method=:disparity)
    @assert method in (:disparity, :fisher)
    pvalues = ones(ne(interactions.graph))
    all_interactions = sum(interactions.edges[!, :nb_ints])+1

    if method === :fisher
        ints_between = interactions.edges[!,:nb_ints]
        ints_other_source = interactions.nodes[interactions.edges[!, :src], :nb_ints] .- ints_between
        ints_other_target = interactions.nodes[interactions.edges[!, :dst], :nb_ints] .- ints_between
        ints_other = all_interactions .- ints_between .- ints_other_source .- ints_other_target
        tests = FisherExactTest.(ints_between, ints_other_target, ints_other_source, ints_other)
        pvalues = pvalue.(tests; tail=:right)
    elseif method === :disparity
        degrees = [degree(interactions.graph, i) - 1 for i in 1:nv(interactions.graph)]
        p_source = (1 .- interactions.edges[!,:nb_ints] ./ interactions.nodes[interactions.edges[!, :src], :nb_ints]).^degrees[interactions.edges[!, :src]]
        p_target = (1 .- interactions.edges[!,:nb_ints] ./ interactions.nodes[interactions.edges[!, :dst], :nb_ints]).^degrees[interactions.edges[!, :dst]]
        pvalues = min.(p_source, p_target)
    end

    adjp = adjust(PValues(pvalues), BenjaminiHochberg())
    interactions.edges[:, :p_value] = pvalues
    interactions.edges[:, :fdr] = adjp
    for colname in (:relmean1, :relmean2, :relmin1, :relmin2, :relmax1, :relmax2)
        interactions.edges[:, colname] = ones(Float64, nrow(interactions.edges)) * -Inf
    end
    return interactions
end

function addrelativepositions!(interactions::Interactions, features::Features; has_merge_type=true, merge_type="merge", cds_type="CDS")
    tus = Dict(name(feature)*type(feature)=>(leftposition(feature), rightposition(feature)) for feature in features)
    for edge_row in eachrow(interactions.edges)
        name_src = interactions.nodes[edge_row[:src], :name]
        type_src = has_merge_type && (interactions.nodes[edge_row[:src], :type] === merge_type) ? cds_type : interactions.nodes[edge_row[:src], :type]
        (feature1_left, feature1_right) = tus[name_src*type_src]
        name_dst = interactions.nodes[edge_row[:dst], :name]
        type_dst = has_merge_type && (interactions.nodes[edge_row[:dst], :type] === merge_type) ? cds_type : interactions.nodes[edge_row[:dst], :type]
        (feature2_left, feature2_right) = tus[name_dst*type_dst]
        isnegative1 = interactions.nodes[edge_row[:src], :strand] === '-'
        isnegative2 = interactions.nodes[edge_row[:dst], :strand] === '-'
        p1 = collect(edge_row[[(isnegative1 ? :meanleft1 : :meanright1), :minleft1, :maxright1]])
        p2 = collect(edge_row[[(isnegative2 ? :meanright2 : :meanleft2), :minleft2, :maxright2]])
        (relpos1, relmin1, relmax1) = min.(1.0, max.(0.0, (p1 .- feature1_left) ./ (feature1_right - feature1_left)))
        (relpos2, relmin2, relmax2) = min.(1.0, max.(0.0, (p2 .- feature2_left) ./ (feature2_right - feature2_left)))
        isnegative1 && ((relpos1, relmin1, relmax1) = (1-relpos1, 1-relmax1, 1-relmin1))
        isnegative2 && ((relpos2, relmin2, relmax2) = (1-relpos2, 1-relmax2, 1-relmin2))
        edge_row[[:relmean1, :relmean2, :relmin1, :relmin2, :relmax1, :relmax2]] =
            round.((relpos1, relpos2, relmin1, relmin2, relmax1, relmax2); digits=4)
    end
    return interactions
end

function asdataframe(interactions::Interactions; output=:edges, min_reads=5, max_fdr=0.05, include_pvalues=true, include_relative_positions=true,
                        annotate_type_from_position=true, merge_type="merge", cds_type="CDS", five_type="5UTR", three_type="3UTR")
    out_df = copy(interactions.edges)
    filter_index = (out_df[!, :nb_ints] .>= min_reads)
    "fdr" in names(out_df) && (filter_index .&= (out_df[!, :fdr] .<= max_fdr))
    out_df = out_df[filter_index, :]
    if output === :edges
        out_df[!, :meanleft1] = Int.(round.(out_df[!, :meanleft1]))
        out_df[!, :meanright1] = Int.(round.(out_df[!, :meanright1]))
        out_df[!, :meanleft2] = Int.(round.(out_df[!, :meanleft2]))
        out_df[!, :meanright2] = Int.(round.(out_df[!, :meanright2]))
        out_df[!, :meanlength1] = Int.(round.(out_df[!, :meanlength1]))
        out_df[!, :meanlength2] = Int.(round.(out_df[!, :meanlength2]))
        out_df[!, :meanmiss1] = round.(out_df[!, :meanmiss1], digits=4)
        out_df[!, :meanmiss2] = round.(out_df[!, :meanmiss2], digits=4)
        out_df[!, :nms1] = round.(out_df[!, :nms1], digits=4)
        out_df[!, :nms2] = round.(out_df[!, :nms2], digits=4)
        out_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        out_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        out_df[:, :ref1] = interactions.nodes[out_df[!,:src], :ref]
        out_df[:, :ref2] = interactions.nodes[out_df[!,:dst], :ref]
        out_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        out_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        out_df[:, :strand1] = interactions.nodes[out_df[!,:src], :strand]
        out_df[:, :strand2] = interactions.nodes[out_df[!,:dst], :strand]
        out_df[:, :in_libs] = sum(eachcol(out_df[!, interactions.replicate_ids]))
        out_columns = [:name1, :type1, :ref1, :name2, :type2, :ref2, :nb_ints, :nb_multi, :in_libs, :strand1,
        :meanleft1, :meanright1, :meanleft2, :meanright2, :strand2, :meanlength1, :meanlength2, :minleft1, :maxright1,
        :minleft2, :maxright2, :nms1, :nms2, :meanmiss1, :meanmiss2]
        include_pvalues && (out_columns = [out_columns[1:9]..., :p_value, :fdr, out_columns[10:end]...])
        include_relative_positions && (out_columns = [out_columns..., :relmean1, :relmean2, :relmin1, :relmax1, :relmin2, :relmax2])
        if annotate_type_from_position
            ("relmean1" in names(out_df)) &&( "relmean2" in names(out_df)) || throw(AssertionError("Cannot annotate type, please run addrelativepositions! first."))
            type1_merge = out_df[:, :type1] .=== merge_type
            type1_5utr = type1_merge .& (out_df[:, :relmean1] .=== 0.0)
            type1_3utr = type1_merge .& (out_df[:, :relmean1] .=== 1.0)
            type1_cds = type1_merge .& .!(type1_3utr .| type1_5utr)
            type2_merge = out_df[:, :type2] .=== merge_type
            type2_5utr = type2_merge .& (out_df[:, :relmean2] .=== 0.0)
            type2_3utr = type2_merge .& (out_df[:, :relmean2] .=== 1.0)
            type2_cds = type2_merge .& .!(type2_3utr .| type2_5utr)
            out_df[type1_5utr, :type1] .= five_type
            out_df[type1_3utr, :type1] .= three_type
            out_df[type1_cds, :type1] .= cds_type
            out_df[type2_5utr, :type2] .= five_type
            out_df[type2_3utr, :type2] .= three_type
            out_df[type2_cds, :type2] .= cds_type
        end
        return sort(out_df[!, out_columns], :nb_ints; rev=true)
    elseif output === :nodes
        out_nodes = copy(interactions.nodes)
        for (i,row) in enumerate(eachrow(out_nodes))
            row[:nb_ints] = sum(out_df[(out_df.src .== i) .| (out_df.dst .== i), :nb_ints])
        end
        return sort(out_nodes[!, [:name, :type, :ref, :nb_single, :nb_ints]], :nb_single; rev=true)
    elseif output === :stats
        stats_df = DataFrame(Matrix(interactions.edgestats)', ["5UTR_1", ["$(i)_1" for i in 1:100]..., "3UTR_1", "5UTR_2", ["$(i)_2" for i in 1:100]..., "3UTR_2"])[filter_index, :]
        stats_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        stats_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        stats_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        stats_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        stats_df[:, :nb_ints] = out_df[:, :nb_ints]
        if annotate_type_from_position
            ("relmean1" in names(out_df)) &&( "relmean2" in names(out_df)) || throw(AssertionError("Cannot annotate type, please run addrelativepositions! first."))
            type1_merge = stats_df[:, :type1] .=== merge_type
            type1_5utr = type1_merge .& (out_df[:, :relmean1] .=== 0.0)
            type1_3utr = type1_merge .& (out_df[:, :relmean1] .=== 1.0)
            type1_cds = type1_merge .& .!(type1_3utr .| type1_5utr)
            type2_merge = stats_df[:, :type2] .=== merge_type
            type2_5utr = type2_merge .& (out_df[:, :relmean2] .=== 0.0)
            type2_3utr = type2_merge .& (out_df[:, :relmean2] .=== 1.0)
            type2_cds = type2_merge .& .!(type2_3utr .| type2_5utr)
            stats_df[type1_5utr, :type1] .= five_type
            stats_df[type1_3utr, :type1] .= three_type
            stats_df[type1_cds, :type1] .= cds_type
            stats_df[type2_5utr, :type2] .= five_type
            stats_df[type2_3utr, :type2] .= three_type
            stats_df[type2_cds, :type2] .= cds_type
        end
        return sort(stats_df, :nb_ints; rev=true)
    else
        throw(AssertionError("output has to be one of :edges, :nodes, :stats"))
    end
end