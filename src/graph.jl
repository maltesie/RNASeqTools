mutable struct Interactions <: InteractionContainer
    graph::SimpleDiGraph
    nodes::DataFrame
    edges::DataFrame
    replicate_ids::Vector{Symbol}
end

"""
<<<<<<< HEAD
Method of write function which saves the Interactions struct in a jld2 file.
=======
Dispatch of Base write function which saves the Interactions struct in a jld2 file.
>>>>>>> 24aa817496d9d29f5ef6d9dfc8ab808dd49f39aa
"""
function Base.write(filepath::String, interactions::Interactions)
    if !endswith(filepath, ".jld2")
        throw(ArgumentError("Append '.jld2' to filepath"))
    else
        save(filepath, "interactions", interactions)
    end
<<<<<<< HEAD
=======
    save(filepath, "interactions", interactions)
>>>>>>> 24aa817496d9d29f5ef6d9dfc8ab808dd49f39aa
end

function leftpos(alnpart::AlignedPart, alnread::AlignedRead)
    minimum(leftposition(p) for p in alnread if sameannotation(p, alnpart))
end

function rightpos(alnpart::AlignedPart, alnread::AlignedRead)
    maximum(rightposition(p) for p in alnread if sameannotation(p, alnpart))
end

function Base.append!(interactions::Interactions, alignments::Alignments, replicate_id::Symbol; min_distance=1000, filter_types=[])
    myhash(part::AlignedPart) = hash(name(part))
    interactions.edges[:, replicate_id] = repeat([0], nrow(interactions.edges))
    push!(interactions.replicate_ids, replicate_id)
    trans = Dict{UInt, Int}(interactions.nodes[i, :hash]=>i for i in 1:nrow(interactions.nodes))
    trans_edges = Dict{Tuple{Int,Int},Int}((interactions.edges[i, :src],interactions.edges[i, :dst])=>i for i in 1:nrow(interactions.edges))
    for alignment in alignments
        !isempty(filter_types) && typein(alignment, filter_types) && continue
        is_chimeric = ischimeric(alignment; min_distance=min_distance)
        is_multi = is_chimeric ? ismulti(alignment) : false

        for (i,part) in enumerate(alignment)
            hasannotation(part) || continue
            any(samename(part, formerpart) for formerpart in alignment.alns[1:i-1]) && continue
            h = myhash(part)
            if !(h in keys(trans))
                trans[h] = length(trans) + 1
                add_vertex!(interactions.graph)
                push!(interactions.nodes, (name(part), type(part) in ("5UTR", "3UTR") ? "CDS" : type(part), refname(part), 0, 0, strand(part), h))
            end
            is_chimeric || (interactions.nodes[trans[h], :nb_single] += 1)
        end

        for (part1, part2) in combinations(alignment.alns[collect(!any(samename(alignment[i], formerpart) for formerpart in alignment.alns[1:i-1]) for i in 1:length(alignment))], 2)
            (hasannotation(part1) && hasannotation(part2)) || continue
            ischimeric(part1, part2; min_distance=min_distance) || continue
            a, b = trans[myhash(part1)], trans[myhash(part2)]
            interactions.nodes[a, :nb_ints] += 1
            interactions.nodes[b, :nb_ints] += 1
            has_edge(interactions.graph, a, b) || (add_edge!(interactions.graph, a, b); trans_edges[(a,b)] = ne(interactions.graph))
            iindex = trans_edges[(a, b)]
            left1, right1 = leftpos(part1, alignment), rightpos(part1, alignment)
            left2, right2 = leftpos(part2, alignment), rightpos(part2, alignment)
            iindex > nrow(interactions.edges) &&
                push!(interactions.edges, (a, b, 0, 0, left1, right1, left2, right2, 0, 0, 0, 0, 0, 0, (0 for i in 1:length(interactions.replicate_ids))...))
            interactions.edges[iindex, :nb_ints] += 1
            is_multi && (interactions.edges[iindex, :nb_multi] += 1)
            interactions.edges[iindex, :minleft1] = min(interactions.edges[iindex, :minleft1], left1)
            interactions.edges[iindex, :maxright1] = max(interactions.edges[iindex, :maxright1], right1)
            interactions.edges[iindex, :minleft2] = min(interactions.edges[iindex, :minleft2], left2)
            interactions.edges[iindex, :maxright2] = max(interactions.edges[iindex, :maxright2], right2)
            for (s,v) in zip((:meanleft1, :meanright1, :meanleft2, :meanright2, :length1, :length2),
                            (left1, right1, left2, right2, right1 - left1 + 1, right2 - left2 + 1))
                interactions.edges[iindex, s] =
                (interactions.edges[iindex, s] * (interactions.edges[iindex, :nb_ints] - 1) + v) / interactions.edges[iindex, :nb_ints]
            end
            interactions.edges[iindex, replicate_id] = 1
        end
    end

    return interactions
end

function Interactions()
    nodes = DataFrame(:name=>String[], :type=>String[], :ref=>String[], :nb_single=>Int[], :nb_ints=>Int[], :strand=>Char[], :hash=>UInt[])
    edges = DataFrame(:src=>Int[], :dst=>Int[], :nb_ints=>Int[], :nb_multi=>Int[], :minleft1=>Int[], :maxright1=>Int[], :minleft2=>Int[],
                        :maxright2=>Int[], :length1=>Float64[], :length2=>Float64[], :meanleft1=>Float64[], :meanleft2=>Float64[],
                        :meanright1=>Float64[], :meanright2=>Float64[])
    return Interactions(SimpleDiGraph(), nodes, edges, [])
end

function Interactions(alignments::Alignments; replicate_id=:first, min_distance=1000, filter_types=[])
    append!(Interactions(), alignments, replicate_id, min_distance=min_distance, filter_types=filter_types)
end

"""
Load Interactions struct from jld2 file.
"""
Interactions(filepath::String) = load(filepath, "interactions")

function annotate!(interactions::Interactions, features::Features; method=:disparity)
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
    interactions.edges = hcat(interactions.edges, DataFrame(repeat([-Inf -Inf -Inf -Inf -Inf -Inf], nrow(interactions.edges)),
                                                            [:relmean1, :relmean2, :relmin1, :relmin2, :relmax1, :relmax2]))
    tus = Dict(name(feature)=>(leftposition(feature), rightposition(feature)) for feature in features if !(type(feature) in ("5UTR", "3UTR")))
    for edge_row in eachrow(interactions.edges)
        (feature1_left, feature1_right) = tus[interactions.nodes[edge_row[:src], :name]]
        (feature2_left, feature2_right) = tus[interactions.nodes[edge_row[:dst], :name]]
        ispositive1 = interactions.nodes[edge_row[:src], :strand] === '+'
        ispositive2 = interactions.nodes[edge_row[:dst], :strand] === '+'
        p1 = collect(edge_row[[:meanright1, :minleft1, :maxright1]])
        p2 = collect(edge_row[[:meanleft2, :minleft2, :maxright2]])
        (relpos1, relmin1, relmax1) = min.(1.0, max.(0.0, (ispositive1 ? (p1 .- feature1_left) : (feature1_right .- p1)) ./ (feature1_right .- feature1_left)))
        (relpos2, relmin2, relmax2) = min.(1.0, max.(0.0, (ispositive2 ? (p2 .- feature2_left) : (feature2_right .- p2)) ./ (feature2_right .- feature2_left)))
        edge_row[[:relmean1, :relmean2, :relmin1, :relmin2, :relmax1, :relmax2]] = round.((relpos1,relpos2, relmin1, relmin2, relmax1, relmax2); digits=4)
    end
    return interactions
end

function asdataframe(interactions::Interactions; output=:edges, min_interactions=5, max_fdr=0.05)
    out_df = copy(interactions.edges)
    "fdr" in names(out_df) && filter!(:fdr => <=(max_fdr), out_df)
    filter!(:nb_ints => >=(min_interactions), out_df)
    if output === :edges
        out_df[!, :meanleft1] = Int.(floor.(out_df[!, :meanleft1]))
        out_df[!, :meanright1] = Int.(floor.(out_df[!, :meanright1]))
        out_df[!, :meanleft2] = Int.(floor.(out_df[!, :meanleft2]))
        out_df[!, :meanright2] = Int.(floor.(out_df[!, :meanright2]))
        out_df[!, :length1] = Int.(floor.(out_df[!, :length1]))
        out_df[!, :length2] = Int.(floor.(out_df[!, :length2]))
        out_df[:, :name1] = interactions.nodes[out_df[!,:src], :name]
        out_df[:, :name2] = interactions.nodes[out_df[!,:dst], :name]
        out_df[:, :ref1] = interactions.nodes[out_df[!,:src], :ref]
        out_df[:, :ref2] = interactions.nodes[out_df[!,:dst], :ref]
        out_df[:, :type1] = interactions.nodes[out_df[!,:src], :type]
        out_df[:, :type2] = interactions.nodes[out_df[!,:dst], :type]
        out_df[:, :strand1] = interactions.nodes[out_df[!,:src], :strand]
        out_df[:, :strand2] = interactions.nodes[out_df[!,:dst], :strand]
        out_df[:, :in_libs] = sum(eachcol(out_df[!, interactions.replicate_ids]))
        return sort(out_df[!, [:name1, :type1, :ref1, :name2, :type2, :ref2, :nb_ints, :nb_multi, :p_value, :fdr, :in_libs, :strand1,
                            :meanleft1, :meanright1, :meanleft2, :meanright2, :strand2, :length1, :length2, :minleft1, :maxright1,
                            :minleft2, :maxright2, :relmean1, :relmean2, :relmin1, :relmin2, :relmax1, :relmax2]], :nb_ints; rev=true)
    elseif output === :nodes
        out_nodes = copy(interactions.nodes)
        for (i,row) in enumerate(eachrow(out_nodes))
            row[:nb_ints] = sum(out_df[(out_df.src .== i) .| (out_df.dst .== i), :nb_ints])
        end
        return sort(out_nodes[!, [:name, :type, :ref, :nb_single, :nb_ints]], :nb_single; rev=true)
    end
end
