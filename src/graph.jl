mutable struct Interactions <: InteractionContainer
    graph::MetaGraph
    replicate_ids::Vector{Symbol}
end

function integrate!(alignments::PairedAlignments, graph::MetaGraph; replicate_id=:first, min_distance=200, filter_types=[])
    trans = Dict{Tuple{String,String}, Int}()
    for (alignment1, alignment2) in alignments
        if hasannotation(alignment1) && hasannotation(alignment2)
            !isempty(filter_types) && typein(alignment1, alignment2, filter_types) && continue
            chimeric = ischimeric(alignment1, alignment2; min_distance=min_distance)
            
            tus = union(Set((name(part),type(part)) for part in alignment1), Set((name(part),type(part)) for part in alignment2))
            nb_diff = length(Set(tu[1] for tu in tus))
            for tu in tus
                if !(tu in keys(trans))
                    trans[tu] = length(trans) + 1
                    add_vertex!(graph)
                    set_props!(graph, nv(graph), Dict(:name=>first(tu), :type=>last(tu), :nb_ints=>0, :nb_single=>0))
                end
                chimeric ?
                set_prop!(graph, trans[tu], :nb_ints, get_prop(graph, trans[tu], :nb_ints)+nb_diff-1) :
                set_prop!(graph, trans[tu], :nb_single, get_prop(graph, trans[tu], :nb_single)+1)
            end

            if chimeric
                multi = nb_diff > 2
                for (a, b) in combinations([trans[tu] for tu in tus], 2)
                    get_prop(graph, a, :name) == get_prop(graph, b, :name) && continue
                    if add_edge!(graph, a, b)
                        set_props!(graph, a, b, Dict(:nb_ints=>1, :nb_multi=> (multi ? 1 : 0), replicate_id=>1))
                    else
                        set_prop!(graph, a, b, :nb_ints , get_prop(graph, a, b, :nb_ints) + 1)
                        set_prop!(graph, a, b, :nb_multi , get_prop(graph, a, b, :nb_multi) + (multi ? 1 : 0))
                    end
                    has_prop(graph, a, b, replicate_id) || set_prop(graph, a, b, replicate_id, 1)
                end
            end
        end
    end
end

function Interactions()
    return Interactions(MetaGraph(), [])
end

function Interactions(alignments::PairedAlignments; replicate_id=:first, min_distance=200, filter_types=[])
    graph = MetaGraph()
    integrate!(alignments, graph; replicate_id=replicate_id, min_distance=min_distance, filter_types=filter_types)
    return Interactions(graph, [replicate_id])
end

function Base.append!(interactions::Interactions, alignments::PairedAlignments; replicate_id=:second, min_distance=200, filter_types=[])
    integrate!(alignments, interactions.graph; replicate_id=replicate_id, min_distance=min_distance, filter_types=filter_types)
    push!(interactions.replicate_ids, replicate_id)
end

function Base.filter!(interactions::Interactions; method=:disparity)
    @assert method in (:disparity, :fisher)
    if method === :fisher
    elseif method === :disparity
    end
    return interactions
end

function asdataframe(interactions::Interactions; output=:edges)
    if output === :edges
        frame = DataFrame(name1=repeat([""], ne(interactions.graph)), type1=repeat([""], ne(interactions.graph)), name2=repeat([""], ne(interactions.graph)), 
                        type2=repeat([""], ne(interactions.graph)), nb_chimeras=repeat([0], ne(interactions.graph)), nb_multi=repeat([0], ne(interactions.graph)))
        for (i,edge) in enumerate(edges(interactions.graph))
            frame[i, :name1] = get_prop(interactions.graph, src(edge), :name)
            frame[i, :type1] = get_prop(interactions.graph, src(edge), :type)
            frame[i, :name2] = get_prop(interactions.graph, dst(edge), :name)
            frame[i, :type2] = get_prop(interactions.graph, dst(edge), :type)
            frame[i, :nb_chimeras] = get_prop(interactions.graph, edge, :nb_ints)
            frame[i, :nb_multi] = get_prop(interactions.graph, edge, :nb_multi)
        end
        return sort(frame, :nb_chimeras; rev=true)
    elseif output === :nodes
        frame = DataFrame(name=repeat([""], nv(interactions.graph)), type=repeat([""], nv(interactions.graph)), nb_single=repeat([0], nv(interactions.graph)))
        for (i,node) in enumerate(vertices(interactions.graph))
            frame[i, :name] = get_prop(interactions.graph, node, :name)
            frame[i, :type] = get_prop(interactions.graph, node, :type)
            frame[i, :nb_single] = get_prop(interactions.graph, node, :nb_single)
        end
        return sort(frame, :nb_single; rev=true)
    end         
end