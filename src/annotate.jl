function annotate_utrs!(annotations::Dict{String, DataFrame}, tss::Dict{String, DataFrame}, terms::Dict{String, DataFrame}; 
    max_distance=150, add_distance=150, guess_distance=150)
    for (chr, annotation) in annotations
        annotation[!, :fiveUTR] = annotation[!, :start] .- guess_distance
        annotation[!, :fiveType] = fill("guess", nrow(annotation))
        annotation[!, :threeUTR] = annotation[!, :stop] .+ guess_distance
        annotation[!, :threeType] = fill("guess", nrow(annotation))
        for row in eachrow(annotation)
            five_hits = tss[chr][row[:start]-max_distance .<= tss[chr][!,:pos] .<= row[:start], :]
            three_hits = terms[chr][row[:stop] .<= terms[chr][!,:pos] .<= row[:stop]+max_distance, :]
            isempty(five_hits) && (five_hits = tss[chr][row[:start]-max_distance-add_distance .<= tss[chr][!,:pos] .<= row[:start]-max_distance, :])
            isempty(three_hits) && (three_hits = terms[chr][row[:stop]+max_distance .<= terms[chr][!,:pos] .<= row[:stop]+max_distance+add_distance, :])
            isempty(five_hits) || (row[:fiveUTR]=five_hits[argmax(five_hits[!, :val]), :pos]; row[:fiveType]="max")
            isempty(three_hits) || (row[:threeUTR]=three_hits[argmax(three_hits[!, :val]), :pos]; row[:threeType]="max")
            #(row[:stop]==-372) && (println(three_hits))
        end
    end
end

strand_filter(a::Interval, b::Interval) = strand(a) == strand(b)

function annotate!(alns::Alignments, features::Features)
    for alignment in alns
        for part in alignment
            for feature_interval in eachoverlap(features.list, part.ref, filter=strand_filter)
                push!(part.ref.metadata, feature_interval.metadata)
            end
        end
    end
end

function annotate!(alns::PairedAlignments, features::Features)
    for (alignment1, alignment2) in alns
        for part in alignment1
            for feature in eachoverlap(features.list, part.ref; filter=strand_filter)
                push!(part.ref.metadata, feature.metadata)
            end
        end
        for part in alignment2
            for feature in eachoverlap(features.list, part.ref; filter=strand_filter)
                push!(part.ref.metadata, feature.metadata)
            end
        end
    end
end