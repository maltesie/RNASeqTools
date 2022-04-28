function dashboard(reads::Sequences)

    app = dash(assets_folder=joinpath(@__DIR__, "assets"))

    app.layout = html_div(id="root") do
        html_div(id="left-column") do
            html_h1("Conserved UTRs"),
            dcc_input(id = "my-id1", value="initial value", type = "text"),
            html_div(id = "my-div1")
        end,
        html_div(id="app-container") do
            dcc_tabs(id ="tabs", value="data") do
                dcc_tab(id="data", label="Data", value="data") do
                    dash_datatable(id="my-data-table", columns=["test"], data=datatabele_data(show_table), page_action="none",
                                    style_table=Dict("height"=>"600px", "overflowY"=>"auto")),
                    html_p("asdasdasdasdasdasdasdasda")
                end,
                dcc_tab(id="about", label="About", value="about") do
                    html_p("asdasdasdasd")
                end
            end
        end
    end

    run_server(app, "0.0.0.0", 8083)
end

function lengthhist(reads::Sequences)
    lengths = [length(read) for read in reads]
    histogram(lengths, legend=false, title=title)
end

function nucleotidedist(reads::Sequences; align=:left, normalize=true, title="nucleotides of reads")
    count = nucleotide_count(reads; normalize=normalize)
    dna_trans = Dict(DNA_A => "A", DNA_T=>"T", DNA_G=>"G", DNA_C=>"C", DNA_N=>"N")
    label = reshape([dna_trans[key] for key in keys(count)], (1, length(dna_trans)))
    plot(collect(values(count)), label=label, title=title)
end

function expressionpca(features::Features, samples::Vector{Coverage}, conditions::Dict{String, UnitRange{Int64}}; plot_pcs=(1,2), legend=:best, topcut=500)
    averages = normalizedcount(features, samples)
    topcut = min(size(averages, 1), topcut)
    averages = averages[sortperm(vec(var(averages; dims=2)); rev=true)[1:topcut],:]
    averages = log2.(averages .+ 1)
    averages = (averages .- mean(averages; dims=1)) ./ std(averages, dims=1)
    M = fit(PCA, averages)
    atrans = MultivariateStats.transform(M, averages)
    p = plot()
    for (n, r) in sort(conditions, by=x->first(last(x)))
        scatter!(atrans[first(plot_pcs),r], atrans[last(plot_pcs),r], label=n)
    end
    ratios = principalvars(M) ./ tvar(M)
    ratio1 = round(ratios[first(plot_pcs)] * 100, digits=1)
    ratio2 = round(ratios[last(plot_pcs)] * 100, digits=1)
    plot(p, xlabel="PC$(first(plot_pcs)) ($ratio1%)", ylabel="PC$(last(plot_pcs)) ($ratio2%)", legend=legend)
end

function expressionpca(counts_file::String, conditions::Dict{String, UnitRange{Int64}}; plot_pcs=(1,2), legend=:best, topcut=500)
    averages = CSV.read(counts_file, DataFrame; header=1, delim=',') |> Matrix{Float64}
    topcut = min(nrow(averages), topcut)
    averages .+= 0.1
    (nfeatures, nsamples) = size(averages)
    avg_sample::Vector{Float64} = [geomean(averages[i, :]) for i in 1:nfeatures]
    norm_factors = [median(averages[:, i] ./ avg_sample) for i in 1:nsamples]
    averages ./= norm_factors'
    averages = averages[sortperm(vec(var(averages; dims=2)); rev=true)[1:topcut],:]
    averages = log2.(averages)
    M = fit(PCA, averages)
    atrans = MultivariateStats.transform(M, averages)
    p = plot()
    for (n, r) in sort(conditions, by=x->first(last(x)))
        scatter!(atrans[first(plot_pcs),r], atrans[last(plot_pcs),r], label=n)
    end
    ratios = principalvars(M) ./ tvar(M)
    ratio1 = round(ratios[first(plot_pcs)] * 100, digits=1)
    ratio2 = round(ratios[last(plot_pcs)] * 100, digits=1)
    plot(p, xlabel="PC$(first(plot_pcs)) ($ratio1%)", ylabel="PC$(last(plot_pcs)) ($ratio2%)", legend=legend)
end

"""
KronaTools wrapper function.
Uses report file from kraken2 created by align_kraken2().

Example use (using .report.txt file created by align_kraken2()):
    kronaplot("notex_01_1.report.txt")

Output:
    notex_01_1.error.txt
    notex_01_1.krona.html
"""
function kronaplot(taxonomy_file::String;
        krona_bin="ktImportTaxonomy",
    )

    output_file = split(taxonomy_file, ".")[1] * ".krona.html"
    # error_file = split(taxonomy_file, ".")[1] * ".error.txt"
    params = [
              "-o", output_file,
              "-t", 2, "-m", 1
             ]

    tmp_file = tempname()
    taxonomy_file = readdlm(taxonomy_file, '\t')
    # cut taxonomy file to appropriate columns
    writedlm(tmp_file, hcat(taxonomy_file[:, 3], taxonomy_file[:, 7]), '\t')

    # ktImportTaxonomy krona.in -o krona.html -t 2 -m 1
    cmd = pipeline(`$krona_bin $tmp_file $params`)# , stderr=error_file)
    run(cmd)
end
