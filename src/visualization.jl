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
                    dash_datatable(id="my-data-table", columns=datatable_columns(show_table), data=datatabele_data(show_table), page_action="none",
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

function lengthhist(reads::PairedSequences, title="length of reads")
    lengths = vcat([[length(read1) length(read2)] for (read1, read2) in reads]...)
    histogram(lengths, labels=["read1" "read2"], title=title)
end

function nucleotidedist(reads::Sequences; align=:left, normalize=true, title="nucleotides of reads")
    count = nucleotide_count(reads; normalize=normalize)
    dna_trans = Dict(DNA_A => "A", DNA_T=>"T", DNA_G=>"G", DNA_C=>"C", DNA_N=>"N")
    label = reshape([dna_trans[key] for key in keys(count)], (1, length(dna_trans)))
    plot(collect(values(count)), label=label, title=title)
end

function nucleotidedist(reads::PairedSequences; align=:left, normalize=true, title1="nucleotides of read1", title2="nucleotides of read2")
    (count1, count2) = nucleotide_count(reads; normalize=normalize)
    dna_trans = Dict(DNA_A => "A", DNA_T=>"T", DNA_G=>"G", DNA_C=>"C", DNA_N=>"N")
    label = reshape([dna_trans[key] for key in keys(count1)], (1, length(dna_trans)))
    p1 = plot(collect(values(count1)), label=label, title=title1)
    p2 = plot(collect(values(count2)), legend=false, title=title2)
    plot(p1, p2, layout=(2,1))
end

function similarityhist(reads::PairedSequences; window_size=10, step_size=5, title="similarity of reads")
    sim = similarity(reads; window_size=window_size, step_size=step_size)
    histogram(collect(values(sim)), title=title, legend=false)
end