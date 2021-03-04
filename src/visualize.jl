function dashboard(reads::Reads)

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

function hist_length_distribution(reads::Reads)
end

function hist_length_distribution(reads::PairedReads, title="length of reads")
    lengths = vcat([[length(read1) length(read2)] for (read1, read2) in values(reads)]...)
    histogram(lengths, labels=["read1" "read2"], title=title)
end

function line_nucleotide_distribution(reads::Reads)
end

function line_nucleotide_distribution(reads::PairedReads; align=:left, title1="nucleotides of read1", title2="nucleotides of read2")
    max_length = maximum(vcat([[length(read1) length(read2)] for (read1, read2) in values(reads)]...))
    count1 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    count2 = Dict(DNA_A => zeros(max_length), DNA_T=>zeros(max_length), DNA_G=>zeros(max_length), DNA_C=>zeros(max_length), DNA_N=>zeros(max_length))
    nb_reads = length(reads)
    for (read1, read2) in reads
        (align==:left) ? 
        (index1 = 1:length(read1); index2 = 1:length(read2)) : 
        (index1 = (max_length - length(read1) + 1):max_length; index2 = (max_length - length(read2) + 1):max_length)
        for ((i1, n1),(i2, n2)) in zip(zip(index1, read1), zip(index2, read2))
            count1[n1][i1] += 1
            count2[n2][i2] += 1
        end
    end
    for ((key1, c1), (key2, c2)) in zip(count1, count2)
        c1 ./= length(reads)
        c2 ./= length(reads)
    end
    dna_trans = Dict(DNA_A => "A", DNA_T=>"T", DNA_G=>"G", DNA_C=>"C", DNA_N=>"N")
    label = reshape([dna_trans[key] for key in keys(count1)], (1, length(dna_trans)))
    p1 = plot(collect(values(count1)), label=label, title=title1)
    p2 = plot(collect(values(count2)), legend=false, title=title2)
    plot(p1, p2, layout=(2,1))
end

function hist_similarity(reads::PairedReads; window_size=10, step_size=5, title="similarity of reads")
    sim = similarity(reads; window_size=window_size, step_size=step_size)
    histogram(collect(values(sim)), title=title, legend=false)
end