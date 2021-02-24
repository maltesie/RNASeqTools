using Dash, DashHtmlComponents, DashCoreComponents, DashTable
using DataFrames, CSV

function datatable_columns(df::SubDataFrame)
    [(id=column, name=column) for column in names(df)[1:end-3]]
end

function datatabele_data(df::SubDataFrame)
    [Dict(name=>row[Symbol(name)] for name in names(df)[1:end-3]) for row in eachrow(df)]
end

function visualize_conserved_utrs(conservation_file::String)

    conservation_table = CSV.read(conservation_file, DataFrame)
    show_table = @view(conservation_table[!, vcat([:name, :length], [Symbol(name) for name in names(conservation_table) if endswith(name, "score")])])

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

    callback!(app, Output("my-div1", "children"), Input("my-id1", "value")) do input_value
        "You've entered kot 2 $(input_value)"
    end

    run_server(app, "0.0.0.0", 8083)
end

#visualize_conserved_utrs("/home/abc/Workspace/ConservedUTRs/alignments/three_alignment_table (copy).csv")