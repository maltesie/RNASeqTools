using Dash
using DashHtmlComponents
using DashCoreComponents

function visualize_conserved_utrs()

    app = dash()

    app.layout = html_div((dcc_input(id = "my-id", value="initial value", type = "text"),
    html_div(id = "my-div")))

    callback!(app, Output("my-div", "children"), Input("my-id", "value")) do input_value
        "You've entered kot $(input_value)"
    end

    run_server(app, "0.0.0.0", 8083)
end

visualize_conserved_utrs()