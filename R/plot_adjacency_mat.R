
#' Visualize graph based on adjacency matrix
#' @description calls DiagrammeR
#' @param adj_mat (numerical matrix) from: col -> to: row
#' @param title   (character) plot title
#' @param draw_edge_labels (TRUE/FALSE) if FALSE, remove edge labels
#' @param node_labels (numeric or character vector) node labels
#' @param round_digit (integer) maximum digits
#' @import DiagrammeR
#' @export
plot_adjacency_mat <- function(adj_mat,
                               title = NULL,
                               draw_edge_labels = TRUE,
                               node_labels = NULL,
                               round_digit = 2) {
    nodes <- DiagrammeR::create_node_df(ncol(adj_mat), label = node_labels)
    graph <- DiagrammeR::create_graph(nodes, directed = TRUE)
    adj_mat <- round(adj_mat, round_digit)
    for (i in 1:ncol(adj_mat)) {
        for (j in which(abs(adj_mat[, i]) > 0)) {
            if (draw_edge_labels) {
                graph <- DiagrammeR::add_edge(graph, from = i, to = j, edge_aes = DiagrammeR::edge_aes(label = adj_mat[j, i]))
            }
            else {
                graph <- DiagrammeR::add_edge(graph, from = i, to = j)
            }

        }
    }
    dot <- DiagrammeR::render_graph(graph, layout = "dot", title = title)
    dot$x$diagram <- gsub("neato", "dot", dot$x$diagram)
    print(DiagrammeR::grViz(dot$x$diagram))
}
