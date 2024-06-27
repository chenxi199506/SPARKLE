#' Cell Phenotype Network Plot
#'
#' This function generates a network plot to visualize cell phenotype interactions based on the provided interaction data.
#' The plot highlights significant interactions with edges colored and weighted according to the p-values and odds ratios (OR).
#'
#' @param interaction.best.info Data frame containing the best interaction information. It should have columns `celltype`, `Pvalue`, and `Beta`.
#' @param cell.color Character vector specifying the colors for the cells. Default is `c("darkred","#b91e45","#f9bcbd", "gray","#c4e5ef","#4198b9","#2d4e76")`.
#' @param cell.size Numeric vector specifying the sizes for the cells. Default is `c(25,20,15,10)`.
#' @param edge.width Numeric value specifying the width of the edges. Default is `3`.
#' @param arrow.size Numeric value specifying the size of the arrows. Default is `0.75`.
#' @param arrow.width Numeric value specifying the width of the arrows. Default is `1.5`.
#' @param edge.color Character vector specifying the colors for the edges based on p-values and OR. Default is `c("darkred","#b91e45","#f9bcbd",  "#c4e5ef","#4198b9","#2d4e76")`.
#'
#' @return A plot of the cell phenotype network.
#' @export
#'

Cell_phenotype_network_plot2  <- function(cwas.test=NULL,best.model=NULL,interaction.best.model=NULL,cell.color=c("#cd3023","#b91e45","#f9bcbd", "gray","#9BBBE1","#90C4E9","#8095CE"),cell.size=c(25,20,15,14),edge.width=c(2,3,4),arrow.size=0.75,arrow.width=1.5,edge.color=c("darkred","#b91e45","#f9bcbd",  "#9BBBE1","#4198b9","#2d4e76")){


  if(is.null(best.model)){
    if(is.null(cwas.test)){stop("Please input caws data  data!",  call. = FALSE)}
    best.all <- cwas_allmodel_cal(cwas.test)
    best.model <-   cwas_autoselected_model(best.all)
    best.info <- best.model[["Chosen_model_info"]]
  }else{
    best.info <- best.model[["Chosen_model_info"]]
  }

  if(is.null(interaction.best.model)){
    interaction <- cwas_2celltype_allmodel_cal(cwas.test,interaction = T)
    interaction.best.model <- cwas_autoselected_model(interaction)
    interaction.best.info <- interaction.best.model[["Chosen_model_info"]]
  }else{
    interaction.best.info <- interaction.best.model[["Chosen_model_info"]]
  }

  # 分割 celltype 列以获取节点
  edges <- strsplit(as.character(interaction.best.info$celltype), " ")
  edges <- do.call(rbind, edges)

  edges <- as.data.frame(edges)
  # 创建一个边数据框

  edges_df <- data.frame(from = edges[,2], to = edges[,1], Pvalue = interaction.best.info$Pvalue, OR = interaction.best.info$Beta * best.info$Beta)

  # 根据 Pvalue 分配权重和颜色
  edges_df$weight <- ifelse(edges_df$Pvalue < 0.001, edge.width[1],
                            ifelse(edges_df$Pvalue < 0.01, edge.width[2],
                                   ifelse(edges_df$Pvalue < 0.05, edge.width[3], 0.1)))

  edges_df$color <- ifelse((edges_df$Pvalue < 0.001) & (edges_df$OR > 0), edge.color[1],
                           ifelse((edges_df$Pvalue < 0.01) & (edges_df$OR > 0), edge.color[2],
                                  ifelse((edges_df$Pvalue < 0.05) & (edges_df$OR > 0), edge.color[3],
                                         ifelse((edges_df$Pvalue < 0.001) & (edges_df$OR < 0), edge.color[4],
                                                ifelse((edges_df$Pvalue < 0.01) & (edges_df$OR < 0), edge.color[5],
                                                       ifelse((edges_df$Pvalue < 0.05) & (edges_df$OR < 0), edge.color[6], NA))))))

  # 过滤掉权重为0的边
  #edges_df <- edges_df[edges_df$weight > 0, ]


  cellall <- unique(edges_df[,1],edges_df[,2])


  # 创建有向图对象
  g <- igraph::graph_from_data_frame(d = edges_df, directed = TRUE)

  # 设置边的宽度和颜色
  igraph::E(g)$width <- edges_df$weight
  igraph::E(g)$color <- edges_df$color
  # 设置节点的颜色和大小

  igraph::V(g)$color <-   ifelse((best.info$Pvalue < 0.001) & (best.info$Beta > 0), cell.color[1],
                                 ifelse((best.info$Pvalue < 0.01) & (best.info$Beta > 0), cell.color[2],
                                        ifelse((best.info$Pvalue < 0.05) & (best.info$Beta > 0), cell.color[3],
                                               ifelse((best.info$Pvalue < 0.001) & (best.info$Beta < 0), cell.color[7],
                                                      ifelse((best.info$Pvalue < 0.01) & (best.info$Beta < 0), cell.color[6],
                                                             ifelse((best.info$Pvalue < 0.05) & (best.info$Beta < 0), cell.color[5], cell.color[4]))))))

  igraph::V(g)$size <- ifelse(best.info$Pvalue < 0.001, cell.size[1],
                              ifelse(best.info$Pvalue < 0.01, cell.size[2],
                                     ifelse(best.info$Pvalue < 0.05, cell.size[3], cell.size[4])))  # 用 Pvalue 设置大小

  # 标准化颜色



  igraph::V(g)$shape <- "circle"
  igraph::V(g)$label.cex <- 1.2
  igraph::V(g)$label.color <- "black"
    igraph::V(g)$frame.color <- "white"

    # 绘制网络图
    # 绘制网络图
    p1 <- plot(g,
               edge.width = igraph::E(g)$width,
               edge.color = adjustcolor(igraph::E(g)$color, alpha.f = 0.7),
               edge.arrow.size = arrow.size,  # 增大箭头的大小
               edge.arrow.width = arrow.width,  # 增大箭头的宽度
               vertex.label.cex = igraph::V(g)$label.cex,
               vertex.label.color = igraph::V(g)$label.color,
               vertex.size = igraph::V(g)$size,
               vertex.color = adjustcolor(igraph::V(g)$color, alpha.f = 0.9),
               vertex.shape = igraph::V(g)$shape,
               vertex.frame.color = igraph::V(g)$frame.color,
               main = "Cell phenotype network plot",
               layout = igraph::layout_with_fr)  # 使用 Fruchterman-Reingold 布局算法

    # 添加图例
    graphics::legend("topright",
                     legend = c("Synergistic effect (p < 0.001)", "Synergistic effect (p < 0.01)", "Synergistic effect (p < 0.05)",
                                "Antagonistic effect (p < 0.05)", "Antagonistic effect (p < 0.01)", "Antagonistic effect (p < 0.001)"),
                     col = edge.color,
                     lty = 1,
                     lwd = 2,
                     cex = 0.8,
                     bg = "#e8e9ec")
    #

    graphics::legend("bottomright",
                     legend = c("Cell rate ↑↑↑ (p < 0.001)", "Cell rate ↑↑ (p < 0.01)", "Cell rate ↑ (p < 0.05)",
                                "No significant change",
                                "Cell rate ↓ (p < 0.05)", "Cell rate ↓↓ (p < 0.01)", "Cell rate ↓↓↓ (p < 0.001)"),
                     fill = cell.color,
                     title = "Cell proportion changes with phenotype",
                     cex = 0.8,
                     bg = "#e8e8ea")

    return(g)

}

