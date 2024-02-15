#' Plot deconstructed difference score correlation
#'
#' Plots the slopes for y1 and y2 by x, and a slope for y1-y2 by x for comparison.
#'
#' @param ddsc_object An object produced by ddsc_sem function.
#' @param diff_color Character. Color for difference score (y1-y2). Default "black".
#' @param y1_color Character. Color for difference score component y1. Default "turquoise".
#' @param y2_color Character. Color for difference score component y2. Default "orange".
#' @param x_label Character. Label for variable X. If NULL (default), variable name is used.
#' @param y_labels Character vector. Labels for variable y1 and y2. If NULL (default), variable names are used.
#' @param densities Logical. Are y-variable densities plotted? Default TRUE.
#' @param point_alpha Numeric. Opacity for data points (default 0.50)
#' @param dens_alpha Numeric. Opacity for density distributions (default 0.75)
#' @param col_widths Numeric vector. Widths of the plot columns: slope figures and density figures; default c(3, 1).
#' @param row_heights Numeric vector. Heights of the plot rows: components, difference score, slope coefs; default c(2, 1, 0.5).
#' @param coef_locations Numeric vector. Locations for printed coefficients. Quantiles of the range of x-variable. Default c(0, 1/3, 2/3).
#' @param coef_names Character vector. Names of the printed coefficients. Default c("b_11", "b_21", "r_x_y1-y2").
#'
#' @export
#'
#' @examples
#' set.seed(342356)
#' d <- data.frame(
#'   y1 = rnorm(50),
#'   y2 = rnorm(50),
#'   x = rnorm(50)
#' )
#' fit<-ddsc_sem(
#'   data = d, y1 = "y1", y2 = "y2",
#'   x = "x"
#' )
#'
#' plot_ddsc(fit,x_label = "X",
#'           y_labels=c("Y1","Y2"))
#'
plot_ddsc <- function(ddsc_object,
                      diff_color = "black",
                      y1_color = "turquoise",
                      y2_color = "orange",
                      x_label = NULL,
                      y_labels = NULL,
                      densities = TRUE,
                      point_alpha = 0.5,
                      dens_alpha = 0.75,
                      col_widths = c(3, 1),
                      row_heights = c(2, 1, 0.5),
                      coef_locations = c(0/3, 1/3, 2/3),
                      coef_names = c("b_11", "b_21", "r_x_y1-y2")) {

  comb_plot <- NULL
  x <- NULL
  y <- NULL
  DV_type <- NULL
  phrase <- NULL

  results <- ddsc_object$results
  data <- ddsc_object$data
  x_scaled <- rownames(ddsc_object$descriptives)[6]
  y1_scaled <- rownames(ddsc_object$descriptives)[2]
  y2_scaled <- rownames(ddsc_object$descriptives)[4]

  if (is.null(x_label)) {
    x_label <- paste0("X: ", x_scaled)
  }

  if (is.null(y_labels)) {
    y1_label <- paste0("Y1: ", y1_scaled)
    y2_label <- paste0("Y2: ", y2_scaled)
  } else {
    y1_label <- y_labels[1]
    y2_label <- y_labels[2]
  }
  # paste0("rh[", expression(x~","~y[1]), "] = ")
  # paste0("hr[", expression(xy1), "] = ")
  # paste0("hr[x,y1] = ")

  # coefs<-data.frame(
  #  DV_type=c("Y1","Y2","Y1-Y2"),
  #  phrase=c(
  #    paste0("hr[\"x,y1\"]"),
  #    "hr2",
  #    "dsc"
  #  ))

  # coefs<-data.frame(
  #  DV_type=c("Y1","Y2","Y1-Y2"),
  #  phrase=c(
  #    paste0("hr[x_y1]",deparse(" = 43")),
  #    "hr2",
  #    "dsc"
  #  ))


  # obtain the coefficients for the plot
  coefs <- data.frame(
    DV_type = c("Y1", "Y2", "Y1-Y2"),
    phrase = c(
      paste0(
        coef_names[1]," = ", round(results[which(rownames(results) == "b_11"), "est"], 2),
        ", ",
        ifelse(results[which(rownames(results) == "b_11"), "pvalue"] < .001, "p < .001",
          paste0("p = ", round(results[which(rownames(results) == "b_11"), "pvalue"], 3))
        )
      ),
      paste0(
        coef_names[2]," = ", round(results[which(rownames(results) == "b_21"), "est"], 2),
        ", ",
        ifelse(results[which(rownames(results) == "b_21"), "pvalue"] < .001, "p < .001",
          paste0("p = ", round(results[which(rownames(results) == "b_21"), "pvalue"], 3))
        )
      ),
      paste0(
        coef_names[3]," = ", round(results[which(rownames(results) == "r_xy1_y2"), "est"], 2),
        ", ",
        ifelse(results[which(rownames(results) == "r_xy1_y2"), "pvalue"] < .001, "p < .001",
          paste0("p = ", round(results[which(rownames(results) == "r_xy1_y2"), "pvalue"], 3))
        )
      )
    )
  )

  # transform to long format

  pd_long <-
    data.frame(
      y = c(
        data[, y1_scaled],
        data[, y2_scaled],
        data[, "diff_score_scaled"]
      ),
      DV_type = rep(c("Y1", "Y2", "Y1-Y2"), each = nrow(data)),
      x = rep(data[, x_scaled], times = 3)
    )
  # plots

  p1 <-
    ggplot2::ggplot(
      data = pd_long[pd_long$DV_type != "Y1-Y2", ],
      ggplot2::aes(x = x, y = y, color = DV_type)
    ) +
    ggplot2::geom_smooth(method = "lm", formula = "y~x") +
    ggplot2::scale_color_manual(
      values = c(y1_color, y2_color),
      labels = c(
        paste0(y1_label),
        paste0(y2_label)
      )
    ) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::ylab(paste0(
      "Component score\n",
      "(", y1_label, " or ", y2_label, ")"
    )) +
    # xlab(x_label)+
    ggplot2::xlab(NULL) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "top",
      panel.background = ggplot2::element_blank()
    )

  p2 <- ggplot2::ggplot(
    data = pd_long[pd_long$DV_type == "Y1-Y2", ],
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_smooth(method = "lm", formula = "y~x", color = diff_color) +
    ggplot2::geom_point(alpha = point_alpha, color = diff_color) +
    ggplot2::ylab(paste0(
      "Difference score\n",
      "(", y1_label, "\U2212", y2_label, ")"
    )) +
    ggplot2::xlab(x_label) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "top"
    )

  # obtain locations for printed coefficients
  x_range <-
    max(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x.range) -
    min(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x.range)

  coef_x_locations <-
    c(
      min(data[, x_scaled]) + x_range * coef_locations[1],
      min(data[, x_scaled]) + x_range * coef_locations[2],
      min(data[, x_scaled]) + x_range * coef_locations[3]
    )

  p3 <- ggplot2::ggplot(
    data = coefs,
    ggplot2::aes(
      label = phrase,
      y = 0,
      x = coef_x_locations,
      hjust = "left"
    )
  ) +
    ggplot2::geom_text(
      color =
        c(
          unique(ggplot2::ggplot_build(p1)$data[[1]]$colour),
          unique(ggplot2::ggplot_build(p2)$data[[1]]$colour)
        ), parse = FALSE
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::xlim(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x.range)

  # obtain the limits for density distribution
  y1y2_min_max <-
    ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$y.range

  # plot the densities for Y1 and Y2
  suppressMessages(
    p1_dens <-
      ggplot2::ggplot(
        data = pd_long[pd_long$DV_type != "Y1-Y2", ],
        ggplot2::aes(x = y, fill = DV_type)
      ) +
      ggplot2::geom_density(alpha = dens_alpha) +
      ggplot2::scale_fill_manual(
        values = c(y1_color, y2_color),
        labels = c(
          paste0(y1_label),
          paste0(y2_label)
        )
      ) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      # ylab("Density")+
      ggplot2::coord_cartesian(xlim = y1y2_min_max) +
      ggplot2::theme(
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        # text=element_text(size=16,  family="sans"),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      ) +
      ggplot2::coord_flip()
  )


  # ggpubr::ggarrange(p1,p1_dens,ncol=2,
  #                  widths=c(3,1),
  #                  common.legend = TRUE,align = "h")



  # obtain the limits for density distribution
  diff_min_max <-
    ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$y.range

  # plot the densities for Y1 and Y2
  suppressMessages(
    p2_dens <-
      ggplot2::ggplot(
        data = pd_long[pd_long$DV_type == "Y1-Y2", ],
        ggplot2::aes(x = y, fill = diff_color)
      ) +
      ggplot2::geom_density(alpha = dens_alpha, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = diff_color) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Density") +
      ggplot2::coord_cartesian(xlim = diff_min_max) +
      ggplot2::theme(
        legend.position = "none",
        legend.title = ggplot2::element_blank(),
        # text=element_text(size=16,  family="sans"),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      ) +
      ggplot2::coord_flip()
  )


  # ggpubr::ggarrange(p2,p2_dens,ncol=2,
  #                  widths=c(3,1),
  #                  common.legend = TRUE,align = "h")


  # check which of the two density plots has higher maximum
  dens_minmax <-
    c(
      min(
        ggplot2::ggplot_build(p2_dens)$layout$panel_params[[1]]$x.range,
        ggplot2::ggplot_build(p1_dens)$layout$panel_params[[1]]$x.range
      ),
      max(
        ggplot2::ggplot_build(p2_dens)$layout$panel_params[[1]]$x.range,
        ggplot2::ggplot_build(p1_dens)$layout$panel_params[[1]]$x.range
      )
    )
  # arranging the plots
  suppressMessages(
    if (densities) {
      comb_plot <- ggpubr::ggarrange(p1, p1_dens + ggplot2::coord_flip(ylim = dens_minmax),
        p2, p2_dens + ggplot2::coord_flip(ylim = dens_minmax),
        p3,
        ncol = 2, nrow = 3,
        widths = col_widths,
        heights = row_heights,
        common.legend = TRUE, align = "hv"
      )
    } else {
      comb_plot <-
        ggpubr::ggarrange(p1, p2, p3,
          ncol = 1,
          heights = row_heights,
          common.legend = TRUE, align = "hv"
        )
    }
  )


  return(comb_plot)
}
