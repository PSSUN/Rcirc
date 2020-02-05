##' Draw the length distribution of the circRNA according to the given bed file
##' User can customize the upper limit of the length
##'
##' @import ggplot2
##' @export
##' @param bed bed file
##' @param max_length the max length of circRNA showed in plot

showLength <- function(bed, max_length) {
  bed <- read.csv(bed)
  len <- bed$End - bed$Start
  less_max <- len[len < max_length]
  m <- mean(less_max)
  ave <- gsub("%", floor(m), "Average length:
% bp")
  ti <- gsub("%", max_length, "length\n(less than %)")

  g1 <- ggplot(as.data.frame(less_max), aes(x = less_max, y = ..density..)) +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "lightblue") +
    geom_line(stat = "density", colour = "#660066") +
    expand_limits(y = 0) +
    geom_vline(aes(xintercept = m), linetype = 5, colour = "darkred") +
    labs(x = ti, title = "Distribution of circRNA length") +
    theme(plot.title = element_text(hjust = 0.5))

  g2 <- ggplot(as.data.frame(less_max), aes(x = less_max)) +
    labs(x = ti, title = "Distribution of circRNA length") +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "lightblue") +
    theme(plot.title = element_text(hjust = 0.5))
  ggplot2.two_y_axis <- function(g1, g2) {
    g1 <- ggplotGrob(g1)
    g2 <- ggplotGrob(g2)

    # Get the location of the plot panel in g1.
    # These are used later when transformed elements of g2 are put back into g1
    pp <- c(subset(g1$layout, name == "panel", se = t:r))

    # Overlap panel for second plot on that of the first plot
    g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

    # Taken from the cowplot package:
    # https://github.com/wilkelab/cowplot/blob/master/R/

    hinvert_title_grob <- function(grob) {

      # Swap the widths
      widths <- grob$widths
      grob$widths[1] <- widths[3]
      grob$widths[3] <- widths[1]
      grob$vp[[1]]$layout$widths[1] <- widths[3]
      grob$vp[[1]]$layout$widths[3] <- widths[1]

      # Fix the justification
      grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust
      grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust
      grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
      grob
    }

    # Get the y axis title from g2
    index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?
    ylab <- g2$grobs[[index]] # Extract that grob
    ylab <- hinvert_title_grob(ylab) # Swap margins and fix justifications

    # Put the transformed label on the right side of g1
    g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
    g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "ylab-r")

    # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
    index <- which(g2$layout$name == "axis-l") # Which grob
    yaxis <- g2$grobs[[index]] # Extract the grob

    # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
    # The relevant grobs are contained in axis$children:
    #   axis$children[[1]] contains the axis line;
    #   axis$children[[2]] contains the tick marks and tick mark labels.

    # First, move the axis line to the left
    yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))

    # Second, swap tick marks and tick mark labels
    ticks <- yaxis$children[[2]]
    ticks$widths <- rev(ticks$widths)
    ticks$grobs <- rev(ticks$grobs)

    # Third, move the tick marks
    ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")

    # Fourth, swap margins and fix justifications for the tick mark labels
    ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])

    # Fifth, put ticks back into yaxis
    yaxis$children[[2]] <- ticks

    # Put the transformed yaxis on the right side of g1
    g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
    g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "axis-r")
    grid.newpage()
    grid.draw(g1)
  }
  ggplot2.two_y_axis(g2, g1)
}

