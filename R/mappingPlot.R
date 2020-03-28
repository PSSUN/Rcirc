# ---
# *title*
# plot circular RNA

# *param*
# summary: result of function circ_summary
# circ_index: choose which circ RNA to plot, input index
# genomefile: virtual genome path
# upstream: visualize upstream distance of junction site
# downstream: visualize downstream distance of junction site
# style: visualization style: track or aside

# *return*
# a plot of target circular RNA

# *details*
# users can get a target circular RNA visualization by passing summary parameters to our function. To be frindly,
# we supply function of zoom in and zoom out if users specify the upstream and downstream parameters and plot complete
# circular RNA by default

# *example*
#
# *author*
# Wang Haoming https://github.com/AlexWanghaoming
# ---

##' @import circlize
##' @export
##' @author Haoming Wang
showMapping <-
  function(summary, circ_index, genomefile, upstream = NULL, downstream = NULL, style = "track", reads_col = "darkgreen") {
    library(circlize)
    trimmed.gr <- summary$trimmed_circRNA
    trimmed.reads <- summary$trimmed_reads
    circ_view <- summary$circ_view

    test_circ <- trimmed.gr[circ_index]
    # extract reads overlap with target circular RNA
    ovl.reads <-
      subsetByOverlaps(trimmed.reads, test_circ, ignore.strand = TRUE)

    # extract coverage info of target circular RNA regions
    cir.cov <-
      as.vector(circ_view[[circ_index]])

    # junction sites on virtual genome
    js <- test_circ$junction_site
    circ.st <- start(test_circ)
    js.pos <- js - circ.st

    if (is.null(upstream) && is.null(downstream)) {
      upstream <- js - start(test_circ)
      downstream <- end(test_circ) - js + 1
    } else {
      if (upstream > js - start(test_circ)) {
        warning(
          "Upstream is beyond the boundary of target circular RNA, adjust automaticly!",
          immediate. = T
        )
      }
      if (downstream > end(test_circ) - js + 1) {
        warning(
          "Downstream is beyond the boundary of target circular RNA, adjust automaticly!",
          immediate. = T
        )
      }
      if (upstream <= 0 || downstream <= 0) {
        stop("Please input positive value!")
      }
      upstream <- min(js - start(test_circ), upstream)
      downstream <- min(end(test_circ) - js + 1, downstream)
    }

    circ.length <- width(test_circ)
    zoom_length <- upstream + downstream
    zoom_cir.gr <-
      GRanges(
        1,
        ranges = IRanges(
          start = test_circ$junction_site - upstream,
          end = test_circ$junction_site +
            downstream
        )
      )

    zoom_ovl.reads <-
      subsetByOverlaps(trimmed.reads, zoom_cir.gr, ignore.strand = TRUE)

    ## trim zoom_ovl.reads heads and tails
    zoom_ovl.reads <- GenomicRanges::restrict(zoom_ovl.reads,
      start = start(zoom_cir.gr),
      end = end(zoom_cir.gr)
    )

    zoom_cir.cov <- cir.cov[(js.pos - upstream + 1):(js.pos + downstream)]

    # preparing data for heatmap track and junction site
    heat.dd <-
      data.frame(
        chr = "circ",
        start = 1:length(zoom_cir.cov) - 1,
        end = 1:length(zoom_cir.cov),
        score = zoom_cir.cov
      )
    juc.dd <-
      data.frame(
        chr = "circ",
        start = c(upstream, upstream),
        end = c(upstream, upstream),
        sorce = c(0, 1)
      )

    ## plot
    circos.clear()
    circos.par(
      "start.degree" = 90,
      cell.padding = c(0, 0, 0, 0)
    )
    circos.initialize("circ", xlim = c(0, zoom_length)) # length of circ RNA

    ## Track1: heatmap of reads density
    f1 <-
      colorRamp2(
        breaks = c(0, max(zoom_cir.cov)),
        colors = c("#A2E4F2", "#FF240A")
      )
    circos.genomicTrackPlotRegion(
      list(heat.dd, juc.dd),
      track.height = 0.1,
      ylim = c(0, 1),
      bg.border = NA,
      panel.fun = function(region, value, ...) {
        i <- getI(...)
        if (i == 1) {
          circos.genomicRect(region,
            value,
            col = f1(value[[1]]),
            border = NA,
            ...
          )
        } else {
          circos.genomicLines(region, value, lwd = 5, col = "#26532B") # plot junction sites
        }
      }
    )

    ## Track2: sequence info and start codon highlight
    # 1. sequence annotation
    genome <-
      readDNAStringSet(genomefile, format = "fasta")
    sequence <-
      as.character(DNAStringSet(genome, start(zoom_cir.gr), end(zoom_cir.gr)))
    seq <-
      substring(sequence, 1:width(zoom_cir.gr), 1:width(zoom_cir.gr))
    seq.color <-
      ifelse(seq == "A", "green", ifelse(seq == "G", "brown", ifelse(seq == "T", "red", "blue")))
    # 2. start,stop codon highlight
    start_codon.pos <- as.vector(gregexpr("ATG", sequence)[[1]])
    start_codon.pos <- start_codon.pos[start_codon.pos != -1]
    stop_codon.pos <- c(
      gregexpr2("TAA", sequence)[[1]],
      gregexpr2("TAG", sequence)[[1]],
      gregexpr2("TGA", sequence)[[1]]
    )
    stop_codon.pos <- stop_codon.pos[stop_codon.pos != -1]

    circos.track(
      ylim = c(0, 1),
      track.height = 0.05,
      bg.border = NA,
      panel.fun = function(x, y) {
        circos.text(
          1:zoom_length - 0.5,
          rep(0.2, zoom_length),
          seq,
          facing = "inside",
          cex = 0.5,
          col = seq.color
        )
        if (length(start_codon.pos) >= 1) {
          circos.rect(
            start_codon.pos - 1,
            0.5,
            start_codon.pos + 3 - 1,
            1,
            col = "green",
            border = "white"
          )
          circos.text(
            start_codon.pos + 0.5,
            rep(0.75, length(start_codon.pos)),
            rep("M", length(start_codon.pos)),
            facing = "inside",
            col = "white",
            cex = 0.5
          )
        }

        if (length(stop_codon.pos) >= 1) {
          circos.rect(
            stop_codon.pos - 1,
            0.5,
            stop_codon.pos + 3 - 1,
            1,
            col = "red",
            border = "white"
          )
          circos.text(
            stop_codon.pos + 0.5,
            rep(0.6, length(stop_codon.pos)),
            rep("*", length(stop_codon.pos)),
            facing = "inside",
            col = "white"
          )
        }
        circos.arrow(
          0,
          1 / 32 * zoom_length,
          y = 0.8,
          col = "black",
          arrow.head.length = 1 / 90 * zoom_length,
          width = 0.2
        )
      }
    )

    # Track3: plot reads coverage region
    zoom_circ.st <-
      start(zoom_cir.gr) # circ RNA start site on virtual genome
    reads.start <- start(zoom_ovl.reads) - zoom_circ.st
    reads.end <- end(zoom_ovl.reads) - zoom_circ.st
    if (style == "track") {
      y1 <- seq(99, 2, length.out = NROW(zoom_ovl.reads))
    }
    if (style == "aside") {
      track.gp <- list(1)
      if (NROW(zoom_ovl.reads) > 1) {
        for (idx in 2:NROW(zoom_ovl.reads)) {
          # print(idx)
          for (i in 1:length(track.gp)) {
            # print(track.gp)
            track <- track.gp[[i]]
            if (start(zoom_ovl.reads[idx]) > end(zoom_ovl.reads[rev(track)[1]])) {
              track.gp[[i]] <- c(track, idx)
              break
            }
            if (i == length(track.gp)) {
              track.gp[length(track.gp) + 1] <- idx
            }
          }
        }
      }

      y1 <- seq(99, 2, length.out = length(track.gp))
      for (x in 1:length(track.gp)) {
        le <- length(track.gp[[x]])
        if (le >= 2) {
          a <- track.gp[[x]][1]
          y1 <- c(y1, rep(y1[a], le - 1))
        }
      }
      y1 <- sort(y1, decreasing = T)[order(unlist(track.gp))]
    }
    y2 <- y1 + 0.5
    circos.track(
      ylim = c(0, 100),
      track.height = 0.70,
      bg.border = NA,
      panel.fun = function(x, y) {
        xlim <- CELL_META$xlim
        ylim <- CELL_META$ylim
        # plot reads coverage region
        circos.rect(reads.start,
          y1,
          reads.end,
          y2,
          col = reads_col,
          border = "white"
        )
      }
    )
    # plot reads region border
    radius <- 0.87
    angle1 <- (min(reads.start) / zoom_length) * 360
    lines(
      c(0, sin(angle1 / 180 * pi) * radius),
      c(0, cos(angle1 / 180 * pi) * radius),
      lwd = 1,
      col = "#021D56",
      lty = 2
    )
    angle2 <- (max(reads.end) / zoom_length) * 360
    lines(
      c(0, sin(angle2 / 180 * pi) * radius),
      c(0, cos(angle2 / 180 * pi) * radius),
      lwd = 1,
      col = "#021D56",
      lty = 2
    )

    p <- recordPlot()
    return(p)
  }

