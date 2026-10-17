# ============================================================================ #
#   THEME_THESIS.R — Shared figure style of the thesis (Chapters 1-4):         #
#   page widths, type sizes, greys, palettes, the ggplot2 themes               #
#   theme_thesis() and theme_thesis_map(), the colour scales and the save      #
#   helper. The rules behind it are in STYLE.md.                               #
#                                                                              #
#   Canonical copy: Thesis/figure-style/theme_thesis.R. Each chapter keeps a   #
#   verbatim copy next to its figure code (sync.sh copies and checks them),    #
#   so a chapter's replication package needs nothing outside it. Edit the      #
#   canonical copy only, then run sync.sh.                                     #
# ============================================================================ #

THESIS_STYLE_VERSION <- "1 (2026-10-17)"

suppressPackageStartupMessages(library(ggplot2))

# ---- Page geometry -------------------------------------------------------------
# Every chapter manuscript (12pt article, A4, 2.75 cm margins) and the thesis
# (A4, 3 + 2.5 cm margins) have a text width of 155 mm. A figure is saved at
# the width it is printed at and included at that width, so LaTeX never
# rescales it and the type sizes below are the printed sizes.

FIG_W_FULL   <- 155   # \includegraphics[width=\linewidth]
FIG_W_NARROW <- 116   # \includegraphics[width=0.75\linewidth]

# ---- Type ------------------------------------------------------------------------
# Printed sizes in pt. Font: the device's sans (Helvetica or Arial).

FONT_PT       <- 8    # axis and legend titles, strip text, panel titles
FONT_SMALL_PT <- 7    # tick labels, legend text, labels inside a panel
FONT_TAG_PT   <- 9    # panel letters (bold a, b, c)

# geom_text()/annotate() sizes are in mm: FONT_LABEL prints at FONT_SMALL_PT
FONT_LABEL <- FONT_SMALL_PT / .pt

# ---- Greys -----------------------------------------------------------------------

INK        <- "#4A4A4A"   # all text
INK_MUTED  <- "#787878"   # secondary text, reference lines and their labels
GREY_LINE  <- "#C0C0C0"   # axis lines, ticks, connectors
GREY_GRID  <- "#E8E8E8"   # grid lines
GREY_STRIP <- "#F0F0F0"   # facet strips
GREY_FILL  <- "#5D6D7E"   # bars and points of a single series

# ---- Colour ----------------------------------------------------------------------
# One main contrast pair for the thesis. Each chapter gives it one meaning and
# keeps that meaning in every figure (STYLE.md lists them).

BLUE <- "#2C7FB8"
RED  <- "#C0392B"

# Categorical order for more than two series: assign in this order, never
# skip or cycle. Passes the colour-vision checks for all pairs on white
# (deutan/protan dE >= 9, normal vision dE >= 15). Orange has 2.8:1 contrast
# on white, so series in orange also carry a marker or line type.
PAL_CAT <- c(blue = BLUE, red = RED, purple = "#4A3AA7",
             green = "#1B9E77", orange = "#E08214")

# Ordinal ramp in slate, dark to light (e.g. k = 2, 4, 5). Its two ends are
# the neutral pair for a grouping that must not look like the main pair.
PAL_ORDINAL <- c("#1D3C47", "#4D7A8A", "#93B5C2")
PAL_PAIR_NEUTRAL <- c("#1D3C47", "#6F95A6")

# Diverging: negative red, zero light grey, positive blue
DIV_MID <- "#F2F2F2"

# Sequential: light grey to one colour of the main pair
pal_sequential <- function(n, to = RED) colorRampPalette(c("#E6E6E6", to))(n)

scale_fill_diverging <- function(..., limits = c(-1, 1)) {
  scale_fill_gradient2(low = RED, mid = DIV_MID, high = BLUE, midpoint = 0,
                       limits = limits, ...)
}

# ---- Themes ----------------------------------------------------------------------
# grid: the value axis carries the grid lines ("y" for vertical bars and
# lines, "x" for horizontal bars and dot plots, "none" for heatmaps).

theme_thesis <- function(base_size = FONT_PT, base_family = "",
                         grid = c("y", "x", "none")) {
  grid  <- match.arg(grid)
  small <- FONT_SMALL_PT / FONT_PT
  grid_line <- element_line(colour = GREY_GRID, linewidth = 0.25)
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title      = element_text(size = rel(1), colour = INK, hjust = 0,
                                      margin = margin(b = 4)),
      plot.subtitle   = element_text(size = rel(small), colour = INK_MUTED,
                                      hjust = 0, margin = margin(b = 4)),
      plot.caption     = element_text(size = rel(small), colour = INK_MUTED,
                                      hjust = 1, margin = margin(t = 4)),
      plot.tag         = element_text(size = rel(FONT_TAG_PT / FONT_PT),
                                      face = "bold", colour = INK,
                                      hjust = 0, vjust = 1),
      plot.tag.position = "topleft",

      axis.title       = element_text(size = rel(1), colour = INK),
      axis.title.x     = element_text(margin = margin(t = 4)),
      axis.title.y     = element_text(margin = margin(r = 4), angle = 90),
      axis.text        = element_text(size = rel(small), colour = INK),
      axis.line        = element_line(colour = GREY_LINE, linewidth = 0.25),
      axis.ticks       = element_line(colour = GREY_LINE, linewidth = 0.25),
      axis.ticks.length = unit(2, "pt"),

      legend.position   = "bottom",
      legend.title      = element_text(size = rel(1), colour = INK),
      legend.text       = element_text(size = rel(small), colour = INK),
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.key.size   = unit(8, "pt"),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.margin     = margin(t = 2),

      panel.grid.major.x = if (grid == "x") grid_line else element_blank(),
      panel.grid.major.y = if (grid == "y") grid_line else element_blank(),
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = NA, colour = NA),
      panel.border       = element_blank(),
      panel.spacing      = unit(8, "pt"),
      plot.background    = element_rect(fill = "white", colour = NA),

      strip.text       = element_text(size = rel(1), colour = INK,
                                      margin = margin(2, 2, 2, 2)),
      strip.background = element_rect(fill = GREY_STRIP, colour = NA),

      plot.margin      = margin(4, 4, 4, 4)
    )
}

theme_thesis_map <- function(base_size = FONT_PT, base_family = "") {
  theme_thesis(base_size, base_family, grid = "none") %+replace%
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank())
}

# Panel letters of cowplot::plot_grid(): bold lower case at the top left
TAG_ARGS <- list(label_size = FONT_TAG_PT, label_fontface = "bold",
                 label_colour = INK, label_x = 0, hjust = 0, vjust = 1)

# ---- Saving ----------------------------------------------------------------------
# Vector PDF through Cairo (embeds the fonts and draws any Unicode); a .png
# name gives a 400 dpi raster for slides. Width = the printed width.

save_thesis_fig <- function(plot, file, width_mm = FIG_W_FULL, height_mm) {
  png_out <- grepl("\\.png$", file, ignore.case = TRUE)
  ggsave(filename = file, plot = plot, width = width_mm, height = height_mm,
         units = "mm", bg = "white",
         device = if (png_out) "png" else cairo_pdf,
         dpi = if (png_out) 400 else 300)
}
