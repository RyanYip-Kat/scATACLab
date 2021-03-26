library(ggplot2)
archr_theme<-function (color = "black", textFamily = "sans", baseSize = 10,
    baseLineSize = 0.5, baseRectSize = 0.5, plotMarginCm = 1,
    legendPosition = "bottom", legendTextSize = 5, axisTickCm = 0.1,
    xText90 = FALSE, yText90 = FALSE)
{
    theme <- theme_bw() + theme(text = element_text(family = textFamily),
        axis.text = element_text(color = color, size = baseSize),
        axis.title = element_text(color = color, size = baseSize),
        title = element_text(color = color, size = baseSize),
        plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm,
            plotMarginCm), "cm"), panel.background = element_rect(fill = "transparent",
            colour = NA), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,
            color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1,
                "points"), "mm"))), axis.ticks.length = unit(axisTickCm,
            "cm"), axis.ticks = element_line(color = color, size = baseLineSize *
            (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"),
            "mm"))), legend.key = element_rect(fill = "transparent",
            colour = NA), legend.text = element_text(color = color,
            size = legendTextSize), legend.box.background = element_rect(color = NA),
        legend.position = legendPosition, strip.text = element_text(size = baseSize,
            color = "black"))
    if (xText90) {
        theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90,
            hjust = 1))
    }
    if (yText90) {
        theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90,
            vjust = 1))
    }
    return(theme)
}

