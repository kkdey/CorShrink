
data(mtcars)
corr <- round(cor(mtcars), 1)
corr


method = "square"
type = "full"
ggtheme = ggplot2::theme_minimal
title = ""
show.legend = TRUE
legend.title = "Corr"
show.diag = FALSE 
colors = c("blue", "white", "red")
outline.color = "gray" 
hc.order = FALSE
hc.method = "complete"
lab = FALSE
lab_col = "black"
lab_size = 4
p.mat = NULL
sig.level = 0.05
insig = "pch"
pch = 4
pch.col = "black"
pch.cex = 5
tl.cex = 12
tl.col = "black"
tl.srt = 45