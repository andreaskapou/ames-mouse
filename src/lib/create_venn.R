create_quad_venn <- function(n1, n2, n3, n4, filename = "venn_diagram.png"){
  library(VennDiagram)

  n1 <- read.table(file = n1)
  n2 <- read.table(file = n2)
  n3 <- read.table(file = n3)
  n4 <- read.table(file = n4)

  venn.diagram(
    x = list(
      DO = n1$V1,
      DY = n2$V1,
      NO = n3$V1,
      NY = n4$V1
    ),
    filename = filename,
    col = "black",
    imagetype = "png",
    lty = "dotted",
    lwd = 4,
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", "white",
                  "white", "white", "white", "darkblue", "white",
                  "white", "white", "white", "darkgreen", "white"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 2.5,
    cat.fontfamily = "serif"
  )
}

create_triple_venn <- function(n1, n2, n3, filename = "venn_diagram.png"){
  library(VennDiagram)
  
  n1 <- read.table(file = n1)
  n2 <- read.table(file = n2)
  n3 <- read.table(file = n3)
  
  venn.diagram(
    x = list(
      K562 = n1$V1,
      GM12878 = n2$V1,
      H1_hESC = n3$V1
    ),
    filename = filename,
    col = "black",
    lty = "dotted",
    fill = c("red2", "blue", "green"),
    imagetype = "png",
    alpha = 0.5,
    label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "text",
    cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 2.5,
    cat.fontfamily = "serif",
    cat.dist = c(0.172, 0.182, -0.14),
    cat.pos = c(350, 10, 0)
  )
}