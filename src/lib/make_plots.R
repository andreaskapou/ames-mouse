plot_scatt_regr_test <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.2, 8.01)
    xlim=c(-3.2, 8.01)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }

  mycorr = round(cor(output$test_pred, output$test$y),3)  	# correlation
  plot(output$test$y,
       output$test_pred,
       ylab="predicted expression (log2)",
       xlab="measured expression (log2)",
       main=main_lab,
       ylim=ylim,
       xlim=xlim,
       cex.main=1.4,
       cex.lab=1.1,
       cex.axis=1.2,
       col="#0000ff22",
       pch=16,
       cex=1.4)
  #my_lm = lm(output$test$Y ~ output$test_pred)
  my_lm = lm(output$test_pred ~ output$test$y)
  a = coefficients(my_lm)
  abline(a[1], a[2], lwd=2, col="red", lty=2)
  l=legend("topleft", paste0("R=", mycorr), bty="n", cex=2)  #no border
  #text(x=l$text$x, y=l$text$y-l$rect$h/2, adj=c(0,1), cex=1.2)
}


plot_scatt_regr_train <- function(output, main_lab = "Methylation Profile", is_margins = FALSE){
  if (is_margins){
    ylim=c(-3.2, 8.01)
    xlim=c(-3.2, 8.01)
  }else{
    max_data <- max(output$test_pred, output$test$y)
    min_data <- min(output$test_pred, output$test$y)
    ylim=c(min_data, max_data)
    xlim=c(min_data, max_data)
  }
  mycorr = round(cor(output$train_pred, output$train$y),3)  	# correlation
  plot(output$train$y,
       output$train_pred,
       ylab="predicted expression (log2)",
       xlab="measured expression (log2)",
       main=main_lab,
       ylim=ylim,
       xlim=xlim,
       cex.main=1.4,
       cex.lab=1.1,
       cex.axis=1.2,
       col="#0000ff22",
       pch=16,
       cex=1.4)
  #my_lm = lm(output$train$Y ~ output$train_pred)
  my_lm = lm(output$train_pred ~ output$train$y)
  a = coefficients(my_lm)
  abline(a[1], a[2], lwd=2, col="red", lty=2)
  l=legend("topleft", paste0("R=", mycorr), bty="n", cex=2)  #no border
  #text(x=l$text$x, y=l$text$y-l$rect$h/2, adj=c(0,1), cex=1.2)
}




plot_cluster_prof <- function(model, basis, add_clust = FALSE, main_lab = "Clustered methylation profiles"){
  xs <- seq(-1,1,len=2000) # create some values
  plot(x=xs, y=eval_probit_function(basis, xs, model$w[,1]), xlim=c(-1,1), ylim=c(0,1), 
       type="l", col="salmon3", lwd=3, xlab="promoter region", ylab="methylation level",
       main=main_lab)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,2]),
        col="black", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,3]),
        col="blue", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,4]),
        col="red3", lwd=3)
  lines(x=xs, y=eval_probit_function(basis, xs, model$w[,5]),
        col="darkgreen", lwd=3)
  
  if (add_clust){
    lines(x=xs, y=eval_probit_function(basis, xs, model$w[,6]),
          col="cyan", lwd=3)
    lines(x=xs, y=eval_probit_function(basis, xs, model$w[,7]),
          col="darkgoldenrod1", lwd=3)
  }
}


plot_cluster_box <- function(gene_expr, add_clust = FALSE, main_lab = "Gene expression levels"){
  col <- c("salmon3", "black", "blue", "red3", "darkgreen")
  if (add_clust){
    col <- c(col, "cyan", "darkgoldenrod1")
  }
  boxplot(gene_expr, col=col, notch=T, xlab="Cluster K", 
          ylab="expression level", main=main_lab)
}
