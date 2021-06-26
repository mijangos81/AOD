# allele frequency difference (AFD) from Berner 2019
AFD_fun <- function(allele_matrix) {
    total_count <-  colSums(allele_matrix,na.rm = T)
    allele_freq <- allele_matrix/total_count
    AFD <- sum(abs(allele_freq[,1] - allele_freq[,2]))/2
  return(AFD)
}

shannon <- function(allele_matrix) {
  column_a <- allele_matrix[, 1]
  column_a <- column_a[column_a > 0]
  U_a <- sum(column_a)
  pi_a <- column_a / U_a
  est_a <- -sum(pi_a * log(pi_a))
  column_b <- allele_matrix[, 2]
  column_b <- column_b[column_b > 0]
  U_b <- sum(column_b)
  pi_b <- column_b / U_b
  est_b <- -sum(pi_b * log(pi_b))
  return(list(est_a, est_b))
}

mutual_information <- function(mat) {
  EstMLEFun <- function(mat) {
    # MLE
    entropyFun <- function(p) {
      p <- p[p > 0]
      out <- -sum(p * log(p))
      return(out)
    }
    n <- sum(mat)
    prob.hat <- mat / n
    px.hat <- apply(prob.hat, 1, sum)
    py.hat <- apply(prob.hat, 2, sum)
    I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
    # MLE of Mutual Information!
    return(I.hat)
  }
  mydata <- as.matrix(mat)
  est <- EstMLEFun(mydata)
  return(est)
}

shannon_pops <- function(alleles) {
  vector_alleles <- alleles
  vector_alleles <- vector_alleles[vector_alleles > 0]
  U <- sum(vector_alleles)
  pi <- vector_alleles / U
  est <- -sum(pi * log(pi))
  return(est)
}
