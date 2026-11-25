Fst_global <- function (x, pop = NULL, quiet = TRUE, na.alleles = "", global = TRUE) 
{
  # Use pegas::: for internal functions when running outside the package namespace
  if (any(pegas:::.checkPloidy(x) != 2)) 
    stop("Fst() requires diploid data")
  NAMESX <- names(x)
  if (is.null(pop)) {
    ipop <- which(NAMESX == "population")
    if (!length(ipop)) 
      stop("no 'population' column in x")
  }
  else {
    if (is.numeric(pop) && length(pop) == 1) {
      ipop <- pop
    }
    else {
      x$populationforthisanalysis <- factor(pop)
      ipop <- length(x)
    }
  }
  LOCI <- attr(x, "locicol")
  nloci <- length(LOCI)
  
  # Store per-locus sums of a, b, and c instead of dividing immediately
  abc <- matrix(NA_real_, nloci, 3)
  dimnames(abc) <- list(NAMESX[LOCI], c("sum_a", "sum_b", "sum_c"))
  
  for (j in 1:nloci) {
    if (!quiet) 
      cat("\rAnalyzing locus", j, "/", nloci)
    Z <- x[, c(LOCI[j], ipop)]
    Z <- na.omit(Z, na.alleles = na.alleles)
    N <- nrow(Z)
    nBYpop <- tabulate(Z[[2]]) # Safely reference the population column
    r <- length(nBYpop)
    
    # Skip locus if fewer than 2 populations have non-missing data
    if (r < 2 || any(nBYpop == 0)) next
    
    nbar <- N/r
    nC <- (N - sum(nBYpop^2)/N)/(r - 1)
    ALLELES <- pegas::getAlleles(Z)[[1]]
    h <- p <- matrix(0, r, length(ALLELES))
    for (i in 1:r) {
      s <- summary(Z[as.integer(Z[[2]]) == i, ])[[1]]
      allel <- names(s$allele)
      genot <- names(s$genotype)
      p[i, ] <- s$allele
      for (k in seq_along(allel)) {
        for (l in seq_along(genot)) {
          ag <- unlist(strsplit(genot[l], "/"))
          if (sum(ag %in% allel[k]) == 1) 
            h[i, k] <- h[i, k] + s$genotype[l]
        }
      }
    }
    ptild <- p/(2 * nBYpop)
    pbar <- colSums(p)/(2 * N)
    s2 <- colSums(nBYpop * (ptild - rep(pbar, each = r))^2)/((r - 1) * nbar)
    hbar <- colSums(h)/N
    A <- pbar * (1 - pbar) - (r - 1) * s2/r
    a <- nbar * (s2 - (A - hbar/4)/(nbar - 1))/nC
    b <- nbar * (A - (2 * nbar - 1) * hbar/(4 * nbar))/(nbar - 1)
    c <- hbar/2
    
    # Save summed variance components for locus j
    abc[j, ] <- c(sum(a), sum(b), sum(c))
  }
  if (!quiet) 
    cat("... Done.\n")
  
  if (global) {
    # Ratio of averages across all valid loci
    tot_a <- sum(abc[, "sum_a"], na.rm = TRUE)
    tot_b <- sum(abc[, "sum_b"], na.rm = TRUE)
    tot_c <- sum(abc[, "sum_c"], na.rm = TRUE)
    
    out <- c(
      Fit = 1 - tot_c / (tot_a + tot_b + tot_c),
      Fst = tot_a / (tot_a + tot_b + tot_c),
      Fis = 1 - tot_c / (tot_b + tot_c)
    )
    return(out)
  } else {
    # Original per-locus output + raw variance components
    res <- cbind(
      Fit = 1 - abc[, "sum_c"] / rowSums(abc),
      Fst = abc[, "sum_a"] / rowSums(abc),
      Fis = 1 - abc[, "sum_c"] / (abc[, "sum_b"] + abc[, "sum_c"]),
      abc
    )
    return(res)
  }
}
