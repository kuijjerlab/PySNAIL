library(NormExpression)
library(preprocessCore)

qstats <- function(object, group_factor, window = 0.05)
  {
  # Compute sample quantiles
  Q = apply(object, 2, sort) 
  
  # Compute quantile reference
  Qref = rowMeans(Q)  
  
  # Compute SST
  SST = rowSums((Q - Qref)^2)
  
  # Compute SSB
  if(is.factor(group_factor)){
    X = model.matrix(~ 0 + group_factor)
  } else {
    X = model.matrix(~ group_factor)
  }
  
  QBETAS = t(solve(t(X) %*% X) %*% t(X) %*% t(Q))
  Qhat = QBETAS %*% t(X)
  SSB = rowSums((Qhat - Qref)^2)
  
  # Compute weights
  roughWeights = 1 - SSB / SST
  roughWeights[SST < 1e-6] = 1
  
  # Compute smooth weights
  k = floor(window * nrow(Q))
  if (k %% 2 == 0) k = k + 1
  smoothWeights = runmed(roughWeights, k = k, endrule="constant")
  
  list(Q=Q, Qref=Qref, Qhat=Qhat, SST=SST, SSB=SSB, SSE=SST-SSB, 
       roughWeights=roughWeights, smoothWeights=smoothWeights)
}

qsmooth <- function(object, group_factor,
                    batch = NULL, norm_factors = NULL, 
                    window = 0.05)
{

    object <- as.matrix(object)

    if(ncol(object) != length(group_factor)){
        stop("Number of columns in object does not match length of group_factor.")
    } 

    # Scale normalization step
    if(!is.null(norm_factors)) { 
        object <- sweep(object, 2, norm_factors, FUN = "/") 
    }

    # Compute quantile statistics
    qs <- qstats(object=object, group_factor=group_factor, window=window)
    Qref <- qs$Qref 
    Qhat <- qs$Qhat  
    w <- qs$smoothWeights
    
    # Weighted quantiles
    objectNorm = w * Qref + (1 - w) * Qhat

    # Re-order objectNorm by rank of object (columnwise)
    for (i in seq_len(ncol(objectNorm))) {
        # Grab ref. i
        ref = objectNorm[,i]
    
        # Grab object column i
        x = object[,i]
    
        # Grab ranks of x (using min rank for ties)
        rmin = rank(x, ties.method="min")
    
        # If x has rank ties then average the values of ref at 
        # those ranks
        dups = duplicated(rmin)
    
        if (any(dups)) {
            # Grab ranks of x (using random ranks for ties) 
            # (needed to uniquely identify the indices of tied ranks)
            rrand = rank(x, ties.method="random")
            # Grab tied ranks
            tied.ranks = unique(rmin[dups])
            for (k in tied.ranks) {
                # Select the indices of tied ranks 
                sel = rrand[rmin == k] 
                ref[sel] = ave(ref[sel])
            }
        }
    
    # Re-order ref and replace in objectNorm
    objectNorm[,i] = ref[rmin]
    }

    rownames(objectNorm) = rownames(object)
    colnames(objectNorm) = colnames(object)
    return(objectNorm)
}

args <- commandArgs(trailingOnly=TRUE) 
count <- args[1]
group <- args[2]

prefix <- substring(count, 1, nchar(count) - 9)
count <- read.table(count, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
group  <- read.table(group, sep='\t', header=FALSE, stringsAsFactors=FALSE)
rownames(group) <- group[, 1]
group <- group[colnames(count), ]

# DESeq
filename <- paste0(prefix, 'rle.tsv')
sf <- getFactors(count, method = "DESeq")
norm <- t(t(count) * sf)
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)
write.table(
    norm, 
    filename,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

# TMM
filename <- paste0(prefix, 'tmm.tsv')
sf <- getFactors(count, method = "TMM")
norm <- t(t(count) * sf)
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)
write.table(
    norm, 
    filename,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

# Qsmooth
filename <- paste0(prefix, 'qsmooth.tsv')
norm <- qsmooth(count, group[, 2])
#norm <- qsmoothData(qsmooth_model)
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)
write.table(
    norm, 
    filename,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)

# Quantile
filename <- paste0(prefix, 'quantile.tsv')
norm <- normalize.quantiles(as.matrix(count))
rownames(norm) <- rownames(count)
colnames(norm) <- colnames(count)
write.table(
    norm, 
    filename,
    sep='\t', quote=FALSE, row.names=TRUE, col.names=NA
)