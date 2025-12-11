## -----------------------------------------------------------------------------
# install.packages(c("rrcov", "robustbase"), repos = "https://cloud.r-project.org")  # install if missing
library(rrcov)          # robust multivariate methods (PCA, covariance)
library(robustbase)     # robust basics (MCD, etc.)

data(machines)          # load the machines dataset
machines_sub <- machines[71:111, ]                  # slice rows hp-3000/64 ... ibm-4331-2
machines_sub$machine <- rownames(machines_sub)      # store machine names for labeling
rownames(machines_sub) <- NULL                      # drop row names to avoid confusion

n_obs <- nrow(machines_sub)                         # count observations
n_vars <- ncol(machines_sub) - 1                    # count numeric variables (exclude machine names)

head(machines_sub, 3)                               # peek at the first rows


## ----results='hide'-----------------------------------------------------------
# summarize each variable (min/mean/median/max/quartiles)
summary(machines_sub)


## -----------------------------------------------------------------------------
head(machines_sub$machine, 5)                       # show first few machine IDs


## -----------------------------------------------------------------------------
winsor_mean <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs, names = FALSE)           # lower/upper cutoffs
  mean(pmin(pmax(x, qs[1]), qs[2]))                 # clamp extremes then average
}

total_variance <- function(S) sum(diag(S))          # sum of variances (trace)
generalized_variance <- function(S) determinant(S, logarithm = TRUE)$modulus  # log-determinant


## ----results='hide'-----------------------------------------------------------
num_vars <- machines_sub[ , setdiff(names(machines_sub), "machine")]  # numeric-only data

stat_table <- data.frame(                              # assemble summary table
  variable = names(num_vars),                     # variable name
  mean = sapply(num_vars, mean),                  # arithmetic mean
  median = sapply(num_vars, median),              # median
  trimmed_mean = sapply(num_vars, mean, trim = 0.1),  # 10% trimmed mean
  winsor_mean = sapply(num_vars, winsor_mean),    # 5–95 winsorized mean
  sd = sapply(num_vars, sd),                      # standard deviation
  var = sapply(num_vars, var),                    # variance
  mad = sapply(num_vars, mad)                     # median absolute deviation
)

stat_table                                        # show table


## ----results='hide'-----------------------------------------------------------
S_classic <- cov(num_vars)                          # classical covariance matrix
total_var <- total_variance(S_classic)              # total variance (trace)
gen_var_log <- generalized_variance(S_classic)      # log generalized variance

list(
  covariance_matrix = S_classic,                  # covariance matrix
  total_variance = total_var,                     # total variance
  generalized_variance_log = gen_var_log          # log generalized variance
)                                                  # print results


## ----results='hide'-----------------------------------------------------------
md_classic <- mahalanobis(
  num_vars,                              # data matrix
  center = colMeans(num_vars),           # classic mean vector
  cov = S_classic                        # classic covariance matrix
)
cutoff <- qchisq(0.975, df = ncol(num_vars))  # 97.5% chi-square threshold

# Robust (MCD) Mahalanobis distances
cmcd <- robustbase::covMcd(num_vars)     # robust center/covariance via MCD
md_robust <- mahalanobis(
  num_vars,                              # data matrix
  center = cmcd$center,                  # robust center
  cov = cmcd$cov                         # robust covariance
)

par(mfrow = c(1, 2))                     # two plots side by side
plot(md_classic, pch = 19, main = "Mahalanobis (classic)", ylab = "Distance")  # classic distances
abline(h = cutoff, col = "red", lty = 2)                                       # cutoff line

plot(md_robust, pch = 19, main = "Mahalanobis (robust MCD)", ylab = "Distance") # robust distances
abline(h = cutoff, col = "red", lty = 2)                                       # cutoff line
par(mfrow = c(1, 1))                     # reset layout

flag_classic <- which(md_classic > cutoff)    # indices flagged by classic distance
flag_robust  <- which(md_robust  > cutoff)    # indices flagged by robust distance

out_table <- list(                            # collect flagged machines
  classical = data.frame(                                       # classic flagged table
    machine = machines_sub$machine[flag_classic],       # names flagged (classic)
    distance = round(md_classic[flag_classic], 2)       # classic distances
  ),
  robust = data.frame(                                          # robust flagged table
    machine = machines_sub$machine[flag_robust],        # names flagged (robust)
    distance = round(md_robust[flag_robust], 2)         # robust distances
  )
)
out_table                                   # show both tables


## -----------------------------------------------------------------------------
pca_raw <- prcomp(num_vars, center = TRUE, scale. = FALSE)   # PCA on raw scale
pca_std <- prcomp(num_vars, center = TRUE, scale. = TRUE)    # PCA on standardized data

# ensure numeric to avoid class quirks on sdev
sdev_raw <- as.numeric(pca_raw$sdev)                         # SDs of PCs (raw)
sdev_std <- as.numeric(pca_std$sdev)                         # SDs of PCs (std)

pve_raw <- sdev_raw^2 / sum(sdev_raw^2)                      # proportion variance explained (raw)
pve_std <- sdev_std^2 / sum(sdev_std^2)                      # proportion variance explained (std)


## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))                                          # side-by-side scree plots
plot(pve_raw * 100, type = "b", pch = 19, xlab = "PC", ylab = "% variance",
     main = "Raw scale")                                      # raw PVE
lines(cumsum(pve_raw) * 100, type = "b", col = "blue")        # cumulative PVE (raw)
abline(h = 95, col = "red", lty = 2)                          # 95% line

plot(pve_std * 100, type = "b", pch = 19, xlab = "PC", ylab = "% variance",
     main = "Standardized")                                   # standardized PVE
lines(cumsum(pve_std) * 100, type = "b", col = "blue")        # cumulative PVE (std)
abline(h = 95, col = "red", lty = 2)                          # 95% line
par(mfrow = c(1, 1))                                          # reset layout


## -----------------------------------------------------------------------------
k_raw <- which(cumsum(pve_raw) >= 0.95)[1]                    # PCs needed (raw) for 95%
k_std <- which(cumsum(pve_std) >= 0.95)[1]                    # PCs needed (std) for 95%

data.frame(
  scale = c("raw", "standardized"),                           # scale type
  pcs_needed_for_95pct = c(k_raw, k_std),                     # count of PCs to reach 95%
  cumulative_variance = c(cumsum(pve_raw)[k_raw], cumsum(pve_std)[k_std])  # achieved cum PVE
)                                                           # show summary


## -----------------------------------------------------------------------------
loadings_std <- round(pca_std$rotation[, 1:2], 3)   # loadings for PCs 1 and 2 (standardized PCA)
loadings_std


## -----------------------------------------------------------------------------
scores_std <- as.data.frame(pca_std$x)                    # PCA scores (std)
scores_std$machine <- machines_sub$machine               # add machine names

plot(scores_std$PC1, scores_std$PC2, pch = 19,           # scatter of PC1 vs PC2
     xlab = "PC1 (std)", ylab = "PC2 (std)",
     main = "Scores on standardized data")
text(scores_std$PC1, scores_std$PC2, labels = scores_std$machine,
     pos = 3, cex = 0.6)                                 # label points


## -----------------------------------------------------------------------------
xnew <- c(75, 2000, 0.8, 80000, 300, 24, 62, 47)          # injected outlier values
machines_out <- machines_sub                              # start from subset
machines_out[machines_out$machine == "hp-3000/64", names(num_vars)] <- xnew  # replace row

pca_out_classic <- prcomp(machines_out[names(num_vars)], center = TRUE, scale. = FALSE)   # classic PCA with outlier
pca_out_robust <- PcaCov(machines_out[names(num_vars)], cov.control = CovControlMcd(), scale = FALSE)  # robust PCA

pve_out_classic <- pca_out_classic$sdev^2 / sum(pca_out_classic$sdev^2)  # PVE classic
pve_out_robust <- pca_out_robust@eigenvalues / sum(pca_out_robust@eigenvalues)  # PVE robust


## ----results='hide'-----------------------------------------------------------
data.frame(
  component = seq_along(pve_out_classic),                 # component index
  classic_pct = round(pve_out_classic * 100, 2),          # classic PVE (%)
  robust_pct = round(pve_out_robust * 100, 2)             # robust PVE (%)
)


## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))                                       # two plots
plot(pca_out_classic$x[,1], pca_out_classic$x[,2], pch = 19,
     main = "Classic PCA with outlier", xlab = "PC1", ylab = "PC2")  # classic scores
text(pca_out_classic$x[,1], pca_out_classic$x[,2],
     labels = machines_out$machine, pos = 3, cex = 0.5)              # labels

plot(pca_out_robust@scores[,1], pca_out_robust@scores[,2], pch = 19,
     main = "Robust PCA (MCD)", xlab = "PC1", ylab = "PC2")          # robust scores
text(pca_out_robust@scores[,1], pca_out_robust@scores[,2],
     labels = machines_out$machine, pos = 3, cex = 0.5)              # labels
par(mfrow = c(1, 1))                                       # reset layout


## -----------------------------------------------------------------------------
plot(pca_out_robust)  # outlier map: orthogonal vs score distances

