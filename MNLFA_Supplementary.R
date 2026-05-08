# ============================================================================
# MNLFA_Supplementary.R — v2.0
# Supplementary Analyses for Manuscript Revision
# ----------------------------------------------------------------------------
# CHANGES v2.0:
#   Block 4 & Block 6: sửa hc3_vcov() — thay diag(X %*% solve() %*% t(X))
#   bằng hatvalues(m) và vector W thay cho diagonal matrix.
#   Lý do: code cũ tạo ma trận N×N trong RAM, gây NaN ở leverage values
#   khi N=537, khiến toàn bộ slope/SE trong Block 4 ra NA.
# ============================================================================
#
# Mục đích: Bổ sung các phân tích còn thiếu trong bản thảo, KHÔNG chỉnh sửa
#           file MNLFA_lavaan.R gốc để đảm bảo reproducibility.
#
# Các khối phân tích:
#   Block 1 — Descriptive Statistics + Correlation Matrix
#   Block 2 — VIF Diagnostics (Multicollinearity)
#   Block 3 — Common Method Bias (Harman's Single Factor Test)
#   Block 4 — Simple Slopes Table (Publication-Ready)
#   Block 5 — Fornell-Larcker Matrix (Discriminant Validity)
#   Block 6 — Robustness Comparison Table (3-column: Main / Outlier / MM)
#
# Yêu cầu:
#   Chạy SAU khi MNLFA_lavaan.R đã chạy xong (các file Excel output đã có)
#   R packages: lavaan, writexl, car, readxl
#
# Cách chạy:
#   Rscript MNLFA_Supplementary.R <CSV_PATH> <OUT_DIR>
#
#   Ví dụ (Windows):
#   Rscript MNLFA_Supplementary.R "D:/CANH/.../Data/_mnlfa_data.csv" "D:/CANH/.../"
# ============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(writexl)
  library(readxl)
})

# ── Đọc arguments ─────────────────────────────────────────────────────────────
args     <- commandArgs(trailingOnly = TRUE)
CSV_PATH <- if (length(args) >= 1) args[1] else stop("Cần cung cấp CSV_PATH")
OUT_DIR  <- if (length(args) >= 2) args[2] else stop("Cần cung cấp OUT_DIR")

cat(rep("=", 70), "\n", sep = "")
cat("MNLFA_Supplementary.R — v1.0\n")
cat(sprintf("CSV: %s\nOUT: %s\n", CSV_PATH, OUT_DIR))
cat(rep("=", 70), "\n", sep = "")

# ── Helper: lưu Excel ─────────────────────────────────────────────────────────
save_xl <- function(obj, fname) {
  path <- file.path(OUT_DIR, fname)
  tryCatch({
    write_xlsx(as.data.frame(obj), path)
    cat(sprintf("  -> Saved: %s\n", fname))
  }, error = function(e)
    cat(sprintf("  [WARN] %s: %s\n", fname, conditionMessage(e))))
}

# ── Helper: significance stars ────────────────────────────────────────────────
sig <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
    ifelse(p < 0.05, "*", "ns"))))

# ── Đọc dữ liệu ──────────────────────────────────────────────────────────────
dat <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
N   <- nrow(dat)
cat(sprintf("\nDữ liệu: N = %d\n", N))

items_map <- list(
  AIP = paste0("AIP", 1:5),
  PA  = paste0("PA",  1:5),
  DR  = paste0("DR",  1:5),
  FSP = paste0("FSP", 1:5)
)
all_items <- unlist(items_map)
factors   <- c("AIP", "PA", "DR", "FSP")

# ── Tái tạo CFA để lấy factor scores (giống MNLFA_lavaan.R) ──────────────────
cat("\n[SETUP] Tái tạo CFA từ Stage 1...\n")
cfa_spec <- '
  AIP =~ AIP1 + AIP2 + AIP3 + AIP4 + AIP5
  PA  =~ PA1  + PA2  + PA3  + PA4  + PA5
  DR  =~ DR1  + DR2  + DR3  + DR4  + DR5
  FSP =~ FSP1 + FSP2 + FSP3 + FSP4 + FSP5
'
fit_cfa <- cfa(cfa_spec, data = dat, estimator = "MLR",
               std.lv = FALSE, missing = "listwise")
cat(sprintf("  CFA converged: %s\n", lavInspect(fit_cfa, "converged")))

# Bartlett factor scores
fs   <- as.data.frame(lavPredict(fit_cfa, type = "lv", method = "Bartlett"))
colnames(fs) <- factors

# Mean-center
fs_c <- as.data.frame(scale(fs, center = TRUE, scale = FALSE))
colnames(fs_c) <- paste0(factors, "_c")

# Tạo các terms phi tuyến và tương tác (giống file gốc)
fs_c$AIP_sq    <- fs_c$AIP_c^2
fs_c$PA_sq     <- fs_c$PA_c^2
fs_c$AIP_x_DR  <- fs_c$AIP_c  * fs_c$DR_c
fs_c$AIP_sq_DR <- fs_c$AIP_sq * fs_c$DR_c
fs_c$PA_x_DR   <- fs_c$PA_c   * fs_c$DR_c
fs_c$PA_sq_DR  <- fs_c$PA_sq  * fs_c$DR_c

# Tái tạo Model C
pred_cols <- c("AIP_c", "AIP_sq", "PA_c", "PA_sq", "DR_c",
               "AIP_x_DR", "AIP_sq_DR", "PA_x_DR", "PA_sq_DR")
modC <- lm(FSP_c ~ AIP_c + AIP_sq + PA_c + PA_sq + DR_c +
             AIP_x_DR + AIP_sq_DR + PA_x_DR + PA_sq_DR, data = fs_c)
fsp_sd <- sd(fs_c$FSP_c)

# DR quantiles
dr_q25 <- quantile(fs_c$DR_c, 0.25)
dr_q50 <- quantile(fs_c$DR_c, 0.50)
dr_q75 <- quantile(fs_c$DR_c, 0.75)
aip_sd <- sd(fs_c$AIP_c)
pa_sd  <- sd(fs_c$PA_c)

# Hệ số từ Model C
b1 <- coef(modC)["AIP_c"];    b2 <- coef(modC)["AIP_sq"]
b3 <- coef(modC)["PA_c"];     b4 <- coef(modC)["PA_sq"]
b5 <- coef(modC)["DR_c"]
b6 <- coef(modC)["AIP_x_DR"]; b7 <- coef(modC)["AIP_sq_DR"]
b8 <- coef(modC)["PA_x_DR"];  b9 <- coef(modC)["PA_sq_DR"]

cat("  Model C tái tạo thành công.\n")


# ============================================================================
# BLOCK 1: DESCRIPTIVE STATISTICS + CORRELATION MATRIX
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 1] DESCRIPTIVE STATISTICS + CORRELATION MATRIX\n")
cat(rep("=", 70), "\n", sep = "")

# 1a. Item-level descriptives
cat("\n  [1a] Item-level descriptives:\n")
desc_items <- data.frame(
  Item  = all_items,
  Mean  = sapply(all_items, function(x) mean(dat[[x]], na.rm = TRUE)),
  SD    = sapply(all_items, function(x) sd(dat[[x]],   na.rm = TRUE)),
  Min   = sapply(all_items, function(x) min(dat[[x]],  na.rm = TRUE)),
  Max   = sapply(all_items, function(x) max(dat[[x]],  na.rm = TRUE)),
  Skew  = sapply(all_items, function(x) {
    xv <- dat[[x]]; n <- length(xv)
    m3 <- mean((xv - mean(xv))^3); s <- sd(xv)
    m3 / s^3
  }),
  Kurt  = sapply(all_items, function(x) {
    xv <- dat[[x]]; n <- length(xv)
    m4 <- mean((xv - mean(xv))^4); s <- sd(xv)
    m4 / s^4 - 3   # excess kurtosis
  }),
  stringsAsFactors = FALSE
)
print(round(desc_items[, -1], 3))

# 1b. Factor score-level descriptives
cat("\n  [1b] Factor score descriptives (centered):\n")
desc_fs <- data.frame(
  Factor = factors,
  Mean   = sapply(factors, function(f) round(mean(fs[[f]]), 4)),
  SD     = sapply(factors, function(f) round(sd(fs[[f]]),   4)),
  Min    = sapply(factors, function(f) round(min(fs[[f]]),  4)),
  Max    = sapply(factors, function(f) round(max(fs[[f]]),  4)),
  stringsAsFactors = FALSE
)
print(desc_fs)

# 1c. Latent factor correlation matrix + AVE diagonal
cat("\n  [1c] Factor score correlation matrix:\n")
cor_fs <- cor(fs)
cat("  [Note: These are factor score correlations, approximating latent r]\n")
print(round(cor_fs, 3))

# 1d. Combined descriptive table (factor level) — paper-ready
#     Row: factor; Cols: Mean, SD, + correlations
cat("\n  [1d] Combined table (Mean, SD, r matrix) — paper-ready:\n")
sl_all  <- standardizedSolution(fit_cfa)
sl_all  <- sl_all[sl_all$op == "=~", ]
ave_vals <- sapply(factors, function(f)
  mean(sl_all[sl_all$lhs == f, "est.std"]^2))

#   Cronbach alpha
alpha_v <- sapply(factors, function(f) {
  cols <- items_map[[f]]; k <- length(cols)
  X    <- dat[, cols, drop = FALSE]
  iv   <- sum(apply(X, 2, var, na.rm = TRUE))
  tv   <- var(rowSums(X, na.rm = TRUE), na.rm = TRUE)
  k / (k - 1) * (1 - iv / tv)
})

desc_combined <- data.frame(
  Factor = factors,
  Mean   = round(sapply(factors, function(f) mean(rowMeans(dat[, items_map[[f]]]))), 2),
  SD     = round(sapply(factors, function(f) sd(rowMeans(dat[, items_map[[f]]]))),   2),
  Alpha  = round(alpha_v, 3),
  AVE    = round(ave_vals, 3),
  r_AIP  = round(cor_fs[, "AIP"], 3),
  r_PA   = round(cor_fs[, "PA"],  3),
  r_DR   = round(cor_fs[, "DR"],  3),
  r_FSP  = round(cor_fs[, "FSP"], 3),
  stringsAsFactors = FALSE
)
print(desc_combined)

# Lưu kết quả
save_xl(desc_items,    "Supp_Desc_Items.xlsx")
save_xl(desc_fs,       "Supp_Desc_FactorScores.xlsx")
save_xl(desc_combined, "Supp_Desc_Combined.xlsx")


# ============================================================================
# BLOCK 2: VIF DIAGNOSTICS
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 2] VIF DIAGNOSTICS (MULTICOLLINEARITY)\n")
cat(rep("=", 70), "\n", sep = "")

# Kiểm tra xem package 'car' có sẵn không; nếu không thì tính thủ công
vif_manual <- function(model) {
  X <- model.matrix(model)[, -1, drop = FALSE]  # bỏ intercept
  k <- ncol(X)
  vif_vals <- numeric(k)
  names(vif_vals) <- colnames(X)
  for (j in seq_len(k)) {
    m_j <- lm(X[, j] ~ X[, -j])
    vif_vals[j] <- 1 / (1 - summary(m_j)$r.squared)
  }
  vif_vals
}

cat("\n  [2a] VIF for Model C (all 9 predictors):\n")
vif_vals <- tryCatch({
  if (requireNamespace("car", quietly = TRUE)) {
    car::vif(modC)
  } else {
    cat("  [INFO] Package 'car' không có sẵn. Dùng VIF thủ công.\n")
    vif_manual(modC)
  }
}, error = function(e) {
  cat(sprintf("  [WARN] VIF error: %s. Dùng VIF thủ công.\n", conditionMessage(e)))
  vif_manual(modC)
})

vif_df <- data.frame(
  Term      = names(vif_vals),
  VIF       = round(as.numeric(vif_vals), 3),
  Tolerance = round(1 / as.numeric(vif_vals), 4),
  Flag      = ifelse(as.numeric(vif_vals) > 10, "[HIGH]",
                ifelse(as.numeric(vif_vals) > 5, "[MODERATE]", "[OK]")),
  stringsAsFactors = FALSE
)
cat(sprintf("  %-20s  %8s  %10s  %s\n", "Term", "VIF", "Tolerance", "Flag"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in seq_len(nrow(vif_df)))
  cat(sprintf("  %-20s  %8.3f  %10.4f  %s\n",
              vif_df$Term[i], vif_df$VIF[i], vif_df$Tolerance[i], vif_df$Flag[i]))

max_vif <- max(vif_df$VIF)
cat(sprintf("\n  Max VIF = %.3f\n", max_vif))
if (max_vif > 10) {
  cat("  [WARN] VIF > 10: Severe multicollinearity. Xem xét lại mô hình.\n")
} else if (max_vif > 5) {
  cat("  [NOTE] VIF 5-10: Moderate multicollinearity. Phổ biến với polynomial terms.\n")
  cat("         Trích dẫn: Aiken & West (1991) — mean-centering giúp giảm VIF cho\n")
  cat("         interaction terms nhưng không loại bỏ hoàn toàn với quadratic terms.\n")
} else {
  cat("  [OK] All VIF < 5: Multicollinearity acceptable.\n")
}

# [2b] VIF cho simplified models để so sánh
cat("\n  [2b] VIF cho Model A (linear only — baseline):\n")
modA <- lm(FSP_c ~ AIP_c + PA_c + DR_c, data = fs_c)
vif_A <- tryCatch(
  if (requireNamespace("car", quietly = TRUE)) car::vif(modA) else vif_manual(modA),
  error = function(e) vif_manual(modA)
)
for (nm in names(vif_A))
  cat(sprintf("  %-15s  VIF = %.3f\n", nm, vif_A[nm]))

save_xl(vif_df, "Supp_VIF_ModelC.xlsx")


# ============================================================================
# BLOCK 3: COMMON METHOD BIAS — HARMAN'S SINGLE FACTOR TEST
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 3] COMMON METHOD BIAS — HARMAN'S SINGLE FACTOR TEST\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n  Phương pháp: EFA 1 nhân tố trên toàn bộ 20 items\n")
cat("  Tiêu chí: Nếu 1 nhân tố giải thích > 50% tổng variance -> CMB là vấn đề\n\n")

# EFA với 1 nhân tố
X_efa <- as.matrix(dat[, all_items])

# Dùng principal component (không phải common factor) — đúng theo Harman (1967)
pca_result <- tryCatch({
  prcomp(X_efa, scale. = TRUE)
}, error = function(e) {
  cat(sprintf("  [WARN] PCA error: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(pca_result)) {
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cat(sprintf("  PC1 variance explained: %.4f (%.2f%%)\n",
              var_explained[1], var_explained[1] * 100))
  cat(sprintf("  PC2 variance explained: %.4f (%.2f%%)\n",
              var_explained[2], var_explained[2] * 100))
  cat(sprintf("  PC3 variance explained: %.4f (%.2f%%)\n",
              var_explained[3], var_explained[3] * 100))

  pct_pc1 <- var_explained[1] * 100
  if (pct_pc1 > 50) {
    cat(sprintf("\n  [WARN] PC1 = %.1f%% > 50%% -> CMB có thể là vấn đề.\n", pct_pc1))
    cat("         Cần báo cáo và thảo luận trong Limitations.\n")
  } else {
    cat(sprintf("\n  [OK] PC1 = %.1f%% < 50%% -> Không có bằng chứng CMB nghiêm trọng.\n",
                pct_pc1))
    cat("         Có thể báo cáo trong bản thảo như sau:\n")
    cat(sprintf("         'Harman's single factor test indicated that the first\n"))
    cat(sprintf("          unrotated factor accounted for %.1f%% of the total\n", pct_pc1))
    cat(sprintf("          variance (threshold: 50%%), suggesting that common\n"))
    cat(sprintf("          method bias does not pose a serious threat to the\n"))
    cat(sprintf("          validity of this study (Podsakoff et al., 2003).'\n"))
  }

  # Cumulative variance cho 20 components
  cum_var <- cumsum(var_explained)
  harman_df <- data.frame(
    PC            = seq_along(var_explained),
    Eigenvalue    = round(pca_result$sdev^2, 4),
    Var_Explained = round(var_explained * 100, 3),
    Cum_Var       = round(cum_var * 100, 3),
    stringsAsFactors = FALSE
  )
  cat("\n  Eigenvalue summary (top 10 PCs):\n")
  print(harman_df[1:10, ])

  save_xl(harman_df, "Supp_Harman_CMB.xlsx")
}

cat("\n  [NOTE] CFA-based CMB test đã được chạy trong MNLFA_lavaan.R.\n")
cat("  Kết quả xem tại: Step1b_CMB_CFA_Test.xlsx\n")
cat("  Báo cáo cả hai kết quả (Harman + CFA-CMF) trong Section 4.1.\n")

# ============================================================================
# BLOCK 4: SIMPLE SLOPES TABLE — PUBLICATION READY
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 4] SIMPLE SLOPES TABLE — PUBLICATION READY\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n  Simple slope = marginal effect of predictor at given X and DR level\n")
cat("  For quadratic: slope at x0 = (b_lin + b_int_lin*DR) + 2*(b_quad + b_int_quad*DR)*x0\n\n")

# Helper: HC3 variance-covariance — vectorised (safe for large N)
# Fix v2: hatvalues() thay cho diag(X %*% solve() %*% t(X)) để tránh
# build ma trận N×N gây NaN/NA khi N lớn (e.g. N=537).
# W là vector thay vì diagonal matrix, dùng broadcasting (t(X) * w).
hc3_vcov <- function(m) {
  X  <- model.matrix(m)
  e  <- residuals(m)
  h  <- hatvalues(m)                   # hat-diagonal, không cần ma trận N×N
  w  <- (e / (1 - h))^2               # HC3 weights — vector
  Xi <- solve(crossprod(X))
  Xi %*% (t(X) * w) %*% X %*% Xi     # broadcasting thay cho diag(W)
}
vc <- hc3_vcov(modC)
co <- coef(modC)

# Slope cho AIP tại x0 với DR = dr_val:
#   d(FSP)/d(AIP) = (b1 + b6*DR) + 2*(b2 + b7*DR)*AIP
# Gradient vector đối với tất cả params:
#   [0, x0, 0, 0, 0, dr_val, x0^2? không — đây là đạo hàm, phức tạp hơn]
# Dùng delta method: slope_aip(x0,dr) = b1 + b6*dr + 2*b2*x0 + 2*b7*x0*dr
# Gradient: d/d_params = [0, 1, 2*x0, 0, 0, 0, dr, 2*x0*dr, 0, 0] (theo thứ tự coef)
# Thứ tự coef: (Intercept) AIP_c AIP_sq PA_c PA_sq DR_c AIP_x_DR AIP_sq_DR PA_x_DR PA_sq_DR

slope_se_aip <- function(x0, dr) {
  # slope = b1 + b6*dr + 2*b2*x0 + 2*b7*x0*dr
  # grad: position của từng coef trong vector co
  g <- setNames(rep(0, length(co)), names(co))
  g["AIP_c"]     <- 1
  g["AIP_sq"]    <- 2 * x0
  g["AIP_x_DR"]  <- dr
  g["AIP_sq_DR"] <- 2 * x0 * dr
  slope_val <- co["AIP_c"] + co["AIP_x_DR"] * dr +
               2 * co["AIP_sq"] * x0 + 2 * co["AIP_sq_DR"] * x0 * dr
  se_val    <- sqrt(as.numeric(t(g) %*% vc %*% g))
  t_val     <- slope_val / se_val
  p_val     <- 2 * pt(-abs(t_val), df = N - length(co))
  c(slope = slope_val, se = se_val, t = t_val, p = p_val)
}

slope_se_pa <- function(x0, dr) {
  g <- setNames(rep(0, length(co)), names(co))
  g["PA_c"]     <- 1
  g["PA_sq"]    <- 2 * x0
  g["PA_x_DR"]  <- dr
  g["PA_sq_DR"] <- 2 * x0 * dr
  slope_val <- co["PA_c"] + co["PA_x_DR"] * dr +
               2 * co["PA_sq"] * x0 + 2 * co["PA_sq_DR"] * x0 * dr
  se_val    <- sqrt(as.numeric(t(g) %*% vc %*% g))
  t_val     <- slope_val / se_val
  p_val     <- 2 * pt(-abs(t_val), df = N - length(co))
  c(slope = slope_val, se = se_val, t = t_val, p = p_val)
}

# DR raw scale grand mean (để báo cáo)
dr_gm_raw  <- mean(rowMeans(dat[, paste0("DR", 1:5)]))
aip_gm_raw <- mean(rowMeans(dat[, paste0("AIP", 1:5)]))
pa_gm_raw  <- mean(rowMeans(dat[, paste0("PA",  1:5)]))

dr_levels  <- list(
  list(lbl = "Low DR (Q25)",  val = dr_q25, raw = dr_q25  + dr_gm_raw),
  list(lbl = "Med DR (Q50)",  val = dr_q50, raw = dr_q50  + dr_gm_raw),
  list(lbl = "High DR (Q75)", val = dr_q75, raw = dr_q75  + dr_gm_raw)
)
aip_levels <- list(
  list(lbl = "Low AIP (−1SD)",  val = -aip_sd, raw = -aip_sd + aip_gm_raw),
  list(lbl = "High AIP (+1SD)", val = +aip_sd, raw = +aip_sd + aip_gm_raw)
)
pa_levels  <- list(
  list(lbl = "Low PA (−1SD)",   val = -pa_sd,  raw = -pa_sd  + pa_gm_raw),
  list(lbl = "High PA (+1SD)",  val = +pa_sd,  raw = +pa_sd  + pa_gm_raw)
)

ss_rows <- list()

cat(sprintf("  %-22s  %-16s  %8s  %6s  %6s  %s\n",
            "AIP Level", "DR Level", "Slope", "SE", "p", "sig"))
cat(paste(rep("-", 68), collapse = ""), "\n")
for (xl in aip_levels) {
  for (dl in dr_levels) {
    res <- slope_se_aip(xl$val, dl$val)
    cat(sprintf("  %-22s  %-16s  %+8.4f  %6.4f  %6.4f  %s\n",
                xl$lbl, dl$lbl, res["slope"], res["se"], res["p"], sig(res["p"])))
    ss_rows[[length(ss_rows) + 1]] <- data.frame(
      Predictor  = "AIP",
      X_Level    = xl$lbl,
      DR_Level   = dl$lbl,
      DR_raw     = round(dl$raw, 3),
      Slope      = round(res["slope"], 4),
      SE         = round(res["se"],    4),
      t_value    = round(res["t"],     3),
      p_value    = round(res["p"],     4),
      sig        = sig(res["p"]),
      stringsAsFactors = FALSE
    )
  }
}

cat(sprintf("\n  %-22s  %-16s  %8s  %6s  %6s  %s\n",
            "PA Level", "DR Level", "Slope", "SE", "p", "sig"))
cat(paste(rep("-", 68), collapse = ""), "\n")
for (xl in pa_levels) {
  for (dl in dr_levels) {
    res <- slope_se_pa(xl$val, dl$val)
    cat(sprintf("  %-22s  %-16s  %+8.4f  %6.4f  %6.4f  %s\n",
                xl$lbl, dl$lbl, res["slope"], res["se"], res["p"], sig(res["p"])))
    ss_rows[[length(ss_rows) + 1]] <- data.frame(
      Predictor  = "PA",
      X_Level    = xl$lbl,
      DR_Level   = dl$lbl,
      DR_raw     = round(dl$raw, 3),
      Slope      = round(res["slope"], 4),
      SE         = round(res["se"],    4),
      t_value    = round(res["t"],     3),
      p_value    = round(res["p"],     4),
      sig        = sig(res["p"]),
      stringsAsFactors = FALSE
    )
  }
}

ss_df <- do.call(rbind, ss_rows)
save_xl(ss_df, "Supp_Simple_Slopes_Full.xlsx")


# ============================================================================
# BLOCK 5: FORNELL-LARCKER MATRIX (DISCRIMINANT VALIDITY)
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 5] FORNELL-LARCKER MATRIX\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n  [Diagonal = sqrt(AVE); Off-diagonal = latent factor correlations]\n\n")

cor_lv   <- lavInspect(fit_cfa, "cor.lv")
sl_std   <- standardizedSolution(fit_cfa)
sl_std   <- sl_std[sl_std$op == "=~", ]

# Tính AVE theo thứ tự của vector factors (dùng index, không dùng tên)
ave_v    <- sapply(factors, function(f) {
  loadings <- sl_std[sl_std$lhs == f, "est.std"]
  if (length(loadings) == 0) return(NA_real_)
  mean(loadings^2)
})
ave_v    <- as.numeric(ave_v)          # bỏ names để tránh mismatch
sqrt_ave <- sqrt(ave_v)

# Sắp xếp cor_lv theo đúng thứ tự factors
lv_names <- rownames(cor_lv)
idx      <- match(factors, lv_names)   # vị trí của mỗi factor trong cor_lv
if (any(is.na(idx))) {
  cat(sprintf("  [WARN] Latent variable names in cor.lv: %s\n",
              paste(lv_names, collapse = ", ")))
  cat(sprintf("  [WARN] Expected factors: %s\n", paste(factors, collapse = ", ")))
  idx <- seq_along(factors)            # fallback: dùng thứ tự hiện có
}
cor_lv_ord <- cor_lv[idx, idx]        # reorder theo factors

fl_mat <- cor_lv_ord
for (i in seq_along(factors)) fl_mat[i, i] <- sqrt_ave[i]  # gán diagonal theo index

cat(sprintf("  %-8s", ""))
for (f in factors) cat(sprintf("  %8s", f))
cat("\n")
cat(paste(rep("-", 8 + 9 * length(factors)), collapse = ""), "\n")
for (i in seq_along(factors)) {
  cat(sprintf("  %-8s", factors[i]))
  for (j in seq_along(factors)) {
    val <- fl_mat[i, j]
    if (i == j)
      cat(sprintf("  [%6.3f]", val))
    else
      cat(sprintf("  %8.3f",   val))
  }
  cat("\n")
}

# Chuẩn Fornell-Larcker: sqrt(AVE_i) > |r_ij| cho mọi j ≠ i
# Dùng index số thay vì tên để tránh NA từ name mismatch
cat("\n  Fornell-Larcker criterion (sqrt(AVE) > max|r|):\n")
fl_pass <- logical(length(factors))
names(fl_pass) <- factors
for (fi in seq_along(factors)) {
  f        <- factors[fi]
  sq_ave_f <- as.numeric(sqrt_ave[fi])
  others   <- seq_along(factors)[-fi]
  max_r    <- max(abs(as.numeric(cor_lv_ord[fi, others])))
  pass     <- isTRUE(sq_ave_f > max_r)
  fl_pass[fi] <- pass
  cat(sprintf("  %s: sqrt(AVE)=%.3f > max|r|=%.3f -> %s\n",
              f, sq_ave_f, max_r, if (pass) "PASS" else "FAIL"))
}
cat(sprintf("  Overall: %s\n",
            if (all(fl_pass, na.rm = TRUE)) "ALL PASS" else "SOME FAIL"))

fl_df <- as.data.frame(round(fl_mat, 3))
colnames(fl_df) <- factors
fl_df$Factor <- factors
fl_df <- fl_df[, c("Factor", factors)]
save_xl(fl_df, "Supp_FornellLarcker.xlsx")


# ============================================================================
# BLOCK 6: ROBUSTNESS COMPARISON TABLE (EXTENDED)
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("[BLOCK 6] ROBUSTNESS COMPARISON TABLE — EXTENDED\n")
cat(rep("=", 70), "\n", sep = "")

cat("\n  Đọc kết quả outlier và MM từ các file Excel đã tạo...\n")

# Đọc file robustness đã có từ file gốc
rob_outlier <- tryCatch(
  as.data.frame(read_xlsx(file.path(OUT_DIR, "Step10a_Robustness_Outlier.xlsx"))),
  error = function(e) { cat(sprintf("  [WARN] %s\n", conditionMessage(e))); NULL }
)
rob_mm <- tryCatch(
  as.data.frame(read_xlsx(file.path(OUT_DIR, "Step10d_Robustness_MM.xlsx"))),
  error = function(e) { cat(sprintf("  [WARN] %s\n", conditionMessage(e))); NULL }
)

# Helper: lấy b và p theo tên term
get_bp <- function(df, term, b_col = "b", p_col = "p") {
  if (is.null(df)) return(c(b = NA, p = NA))
  row <- df[df[[1]] == term | df$Term == term, , drop = FALSE]
  if (nrow(row) == 0) return(c(b = NA, p = NA))
  b_val <- as.numeric(row[[b_col]][1])
  # Cột p có thể tên khác nhau
  p_col_found <- intersect(c(p_col, "p", "p_MM", "p_value"), colnames(df))[1]
  p_val <- if (!is.na(p_col_found)) as.numeric(row[[p_col_found]][1]) else NA
  c(b = b_val, p = p_val)
}

# Tái chạy HC3 cho Model C (main) — vectorised, safe for large N
# Fix v2: same approach as Block 4
hc3_vcov_fn <- function(m) {
  X  <- model.matrix(m)
  e  <- residuals(m)
  h  <- hatvalues(m)
  w  <- (e / (1 - h))^2
  Xi <- solve(crossprod(X))
  Xi %*% (t(X) * w) %*% X %*% Xi
}
vc_main <- hc3_vcov_fn(modC)
se_main <- sqrt(diag(vc_main))
p_main  <- 2 * pt(-abs(coef(modC) / se_main), df = N - length(coef(modC)))

key_terms <- c("AIP_sq", "PA_sq", "AIP_sq_DR", "PA_sq_DR", "DR_c")
term_labels <- c(
  AIP_sq     = "AIP² (H1)",
  PA_sq      = "PA² (H2)",
  AIP_sq_DR  = "AIP²×DR (H3)",
  PA_sq_DR   = "PA²×DR (H4)",
  DR_c       = "DR direct (H5)"
)

cat(sprintf("\n  %-18s  %10s  %10s  %10s  %10s  %10s  %10s\n",
            "Term",
            "Main_b", "Main_p",
            "Outlier_b", "Outlier_p",
            "MM_b", "MM_p"))
cat(paste(rep("-", 92), collapse = ""), "\n")

rob_rows <- list()
for (tm in key_terms) {
  b_main <- coef(modC)[tm]
  p_m    <- p_main[tm]
  bp_out <- get_bp(rob_outlier, tm, b_col = "b", p_col = "p")
  bp_mm  <- get_bp(rob_mm, tm, b_col = "b_MM", p_col = "p_MM")

  cat(sprintf("  %-18s  %+10.4f  %10.4f  %+10.4f  %10.4f  %+10.4f  %10.4f\n",
              term_labels[tm],
              b_main, p_m,
              bp_out["b"], bp_out["p"],
              bp_mm["b"],  bp_mm["p"]))

  rob_rows[[tm]] <- data.frame(
    Term      = term_labels[tm],
    Main_b    = round(b_main,       4),
    Main_p    = round(p_m,          4),
    Main_sig  = sig(p_m),
    Outlier_b = round(bp_out["b"],  4),
    Outlier_p = round(bp_out["p"],  4),
    Outlier_sig = sig(bp_out["p"]),
    MM_b      = round(bp_mm["b"],   4),
    MM_p      = round(bp_mm["p"],   4),
    MM_sig    = sig(bp_mm["p"]),
    stringsAsFactors = FALSE
  )
}
cat(sprintf("  Main model N = %d\n", N))

rob_ext_df <- do.call(rbind, rob_rows)
rownames(rob_ext_df) <- NULL

# Thêm R2
cat(sprintf("\n  R² (Main): %.4f | Adj-R²: %.4f\n",
            summary(modC)$r.squared, summary(modC)$adj.r.squared))

save_xl(rob_ext_df, "Supp_Robustness_Extended.xlsx")


# ============================================================================
# TỔNG KẾT
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("SUPPLEMENTARY ANALYSES HOÀN THÀNH\n")
cat(rep("=", 70), "\n", sep = "")

expected_outputs <- c(
  "Supp_Desc_Items.xlsx",
  "Supp_Desc_FactorScores.xlsx",
  "Supp_Desc_Combined.xlsx",
  "Supp_VIF_ModelC.xlsx",
  "Supp_Harman_CMB.xlsx",
  "Supp_Simple_Slopes_Full.xlsx",
  "Supp_FornellLarcker.xlsx",
  "Supp_Robustness_Extended.xlsx"
)

cat("\nFiles được tạo:\n")
for (f in expected_outputs) {
  fp <- file.path(OUT_DIR, f)
  cat(sprintf("  %s  %s\n", if (file.exists(fp)) "[OK]" else "[NOT FOUND]", f))
}

cat("\nGhi chú cho bài báo:\n")
cat("  Block 1 -> Table 1 (Descriptive Statistics & Correlations)\n")
cat("  Block 2 -> Footnote trong Section 3.4 hoặc Appendix\n")
cat("  Block 3 -> Paragraph trong Section 3.4 (CMB)\n")
cat("  Block 4 -> Table bổ sung hoặc Appendix (Simple Slopes)\n")
cat("  Block 5 -> Table 2 (Measurement Model) — thêm vào FL matrix\n")
cat("  Block 6 -> Table Robustness Summary\n")
cat(rep("=", 70), "\n", sep = "")