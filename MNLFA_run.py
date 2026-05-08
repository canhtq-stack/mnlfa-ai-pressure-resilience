#!/usr/bin/env python3
"""
MNLFA_run.py — v5.1  PYTHON WRAPPER
=====================================
Buoc 1: Doc Excel -> kiem tra data -> export CSV
Buoc 2: Goi MNLFA_lavaan.R qua Rscript (phan tich chinh)
Buoc 3: Goi MNLFA_Supplementary.R qua Rscript (phan tich bo sung)
Buoc 4: In tom tat ket qua + kiem tra output files

Su dung:
    python MNLFA_run.py

Truoc khi chay:
    pip install pandas openpyxl scipy numpy
    R: install.packages(c('lavaan','semTools','writexl','MASS','readxl'))
    R: install.packages('car')  # tuy chon — de tinh VIF; neu khong co thi dung thu cong

Thay doi tu v5.0:
    - Them Buoc 3: chay MNLFA_Supplementary.R sau khi main script hoan thanh
    - Them kiem tra output file supplementary
    - Them huong dan su dung ket qua vao ban thao
"""

import os, sys, subprocess
import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════════
# CẤU HÌNH — CHỈNH ĐƯỜNG DẪN Ở ĐÂY
# ═══════════════════════════════════════════════════════════════════════════
DATA_PATH   = r"D:\CANH\Nghiên cứu khoa học\VIẾT BÁO - CHẠY SỐ LIỆU\Cân bằng giữa Đổi mới và Sức khỏe Tinh thần\Data.xlsx"
OUT_DIR     = os.path.dirname(DATA_PATH)       # Lưu kết quả cùng thư mục data
N_BOOTSTRAP = 2000                             # Giảm xuống 500 nếu chậm
SEED        = 42

# Đường dẫn Rscript — R 4.4.1
RSCRIPT_PATH = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"

# R script chính (cùng thư mục với file này hoặc chỉ định đường dẫn đầy đủ)
R_SCRIPT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "MNLFA_lavaan.R")

# R script supplementary (phân tích bổ sung cho bài báo)
R_SUPP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "MNLFA_Supplementary.R")

# Bật/tắt chạy supplementary (True = chạy sau main script)
RUN_SUPPLEMENTARY = True
# ═══════════════════════════════════════════════════════════════════════════

CSV_PATH = os.path.join(OUT_DIR, "_mnlfa_data.csv")
LOG_PATH = os.path.join(OUT_DIR, "_mnlfa_lavaan.log")

print("=" * 78)
print("MNLFA v5.0 — PYTHON WRAPPER + LAVAAN R BACKEND")
print("=" * 78)


# ─────────────────────────────────────────────────────────────────────────────
def cronbach_alpha(df_items):
    k = df_items.shape[1]
    iv = df_items.var(axis=0, ddof=1).sum()
    tv = df_items.sum(axis=1).var(ddof=1)
    return k / (k - 1) * (1 - iv / tv)


# ═══════════════════════════════════════════════════════════════════════════
# BƯỚC 1: ĐỌC & KIỂM TRA DỮ LIỆU
# ═══════════════════════════════════════════════════════════════════════════
print("\n[STEP 1] ĐỌC VÀ KIỂM TRA DỮ LIỆU")

df = pd.read_excel(DATA_PATH)
print(f"    Raw: {df.shape[0]} obs × {df.shape[1]} vars")

items = {
    'AIP': [f'AIP{i}' for i in range(1, 6)],
    'PA':  [f'PA{i}'  for i in range(1, 6)],
    'DR':  [f'DR{i}'  for i in range(1, 6)],
    'FSP': [f'FSP{i}' for i in range(1, 6)],
}
all_items = [it for lst in items.values() for it in lst]

# Kiểm tra items có tồn tại không
missing_cols = [c for c in all_items if c not in df.columns]
if missing_cols:
    print(f"    [ERROR] Thiếu cột: {missing_cols}")
    sys.exit(1)

# Loại bỏ missing
df = df.dropna(subset=all_items).reset_index(drop=True)
print(f"    Sau listwise deletion: {len(df)} obs")

# Loại bỏ straight-liners
for name, cols in items.items():
    mask = df[cols].std(axis=1) < 1e-6
    if mask.any():
        print(f"    Removed {mask.sum()} straight-liners ({name})")
        df = df[~mask].reset_index(drop=True)

N = len(df)
print(f"    Final N = {N}")

# Descriptive statistics
print("\n    Descriptive Statistics:")
desc = df[all_items].agg(['mean', 'std', 'min', 'max']).T
desc['skew'] = df[all_items].apply(stats.skew)
desc['kurt'] = df[all_items].apply(lambda x: stats.kurtosis(x, fisher=True))
print(desc.round(3).to_string())

# Cronbach alpha
print("\n    Cronbach Alpha:")
for name, cols in items.items():
    a = cronbach_alpha(df[cols])
    print(f"    {name}: α = {a:.3f}")

# Kiểm tra phân phối (Mardia multivariate skewness — simplified)
X = df[all_items].values
Xc = X - X.mean(axis=0)
try:
    from numpy.linalg import inv, LinAlgError
    S_inv = inv(np.cov(X.T))
    D2 = np.einsum('ij,jk,ik->i', Xc, S_inv, Xc)
    mardia_b1 = float((D2**3).mean() / 6)
    df_mardia = len(all_items) * (len(all_items)+1) * (len(all_items)+2) / 6
    p_mardia = 1 - stats.chi2.cdf(mardia_b1, df=df_mardia)
    print(f"\n    Mardia multivariate skewness: b1p = {mardia_b1:.3f} (p = {p_mardia:.4f})")
    if p_mardia < 0.05:
        print("    → Non-normality detected: MLR estimator được dùng trong R (robust SE)")
    else:
        print("    → Multivariate normality OK")
except Exception:
    print("    → Mardia test skipped")

# ═══════════════════════════════════════════════════════════════════════════
# BƯỚC 2: EXPORT CSV CHO R
# ═══════════════════════════════════════════════════════════════════════════
print(f"\n[STEP 2] EXPORT CSV → {CSV_PATH}")
df[all_items].to_csv(CSV_PATH, index=False)
print(f"    Exported {N} rows × {len(all_items)} item columns")

# ═══════════════════════════════════════════════════════════════════════════
# BƯỚC 3: CHẠY R SCRIPT
# ═══════════════════════════════════════════════════════════════════════════
print(f"\n[STEP 3] CHẠY MNLFA_lavaan.R (N_BOOT={N_BOOTSTRAP}, seed={SEED})")
print(f"    R script: {R_SCRIPT_PATH}")

# Kiểm tra R script tồn tại
if not os.path.exists(R_SCRIPT_PATH):
    print(f"    [ERROR] R script không tìm thấy: {R_SCRIPT_PATH}")
    print("    Đảm bảo MNLFA_lavaan.R nằm cùng thư mục với file này.")
    sys.exit(1)

cmd = [
    RSCRIPT_PATH, '--vanilla',
    R_SCRIPT_PATH,
    CSV_PATH.replace('\\', '/'),
    OUT_DIR.replace('\\', '/'),
    str(N_BOOTSTRAP),
    str(SEED)
]

print(f"    Command: {' '.join(cmd[:2])} ... (4 args)")
print("    Đang chạy... (bootstrap có thể mất 5-15 phút)\n")
print("=" * 78)

try:
    proc = subprocess.run(
        cmd,
        capture_output=True, text=True, encoding='utf-8',
        timeout=3600   # 60 min max
    )

    # In stdout ngay (R output)
    stdout = proc.stdout
    if stdout:
        print(stdout)

    # Ghi log đầy đủ
    with open(LOG_PATH, 'w', encoding='utf-8') as f:
        f.write("=== COMMAND ===\n")
        f.write(' '.join(cmd) + '\n\n')
        f.write("=== STDOUT ===\n")
        f.write(stdout)
        f.write("\n=== STDERR ===\n")
        f.write(proc.stderr)

    if proc.returncode != 0:
        print("\n" + "=" * 78)
        print("[WARN] R script kết thúc với lỗi (returncode ≠ 0)")
        print("Stderr (last 2000 chars):")
        print(proc.stderr[-2000:])
    else:
        print("\n" + "=" * 78)
        print("[OK] R script hoàn thành thành công")
        print(f"    Full log: {LOG_PATH}")

        print("\n    Output files (cùng thư mục dữ liệu):")
        expected_files = [
            "Step1_CFA_Results.xlsx",
            "Step0_DIF_Test.xlsx",
            "ModelA_Linear_Results.xlsx",
            "ModelB_Quadratic_Results.xlsx",
            "ModelC_Full_Results.xlsx",
            "Step6_Bootstrap_Results.xlsx",
            "Step8_Inflection_Points.xlsx",
            "Step9_Simple_Slopes.xlsx",
            "Table3_Fit_Indices.xlsx",
            "Table4_Structural_Coefficients.xlsx",
            "Table_Inflection_Raw_Scale.xlsx",
            "Table_Robustness_Summary.xlsx",
            "Fig5_Simple_Slopes_Publication.pdf",
            "Fig1_GAM_Nonlinear_Shapes.pdf",
            "Fig4_Overlay_Interactions.pdf",
            "Step10a_Robustness_Outlier.xlsx",
            "Step10d_Robustness_MM.xlsx",
            "Table_Power_Analysis.xlsx",
        ]
        for fname in expected_files:
            fpath = os.path.join(OUT_DIR, fname)
            status = "✓" if os.path.exists(fpath) else "✗ (not created)"
            print(f"    {status}  {fname}")

except FileNotFoundError:
    print(f"\n[ERROR] Không tìm thấy '{RSCRIPT_PATH}'")
    print("\nCách khắc phục trên Windows:")
    print("    1. Đảm bảo R đã được cài (https://cran.r-project.org/)")
    print("    2. Chỉnh biến RSCRIPT_PATH trong file này:")
    print(r'       RSCRIPT_PATH = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"')

except subprocess.TimeoutExpired:
    print("\n[ERROR] R script timeout (>60 min)")
    print("    Thử giảm N_BOOTSTRAP xuống 500 và chạy lại")

# ═══════════════════════════════════════════════════════════════════════════
# BƯỚC 4: CHẠY MNLFA_Supplementary.R (phân tích bổ sung)
# ═══════════════════════════════════════════════════════════════════════════
if RUN_SUPPLEMENTARY:
    print(f"\n[STEP 4] CHẠY MNLFA_Supplementary.R")
    print(f"    R script: {R_SUPP_PATH}")

    if not os.path.exists(R_SUPP_PATH):
        print(f"    [WARN] MNLFA_Supplementary.R không tìm thấy: {R_SUPP_PATH}")
        print("    Bỏ qua bước này. Đảm bảo file nằm cùng thư mục với MNLFA_run.py")
    else:
        cmd_supp = [
            RSCRIPT_PATH, '--vanilla',
            R_SUPP_PATH,
            CSV_PATH.replace('\\', '/'),
            OUT_DIR.replace('\\', '/')
        ]
        print("    Đang chạy MNLFA_Supplementary.R...\n")
        print("=" * 78)

        try:
            proc_supp = subprocess.run(
                cmd_supp,
                capture_output=True, text=True, encoding='utf-8',
                timeout=600   # 10 min max cho supplementary
            )

            if proc_supp.stdout:
                print(proc_supp.stdout)

            # Ghi log supplementary
            supp_log = os.path.join(OUT_DIR, "_mnlfa_supplementary.log")
            with open(supp_log, 'w', encoding='utf-8') as f:
                f.write("=== SUPPLEMENTARY ANALYSES LOG ===\n")
                f.write(' '.join(cmd_supp) + '\n\n')
                f.write("=== STDOUT ===\n")
                f.write(proc_supp.stdout)
                f.write("\n=== STDERR ===\n")
                f.write(proc_supp.stderr)

            if proc_supp.returncode != 0:
                print("\n" + "=" * 78)
                print("[WARN] Supplementary script kết thúc với lỗi")
                print("Stderr (last 1000 chars):")
                print(proc_supp.stderr[-1000:])
            else:
                print("\n" + "=" * 78)
                print("[OK] MNLFA_Supplementary.R hoàn thành")
                print(f"    Log: {supp_log}")

                print("\n    Supplementary output files:")
                supp_files = [
                    "Supp_Desc_Items.xlsx",
                    "Supp_Desc_FactorScores.xlsx",
                    "Supp_Desc_Combined.xlsx",
                    "Supp_VIF_ModelC.xlsx",
                    "Supp_Harman_CMB.xlsx",
                    "Supp_Simple_Slopes_Full.xlsx",
                    "Supp_FornellLarcker.xlsx",
                    "Supp_Robustness_Extended.xlsx",
                ]
                for fname in supp_files:
                    fpath = os.path.join(OUT_DIR, fname)
                    status = "✓" if os.path.exists(fpath) else "✗ (not created)"
                    print(f"    {status}  {fname}")

                print("\n    Hướng dẫn sử dụng kết quả vào bản thảo:")
                print("    Supp_Desc_Combined.xlsx   -> Table 1 (Descriptive Stats & Correlations)")
                print("    Supp_VIF_ModelC.xlsx      -> Footnote Section 3.4 (VIF values)")
                print("    Supp_Harman_CMB.xlsx      -> Paragraph Section 3.4 (CMB)")
                print("    Supp_Simple_Slopes_Full.xlsx -> Appendix Table / Supplementary")
                print("    Supp_FornellLarcker.xlsx  -> Table 2 (Measurement Model)")
                print("    Supp_Robustness_Extended.xlsx -> Table 5 (Robustness)")

        except subprocess.TimeoutExpired:
            print("\n[ERROR] Supplementary script timeout (>10 min)")
        except Exception as ex:
            print(f"\n[ERROR] Supplementary: {ex}")

print("\n" + "=" * 78)
print("PIPELINE v5.1 HOÀN THÀNH")
print("=" * 78)
