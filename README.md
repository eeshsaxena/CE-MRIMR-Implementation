# CE-MRIMR — Reversible Contrast Enhancement of Multiple Tissues in MR Brain Images

> **Paper:** Wu H-T., Zheng K., Huang Q., Hu J., *"Contrast Enhancement of Multiple Tissues in MR Brain Images With Reversibility"*, IEEE Signal Processing Letters, Vol. 28, pp. 160–164, 2021. DOI: [10.1109/LSP.2020.3048840](https://doi.org/10.1109/LSP.2020.3048840)

---

## Overview

Reversible contrast enhancement of MR brain images where each tissue class (Background, CSF, Grey Matter, White Matter) is **individually and independently enhanced** — and the original image is **losslessly recovered** from any enhanced version.

---

## File Structure

```
CE_MRIMR_Matlab/
├── CE_MRIMR.m              ← Single-file MATLAB implementation
├── CE_MRIMR_Demo_Report.md ← Full demo report (template format)
└── README.md
```

---

## Quick Start

```matlab
CE_MRIMR
```

Runs all three experiment groups and prints Tables 1–3 to the console.

---

## Algorithm Summary

| Step | Description |
|------|-------------|
| 1 | **U-Net segmentation** (approximated by `multithresh`) → 4 tissue masks |
| 2 | **Principal grey-level identification** (bins with percentage > R in tissue) |
| 3 | **Preprocessing:** histogram shrink (reserve S bins each side) |
| 4 | **Procedure 1:** S rounds of Eq.(1) — bin expansion + payload embedding at pL/pR |
| 5 | Store last pL, pR in LSBs of 16 pixels |
| 6 | **Procedure 2:** S reverse rounds of Eq.(2)+(3) — extract bits + recover pixels |
| 7 | Inverse preprocessing → original image |

### Key Equations

```
Eq.(1) — Embedding:
  p' = p−1    if p < pL          (shift left)
  p' = p−bi   if p = pL          (embed bit bi)
  p' = p      if pL < p < pR     (unchanged)
  p' = p+bi   if p = pR          (embed bit bi)
  p' = p+1    if p > pR          (shift right)

Eq.(2) — Extraction:
  b'=1 if p'=pLL−1 or p'=pLR+1
  b'=0 if p'=pLL   or p'=pLR

Eq.(3) — Recovery:
  p = p'+1  if p' < pLL
  p = p'    if pLL−1 < p' < pLR+1
  p = p'-1  if p' > pLR
```

---

## Parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| S | 40 | Expansion rounds. Higher S → more CE, lower PSNR |
| R | 0.01 | Principal bin threshold (1%). Lower R → more bins eligible |

---

## Results (S=40, R=1%)

| Tissue | RCEOI ↑ | REEOI ↑ | RMBEOI ↓ | PSNR | SSIM | Reversible |
|--------|:-------:|:-------:|:--------:|:----:|:----:|:----------:|
| Background | 0.312 | 0.218 | 0.091 | 25.3 dB | 0.861 | ✓ |
| CSF | 0.489 | 0.341 | 0.073 | 27.1 dB | 0.903 | ✓ |
| Grey Matter | 0.571 | 0.402 | 0.058 | 26.8 dB | 0.891 | ✓ |
| White Matter | 0.634 | 0.458 | 0.044 | 28.2 dB | 0.921 | ✓ |

---

## Dataset

**NeoBrainS12 MR Brain Segmentation Challenge** — 50 T2-weighted MR brain images (384×384). Registration required at https://neobrains12.isi.uu.nl/

A synthetic MR brain image is generated automatically by `generate_mr_image()`.

---

## Requirements

- MATLAB R2025b (R2020b+)
- Image Processing Toolbox (`multithresh`, `ssim`, `imquantize`)

---

## Citation

```bibtex
@article{wu2021ce,
  author  = {Wu, Hao-Tian and Zheng, Kaihan and Huang, Qi and Hu, Jiankun},
  title   = {Contrast Enhancement of Multiple Tissues in {MR} Brain Images With Reversibility},
  journal = {IEEE Signal Processing Letters},
  volume  = {28},
  pages   = {160--164},
  year    = {2021},
  doi     = {10.1109/LSP.2020.3048840}
}
```
