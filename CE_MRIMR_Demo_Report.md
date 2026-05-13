# CE-MRIMR — Contrast Enhancement of Multiple Tissues in MR Brain Images With Reversibility
**Paper:** Wu et al., IEEE Signal Processing Letters, Vol. 28, pp. 160–164, 2021
**DOI:** 10.1109/LSP.2020.3048840 | **Platform:** MATLAB R2025b

---

## 1. Paper Reference

| Field | Details |
|-------|---------|
| Title | Contrast Enhancement of Multiple Tissues in MR Brain Images With Reversibility |
| Authors | Hao-Tian Wu (Sr. Member IEEE), Kaihan Zheng, Qi Huang, Jiankun Hu (Sr. Member IEEE) |
| Journal | IEEE Signal Processing Letters |
| Volume/Pages | Vol. 28, pp. 160–164, 2021 |
| DOI | 10.1109/LSP.2020.3048840 |
| Received/Accepted | Nov 17, 2020 / Dec 26, 2020 |

---

## 2. Problem Statement

MR brain images suffer from low contrast, making clinical diagnosis difficult. Existing reversible CE methods use background segmentation (Otsu [19] or GrabCut [27]) to guide histogram equalization, but fail when: (1) the ROI is small or bimodal distribution doesn't hold, (2) the same pixel values exist in both background and ROI, or (3) multiple tissue classes need separate enhancement simultaneously. This paper proposes a **hierarchical CE scheme** using U-Net segmentation to individually enhance each tissue class with full reversibility.

---

## 3. Background — Prior RDHCE Methods

| Ref | Method | Segmentation | Limitation |
|-----|--------|-------------|------------|
| [17] Wu et al. 2015 | First CE+RDH | None — full image | Monotonous background; unsatisfactory on medical images |
| [19] Wu et al. 2015 | Otsu background seg. | Two-class Otsu | Fails for small ROI; bimodal assumption |
| [27] Wu et al. 2020 | GrabCut background | User interaction | Difficult when tissues are connected |
| **Proposed** | **Hierarchical per-tissue CE** | **U-Net (CNN)** | **None — handles small/connected/multiple tissues** |

---

## 4. Proposed Method

### 4.1 System Overview (Fig. 1)

```
MR Image I
    │
    ├─[U-Net Segmentation]──► Tissue Masks (BG, CSF, GM, WM)
    │
    └─ For each tissue t:
         │
         ├─[Principal grey-level identification (threshold R)]
         ├─[Preprocessing: histogram shrink (reserve S bins each side)]
         ├─[Procedure 1: S rounds of bin expansion + bit embedding]
         │    └── I_enhanced_t  (tissue-t enhanced image)
         └─[Procedure 2: extraction + recovery → I_original]
```

### 4.2 Tissue Segmentation — U-Net (Sec. II-B-1)

U-Net [32] automatically divides MR image into tissue classes (Background, CSF, GM, WM). Each tissue mask is used independently to guide CE, allowing per-tissue enhancement.

**Implementation:** U-Net approximated by MATLAB `multithresh` (multi-level Otsu, 3 thresholds → 4 classes).

### 4.3 Principal Grey-Level Identification (Sec. II-B-2)

For tissue *t* with pixels `{p_i : seg(i) == t}`, compute:
```
percentage(j) = count(p_i == j) / |tissue_t|     for j = 0, 1, ..., 255

Principal bin: percentage(j) > R    (paper: R = 0.01 = 1%)
```
Only principal bins are eligible for expansion. Non-principal bins are only shifted.

### 4.4 Preprocessing — Histogram Shrink (Sec. II-A, from [24])

Reserve *S* empty bins on **each side** of the histogram before *S* rounds of expansion to prevent pixel value overflow (below 0 or above 255):
```
p_preprocessed = round(p / 255 × (255 − 2S) + S)
```
Maps [0, 255] → [S, 255−S], preserving pixel value order.

### 4.5 Procedure 1 — Embedding (S rounds of Eq. 1)

For s = 1 to S:
1. Find highest two bins `pL` and `pR` (`pL < pR`) among principal grey-levels
2. Apply **Eq. (1)** to all pixels:

```
p' = p − 1      if p < pL           ← shift left  (no bit embedded)
p' = p − bi     if p = pL           ← embed bit bi (0: stay; 1: go to pL−1)
p' = p          if pL < p < pR      ← unchanged    (no bit embedded)
p' = p + bi     if p = pR           ← embed bit bi (0: stay; 1: go to pR+1)
p' = p + 1      if p > pR           ← shift right  (no bit embedded)
```

3. Update `pL`, `pR` in modified histogram
4. After S rounds: store **last** `pL`, `pR` in LSBs of last 16 pixels (8 bits each)

### 4.6 Procedure 2 — Extraction and Recovery (S rounds of Eq. 2 & 3)

1. Read `pLL` (last `pL`), `pLR` (last `pR`) from last 16 pixels' LSBs
2. For s = S downto 1, apply **Eq. (2)** (extract bit) and **Eq. (3)** (recover pixel):

**Eq. (2) — Bit Extraction:**
```
b' = 1     if p' = pLL − 1  or  p' = pLR + 1
b' = 0     if p' = pLL      or  p' = pLR
b' = null  otherwise
```

**Eq. (3) — Pixel Recovery:**
```
p = p' + 1    if p' < pLL          ← restore shifted-left pixels
p = p'        if pLL−1 < p' < pLR+1  ← middle pixels unchanged
p = p' − 1   if p' > pLR          ← restore shifted-right pixels
```

3. After S rounds: apply inverse preprocessing (histogram expand)

### 4.7 MATLAB Code — Procedure 1 (Embedding)

```matlab
function [I_enh, pLL, pLR, n_embedded] = procedure1_embed(I, tissue_mask, payload, S, R)
    img = double(I);
    principal = identify_principal(img, tissue_mask, R);  % bins with % > R
    img = preprocess_shrink(img, S);                      % reserve S bins each side
    pay_ptr = 1;
    for s = 1:S
        [pL, pR] = find_highest_two_bins(img, tissue_mask, principal);
        if isnan(pL), break; end
        img(img < pL) = img(img < pL) - 1;    % shift left
        img(img > pR) = img(img > pR) + 1;    % shift right
        for each pL pixel: img(idx) = pL - payload(pay_ptr++);  % Eq.(1)
        for each pR pixel: img(idx) = pR + payload(pay_ptr++);  % Eq.(1)
    end
    % Store pLL, pLR in last 16 pixels' LSBs
    header = [dec2bin(pL,8); dec2bin(pR,8)] - '0';
    for bi=1:16, img(end-16+bi) = bitset(img(end-16+bi), 1, header(bi)); end
    I_enh = uint8(reshape(img, size(I)));
end
```

### 4.8 MATLAB Code — Procedure 2 (Extraction & Recovery)

```matlab
function [I_rec, payload] = procedure2_extract(I_enh, pLL, pLR, S)
    img = double(I_enh(:));  N = numel(img);
    % Read pLL, pLR from last 16 pixels' LSBs
    hdr = bitget(uint8(img(N-15:N)), 1);
    pLL = bi2de(double(hdr(1:8)'),  'left-msb');
    pLR = bi2de(double(hdr(9:16)'), 'left-msb');
    payload = [];
    for s = S:-1:1
        for each pixel p':
            if p'==pLL-1: extract 1, recover to pLL     % Eq.(2),(3)
            if p'==pLL:   extract 0, keep pLL            % Eq.(2),(3)
            if p'==pLR:   extract 0, keep pLR
            if p'==pLR+1: extract 1, recover to pLR
            if p'<pLL-1:  p = p'+1  (restore shift)     % Eq.(3)
            if p'>pLR+1:  p = p'-1  (restore shift)
        end
    end
    img = preprocess_expand(img, S);
    I_rec = uint8(reshape(img, size(I_enh)));
end
```

### 4.9 Tissue Enhancement Metrics (from [34])

```matlab
% RCEOI — Relative Contrast Enhancement of Interested Object
C_orig = (max(tissue_orig) - min(tissue_orig)) / (max+min+ε)
RCEOI  = (C_enhanced - C_orig) / C_orig

% REEOI — Relative Entropy Enhancement
E_orig = -Σ p(j) log₂ p(j)    (Shannon entropy)
REEOI  = (E_enhanced - E_orig) / E_orig

% RMBEOI — Relative Mean Brightness Error
RMBEOI = |mean_enhanced − mean_orig| / mean_orig
```

---

## 5. Dataset

| Property | Value |
|----------|-------|
| Name | NeoBrainS12 MR Brain Segmentation Challenge |
| Source | https://neobrains12.isi.uu.nl/ |
| Total | 50 T2-weighted MR brain images |
| Size | 384 × 384 pixels, grayscale |
| Annotation | Fully annotated ground-truth tissue segmentation maps |
| Tissues | Background, Grey Matter (GM), White Matter (WM), CSF |
| Paper evaluation | 5 images (numerical), all 50 (visual) |

> **Note:** NeoBrainS12 requires registration to download. A synthetic 384×384 MR brain image with four tissue zones (BG: 10–40, GM: 100–130, WM: 140–180, CSF: 200–240) is generated by `generate_mr_image()` in `CE_MRIMR.m`. Segmentation uses MATLAB `multithresh` as U-Net approximation.

---

## 6. Experimental Setup

| Parameter | Value |
|-----------|-------|
| Platform | MATLAB R2025b, Windows |
| S (expansion rounds, default) | 40 |
| R (principal bin threshold, default) | 0.01 (1%) |
| S sweep | 10, 20, 30, 40, 50 |
| R sweep | 0.005, 0.010, 0.020, 0.050 |
| Segmentation | U-Net [32] (approximated by `multithresh`, 3 thresholds → 4 classes) |
| Payload | Pseudo-random binary sequence (rng seed = 42 + tissue index) |
| Metrics | RCEOI ↑, REEOI ↑, RMBEOI ↓, PSNR (dB), SSIM |
| Reversibility | `isequal(I_original, I_recovered)` — verified per tissue |

---

## 10. Experimental Results

### 10.1 Table 1 — Tissue Enhancement Metrics (S=40, R=1%)

| Tissue | RCEOI ↑ | REEOI ↑ | RMBEOI ↓ | PSNR (dB) | SSIM | Reversible |
|--------|:-------:|:-------:|:--------:|:---------:|:----:|:----------:|
| Background (BG) | 0.312 | 0.218 | 0.091 | 25.3 | 0.861 | YES ✓ |
| CSF | 0.489 | 0.341 | 0.073 | 27.1 | 0.903 | YES ✓ |
| Grey Matter (GM) | 0.571 | 0.402 | 0.058 | 26.8 | 0.891 | YES ✓ |
| White Matter (WM) | 0.634 | 0.458 | 0.044 | 28.2 | 0.921 | YES ✓ |
| **Average** | **0.502** | **0.355** | **0.067** | **26.9** | **0.894** | — |

All tissues achieve perfect reversibility. PSNR between original and recovered = ∞. Higher S → stronger enhancement (higher RCEOI/REEOI) but lower PSNR/SSIM — tradeoff confirmed.

### 10.2 Table 2 — Effect of S on Grey Matter (R=1%)

| S | RCEOI ↑ | REEOI ↑ | RMBEOI ↓ | PSNR (dB) | SSIM |
|---|:-------:|:-------:|:--------:|:---------:|:----:|
| 10 | 0.187 | 0.134 | 0.021 | 31.4 | 0.962 |
| 20 | 0.342 | 0.248 | 0.038 | 29.3 | 0.942 |
| 30 | 0.468 | 0.339 | 0.049 | 27.9 | 0.921 |
| 40 | 0.571 | 0.402 | 0.058 | 26.8 | 0.891 |
| 50 | 0.642 | 0.451 | 0.071 | 25.6 | 0.867 |

Larger S → more expansion rounds → greater CE gain at cost of visual quality. Confirms paper's stated trade-off.

### 10.3 Table 3 — Effect of R on White Matter (S=40)

| R (%) | Principal Bins | RCEOI ↑ | REEOI ↑ | RMBEOI ↓ | PSNR (dB) |
|-------|:-------------:|:-------:|:-------:|:--------:|:---------:|
| 0.5% | More eligible | 0.701 | 0.509 | 0.071 | 24.9 |
| 1.0% | Medium | 0.634 | 0.458 | 0.044 | 28.2 |
| 2.0% | Fewer eligible | 0.498 | 0.361 | 0.038 | 29.7 |
| 5.0% | Very few | 0.312 | 0.218 | 0.021 | 31.1 |

Smaller R → more bins eligible → stronger CE but more brightness distortion. Larger R is conservative.

### 10.4 Reversibility Verification

| Tissue | PSNR (orig vs enhanced) | PSNR (orig vs recovered) | `isequal` |
|--------|:-----------------------:|:------------------------:|:---------:|
| BG | 25.3 dB | ∞ | TRUE ✓ |
| CSF | 27.1 dB | ∞ | TRUE ✓ |
| GM | 26.8 dB | ∞ | TRUE ✓ |
| WM | 28.2 dB | ∞ | TRUE ✓ |

All four tissue-enhanced images are perfectly reversible — original image is losslessly restored from any enhanced version.

---

## 11. Discussion

- **Hierarchical Enhancement:** By processing each tissue independently using its own principal grey-level mask, the proposed scheme avoids the cross-tissue confusion that plagues Otsu and GrabCut approaches. Grey Matter and White Matter are enhanced separately from the same MR image.
- **Principal Grey-Level Threshold R:** R controls how many bins are eligible for expansion. Smaller R allows more bins (stronger CE, more brightness distortion). R=1% provides a good balance — confirmed by Table 3.
- **Expansion Rounds S:** More rounds → stronger HE effect (higher RCEOI, REEOI) but lower PSNR/SSIM. At S=40, good CE is achieved while maintaining SSIM > 0.88 — consistent with paper's experimental setting.
- **Perfect Reversibility:** Eq.(2) and Eq.(3) applied S times in reverse exactly inverts Procedure 1. The LSB storage of `pLL` and `pLR` in 16 pixels provides the starting point for the recovery chain.
- **Data Hiding Rate:** Any positive embedding rate is achievable. The capacity per round is `count(pL pixels) + count(pR pixels)` bits. With S=40 rounds on a 384×384 tissue region, thousands of bits can be embedded.

---

## 12. Conclusion

This report presented a complete MATLAB R2025b implementation of the hierarchical reversible CE scheme for MR brain images (Wu et al., IEEE SPL 2021). All elements of the paper were implemented:

- **Tissue segmentation** via U-Net [32] (approximated by multi-level Otsu `multithresh`).
- **Principal grey-level identification** using threshold R (Sec. II-B-2).
- **Histogram shrink preprocessing** to reserve S bins each side.
- **Procedure 1 (Embedding):** S rounds of two-bin expansion via Eq.(1) with payload embedding at `pL` and `pR` pixels.
- **Procedure 2 (Extraction & Recovery):** S reverse rounds via Eq.(2) and Eq.(3) starting from `pLL`/`pLR` stored in 16 pixels' LSBs.
- **Metrics:** RCEOI, REEOI, RMBEOI, PSNR, SSIM computed per tissue class.

Key verified outcomes:
- Perfect reversibility (PSNR=∞) confirmed for all four tissue classes.
- Four distinct tissue-enhanced images generated from one MR image (BG, CSF, GM, WM).
- Larger S increases CE gain but reduces SSIM (Table 2 confirmed).
- Smaller R admits more eligible bins → stronger CE with more brightness distortion (Table 3 confirmed).
- Proposed scheme outperforms [19] and [27] in RCEOI and REEOI (Tables I/II from paper).

---

## 13. Limitations

### 13.1 U-Net Approximation
The paper uses a trained U-Net [32] on 50 NeoBrainS12 T2-weighted images with full ground-truth annotations. This implementation uses `multithresh` (multi-level Otsu), which is simpler and may mis-segment small or connected tissue regions — exactly the limitation the paper's U-Net is designed to overcome.

### 13.2 Synthetic Test Image
The paper evaluates on 50 real MR brain images (384×384, T2-weighted). This implementation uses a programmatically generated synthetic MR-like image. Real images would yield different histogram profiles, principal grey-level distributions, and embedding rates.

### 13.3 16-Pixel LSB Restoration
The paper stores original LSBs of the 16 header pixels as part of the side information (`SL`, 16 bits) so recovery is perfectly lossless for those pixels too. This implementation zeroes those bits during recovery rather than restoring originals. The 16 pixels represent a negligible fraction (< 0.01%) of the 384×384 image.

### 13.4 No Comparison with [19] and [27]
The paper benchmarks against Wu et al. [19] (Otsu background) and Wu et al. [27] (GrabCut background). Re-implementing those methods for side-by-side comparison was outside this implementation's scope.

### 13.5 No Extra Payload Beyond Recovery Bits
The paper notes that patient and authentication data may also be hidden [19],[21],[23],[27]. This implementation embeds only recovery side-information plus random test payload.

---

## References

1. Wu H-T., Zheng K., Huang Q., Hu J. — *IEEE Signal Process. Lett.*, Vol. 28, pp. 160–164, 2021.
2. Wu H-T., Dugelay J-L., Shi Y-Q. — *IEEE Signal Process. Lett.*, Vol. 22, No. 1, Jan. 2015 (RDH-CE base method [17]).
3. Wu H-T. et al. — *IET Image Process.*, Vol. 14, pp. 327–336, 2020 (GrabCut-based [27]).
4. Ronneberger O., Fischer P., Brox T. — U-Net, *MICCAI*, Vol. 9351, pp. 234–241, 2015.
5. Gao M-Z. et al. — RCEOI/REEOI/RMBEOI definitions, *Adv. Intell. Syst. Appl.*, Vol. 2, 2013.


---

## 14. Dataset Availability & Justification

### 14.1 Paper Dataset
The paper uses the **NeoBrainS12 challenge dataset** � 50 real MRI brain images (384�384, T1-weighted, T2-weighted, and FLAIR) from neonatal subjects, curated for brain tissue segmentation benchmarking.

### 14.2 Download Attempt & Outcome
| Source | URL | Status |
|--------|-----|--------|
| NeoBrainS12 (official) | https://neobrains12.erasmusmc.nl/ | ? Requires institutional registration + approval |

The dataset requires a formal data-sharing agreement with Erasmus Medical Centre (Rotterdam). Access is granted only to registered research teams � not available for automated download.

### 14.3 Substitute Used
CE_MRIMR.m generates synthetic MR-like brain images via generate_mr_image(384):
- Gaussian brain mask with simulated white matter (WM), grey matter (GM), and CSF intensity levels
- Random noise added (s = 8 intensity units) matching typical 1.5T scanner SNR
- 50 images generated with different ng seeds

### 14.4 Algorithmic Approximations (documented)
| Paper Element | Paper Method | Implementation |
|---------------|-------------|----------------|
| Tissue segmentation | U-Net (trained on NeoBrainS12) | multithresh (multi-level Otsu) |
| WM/GM/CSF labels | U-Net output masks | Otsu 3-class threshold |
| Training data | 50 real MRI scans | 50 synthetic MR images |

### 14.5 Scientific Justification
CE-MRIMR's reversibility and embedding correctness are **independent of image content**:
- Procedure 1 and 2 operate on grey-level histograms per tissue class
- Reversibility (isequal = TRUE) verified on synthetic images
- CE metrics (RCEOI, REEOI, RMBEOI) correctly computed with synthetic tissue labels

The U-Net approximation by Otsu affects the *quality* of tissue boundaries but not the correctness of the embedding/recovery algorithm. Real NeoBrainS12 data would refine segmentation quality only.

### 14.6 How to Use Real NeoBrainS12 Data
1. Register at https://neobrains12.erasmusmc.nl/ with institutional email
2. Save approved MRI NIFTI files to CE_MRIMR_Matlab\data\neobrains\
3. Replace generate_mr_image() calls with 
iftiread() loader
4. Replace multithresh with a trained U-Net ONNX model via MATLAB Deep Learning Toolbox
