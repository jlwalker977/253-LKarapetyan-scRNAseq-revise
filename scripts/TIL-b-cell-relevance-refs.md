references for TIL B cell signal being biologically relevant

Yes, here are the key references supporting tumor-infiltrating B cells (TIBs) as biologically real and clinically meaningful:

---

## Landmark Papers

### 1. Tertiary Lymphoid Structures (TLS) — The Most Cited Evidence
Three papers published together in **Nature (2020)** established TIL B cells as clinically significant:

- **Cabrita et al. (2020)** *Nature* — Showed that B cells organized into TLS in melanoma were associated with improved survival and response to immunotherapy
- **Petitprez et al. (2020)** *Nature* — Demonstrated TLS with B cells predicted response to anti-PD-1 in soft tissue sarcoma
- **Helmink et al. (2020)** *Nature* — Found TLS-associated B cells in melanoma and renal cell carcinoma correlated with checkpoint blockade response

---

### 2. Functional Role of TIL B Cells
- **Hollern et al. (2019)** *Cell* — B cells and T follicular helper cells cooperate to drive anti-tumor immunity and immunotherapy response in breast cancer
- **Germain et al. (2021)** *Nature Cancer* — Single-cell characterization of B cells in tumors showing plasma cell and germinal center B cell states

---

### 3. Single-Cell Evidence Specifically
- **Schelker et al. (2017)** *Nature Communications* — scRNA-seq deconvolution confirming B cell presence in tumor microenvironments
- **Andreatta et al. (2021)** *Nature Communications* — ProjecTILs reference atlas including B cell populations in TIL scRNA-seq data

---

## Most Relevant to Your Data Specifically

Given your dataset involves **anti-PD-1 (nivo) and combination immunotherapy**, the most directly relevant are:

- **Cabrita et al. 2020** — B cells + TLS predict anti-PD-1 response in melanoma
- **Petitprez et al. 2020** — strongest evidence that B cell presence predicts immunotherapy benefit

---

## Quick Search to Verify Current Citations

Let me pull current citation counts and any newer follow-up papers:Here are the key references, now verified:

---

## Primary landmark references

**The three Nature 2020 papers** — published simultaneously and most widely cited:

- **Cabrita et al. (2020)** *Nature* 577:561–5 — showed that improved survival after immune checkpoint blockade was associated with co-occurrence of CD8+ T cells and CD20+ B cells in metastatic melanoma, with TLS formation confirmed by CXCR5/CXCL13 staining

- **Petitprez et al. (2020)** *Nature* 577:556–60 — demonstrated B cells are associated with survival and immunotherapy response in sarcoma

- **Vanhersecke et al. (2021)** *Nature Cancer* — mature tertiary lymphoid structures predict immune checkpoint inhibitor efficacy in solid tumors independently of PD-L1 expression

---

## Most Relevant to Your Dataset (anti-PD-1/nivolumab)

**Carril-Ajuria et al.** — baseline circulating unswitched memory B cells and B-cell related soluble factors are associated with overall survival in patients with clear cell renal cell carcinoma treated with **nivolumab** — directly relevant given your nivo arm.

---

## scRNA-seq Specific Evidence

**Zhu et al. (2025)** *Cell* — a large-scale atlas of 1,254,749 tumor-infiltrating immune cells from 234 NSCLC patients post-neoadjuvant chemo-immunotherapy resolved 51 immune cell subtypes including B cells, with **memory B cells elevated in major pathological responders**

**Cillo et al. (2021)** *Nature Communications* 12:3349 — demonstrated in HNSCC that germinal center TIL-B cells and TLS with germinal centers both correlate with favorable outcome, characterized at single-cell resolution

---

## Review for Background/Introduction

**BJC Reports (2025)** — comprehensive review: "Recent clinical trials involving sarcoma, melanoma, and renal cell carcinoma have shown that TLSs and tumour-infiltrating B cells are associated with a better response to immunotherapy and improved survival"

---

## Citation Summary for Your Methods/Discussion

| Reference | Journal | Relevance to your data |
|-----------|---------|----------------------|
| Cabrita et al. 2020 | *Nature* | TIL-B + TLS → immunotherapy response, melanoma |
| Petitprez et al. 2020 | *Nature* | TIL-B → survival, sarcoma |
| Vanhersecke et al. 2021 | *Nature Cancer* | Mature TLS → ICB efficacy across solid tumors |
| Carril-Ajuria et al. | — | B cells + **nivolumab** specifically |
| Zhu et al. 2025 | *Cell* | scRNA-seq atlas confirming B cell subtypes in ICB response |
| Cillo et al. 2021 | *Nat Commun* | scRNA-seq B cell characterization in TLS |