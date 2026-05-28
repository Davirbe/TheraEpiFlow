# TheraEpiFlow — Validation Suite Report

> Timings in this report were measured on a Lenovo IdeaPad with AMD Ryzen 5 5500U, AMD Radeon iGPU, 8 GB RAM, running WSL2 (Ubuntu) on Windows. All replicates ran one-at-a-time on a quiet machine (no Chrome / IDE / Teams). Absolute seconds do not generalize to other hardware; compare _ratios_ (tr3 : tr1, SARS_NCAP : HPV16_E7) and the _shape_ of the loss curves.

### Execution environment

| Field | Value |
|---|---|
| Captured at | `2026-05-27T14:12:23` |
| Platform | `Linux-6.6.87.2-microsoft-standard-WSL2-x86_64-with-glibc2.39` |
| WSL | `True` |
| Kernel | `Linux DESKTOP-E5VJ9RK 6.6.87.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun  5 18:30:46 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux` |
| Python | `3.10.20` |
| netMHCpan | `?` |

### Package versions

| Package | Version |
|---|---|
| `biopython` | `1.87` |
| `jinja2` | `3.1.6` |
| `joblib` | `1.5.3` |
| `matplotlib` | `3.10.9` |
| `matplotlib-venn` | `1.1.2` |
| `mhcflurry` | `2.0.6` |
| `networkx` | `3.4.2` |
| `numpy` | `1.26.4` |
| `openpyxl` | `3.1.5` |
| `pandas` | `2.3.3` |
| `requests` | `2.33.1` |
| `rich` | `15.0.0` |
| `scikit-learn` | `1.2.2` |
| `seaborn` | `0.13.2` |
| `tensorflow` | `2.15.0` |
| `toxinpred3` | `1.4` |
| `xlsxwriter` | `3.2.9` |

## Experiment 1 — Execution time + reproducibility

- Started: `2026-05-27T14:12:23`
- Ended:   `2026-05-27T16:03:28`
- Total wall: **108.1 min**

### Replicate timing summary

| Preset | n | median | IQR | min | max |
|---|---|---|---|---|---|
| `hpv16_e7_tr1` | 10 | 1.2 min | 1.2 min – 1.3 min | 38.3s | 2.1 min |
| `hpv16_e7_tr3` | 10 | 3.6 min | 3.6 min – 3.7 min | 3.2 min | 3.7 min |
| `sars_nucleo_tr1` | 10 | 1.3 min | 1.3 min – 1.3 min | 39.2s | 1.4 min |
| `sars_nucleo_tr3` | 10 | 4.4 min | 4.3 min – 4.5 min | 3.9 min | 4.5 min |

### Figures

![exp1 loss_by_stage hpv16_e7_tr1](exp1_loss_by_stage_hpv16_e7_tr1.png)

![exp1 time_by_stage hpv16_e7_tr1](exp1_time_by_stage_hpv16_e7_tr1.png)

![exp1 consistency_matrix hpv16_e7_tr1](exp1_consistency_matrix_hpv16_e7_tr1.png)

![exp1 venn hpv16_e7_tr1](exp1_venn_hpv16_e7_tr1.png)

![exp1 loss_by_stage hpv16_e7_tr3](exp1_loss_by_stage_hpv16_e7_tr3.png)

![exp1 time_by_stage hpv16_e7_tr3](exp1_time_by_stage_hpv16_e7_tr3.png)

![exp1 consistency_matrix hpv16_e7_tr3](exp1_consistency_matrix_hpv16_e7_tr3.png)

![exp1 venn hpv16_e7_tr3](exp1_venn_hpv16_e7_tr3.png)

![exp1 loss_by_stage sars_nucleo_tr1](exp1_loss_by_stage_sars_nucleo_tr1.png)

![exp1 time_by_stage sars_nucleo_tr1](exp1_time_by_stage_sars_nucleo_tr1.png)

![exp1 consistency_matrix sars_nucleo_tr1](exp1_consistency_matrix_sars_nucleo_tr1.png)

![exp1 venn sars_nucleo_tr1](exp1_venn_sars_nucleo_tr1.png)

![exp1 loss_by_stage sars_nucleo_tr3](exp1_loss_by_stage_sars_nucleo_tr3.png)

![exp1 time_by_stage sars_nucleo_tr3](exp1_time_by_stage_sars_nucleo_tr3.png)

![exp1 consistency_matrix sars_nucleo_tr3](exp1_consistency_matrix_sars_nucleo_tr3.png)

![exp1 venn sars_nucleo_tr3](exp1_venn_sars_nucleo_tr3.png)

![exp1_time_tr1_vs_tr3.png](exp1_time_tr1_vs_tr3.png)

![exp1_time_protein_size.png](exp1_time_protein_size.png)

## Experiment 2 — Clinical validation (IEDB IQ-API)

- Started: `2026-05-27T16:13:15`
- Ended:   `2026-05-27T16:57:35`
- Wall: **43.0 min**
- Pipeline replicates: 10
- Median pipeline wall per replicate: 2.8 min

### Literature support summary

| Track | ★ unique | ★ with hits | neg.control unique | neg.control with hits |
|---|---|---|---|---|
| `H16_E7` | 7 | 6 | 10 | 4 |
| `SARS2_N` | 30 | 20 | 10 | 0 |

### Figures

![exp2 IEDB counts H16_E7](exp2_iedb_counts_H16_E7.png)

![exp2 IEDB counts SARS2_N](exp2_iedb_counts_SARS2_N.png)

![exp2 negative-control comparison](exp2_neg_control_comparison.png)

_Supplementary publications table (PMID + title): `validation/results/exp2/iedb_pmid_titles.csv`_
