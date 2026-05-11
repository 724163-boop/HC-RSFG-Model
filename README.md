# HC-RSFG: Hierarchical Core Red Star-Forming Galaxy Theory

**English | 日本語 | 中文**

This repository contains the numerical simulation codes for the **HC‑RSFG (Hierarchical Core Red Star‑Forming Galaxy)** theory.  
This theory proposes that the extreme redness of galaxies — from the early Universe's Little Red Dots (LRDs) to present‑day relic galaxies — can be explained **entirely by pure stellar formation** with an extremely bottom‑heavy Initial Mass Function (IMF), **without requiring dust or an Active Galactic Nucleus (AGN)**.

🌐 **Project Website & Interactive Visualizations**: *[URL to be added]*

---

## 1. Project Overview | プロジェクト概要 | 项目概述

### English
The **HC‑RSFG theory** demonstrates that the red colours of JWST's Little Red Dots (LRDs) can be fully reproduced by pure star formation with an extremely bottom‑heavy IMF (few O/B/A stars), **without the need for dust extinction or AGN**.  
The theory provides a unified physical scenario that connects all evolutionary stages of massive galaxy cores:
- **Dark‑RSFG** (z≈10, currently undetected)
- **Little Red Dots (LRDs)** (z≈5–7)
- **Red Misfits** (z≈0.5)
- **Relic Galaxies** (e.g., NGC 1277, z≈0)

### 日本語
**HC‑RSFG（段階的コア形成RSFG）理論**は、JWSTが発見したリトル・レッド・ドット（LRD）の赤い色を、ダストやAGNに頼らず、極端なボトムヘビーな初期質量関数（O型・B型・A型星がほとんど作られない）のみで説明するものです。  
本理論は初期宇宙から現在に至る銀河の異常な赤みを説明可能とします。具体的には、**Dark‑RSFG**（z≈10, 未検出）、**LRD**（z≈5–7）、**Red Misfits**（z≈0.5）、そして現在の**レリック銀河**（NGC 1277）までを統一的に説明し、一貫した物理シナリオである可能性を論じます。

### 中文
**HC‑RSFG（分层核心红星星系形成）理论**表明，JWST观测到的小红点（LRD）的红色可以用纯星族合成来解释——其初始质量函数极度“底重”（缺乏大质量O/B/A星），**不需要尘埃消光或活动星系核**。  
该理论统一描述了从早期宇宙到现在的星系核心演化，包括：**Dark‑RSFG**（z≈10，不可见）、**LRD**（z≈5–7）、**Red Misfits**（z≈0.5）和**遗迹星系**（例如 NGC 1277）。

---

## 2. The Mystery of Little Red Dots & Existing Theories | LRDの謎と既存理論 | LRD的谜题与现有理论

Shortly after JWST began operations, astronomers noticed something unexpected: red, compact objects dubbed **“Little Red Dots” (LRDs)** appearing in large numbers in the very early Universe. Their extreme redness has sparked intense debate:

- **AGN (Active Galactic Nucleus) hypothesis**: Explains strong emission lines, but many LRDs are “X‑ray quiet” and some lack broad lines.
- **Dust‑obscured star formation hypothesis**: Explains the redness, but several LRDs show no dust at all (`Xiao et al. 2025`).

**The HC‑RSFG theory proposes a third, standalone solution**: the stars themselves are red because **massive stars simply never form**.

---

## 3. Theoretical Solution | 理論的解決策 | 理论解决方案

### Core Hypothesis | 核心仮説 | 核心假说
In the earliest stages of galaxy formation (z ≈ 7.5), the central core reaches extreme densities (`n ≈ 10⁷ cm⁻³`) and cools **30,000 times faster than it collapses**.  
This keeps the Jeans mass below ~1.5 M⊙, preventing the formation of O-, B-, and most A‑type stars.  
The resulting stellar population — dominated by K‑ and M‑dwarfs — is **intrinsically born red** without any dust.

### An example of the course of evolution | 進化の流れの一例 | 演化序列
| Phase | Redshift | Observable Counterpart | Status |
|:---|:---|:---|:---|
| **Dark‑RSFG** | z ≈ 7–10 | — (extremely faint) | **Theoretically predicted** |
| **LRD** | z ≈ 5–7 | Narrow‑line Little Red Dots | **Observed (JWST)** |
| **Red Misfit** | z ≈ 0.5 | Star‑forming red galaxies | **Observed (SDSS)** |
| **Relic Galaxy** | z ≈ 0 | NGC 1277, etc. | **Observed (Local)** |

---

## 4. Repository Contents | リポジトリの内容 | 仓库内容

All simulation codes are written in **Python 3** and require the **FSPS** library. They are designed for **Google Colab**.

| File | Description (EN) | 説明 (JP) | 说明 (ZH) |
|------|----------------|-----------|-----------|
| `HC_RSFG_Code1_LRD.py` | Phase 4: LRD state simulation (z=6, tage=0.9 Gyr) | Phase 4: LRD状態のシミュレーション | 阶段4：LRD态模拟 |
| `HC_RSFG_Code2_RedMisfit.py` | Phase 5: Red Misfit simulation (z=0.5, two‑component model) | Phase 5: Red Misfit状態のシミュレーション | 阶段5：Red Misfit态模拟 |
| `HC_RSFG_Code3_RelicGalaxy.py` | Phase 6: Relic galaxy simulation (NGC 1277 analog) | Phase 6: レリック銀河シミュレーション | 阶段6：遗迹星系模拟 |
| `HC_RSFG_Code4_DarkRSFG.py` | Phase 2–3: Dark‑RSFG simulation (undetected stage) | Phase 2–3: Dark‑RSFGシミュレーション | 阶段2–3：Dark‑RSFG模拟 |
| `hc_rsfg_hand_calculations.py` | Analytical verification (no FSPS needed) | 手計算による検証 | 手算验证 |

---

## 5. Key Numerical Results | 主要数値結果 | 关键数值结果

| Observable | HC‑RSFG Prediction | Observation | Source |
|:---|:---|:---|:---|
| F200W−F444W (LRD) | **+1.97 mag** (dust‑free) | > 1.0 mag | Hviding+2025 |
| EW(Hα) (LRD core) | < 1 Å | — (narrow lines only) | — |
| u−r (Red Misfit) | **+3.55 mag** | > 2.22 mag | Evans+2018 |
| M/L ratio (NGC 1277 centre) | **1.56 × Kroupa** (lower limit) | ~2.0 × Kroupa | Martín‑Navarro+2015 |

---

## 6. Setup & Requirements | 環境構築と必要条件 | 环境配置与依赖

```bash
# Install required Python packages
!pip install numpy fsps

The codes require FSPS (Flexible Stellar Population Synthesis), which is automatically cloned and compiled in the Google Colab environment cells provided in each script.
