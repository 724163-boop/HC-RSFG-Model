# HC-RSFG-Model
Source code for numerical simulations and theoretical modeling of hierarchical core red star-forming galaxies (HC-RSFG) to explain the origin of Little Red Dots (LRDs) in the early universe.
# HC-RSFG: Hierarchical Core Red Star-Forming Galaxy Theory

**English | 日本語 | 中文**

This repository contains the numerical simulation codes for the **HC-RSFG (Hierarchical Core Red Star-Forming Galaxy)** theory, which proposes a new galaxy formation scenario in which massive elliptical galaxy cores are built by a "bottom-heavy" initial mass function (IMF) in the early Universe without requiring dust or AGN.

🌐 **Live Interactive Website**: [https://your-username.github.io/HC-RSFG](https://your-username.github.io/HC-RSFG)  
*(Replace with your actual URL)*

---

## 1. Project Overview | プロジェクト概要 | 项目概述

### English
The **HC-RSFG theory** demonstrates that the red colors of Little Red Dots (LRDs) observed by JWST can be fully reproduced by pure star formation with an extremely bottom-heavy IMF (few O/B/A stars), without the need for dust extinction or active galactic nuclei (AGN).  
The theory spans the whole evolution of massive galaxy cores from the early Universe (Phase 0–6), linking **Dark-RSFG** (z≈10, invisible), **LRDs** (z≈5–7), **Red Misfits** (z≈0.5), and present-day **relic galaxies** (NGC 1277).

### 日本語
**HC-RSFG（段階的コア形成RSFG）理論**は、JWSTが発見したリトル・レッド・ドット（LRD）の赤い色を、ダストやAGNに頼らず、極端なボトムヘビーな初期質量関数（O型・B型・A型星がほとんど作られない）のみで説明するものです。  
本理論は初期宇宙から現在に至る銀河中心核の進化を統一的に記述し、**Dark-RSFG**（z≈10, 未検出）、**LRD**（z≈5–7）、**Red Misfits**（z≈0.5）、そして現在の**レリック銀河**（NGC 1277）まで、一貫した物理シナリオを提供します。

### 中文
**HC-RSFG（分层核心红星星系形成）理论**表明，JWST观测到的小红点（LRD）的红色可以用纯星族合成来解释——其初始质量函数极度“底重”（缺乏大质量O/B/A星），不需要尘埃消光或活动星系核。  
该理论统一描述了从早期宇宙到现在的星系核心演化：**Dark-RSFG**（z≈10，不可见）→ **LRD**（z≈5–7）→ **Red Misfits**（z≈0.5）→ **遗迹星系**（NGC 1277）。

---

## 2. Repository Contents | リポジトリの内容 | 仓库内容

All simulation codes are written in **Python 3** and require the **FSPS** (Flexible Stellar Population Synthesis) library. They are designed to run on **Google Colab** (free GPU environment).  
The **`/docs`** directory contains the source files for the interactive website (HTML, CSS, JavaScript) that visualizes the simulation results.

| File | Description (EN) | 説明 (JP) | 说明 (ZH) |
|------|----------------|-----------|-----------|
| `HC_RSFG_Code1_LRD.py` | Phase 4: LRD state simulation (z=6, tage=0.9 Gyr) | Phase 4: LRD状態のシミュレーション | 阶段4：LRD态模拟 |
| `HC_RSFG_Code2_RedMisfit.py` | Phase 5: Red Misfit simulation (z=0.5, two-component model) | Phase 5: Red Misfit状態のシミュレーション（コア＋円盤の2成分モデル） | 阶段5：Red Misfit态模拟（核+盘双成分） |
| `HC_RSFG_Code3_RelicGalaxy.py` | Phase 6: Relic galaxy simulation (NGC 1277 analog, radial IMF gradient) | Phase 6: レリック銀河シミュレーション（NGC 1277類似、半径方向IMF勾配） | 阶段6：遗迹星系模拟（NGC 1277 类似，径向IMF梯度） |
| `HC_RSFG_Code4_DarkRSFG.py` | Phase 2–3: Dark-RSFG simulation (z≈7.5–10, detection feasibility) | Phase 2–3: Dark-RSFGシミュレーション（z≈7.5–10、JWST検出可能性） | 阶段2–3：Dark-RSFG模拟（z≈7.5–10，JWST可探测性） |
| `/docs` | Interactive website (HTML, CSS, JS) | インタラクティブウェブサイト | 交互式网站 |

---

## 3. Interactive Website | ウェブサイト | 交互式网站

🌐 **URL**: [https://your-username.github.io/HC-RSFG](https://your-username.github.io/HC-RSFG)  
*(Replace with your actual URL)*

The website provides an interactive exploration of the HC-RSFG simulation results without the need to install any software. You can:
- View the predicted spectra and color evolution for each evolutionary phase.
- Compare the HC-RSFG model with standard IMFs (Chabrier, Kroupa).
- Examine the emission line equivalent widths and their time evolution.
- Explore the radial IMF gradient in the relic galaxy model.

The website is hosted using GitHub Pages from the `/docs` folder of this repository.

---

## 4. Quick Start (Google Colab) | クイックスタート | 快速启动

1. Open a new Colab notebook.
2. Copy the entire content of one of the code files into a cell.
3. Run the cell. The first execution will:
   - Clone FSPS from GitHub
   - Set the `SPS_HOME` environment variable
   - Install `fsps` Python wrapper
   - Initialize the stellar population with the HC-RSFG IMF
   - Produce the results and save a PNG figure
4. The output PNG will be saved to the Colab runtime filesystem (download from the left panel).

**※ Note**: The first run may take a few minutes because FSPS needs to compile its Fortran SSP generators.

---

## 5. Key Parameters | 主要パラメータ | 关键参数

All codes share the same **broken power-law IMF** to realize the extreme bottom-heavy star formation:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `imf_type` | 0 | Broken power-law IMF |
| `imf1` | 1.8 | Slope for 0.08–0.5 Msun |
| `imf2` | 1.8 | Slope for 0.5–1.5 Msun |
| `imf3` | 3.0–3.5 | Slope for >1.5 Msun (steep to suppress O/B/A stars) |
| `imf_upper_limit` | 1.5 Msun | Effective upper mass cutoff (in some codes) |
| `logzsol` | –2.0 to 0.0 | Metallicity (depends on phase) |
| `dust2` | 0.0–0.08 | Dust extinction (usually very low or zero) |

---

## 6. Citation & References | 引用文献 | 参考文献

If you use these codes or the HC-RSFG theory in your research, please cite the following papers:

- **Steinhardt (2025)**, *ApJ*, 982, 189 – original RSFG concept.
- **Hviding et al. (2025)**, *A&A*, 702, A57 – RUBIES survey of LRDs.
- **Martín-Navarro et al. (2015)**, *MNRAS*, 451, 1081 – NGC 1277 bottom-heavy IMF.
- **Evans et al. (2018)**, *MNRAS*, 476, 5284 – Red Misfits in SDSS.

(You can also cite this repository directly via a Zenodo DOI once archived.)

---

## 7. License | ライセンス | 许可证

This project is released under the **MIT License**. See the `LICENSE` file for details.  
You are free to use, modify, and distribute the code, provided that proper attribution is given.

---

## 8. Contact | お問い合わせ | 联系方式

For questions or collaboration requests, please open an issue on this GitHub repository or contact the author (724163@kng.ed.jp).
