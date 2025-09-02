TFPT — Topological Fixed Point Theory

# TFPT — Topological Fixed Point Theory

> **From topology to dynamics to observables**

[![Website](https://img.shields.io/badge/Website-fixpoint--theory.com-blue)](https://fixpoint-theory.com/)
[![Paper](https://img.shields.io/badge/Paper-v1.06-green)](https://claude.ai/chat/Paper%20V1.06%20-%2001.09.2025.pdf)
[![License](https://img.shields.io/badge/License-MIT-yellow)](https://claude.ai/chat/LICENSE)

## 🌐 Live Site & Playground

Explore the theory, visuals, and the interactive alpha playground here:

* **Main Site** :[fixpoint-theory.com](https://fixpoint-theory.com/)
* **Alpha Playground** : Interactive cubic solver and plots

The notebooks in this repo mirror what you see on the site and let you reproduce the numbers locally.

---

## 📦 What This Repository Contains

### 📓 Notebooks

* **`TFPT_alpha_only_v1.ipynb`**

  A focused, narrative notebook that derives the fine structure constant from the cubic fixed point in one sitting. It includes the one parameter normal form in $c_3$, exact and approximate solvers, error analysis, and side-by-side checks against CODATA 2022.
* **`TFPT_Playground_v1_06.ipynb`**

  A richer, explorable notebook: slider-style inputs for $c_3$ and $\varphi_0$, plots for the cubic, sensitivity panels, and quick links to the E₈ ladder, inflation tiles, and two loop fingerprints.

### 🛠 Tools

* **`Tools/E8ChainSolver/e8_orbit_engine`**

  A compact engine to build and traverse the E₈ nilpotent-orbit chain and to compute the log-exact ladder $\varphi_{n+1}=\varphi_n,\mathrm{e}^{-\gamma(n)}$.
  **Folders:**

  * `src/` — core routines for the chain, damping $\gamma(n)$, and ladder values
  * `data/` — source tables for orbit dimensions and labels
  * `tables/` — precomputed ladder snapshots
  * `plots/` — quick figures for reviews
  * `results/` — cached numeric outputs used by the notebooks
  * `reports/` — short, human-readable summaries of a run
  * `tests/` — minimal checks for chain uniqueness and ratio invariants
  * `notebooks/` — small dev notebooks for the engine
* **`GeneticAlgorithm/`**

  The search code that originally discovered the invariant patterns in candidate Lagrangians. Useful if you want to re-run the evolutionary search that converges on the same structural constants.
* **`Pyr@ate/`**

  Model files and exports for the two loop RGE runs with optional three loop SU(3) check. The YAML reflects the field content and thresholds used in the paper.

### 📄 Papers & Documentation

* **`Paper V1.06 - 01.09.2025.pdf`** and **`.md`**

  The complete reference: geometric derivation of $c_3$ and $\varphi_0$, the cubic fixed point for $\alpha$, the E₈ ladder, two loop fingerprints, inflation chapter, the chirality appendix, and a detailed changelog.
* **`TFPT Paper short.pdf`**

  A compressed narrative for quick reading and press contacts.

### 🔧 Build & Usage Files

* `pyproject.toml`,`setup.py` — package metadata for editable installs
* `requirements.txt` — minimal environment
* `run_demo.py` — command line demo for the orbit engine and ladder ratios
* `README_USAGE.md` — short task oriented recipes
* `README.md` — this file

---

## 💡 The Idea in One Page

TFPT claims that a handful of constants are not inputs but forced solutions of a stationary action across topology, geometry, and symmetry. Three pillars:

### 1. **Topology fixes the number**

$$
c_3 = \frac{1}{8\pi}
$$

via 11D Chern–Simons reduction with integer intersection on $Y_7$. The same scale appears in ABJ-type terms. See Sections 3.2.1 and Appendix E.

### 2. **Geometry fixes the length**

$$
\varphi_0 = \frac{1}{6\pi} + \frac{3}{256\pi^4}

$$

from Gauss–Bonnet on the orientable double cover of the Möbius fiber, with a canonical $6\pi$ boundary coefficient and a universal topological add-on. See Sections 3.2.2 and Appendix D.

### 3. **Symmetry orders the scales**

The E₈ nilpotent-orbit chain $D_n=60-2n$ yields a log-exact damping $\gamma(n)$ and the ladder $\varphi_n$. Its initial value $\gamma(0)=0.834$ is set by an extremal smoothing principle on the chain and makes the ratio laws calibration-free. See Section 4.1 to 4.5 and Appendix B.

From these three, a cubic fixed point for the electromagnetic coupling follows:

$$
\alpha^3 - 2c_3^3\alpha^2 - 8b_1c_3^6\ln\frac{1}{\varphi_0} = 0, \quad b_1 = \frac{41}{10}

$$

which reproduces $\alpha^{-1} \approx 137.03650146$ within a few parts per million, with **no fitting** . See Sections 7.2 to 7.6.

Two loop RGE runs show fingerprints of the fixed points on the QCD flow: $\alpha_3$ meets $\varphi_0$ near one PeV and $c_3$ near $2.5 \times 10^8$ GeV, and the three inverse couplings form a narrow corridor around $10^{15}$ GeV. See Section 5.2 and the plots on pages 23 to 24.

There is also a compact inflation chapter that ties the same invariants to an alpha-attractor plateau with $r \sim 2.5 \times 10^{-3}$ at $N \approx 57$. See Section 6.

---

## 🚀 Quick Start

```bash
# 1) Create a fresh environment
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# 2) Install
pip install -U pip
pip install -r requirements.txt
pip install -e .

# 3) Sanity check
python run_demo.py  # prints the E8 chain, ratio tests, and a few ladder values
```

### Open the notebooks:

* **Alpha only** : Start with `TFPT_alpha_only_v1.ipynb`

  You will compute $\alpha$ from the cubic, compare to CODATA 2022, and try the high accuracy Newton update.
* **Playground** : Open `TFPT_Playground_v1_06.ipynb`

  Change $\varphi_0$ within sensible bounds and watch the self-consistency loop and sensitivity bars. The E₈ ladder and the cubic plots update together.

---

## 🔬 Reproduce Key Results

### 1. **One parameter normal form**

In the alpha notebook, evaluate the closed form with $c_3=1/(8\pi)$ and $\varphi_0=\frac{1}{6\pi}+\frac{3}{256\pi^4}$. The Cardano branch and the fast approximation both land in the ppm window. See Section 7.3 and 7.5.

### 2. **E₈ ladder ratios without any unit choice**

Run `run_demo.py` or the orbit engine to confirm:

```
φ₁₂/φ₁₀ = (36/40)^λ
φ₁₅/φ₁₂ = (30/36)^λ
φ₂₅/φ₁₅ = (10/30)^λ
```

with $\lambda = \gamma(0)/(\ln 248 - \ln 60)$. See Section 4.3.

### 3. **Two loop fingerprints**

Use the Pyr@ate model to run the flow and replicate the markers at one PeV and $2.5 \times 10^8$ GeV and the corridor near $10^{15}$ GeV. See Section 5.2.

---

## 📁 Repository Map

```
TFPT/
├── 📓 Notebooks/
│   ├── TFPT_alpha_only_v1.ipynb
│   └── TFPT_Playground_v1_06.ipynb
├── 🛠 Tools/
│   └── E8ChainSolver/e8_orbit_engine/
│       ├── src/         # chain construction, damping, ladder
│       ├── data/        # orbit tables and labels
│       ├── tables/      # precomputed ladders
│       ├── results/     # cached outputs
│       ├── plots/       # quick figures
│       ├── reports/     # text summaries
│       ├── tests/       # unit checks and invariants
│       └── notebooks/   # small dev notebooks
├── 🧬 GeneticAlgorithm/    # GA that unearthed the invariant pattern
├── 📊 Pyr@ate/             # model and export for two loop runs
├── 🎯 run_demo.py          # quick demo for the ladder and ratios
├── 📦 pyproject.toml       # packaging
├── 📦 setup.py             # legacy installer
├── 📋 requirements.txt     # pinned versions
├── 📖 README_USAGE.md      # task recipes
├── 📄 Paper V1.06 ...md    # paper in markdown
├── 📄 Paper V1.06 ...pdf   # full paper
└── 📄 TFPT Paper short.pdf # compressed version
```

---

## 🎯 Design Choices

* **Everything is testable**

  The orbit engine exposes pure functions. Ratio laws and uniqueness of the chain are covered by unit tests.
* **Separation of concerns**

  Notebooks tell the story. Tools do the math. Papers justify the steps with references and proofs.
* **No fitting**

  Numbers come from topology, geometry, and group data. Where a unit is inevitable, it is a block-level calibration rather than a knob. See Sections 8.1 to 8.5.

---

## 🔧 How to Extend

### Add a new observable block

Define $(r_B, k_B, n_B)$, compute $\zeta_B = (\pi c_3),\mathrm{e}^{-\beta_B\pi c_3}\mathrm{e}^{-k_B/c_3}$ with $\beta_B = (8-r_B)/8$, and set $X_B = \zeta_B M_{\mathrm{Pl}}\varphi_{n_B}$. See Sections 8.1.1 and 8.1.2.

### Add a chain diagnostic

Use the engine's exported $\ln D_n$ to compute second differences and the extremal functional from Section 4.5, then compare $\gamma(0)$ against 0.834.

### Add a new flow scenario

Duplicate the Pyr@ate YAML, change threshold scales in one decade both ways, and record the corridor stability and the two fingerprints from Section 5.2.

---

## 🌐 How the Website Relates to This Code

The website uses the same formulas and tables:

* The**Alpha page** evaluates the cubic and shows the exact and approximate roots for the same $c_3$ and $\varphi_0$ as the paper.
* The**E₈ ladder widget** reads the same $D_n=60-2n$ chain and the identical log-exact $\gamma(n)$.
* The**RGE view** plots the corridor and the two fingerprints from the model in Pyr@ate.

If you want to align the site and the notebooks perfectly, run the orbit engine once and export `tables/ladder_v1_06.csv` for the site build.

---

## 📚 Citation

If this work helps your research or product, please cite the paper:

> Hamann, S. *Von Topologie zu Dynamik: Die Ordnung hinter α und den Naturkonstanten* . Version 1.0.6, 01.09.2025.

### BibLaTeX

```bibtex
@misc{hamann2025tfpt,
  title  = {Von Topologie zu Dynamik: Die Ordnung hinter α und den Naturkonstanten},
  author = {Hamann, Stefan},
  year   = {2025},
  note   = {Version 1.0.6, 01.09.2025},
  url    = {https://fixpoint-theory.com}
}
```

---

## 🤝 Contributing

Issues and pull requests are welcome. High value items right now:

* Formal algebraic proof for the closed $\gamma(n)$ form from E₈ orbit theory
* Full two loop electroweak sector with consistent threshold matching
* Cleaner tests and fixtures for the chain uniqueness algorithm

Please keep the tone scientific, add checks, and reference the exact section in the paper that your change touches.

---

## 📜 License

Add your license of choice here if not already present. A permissive license makes reproduction and independent checks much easier.

---

## ✨ Final Note

> **If you only remember one thing, remember this: constants are not knobs. They are what remains when topology, geometry, and symmetry agree. The rest of this repo is a set of tools and proofs to make that statement reproducible.**
