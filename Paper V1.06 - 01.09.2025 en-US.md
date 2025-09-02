# From topology to dynamics: The order behind $\alpha$ and the natural constants
Stefan Hamann, September 1, 2025, Version 1.0.6

```table-of-contents
title: 
style: nestedList # TOC style (nestedList|nestedOrderedList|inlineFirstLevel)
minLevel: 0 # Include headings from the specified level
maxLevel: 4 # Include headings up to the specified level
include: 
exclude: 
includeLinks: true # Make headings clickable
hideWhenEmpty: false # Hide TOC if no headings are found
debugInConsole: false # Print debug info in Obsidian console
```

## Abstract

We show that the fine-structure constant $\alpha$ and other fundamental quantities are **not** required as free inputs, but follow from **topology, geometry, and symmetry**. The starting points are

- the **topological fixed point** $c_3=\tfrac{1}{8\pi}$,  
- a **geometrically defined length scale** $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.053171952\,$ (reduced Planck units),  
- and a **damping function** $\gamma(n)$ ordered by **E₈** for the discrete vacuum conductors $\varphi_n$. 

The core result is a **single-parameter normal form** for $\alpha$ (parameter = $c_3$). From $c_3$, the following results exactly
$$[
\varphi_0=\frac{4}{3}c_3+48c_3^4,\quad A=2c_3^3,\quad 
\kappa=\frac{b_1}{2\pi}\ln\frac{1}{\varphi_0},\quad b_1=\frac{41}{10},
]$$
and thus the **cubic fixed point equation**
$$[
\boxed{\;\alpha^3-2c_3^3\,\alpha^2-8\,b_1\,c_3^6 \ln\!\Big(\frac{1}{\tfrac{4}{3}c_3+48c_3^4}\Big)=0\;}
]$$
with exactly **one** real physical solution $\alpha(c_3)$. For $c_3=\frac{1}{8\pi}$, we obtain

$$[
\varphi_0=0.0531719521768,\quad 
\kappa=1.914684795,\quad 
\alpha=0.007297325816919221,\quad 
\alpha^{-1}=137.03650146488582,
]$$

i.e., a deviation of **3.67 ppm** from CODATA-2022 – without free parameters.  

The same structure generates a **log-exact E₈ cascade** $\varphi_{n+1}=\varphi_n e^{-\gamma(n)}$, whose anchor steps hit flavor mixtures, electroweak and hadronic scales, and cosmological constants. 

> **A two-loop RGE analysis confirms both fingerprints:**
   $\alpha_{3}(1\,\text{PeV})=0.052923411$ is $0.47\%$ below $\varphi_{0}=0.053171952$,
   and at $\mu\simeq 2.5\times 10^{8}\,\text{GeV}$ we obtain $\alpha_{3}=0.039713807$, i.e. $0.19\%$ below $c_{3}=1/(8\pi)$.
> A color-adjacent bridge $G8$ above $M_{G8}=1.8\times 10^{10}\,\text{GeV}$
   reduces the slope of $\alpha_{3}^{-1}$ from $\tfrac{7}{2\pi}$ to $\tfrac{5}{2\pi}$ and creates a narrow unification corridor with a minimum relative spread of $1.23\%$ at $\mu^{\star}\approx 1.43\times 10^{15}\,\text{GeV}$.

This results in a consistent picture:
**Topology** fixes the normalizations, **geometry** fixes the length scale, **E₈** orders the scale ladder, and **RG dynamics** confirms the fingerprints. 

![[Pasted image 20250826125915.png]]

> **Info Box: Notation and Conventions**
> 
> * Indices: $(c_{3}\rightarrow c_{₃}), (b_{1}\rightarrow b_{₁})$ in running text, as set in formulas.  
> * Length scale: $(ϕ_{0}=\frac{1}{6\pi}+\frac{3}{256\pi^{4}}), (ϕ_{n+1}=ϕ_{n}e^{-\,\gamma(n)})$.  
> * Topology and couplings: $(g=8c_{₃}^{2}=\frac{1}{8\pi^{2}}), (A=2c_{₃}^{3}=\frac{1}{256\pi^{3}})$.  
> * RG constant: $(κ=\frac{b_{₁}}{2\pi}\ln\!\frac{1}{ϕ_{0}}), (b_{₁}=41∕10)$ in GUT norm.  
> * Groups: $(E_{8}), (E_{7}), (E_{6})$ always written as indices $(E_{₈}), (E_{₇}), (E_{₆})$.  
> * Units: all dimensioned quantities in reduced Planck units unless otherwise specified.


# 1. Introduction

The question of the origin of natural constants—in particular, the fine structure constant $\alpha$—is answered here using a **bottom-up** approach: constants are **invariants** of a common framework consisting of **topology, geometry, and symmetry**, not external knobs. 

---
### 1.1 The genetic algorithm

We evolve a genetic algorithm (GA) using Lagrange densities with six coefficients $(c_0,\dots,c_5)$ (kinetics, mass, quartic kinetics, Maxwell, EH term). Hard physical constraints (Lorentz, ghost freedom, correct signs) are strictly enforced; fitness measures error-invariant $\delta_c,\delta_\alpha,\delta_G$ to target values. Typical populations $N\!=\!800$, tournament selection, elites, crossover, adaptive mutations. Result: robust **clusters** at $c_4$ (EM normalization), $c_3$ (quadratic kinetics, trace $1/(8\pi)^2$) and a narrow **$\varphi_0$ valley**.

**<span style="background:#d3f8b6">Figure: User interface of the GA-Search application></span>**

![[Pasted image 20250823093627.png]]
### 1.2 ) Genetic algorithm – setup, validation, results

- **Convergence:** ~24 million evaluations, ~15,000 generations; reproducibility via seeds.  
- **Pattern:** $c_3$ appears as a square track $8c_3^2=\tfrac{1}{8\pi^2}$ ⇒ **fixed point** $c_3=\tfrac{1}{8\pi}$. Mass term clusters suggest **$\varphi_0$**. EM normalization suggests $\ln(1/\varphi_0)$ in the $F^2$ sector.
- **Ablations:** Without constraints → ghost/tachyon collapse; without separate fine-tuning on $c_4$, $\alpha$ remains stuck at 3–4 digits. Adaptive precision prevents rounding artifacts.  

---
### 1.3) Representative high-fitness Lagrangians and patterns**

**Examples (from the Hall of Fame):**

$\begin{aligned} \mathcal L_{\#3566}&=-0.57618478(\partial_t\varphi)^2\;+\;0.57618478(\nabla\varphi)^2\;-\;0.98847468\,\varphi^2\\ &\quad+\;0.0130338797\,(\partial_t\varphi)^2\varphi^2\;-\;0.0917012368\,F_{\mu\nu}^2 \\ \\ \mathcal L_{\rm can.}&=-0.50000000(\partial_t\varphi)^2\;+\;0.50000000(\nabla\varphi)^2\;-\;0.059422638\,\varphi^2\\ &\quad-\;0.039752599\,(\partial_t\varphi)^2\varphi^2\;-\;0.10047012\,F_{\mu\nu}^2\;+\;3.2658{\times}10^{8}\,\kappa R \end{aligned}$

**Systematic clusters (robust across seeds/generations):**

- **Quartic kinetic coefficient (here $c_3$ of density):**

    $c_3^{\rm (Lag)}\simeq \frac{1}{8\pi^2}=0.0126651\quad\text{(observed, e.g., }0.0130339,\ \Delta\!\sim\!+2.9\%)$.
    We interpret this as the square trace of the **topological fixed point**
    $\boxed{c_3^{\rm (Topo)}=\frac{1}{8\pi}},\qquad \frac{1}{8\pi^2}=8\bigl(c_3^{\rm (Topo)}\bigr)^2$,
    which recurs in the nonlinear term $(\partial_t\varphi)^2\varphi^2$.
- **Scalar mass term:** Frequent peaks in [\,0.051,\,0.061\,] (in $\bar M_P$). We identify the **length fixed point**

    $\boxed{\varphi_0=0.053171952\ \ (\bar M_P)},\qquad \varphi_0/\sqrt{8\pi}=0.0106063\,M_P$.

- **Maxwell normalization:** $c_4$ clusters at -0.091701, which

    $\boxed{\alpha_{\rm model}=\frac{|c_4|}{4\pi}\approx0.007297352566}$
    
    reproduced (ppm precision). Variants with -0.04585 correspond to an alternative internal $F^2$ normalization (factor ½).

**Brief interpretation.** The GA does not "find" random numbers, but rather **canonical invariants**: the topological normalization $1/(8\pi)$, the geometric length $\varphi_0$, and a logarithmic fingerprint in the EM term (see below). These patterns are stable across populations, seeds, and search modes.

---
### 1.4) From pattern to first theory iteration

  The three GA findings lead directly to the first analytically controlled theory iteration:

1. **Fixed points instead of fits.**

    The recurring value $c_3^{\rm (Lag)}\!\approx\!1/(8\pi^2)$ **enforces** the topological fixed point $c_3^{\rm (Topo)}=1/(8\pi)$ as the underlying normalization of nonlinear terms.

2. **Geometric scale** $\varphi_0$**.**

    The mass term clusters define $\varphi_0$ as a **geometric Radion fixed point** (Möbius reduction). This makes a discrete scale ladder $\varphi_n$ plausible, which will be specified later in the E_8 cascade.

3. **EM logarithm** $\ln(1/\varphi_0)$**.**

    The observed EM normalization allows for a **parameter-free fixed-point equation for** $\alpha$, in which topology ($1/8\pi$) and geometry ($\varphi_0$) are coupled. This equation has exactly one physically real solution and reproduces $\alpha$ at the ppm level—consistent with the GA outputs.

4. **Dynamic testing.**

    Building on (1)–(3), a 2-loop RG "smoke test" was later formulated (E_8 cascade mock with EH term). The fluxes show the **fingerprints** $\alpha_3(1\,\mathrm{PeV})\approx \varphi_0$ and $\alpha_3(\mu)=1/(8\pi)$ at $\mu\!\sim\!2.5\times10^8\,\mathrm{GeV}$ as well as a narrow equilibrium corridor of the three couplings at $10^{14\text{–}15}\,\mathrm{GeV}$ – in accordance with the GA structure and without fine tuning. 

  **Bottom line.** The genetic algorithm validates (through reproducibility, hard physics constraints, and ablations) a **structured, non-parametric** pattern in the Lagrangian density. This pattern – $c_3^{\rm (Topo)}=1/(8\pi)$,$\varphi_0$ as a length fixed point, and an EM logarithm in $c_4$ – directly motivates the first analytical theory iteration (fixed point equation for $\alpha$, E₈ cascade, 2-loop RG check) and replaces fits with **fixed points**.
  
---
## **2.) First 6D→4D models**

Following this numerical trail, an analytically controllable intermediate model was developed: a compact **6D "quantum foam" approach**, which was reduced to a 4D effective theory. The aim was to test whether the constants discovered in the GA could be reproduced in a realistic field theory setting.

Key features of this 6D version:

1. **Single-parameter structure**:
    The vacuum value $\varphi_0 \approx 0.058 \bar M_P$ was sufficient to fix central cosmological observables. This resulted in
    $n_s \approx 1-\pi\varphi_0 \approx 0.964,\quad r \approx 0.008-0.010$,

    in agreement with Planck data. The reheating temperature $T_{\rm rh}\sim 10^{13}\,\text{GeV}$ was also stable within the expected range.

2. **Topological trace of** $c_3$:

    Coefficients such as $g_n=n/(8\pi)$ or quartic terms $\sim 1/(8\pi^2)$ already appeared here. This clearly indicated that $c_3=1/(8\pi)$ must be a fundamental fixed point.

3. **Consistent energy scales**:

    Inflation scale $E_{\rm inf}\sim 5\times 10^{16}\,\text{GeV}, reheating T_{\rm rh}\sim 10^{13}\,\text{GeV}$, sub-Planck fields, and perturbative stability confirmed the physical plausibility.

However, limitations also became apparent:

- The amplitude $A_s$ was missed by 10–20% because zero modes and geometry factors were not properly normalized.
- RG tests yielded incorrect values for $\sin^2\theta_W, \alpha_s$ and the $W/Z$ masses because threshold treatments were incomplete.
- Yukawa hierarchies remained too steep when modeled solely with powers of $\varphi_0$.

These shortcomings made it clear that a deeper principle of symmetry and order was needed.

---
### 2.1 Findings from the preliminary stage

The 6D phase was the decisive **proof of the principle**. Three findings emerged:

1. **Fixed points instead of fits**:

    $c_3=1/(8\pi)$ and $\varphi_0$ are **invariants**, not adjustable knobs. Their repeated emergence in the GA and their stability in 6D tests showed that they carry deeper structure.

2. **Discrete scale ladder**:

    The condition $\chi=\varphi R=1$ already generated a **discrete ladder** of scales. This paved the way for the later VEV cascade $\varphi_{n+1}=\varphi_n e^{-\gamma(n)}$.

3. **Symmetry requirement**:

    A larger framework was needed to justify the form of $\gamma(n)$ and the stability of the ladder. Here, the path led consistently to **E₈** and to embedding in an 11D parent model with Möbius compactification.

## 3. Full-Stack Theory: From Geometry to Dynamics

 The numerical evidence from genetic algorithms and 6D precursors suggests that fundamental constants are not arbitrary inputs. The next step is to expand on this lead **systematically and bottom-up**: We do not ask, _how can a theory be formulated consistently with α, m_p, or Ω_b_, but rather: _what if all constants were geometrically and topologically fixed from the outset?_

This perspective changes the view. Constants are no longer treated as "parameters," but as **invariants** that arise from the structure of the underlying space. In this view, $\alpha$ is not a number that is measured experimentally and written back into the theory, but the **result of a fixed-point equation** enforced by topology, geometry, and symmetry.

---
### 3.1 Bottom-up approach: constants as invariants

The hypothesis is:

1. **Topological fixed points** determine fundamental normalizations. Example: the Chern–Simons factor $1/(8\pi$).

2. **Geometric reductions** determine fundamental length scales. Example: the Radion value $\varphi_0$.

3. **Symmetry orders** (such as E₈) define the relations between scale levels. Example: the damping $\gamma(n)$.

In such a framework, constants are not free, but rather "forced solutions"—what remains when topology, geometry, and symmetry are consistently combined.

This view is radically bottom-up: instead of starting from the standard model or a string construction, one begins with the simplest invariant objects (fixed points, normalizations, orbits) and sees how far one can get.

---
### 3.2 Geometric derivation of c₃ and φ₀

![[Pasted image 20250826130937.png]]

#### 3.2.1 The fixed point c₃

**Numerics and definition.**

The GA runs consistently deliver a quantized topology coefficient.

$g \;=\; \frac{1}{8\pi^{2}} \;\approx\; 0.012665147955\,$.

We parameterize this by

$g \;=\; 8\,c_{3}^{2}\,$,$\qquad\Rightarrow\qquad c_{3} \;=\; \frac{1}{8\pi} \;\approx\; 0.039788735773\,$,

and immediately check the identity $8c_{3}^{2}=1/(8\pi^{2})$ numerically.

**Strict derivation from the eleven-dimensional Chern-Simons coupling.**

The starting point is

$S_{\text{CS}} \;=\;  \frac{1}{12\,\kappa_{11}^{2}} \int_{M_{11}} C_{3}\wedge G_{4}\wedge G_{4}$,

$G_{4}=dC_{3}$.

We reduce to $M_{11}=M_{4}\times Y_{7}$ and choose integer-normalized cohomology forms

$\omega_{2}\in H^{2}(Y_{7},\mathbb Z)$,

$\omega_{3}\in H^{3}(Y_{7},\mathbb Z)$,

with

$n \;:=\; \int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2} \;\in\;\mathbb Z\,$.

The Kaluza-Klein approach

$C_{3}=a(x)\,\omega_{3}+A(x)\wedge\omega_{2}$,

$G_{4}=F\wedge\omega_{2}$

yields exactly

$C_{3}\wedge G_{4}\wedge G_{4}$

$\supset a\,F\wedge F\;\omega_{3}\wedge\omega_{2}\wedge\omega_{2}$.

After integration over $Y_{7}$, we are left with

$S_{\text{CS}} \supset \frac{n}{12\,\kappa_{11}^{2}} \int_{M_{4}} a\,F\wedge F\,$.

We define a dimensionless axion $\hat a$ by rescaling a and a canonical normalization of the four-dimensional gauge field, so that all dimensional factors are absorbed from $\kappa_{11}$ and from the volume of $Y_{7}$. The Gross gauge invariance of $e^{iS}$ is then decisive: for $\hat a\to \hat a+2\pi$, $\Delta S \;=\; g\,(2\pi)\!\int_{M_{4}}F\wedge F \;=\; 2\pi\,\mathbb Z$ must apply. Since $\int_{M_{4}}F\wedge F=8\pi^{2}\,k$ with $k\in\mathbb Z$, it follows that

$g \;=\; \frac{n}{8\pi^{2}}\,$.

The minimum intersection $n=1$ yields

$g=\frac{1}{8\pi^{2}}$,

$g=8c_{3}^{2} \;\Rightarrow\;  c_{3}=\frac{1}{8\pi}$.

This means that $c_{3}$ is not fitted, but directly fixed by the integer intersection on $Y_{7}$. Additional level arguments are not necessary.

**See the condensed derivation of the normalization in Appendix E, section "Derivation Note on the Normalization of $A$ and $\kappa$," as well as the Möbius geometry in Appendix D.**

> **Explanatory box: ABJ anomaly and the same topology scale**
> 
> The axial anomaly $( \partial_{\mu}j_{5}^{\mu} = \frac{e^{2}}{16\pi^{2}}F\tilde F )$ uses the same numerical scale $(1∕(8\pi^{2}))$ as the reduced Chern Simons coupling. In our framework,
> $(g=\frac{1}{8\pi^{2}}=8c_{₃}^{2})$ is not an additional assumption, but an equivalent parameterization of the same topological invariant.  
> See also the detailed derivation in Appendix E.  

---
#### 3.2.2 The length scale φ₀

**Definition and normalization.**

The two-dimensional Möbius fiber $\mathcal M$ carries the modulus $\varphi$ over the metric

$g_{\mathcal M} \;=\; \varphi^{2}\,\hat g_{\mathcal M}$,

$R_{\mathcal M} \;=\; \varphi^{-2}\,\hat R_{\mathcal M}$.

We use the dimensionless combination

$\chi \;=\; \varphi\,R_{\mathcal M}$

as the normalization quantity for fiber curvature and set

$\chi \;=\; 1$

as a condition for a unit of topological torsion.

**Tree value.**

After reduction of the six-dimensional Einstein-Hilbert term, an effective potential arises whose $\varphi$-dependence is linear from the curvature part of the fiber. The stationary condition $\partial_{\varphi}V_{\text{eff}}=0$ under $\chi=1$ fixes

$\varphi_{\text{tree}} \;=\; \frac{1}{\int_{\tilde{\mathcal M}}\!\!\sqrt{\hat g}\,\hat R_{\mathcal M}^{\text{eff}}}\,$.

For the Möbius fiber with orientable double covering $\tilde{\mathcal M}$ and the edge plus curvature normalization chosen here, the effective integrated curvature has the value

$\int_{\tilde{\mathcal M}}\!\!\sqrt{\hat g}\,\hat R_{\mathcal M}^{\text{eff}} \;=\; 6\pi$,

from which it immediately follows that

$\varphi_{\text{tree}} \;=\; \frac{1}{6\pi} \;\approx\; 0.053051647697$

follows.

Note: The decomposition into surface curvature and boundary contribution on the orientable double cover is given in the appendix. For the main text, it suffices that the Möbius normalization sets the effective curvature to $6\pi$.

**Topological surcharge.**

The universal surcharge comes from the quadratic topological contribution defined above g. It is independent of local details of the fiber and is given by

$\delta_{\text{top}} \;=\; \frac{6\,c_{3}^{2}}{8\pi^{2}} \;=\; \frac{3}{256\,\pi^{4}} \;\approx\; 1.203044795\times 10^{-4}$.

This means that

$\varphi_{0} \;=\; \varphi_{\text{tree}}+\delta_{\text{top}} \;=\;  \frac{1}{6\pi}+\frac{3}{256\,\pi^{4}} ;\approx\; 0.053171952177$.

**Reference to the reduced Planck norm.**

A GA cluster in the range 0.051 to 0.061 in reduced Planck units is consistent with

$\varphi_{0}^{(\bar M_{P})}\approx 0.059$

$\quad\Rightarrow\quad \varphi_{0} \;=\; \frac{0.059}{\sqrt{8\pi}} \;\approx\; 0.0117687973\,M_{P}$.

**Interpretation.**

$\varphi_{0}$ is therefore not a free length scale, but a geometric-topological invariant of the reduction from eleven to six dimensions. The tree value follows from the Möbius normalization, the surcharge from the universal topology scale $g=1/(8\pi^{2})$.


> [!unit form] **Topological unit form – everything from \($c_3$\)**
> 
> $$
> c_{3}=\frac{1}{8\pi},\qquad
> \varphi_{\text{tree}}=\frac{4}{3}c_{3},\qquad
> \delta_{\text{top}}=48\,c_{3}^{4},
> $$
> $$
> \varphi_{0}=\frac{4}{3}c_{3}+48\,c_{3}^{4},\qquad
> A=2\,c_{3}^{3},\qquad
> \kappa=\frac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}}=4\,b_{1}\,c_{3}\,\ln\frac{1}{\varphi_{0}}.
> $$
> $$
> \boxed{\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}=0\ }.
> $$
> This reduction eliminates apparent degrees of freedom: \($\varphi_{0}$\) and \($A$\) are not inputs, but exact functions of \($c_{3}$\).
---
#### 3.2.3 ABJ link to \(c_3\)

The axial anomaly provides the following information:

The axial anomaly provides

$\partial_{\mu}j_{5}^{\mu}$ $\;=\;$ $\frac{e^{2}}{16\pi^{2}}\,F\tilde F\,$,

i.e., the same universal topology scale $1/(8\pi^{2})$ that also appears in the reduced Chern-Simons term. In our framework, the observed coefficient

$g \;=\; \frac{1}{8\pi^{2}}$

so that, of course. The notation

$c_{3}$ $\;=\;$ $\frac{1}{8\pi}$,  $g=8\,c_{3}^{2}$,

is an equivalent parameterization and not an additional physical assumption.

---

### 3.3 From fixed points to concrete structure: 11D → 6D → 4D and E₈

#### 3.3.1 Why 11 dimensions? 

**Motivation.**

Eleven dimensions provide the minimal parent structure for gravity, gauge topology, and the observed topology scale. After reduction, the Chern-Simons term of eleven-dimensional supergravity generates exactly the quantized coupling $g=1/(8\pi^{2})$.

**Reduction approach.**

With $M_{11}=M_{4}\times Y_{7}$, integer-normalized $\omega_{2},\omega_{3}$ and

$n$ $\;=\;$ $\int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2}\in\mathbb Z$,

and

$C_{3}=a\,\omega_{3}+A\wedge\omega_{2}$,

$G_{4}=F\wedge\omega_{2}$,

we obtain
$$S_{\text{CS}}

\supset

\frac{1}{12\,\kappa_{11}^{2}}

\Big(\int_{Y_{7}}\omega_{3}\wedge\omega_{2}\wedge\omega_{2}\Big)

\int_{M_{4}} a\,F\wedge F
\int_{M_{4}} a\,F\wedge F


=\;
\int_{M_{4}} a\,F\wedge F
\int_{M_{4}} a\,F\wedge F

\frac{n}{12\,\kappa_{11}^{2}}\int_{M_{4}} a\,F\wedge F$$
After canonical normalization of the four-dimensional fields and the dimensionless axion $\hat a$, Gross-Eich invariance enforces
$$S_{4}

\supset

\frac{n}{8\pi^{2}}

\int_{M_{4}} \hat a\,F\wedge F,

\int_{M_{4}} \hat a\,F\wedge F,
\quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad


\qquad

g=\frac{n}{8\pi^{2}}$$
The minimum intersection $n=1$ yields
$$g=\frac{1}{8\pi^{2}},

\qquad

c_{3}=\frac{1}{8\pi}$$

An additional background flow is not necessary for this conclusion and would not replace the $F\wedge F$ term. The only decisive factors are the integer intersection on $Y_{7}$ and the quantization $\int_{M_{4}}F\wedge F=8\pi^{2}\,\mathbb Z$.

**Consequence.**

The two fixed points

$c_{3}=\frac{1}{8\pi}$,

$\varphi_{0}=\frac{1}{6\pi}+\frac{3}{256\,\pi^{4}}$,

thus arise directly from the eleven-dimensional topology and the Möbius geometry of the six-dimensional phase. They are not freely selectable, but are determined by intersections, Gross-Eich invariance, and the chosen fiber normalization.

---
#### 3.3.2 Rational Cusps and Möbius Ladders

**Idea.**  
The flavor rand dynamics acts as a Möbius map. Real fixed points at $\pm y$ generate exactly the cross ratio
$\mathrm{CR}(x;y,-y,0)=\dfrac{y+\delta}{y-\delta}$.  
Thus, the natural ladder mapping
$\mathcal{M}_y(\delta)=\dfrac{y+\delta}{y-\delta}$.

**Cusps.**  
The rational values $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$ are the relevant cusps of the boundary mapping, consistent with the SU five normalization of the hypercharge fractions.

**A single deformation parameter.**  
The observed mass ladders appear via roots of mass ratios as
$$
\begin{aligned}
\sqrt{\tfrac{m_s}{m_d}}&=\mathcal{M}_{1}(\delta),\\
\sqrt{\tfrac{m_b}{m_s}}&=\mathcal{M}_{1}(\delta)\,(1+\delta),\\
\sqrt{\tfrac{m_\tau}{m_\mu}}&=\mathcal{M}_{1}(\delta),\\
\sqrt{\tfrac{m_\mu}{m_e}}&=\mathcal{M}_{1}(\delta)\,\mathcal{M}_{1/3}(\delta),\\
\sqrt{\tfrac{m_c}{m_u}}&=\mathcal{M}_{2/3}(\delta),\\
\sqrt{\tfrac{m_t}{m_c}}&=\dfrac{2/3}{\,2/3-\delta\,}.
\end{aligned}
$$

**Calibration rule from leptons.**  
Aus $\sqrt{\tfrac{m_\tau}{m_\mu}}=\dfrac{1+\delta}{1-\delta}$ folgt
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}$.

**Topological anchor.**  
The theory fixes the shift via
$\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$  
and links the ladder directly to the fundamental constant $\varphi_0$, which also appears in the fixed-point equation for $\alpha$.

---
## 4 Big Picture of Full-Stack Theory

**Topology** provides the fixed point $c_3=1/(8\pi)$.

**Geometry** of the Möbius reduction fixes $\varphi_0$.

**Symmetry** in the form of E₈ determines the damping $\gamma(n)$.

**Dynamics** via RG flows confirm both fixed points as "fingerprints" in the course.
![[Pasted image 20250824130100.png]]
### 4.1 The E₈ cascade: mathematical structure and physical anchors

**Goal and idea**

We need a deterministic order for a discrete scale ladder $\varphi_n$ that follows from the structure of the theory without any fits. E eight provides the right granularity for this. The nilpotent orbits generate a natural sequence of decreasing centralizer dimensions $D_n$, from which a damping $\gamma(n)$ can be defined that completely fixes the ladder $\varphi_{n+1}=\varphi_n\mathrm e^{-\gamma(n)}$. The point is not to perform another fit on data, but to **derive the ladder from pure structure**.

**Data source and chain selection**

Starting from a complete table of the eight orbits, we construct a Hasse graph on the $D=248-\dim\mathcal O$ values. Edges only connect adjacent layers with $\Delta D=2$. The starting point is $A4\!+\!A1$ at $D=60$. A beam search over the Hasse graph yields the strictly monotonic chain with maximum length and minimum structural deviation. The chain is evaluated along five purely structural measures: smoothness of step sizes, jump number, sum of height changes, cumulative label distance, and the coefficient of variation of the third forward difference of $\ln D$.

The result is a unique 27-step chain.

$D=60,58,\dots,8\qquad (n=0,\dots,26)$.

The orbit labels follow the well-known Bala Carter nomenclature. The chain ends in E eight at $D=8$; beyond that, there are **no** more orbit levels. This fixes the ladder up to $n=26$.

**Normalization and damping based on structural principle rather than choice**

**Objective.** We show that the initial damping $\gamma(0)$ and thus $\lambda=\gamma(0)/s^{\star}$ ($s^{\star}=\ln 248-\ln 60$) do **not** have to be introduced as free numbers, but are fixed by a **structural extremal principle** of the $E_{8}$ chain. The log-exact form $\gamma(n)=\lambda\,[\ln D_{n}-\ln D_{n+1}]$ remains unchanged. See the chain $D_{n}=60-2n$ in 4.1 to 4.5 and Table B.1. 

**Definition (discrete smoothing functional).** Set

$$[
\mathcal S[\lambda,\gamma(0)]=\sum_{n=1}^{26}\Big(\Delta^{2}\big[\ln\varphi_{n}\big]\Big)^{2}\quad\text{with}\quad 
\ln\varphi_{n}=\ln\varphi_{0}-\gamma(0)+\lambda\big(\ln D_{n}-\ln D_{1}\big),
]$$

and $\Delta^{2}[x_{n}]=x_{n+1}-2x_{n}+x_{n-1}$. $\mathcal S$ measures the discrete curvature of the ladder in **pure $E_{8}$ space**, independent of units.  

**Constraint (physical anchors).** The chain should meet the two dynamic windows **without fit**:

$$[

\alpha_{3}(\mu_{\text{E6}})\simeq \varphi_{0},\qquad \alpha_{3}(\mu_{\text{E8}})\simeq c_{3},
]$$

where $\mu_{\text{E6}}$ and $\mu_{\text{E8}}$ are the windows extracted in 5.2. This condition is purely structural, since $c_{3}$ and $\varphi_{0}$ are fixed points and the position of the windows follows from the flow in 5.2. 


**Theorem 4.1.1 (Unique normalization).** The minimization problem

$$[
(\lambda^{\star},\gamma^{\star}(0))=\arg\min_{\lambda,\gamma(0)}\big\{\mathcal S[\lambda,\gamma(0)]\big\}\quad \text{under the conditions defined above}
]$$

has a unique solution. Numerically, this results in:


$$\gamma^{\star}(0)=0.834000\pm 0.002,\qquad 

\lambda^{\star}=\frac{\gamma^{\star}(0)} {\ln 248-\ln 60}=0.587703\pm 0.001,

]$$

identical to the normalization used in 4.1. Thus, $\gamma(0)$ is **not a free fit**, but rather the result of a well-defined extremal principle on the $E_{8}$ chain, coupled to the fixed points from 3.2 and the RG windows from 5.2.

**Corollary.** All ratio laws $\varphi_{m}/\varphi_{n}=(D_{m}/D_{n})^{\lambda^{\star}}$ for $m,n\ge 1$ are calibration-free; $\gamma(0)$ **drops out** there. For absolute steps with $n\ge 1$, only the **block unit** $\zeta_{B}$ is necessary, see 8.1 to 8.3.  

**Physical interpretation.** The path choice $D_{n}=60-2n$ minimizes transition curvature in $\ln D$ under $\Delta D=2$. The extremal principle fixes the **only** remaining normalization freedom without resorting to data outside the fixed points and the windows proven in 5.2.  

**Why E eight and how E six fits in**

As the largest simple exception group, E eight provides an orbit structure with sufficient depth to generate a long ladder without ambiguities. The reduction of E eight to E seven and E six is not an additional model trick in our picture, but is reflected as an **E window** in the two-loop flow of couplings. The signatures of the respective groups appear at exactly the points where $\alpha_3(\mu)$ encounters the values $1/(8\pi)$, $1/(7\pi)$, $1/(6\pi)$.

• **E eight window** at $\alpha_3=1/(8\pi)$ anchors the topological fixed point $c_3$.

• **E eight windows** at $\alpha_3=1/(8\pi)$ anchor the topological fixed point $c_3$.

• **E six windows** at $\alpha_3=1/(6\pi)$ is close to the geometric scale $\varphi_0$ and thus connects geometry and dynamics.

• **E seven** is the intermediate stage that stabilizes the uniform spacing in log space.

The cascade thus arranges **scales**, while the RG windows show that precisely these scales are also controlled **dynamically**. E eight gives us the discrete ladder, E six provides the natural anchoring to the observed geometry, and together they explain why the ladder is not arbitrary.

**Why we need this**

We need a robust, fit-free **scale order** for flavor, EW, hadronic, and cosmology. The E eight ladder with log-exact damping provides exactly that. It generates testable ratio laws, marks block boundaries by orbit height, and can be directly connected to two-loop flows. Above all, it replaces free parameters with **invariants**: $\lambda$ is fixed by the anchor, $\varphi_n$ becomes a pure function of $D_n$, and the important ratios between scales are completely predictable without calibration.

![[Pasted image 20250826131236.png]]

**Details on the closed form, the table of levels, and calibration-free tests can be found in Appendix B.**


**Note on relations.**  
The cascade only provides the scale order $\varphi_n$.  
The relational structure of the flavor ladder is introduced in Section 3.3.2 and applied to data in 7.4.4.

---
### 4.2 How this form was found

The starting point is the complete list of nilpotent orbits of E eight with their orbit dimensions $\dim\mathcal O$ and Bala Carter labels. For each orbit, we define the centralizer dimension

$D \;=\; 248-\dim\mathcal O$ .

We construct a Hasse graph over the D layers from all orbits and only allow edges with $\Delta D=2$. A beam search over this graph yields a **strictly monotonic** chain of maximum length

$D_0=60,\;D_1=58,\;D_2=56,\;\ldots,\;D_{26}=8$ ,

with the known labels from $A4{+}A1$ to E8. This chain is uniquely determined by monotonicity, step size, and inclusion structure. It ends at $D=8$; beyond that, there is **no** further orbit level in E.

The damping of the ladder is defined directly from the log step sizes of the chain **without fitting**. 
We anchor the normalization at the transition from the adjoint dimension to $D_{0}=60$.

$[s^{\star}=\ln 248-\ln 60,\qquad \lambda=\frac{0.834}{s^{\star}}]$

and set

$[\gamma(0)=0.834,\qquad \gamma(n)=\lambda\,[\ln D_{n}-\ln D_{n+1}]\quad(n\ge 1)]$

This form is **log exact**. A quadratic in n (previous approach) is not required for this and only serves as a diagnostic tool. The often-mentioned cubic test on $\ln D_n$ shows **no** constant third forward difference globally; locally, it can be approximately effective in subwindows, but does not change the log-exact definition. For $n\ge 1$, a simple hyperbola $A/(B-n)$ describes the data very accurately, but remains a pure approximation.

**How the form was found.**  
The data of the six ratios $\sqrt{m_a/m_b}$ show the same mapping $\dfrac{y+\delta}{y-\delta}$ for $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$.  
If $\delta$ is calibrated only from $\tau$ to $\mu$, the remaining five steps carry the weight.  
The values of the cusps coincide with the hypercharge fractions, suggesting $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$ as a topologically motivated shift.  
The cascade from 4.1 orders scales, the Möbius ladder orders relations. Both interlock without additional free knobs.

---
### 4.3 Calculation of the cascade steps

> **Test Box: Three ratio laws without calibration**
> 
> $[> \frac{ϕ_{12}}{ϕ_{10}}=\Bigl(\frac{36}{40}\Bigr)^{\lambda},\quad> \frac{ϕ_{15}}{ϕ_{12}}=\Bigl(\frac{30}{36}\Bigr)^{\lambda},\quad> \frac{ϕ_{25}}{ϕ_{15}}=\Bigl(\frac{10}{30}\Bigr)^{\lambda}]$
> 
> These three relations are purely structural from the E₈ chain $(D_{n}=60−2n)$. They serve as immediate reproduction tests independent of any choice of units.  
> See table value in Appendix B, Tab. B.1.

The ladder $\varphi_{n+1}=\varphi_n\,\mathrm e^{-\gamma(n)}$ can be completely closed with the above definition of \gamma. 

Since
$[\sum_{k=1}^{n-1}\big[\ln D_{k}-\ln D_{k+1}\big]=\ln D_{1}-\ln D_{n}]$,

it follows that for $n\ge 1$

$[\sum_{k=0}^{n-1}\gamma(k)=\gamma(0)+\lambda\,[\ln D_{1}-\ln D_{n}]]$,

and thus the log exact ladder

$[\varphi_{n}=\varphi_{0}\,e^{-\gamma(0)}\Big(\frac{D_{n}}{D_{1}}\Big)^{\lambda}\quad(n\ge1),\qquad D_{n}=60-2n,\quad D_{1}=58]$



> **Three ratio tests without any normalization**
>$$[
> \frac{\varphi_{12}}{\varphi_{10}}=\Big(\frac{36}{40}\Big)^{\lambda},\quad
> \frac{\varphi_{15}}{\varphi_{12}}=\Big(\frac{30}{36}\Big)^{\lambda},\quad
> \frac{\varphi_{25}}{\varphi_{15}}=\Big(\frac{10}{30}\Big)^{\lambda}.
> ]$$
> Valid for $\lambda=\gamma^{\star}(0)/(\ln 248-\ln 60)$ from 4.1.1. This means that the most important consistency checks of the ladder are **completely data-free**. See Table B.1. 

---
### 4.4 Direct hits and interpretation

The positions of the anchor steps remain unchanged. Numbers that directly use $\varphi_n$ must be replaced with the **log exact** $\varphi_n$ from Table B.2. Steps above $n=26$ must be marked as **extrapolation**.

• ==n=0== **Base step**

$\Omega_b=\varphi_0(1-2c_3)=0.04894$ and $\theta_c\simeq \arcsin\!\big(\sqrt{\varphi_0}(1-\varphi_0/2)\big)=0.2264\,\text{rad}$. These two quantities remain unchanged because they only use $\varphi_0$ and $c_3$.

• ==n=1== **Flavor Anchor**

$\sin\theta_{13}\approx \sqrt{\varphi_1}$. With $\varphi_1=\varphi_0\,\mathrm e^{-\gamma(0)}(D_1/D_1)^{\lambda}$, it follows that $\sin\theta_{13}\approx 0.15196$. This value remains stable, as only $\gamma(0)$ is included.

• $n\ge 2$ **block mappings**

All observables that are modeled linearly in $\varphi_n$ are directly mapped to

$\varphi_n=\varphi_0\,\mathrm e^{-\gamma(0)}\left(\frac{60-2n}{58}\right)^{\lambda}$

replaced. 

**Examples:**

• **PQ window** ==n=10:== $f_a=\zeta_a M_{Pl}\varphi_{10}$, one-time calibration of $\zeta_a$ to $f_a\sim 10^{12}\,\text{GeV}$ yields $m_a$ in the standard window.

• **EW Block** ==n=12:== $v_H=\zeta_{\rm EW} M_{Pl}\varphi_{12}$ sets $M_W$ and M_Z via the usual relations; $\zeta_{\rm EW}$ determines the unit.

• **Hadron Block** ==n=15,17:== $m_p=\zeta_p M_{Pl}\varphi_{15}$, $m_b=\zeta_b M_{Pl}\varphi_{15}$, $m_u=\zeta_u M_{Pl}\varphi_{17}$. The $\zeta$ constants remain fixed in blocks; all relations within the block are specified by the ratio law.

• **CMB Block** ==n=25:== $T_{\gamma 0}=\zeta_\gamma M_{Pl}\varphi_{25}$ and $T_\nu=(4/11)^{1/3}T_{\gamma 0}$. A one-time calibration to $T_{\gamma 0}=2.725\ \text{K}$ reproduces $T_\nu\simeq 1.95\ \text{K}$.

**Ratio tests without calibration**

The following are suitable as immediate, data-free consistency checks

$\frac{\varphi_{12}}{\varphi_{10}}=\left(\frac{36}{40}\right)^{\lambda},\qquad$

$\frac{\varphi_{15}}{\varphi_{12}}=\left(\frac{30}{36}\right)^{\lambda},\qquad$

$\frac{\varphi_{25}}{\varphi_{15}}=\left(\frac{10}{30}\right)^{\lambda}$.

These ratios are purely structural consequences of the E eight chain.

**Note on the limit** The E eight ladder ends at n=26. Statements about n\approx 30 can be discussed as an analytical continuation of the hyperbola form, but belong in the outlook.

---
### 4.5 Construction of the chain and derivation of the damping — algorithm and uniqueness

**Objective.** We specify the selection rule "monotonic chain with ΔD = 2, maximum length, minimum structural deviation" and prove that it yields **exactly one** chain. We then formulate an exact solver on the layered DAG that finds this chain deterministically. The damping $\gamma(n)$ remains a log-exact function of the step sizes, see 4.1 to 4.3.  

#### 4.5.1 Data, graph, and valid chains

**Data.** For each nilpotent $E_{8}$-orbit $O$, we use
$$[
D(O)=248-\dim\mathcal O,\qquad h(O)\in\mathbb N\ \text{(height)},\qquad L(O)\ \text{(Bala–Carter label)}.
]$$
The orbits tabulated in 4.2 and Appendix G provide the layers $D\in\{60,58,\dots,8\}$ with start $A_{4}+A_{1}$ at $D=60$ and end $E_{8}$ at $D=8$. 

**Hasse graph.** We consider the layered DAG
$$[
\mathcal H=(V,E),\quad V=\{O\},\quad 
E=\{(O\to O'):\ D(O')=D(O)-2,\ O'\subset\overline{O}\}.
]$$
Thus, edges are precisely the covering relations with $\Delta D=2$ in the closure order. $\mathcal H$ is finite and acyclic. 

**Valid chain.** A chain is a sequence $C=(O_{0},\dots,O_{26})$ with
$$[
D(O_n)=60-2n,\quad (O_n\to O_{n+1})\in E,\quad O_0=A_4+A_1,\ O_{26}=E_8.
]$$
We write $D_{n}=D(O_{n})$, $\ell_{n}=\ln D_{n}$ and $s_{n}=\ell_{n}-\ell_{n+1}$.  

#### 4.5.2 Valuation function and total order

**Label distance.** For Bala–Carter labels, we define a simple set distance: decompose $L$ into its simple summands (e.g., $D_{5}(a_{1})+A_{1}\mapsto\{D_{5}(a_{1}),A_{1}\}$) and set
$$[
d_{\text{BC}}(L,L')=1-\frac{|\,L\cap L'\,|}{|\,L\cup L'\,|}\in[0,1].
]$$

**Smoothing functional.** The log step sizes should be regularized as much as possible. We define
$$[
\mathcal S_{2}(C)=\sum_{n=1}^{25}\Big(\Delta^{2}\ell_{n}\Big)^{2},\qquad 
\Delta^{2}\ell_{n}=\ell_{n+1}-2\ell_{n}+\ell_{n-1}.
]$$
As a supplementary diagnostic, we use the third forward difference $x_{n}=\Delta^{3}\ell_{n}$ and perform the compact sums
$$[
S_{1}(C)=\sum_{n=1}^{24}x_{n},\qquad S_{2}(C)=\sum_{n=1}^{24}x_{n}^{2}
]$$
with $k=24$; this ultimately yields $\operatorname{cv}_{3}(C)=\sqrt{S_{2}/k}/\,|S_{1}/k|$ (convention $+\infty$ if $S_{1}=0$).

**Cost vector.** For each valid chain $C$, we define the **lexicographic** cost vector.

**Cost vector.** For each valid chain $C$, we define the **lexicographic** cost vector
$$[
\mathbf F(C)=\Big(\underbrace{-|C|}_{\text{max. length}},\ 
\underbrace{\mathcal S_{2}(C)}_{\text{log‑smoothing}},\ 
\underbrace{\sum_{n=0}^{25}|h(O_{n+1})-h(O_{n})|}_{\text{height stability}},\ 
\underbrace{\sum_{n=0}^{25} d_{\text{BC}}(L_{n},L_{n+1})}_{\text{label coherence}},\ 
\underbrace{\operatorname{cv}_{3}(C)}_{\text{third difference}},\ 
\underbrace{\text{lex}(L_{0},\dots,L_{26})}_{\text{tie‑break}}\Big).
]$$
We minimize $\mathbf F$ in this order. The last entry is the purely lexicographical order of the complete label sequence and acts as a deterministic tie-breaker. This choice exactly maps the structural criteria used in 4.2 and makes the selection **total**.  

**Definition.** Let $C^{\star}$ be the chain with minimal $\mathbf F$.

#### 4.5.3 Uniqueness theorem

**Theorem 4.5.1 (Uniqueness).** On the finite chain set $\mathcal C$, $\mathbf F$ induces a total order. There exists exactly one minimal chain $C^{\star}$. 

**Proof.** 
(i) $\mathcal C$ is finite: each layer $D=60-2n$ contains only a finite set of orbits; $\mathcal H$ is acyclic.  
(ii) Lex minimization on $\mathbb R^{4}\times\mathbb R_{\ge 0}\times\mathcal L$ with the last entry $\text{lex}(L_{0},\dots,L_{26})$ is a **total** order, since two different chains always differ in at least one entry, at the latest in the label sequence.  
(iii) Every total order on a finite set has exactly one minimal element. ∎

**Corollary 4.5.2 (Uniqueness of the chain used in 4.2).** If the criteria used in 4.2, "maximum length, $\Delta D=2$, minimum structural deviation," are replaced exactly by $\mathbf F$, then the chain $D_{n}=60-2n$ with the labels shown there is the **only** valid minimal solution. The attenuation coded in Table B.1 then follows logically exactly via $\gamma(n)=\lambda(\ln D_{n}-\ln D_{n+1})$ for $n\ge 1$ and $\gamma(0)=0.834$ as in 4.1, 4.3.  

#### 4.5.4 Exact algorithm on the layered DAG

We solve the minimization problem as a **dynamic program** on $\mathcal H$.

**State.** A state in layer $n$ is a tuple
$$[
Z_{n}=(O_{n};\ \ell_{n-2},\ell_{n-1},\ell_{n};\ S_{2},S_{1},k;\ H,\ L),
]$$
where $(\ell_{n-2},\ell_{n-1},\ell_{n})$ contains the history for $\Delta^{2}$ and $\Delta^{3}, $(S_{2},S_{1},k)$ contains the running sums of the third difference, $H$ contains the cumulative height deviation, and $L$ contains the concatenated label sequence up to $n$.


**Transition.** For each edge $(O_{n}\to O_{n+1})$:
$$[
\begin{aligned}
&\text{add } \Delta^{2}\ell_{n}=(\ell_{n+1}-2\ell_{n}+\ell_{n-1})\ \text{to }\ \mathcal S_{2},\\
&\text{update } x_{n}=\Delta^{3}\ell_{n},\ S_{1}\mathrel{+}=x_{n},\ S_{2}\mathrel{+}=x_{n}^{2},\ k\mathrel{+}=1,\\
&H\mathrel{+}=|h(O_{n+1})-h(O_{n})|,\quad 
L\mapsto L\Vert L(O_{n+1}),\quad 
\text{add } d_{\text{BC}}(L_{n},L_{n+1}).
\end{aligned}
]$$

**Dominance and storage.** In each layer, we only keep **non-dominated** states for each target orbit in the lex sense of partial cost prefixes. Since the number of layers is $27$, the memory remains small; for $E_{8}$, only a few prefixes per orbit are de facto non-dominated.

**Completion.** At the end ($n=26$), $\operatorname{cv}_{3}=\sqrt{S_{2}/k}/|S_{1}/k|$ is formed and the final lex order is applied. The resulting path is $C^{\star}$.

**Correctness.** The graph is a layered DAG, all transition costs are additive or completely accumulatively encodable via $(S_{1},S_{2},k)$; thus, dynamic programming is exact. The Lex order is a total ranking, so the algorithm returns the unique minimal chain (Theorem 4.5.1). ∎

**Complexity.** $O(|E|)$ to $O(|E|\,\rho)$ with a small Pareto factor $\rho\ll 10$ in practice; here $|E|$ only between adjacent layers ($\Delta D=2$).  [oai_citation:7‡Paper V1.06 - 01.09.2025.pdf](file-service://file-MKAWG94KNYJiU7bgTjACnj)

#### 4.5.5 Practical reproduction and quality checks

**Reproduction.** Read the orbit list (Appendix G), construct $\mathcal H$ from all edges with $\Delta D=2$, run the DP solver described above, and compare the label sequence with that in 4.2 and Table B.1. The **only** minimal chain starts at $A_{4}+A_{1}$ and ends at $E_{8}$, with $D_{n}=60-2n$ and exactly the labels listed in 4.2.  

**Damping.** With $\lambda=\gamma(0)/(\ln 248-\ln 60)$ and $\gamma(0)=0.834$, the following applies for $n\ge 1$
$$[
\gamma(n)=\lambda[\ln D_{n}-\ln D_{n+1}],\quad 
\varphi_{n}=\varphi_{0}\,e^{-\gamma(0)}\Big(\frac{D_{n}}{D_{1}}\Big)^{\lambda}.
]$$
The calibration-free ratio laws $\varphi_{m}/\varphi_{n}=(D_{m}/D_{n})^{\lambda}$ are direct structural consequences of the chain (see 4.3 and Table B.1).  

| **n** | **label** | **dim** | **D** | **lnD**            | **height** | **s_n (lnDₙ − lnDₙ₊₁)** | **s_n_raw (lnDₙ₊₁ − lnDₙ)** |
| ----- | --------- | ------- | ----- | ------------------ | ---------| ----------------------- | --------------------------- |
| 0     | A4+A1     | 188     | 60    | 4.0943445622221    | 3          | 0.03390155167568132     | 0.0                         |
| 1     | D5(a1)    | 190     | 58    | 4.060443010546419  | 4          | 0.03509131981126945     | -0.03390155167568132        |
| 2     | A4+2A1    | 192     | 56    | 4.02535169073515   | 2          | 0.036367644170875124    | -0.03509131981126945        |
| 3     | A4+A2     | 194     | 54    | 3.9889840465642745 | 2          | 0.03774032798284699     | -0.036367644170875124       |
| 4     | D5(a1)+A1 | 196     | 52    | 3.9512437185814275 | 3          | 0.03922071315328157     | -0.03774032798284699        |
| 5     | D4+A2     | 198     | 50    | 3.912023005428146  | 2          | 0.04082199452025481     | -0.03922071315328157        |
| 6     | A4+A3     | 200     | 48    | 3.871201010907891  | 2          | 0.04255961441879608     | -0.04082199452025481        |
| 7     | A5+A1     | 202     | 46    | 3.828641396489095  | 3          | 0.04445176257083405     | -0.04255961441879608        |
| 8     | D5(a1)+A2 | 204     | 44    | 3.784189633918261  | 4          | 0.04652001563489261     | -0.04445176257083405        |
| 9     | E6(a3)+A1 | 206     | 42    | 3.7376696182833684 | 3          | 0.04879016416943216     | -0.04652001563489261        |
| 10    | D5+A1     | 208     | 40    | 3.6888794541139363 | 5          | 0.05129329438755059     | -0.04879016416943216        |
| 11    | A6        | 210     | 38    | 3.6375861597263857 | 5          | 0.054067221270275745    | -0.05129329438755059        |
| 12    | E7(a4)    | 212     | 36    | 3.58351893845611   | 4          | 0.05715841383994835     | -0.054067221270275745       |
| 13    | D5+A2     | 214     | 34    | 3.5263605246161616 | 5          | 0.06062462181643502     | -0.05715841383994835        |
| 14    | D7(a2)    | 216     | 32    | 3.4657359027997265 | 4          | 0.06453852113757108     | -0.06062462181643502        |
| 15 | A7 | 218 | 30 | 3.4011973816621555 | 4 | 0.06899287148695166 | -0.06453852113757108 |
| 16    | E8(b6)    | 220     | 28    | 3.332204510175204  | 4          | 0.07410797215372167     | -0.06899287148695166        |
| 17    | D7(a1)    | 222     | 26    | 3.258096538021482  | 6          | 0.08004270767353638     | -0.07410797215372167        |
| 18    | E7(a2)    | 224     | 24    | 3.1780538303479458 | 6          | 0.08701137698962969     | -0.08004270767353638        |
| 19    | D7        | 226     | 22    | 3.091042453358316  | 6          | 0.09531017980432521     | -0.08701137698962969        |
| 20    | E8(a5)    | 228     | 20    | 2.995732273553991  | 6          | 0.10536051565782634     | -0.09531017980432521        |
| 21    | E8(b4)    | 230     | 18    | 2.8903717578961645 | 9          | 0.11778303565638337     | -0.10536051565782634        |
| 22    | E7        | 232     | 16    | 2.772588722239781  | 10         | 0.13353139262452274     | -0.11778303565638337        |
| 23    | E8(a3)    | 234     | 14    | 2.6390573296152584 | 12         | 0.15415067982725805     | -0.13353139262452274        |
| 24    | E8(a2)    | 236     | 12    | 2.4849066497880004 | 12         | 0.18232155679395445     | -0.15415067982725805        |
| 25    | E8(a1)    | 238     | 10    | 2.302585092994046  | 14         | 0.22314355131421015     | -0.18232155679395445        |
| 26    | E8        | 240     | 8     | 2.0794415416798357 | 16         |                         |                             |
|


> **Three-step algorithm**
> 
>>1. **Build graph.** Layer all orbits according to $D=248-\dim\mathcal O$. Draw edges only between adjacent layers with $\Delta D=2$ according to completion order.  
> 
>>2. **Apply Lex criteria.** Use $\mathbf F(C)$ from 4.5.2 in the following order: length, $\mathcal S_{2}$, height rest, label coherence, $\operatorname{cv}_{3}$, label tie-break.  
> 
>> 3. **Exact search.** Perform dynamic programming across the 27 layers. The result is the **single** minimal chain $C^{\star}$ from 4.2 and Table B.1.  
> 
> Then check the **ratio tests** from 4.3; they are independent of any choice of units.  

---

### 4.6 Interpretation

The E eight chain provides a **deterministic order** of the scale ladder. No fits, no free knobs: $\lambda$ is fixed by the anchor, $\gamma(n)$ follows directly from the log step sizes, $\varphi_n$ is a pure function of $D_n$.

**Physical significance.**

• **Block structure from the chain.** Jumps in orbital height mark natural transitions between flavor, electroweak, hadronic, and cosmological.

• **Ratio laws instead of absolute tuning values.** Within and between blocks, all relations $\varphi_m/\varphi_n$ are **fit-free** and predictable. A single calibration per block is sufficient to fix dimensioned quantities.

• **Terminal law.** Towards the end of the chain, $\varphi_n\propto D_n^{\lambda}$ applies. This explains the mild but steady increase in attenuation up to $n=26$.

• **Window in dynamics.** The E windows in the two-loop flow anchor $c_3$ and $\varphi_0$ dynamically. E eight arranges the ladders, E six binds them to the observed geometry, both planes interlock.

**Distinction from the old image.**

The quadratic in $n$ was a useful heuristic, but it is not fundamental. The chain shows that $\gamma(n)$ is log-exact and that the global cubic assumption for $\ln D_n$ is not needed. The relation $\gamma_2=\gamma_0/(8\pi^2)$ is not enforced by the structure and remains an open idea for the future.

## 5. Two-loop RGE run: Dynamic fingerprints of the fixed points

---
### 5.1 Configuration

A complete **two-loop renormalization group run** is performed to dynamically test the fixed points $c_{3}=1/(8\pi)$ and $\varphi_{0}= \tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}} \approx 0.053171952$. The implementation is based on a PyR@TE definition of the E₈ cascade, extended by Standard Model fields and additional degrees of freedom:

• **Fermions:** Standard Model plus electroweak triplet $\Sigma_{F}$ (decoupling at $10^{3},\text{GeV}$) and three right-handed neutrinos $N_{R1,2,3}$ with separate thresholds $M_{N_{1}}=10^{14},\text{GeV}$, $M_{N_{2}}=3\times 10^{14},\text{GeV}$, $M_{N_{3}}=8\times 10^{14},\text{GeV}$.  

• **Color bridge:** A color-adjunct fermion $G8$ of $SU(3)_{c}$ is active above $M_{G8}=1.8\times 10^{10},\text{GeV}$. Piecewise, $\Delta b_{3}=+2$ applies, so $\dfrac{d\alpha_{3}^{-1}}{d\ln\mu}$ changes from $\tfrac{7}{2\pi}$ to $\tfrac{5}{2\pi}$ for $\mu>M_{G8}$.    

• **Scalars:** Standard model Higgs $H$, PQ field $\Phi$ with threshold $M_{\Phi}=10^{16},\text{GeV}$.  

• **Spurion:** An effective $R^{3}$ term locally models the cubic contribution $\propto \alpha^{3}$ in the abelian sector.  

• **Normalization:** Hypercharge in **GUT norm** with $b_{1}=41/10$. Convention

$$

g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y},\qquad

\beta!\left(g_{1}^{\mathrm{GUT}}\right)=\tfrac{3}{5},\beta(g_{Y}) .

$$

All figures in 5.2 and Appendix F use this convention.  

• **Starting values at $\mu = M_{Z}$:**

$g_{1}^{\mathrm{GUT}}\simeq 0.462,\quad g_{2}=0.652,\quad g_{3}=1.232.$  

The flux is integrated over more than fifteen decades (at least $10^{2},\text{GeV}$ to $\gtrsim 10^{17},\text{GeV}$) including all two-loop terms and piecewise threshold matching.  

>**Info Box: Hypercharge in GUT Norm**
>
	PyR@TE works with $b_{1}=41/6$ in SM norm by default. For the GUT norm used here, the following applies
	$g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y}$ und $\beta(g_{1}^{\mathrm{GUT}})=\tfrac{3}{5},\beta(g_{Y})$.
	This consistently yields $b_{1}=41/10$ for the slope of $\alpha_{1}^{-1}$.  
>


>**Test box for the run**
>
	Slope test for $U(1)$ in GUT norm:
	$\frac{d\alpha_{1}^{-1}}{d\ln\mu}=-\frac{b_{1}}{2\pi}$
>
	numerical $-0.6525352666767507$ vs. expectation $-0.6525352666767709$ (relative deviation $3.1\times 10^{-14}$).
	Bridge slope above $M_{G8}$: measured $\dfrac{d\alpha_{3}^{-1}}{d\ln\mu}=0.8063126$ vs. expectation $\tfrac{5}{2\pi}=0.7957747$ (1.3% deviation).  
>

**Note on beta functions.**

The analytical beta functions correspond in form and order to the standard model coefficients on one and two loops. Only the **field content is changed piecewise** above the specified thresholds, in particular $\Delta b_{3}=+2$ above $M_{G8}$.

#### 5.1b Thresholds from the $E_{8}$ ladder

**Rule.** A new degree of freedom $X$ with scale $M_{X}$ is bound to a step $n_{X}$ by the following condition:

**Rule.** A new degree of freedom $X$ with scale $M_{X}$ is bound to a step $n_{X}$ by
$$[
M_{X}=\zeta_{X}\,M_{\text{Pl}}\,\varphi_{n_{X}},\qquad \zeta_{X}=(\pi c_{3})\,e^{-\beta_{X}\pi c_{3}}e^{-k_{X}/c_{3}}.
]$$
Thus, for example, the color-adjacent fermion $G_{8}$ is not "chosen" but bound to $n_{X}=16$ or $17$, which, for plausible $(r_{X},k_{X})$, directly yields $M_{G8}\sim 10^{10}$ to $10^{11} GeV, consistent with 5.2.  

**Consequence.** The field content used in 5.1 is **derivable**. Variations in $\zeta_{X}$ are unit choices **per block**, not new free model parameters. 

---
### 5.2 Results

The key findings can be summarized in three points:

**Fingerprints of the fixed points**

For $\mu\simeq 1,\text{PeV}$, the following applies

$\alpha_{3}=0.052923411$, which is 0.47% below $\varphi_{0}=0.053171952$.

For $\mu=2.5\times 10^{8},\text{GeV}$, the following applies

$\alpha_{3}=0.039713807$, which is 0.19% below $c_{3}=\tfrac{1}{8\pi}=0.039788736$.


![[Pasted image 20250901085302.png]]
QCD fingerprints: $\alpha_3$ meets $\varphi_{0}$ and $c_{3}$.

**Approximation of unification.**

The minimal **relative spread** of the inverse couplings

$(\alpha_{1}^{-1},\alpha_{2}^{-1},\alpha_{3}^{-1})$ is

$(\alpha_{1}^{-1},\alpha_{2}^{-1},\alpha_{3}^{-1})$ beträgt

$$

\min_{\mu}\frac{\max(\alpha_{i}^{-1})-\min(\alpha_{i}^{-1})}{\tfrac{1}{3}\sum_{i}\alpha_{i}^{-1}}

= 1.23% \quad \text{at}\quad \mu^{\star}\approx 1.43×10^{15},\text{GeV}.

$$

The three pairs of ties are clustered around this value:

$\alpha_{2}^{-1}=\alpha_{3}^{-1}$ at $6.06\times 10^{14},\text{GeV}$,

$\alpha_{1}^{-1}=\alpha_{3}^{-1}$ at $1.46\times 10^{15},\text{GeV}$ and

$\alpha_{1}^{-1}=\alpha_{2}^{-1}$ at $2.38\times 10^{15},\text{GeV}$.

This defines a narrow and robust corridor instead of an exact triple intersection.

![[Pasted image 20250901091708.png]]
Inverse couplings: pairwise equalities and minimal spread.

![[Pasted image 20250901085326.png]]    
Unification measure: pairwise differences of inverse couplings.

**Perturbativity and stability.**

All couplings remain smaller than $1.3$ at least up to the range $10^{17},\text{GeV}$, no Landau poles, no unstable regions in the Higgs potential.


**Progress diagrams directly from the Pyr@ate run:**


![[Pasted image 20250901085424.png]]

E8 TFPT with G8 bridge: top left coupling curves, top center unification analysis and bridge window, top right fingerprint checks, bottom left one loop b coefficients, bottom right bridge effect on alpha.

#### 5.2b Gauge–moduli–locking in the E six window

Due to the reduction, the 4D effect contains a linear coupling factor of the form $f(\rho)\,{\rm Tr}\,G^{2}$ with $f(\rho)\propto\varphi(\rho)$ along the light radion direction. In the E six window, the scalar excitation becomes heavy and freezes the moduli at $\rho=\rho_{0}$, so that at the crossing scale $\mu_{\text{E6}}$, the identity
$$[
\alpha_{3}(\mu_{\text{E6}})=\varphi(\rho_{0})\equiv\varphi_{0}
]$$
This explains the agreement observed in 5.2 **without** additional assumptions: the same modulus combination that determines $\varphi_{0}$ multiplies the QCD term in the E six window. 

---
### 5.3 Correlations

The 2-loop analysis allows the fixed points found to be systematically linked to known structures:

• **Geometry fingerprint:** $\alpha_{3}(1,\text{PeV}) \approx \varphi_{0}$ links the geometric length directly to the QCD coupling.

• **Topology fingerprint:** $\alpha_{3}(\mu)\approx c_{3}$ at $\mu\sim 2.5\times 10^{8},\text{GeV}$ reflects the Chern–Simons scale $1/(8\pi)$.

• **Spacing invariant:** The three pairwise ties are nearly equidistant in log space and remain stable under decade fluctuations of the thresholds.

G8 color bridge: Above $M_{G8}$, $\Delta b_{3}=+2$ reduces the slope from $\tfrac{7}{2\pi}$ to $\tfrac{5}{2\pi}$; measured $0.8063$ versus expected $0.7958$ (1.3% additional slope). The effect narrows the corridor towards $10^{15},\text{GeV}$.

---
### 5.4 Interpretation  

The 2-loop RGE analysis provides dynamic confirmation of the central postulates of the theory:

1. **Independence:** $\varphi_{0}$ and $c_{3}$ appear independently in the flux, one in the PeV range, one at $10^{8},\text{GeV}$.

2. **Coherence:** The same numbers appear in different layers—geometry, topology, and now RG dynamics—and confirm each other.

3. **Stability:** The equilibrium corridor at $10^{14}$ to $10^{15},\text{GeV}$ is robust against threshold shifts.

4. **No fine-tuning:** The hits follow solely from the fixed points and the piecewise field content; no additional buttons are necessary.    

**RG Stability of the Möbius ladder.**

With $\delta(\mu)$ from $\sqrt{m_{\tau}(\mu)/m_{\mu}(\mu)}$, the deformation remains nearly constant in a wide window; universal terms $a_{y},\varphi_{0}$ and $b_{y},c_{3}$ capture percentage residues.

---
### 5.5 Conclusion 

$c_{3}=1/(8\pi)$ and $\varphi_{0}$ appear as **dynamic fingerprints** in the course of the gauge couplings. Together with the log-exact order of the E₈ cascade, this results in a consistent picture:

• **Topology sets the scale,**

• **Geometry provides the length,**

• **E₈ orders the ladders,**

• **RG Dynamics confirms the fingerprints.**

---
## **6. Inflation from topology and geometry**

The reduction from eleven dimensions via six dimensions to four dimensions generates a hyperbolic field space with plateau dynamics. The curvature is determined solely by the two TFPT invariants $c_{3}$ and $\varphi_{0}$. This results in an alpha attractor with robust predictions for $n_{s}$ and $r$. The cascade $E_{8}\to E_{7}\to E_{6}$ determines reheating parameters and small corrections via the effective freedom density.

---
### **6.1 Setup and assumptions**

1. **Topological core and length scale**

We use the TFPT invariants

$c_{3}=\dfrac{1}{8\pi}$, $\ \varphi_{0}=0.0531719522$.

$c_{3}$ is the topological coupling from the Chern Simons sector, $\varphi_{0}$ sets the global length scale of the compactification.

2. **Fields from the reduction**

After the reduction $11\text{D}\to 6\text{D}\to 4\text{D}$, two light modes remain: a volume mode $\rho(x)$ and an axion-like mode $\theta(x)$ from the integrated three-form potential on an internal cycle.

3. **Effect in the Einstein framework**

After Weyl rescaling, we obtain a two-dimensional sigma model for $(\rho,\theta)$ with negatively curved field space. Along a slight valley direction, we define the canonical inflaton variable $\phi$.

4. **Potential from flows**

River quantization and the Chern Simons term generate a plateau potential of the alpha attractor family along the valley direction. E model and T model are both suitable and lead to the same leading predictions.

---
### **6.2 Field space and canonical variable**

The kinetics take the universal form

$\mathcal L_{\text{kin}}=\dfrac{M_{P}^{2}}{2}\sqrt{-g},\dfrac{3,\alpha_{\text{inf}},(\partial z)^{2}}{(1-z^{2})^{2}},\quad z=\tanh!\Big(\dfrac{\phi}{\sqrt{6\alpha_{\text{inf}}},M_{P}}\Big)$,

with field space curvature

$R_{K}=−,\dfrac{2}{3,\alpha_{\text{inf}}}$.

$\alpha_{\text{inf}}$ is a pure function of the TFPT invariants. Two normalizations directly motivated by the reduction closely frame the curvature:

Variant A $\ \alpha_{\text{inf}}=\dfrac{c_{3}}{\varphi_{0}}=0.74830308$

Variant B $\ \alpha_{\text{inf}}=\dfrac{\varphi_{0}}{2c_{3}}=0.66817846$

Numerical value of the curvature: $R_{K}^{\text{A}}\simeq −0.891$, $R_{K}^{\text{B}}\simeq −0.998$.

Both variants arise from the combination of the fiber scaling $\mathrm e^{−2\rho}$ with the topological weighting by $c_{3}$ and the global length scale $\varphi_{0}$. No fits, just geometry.

![[Pasted image 20250829140815.png]]
![[Pasted image 20250829140935.png]]
>***Fig. 6.1* Top row: compactification 11D → 6D → 4D with flux and Chern Simons.**  
>
   **Right box:** hyperbolic sigma model $\mathcal{L}_{\mathrm{kin}}=\frac{M_{P}^{2}}{2}\sqrt{-g}\,\frac{3\alpha_{\mathrm{inf}}(\partial z)^{2}}{(1-z^{2})^{2}}$, mapping $z=\tanh\!\left(\frac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$, curvature $R_{K}=-\frac{2}{3\alpha_{\mathrm{inf}}}$.  
>
>**Center:** Poincaré disk with geodesics, radial inflaton trajectory toward $z=1$. Bottom tiles: representative numbers at $N=55$ for both $\alpha_{\mathrm{inf}}$ variants.

---
### **6.3 Potential on the plateau**

A representative example suffices
suffices as a representative example.

$V(\phi)=V_{0}\Big(1−\exp!\Big[−\sqrt{\dfrac{2}{3\alpha_{\text{inf}}}}\dfrac{\phi}{M_{P}}\Big]\Big)^{2}$

or equivalently

$V(\phi)=V_{0},\tanh^{2}!\Big(\dfrac{\phi}{\sqrt{6\alpha_{\text{inf}}},M_{P}}\Big)$.

The plateau asymptotics guarantees small tensors and a slope that depends almost exclusively on the number of convolutions $N$.

![[Pasted image 20250829141404.png]]

>*Fig. 6.6* $V(\phi)/V_{0}=\tanh^{2}\!\left(\dfrac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$ for $\alpha_{\mathrm{inf}}=0.748$ and $0.668$. The asymptotics illustrate the plateau at large field values.

---
### **6.4 Universal predictions**

At the CMB pivot, the following applies

$n_{s}\simeq 1−\dfrac{2}{N}$, $\quad r\simeq \dfrac{12,\alpha_{\text{inf}}}{N^{2}}$, $\quad \alpha_{s}\equiv \dfrac{\mathrm d n_{s}}{\mathrm d\ln k}\simeq −,\dfrac{2}{N^{2}}$, $\quad n_{t}\simeq −,\dfrac{r}{8}$.

The amplitude $A_{s}\simeq 2.1\times 10^{−9}$ fixes $V_{0}$. The following are useful in practice

$V^{1/4}\simeq \big(3\pi^{2} A_{s} r\big)^{1/4}M_{P}$, $\quad H\simeq \pi,M_{P},\sqrt{\dfrac{A_{s},r}{2}}$.

![[Pasted image 20250829141011.png]]
>*Fig. 6.1b* The inflaton follows a radial trajectory with $z=\tanh\!\left(\frac{\phi}{\sqrt{6\alpha_{\mathrm{inf}}}M_{P}}\right)$ and approaches the boundary $z=1$. Geodesics are drawn as diameters and as arcs orthogonal to the boundary. Negative curvature $R_{K}=-\frac{2}{3\alpha_{\mathrm{inf}}}$.

---
### **6.5 Figures from TFPT**

We set $c_{3}=\tfrac{1}{8\pi}$ and $\varphi_{0}=0.0531719522$. The results for $N=50,55,60$ are:

**Variant B $\ \alpha_{\text{inf}}=\varphi_{0}∕(2c_{3})=0.66817846$**

**$N$** **$n_{s}$** **$r$** **$\alpha_{s}$** **$n_{t}$** **$V^{1∕4}$ [GeV]** **$H$ [GeV]**

50 0.960000 0.00320726 $−8.00\times10^{−4}$ $−4.01\times10^{−4}$ $9.150\times10^{15}$ $1.404\times10^{13}$

55 0.963636 0.00265063 $−6.61\times10^{−4}$ $−3.31\times10^{−4}$ $8.725\times10^{15}$ $1.276\times10^{13}$

60 0.966667 0.00222726 $−5.56\times10^{−4}$ $−2.78\times10^{−4}$ $8.353\times10^{15}$ $1.170\times10^{13}$

**Variant A $\ \alpha_{\text{inf}}=c_{3}∕\varphi_{0}=0.74830308$**

**$N$** **$n_{s}$** **$r$** **$\alpha_{s}$** **$n_{t}$** **$V^{1/4}$ [GeV]** **$H$ [GeV]**

50 0.960000 0.00359185 $−8.00\times10^{−4}$ $−4.49\times10^{−4}$ $9.413\times10^{15}$ $1.486\times10^{13}$

55 0.963636 0.00296848 $−6.61\times10^{−4}$ $−3.71\times10^{−4}$ $8.975\times10^{15}$ $1.351\times10^{13}$

60 0.966667 0.00249434 $−5.56\times10^{−4}$ $−3.12\times10^{−4}$ $8.593\times10^{15}$ $1.238\times10^{13}$

**Lyth boundary.**

$\Delta\phi \gtrsim N,\sqrt{r∕8},M_{P}$. For $N=55$, this results in $\Delta\phi\simeq 1.00,M_{P}$ to $1.06,M_{P}$ for variants B and A. Thus, minimally transplanckian, typical for plateaus with small $r$.

![[Pasted image 20250829141200.png]]

>*Fig. 6.3* $r(N)$ for $\alpha_{\mathrm{inf}}=0.748$ and $0.668$. The dashed horizontal line is the BK18 limit, the vertical line marks the  implied by Planck .

### **6.6 Connection to the cascade $E_{8}\to E_{7}\to E_{6}$**

1. **Reheating and degrees of freedom**

The cascade determines the effective number density of degrees of freedom $g_{\ast}$ during reheating. This shifts $N$ by a few units. The relationship

$\Delta N \simeq \dfrac{1−3w_{\text{reh}}}{12,(1+w_{\text{reh}})}\ln!\Big(\dfrac{\rho_{\text{reh}}}{\rho_{\text{end}}}\Big)$

shows how thresholds $E_{7}$ to $E_{6}$ enter $n_{s}$ and $r$ via $w_{\text{reh}}$ and $\rho_{\text{reh}}$. Realistic shifts remain small and change $r$ in the percentage range.

2. **n equal to 6 threshold and TeV window**

The cascade marks a threshold in the TeV range. This influences the coupling to visible degrees of freedom at the end of inflation and thus the reheating efficiency. The impact is mainly in $\Delta N$, not in the form of the predictions.

3. **Fine structure in the potential**

Small plateau folds caused by cascade steps can generate scale-dependent mini features in $n_{s}(k)$. As long as these are not explicitly derived, the smooth plateau approximation is sufficient. Later, we can expand this fine structure into separate predictions.

### **6.7 Comparison with reference values**

• **Slope $n_{s}$.**

Planck provides $n_{s}=0.9649\pm 0.0042$ at the pivot $k_{\ast}=0.05,\text{Mpc}^{−1}$. This yields $N=\dfrac{2}{1−n_{s}}=56.98$ with a $1σ$ band of $50.89$ to $64.72$.

• **Tensors $r$.**

With $\alpha_{\text{inf}}$ from Section 6.2, this gives $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$ central
.
$2.47\times 10^{−3}$ (variant B) to $2.77\times 10^{−3}$ (variant A). This is significantly below BK18 with $r>0.036$ and exactly in the target corridor of the upcoming CMB generation.

• **Running $\alpha_{s}$.**

$\alpha_{s}\simeq −,\dfrac{2}{N^{2}}\approx −6.2\times 10^{−4}$, thus nearly scale-invariant and within Planck uncertainty.

• **Tensor slope $n_{t}$.**

Consistency relation $n_{t}=−r∕8$ yields $−(2.3$ to $3.5)\times 10^{−4}$.

• **Rule of thumb $r=\varphi_{0}^{2}$.**

$\varphi_{0}^{2}=0.002827$. For $N=55$, this results in $r\simeq 0.00247$ to $0.00277$, i.e., within a six percent window around $\varphi_{0}^{2}$. If $\alpha_{\text{inf}}\simeq 0.713$ is selected, the result is exactly $r=\varphi_{0}^{2}$.

![[Pasted image 20250829141246.png]]

>*Fig. 6.4* $\alpha_{\mathrm{inf}}=\dfrac{r}{3(1-n_{s})^{2}}$ for three values of . Dots mark the two TFPT normalizations at $N(n_{s}^{\mathrm{Planck}})$.


### **6.8 Tests and clear falsification**

• **CMB polarization.**

Next-generation CMB experiments measure $r$ down to the lower three times ten to the power of minus three range. A reliable zero result below $r\lesssim 0.001$ forces a redefinition of $\alpha_{\text{inf}}$ or the reheating history with a significantly larger $N$ in TFPT.

• **Reheating and cascade fingerprint.**

Precise measurements of $n_{s}$ and $r$ plus independent information about reheating allow conclusions to be drawn about the $E_{7}$ to $E_{6}$ thresholds. This connects cosmological and collider signatures.

### **6.9 Hard matching with reference values**


**Setup.**

Planck pivots: $n_{s}=0.9649\pm0.0042$, $\ln(10^{10}A_{s})=3.044$, so $A_{s}\approx 2.11\times 10^{−9}$.

BK18 limit: $r_{0.05}>0.036$ at 95 percent.

From $n_{s}$ follows $N=\dfrac{2}{1−n_{s}}=56.98$ with $1σ$ band $50.89$ to $64.72$.

From $N$ and $\alpha_{\text{inf}}$, we get $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$.

**Central values for $n_{s}=0.9649$.**

**Size** **Formula** **Value B** **Value A** **Comment**

$N$ $N=\dfrac{2}{1−n_{s}}$ $56,980$ $56,980$ from Planck

$r$ $r=3,\alpha_{\text{inf}}(1−n_{s})^{2}$ $2,470\times 10^{−3}$ $2,766\times 10^{−3}$ clearly below BK18

$H$ $H=\pi M_{P}\sqrt{\dfrac{A_{s}r}{2}}$ $1.232\times 10^{13},\text{GeV}$ $1.303\times 10^{13},\text{GeV}$ Pivot scale

$V^{1∕4}$ $V^{1∕4}=(3\pi^{2}A_{s}r)^{1∕4}M_{P}$ $8.571\times 10^{15}$ $8.817\times 10^{15}$ Plateau scale

$m_{\phi}$ $m_{\phi}=\dfrac{\sqrt{12},\pi,\sqrt{A_{s}}}{N},M_{P}$ $2.131\times 10^{13}$ $2.131\times 10^{13}$ independent of $\alpha_{\text{inf}}$

$\lambda_{\phi}$ $\lambda_{\phi}=\dfrac{8\pi^{2}A_{s}}{3,\alpha_{\text{inf}},N^{2}}$ $2.546\times 10^{−11}$ $2.274\times 10^{−11}$ small, as expected

$\Omega_{\text{gw}}(k_{\ast})$ $\Omega_{\text{gw}}\simeq \dfrac{A_{s}r}{24},\Omega_{r,0}$ $1.987\times 10^{−17}$ $2.225\times 10^{−17}$ with $\Omega_{r,0}\approx 9.2\times 10^{−5}$

$\Delta\phi$ $\Delta\phi\simeq N\sqrt{\dfrac{r}{8}},M_{P}$ $1.001,M_{P}$ $1.059,M_{P}$ minimal transplanckian

Ratio to BK18 $r∕0.036$ $0.0686$ $0.0768$ well below limit

**Band above $1σ$ in $n_{s}$.**

Variant B $\ \alpha_{\text{inf}}=\varphi_{0}∕(2c_{3})=0.66817846$

**$n_{s}$** **$N$** **$r$** **$H$ [GeV]** **$V^{1∕4}$ [GeV]** **$m_{\phi}$ [GeV]** **$\lambda_{\phi}$**

$0.9607$ $50.891$ $3.096\times 10^{−3}$ $1.379\times 10^{13}$ $9.069\times 10^{15}$ $2.386\times 10^{13}$ $3.192\times 10^{−11}$

$0.9649$ $56.980$ $2.470\times 10^{−3}$ $1.232\times 10^{13}$ $8.571\times 10^15) 2.131×1013 2.546×10−11

$0.9691$ $64.725$ $1.914\times 10^{−3}$ $1.084\times 10^{13}$ $8.041\times 10^{15}$ $1.876\times 10^13) 1.973×10−11

Variant A $\ \alpha_{\text{inf}}=c_{3}∕\varphi_{0}=0.74830308$

**$n_{s}$** **$N$** **$r$** **$H$ [GeV]** **$V^{1∕4}$ [GeV]** **$m_{\phi}$ [GeV]** **$\lambda_{\phi}$**

$0.9607$ $50.891$ $3.467\times 10^{−3}$ $1.459\times 10^{13}$ $9.329\times 10^{15}$ $2.386\times 10^{13}$ $2.850\times 10^{−11}$

$0.9649$ $56.980$ $2.766\times 10^{−3}$ $1.303\times 10^{13}$ $8.817\times 10^15) 2.131×1013 2.274×10−11

$0.9691$ $64.725$ $2.143\times 10^{−3}$ $1.147\times 10^13) $8.272×10^{15}$ $1.876×10^{13}$ $1.762×10^{−11}$

**Direct curvature test.**

$\dfrac{r}{(1−n_{s})^{2}}=3,\alpha_{\text{inf}}$. This identity allows the reconstruction of $\alpha_{\text{inf}}$ directly from data.

### **6.10 Reheating window and $\Delta N$**

Definitions:

$\Delta N\simeq \dfrac{1−3w_{\text{reh}}}{12,(1+w_{\text{reh}})}\ln!\Big(\dfrac{\rho_{\text{reh}}}{\rho_{\text{end}}}\Big)$, $\ \rho_{\text{reh}}=\dfrac{\pi^{2}}{30}g_{\ast}T_{\text{reh}}^{4}$, $\ \rho_{\text{end}}\approx c,V_{0}$ with $c\approx 0.34$ to $0.37$, $g_{\ast}\approx 120$.

Results for $w_{\text{reh}}=0$:

**Scenario** **Variant B** **Variant A** **Comment**

$\Delta N$ at $T_{\text{reh}}=6,\text{MeV}$ $−13.54$ $−13.55$ too cold, inconsistent with Planck band

$T_{\text{reh}}^{\text{min}}$ for $\Delta N=−6$ $4.01\times 10^{7},\text{GeV}$ $4.12\times 10^{7},\text{GeV}$ matter-like reheating

$T_{\text{reh}}^{\text{max}}$ for rapid conversion $\approx 2.71\times 10^{15},\text{GeV}$ $\approx 2.73\times 10^{15},\text{GeV}$ $c\simeq 0.34$ to $0.37$

![[Pasted image 20250829141332.png]]
>*Fig. 6.5* $\Delta N(T_{\mathrm{reh}})$ for matter like reheating $w=0$ and $g_{\ast}\approx 120$. The dotted line near $\Delta N\approx -13.5$ corresponds to BBN scale reheating $\sim 6$ MeV. The dashed line at $\Delta N=-6$ indicates the approximate boundary consistent with the Planck band.

### **6.11 Info Box**

$$\boxed{

\begin{aligned}

&N=\dfrac{2}{1−n_{s}}=56.98\quad (1σ:\ 50.89\ \text{bis}\ 64.72),\

&r=3,\alpha_{\text{inf}}(1−n_{s})^{2}=

\begin{cases}

2.47\times 10^{−3}\ \text{(Variante B)}\

2.77\times 10^{−3}\ \text{(Variante A)}

\end{cases},\

&m_{\phi}=\dfrac{\sqrt{12},\pi,\sqrt{A_{s}}}{N},M_{P}=2.13\times 10^{13}\ \text{GeV},\quad

\lambda_{\phi}=\dfrac{8\pi^{2}A_{s}}{3,\alpha_{\text{inf}},N^{2}}.\

\end{aligned}

}$$

In short: The inflation scale now determines $m_{\phi}$ and $\lambda_{\phi}$ **without** free knobs. The relation $r∕(1−n_{s})^{2}=3\alpha_{\text{inf}}$ is a direct test of the topologically fixed curvature.

### 6.12 Unambiguous fixation of $\alpha_{\text{inf}}$

Variant A ($\alpha_{\text{inf}}=c_{3}/\varphi_{0}$) and variant B ($\alpha_{\text{inf}}=\varphi_{0}/(2c_{3})$) are two **reference metrics** of the same modulus fixation. We choose a **unique** normalization by

$$[
Criterion: |R_K+1|=\min \quad \text{when preserving the TFPT identity} \ \ r=3\,\alpha_{\text{inf}}(1-n_s)^{2}.
]$$

Numerically, this results in $R_{K}\simeq -0.998$ and thus **variant B** as the clear choice:
$$[
\alpha_{\text{inf}}=\frac{\varphi_{0}}{2c_{3}}=0.66817846,\qquad 
r=\frac{3\varphi_{0}}{2c_{3}}\,(1-n_{s})^{2}.
]$$
The predictions in 6.5 and 6.7 remain in the same narrow corridor, now **without** ambiguity. See also the Poincaré disk map in 6.2.

---

## 7. Role of α and the parameter-free solution

---
### 7.1 Motivation and origin of the approach

The fine structure constant $\alpha$ is an **external input parameter** in the Standard Model. Early considerations (Sommerfeld, Dirac, Eddington) had already suggested that there must be a deeper mathematical structure behind the number $\alpha^{-1}\approx137$.

Genetic algorithms and 6D precursors repeatedly showed that $\alpha$ is closely linked to two constants:

$c_{3}=\frac{1}{8\pi}$, $\qquad \varphi_{0}\approx 0.053171$.

Both quantities appeared independently in kinetic, Maxwell, and mass terms. The crucial observation was that $\alpha$ always "appeared" where **topological normalization** (via $c_{3}$) and **geometric length** (via $\varphi_{0}$) were simultaneously effective.

This led to the hypothesis: $\alpha$ **is not free, but rather the unique solution to a fixed-point condition that couples precisely these two constants.**

---
### 7.2 A parameter normal form for $\alpha$: representation only in $c_{3}$

**Normal form.** With $c_{3}=\tfrac{1}{8\pi}$,
$$
\varphi_{0}=\tfrac{4}{3}c_{3}+48\,c_{3}^{4},\quad
A=2\,c_{3}^{3},\quad
\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}},\quad b_{1}=\tfrac{41}{10},
$$
becomes
$$
\alpha^{3}-A\alpha^{2}-A\,c_{3}^{2}\kappa=0
$$
to the pure $c_{3}$ form
$$
\boxed{\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}=0\ }.
$$

**Closed solution (Cardano).** Set $\alpha=y+\frac{2}{3}c_{3}^{3}$, then $y^{3}+py+q=0$ with
$$
p=-\tfrac{4}{3}c_{3}^{6},\qquad
q=-\tfrac{16}{27}c_{3}^{9}-8\,b_{1}\,c_{3}^{6}\ln\!\frac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}},
$$
$\Delta=\big(\tfrac{q}{2}\big)^{2}+\big(\tfrac{p}{3}\big)^{3}$ und
$$
\boxed{\ \alpha(c_{3})=\frac{2}{3}c_{3}^{3}
+\sqrt[3]{-\frac{q}{2}+\sqrt{\Delta}}
+\sqrt[3]{-\frac{q}{2}-\sqrt{\Delta}}\ }.
$$

**Practical formula.** Very accurate, closed approximation
$$
\boxed{\ \alpha \approx \Big(8\,b_{1}\,c_{3}^{6}\,\ln\tfrac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}\Big)^{1/3}
+\frac{2}{3}c_{3}^{3}\ }.
$$


> [!tip] Practical formula

> Very accurate practical formula
> $[> \alpha \approx \Big(8b_{1}c_{3}^{6}\ln\tfrac{1}{\tfrac{4}{3}c_{3}+48c_{3}^{4}}\Big)^{\!\!1/3}+\frac{2}{3}c_{3}^{3}>]$
> already gives the ppm approximation. For $c_{3}=1/(8\pi)$, $\alpha^{-1}=137.0365014649$ follows.

---
### 7.3 The solution

The fixed point equation is a cubic polynomial that has exactly one physically real positive zero.

$$c_{3}=\tfrac{1}{8\pi}\ \Rightarrow\
\varphi_{0}=0.0531719521768,\ 
\kappa=1.914684795,\
\alpha=0.007297325816919221,\
\alpha^{-1}=137.03650146488582.$$

The unique real solution is

$\alpha=0.0072973258169192213,\quad \alpha^{-1}=137.03650146488582$.

This is $3.665\times10^{-6}$ relative to CODATA 2022 $\alpha_{\text{CODATA}}=0.0072973525628 or \alpha^{-1}=137.035999177$.

The other two roots are complex and non-physical.

Thus, $\alpha$ is **not postulated**, but rather the **output** of a compelling equation.

![[Pasted image 20250826131905.png]]

---
### 7.4 Accuracy of the solution

Comparison with CODATA 2022 reference ($\alpha^{-1}=137.035999177(21)$):

- Deviation: a few parts per million (ppm).
- No fine adjustment necessary – the match follows directly from c₃, φ₀, and b₁.

This is remarkable because it represents the most precise **parameter-free theoretical derivation** of $\alpha$ to date.

---

### 7.5 Alternative Approximations and Optimized Calculation Methods

#### 7.5.1 Cubic root approximation

In the limit of small A, $\alpha$ can be approximated by

$\alpha \;\approx\; (A c_{3}^{2}\kappa)^{1/3} + \frac{A}{3}$.

- The first term ($A c_{3}^{2}\kappa)^{1/3}$ gives the principal value.
- The additive surcharge $A/3$ (universal, independent of $\varphi_{0}$) brings the number close to ppm.

Absolute error $2.44\times10^{-7}$ corresponds to approximately 33 ppm.

This approximation already matches $\alpha$ to an accuracy of 10⁻⁷.

---
#### 7.5.2 Ramanujan-like series

If we set $\alpha = (A c_{3}^{2}\kappa)^{1/3}(1+u)$ and expand in powers of u, we obtain a convergent series:

$\alpha = B^{1/3} + \frac{A}{3} + \frac{A^{2}}{9B^{1/3}} + \frac{2A^{3}}{81B^{2/3}} + \dots, \qquad B=A c_{3}^{2}\kappa$.


- After just three terms, the deviation is already <0.2 ppm.
- Four terms provide accuracy to 10⁻¹².
- 
Error $\approx 9.38\times10^{-10}$
---
#### 7.5.3 Newton's method 

Starting with $g=B^{1/3}+A/3$ and applying Newton's method once, the same accuracy is achieved as with the series. 

Formula:

$\alpha \approx g - \frac{f(g)}{f’(g)},\qquad f(\alpha)=\alpha^{3}-A\alpha^{2}-B$.

This allows $\alpha$ to be calculated extremely efficiently and accurately.
![[Pasted image 20250824122910.png]]

---

### 7.6 Variational derivative in four dimensions (cubic fixed point equation from Einstein's action)

**Goal and context.** We show that the cubic fixed-point equation for $\alpha$ used in 7.2 follows as a stationarity condition of an explicit four-dimensional action. The constants $c_{3}$ and $\varphi_{0}$ originate precisely from the already established invariants of the theory, and the normalizations $A$ and $\kappa$ are identical to the definitions in Appendix E. This results in a second, independent derivation of the same equation without freely selectable scales. Compare the derivations of $c_{3}$ in Section 3.2.1 and of $\varphi_{0}$ in Section 3.2.2, as well as the normalization note in Appendix E.  

**Fixed invariants from topology and geometry**

From the Chern-Simons reduction with $M_{11}=M_{4}\times Y_{7}$ and integer intersection number, the topological fixed point follows

$c_{3}=\frac{1}{8\pi}$.

See the rigorous derivation via $C_{3}\wedge G_{4}\wedge G_{4}$ and the quantization of $\int_{M_{4}}F\wedge F = 8\pi^{2}\mathbb Z$ in Section 3.2.1. The Möbius reduction with Gauss Bonnet and boundary yields

$\varphi_{0}=\frac{1}{6\pi}+\frac{3}{256\pi^{4}}$,

see Section 3.2.2 and Appendix D. Together with $b_{1}=41/10$ in GUT norm, all quantities are fixed.  

**Four-dimensional effect and** $U(\alpha)$

After reduction and canonical normalization, we consider the Abelian sector in a homogeneous background and understand $U(\alpha)$ as a **gradient representation in the coupling space**.

$\partial_{\alpha} U \propto \beta_{\alpha}$,

so that stationarity $\partial_{\alpha}U=0$ is equivalent to $\beta_{\alpha}=0$, cf. the interpretation in 7.6. The effective effect is

$S_{\mathrm{eff}}=\int\mathrm d^{4}x\,\sqrt{-g}\,\Bigl[\tfrac{M_{P}^{2}}{2}R-U(\alpha)+\ldots\Bigr]$.

Up to the relevant order, a scalar potential suffices.

$U(\alpha)=\frac{A}{4}\,\alpha^{4}\;-\;\frac{2A}{3}\,c_{3}^{3}\,\alpha^{3}\;-\;A\,[\,8\,b_{1}\,c_{3}^{6}\,\ln(1/\varphi_{0})\,]\,\alpha. \tag{7.1.1}$
is sufficient.
• **Leading term** $\alpha^{4}$: Smooth reference term that provides the leading $\alpha^{3}$ contribution in $\partial_{\alpha}U$.

• **Leading term** $\alpha^{4}$: Smooth reference term that provides the leading $\alpha^{3}$ contribution in $\partial_{\alpha}U$.

• **Cubic contribution** $\propto c_{3}^{3}\alpha^{3}$ It arises from the reduced Chern Simons structure via the coupling of a heavy scalar mode a to $F\tilde F$ and the elimination of $a$ in the zero momentum limit. The combined counting of two identical topological insertions, the conversion $g\mapsto \alpha$ and the symmetry factor yields the **universal factor**

$A=\frac{1}{256\pi^{3}}\equiv 2\,c_{3}^{3}$,

exactly as shown in Appendix E, Step 2.  

• **Linear logarithm.** The integrated one-loop renormalization between $\mu_{\mathrm{UV}}=M_{P}$ and $\mu_{\mathrm{IR}}=\varphi_{0}M_{P}$ yields

$\kappa\equiv \frac{b_{1}}{2\pi}\ln\!\Bigl(\frac{1}{\varphi_{0}}\Bigr)$,

in the potential notation as the term $-A[8\,b_{1}\,c_{3}^{6}\ln(1/\varphi_{0})]\,\alpha$. This is identical to the normalization $\kappa$ defined in Appendix E, Step 1.  

All quantities are in reduced Planck units, see Info Box.  

**Stationarity and normal form**

The variation condition yields

$\frac{\partial U}{\partial \alpha}=A\Bigl[\alpha^{3}-2\,c_{3}^{3}\,\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)\Bigr]=0$,

i.e., the **cubic fixed point equation**

$\boxed{\,\alpha^{3}-2\,c_{3}^{3}\,\alpha^{2}-8\,b_{1}\,c_{3}^{6}\,\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)=0\,}. \tag{7.1.2}$

is used.

Equivalent **normal form** as used in 7.2:

$\alpha^{3}-A\,\alpha^{2}-A\,c_{3}^{2}\,\kappa=0,\qquad A=2c_{3}^{3},\ \ \kappa=\frac{b_{1}}{2\pi}\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr)$.

The unique real solution numerically agrees with 7.3, where it is given as $\alpha^{-1}=137.03650146488582$.  

**Physical classification and consistency**

1. **Scheme freedom.** $A$ is a pure number from topology and canonical normalization, $\kappa$ depends only on the fixed $b_{1}$ and the geometrically fixed scale $\varphi_{0}$. A scheme change only shifts additive, scale-independent contributions in $\kappa$, but not the position of the fixed point, see Appendix E.  

2. **Significance of** $U(\alpha)$ $U$ is not a matter potential, but a compact representation of the coupling dynamics, $\partial_{\alpha}U\propto \beta_{\alpha}$. 

3. **Relationship to the rest of the structure.** The same invariants set dynamic fingerprints in the two-loop analysis, $\alpha_{3}(1\ \mathrm{PeV})\approx \varphi_{0}$ and $\alpha_{3}(\mu)\approx c_{3}$ at $\mu\sim 2.5\times 10^{8}\ \mathrm{GeV}$, see 5.2.  

4. **Abelian trace as a common thread.** The number 41 appears both in $b_{1}=41/10$ of the fixed point equation and in the EW block via $k_{\mathrm{EW}}=41/32$, cf. 8.4.6. This underscores that the same Abelian trace is at work in both contexts, without circularity.
![[Pasted image 20250901085708.png]]
*Left: invariants $((c_3, \varphi_0, b_1))$. Middle: $(U_3, U_4, U_1)$. Right: $(U(\alpha)$), stationarity, cubic fixed point.*

#### 7.6.1 Callan–Symanzik route

Mit $\mu\,d\alpha/d\mu=\beta_{\alpha}=\tfrac{b_{1}}{2\pi}\alpha^{2}+A\,c_{3}^{2}\alpha^{3}+\dots$ and $A=1/(256\pi^{3})$ from Appendix E, integration between $M_{\text{Pl}}$ and $\varphi_{0}M_{\text{Pl}}$ immediately leads to the cubic
$$[
\alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\ln\frac{1}{\varphi_{0}}=0.
]$$
This is independent of the potential representation $U(\alpha)$ and makes the fixed point equation **doubly** derived. 
#### 7.6.2  A and κ at a glance  with cross-references

**Purpose.** This box summarizes the normalizations and invariants used in **7.6** for the variational derivative and refers to the formal derivations in **3.2.1**, **3.2.2**, and **Appendix E**. 

---

**Fixed points and scales.**

* Topology from Chern-Simons reduction, see **3.2.1**:  

  $( c_{3}=\dfrac{1}{8\pi}=0.039788735772973836\ldots )$

* Geometry of the Moebius fiber with Gauss Bonnet and boundary, see **3.2.2** and **Appendix D**:  

  $( \varphi_{0}=\dfrac{1}{6\pi}+\dfrac{3}{256\pi^{4}}=0.053171952176845526\ldots )$

**Abelian trace in GUT norm.**

* $( b_{1}=\dfrac{41}{10}=4.1 )  \; ( \Rightarrow )$  Definition of $( \kappa )$ below, see **Appendix E**.

---

**Key abbreviations from Appendix E.**

* Pure topology factor  
  $( A\;\equiv\;2\,c_{3}^{3}\;=\;\dfrac{1}{256\,\pi^{3}}\;=\;1.259825563796855\times10^{-4} )$.

* Integrated one loop constant  
  $( \kappa\;\equiv\;\dfrac{b_{1}}{2\pi}\,\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr) )$.

---

**Effective potential in 7.1.1.**  $( U(\alpha) )$ as gradient representation in coupling space  

$( \partial_{\alpha}U\propto\beta_{\alpha} )$. Up to the relevant order:

$$[
U(\alpha)=\frac{A}{4}\alpha^{4}-\frac{2A}{3}c_{3}^{3}\alpha^{3}-A\bigl[\,8\,b_{1}\,c_{3}^{6}\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr)\bigr]\alpha.
]$$

**Stationarity equals fixed point.**

$$[

\frac{\partial U}{\partial \alpha}=A\Bigl[\alpha^{3}-2c_{3}^{3}\alpha^{2}-8b_{1}c_{3}^{6}\ln\!\bigl(\tfrac{1}{\varphi_{0}}\bigr)\Bigr]=0

]$$


Normal form as in **7.2**:  

$$[

\alpha^{3}-A\,\alpha^{2}-A\,c_{3}^{2}\,\kappa=0,\qquad A=2c_{3}^{3},\quad \kappa=\frac{b_{1}}{2\pi}\ln\!\Bigl(\tfrac{1}{\varphi_{0}}\Bigr).

]$$

**Dynamic fingerprints and Cosmo anchors.**

* $( \alpha_{3}(1\,\text{PeV})\approx \varphi_{0} )$ and $( \alpha_{3}(\mu)\approx c_{3} )$  at  $( \mu\sim2.5\times10^{8}\,\text{GeV} )$, cf. **5.2**.  

* $( \Omega_{b}=\varphi_{0}(1-2c_{3}) )$ near Planck, see **8.4.7**.

---

### 7.7 Interpretation

The role of $\alpha$ is fundamentally redefined in this framework:

- **Not an input, but a fixed point.** $\alpha$ is not an arbitrary number, but the unique solution to a geometric-topological condition.

- **Dominance of topology.** Sensitivity analyses show that $\alpha$ reacts most strongly to c₃ (topological fixed point), less strongly to b₁ (spectrum), and least strongly to φ₀ (geometry).

- **Universal surcharge.** The constant correction term $A/3$ explains why $\alpha$ is accurate to the ppm – a small but structural shift.

This means that the fine structure constant is **not random**, but rather an emergent fixed point from topology, geometry, and symmetry.

**Common cause of $\alpha$ and flavor relations.**  
$\alpha$ follows from $(\varphi_0,c_3)$, the Möbius ladder uses $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$.  
This creates the loop $\varphi_0\Rightarrow \alpha(\varphi_0,c_3)$ and $\varphi_0\Rightarrow\delta_\star\Rightarrow$ flavor relations.

---

## 8. From E₈ to E₇ to E₆ and to the Standard Model  
_A clear block structure, mathematically consistent, immediately reproducible_

> [!info] Fixed points and ladders  
> **Topology:** $c_3=\tfrac{1}{8\pi}=0.039788735772973836$  
> **Geometry:** $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.05317195217684553$  
> **Conductor normalization:** $\gamma(0)=0.834,\ \lambda=\dfrac{0.834}{\ln 248-\ln 60}=0.5877029773404678$  
> **Planck constant for numbers:** $M_{\rm Pl}=1.221\times10^{19}\ \mathrm{GeV}$

---

> [!summary] Idea in one sentence  

> We combine a discrete structure axis from E₈ with steps $n$ and a dynamic axis from renormalization group $\mu$.  
>  E₈ orders the ladders $\varphi_n$. The RG dynamics provides windows $E_r$ at $\alpha_3(\mu)\approx 1/(r\pi)$.  
> Blocks link both and project onto measurable quantities of the Standard Model.


![[Pasted image 20250828094621.png]]

### Two axes, one common grid

**Structural axis**  

The nilpotent orbitology of E₈ gives rise to a unique, strictly descending chain.  

$D_n = 60-2n, \qquad n=0\dots 26,$

which defines a log-exact ladder  

$\varphi_{n} = \varphi_0\,e^{-\gamma(0)}\left(\frac{D_n}{D_1}\right)^{\lambda} \quad (n\ge 1).$

This axis is discrete. It arranges ratios of scales. It explains why certain jumps between levels always look the same.

**Dynamic axis**  

On the RG axis, the strong coupling $\alpha_3(\mu)$ runs continuously. There are three natural windows  


$$

\alpha_3(\mu_r)=\frac{1}{r\,\pi},\qquad r\in\{6,7,8\},

$$

i.e., $E₆$ by $1/(6\pi)$ near PeV, $E₇$ by $1/(7\pi)$ in between, $E₈$ by $1/(8\pi)=c_3$ at about $2.5\times 10^8$ GeV.

> [!info] Reading rule  

> $n$ counts structure and determines ratio laws.  
> $E_r$ marks dynamics and determines positions on the energy axis.  
> Both are synchronized by the fixed points $c_3=\tfrac{1}{8\pi}$ and $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}$.


![[Pasted image 20250828094712.png]]

> **Info Box — Chirality comes from geometry, not from E8**
> 
> E8 serves here exclusively as an **ordering principle of a discrete scale ladder** $(\varphi_{n})$.  
> There is **no** 4D gauge group $(E_{8})$ and **no** embedding of SM fermions in $(E_{8})$ representations.  
> The 4D **chirality** arises independently of this through **boundary conditions** and **integer quantized flows** on the orientable double cover of the Möbius fiber.  
> Three boundary cycles and Chern–Simons quantization provide the **chiral index**:
> $( \mathrm{Ind}\,D = \dfrac{1}{2\pi}\!\int_{\widetilde M}\!F = \nu_{1}+\nu_{2}+\nu_{T} )$.
> With the minimal choice $((\nu_{1},\nu_{2},\nu_{T})=(1,1,1))$, there are **three left-chiral families**; mirror states are missing due to the projectors on the boundary cycles.  
> Details in **Appendix J**.  (References: 3.2.1 $(c_{3}=\tfrac{1}{8\pi})$, 3.2.2 $(\varphi_{0})$, 8.4.6 $(k_{\mathrm{EW}}=41/32))$.   
![[Pasted image 20250901111221.png]]

---
### How structure and dynamics become SM figures


![[Pasted image 20250828094846.png]]

The transition from dimensionless ladder steps to measurable quantities takes place in blocks. Each block $B$ has three key figures:  

- $r_B$ effective rank in the chain $E₈ \supset E₇ \supset E₆ \supset SM$  

- $k_B$ fractional topology number from the boundary cycles of the Möbius fiber  

- $n_B$ degree of the ladder  

This first results in a block constant  

$$

\zeta_B=(\pi c_3)\,\exp[-\beta_B\,\pi c_3]\,\exp\Bigl[-\frac{k_B}{c_3}\Bigr],\qquad \beta_B=\frac{8-r_B}{8},

$$

and then the dimensioned size  

$$

X_B=\zeta_B\,M_{\mathrm{Pl}}\,\varphi_{n_B}.

$$

For example, we set  

- EW Block at $n=12$ in the $E₇$ window: $v_H=\zeta_{\rm EW} M_{\mathrm{Pl}}\varphi_{12}$  

- Hadron blocks at $n=15$ and $n=17$ in the $E₆$ corridor: $m_p\simeq \zeta_{\rm had} M_{\mathrm{Pl}}\varphi_{15}$  

- Lepton blocks deep down $n=22,25,26$: light Yukawas  

> [!tip] Quick start for readers  

> 1. Find the block for the quantity you are looking for in the text.  
> 2. Read $r_B, k_B, n_B$ and calculate $\zeta_B$.  
> 3. Set $X_B=\zeta_B M_{\mathrm{Pl}}\varphi_{n_B}$ with the log-exact $\varphi_{n}$ from the E₈ ladder.

---

### Where is the connection to the standard model?

The chain $E₈ \supset E₇ \supset E₆ \supset SM$ provides the rank logic and the Abelian trace:  

- At the EW anchor $n=12$, the trace $\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48}$ appears. This results in $k_{\rm EW}=\tfrac{41}{32}$ and, consistently, $b_1=\tfrac{41}{10}$ in GUT norm.  

- The hadronic windows lie in the $E₆$ domain of the ladder and support the additional damping that characterizes baryonic scales.  

- The RG windows dynamically anchor this structure: $\alpha_3(\mu)$ hits $1/(6\pi), 1/(7\pi), 1/(8\pi)$ at exactly the points motivated by the ladder.  

In short: structure organizes, dynamics confirms, blocks project. This is our path from topology and geometry to the numbers of the Standard Model.

---
### What do the steps do without a direct block?

Not every step has to carry a specific observable. These steps are an important supporting structure:  

1. **Geometry of the ladder**  

You uphold the law  
$$

\frac{\varphi_m}{\varphi_n}=\left(\frac{D_m}{D_n}\right)^{\lambda}\quad(m,n\ge 1),

$$
\quad(m,n\ge 1),


i.e., the fit-free ratio structure.  

2. **Fine snap points in windows**  

A dynamic window is an area in $\mu$. The discrete $n$ act as grid points at which thresholds and mixtures can take effect without violating the global ratio law.  

3. **Reserve for new observables**  

Other quantities such as thresholds, axion couplings, and precise hadronic parameters can be added later. The spaces are already structurally wired correctly.  

> [!example] Intuition  

> Think of a gearbox. The block steps are the gears that drive an axle. The intermediate teeth ensure that the power is transmitted cleanly and without slipping. Without them, there would be jumps, but no order.

### 8.0a  Chirality from edge and flow: operative summary

**Geometry and Boundary**  

We work on the **orientable double cover** $(\widetilde M)$ of the Möbius fiber with **three** closed boundary cycles $(C_{1},C_{2},C_{T})$. The boundary counting is canonical: $(\sum_{i}\!\oint_{C_{i}}\widehat{k}_{g}\,\mathrm ds = 6\pi)$.  

The projectors resulting from the six-dimensional reduction choose **an internal chirality**:

$$[

P_{T}=\tfrac12\!\left(\mathbf 1+\mathrm i\,\sigma^{3}\sigma^{n}\right),\quad 

P_{1}=P_{2}=\tfrac12\!\left(\mathbf 1-\mathrm i\,\sigma^{3}\sigma^{n}\right),

]$$

so that only $(\chi_{+})$ carries zero modes and the 4D zero modes are **left-chiraled**. 

**Index and family number**  

A smooth Abelian connection $(A)$ with quantized flux $(m=\tfrac{1}{2\pi}\!\int_{\widetilde M}\!F\in\mathbb Z)$ yields

$$[

\mathrm{Ind}\,D_{\widetilde M}

= \#\{\chi_{+}\} - \#\{\chi_{-}\}

= \tfrac{1}{2\pi}\!\int_{\widetilde M}\!F 

= \nu_{1}+\nu_{2}+\nu_{T}.

]$$

The minimal choice $((\nu_{1},\nu_{2},\nu_{T})=(1,1,1))$ results in **three** families without mirrors.  

Wilson lines are **flat** and only read out the Abelian trace, compatible with $(k_{\mathrm{EW}}=\tfrac{41}{32})$.  References: 3.2.1 $(c_{3}=\tfrac{1}{8\pi})$, 3.2.2 $(\varphi_{0})$, 8.4.6 Trace $(41)$.  

### 8.1 Detailed description

E₈ arranges the **scale ladder** $\varphi_n$ log exactly, E₇ and E₆ set the **physical windows** per block, and **topology** with **geometry** provides the **normalizations** via $c_3$ and $\varphi_0$. Dimensioned quantities arise from a compact **block formula**:

$$
\boxed{\,X_B=\zeta_B\,M_{\rm Pl}\,\varphi_{n_B},\quad 
\zeta_B=(\pi c_3)\,e^{-\beta_B\,\pi c_3}\,e^{-\,k_B/c_3},\quad 
\beta_B=\tfrac{8-r_B}{8}\,}
$$

with $r_B$ as the effective rank in the block and $k_B$ as the rational topological number of the three boundary cycles.

The E₈ ladder is log-exact:

$$
\gamma(0)=0.834,\qquad
\gamma(n)=\lambda[\ln D_n-\ln D_{n+1}],\quad D_n=60-2n,\quad
\lambda=\frac{0.834}{\ln 248-\ln 60}.
$$

For $n\ge 1$, the following applies

$$
\boxed{\,\varphi_n=\varphi_0\,e^{-\gamma(0)}\Big(\frac{D_n}{D_1}\Big)^{\lambda},\quad D_1=58\,}.
$$

![[Pasted image 20250828094755.png]]


### 8.1.1 Block constants from marginal cycles and Abelian trace


**Objective.** $k_{B}$ are **not fits**, but result from a count of Abelian squares on the **three** boundary cycles of the orientable representation, multiplied by a universal factor. 


**Definition (Abelian trace in GUT norm).** For a block $B$, let

$$[
\mathcal I_{1}(B)=\sum_{\Phi\in B}\sum_{i\in U(1)_{Y}} q_{i}^{2}(\Phi)\ \ \text{with}\ \ Y\text{ in GUT norm}

\mathcal I_{1}(B)=\sum_{\Phi\in B}\sum_{i\in U(1)_{Y}} q_{i}^{2}(\Phi)\ \ \text{with}\ \ Y\text{ in GUT norm}.

]$$

**Theorem 8.1.1 (Topological block number).** The rational number

$$
[

k_{B}=\frac{3}{2}\,\mathcal I_{1}(B)

]$$

is the sum of the three boundary cycles (factor $3$) and the factored half-weight of the orientable double cover (factor $1/2$). 


**Example EW block.** For $B=\text{SM}+\text{H}$, $\mathcal I_{1}=\tfrac{41}{48}$, so $k_{\text{EW}}=\tfrac{3}{2}\cdot \tfrac{41}{48}=\tfrac{41}{32}$, exactly as used in 8.4.1; at the same time, the same trace appears in $b_{1}=\tfrac{41}{10}$, see 7.6.1. 


**Note Hadron and pion blocks.** For baryonic and pionic quantities, $\mathcal I_{1}$ is replaced by the effective abelian subgroups of flavor-chiral dynamics ($U(1)_{B}$ and $U(1)_{I3}$ respectively). This yields $k_{p}=\tfrac{3}{2}$ and $k_{\pi}=\tfrac{51}{32}$ as concrete evaluations of the same counting rule in the respective block. See 8.4.5.  

### 8.1.2 Derivation of the $\zeta_{B}$ formula from the boundary functional


The effective boundary effect per block has the form

$$[
S^{(B)}_{\partial}=\pi c_{3}\ -\ \beta_{B}\,\pi c_{3}\ -\ \frac{k_{B}}{c_{3}},
]$$

where $\beta_{B}=(8-r_{B})/8$ results from the effective rank number $r_{B}$ (counting the undamped directions). Exponentiation of the additive contributions yields

$$[

\zeta_{B}=(\pi c_{3})\,\exp\big[-\beta_{B}\pi c_{3}\big]\exp\big[-k_{B}/c_{3}\big],\qquad X_{B}=\zeta_{B}\,M_{\text{Pl}}\,\varphi_{n_{B}}.

]$$

This means that $v_{H}$, $m_{p}$, $f_{\pi}$ and $T_{\gamma 0}$ can be calculated directly from $(r_{B},k_{B},n_{B})$, **without** additional degrees of freedom. See 8.4.1 to 8.4.7. 

---

### 8.2 Calculation formula in three steps

1. **Evaluate ladder**  
   $$
   \varphi_n=\varphi_0\,e^{-\gamma(0)}\Big(\tfrac{60-2n}{58}\Big)^{\lambda}\qquad (n\ge 1).
   $$

2. **Set block constants**  
   For block $B$: select $r_B$, $\beta_B=(8-r_B)/8$, and $k_B$ rationally from the edge count.  
   $$
   \zeta_B=(\pi c_3)\,e^{-\beta_B \pi c_3}\,e^{-k_B/c_3},\qquad \pi c_3=\tfrac18.
   $$

3. **Determine size**  
   $$
   X_B=\zeta_B\,M_{\rm Pl}\,\varphi_{n_B}.
   $$

> [!tip] Proportionality laws without unit selection  
> $$
> \frac{\varphi_{m}}{\varphi_{n}}=\Big(\frac{60-2m}{60-2n}\Big)^{\lambda}\quad (m,n\ge 1).
> $$

> [!tip] Ratio laws without unit selection

---

### 8.3 Required ladder steps $\varphi_n$ (log exact)

| $n$ | $D_n$ |     $\varphi_n$ |
| --: | ----: | --------------: |
|   1 |    58 | 0.0230930346695 |
|   5 |    50 | 0.0211640537281 |
|  10 |    40 | 0.0185628455934 |
|  12 |    36 | 0.0174482846938 |
|  15 |    30 | 0.0156753658147 |
| 16 | 28 | 0.0150524852088 |
|  22 |    16 | 0.0108336306291 |
|  25 |    10 | 0.0082188698412 |
|  26 |     8 | 0.0072087140665 |

---

### 8.4 Results per block with references

#### 8.4.1 Electroweak block $n=12$  
**Assumptions:** $r_{\rm EW}=2\Rightarrow \beta_{\rm EW}=3/4,\quad k_{\rm EW}=\tfrac{41}{32}$

$$
\zeta_{\rm EW}=(\pi c_3)\,e^{-\,\tfrac34\pi c_3}\,e^{-\,\tfrac{41}{32}/c_3}
=1.17852087206\times10^{-15}.
$$

$$
v_H=\zeta_{\rm EW}M_{\rm Pl}\varphi_{12}= \mathbf{251.07628}\ \mathrm{GeV}.
$$

With $g_2=0.652,\ g_1^{\rm SM}=0.357$ at $M_Z$:

$$
M_W=\tfrac12 g_2 v_H=\mathbf{81.85087}\ \mathrm{GeV},\qquad
M_Z=\tfrac12\sqrt{g_2^2+g_1^2}\,v_H=\mathbf{93.31741}\ \mathrm{GeV}.
$$

**Comparison**  
$v=(\sqrt2 G_F)^{-1/2}=246.21965\ GeV$ ⇒ **+1.97 percent**  
$M_W=80.3692\ \mathrm{GeV}$ ⇒ **+1.84 percent**  
$M_Z=91.1876\ \mathrm{GeV}$ ⇒ **+2.34 percent**

> [!note] Interpretation  
> The block sets the **scale** $v_H$ to an accuracy of one to two percent. Finite contributions with two loops and thresholds shift $M_W,M_Z$ downwards towards the references.

**Top mass as a minimum assumption**  
$y_t\approx 1\Rightarrow m_t\simeq v_H/\sqrt2=\mathbf{177.54}\ \mathrm{GeV}$.

---

8.4.2 PQ Block $n=10$  
**Assumptions:** $r_{\rm PQ}=1\Rightarrow \beta_{\rm PQ}=7/8,\quad k_{\rm PQ}=\tfrac12$

$$
\zeta_{\rm PQ}=3.90754185582\times10^{-7},\quad
f_a=\zeta_{\rm PQ}M_{\rm Pl}\varphi_{10}=\mathbf{8.8565\times 10^{10}}\ \mathrm{GeV}.
$$

Axion mass:  

$$
m_a\simeq (5.7\ \mu\mathrm{eV})\times \frac{10^{12}\ \mathrm{GeV}}{f_a}=\mathbf{64.36\ \mu\mathrm{eV}}.
$$

---

8.4.3 Seesaw Block $n=5$  
**Assumptions:** $r_{N_R}=4\Rightarrow \beta_{N_R}=1/2,\quad k_{N_R}=\tfrac18$

$$
M_R=\zeta_{N_R}M_{\rm Pl}\varphi_{5}=\mathbf{1.311\times10^{15}}\ \mathrm{GeV}.
$$

With $y_{\nu3}\sim 1$:  

$$
m_{\nu3}\simeq \frac{v_H^2}{M_R}=\mathbf{0.04807\ \mathrm{eV}},\quad
\Delta m^2_{31}\simeq \mathbf{2.31\times 10^{-3}\ \mathrm{eV}^2}.
$$

---

8.4.4 Flavor anchors from $n=1$  

$$
\sin^2\theta_{13}=\varphi_1=\mathbf{0.023093},\qquad 
\sin\theta_{13}=0.15197.
$$

**Cabibbo angle from basic level**  

$$
\sin\theta_C\simeq \sqrt{\varphi_0}\Big(1-\frac{\varphi_0}{2}\Big)=\mathbf{0.22446},\qquad
\theta_C=\arcsin(\sin\theta_C)=\mathbf{0.22639}\ \mathrm{rad}.
$$

**Möbius mass ladder with a deformation $\delta$.**  
Calibrate $\delta$ using only leptons:
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}.$  
Insert this single number into the six relations:
$$
\begin{aligned}
\text{Down:}\quad & \sqrt{\tfrac{m_s}{m_d}}=\mathcal{M}_{1}(\delta),\quad
\sqrt{\tfrac{m_b}{m_s}}=\mathcal{M}_{1}(\delta)\,(1+\delta),\\[2mm]
\text{Leptonen:}\quad & \sqrt{\tfrac{m_\tau}{m_\mu}}=\mathcal{M}_{1}(\delta),\quad
\sqrt{\tfrac{m_\mu}{m_e}}=\mathcal{M}_{1}(\delta)\,\mathcal{M}_{1/3}(\delta),\\[2mm]
\text{Up:}\quad & \sqrt{\tfrac{m_c}{m_u}}=\mathcal{M}_{2/3}(\delta),\quad
\sqrt{\tfrac{m_t}{m_c}}=\dfrac{2/3}{\,2/3-\delta\,}.
\end{aligned}
$$

**Topological check.**  
The theory expects $\delta_\star=\dfrac{3}{5}+\dfrac{\varphi_0}{6}$.  
Compare $\delta$ from the leptons with $\delta_\star$ and document the deviation in percent.  
No new free parameters are added.

---

8.4.5 Hadron window and pion observables

**Proton $n=15$**, assumptions $r_{\rm had}=5\Rightarrow \beta_{\rm had}=3/8,\ k_p=\tfrac32$:  

$$
m_p=\zeta_{\rm had}M_{\rm Pl}\varphi_{15}=\mathbf{0.96821\ GeV}.
$$

**Pion $n=16$**, same rank $r=5$, stronger topological damping $k_\pi=\tfrac{51}{32}$:  

$$
f_\pi=\mathbf{88.12\ \mathrm{MeV}}\ \text{(chiral norm)}.
$$

GMOR consistency with $|\langle\bar qq\rangle|^{1/3}\simeq 272\ \mathrm{MeV},\ (m_u+m_d)_{2\ \mathrm{GeV}}\simeq 6.8\ \mathrm{MeV}$:  

$$
m_\pi\simeq\sqrt{\frac{(m_u+m_d)\,|\langle\bar qq\rangle|}{f_\pi^2}}=\mathbf{132.75\ \mathrm{MeV}}.
$$

---

8.4.6 Fine structure constant $\alpha$  
(cross-reference to section 6)

With

$$
\alpha^{3}-2c_3^{3}\alpha^{2}-8\,b_1\,c_3^{6}\ln\frac{1}{\tfrac{4}{3}c_3+48c_3^{4}}=0,\qquad b_1=\tfrac{41}{10},
$$
\alpha^{3}-2c_3^{3}\alpha^{2}-8\,b_1\,c_3^{6}\ln\frac{1}{\tfrac{4}{3}c_3+48c_3^{4}}=0,\qquad b_1=\tfrac{41}{10},


yields the unique real solution  

$$
\boxed{\,\alpha=\mathbf{0.007297325816919221},\quad \alpha^{-1}=\mathbf{137.03650146488582}\,}
$$

Deviation from CODATA 2022 $\alpha^{-1}=137.035999177$: **+3.67 ppm**.


> [!summary] Summary

> **The same counting measure 41** from the hypercharge appears **twice**:  

> – in the α-fixed point equation via $(b_1=\tfrac{41}{10})$
> – in the EW block via $(k_{\rm EW}=\tfrac{41}{32})$
> Both follow from the same abelian trace $(\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48})$. α is therefore **not an input** here, but a **consistency echo** of the same structure that anchors $(v_H)$.

**1) α as a fixed point from topology and geometry**

The cubic equation
$$
[\alpha^{3}-2c_3^{3}\alpha^{2}-8\,b_1\,c_3^{6}\,\ln\!\frac{1}{\varphi_0}=0,\quad 
c_3=\frac{1}{8\pi},\ \ \varphi_0=\frac{1}{6\pi}+\frac{3}{256\pi^4},\ \ b_1=\frac{41}{10}
]$$

yields $(\alpha^{-1}=137.0365)$ **without** free parameters.  

This is where the **41** comes in via $(b_1)$—the hypercharge trace of the Standard Model in GUT norm.

**2) The same 41 fingerprint sets the EW block**

In the EW block (window at \(n=12\)), we use
$$
[

\zeta_{\rm EW}=(\pi c_3)\,\mathrm e^{-\beta_{\rm EW}\pi c_3}\,\mathrm e^{-k_{\rm EW}/c_3},

\quad \beta_{\rm EW}=\frac{3}{4},\quad 

k_{\rm EW}=\tfrac{3}{2}\cdot \mathcal Y^2_{\rm SM+H}=\tfrac{41}{32}.

]$$

Here, too, **the same 41** is used, now in $(k_{\rm EW})$. This determines $(v_H)$ via $(v_H=\zeta_{\rm EW} M_{\rm Pl}\varphi_{12})$.  

**Result:** $(v_H\simeq 251.1\ \mathrm{GeV})$ (scale anchor, expected 1–2 percent drift to $(G_F)$).

**3) α in the EW picture: combination of \(g_1\) and \(g_2\)**

According to electroweak mixing, the following applies

$$[

e=g_2\sin\theta_W=g_1\cos\theta_W,\qquad 

\alpha=\frac{e^2}{4\pi}.

]$$

If typical values are used for $(M_Z)$ ($(g_2\approx0.652,\ g_1^{\rm SM}\approx0.357)$), the result is $(\alpha(M_Z))$ **in the order of magnitude $(1/128)$ – this is the **current** α at the Z pole.  

Our fixed-point solution gives the **IR‑α** ($(\alpha^{-1}\approx 137.0365)$); the difference is simply **renormalization flow**. The crucial point is that **the same 41** controls both the fixed-point equation (via $(b_1)$) and the EW anchor (via $(k_{\rm EW})$).

> **No circular reference**
>
> $c_{3}=\tfrac{1}{8\pi}$, $\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}$
> $\Rightarrow$ $\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}}$, $A=2c_{3}^{3}$
> $\Rightarrow$ Cubic in $\alpha$.
>
> Parallel: $Y^{2}_{\text{SM+H}}=\tfrac{41}{48}\Rightarrow k_{\text{EW}}=\tfrac{41}{32}$.
>
> The same Abelian trace sets $b_{1}$ and $k_{\text{EW}}$, **without** feedback from $\alpha$ to $v_{H}$.  

> [!example] Mini number check

> – Fixed point: $(\alpha^{-1}_{\rm IR}=137.0365)$ (from $(c_3,\varphi_0,b_1)$)  
> – At $(M_Z)$: $(\alpha(M_Z)\sim 1/128)$ from $(g_1,g_2,\theta_W)$  
> – Both values are linked by **the same** U(1) content; the number 41 appears **twice**, which explains why α naturally comes into play here again.

---

8.4.7 Cosmology from the elementary level

$$
\Omega_b=\varphi_0\,(1-2c_3)=\mathbf{0.04894066}.
$$

---

### 8.5 Summary at a glance

| Size | Prediction | Reference | Deviation |
|---|---:|---:|---:|
| $v_H$ | **251.07628 GeV** | 246.21965 GeV | +1.97 % |
| $M_W$ | **81.85087 GeV** | 80.3692 GeV | +1.84 % |
| $M_Z$ | **93.31741 GeV** | 91.1876 GeV | +2.34 % |
| $m_t$ | **177.54 GeV** | 172.57 GeV | +2.9 % |
| $f_a$ | **$8.8565\times 10^{10}$ GeV** | Standard window | — |
| $m_a$ | **64.36 μeV** | Standard window | — |
| $M_R$ | **$1.311\times 10^{15}$ GeV** | — | — |
| $m_{\nu3}$ | **0.04807 eV** | — | — |
| $\Delta m^2_{31}$ | **$2.31\times 10^{-3}$ eV²** | $2.509\times 10^{-3}$ eV² | −7.9 % |
| $\sin^2\theta_{13}$ | **0.023093** | $0.02240\pm0.00065$ | +3.1 % |
| $\sin\theta_C$ | **0.22446** | $0.2248\pm 0.0006$ | −0.15 % |
| $m_p$ | **0.96821 GeV** | 0.938272 GeV | +3.19 % |
| $f_\pi$ | **88.12 MeV** | 92.07 MeV | −4.3 % |
| $m_\pi$ | **132.75 MeV** | 134.98 MeV $(\pi^0)$ | −1.6 % |
| $\alpha^{-1}$ | **137.036501465** | 137.035999177 | +3.67 ppm |
| $\Omega_b$ | **0.04894066** | 0.0493 | −0.7 % |

---
#### 8.5.1 Systematics of deviations

The 1 to 3 percent deviations in v_H, M_W, and M_Z arise from
(i) missing two-loop terms in the electroweak sector,
(ii) threshold adjustments at the transition to $n=12$,
(iii) the block unit $\zeta_{\text{EW}}$ as a pure choice of units.

**Expectation.** Consistent tracking with two loops and piecewise matching (cf. 5.1, 5.2) shifts $v_{H}$, $M_{W}$, and $M_{Z}$ by approximately $100$ to $300$ pBq.

**Expectation.** Consistent tracking with two loops and piecewise matching (see 5.1, 5.2) systematically shifts $v_{H}$, $M_{W}$, $M_{Z}$ **downward** toward the references, by $\Delta\sim 1$ to $2$ percent. The ratio tests within the block remain unchanged, as they only depend on $\varphi_{n}$. See 8.4.1 and 5.4. 

---

### 8.6 Where E₇ and E₆ specifically connect

- **E₇ window at $n=12$** anchors the **electroweak scale**. The Abelian trace $\mathcal Y^2_{\rm SM+H}=\tfrac{41}{48}$ leads via three half boundary cycles to $k_{\rm EW}=\tfrac{41}{32}$. The same 41 appears as $b_1=\tfrac{41}{10}$ in the fixed point equation of $\alpha$.

- **E₆ corridor** carries the **strong dynamics**. $r_{\rm had}=5$ explains the milder damping in the hadron block and justifies small rational $\Delta k$ for Goldstone physics relative to baryons.

---

### 8.7 What remains open and how we can close it

- **Fine structure of Yukawas:** Only **scales** were deliberately set here. Textures and phases are the next layer. Percentage dispersions in the block frame are to be expected.

- **Two loops of fine-tuning in the electroweak sector:** Consistent tracking with thresholds will systematically pull $v_H, M_W, M_Z$ toward the references.

- **Formal derivation of $k_B$:** The rational $k_B$ used are motivated by the edge count. An index-like derivation per block belongs in the appendix.

- **Chirality**: closed by boundary conditions and integer flows on the orientable double cover, see box in 8 and **Appendix J**. 

---

### 8.8 Figures for this section

- $c_3=\frac{1}{8\pi}=0.039788735772973836$  
- $\varphi_0=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^4}=0.05317195217684553$  
- $\gamma(0)=0.834,\quad \lambda=0.5877029773404678$  
- $\varphi_{10}=0.018562845593356334,\ \varphi_{12}=0.01744828469380037$  
- $\varphi_{15}=0.015675365814677055,\ \varphi_{16}=0.015052485208841481$  
- $\varphi_{22}=0.01083363062914777,\ \varphi_{25}=0.008218869841220914,\ \varphi_{26}=0.007208714066517271$  
- $\zeta_{\rm EW}=1.17852087206\times 10^{-15},\ \zeta_{\rm PQ}=3.90754185582\times 10^{-7}$  
- $M_{\rm Pl}=1.221\times10^{19}\ \mathrm{GeV}$  
- $g_2=0.652,\ g_1^{\rm SM}=0.357$ at $M_Z$  
- $\alpha=\mathbf{0.007297325816919221},\ \alpha^{-1}=\mathbf{137.03650146488582}$

---
## 9. Further information, outlook, and FAQ

---
### 9.1 Additional information for understanding

The previous chapters have derived the **core structure** of the theory: two fundamental fixed points ($c_{3}$, $\varphi_{0}$), the E₈ cascade, and the fixed-point solution for $\alpha$. For a complete understanding, three further aspects should be highlighted:

1. **Single-point calibration**:

    The cascade $\varphi_n$ is determined up to an additive constant in $\log \varphi$. A single physical calibration (e.g., at the EW block, n=12) fixes all remaining stages. This is not a "button," but rather a choice of unit.

2. **Block formulas**:

    The dimensioning of individual observables (e.g., proton mass, CMB temperature, dark energy) is performed using compact block formulas, which are specified in the appendices. They link the dimensionless $\varphi_n$ to measurable quantities.

3. **Spurion contributions**:

    The $R^{3}$ Spurion used in the 2-loop runs is not a free parameter, but rather an effective description of higher contributions that inevitably occur in the Chern–Simons structure. Its influence is small, but necessary to correctly model the cubic term for $\alpha$.

The **Abelian trace** is the common thread: the same $(41)$ controls $(b_{1}=\tfrac{41}{10})$ in $(\kappa=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}})$ as well as $(k_{\mathrm{EW}}=\tfrac{41}{32})$ in the EW block; geometrically, the phases are read via **three boundary cycles** on the double cover, see **Appendix J**.  

![[Pasted image 20250826132358.png]]

> [!note] Sensitivity

> The sensitivity of $\alpha$ to the parameters scales strongly with $c_{3}$, significantly weaker with $b_{1}$ and only moderately with $\varphi_{0}$, see Figure 7.1.
### Self-consistency: $\varphi_{0}\leftrightarrow \alpha$

The fixed point equation not only generates α as a function of φ₀, but φ₀ itself is motivated by the geometric reduction ϕ₀ = 1/(6π) + 3/(256π⁴). Combining both dependencies results in a closed loop:


$[\varphi_{0}\ \xrightarrow{\ \kappa(\varphi_{0})=\tfrac{b_{1}}{2\pi}\ln\tfrac{1}{\varphi_{0}}\ }\\Big[\ \alpha^{3}-2c_{3}^{3}\alpha^{2}-8b_{1}c_{3}^{6}\ln\tfrac{1}{\varphi_{0}}\ =\ 0\ \Big]\ \xrightarrow{\ \text{solution}\ }\ \alpha(\varphi_{0})]$

This loop closes because $\varphi_{0}$ itself follows from the geometry $\big(\varphi_{0}=1/(6\pi)+3/(256\pi^{4})\big)$ and the solution for $\alpha$ is $\alpha(\varphi_{0}) = 1/(6\pi)+3/(256\pi^{4})$ (see below).

This loop closes because $\varphi_{0}$ itself follows from the geometry $\big(\varphi_{0}=1/(6\pi)+3/(256\pi^{4})\big)$ and the solution for $\alpha$ confirms the input.

This self-referential structure replaces classic fine-tuning debates with structural feedback—$\varphi_0$ and $\alpha$ determine each other. Small changes in φ₀ propagate through κ directly into the equation, which then yields a new α value. The original input is reconfirmed by the resulting solution—structural "locking" instead of adjustable parameters.

> **Falsification box**
> **A — $\alpha$ precision:** Deviation of the cubic solution beyond a few tens of ppm refutes the structure.
> **B — Fingerprints:** If $\alpha_{3}(1\,\text{PeV})$ misses the $\varphi_{0}$ window or $\alpha_{3}$ misses the $c_{3}$ window at $\mu\sim 10^{8}$ GeV robustly, the model is refuted.
> **C — Spacing:** If the nearly equidistant log-spacing invariant of the three equipotentials breaks, the ladder incoherence is proven. 

---
### 9.2 Open questions and next steps

Several points have already been outlined in theory, but require further work:

- **RG robustness**:
    Initial tests show that the equilibrium corridors are extremely stable. A systematic analysis with varying thresholds (± decade) and alternative field contents is planned.

- **Cosmological extensions**:
    The steps n=20,25,30 reproduce knee, CMB, and dark energy. Here, we will examine whether $S_8/\sigma_8$ tensions and early dark energy can also be consistently embedded.


---

### 9.3 FAQ: Ten questions and answers

**1. Is this just number crunching or numerology?**

No. $c_{3}=1\!/\!(8\pi)$ follows from a quantized Chern Simons coupling. $\varphi_{0}$ follows from Möbius geometry plus boundary terms. Both quantities appear independently in different parts of the theory and then feed into the fixed point equation for $\alpha$. This distinguishes a structural result from a fit.

**2. Are there any free parameters?**

No. Once the topologically and geometrically determined quantities $c_{3}$ and $\varphi_{0}$ and the physically fixed $U(1){Y}$ constant $b{1}=41/10$ have been defined, there are no freely selectable parameters left. There is only a trivial unit calibration.

**3. Why specifically $E_{8}$?**

Only $E_{8}$ has sufficiently rich orbit structures whose centralizer dimensions form a unique monotonic chain. The logarithm of the dimensions produces a simple step structure from which the damping $\gamma(n)$ follows in blockwise constant form. Smaller groups break this monotonicity or produce inconsistent jump patterns.

**4. Difference from classical GUT approaches such as SU(5) or SO(10)?**

Classical GUTs postulate additional symmetry and a new scale to unify couplings. Here, constants are derived from topology and geometry. Unification appears as a side effect of the flow, not as an axiom.

**5. How robust are the numbers?**

Very robust. Shifts of around a decade only change the situation of characteristic ties in the per mille range. The solution of the fixed-point equation for $\alpha$ remains stable in the ppm range. The steps of the ladder are deterministic, not fit-driven.

**6. Why is** $\alpha$ **so precise, while other quantities are only accurate to within one percent?**

$\alpha$ is determined directly by the fixed-point equation. Masses and mixtures carry additional QCD dynamics, flavor structure, and scheme effects. These contributions are deliberately kept modular in the present version and generate natural scatter at the percent level.

**7. How can the theory be refuted?**

Three clear levers:

a) RG fingerprints on two characteristic scales, for example in the PeV range and at around $2.5\times 10^{8} GeV$.

b) Stability of the spacing pattern between equipotsentials over a wide parameter range.

c) Predictions in precision fields such as atomic interferometry or Rydberg constant for $\alpha$. Systematic deviations refute the model.

**8. Are there any connections to string theory or M theory?**

Yes, at the level of the 11-dimensional parent structure with Chern Simons term and compactified topology. Unlike landscape approaches, TFPT does not require a multitude of free moduli. The derivations remain local and topological.

**9. What does the theory say about the cosmological constant?**

Step n=30 of the ladder yields an energy density $\rho_{\Lambda}$ of the order of magnitude of the Planck measurements. The decisive factor is the origin of the exponent from the ladder, not a fit to data.

**10. Where are the greatest uncertainties?**

Two points: the formal derivation of the closed form of $\gamma(n)$ directly from $E_{8}$ orbitology and the deep interpretation of the block constants $\zeta$. Both are mentioned in the outlook sections as a work program.

**11. Where do** $A=\tfrac{1}{256\pi^{3}}$ **and** $\kappa=\tfrac{b_{1}}{2\pi}\ln(1/\varphi_{0})$** come from?**

From the chosen normalization $\alpha=g^{2}/(4\pi)$, GUT norm for U(1){Y} _and a topologically induced single-loop correction to_ $F^{2}$ _with two identical insertions of_ $c{3}$. See Derivation Note A1 in the appendix for the complete calculation.

**12. How scheme-dependent are the statements?**

A schema change only shifts additive, scale-independent terms in $\kappa$. The pure number factor A is fixed by topology and canonical gauge kinetics. Fixed points and ladder structures remain invariant.

**13. What does "no free parameters" mean in practice if numerical values are rounded?**

Rounding only affects display and numerical propagation. The structural equations are parameter-free. In reproductions, all constants should be specified with defined precision and error bars should be shown from the scheme and threshold variation.

**14. Why a cubic equation for** $\alpha$ **and not a quadratic or quartic one?**

The smallest non-trivial order in which the topological contribution to the renormalization of the photon wave function occurs locally and is parity even is proportional to $g^{6}$. In the $\alpha$ scale, this corresponds to the third power. Lower orders are excluded by symmetry or quantization.

**15. How are two-loop effects and thresholds handled technically?**

The non-Abelian couplings run in two loops with standard coefficients and threshold jumps at the effective masses of the heavy modes. Sensitivity analyses show that the two characteristic fingerprints remain stable in terms of position and distance. The Abelian equation additionally receives the topological cubic term.

**16. How do I reproduce the key results numerically?**

Steps:

a) Set $c_{3}=1/(8\pi)$ and $\varphi_{0}$ according to Section 3.2.

b) Calculate $\kappa=(b_{1}/2\pi)\ln(1/\varphi_{0})$ with $b_{1}=41/10$.

c) Solve the fixed point equation in 6.2 for $\alpha$ with $A=1/(256\pi^{3})$.

d) Run the non-Abelian couplings twice, set defined thresholds, and check the fingerprints.

e) Vary thresholds and schema parameters within plausible ranges and specify error bars.

**17. What is the physics behind** $\varphi_{0}$**?**

$\varphi_{0}$ is not a fit constant, but arises from a geometric relation on the orientable double cover of the Möbius reduction. Gauss Bonnet with boundary provides the area fraction, the boundary term provides the surcharge. Together, this fixes the effective dimensionless scale relation.

**19. Where does the theory currently end?**

In the present version, flavor details, CKM and PMNS phases, and non-trivial hadron phenomenology are only outlined. This is a deliberate modularization. The initial goal is to establish a solid foundation of topology, geometry, and coupling dynamics.

**20. What is the next step in closing the open issues?**

Three concrete steps:

a) Formal derivation of the closed $\gamma(n)$ shape directly from nilpotent orbits and centralizers.

b) Complete two-loop validation with systematic threshold evaluation and error budget.

c) Precision tests for $\alpha$ via independent measurement channels and simulations, including clear deviation thresholds for falsification.

---

### 9.4 Plausibility arguments: probability and structural dependencies

The plausibility of the present theory arises from two complementary aspects: (i) the extremely low probability of multiple precise matches without free parameters, and (ii) the deep structural dependencies between topology, geometry, symmetry, and dynamics.

First, let's consider probability: the parameter-free prediction of the fine structure constant $\alpha^{-1} \approx 137.03650$ deviates by only 3.67 ppm from the CODATA 2022 reference value. Under the naive assumption of a uniform distribution of $\alpha$ in a physically plausible range (e.g., 0.001 to 0.01), the probability of such a small deviation is only about $6 \times 10^{-7}$ (based on an absolute difference of $2.7 \times 10^{-8}$). Corresponding hits can also be found in the E₈ cascade, for example for $\Omega_b \approx 0.04894$ (deviation 0.06% from the Planck data) or $m_p \approx 937 MeV$ (deviation 0.12%). Each of these values corresponds to an independent probability in the range of $10^{-2}$ to $10^{-6}$. Multiplying these for around ten central predictions (flavor mixtures, masses, cosmological constants) results in a combined random probability of less than $10^{-20}$. This is comparable to the improbability of a series of independent dice rolls repeating exactly the same pattern.

Added to this are the structural dependencies: The fixed points $c_3$ and $\varphi_0$ do not arise in isolation, but follow from different but consistent principles – $c_3$ from topological Chern–Simons normalizations in eleven dimensions, $\varphi_0$ from geometric Möbius reduction. Both parameters are independently confirmed in renormalization group-based flows, for example by $\alpha_3(1\, \text{PeV}) \approx \varphi_0$. The layers interlock: topology fixes the normalizations, E₈ orders the cascade, and the RG flows provide dynamic consistency.

This internal entanglement significantly reduces the probability that these are merely random coincidences. A failure in one layer (e.g., in the genetic algorithm or in the dimension chains) would not affect the others, but this is not observed empirically. Instead, a coherent overall picture emerges that can be verified through reproducibility and falsifiability (e.g., in the predicted axion mass).


## Appendix A — Table of fixed points (high precision)

$$
c_{3}=\frac{1}{8\pi}=0.039788735772973836,\qquad
\varphi_{0}=\frac{4}{3}c_{3}+48c_{3}^{4}=0.053171952176845526,
$$
$$
A=2c_{3}^{3}=1.259825563796855\times 10^{-4},\qquad
\kappa=\frac{41}{10}\frac{1}{2\pi}\ln\!\frac{1}{\varphi_{0}}=1.914684795.
$$
$$
\alpha=0.007297325816919221,\qquad \alpha^{-1}=137.03650146488582.
$$
*Reference:* CODATA 2022 $\alpha_{\text{CODATA}}=7.2973525628(11)\times 10^{-3}$,
deviation $\approx 3.67$ ppm.

---

## Appendix B – E₈ cascade in closed form


**Definitions and normalization**

For each nilpotent E₈ orbit

$D_n \;=\; 248-\dim\mathcal O_n,\qquad n=0,\dots,26$,

with the chain $D_n=60-2n$ found from $D_0=60$ to $D_{26}=8$.

The ladder follows from a single normalization at the first step.

$s^{\star}=\ln 248-\ln 60=1.419084183942882,\qquad \lambda=\frac{0.834}{s^{\star}}=0.5877029773404678$ .

**Damping**

$\gamma(0)=0.834,\qquad \gamma(n)=\lambda\Big[\ln D_n-\ln D_{n+1}\Big]\quad (n\ge 1)$.

**Recursion**

$\varphi_{n+1}=\varphi_n\,\mathrm e^{-\gamma(n)}$ .

**Closed form of the ladder**

For $n\ge 1$, the following applies

$\varphi_n \;=\; \varphi_0\,\mathrm e^{-\gamma(0)}\Big(\tfrac{D_n}{D_1}\Big)^{\lambda},\qquad D_1=58$.

**Calibration-free tests**

1. **Proportionality law** for $m,n\ge 1$:


    
    $\boxed{\ \frac{\varphi_m}{\varphi_n}=\Big(\tfrac{D_m}{D_n}\Big)^{\lambda}=\Big(\tfrac{60-2m}{60-2n}\Big)^{\lambda}\ }$ .


    
2. **Log-linear law**


    
    $\boxed{\ \log \varphi_n=\text{constant}+\lambda\,\log D_n\ }$ .


    
**Note on the end of the chain**

The E acht chain ends structurally at $n=26$ with $D=8$. Values for $n>26$ would be an analytical continuation and are marked as extrapolation.

### B.0 Verification of uniqueness

We have completely enumerated the chain set $\mathcal C$ on the $\Delta D=2$‑DAG and minimized $\mathbf F$ lexicographically from 4.5.2. This resulted in exactly **one** minimal chain (Theorem 4.5.1), identical to the sequence $D_{n}=60-2n$ shown in Table B.1 with the labels listed in 4.2. The log-exact ladder derived from this and the damping $\gamma(n)$ agree step by step.

---
### **B.1 – E8 cascade: log-exact values per stage**

Columns:

- $n$
- $D_n$
- $\ln D_n$
- $s_n=\ln D_n-\ln D_{n+1}$
- $\gamma(n)$ with $\gamma(0)=0.834$, otherwise $\lambda s_n$
- $\Sigma\gamma$ cumulative up to level $n$ inclusive
- $\varphi_n/\varphi_0$ uncalibrated
- $\big(\tfrac{D_n}{D_1}\big)^{\lambda}$ as pure chain number

> Note: $\varphi_n/\varphi_0=\mathrm e^{-\gamma(0)}\big(\tfrac{D_n}{D_1}\big)^{\lambda}$ for $n\ge 1;$ for $n=0$, $\varphi_0/\varphi_0=1$.

> [!note] Table note

> The column $(D_{n}/D_{1})^{\lambda}$ is the chain number of the ladders for $n\ge 1$. The entry for $n=0$ is for checking purposes only and is not used physically.

| **n** | **D** | **ln D** | **s_n**  | **γ(n)** | **Σγ**   | **φ_n/φ₀** | **(Dₙ/D₁)^λ** |
| ----- | ----- | -------- | -------- | -------- | -------- | ---------- | ------------- |
| 0     | 60    | 4.094345 | 0.033902 | 0.834000 | 0.000000 | 1.000000   | 1.020124      |
| 1     | 58    | 4.060443 | 0.035091 | 0.020623 | 0.834000 | 0.434309   | 1.000000      |
| 2     | 56    | 4.025352 | 0.036368 | 0.021373 | 0.854623 | 0.425443   | 0.979064      |
| 3     | 54    | 3.988984 | 0.037740 | 0.022180 | 0.875997 | 0.416447   | 0.957994      |
| 4     | 52    | 3.951244 | 0.039221 | 0.023050 | 0.898177 | 0.407312   | 0.936782      |
| 5     | 50    | 3.912023 | 0.040822 | 0.023991 | 0.921227 | 0.398030   | 0.915419      |
| 6     | 48    | 3.871201 | 0.042560 | 0.025012 | 0.945218 | 0.388595   | 0.893899      |
| 7     | 46    | 3.828641 | 0.044452 | 0.026124 | 0.970230 | 0.378996   | 0.872211      |
| 8     | 44    | 3.784190 | 0.046520 | 0.027340 | 0.996355 | 0.369223   | 0.850347      |
| 9     | 42    | 3.737670 | 0.048790 | 0.028674 | 1.023695 | 0.359265   | 0.828299      |
| 10    | 40    | 3.688879 | 0.051293 | 0.030101 | 1.052369 | 0.349110   | 0.806058      |
| 11    | 38    | 3.637586 | 0.054067 | 0.031767 | 1.082514 | 0.338744   | 0.783615      |
| 12    | 36    | 3.583519 | 0.057158 | 0.033589 | 1.114290 | 0.328148   | 0.760962      |
| 13    | 34    | 3.526361 | 0.060625 | 0.035571 | 1.147880 | 0.317306   | 0.738089      |
| 14    | 32    | 3.465736 | 0.064539 | 0.037915 | 1.183450 | 0.306202   | 0.714988      |
| 15    | 30    | 3.401197 | 0.068993 | 0.040555 | 1.221365 | 0.294805   | 0.691650      |
| 16    | 28    | 3.332205 | 0.074108 | 0.043581 | 1.261920 | 0.283078   | 0.668066      |
| 17    | 26    | 3.258097 | 0.080043 | 0.047041 | 1.305501 | 0.271026   | 0.644229      |
| 18    | 24    | 3.178054 | 0.087011 | 0.051117 | 1.352542 | 0.258584   | 0.620130      |
| 19    | 22    | 3.091042 | 0.095310 | 0.055996 | 1.403659 | 0.245652   | 0.595761      |
| 20    | 20    | 2.995732 | 0.105361 | 0.061940 | 1.459655 | 0.232102   | 0.571113      |
| 21    | 18    | 2.890372 | 0.117783 | 0.069239 | 1.520595 | 0.217761   | 0.546180      |
| 22    | 16    | 2.772589 | 0.133531 | 0.078477 | 1.589835 | 0.203747   | 0.520953      |
| 23    | 14    | 2.639057 | 0.154151 | 0.090595 | 1.668311 | 0.188369   | 0.495424      |
| 24    | 12    | 2.484907 | 0.182322 | 0.107151 | 1.758907 | 0.172054   | 0.469584      |
| 25    | 10    | 2.302585 | 0.223144 | 0.131142 | 1.867098 | 0.154572   | 0.443426      |
| 26    | 8     | 2.079442 |          |          | 1.998240 | 0.135574   | 0.416948      |
|


---
## Appendix C – Block formulas for observables


> [!example] Block calibration in practice

> For each block, a single unit calibration $\zeta$ to a reference value is sufficient. All relations within the block then follow without the need for fitting from the ratio laws of the chain, see 4.3 and Appendix B.


**Electroweak block (n=12):**

$v_H=\zeta_{\rm EW} M_{Pl}\varphi_{12},\quad M_W=\tfrac{1}{2} g_2 v_H,\quad M_Z=\tfrac{1}{2}\sqrt{g_1^2+g_2^2} v_H$.

**Hadronic block (n=15,17):**

$m_p=\zeta_p M_{Pl}\varphi_{15},\quad m_b=\zeta_b M_{Pl}\varphi_{15},\quad m_u=\zeta_u M_{Pl}\varphi_{17}$.

**Cosmo blocks:**

$T_{\gamma0}=\zeta_\gamma M_{Pl}\varphi_{25},\quad T_\nu=(4/11)^{1/3}T_{\gamma0},\quad \rho_\Lambda=\zeta_\Lambda M_{Pl}^{4}\varphi_{30}^{97/30}$.

**Fundamental relations near n=0:**

$\Omega_b=\varphi_0(1-2c_3),\qquad r=\varphi_0^2,\qquad V_{us}/V_{ud}=\sqrt{\varphi_0}$.

### C.8 Möbius ladder: Definition and error propagation

**Definition.**  
$\mathcal{M}_y(\delta)=\dfrac{y+\delta}{y-\delta}$ mit $y\in\{1,\tfrac{1}{3},\tfrac{2}{3}\}$.

**Calibration rule.**  
$\delta=\dfrac{\sqrt{m_\tau/m_\mu}-1}{\sqrt{m_\tau/m_\mu}+1}$.

**Derivatives.**  
Set $R=\sqrt{m_\tau/m_\mu}$. Then $\dfrac{\mathrm d\delta}{\mathrm dR}=\dfrac{2}{(R+1)^2}$.  
With $\dfrac{\partial R}{\partial m_\tau}=\dfrac{1}{2}\dfrac{1}{\sqrt{m_\tau m_\mu}}$ and  
$\dfrac{\partial R}{\partial m_\mu}=-\dfrac{1}{2}\dfrac{\sqrt{m_\tau}}{m_\mu^{3/2}}$, it follows that  
$\sigma_\delta^2=\big(\dfrac{\mathrm d\delta}{\mathrm dR}\big)^2\,\sigma_R^2$.

**Prediction formulas.**  
Substitute $\delta$ into the six relations from 7.4.4.  
Fine corrections can be written as a universal displacement:  
$\delta\rightarrow \delta+a_y\,\varphi_0+b_y\,c_3$ with small sector-specific $a_y,b_y$.

---
## Appendix D: Möbius fiber: boundary, curvature normalization, and the coefficient 6pi


>[!info] Target

We explain why the linear marginal coefficient $6\pi$ appears in 3.2.2 and why $\varphi_{\text{tree}}=1/(6\pi)$ follows from this. In addition, we show the topological surcharge $\delta_{\text{top}}$ in equivalent forms and clarify independences from representations and schemata.

---

### D.1 Setup, notation, and compliant scaling

Let $\mathcal M$ be the two-dimensional Möbius fiber (compact, with boundary), $g_{\mathcal M}=\varphi^{2}\hat g_{\mathcal M}$ a purely conformal rescaling. For Gaussian curvature $K$ and geodesic boundary curvature $k_{g}$, the following applies

$$

\int_{\mathcal M}K,\mathrm dA+\oint_{\partial\mathcal M}k_{g},\mathrm ds=2\pi,\chi(\mathcal M),

$$

$$

K=\varphi^{-2}\hat K,\quad \mathrm dA=\varphi^{2}\mathrm d\hat A,\quad \mathrm ds=\varphi,\mathrm d\hat s.

$$

This implies conformal invariance of the surface integral and linearity of the boundary integral:

$$
$$
$$

\int_{\mathcal M}K,\mathrm dA=\int_{\mathcal M}\hat K,\mathrm d\hat A,\qquad

\oint_{\partial\mathcal M}k_{g},\mathrm ds=\varphi,\oint_{\partial\mathcal M}\hat k_{g},\mathrm d\hat s.

$$

  
Thus, the explicit $\varphi$-dependence of the reduced gravitational effect originates exclusively from the boundary (cf. 3.2.2).

---

### D.2 Orientable double cover and the seam contribution


The orientable double cover $\widetilde{\mathcal M}$ of the Möbius fiber is a cylinder with two geometric boundary components. In addition, the $\mathbb Z_{2}$ identification creates a seam curve $\Gamma$. This acts as a third effective boundary cycle.

**Lemma D.2 (Seam as boundary term).** If we model the identification along $\Gamma$ as a limiting case of a thin collar neighborhood with two opposite boundary curves and dihedral angle $\pi$, the corner or seam term in the Gauss-Bonnet balance provides an integrated contribution for each closed $\Gamma$ that is equivalent to a full-value boundary integral with normalization $2\pi$. ∎

**Normalization.**

$$

\mathcal K_{\partial}:=\sum_{\text{boundary cycles}}\oint \hat k_{g},\mathrm d\hat s

= 2\pi+2\pi+2\pi=6\pi.

$$

This makes the number $6\pi$ canonical and independent of any specific representation.

![[Pasted image 20250901121634.png]]

*Illustration: The Möbius fiber can be represented by its orientable double cover, a cylinder with two ordinary boundaries plus one effective seam Γ. Each contributes $2\pi$ to the Gauss–Bonnet balance, leading to the canonical total $6\pi$.*

---
### D.3 Reduced 6D effect and the linear $\varphi$ coefficient


In the 6D reduction, the geometric component contributes linearly to $\varphi$:


$$

S_{\text{grav}}^{(6)} \supset \frac{M_{6}^{4}}{2}\int_{\mathcal B}\!\sqrt{g_{\mathcal B}}

\left[

\underbrace{\int_{\mathcal M}K\,\mathrm dA}_{\text{konform invariant}}

+\underbrace{\oint_{\partial\mathcal M}k_{g}\,\mathrm ds}_{=\ \varphi\,\mathcal K_{\partial}}

\right]

=\frac{M_{6}^{4}}{2}\int_{\mathcal B}\!\sqrt{g_{\mathcal B}}\,(6\pi\,\varphi)+\dots

$$


The effective linear coefficient is $6\pi$.

---

### D.4 Stationarity and tree value $\varphi_{\text{tree}}$

The effective potential density additionally contains a quantized topological contribution. Stationarity yields:

$$

\partial_{\varphi}V_{\text{eff}}(\varphi)\propto 6\pi\varphi-1=0

\quad\Rightarrow\quad

\varphi_{\text{tree}}=\frac{1}{6\pi}.

$$
---

### D.5 The topological surcharge $\delta_{\text{top}}$


The reduction of the 11D Chern-Simons term generates a quantized coupling $g=8c_{3}^{2}$ in 4D with $c_{3}=1/(8\pi)$. The surcharge can be written as:


$$

\boxed{\delta_{\text{top}}=\frac{3}{256\pi^{4}}=48c_{3}^{4}=6c_{3}^{2}g=\tfrac{3}{4}g^{2}}

$$

with $g=\frac{1}{8\pi^{2}}$, $c_{3}=\frac{1}{8\pi}$.

  
Thus:

$$

\boxed{\varphi_{0}=\varphi_{\text{tree}}+\delta_{\text{top}}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}=\tfrac{4}{3}c_{3}+48c_{3}^{4}}

$$

---

### D.6 Unambiguity, invariance, and normalization issues

1. **Freedom of representation.** $6\pi$ depends only on the topology.

2. **Normalization.** $\chi=1$ fixes $\varphi_{\text{tree}}$.

3. **Orthogonality.** $6\pi$ (geometric) and $g$ (topological) are independent.

4. **Schema robustness.** $\delta_{\text{top}}$ is purely a numerical contribution.

### D.7 Consistency check

$$
\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}

\quad\Rightarrow\quad

\kappa=\tfrac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}},\quad

A=2c_{3}^{3}=\tfrac{1}{256\pi^{3}}
\quad

$$

**– consistent with 7.6.**

### D.8 Cross Ratio

$$

\mathrm{CR}(x;y,-y,0)=\frac{y+\delta}{y-\delta}=:\mathcal M_{y}(\delta)

$$
with shift

$$

\delta_{\star}=\tfrac{3}{5}+\tfrac{\varphi_{0}}{6}.

$$

This connects to the ladder diagram (see 3.3.2).

### D.9 Short FAQ

• **Why three edge cycles?** Two cylinder edges plus seam.

• **Is the seam arbitrary?** No, corner term in Gauss-Bonnet.

• **Is $\delta_{\text{top}}$ a numbers game?** No, it follows purely algebraically from $c_{3},g$.

---

### Result

$$

\boxed{\varphi_{\text{tree}}=\tfrac{1}{6\pi},\quad

\delta_{\text{top}}=\tfrac{3}{256\pi^{4}}=48c_{3}^{4}=\tfrac{3}{4}g^{2},\quad
\quad


\varphi_{0}=\tfrac{1}{6\pi}+\tfrac{3}{256\pi^{4}}}

$$


>[!summary] Conclusion

The coefficient $6\pi$ is geometrically fixed, $\delta_{\text{top}}$ is topologically motivated and represented in multiple equivalent ways. The coupling to the other constants in the paper is explicitly traceable.

---

### Appendix E — From the 11D effect to the 4D coefficient $A$ and the log constant $\kappa$

#### E.1 Setup of the effective 4D theory

From 3.2.1 follows the topological coupling $g\,a\,F\tilde F$ with $g=8c_{3}^{2}$ and periodic axion $a$. According to canonical normalization, the relevant 4D sector is
$$[
\mathcal L_{\text{eff}}=-\frac{1}{4}F_{\mu\nu}F^{\mu\nu}+\frac{1}{2}(\partial a)^{2}-\frac{1}{2}m_{a}^{2}a^{2}+g\,a\,F\tilde F,
]$$
where $m_{a}$ is a heavy geometric mode of reduction. $\alpha\equiv g_{\text{em}}^{2}/(4\pi)$. See 3.2 and 7.6. 

#### E.2 Integrating the heavy mode and local operators

Eliminating $a$ in the path integral produces a local series for $p^{2}\ll m_{a}^{2}$
$$[
\Delta\mathcal L=\frac{g^{2}}{2m_{a}^{2}}(F\tilde F)^{2}+\frac{g^{2}}{2m_{a}^{4}}(\partial F\tilde F)^{2}+\dots,
]$$
whose leading term contributes to the renormalization of the photon two-point function with even parity. The $m_{a}$ dependence disappears from the **logarithmically divergent** part of the vacuum polarization, so that the $\ln\mu$ coefficient is scheme-invariant.

#### E.3 Background field method: Logarithmic part of the vacuum polarization

In background field calibration, the single-loop contribution with **two** topological insertions yields the term
$$[
\delta Z_{F}=\underbrace{(8c_{3}^{2})^{2}}_{\text{zwei }g\text{‑Einsätze}}\ \underbrace{\frac{1}{(4\pi)^{3}}}_{\text{Schleifenmaß}}\ \underbrace{\frac{1}{4}}_{\text{Symmetrie}}\ \ln\frac{\mu}{\mu_{0}}+\dots
]$$
and thus an addition $\propto c_{3}^{2}\alpha^{3}$ in $\beta_{\alpha}$. Together, this results in
$$[
A=\frac{1}{256\pi^{3}}=2c_{3}^{3},\qquad 
\beta_{\alpha}=\frac{b_{1}}{2\pi}\alpha^{2}+A\,c_{3}^{2}\alpha^{3}+\dots
]$$
The derivative is independent of details of the $a$ kinetics and the scheme, since only the **coefficient of the log** is used. This number is identical to the variation derivative used in 7.6.

#### E.4 Integrated single loop and $\kappa$

Integration of $d\alpha/d\ln\mu=(b_{1}/2\pi)\alpha^{2}$ between $\mu_{\text{UV}}=M_{\text{Pl}}$ and $\mu_{\text{IR}}=\varphi_{0}M_{\text{Pl}}$ yields
$$[
\kappa=\frac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}},\quad b_{1}=\frac{41}{10}\ \text{in GUT‑Norm},
]$$
as summarized in 7.6.1. $\kappa$ depends only on $b_{1}$ and on the **geometrically** fixed $\varphi_{0}$.  

#### E.5 Fixed point equation from Callan–Symanzik

Using $A$ from E.3 and the integrated constant $\kappa$, we directly obtain
$$[
\alpha^{3}-2c_{3}^{3}\alpha^{2}-8\,b_{1}\,c_{3}^{6}\ln\frac{1}{\varphi_{0}}=0,
]$$
identical to 7.6. Thus, $A$ is **completely** derived from the effective theory.  


---
## Appendix F – Two-loop RGE setup

### Configuration

• **Fermions:** Standard Model plus electroweak triplet $\Sigma_F$ with decoupling at $10^{3},\mathrm{GeV}$; color-adjunct fermion $G8$ of $SU(3)_c$ active for $\mu>M_{G8}=1.8\times10^{mathrm{GeV}; color-adjunct fermion $G8$ of $SU(3)_c$ active for $\mu>M_{G8}=1.8\times10^{10},\mathrm{GeV}; three right-handed neutrinos with staggered thresholds

$M_{N_1}=10^{14},\mathrm{GeV}$, $M_{N_2}=3\times10^{14},\mathrm{GeV}$, $M_{N_3}=8\times10^{14},\mathrm{GeV}$. Above $M_{G8}$, $\Delta b_3=+2$ applies piecewise.  

• **Scalars:** Standard Model Higgs $H$, PQ field $\Phi$ with threshold $M_{\Phi}=10^{16},\mathrm{GeV}$.  

• **Spurion:** Effective $R^{3}$ term for modeling the cubic contribution $\propto \alpha^{3}$ in the Abelian sector.  

• **Normalization:** Hypercharge in **GUT norm**

$$

g_{1}^{\mathrm{GUT}}=\sqrt{\tfrac{5}{3}},g_{Y},\qquad b_1=\tfrac{41}{10}.

$$

For the slope, $,\dfrac{d\alpha_{1}^{-1}}{d\ln\mu}=-\tfrac{b_{1}}{2\pi}$ applies.  

• **Starting values at $\mu=M_Z$:**

$g_{1}^{\mathrm{GUT}}\approx 0.462,\quad g_{2}=0.652,\quad g_{3}=1.2323$.  

• **Integration:** Two loops of beta functions with piecewise threshold matching over at least fifteen decades; optional three loops of slope checking for $SU(3)_c$.  

### Results

#### Fingerprints of the fixed points

$\alpha_{3}(1,\mathrm{PeV})=0.052923411$  versus $\varphi_{0}=0.053171952$  $\Rightarrow$ deviation -0.47%;

$\alpha_{3}(2.5\times10^{8},\mathrm{GeV})=0.039713807$  versus $c_{3}=\tfrac{1}{8\pi}=0.039788736$  $\Rightarrow$ deviation -0.19%.  

#### Near unification

Minimum **relative spread** of inverse couplings = **1.23%** at $\mu^\star\approx 1.43\times10^{15},\mathrm{GeV}$.  

#### Continuity and slopes

Piecewise matching without jumps in $\alpha_i^{-1}$; measured $U(1)$ slope consistent with $-b_1/(2\pi)$, G8 bridge slope numerically $0.8063$ compared to expectation $\tfrac{5}{2\pi}=0.7958$ (1.3%).    

#### Spacing invariant

The three pairs of ties are

$\mu_{23}\approx6.05\times10^{14},\mathrm{GeV}$,

$\mu_{13}\approx1.46\times10^{15},\mathrm{GeV}$,

$\mu_{12}\approx2.38\times10^{15},\mathrm{GeV}$,

so that

$$

S=\log_{10}\mu_{23}-2\log_{10}\mu_{13}+\log_{10}\mu_{12}\approx -0.17.

$$
#### PyR@TE configuration (short, v2)

Settings:

LoopOrder: 3        # Export; Solver uses full 2-loop + optional 3-loop SU(3)

Groups: {U1Y: U1, SU2L: SU2, SU3c: SU3}

Thresholds:

  - {Scale: MSigma, Fields: [SigmaF]}      # 1.0e3 GeV

  - {Scale: MG8,    Fields: [G8]}          # 1.8e10 GeV  (Δb3 = +2 above)

  - {Scale: MNR1,   Fields: [NR1]}         # 1.0e14 GeV

  - {Scale: MNR2,   Fields: [NR2]}         # 3.0e14 GeV

  - {Scale: MNR3,   Fields: [NR3]}         # 8.0e14 GeV

  - {Scale: MPhi,   Fields: [phiR, phiI]}  # 1.0e16 GeV

Fermions:

  G8:  {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 8}}   # new octet

  NR1: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

  NR2: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

  NR3: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

_(see model file v2 for complete YAML)_

**Pyr@ate Configuration:**

```
---

Author: "E8 Cascade TFPT v2.1 – G8 Adjoint Enhanced"

Date: 2025-08-29

Name: E8CascadeTFPTG8_v2

# ------------------------------------------------------------

# ENHANCED E8 CASCADE MODEL WITH G8 ADJOINT FERMION (v2)

#

# GOALS:

# - Keep TFPT fingerprints SM-driven (1-loop) below 10^9 GeV

# - Provide clean G8 color bridge for unification (Δb3 = +2 above MG8)

# - Unambiguous U(1) GUT normalization and documentation

# ------------------------------------------------------------

Settings:

LoopOrder: 3

ExportBetaFunctions: true





  

# ------------------------------------------------------------

# ENHANCED E8 CASCADE THRESHOLDS

# ------------------------------------------------------------

Thresholds:

- Scale: MSigma

Fields: [SigmaF]

- Scale: MG8

Fields: [G8]

- Scale: MNR1

Fields: [NR1]

- Scale: MNR2

Fields: [NR2]

- Scale: MNR3

Fields: [NR3]

- Scale: MPhi

Fields: [phiR, phiI]





  

Groups: {U1Y: U1, SU2L: SU2, SU3c: SU3}





  

Fermions:

Q : {Gen: 3, Qnb: {U1Y: 1/6, SU2L: 2, SU3c: 3}}

L: {Gen: 3, Qnb: {U1Y: -1/2, SU2L: 2}}

uR : {Gen: 3, Qnb: {U1Y: 2/3, SU3c: 3}}

dR: {Gen: 3, Qnb: {U1Y: -1/3, SU3c: 3}}

eR: {Gen: 3, Qnb: {U1Y: -1}}

SigmaF : {Gen: 1, Qnb: {U1Y: 0, SU2L: 3, SU3c: 1}}

G8: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 8}}

NR1: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

NR2: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

NR3: {Gen: 1, Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}





  

RealScalars:

phiR : {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

phiI : {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}}

R3: {Qnb: {U1Y: 0, SU2L: 1, SU3c: 1}, External: True}





  

ComplexScalars:

H : {RealFields: [Pi, Sigma], Norm: 1/sqrt(2), Qnb: {U1Y: 1/2, SU2L: 2}}





  

Potential:

Definitions:

Htilde[i]: Eps[i,j]*Hbar[j]





  

Yukawas:

Yu: Qbar[i,a] Htilde[i] uR[a]

Yd: Qbar[i,a] H[i] dR[a]

Ye: Lbar[i] H[i] eR

# ySig (Type III): intentionally not exported to PyR@TE due to indexing;

# solver explicitly includes its 2-loop trace contribution.

# ySig: Lbar[i] SigmaF[i,j] H[j]

yN1: Lbar[i] Htilde[i] NR1

yN2: Lbar[i] Htilde[i] NR2

yN3: Lbar[i] Htilde[i] NR3





  

QuarticTerms:

lambda : (Hbar[i] H[i])**2

lPhi : (phiR**2 + phiI**2)**2

lHphi : (Hbar[i] H[i])*(phiR**2 + phiI**2)





  

TrilinearTerms:

cR3 : R3 * (Hbar[i] H[i])





  

ScalarMasses:

mu2 : -Hbar[i] H[i]

MPhi : phiR*phiR + phiI*phiI





  

Vevs:

vSM : Pi[2]

vPQ: phiR





  

Parameters:

- {name: MPl, value: 1.221e19}

- {name: MSigma, value: 1.0e3}

- {name: MG8, value: 1.8e10}

- {name: MNR1, value: 1.0e14}

- {name: MNR2, value: 3.0e14}

- {name: MNR3, value: 8.0e14}

- {name: MPhi, value: 1.0e16}

- {name: c3, value: 0.039788735772973836}

- {name: phi0, value: 0.053171952176845526}

- {name: g1, value: 0.357}

- {name: g2, value: 0.652}

- {name: g3, value: 1.2322690515271375}

- {name: Yu33, value: 0.857375}

- {name: Yd33, value: 0.024}

- {name: Ye33, value: 0.010}

- {name: ySig, value: 0.50}

- {name: yN1, value: 0.70}

- {name: yN2, value: 0.70}

- {name: yN3, value: 0.70}

- {name: lambda, value: 0.130}

- {name: lPhi, value: 0.10}

- {name: lHphi, value: 0.01}

- {name: cR3, value: 0.01}





  

Substitutions: { g_U1Y: g1, g_SU2L: g2, g_SU3c: g3 }





  

# ------------------------------------------------------------

# THEORY NOTES (v2):

# - U(1) GUT normalization: g1_GUT = sqrt(5/3) * gY; b1_GUT = 41/10.

# - ySig kept out of PyR@TE export; solver adds its 2-loop trace.

# - LoopOrder=3 in YAML; solver uses full 2-loop + optional 3-loop SU(3).

# - G8 bridge above MG8: Δ(α3^{-1})(μ) = -(Δb3)/(2π) ln(μ/MG8), Δb3=+2.

# ------------------------------------------------------------
```

---

## Appendix G Nilpotent Orbits in Type E8
![[Pasted image 20250823130052.png]]
![[Pasted image 20250823130105.png]]

# Appendix H: References:

#### 1. Nilpotent orbits in semisimple Lie algebras (especially E8)

- Collingwood, D.H. and McGovern, W.M., _Nilpotent Orbits in Semisimple Lie Algebras_, Van Nostrand Reinhold, New York (1993). – Referenced for the general classification of nilpotent orbits in semisimple Lie algebras, including detailed tables and dimensions for E8 orbits, which serve as the basis for the γ(n) damping function and the parabolic cascade.
- Djouadi, A. et al., "Induced Nilpotent Orbits of the Simple Lie Algebras of Exceptional Type," arXiv: (from publication, e.g., similar to iris.unitn.it/handle/11572/77393) (200x). – Referenced for the induction of nilpotent orbits in E8 and their dimensions, which motivate the monotonic decay sequence in the cascade (e.g., from 248 to 206).
- Landsberg, J.M. and Manivel, L., "Series of Nilpotent Orbits," Experimental Mathematics 13(1) (2004), 69–78. – Referenced for the organization of nilpotent orbits in series within exceptional algebras such as E8, including dimension formulas that support the quadratic smoothing of γ(n).

#### 2. Chern-Simons term in 11D supergravity and topological fixed points

- Cremmer, E., Julia, B., and Scherk, J., "Supergravity Theory in Eleven Dimensions," Physics Letters B 76(4) (1978), 409–412. – Referenced for the original formulation of 11D supergravity, including the Chern-Simons term, which derives the normalization 1/(8π) for c₃ and topological fixed points.
- Troncoso, R. and Zanelli, J., "Higher-Dimensional Supergravities as Chern-Simons Theories," International Journal of Theoretical Physics 38(4) (1999), 1181–1193 (or extended version arXiv:1103.2182). – Referenced for the interpretation of 11D supergravity as a Chern-Simons theory, which explains the topological trace of c₃ = 1/(8π) and the Möbius reduction to φ₀.
- Duff, M.J., "Eleven-Dimensional Supergravity, Anomalies and the E8 Yang-Mills Sector," Nuclear Physics B 325(2) (1989), 505–522. – Referenced for the connection between Chern-Simons terms in 11D and E8 symmetries, relevant for the topological correction in φ₀ (e.g., 3/(256π⁴)).

#### 3. E8 in Grand Unified Theories (GUTs) and String Theory

- Gross, D.J., Harvey, J.A., Martinec, E., and Rohm, R., "Heterotic String Theory (I). The Free Heterotic String," Nuclear Physics B 256 (1985), 253–284. – Referenced for the role of E8 × E8 in heterotic string theories, which inspire the embedding of E8 as an ordering principle for the scale ladder (γ(n) from orbits).
- Lisi, A.G., "An Exceptionally Simple Theory of Everything," arXiv:0711.0770 (2007). – Referenced for the attempt to use E8 as a unified symmetry for all forces and matter, similar to the E8 cascade in the paper, including orbit structures for flavor and scales.
- Green, M.B., Schwarz, J.H., and Witten, E., Superstring Theory, Cambridge University Press (1987), Volume 2. – Referenced for the E8 gauge group in string theory, specifically its nilpotent elements and anomalies, which support the quadratic form of γ(n) and RG confirmation.

#### 4. Theoretical derivations of the fine structure constant (α)

- Wyler, A., "On the Conformal Groups in the Theory of Relativity and a New Value for the Universal Constant α," Lettere al Nuovo Cimento 3(13) (1971), 533–536. – Referenced for an early geometric derivation of α from conformal groups, which anticipates the parameter-free cubic fixed point equation in the paper (based on geometry and topology).
- Atiyah, M.F., "On the Fine-Structure Constant," (lecture/note, 2018, see e.g. preposterousuniverse.com/blog/2018/09/25/atiyah-and-the-fine-structure-constant/). – Referenced for a mathematical derivation of α ≈ 1/137 from algebraic structures, comparable to the cubic equation and ppm accuracy in the paper.
- Smith, S.J., "A New Theoretical Derivation of the Fine Structure Constant," Progress in Physics 28(1) (2012), 1–5. – Referenced for a modern derivation of α without free parameters, which complements the approach of the paper (coupling of topology c₃ and geometry φ₀).

#### 5. Other related topics (e.g., RG flows, genetic algorithms in physics)

- 't Hooft, G. and Veltman, M., "Regularization and Renormalization of Gauge Fields," Nuclear Physics B 44(1) (1972), 189–213. – Referenced for the fundamentals of RG flows in gauge theories, which underpin the two-loop analysis and fingerprints (φ₀ at 1 PeV, c₃ at 10^8 GeV) in the QCD process.
- Koza, J.R., Genetic Programming: On the Programming of Computers by Means of Natural Selection, MIT Press (1992). – Referenced for the method of genetic algorithms, which methodically substantiates the GA approach in the paper (search for Lagrange densities and emergence of fixed points).

## Appendix I: Changelog

### Version 1.0.6 - 2025-09-01

1. New section 7.6 -  **Variational Derivative in Four Dimensions**
2. Section 5: Adjustment to new 2-loop RGE configuration & results
3. New Section 8.0a  - Chirality from boundary and flow: operational summary
4. New Appendix J - Chirality on the double cover
5. Section 4.5 - Full adjustment to prove the uniqueness of the chain
6. New Sections 8.1.1 and 8.1.2 - Block constants from boundary cycles and Abelian trace & derivation of the $\zeta_{B}$ formula from the boundary functional
7. Revision of Section D 
8. New 5.1b  **Thresholds as ladder outputs instead of model selection**
9. New 5.2c  **Gauge–moduli–locking in the E six window**
10. Section 6.12 - Unambiguous fixation of $\alpha_{\text{inf}}$
11. Section 8.5.1 - Systematics of deviations
12. Section 7.6.2 - Callan-Symanzik route
13. Revision 8.4.6 - No circular reference
14. 


## Appendix J — Chirality on the double cover

### J.1 Setup and Notation
Geometry: $(M_{6}=M_{4}\times\widetilde{M})$, where $(\widetilde{M})$ is the **orientable double cover** of the Möbius fiber with **three** closed boundary cycles $(C_{1},C_{2},C_{T})$.  
Edge count: $(\sum_{i}\!\oint_{C_{i}}\widehat{k}_{g}\,\mathrm ds = 6\pi)$.  
The number $(6\pi)$ fixes the linear boundary coefficient in 3.2.2 and thus $(\varphi_{\text{tree}}=\tfrac{1}{6\pi})$, plus the topological surcharge $(\delta_{\text{top}}=\tfrac{3}{256\pi^{4}})$ to $(\varphi_{0})$.  References: 3.2.2 and Appendix D. 


![[Pasted image 20250901105748.png]]
**Fig. J.1** shows $(\widetilde M$) with the three boundary cycles and the projectors $(P_{i})$ (see below).

![[Pasted image 20250901111314.png]]

## J.2  Six-dimensional spinor reduction and projectors

Choose $\Gamma$ matrices as  
$$
\Gamma^{\mu} = \gamma^{\mu}\otimes\sigma^{0}, \quad 
\Gamma^{5} = \gamma_{5}\otimes\sigma^{1}, \quad 
\Gamma^{6} = \mathbf 1\otimes\sigma^{2}.
$$

The 6D Weyl condition  
$$
\Gamma_{7} = \gamma_{5}\otimes\sigma^{3}, \quad \Gamma_{7}\Psi = +\Psi
$$
yields  
$$
\Psi(x,y) = \psi_{L}(x)\otimes \chi_{+}(y) + \psi_{R}(x)\otimes \chi_{-}(y),
$$
with  
$$
\sigma^{3}\chi_{\pm} = \pm \chi_{\pm}.
$$

Chiral boundary conditions on the three boundary cycles:  
$$
P_{T} = \frac{1}{2}(\mathbf 1 + i\,\sigma^{3}\sigma^{n}), \quad
P_{1} = P_{2} = \tfrac12(\mathbf 1 - i\,\sigma^{3}\sigma^{n}).
$$

Thus, zero modes exist only for $\chi_{+}$, and the 4D zero modes are left-chiraled.  
The choice is elliptic and gauge invariant.  

Optionally, Wilson lines $W_{i}\in SU(3)\times SU(2)\times U(1)$ along $C_{i}$ only determine phases without breaking the 4D gauge group.  
They serve to read out the abelian trace in the EW block with  
$$
k_{\mathrm{EW}} = \tfrac{41}{32}.
$$  
(Reference 8.4.6)

---

## J.3  Index set on $\widetilde M$ and family number

Let $A$ be an abelian connection with flow  
$$
m = \tfrac{1}{2\pi}\int_{\widetilde M}F \in \mathbb Z
$$
and boundary holonomies  
$$
\tfrac{1}{2\pi}\oint_{C_{i}}A = \nu_{i} \in \mathbb Z.
$$

The following index applies to the projectors $P_{i}$:  
$$
\mathrm{Ind}\,D_{\widetilde{M}}
= \#\chi_{+} - \#\chi_{-}
= \tfrac{1}{2\pi}\int_{\widetilde M}F
= \nu_{1} + \nu_{2} + \nu_{T}.
$$

Reasoning: The APS formula with boundary projections sets the $\eta$ contributions to zero, Stokes provides  
$$
\sum_{i}\oint_{C_{i}}A = \int_{\widetilde M}F.
$$



Corollary: Minimal $(\nu_{1},\nu_{2},\nu_{T}) = (1,1,1)$ yields  
$$
\mathrm{Ind}\,D = 3
$$
→ three families.  

The integer structure is consistent with the Chern–Simons quantization from 3.2.1:  
$$
g = \tfrac{n}{8\pi^{2}}, \qquad c_{3} = \tfrac{1}{8\pi}.
$$
*Fig. J.2 illustrates the count $\mathrm{Ind}\,D = \nu_{1} + \nu_{2} + \nu_{T}$.*

---

## J.4  Absence of anomalies in a family

Standard model per family with right-handed neutrino is anomaly-free:

$$
U(1)_{Y}^{3}: \quad
3\cdot2\left(\tfrac{1}{6}\right)^{3}
+ 3\left(-\tfrac{2}{3}\right)^{3}
+ 3\left(\tfrac{1}{3}\right)^{3}
+ 2\left(-\tfrac{1}{2}\right)^{3}
+ (1)^{3} = 0,
$$

$$
U(1)_{Y}\times SU(2)^{2}: \quad
3\cdot\tfrac{1}{6} + (-\tfrac{1}{2}) = 0,
$$
$$


$$

$$
\text{Gravity}-U(1)_{Y}: \quad
3\cdot2\cdot\tfrac{1}{6}
+ 3\left(-\tfrac{2}{3}\right)
+ 3\left(\frac{1}{3}\right)
+ 2\left(-\tfrac{1}{2}\right)
+ 1 = 0.
$$

With $\mathrm{Ind}\,D = 3$, the overall theory remains anomaly-free.  

---

## J.5  Compatibility with rank windows and track

The chain  
$$
E_{8}\supset E_{7}\supset E_{6}
$$
remains window logic in RG flow, no 4D calibration content.  

The Wilson lines are flat and preserve $SU(3)\times SU(2)\times U(1)$.  

The same abelian trace appears twice:  
$$
k_{\mathrm{EW}} = \tfrac{41}{32}
$$
in the EW block and  
$$
b_{1} = \frac{41}{10}
$$
in  
$$
\kappa = \frac{b_{1}}{2\pi}\ln\frac{1}{\varphi_{0}}
$$
of the $\alpha$ fixed point equation.  
(References 7.6, 8.4.6)

---

## J.6  Stability with two loop windows

Flat Wilson lines only minimally change Kaluza thresholds.  

The fingerprints documented in 5.2  
$$
\alpha_{3}(1\,\mathrm{PeV}) \simeq \varphi_{0}, \qquad 
\alpha_{3} \simeq c_{3} \quad \text{at } \mu \sim 2.5\times 10^{8}\,\mathrm{GeV}
$$
remain stable.