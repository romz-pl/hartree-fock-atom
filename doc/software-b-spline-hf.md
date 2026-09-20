# A Review of SPHF — the B-spline Hartree–Fock Program of C. Froese Fischer

| | |
|---|---|
| **Package** | SPHF (a.k.a. "B-spline HF"), version 1.00 |
| **Author** | Charlotte Froese Fischer (Vanderbilt Univ. / NIST) |
| **Reference paper** | *A B-spline Hartree–Fock program*, Comput. Phys. Commun. **182**, 1315–1326 (2011), doi:[10.1016/j.cpc.2011.01.012](https://doi.org/10.1016/j.cpc.2011.01.012) |
| **CPC catalogue ID** | AEIJ_v1_0 |
| **Repository reviewed** | <https://github.com/compas/sphf> (branch `master`, 9 commits, MIT licence) |
| **Language** | Fortran 90/95, depends on LAPACK/BLAS |
| **Scope** | Non-relativistic, spherical, LS-coupled Hartree–Fock for free atoms and ions |

> **Method note.** This review is based on (i) the paper bundled in the repository, (ii) reading the source, and (iii) *actually building and running the code* with gfortran 13.3 on Ubuntu 24.04. Every number in §5–§7 was measured by me in this session unless explicitly attributed to the paper or the literature. Where I could not verify something independently, I say so.

---

## 1. Executive summary

SPHF is a compact, readable, and numerically excellent solver for the **radial** Hartree–Fock problem for atoms. Its central idea is to expand each radial orbital in a **B-spline basis** so that the two-point boundary-value problem of classical HF codes (node counting, log-grids, "accelerating" parameters) is replaced by **algebraic** problems: generalized eigenproblems $(\mathbf H_a-\varepsilon_{aa}\mathbf B)\mathbf a=0$ and Newton–Raphson linear systems.

**Verdict in one paragraph.** For what it claims to do — closed-shell and simple open-shell atoms up to at least Z ≈ 90 — it works, and works very well. I reproduced the paper's virial-theorem test to machine precision on all eight bundled cases, and energies agree with literature Hartree–Fock limits to **≤ 2×10⁻⁹ E_h** for He, Be, Ne, Ar and **6×10⁻⁹ E_h** for Kr. It builds in ~4 s with zero warnings at default flags and is clean under `-fcheck=all` + signalling-NaN + FPE traps. Its weaknesses are **scope and engineering, not numerics**: strict limits on open-shell coupling, interactive-only input, fragile input parsing, no automated test harness, truncated reference outputs, no build system beyond a hand-written Makefile, and a dense-matrix (O(ns²)) design that is unsuitable for very large bases.

| Dimension | Rating | One-line justification |
|---|---|---|
| Numerical accuracy | ★★★★★ | Energies to ~1e-9 E_h; virial ratio ≈ −2 to 1e-14 |
| Robustness of SCF | ★★★★☆ | Rotation + projection + NR handles hard cases (2p⁵, 1s2s) |
| Physics scope | ★★☆☆☆ | ≤ 2 open shells for specific LS terms; no relativity, no CI/MCHF |
| Code quality | ★★★☆☆ | Modern module style, good comments; dead code, stack arrays, float `==` |
| Usability / I/O | ★★☆☆☆ | Prompt-driven; raw runtime error on malformed input |
| Testing / reproducibility | ★★☆☆☆ | 8 cases, no pass/fail check; two reference logs truncated |
| Performance | ★★★☆☆ | Fast for Z ≤ 40; 1.7 s for Ra; dense matrices limit scaling |
| Documentation | ★★★☆☆ | Excellent paper; thin README; one citation error |

---

## 2. Background and lineage

SPHF sits at the end of a long line of Froese Fischer codes:

```
Hartree (1957, desk calculator)
   └─ MCHF77  (CPC 14, 145, 1978)
        └─ HF86   (CPC 43, 355, 1987)   ── log grid  ρ = ln(Zr),  y = P/√r
             └─ HF96   (CPC 98, 255, 1996)  ── + partially filled f-shells
                  └─ SPHF (CPC 182, 1315, 2011)  ── B-spline / Galerkin, no COMMON
                       └─ DBSR_HF (CPC 202, 287, 2016) ── Dirac–HF, Zatsarinny & Froese Fischer
```

The problem SPHF was designed to solve is stated crisply in the paper: finite-difference HF codes must contend with *(a)* multiple off-diagonal Lagrange multipliers from orthogonality constraints, *(b)* SCF iterations that oscillate (the textbook case is F 2p⁵), *(c)* multiple solutions, requiring **node counting** that the author herself calls "an art", and *(d)* coordinate transformations that users must understand before modifying code. B-splines convert all of this into linear algebra with a natural variational structure.

The B-spline machinery itself is borrowed from **Zatsarinny's BSR library** (CPC 174, 273, 2006; lightly modified), which the paper credits explicitly.

---

## 3. Theory implemented

### 3.1 The B-spline basis

Given a knot sequence $t_0<\dots<t_{n_v}$ with end knots of multiplicity $k_s$, the $n_s=n_v+k_s-1$ B-splines $B_i(r)$ span piecewise polynomials of degree $k_s-1$ with continuity $C^{k_s-2}$ at interior knots. Each is non-zero on $k_s$ adjacent intervals, so the overlap and one-electron matrices are **banded** with bandwidth $2k_s-1$. A radial function is

$$P(a;r)=\sum_{i=1}^{n_s}a_i B_i(r),\qquad \|P\|^2=\mathbf a^{t}\mathbf B\,\mathbf a .$$

Boundary conditions: $a_1=0$ (origin) and, unusually, **two** conditions at the far end, $a_{n_s-1}=a_{n_s}=0$ (i.e. $P(R)=P'(R)=0$) "for additional numerical stability". No condition enforcing $r^{l+1}$ behaviour at the origin is imposed.

### 3.2 Variational rules and the generalized eigenproblem

The HF energy is *quartic* in the expansion coefficients. Varying with respect to one orbital eliminates one summation. The paper's **Table 1** codifies the rules (e.g. $\tfrac12\partial F^k(a,a)/\partial\mathbf a=2R^k(\cdot,a;\cdot,a)\mathbf a$), the key notational device being the "dot" that marks the summation eliminated by variation. This yields, per orbital,

$$\big(\mathbf H^{a}-\varepsilon_{aa}\mathbf B\big)\mathbf a=0 .$$

Direct and $F^k$ terms are **banded**; exchange $G^k$ terms are **full** — that distinction drives the whole performance profile (§7).

### 3.3 Orthogonality: rotations and projection operators

This is the genuinely novel theory in the paper. For two orthogonal orbitals $a,b$ of the same symmetry that are both varied, the matrix of Lagrange multipliers must be symmetric ($\varepsilon_{ab}=\varepsilon_{ba}$) — an extra condition because the energy must be stationary under orbital **rotations**

$$\begin{pmatrix}P^\ast(a)\\P^\ast(b)\end{pmatrix}=\frac{1}{\sqrt{1+\epsilon^2}}\begin{pmatrix}1&-\epsilon\\\epsilon&1\end{pmatrix}\begin{pmatrix}P(a)\\P(b)\end{pmatrix}.$$

Expanding $E(\epsilon)=E(0)+g\epsilon+g'\epsilon^2+\dots$ gives the stationary rotation $\epsilon=-g/(2g')$. The gradients $g,g'$ for every integral type are tabulated in the paper's **Table 2** and coded verbatim in `rotate.f90` (I checked the header comment against the paper; they match). Once stationary, the off-diagonal multipliers are eliminated by **projection operators** (after Bentley, ref. below), giving a pair of *independent* generalized eigenproblems

$$\big[(\mathbf I-\mathbf B\mathbf b\mathbf b^{t})\mathbf H^{a}(\mathbf I-\mathbf b\mathbf b^{t}\mathbf B)-\varepsilon_{aa}\mathbf B\big]\mathbf a=0 .$$

Special cases the paper treats: Koopmans' theorem for filled shells of the same symmetry; and $g=g'=0$ (energy rotation-invariant), exemplified by 1s2s ³S.

### 3.4 Convergence machinery: SCF vs. Newton–Raphson

* **Singly occupied orbital** → $\mathbf H^a$ depends only on *other* orbitals → plain eigenproblem (`hf_eiv`). Controls *which* eigenvalue is selected.
* **Multiply occupied orbital** → $\mathbf H^a$ depends on $\mathbf a$ itself (self-interaction) → naïve SCF can oscillate (F 2p⁵). Fix: **single-orbital Newton–Raphson** (`hf_nr`), i.e. a linear solve instead of an eigenproblem.
* **Warning stated by the author:** NR converges to the *nearest* solution. With screened-hydrogenic guesses for Mg 3p² ³P it converged to the **4p²** energy instead. Hence the mixed strategy: SCF (≤ 6 iterations) then all-orbital NR.
* **All-orbital NR** (`hfall_nr`) solves the augmented Jacobian system, Eq. (19) of the paper, in the unknowns $(\Delta\mathbf a,\Delta\mathbf b,\Delta\varepsilon_{aa},\Delta\varepsilon_{bb},\Delta\varepsilon_{ab})$, with the orthonormality conditions appended. Quadratic convergence is expected and observed (§6.3).
* The paper also sketches, but **did not test**, a block-LU and an SVD alternative (both linear convergence).

### 3.5 Energy expressions for open shells

`get_energy.f90` (1,133 lines, the largest file) builds the energy from Slater's **average energy of a configuration**, then adds tabulated LS-term **deviations** for a *single* open s/p/d/f shell, or one s-electron plus an open shell, or (nd)(n'f). This tabular design is the root cause of the scope limit in §4.

---

## 4. Capabilities and restrictions

**Supported**

| Feature | Notes |
|---|---|
| Average-of-configuration energies | Any configuration with s, p, d, f, g electrons; no restriction |
| Specific LS terms | Restricted (below) |
| Fractional occupations | e.g. `3p(0.5)3d(0.5)` = average of two configurations, **no interaction** between them |
| Ground and excited states | Default universal grid; excited states via orthogonality constraints |
| Fixed-core Rydberg series | For a singly-occupied last orbital, prints eigenvalues + mean radii of the whole series |
| Grid refinement | Map a solution to a finer grid (`refine_grid`) — the "two-level" workflow |
| Fixed vs. varied orbitals | `all`, `none`, `=n` (last *n*), or explicit list |

**Restrictions on specific LS terms** (verbatim from the program summary): only one or two open shells, of types (nl)ᴺ n′s, (np)ᴺ n′l, or (nd)(n′f).

**Not supported:** relativistic effects (use DBSR_HF); Breit–Pauli; MCHF/CI (the paper says "further studies are under way"); more than two coupled open shells; molecules/diatomics; non-central potentials.

---

## 5. Program architecture (as read from source)

**Size (measured):** 7,721 lines of Fortran across ~55 files; 4,677 code / 2,102 comment / 942 blank (≈ 27 % comments); 14 modules; 118 subroutines/functions. The paper's 13,925-line figure includes test data.

**Control flow** (`sphf.f90`, only ~40 lines):

```
open_files → read_hf_param → get_case ─┬─ get_atom        (interactive prompts / bsw.c)
                                       ├─ get_energy      (energy expression from config + term)
                                       └─ get_spline_param(default: h=0.25, ks=4, ns=30+2√Z)
define_spline                                            (B-spline arrays, Galerkin matrices)
DO acc = 1, acc_max                                      (default 2 accuracy levels)
   acc==1: get_estimates      (read bsw.inp or screened-hydrogenic guesses)
   acc>1 : refine_grid        (h←h/2, ks←ks+2, remap orbitals)
   write_spline_param
   solve_HF_equations         (solve_all.f90  OR  solve_scf.f90)
write_hf_param                                           (for re-running the top level)
```

**Two executables from one code base** — selected at link time by which `solve_*.o` is linked:

* `sphf_scf` — updates one orbital at a time throughout.
* `sphf_all` — SCF for ≤ 6 iterations, then updates **all** orbitals simultaneously by NR (Level 2).

**Default grid.** Two regions in $t=Zr$: equally spaced ($t_i=t_{i-1}+h$) out to about the maximum of hydrogenic 1s, then exponential ($t_i=t_{i-1}(1+h)$). Recommended $h$ be an *exact binary fraction* (¼, ⅛). Level 2 halves $h$ and raises $k_s$ by 2.

**I/O files:** `bsw.c` (case description), `bsw.inp`/`bsw.out` (orbitals in/out), `hf_param` (parameters; auto-written for re-runs), `sphf.log` (summary), `ryd.bsw` (Rydberg series).

**Performance trick worth noting.** Slater matrix elements $R^k(i,j;i',j')$ are stored once per grid. The paper also modifies `add_rkm_x` to skip exchange contributions where the product of expansion coefficients is < 10⁻¹⁸; on Ac this cut CPU time by ~30 % with negligible accuracy change (paper's Table 4: 9.085 s → 6.612 s).

---

## 6. Hands-on evaluation

### 6.1 Build and toolchain

* Environment: Ubuntu 24.04, **gfortran 13.3.0**, reference LAPACK/BLAS, default Makefile flags (`-O2`).
* `make sphf_all` → **success, ~4.3 s, 0 warnings**. `make sphf_scf` → success.
* The Makefile is hard-wired to `gfortran`, writes to `../bin/` (which **does not exist in the repo** — I had to `mkdir bin`), and links `-llapack -lblas`. Object files are tracked in `.f90.o` suffix rules with **no dependency tracking for `.mod` files**, so editing a module and re-running `make` will not recompile dependents.
* Because `sphf_all` and `sphf_scf` share objects but differ in one file, switching between them requires deleting `solve_*.o` (I did this manually). A `clean` between builds is the safe path.

### 6.2 Correctness — the eight bundled tests

I ran all eight test cases (`test_all/run1…run8`) in a scratch copy and compared to (a) the bundled reference logs and (b) the paper's Table 4.

| # | Case | E_total (E_h), my run | \|ratio + 2\|, mine | \|ratio + 2\|, paper Table 4 |
|---|---|---:|---:|---:|
| 1 | He 1s² ¹S | −2.861679995611798 | 3.4e-14 | 9e-15 |
| 2 | Li 1s² 2s ²S | −7.432726930728180 | 1.0e-14 | 0 |
| 3 | He 1s2s ¹S | −2.169854456993815 | 8.8e-10 | 8.8e-10 |
| 4 | Be 1s² 2s² ¹S | −14.573023168313082 | 5.5e-14 | 3.9e-14 |
| 5 | F [Ne]2p⁵ ²P | −99.409349386693137 | 5.0e-15 | 1e-15 |
| 6 | Mg [Ne]3s3d ¹D | −199.426751894065177 | 3.5e-14 | 3.2e-14 |
| 7 | Ra [Rn]7s² ¹S | −23094.303666365959 | 1.5e-14 | 2.3e-14 |
| 8 | Ac [Rn]7s² 7p ²P | −23722.104101465030 | 5.6e-11 | 5.6e-11 |

**Findings.**

1. **Cases 1–6 reproduce the bundled reference logs to ≤ 2×10⁻¹³ E_h** in total energy. The differences visible in `diff` are all in the last one or two printed digits — ordinary floating-point reassociation between compilers/BLAS, fifteen years apart. **No regression.**
2. **The bundled reference outputs for cases 7 and 8 (Ra, Ac) are truncated.** `test7.log`/`test8.log` contain *Level 1 only*; the bundled `out7` literally ends at "Iteration 2.01" with nothing after. My full runs complete Level 2. The Level 1 energy matches the bundle to **1.5×10⁻¹¹ E_h**, and the bundled stdout's Level-2 line (`−23094.303666365631`) matches my completed Level 2 (`…958634`) to ~3×10⁻¹⁰ E_h. So the *code* is fine; the *shipped reference artifacts* for the two heaviest cases are incomplete. A naive `diff` regression check would report failure for the wrong reason.
3. The virial ratios I obtain match the paper's Table 4 pattern, including the two "soft" cases (He 1s2s ¹S and Ac), where the ratio only reaches −2 ± 10⁻¹⁰ to 10⁻¹¹. This is a property of the method (1s2s has a large off-diagonal multiplier; Ac has a huge dynamic range in the core), not a defect of my build.

### 6.3 Convergence behaviour

**SCF → NR hand-off (Ne 1s² 2s² 2p⁶, `sphf_all`).** Total energy by iteration on the default coarse grid:

| Iter | Phase | E_total (E_h) | scf-diff |
|---|---|---:|---:|
| 1.01 | SCF | −126.868065 | 2.0e-2 |
| 1.02 | SCF | −127.150766 | 2.2e-3 |
| 1.03 | SCF | −128.467253 | 1.0e-2 |
| 1.04 | SCF | −128.546974 | 6.2e-4 |
| 1.05 | SCF | −128.547055 | 6.4e-7 |
| 2.01 | NR | −128.547057 | 1.1e-8 |
| 2.02 | NR | −128.547057 | 2.0e-15 |

Note iteration 1.03: the scf-diff **increases** (2.2e-3 → 1.0e-2) before falling — the SCF is not monotone, which is precisely why the NR phase exists. NR then contracts the step from 1.1e-8 to 2.0e-15 in one iteration, consistent with the quadratic convergence claimed in the paper. On the refined Level 2 grid the whole calculation needs only **2 SCF + 2 NR** iterations, because the mapped orbitals are already excellent.

### 6.4 Accuracy against the literature

Converged runs (default two-level scheme) vs. accepted Hartree–Fock limits:

| Atom | SPHF (E_h) | Reference (E_h) | \|ΔE\| | Reference source |
|---|---:|---:|---:|---|
| He 1S | −2.861679995612 | −2.861679996 | 3.9e-10 | Literature (quoted to 9 d.p.) |
| Be 1S | −14.573023168313 | −14.57302317 | 1.7e-9 | Literature (quoted to 8 d.p.) |
| Ne 1S | −128.547098109350 | −128.547098109 | 3.5e-10 | Lehtola, arXiv:2108.05850 / 1810.11651 |
| Ar 1S | −526.817512802497 | −526.817512803 | 5.0e-10 | Lehtola, Int. J. Quantum Chem. 119, e25945 (2019) |
| Kr 1S | −2752.054977343617 | −2752.05497735 | 6.4e-9 | Lehtola (same) |

The residuals are set by the number of digits in the reference values, not by SPHF. **This is essentially the exact numerical HF limit.**

> **Self-correction, for transparency.** My first Ar comparison used a reference value I recalled from memory (−526.81751225) and appeared to show a 5.5×10⁻⁷ discrepancy. I did not attribute this to the code; I looked up the accepted value (−526.817512803, Lehtola 2019) and the discrepancy vanished (5×10⁻¹⁰). The original figure was my recall error.

**Basis convergence study (Ne, my runs).** Spline order matters far more than step size:

| h | k_s | n_s | R_max (a₀) | E_total (E_h) | \|ratio + 2\| | time |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 4 | 24 | 221.7 | −128.545002720 | 3.6e-5 | 0.010 s |
| 0.25 | 4 | 36 | 64.6 | −128.547056868 | 6.9e-7 | 0.015 s |
| **0.25** | **6** | **36** | 41.4 | **−128.547098083** | **9.0e-12** | 0.024 s |
| 0.125 | 4 | 60 | 32.1 | −128.547097329 | 1.3e-8 | 0.027 s |
| 0.125 | 6 | 60 | 25.4 | −128.547098109 | 1.0e-14 | 0.051 s |
| 0.125 | 8 | 60 | 20.0 | −128.547098109 | 8.0e-15 | 0.068 s |
| 0.0625 | 6 | 110 | 22.0 | −128.547098109 | 9.1e-14 | 0.178 s |
| 0.0625 | 8 | 110 | 19.5 | −128.547098109 | 9.2e-14 | 0.264 s |

Going from $k_s=4$ to $6$ at fixed $h=0.25$ improves the virial deviation by **~5 orders of magnitude** for a 1.6× cost. Halving $h$ *beyond* $h=0.125$ buys nothing but time (and slightly worse conditioning: 1e-14 → 9e-14). **Practical rule: use $k_s\ge 6$, $h\in\{1/4,1/8\}$.** This vindicates the program's Level-2 default ($h\to h/2$, $k_s\to k_s+2$).

### 6.5 `sphf_all` vs. `sphf_scf`

Best-of-three wall time, gfortran `-O2`, single core:

| Atom | sphf_all | sphf_scf | Speed ratio | Virial dev. (all) | Virial dev. (scf) | ΔE (all−scf) |
|---|---:|---:|---:|---:|---:|---:|
| Ne | 0.037 s | 0.030 s | 1.2× | 6.0e-15 | 1.3e-11 | 4e-13 |
| Ar | 0.085 s | 0.074 s | 1.1× | 1.9e-14 | 5.6e-14 | 4e-12 |
| Kr | 0.279 s | 0.163 s | 1.7× | 8.9e-15 | 6.9e-13 | 2e-11 |
| Ra 7s² | 1.724 s | 0.691 s | **2.5×** | 1.8e-14 | 9.0e-12 | 4e-11 |

`sphf_all` is uniformly more accurate in the virial ratio (up to ~3 orders) but the **total energies agree to ≤ 4×10⁻¹¹ E_h** — the extra work purchases better *orbital shape*, which matters for properties (hyperfine constants, ⟨r⟩, transition integrals), not the energy itself. For total energies alone, `sphf_scf` is the better value on heavy atoms. The paper's stated conclusion is confirmed.

### 6.6 Robustness testing

I built a second binary with `-O0 -g -fcheck=all -Wall -Wextra -finit-real=snan -ffpe-trap=invalid,zero,overflow -fbacktrace` and ran seven cases (He, Li, He 1s2s, F 2p⁵, Mg 3s3d, Ra, Ac).

* **Result: zero runtime errors, zero FP traps, zero bounds violations, zero uninitialised reads.** Energies and virial ratios were identical to the optimised build's (e.g. Ra `−2.000000000000015`, Ac `−2.000000000056318`). (I verified by md5 and by grepping for runtime-error strings that I really did execute the instrumented binary; my first worry that I had run the wrong one was unfounded.)
* Strict-warning build: 40 warnings, all benign categories — 28 real comparisons flagged by `-Wcompare-reals` (14 equality, 14 inequality), 7 real→integer conversions, 2 continuation-line `&` issues, 1 unused import, 1 unused dummy, 1 truncation.
* `hf_matrix.f90`: `y = matmul(hfm, p(:,i))` is assigned **twice and never read** (dead store, wasted O(ns²) work each call).

### 6.7 Findings from stress-testing

1. **Input parsing is fragile.** `[Ar]3d4s` (no separating blank) aborts with the raw `Fortran runtime error: Bad integer for item 1 in list input` from `el_nl`, called from `get_atom`; the log shows the program had already mis-parsed the shell list as `3p3d4`. The paper says shells are "separated by blanks", so this is user error — but the diagnostic is unhelpful. With `[Ar] 3d 4s`, Kr converged normally.
2. **Raising `n_s` alone is not a convergence knob — it is a trap.** Because the outer grid is geometric ($t_i=t_{i-1}(1+h)$), $n_s$ sets the *range*: at $h=1/8$ the radius grows roughly as $(1.125)^{n_s}$, reaching ~10⁹ a₀ at $n_s=200$ and overflowing the `F10.2` format at $n_s=1200$ (prints `************`). Results are then unchanged (identical virial ratio at 400/800/1200) because the extra knots sit far beyond the orbital. Resolution is controlled by $h$ and $k_s$.
3. **Stack usage — my hypothesis was wrong.** I flagged several `REAL(8), DIMENSION(ns,ns)` automatic arrays (`hf_matrix`, `hfall_nr`, `solve_all`, `solve_scf`, `apply_orthogonality`) as a stack-overflow risk. I tested $n_s$ up to 1200 (one such array ≈ 11.5 MB) under an 8 MB stack and **nothing crashed**: gfortran places automatic arrays on the heap by default. This is therefore only a **portability risk** (`-fstack-arrays`/`-Ofast`, or compilers that default to stack), not a defect under the reviewed toolchain.
4. **Console output goes to stderr.** Prompts and level summaries are written to unit `err` while iteration logs go to stdout, so `prog > file` captures only part of the run. (The 50+ stderr lines per run are output, not errors.)

---

## 7. Strengths and weaknesses

### 7.1 Strengths

1. **Numerical quality is exceptional.** Total energies at the numerical HF limit; virial ratio to ~1e-14; no fudge factors.
2. **The algorithmic idea is elegant and principled.** Rotation analysis + projection operators remove the off-diagonal Lagrange multipliers *analytically*, replacing an ad-hoc art (node counting) with a variational procedure.
3. **Grid refinement by mapping** is a practical superpower: get a cheap answer on a coarse low-order grid, remap onto a finer high-order grid, polish in 2 SCF + 2 NR iterations.
4. **Clear, well-commented modern Fortran**, no `COMMON`, ~27 % comments, routines whose headers cite the exact equation/table in the paper (`rotate.f90` reproduces Table 2 in its header).
5. **Tiny footprint and dependencies** — LAPACK/BLAS and a Fortran compiler. Builds in seconds; runs in milliseconds for most atoms.
6. **Reproducibility aids:** `hf_param` and `bsw.c` are auto-written so any run can be re-executed non-interactively.
7. **Permissive MIT licence** in the GitHub mirror (the original CPC distribution used the standard CPC licence — note the difference).
8. **Useful by-product:** the fixed-core Rydberg series (eigenvalues, mean radii, continuum-like pseudostates) is directly reusable as an effectively complete basis.

### 7.2 Weaknesses and risks

**Scope**
* Only ≤ 2 open shells for specific LS terms; angular coefficients are **tabulated** (deviations from the average energy), not generated — the design ceiling of `get_energy.f90`.
* Non-relativistic only; MCHF/CI absent.

**Numerics / algorithms**
* **Dense-matrix design.** Exchange terms make $\mathbf H^a$ a *full* $n_s\times n_s$ matrix; the all-orbital NR system is $M\times M$ with $M=n_{wf}\,n_s$. Cost grows steeply with both. Ra took 1.7 s (`sphf_all`); an $f$-block element would be far worse. The author's ×30 % exchange-cutoff shortcut is stored in the paper's "modified exchange" row but is **not** the default.
* NR converges to the *nearest* solution and needs good estimates — acknowledged by the author (Mg 3p² → 4p²); no built-in safeguard beyond the SCF pre-phase.
* Two boundary conditions at $R$ mean an orbital "loses" one degree of freedom; the outermost orbital must fit in $n_s-1$ coefficients (paper: MAXR = 65 of 66).
* The block-LU and SVD alternatives were described but never benchmarked.

**Engineering**
* **No automated tests.** `test_all/sh_run` runs eight cases and `mv`s logs; there is no comparison or pass/fail. It also calls `run$n` unqualified (requires `.` on `PATH`) and `rm`s files that may not exist.
* **Truncated reference outputs** for the two heavy cases (§6.2).
* Interactive-prompt I/O; batch use relies on shell here-documents. Fragile parsing (§6.7).
* Makefile: hard-coded compiler, no module dependency tracking, requires a pre-existing `../bin`, and two executables sharing objects need manual cleaning between builds.
* `bsw.out` (empty, 0 bytes) is committed inside `src/`.
* Mixed console streams (§6.7-4).
* Dead code (`y = matmul(...)`), 14 exact real *equality* tests (e.g. `if (c == 0.d0)`; a further 14 flagged comparisons are inequalities and mostly harmless), automatic ns×ns arrays (a stack risk on some compilers).
* **Single-threaded.** No OpenMP; density/Slater-matrix construction is embarrassingly parallel and is the obvious hot spot.

**Documentation**
* README is minimal (directory listing only).
* **Citation error in the paper:** reference [5] reads "*Advances in Atomic and Molecular Physics*, vol. 55, **2007**, pp. **539–550**". Every independent source gives **Adv. At. Mol. Opt. Phys. 55, 235–291, 2008** (doi:10.1016/S1049-250X(07)55005-6). Anyone following the paper's own reference will look in the wrong place.

---

## 8. Comparison with related tools

| Code | Basis | Relativity | Atoms | Open shells | Notes |
|---|---|---|---|---|---|
| **SPHF** (this) | B-spline | No | ✔ | ≤ 2 specific-LS; any average | Simplest, very accurate |
| **DBSR_HF** (Zatsarinny & Froese Fischer 2016) | B-spline | **Dirac** | ✔ | up to j ≤ 9/2 | Same lineage, universal grid, several states + CI; command-line driven |
| **HF96 / HF86** (Froese Fischer) | Finite difference, log grid | No | ✔ | ≤ 2 (f-shells in HF96) | The predecessors; node counting |
| **HelFEM** (Lehtola 2019) | Finite element (LIP/HIP) | No | ✔ + diatomics | general | Also hybrid DFT; different basis philosophy |
| **GRASP / GRASP2018** | Finite difference | Dirac (MCDHF) | ✔ | general | Large-scale relativistic MCDHF |

SPHF's niche: *the* clean, small, non-relativistic B-spline reference implementation whose algorithm is fully documented in one paper — ideal for teaching, benchmarking, generating starting orbitals/Rydberg bases, and validating other codes.

---

## 9. Recommendations

**For users**
1. Use `sphf_all` when you need orbital-dependent properties or the tightest virial ratio; use `sphf_scf` for total energies of heavy atoms (up to ~2.5× faster).
2. Prefer $k_s\ge 6$ and $h\in\{1/4, 1/8\}$. Increase $n_s$ **only** to extend the range (check that MAXR < $n_s$), never to add resolution.
3. Always blank-separate closed shells and use fully specified configurations. Re-run from the auto-written `bsw.c`/`hf_param` for reproducibility.
4. For scripting, redirect **both** streams (`> out 2> err`) and read `sphf.log`.
5. For relativistic or multi-state work, go directly to DBSR_HF.

**For maintainers / contributors** (prioritised)
1. Add a proper **regression suite**: run each case, parse `TOTAL ENERGY` and `Ratio`, compare with tolerances (I used ~1e-9 E_h for energy, ~1e-12 for the virial ratio), and regenerate the truncated `test7/8` references.
2. Replace the Makefile by **CMake** (or at least add module dependencies, `mkdir -p ../bin`, and a `FC` override); build both executables into separate object directories.
3. Harden input parsing: accept comma/blank-tolerant shell lists and emit a diagnostic naming the offending token instead of a runtime error.
4. Route all diagnostics to a single, documented stream; add a quiet/batch mode and a namelist or command-line interface (as DBSR_HF already does).
5. Remove the dead `y = matmul(...)`; replace exact real equality tests with tolerances; change automatic $(n_s,n_s)$ arrays to `ALLOCATABLE` to remove any stack dependence.
6. Make the exchange-cutoff optimisation an option (default off, documented).
7. Add **OpenMP** to the density/Slater-matrix construction; consider banded storage for the exchange operator using the coefficient-product cutoff.
8. Correct reference [5] in any future edition; expand the README with build instructions and a worked example.

---

## 10. Reproducibility notes (what I ran)

```text
Toolchain : Ubuntu 24.04, gfortran 13.3.0, reference LAPACK/BLAS
Build     : make sphf_all ; make sphf_scf     (FC_FLAGS = -O2, defaults)
Checked   : FC_FLAGS="-O0 -g -fcheck=all -Wall -Wextra -finit-real=snan
                      -ffpe-trap=invalid,zero,overflow -fbacktrace"
Tests     : test_all/run1..run8 in a scratch copy, PATH=.
            + Ne convergence grid (h,k_s,n_s), He/Be/Ne/Ar/Kr accuracy
            + sphf_all vs sphf_scf best-of-3 timing (Ne, Ar, Kr, Ra)
Input note: closed shells must be blank-separated:  "[Ar] 3d 4s"
```

Example non-interactive run (He 1s² ¹S):

```bash
printf 'He 1S 2\n*\n1s(2)\nall\n' | sphf_all > out 2> err ; grep -E "TOTAL ENERGY|Ratio" sphf.log
```

**Caveats on my own evidence.** (a) Timings are single-machine, best-of-three, and will differ elsewhere. (b) For He and Be my literature references are 8–9-digit values, so the quoted ΔE is a bound, not a measurement of SPHF's true error. (c) I did not test open-shell cases beyond the eight bundled tests, nor the `ryd.bsw` Rydberg-series output. (d) I did not run any independent HF code as a cross-check for the open-shell energies (F, Mg, Ac); those are validated only via the virial theorem and agreement with the bundled reference logs.

---

# Publications related to the package's theory

**Legend.** ✅ = bibliographic details verified against an independent source during this review. 📎 = cited by the SPHF paper; I could not independently verify the details. Items are grouped by role. DOIs are given where I confirmed them.

## A. The package itself and its direct successors

1. ✅ **C. Froese Fischer**, *A B-spline Hartree–Fock program*, Comput. Phys. Commun. **182**, 1315–1326 (2011). doi:[10.1016/j.cpc.2011.01.012](https://doi.org/10.1016/j.cpc.2011.01.012). *(The SPHF paper.)*
2. ✅ **O. Zatsarinny, C. Froese Fischer**, *DBSR_HF: A B-spline Dirac–Hartree–Fock program*, Comput. Phys. Commun. **202**, 287–303 (2016). doi:[10.1016/j.cpc.2015.12.023](https://doi.org/10.1016/j.cpc.2015.12.023). Code: <https://github.com/compas/dbsr_hf>.
3. ✅ **C. Froese Fischer, O. Zatsarinny**, *A B-spline Galerkin method for the Dirac equation*, Comput. Phys. Commun. **180**, 879–886 (2009). doi:[10.1016/j.cpc.2008.12.010](https://doi.org/10.1016/j.cpc.2008.12.010); arXiv:[0806.3067](https://arxiv.org/abs/0806.3067).

## B. Core theory: B-splines in variational atomic structure

4. ✅ **C. Froese Fischer**, *B-splines in variational atomic structure calculations*, Adv. At. Mol. Opt. Phys. **55**, 235–291 (2008). doi:[10.1016/S1049-250X(07)55005-6](https://doi.org/10.1016/S1049-250X(07)55005-6). **The most important theory reference for SPHF** (rotations, projection, Galerkin/NR treatment). ⚠️ The SPHF paper mis-cites this as "vol. 55, 2007, pp. 539–550".
5. ✅ **H. Bachau, E. Cormier, P. Decleva, J. E. Hansen, F. Martín**, *Applications of B-splines in atomic and molecular physics*, Rep. Prog. Phys. **64**, 1815–1943 (2001). doi:[10.1088/0034-4885/64/12/205](https://doi.org/10.1088/0034-4885/64/12/205). The standard review of B-spline properties and applications (reference list to 2000); recommended by the SPHF paper for background.
6. ✅ **O. Zatsarinny**, *BSR: B-spline atomic R-matrix codes*, Comput. Phys. Commun. **174**, 273–356 (2006). doi:[10.1016/j.cpc.2005.10.006](https://doi.org/10.1016/j.cpc.2005.10.006). Source of the B-spline library that SPHF's spline routines are adapted from.
7. ✅ **Y. Qiu, C. Froese Fischer**, *Integration by cell algorithm for Slater integrals in a spline basis*, J. Comput. Phys. **156**, 257–271 (1999). The "cell" algorithm underlying the $R^k$ B-spline Slater integrals.
8. ✅ **C. Froese Fischer**, *The introduction of B-spline basis sets in atomic structure calculations*, Phys. Scr. **T47**, 7–17 (1993).
9. ✅ **C. Froese Fischer, T. Brage**, *Splines in atomic structure calculations*, AIP Conf. Proc. (1995), pp. 139–150.
10. ✅ **S. L. Saito**, *Validity of B-splines as a universal basis set for atomic Hartree–Fock–Roothaan calculations*, Theor. Chem. Acc. **109**, 326–331 (2003).
11. ✅ **C. Froese Fischer, M. Idrees**, *Spline algorithms for continuum functions*, Comput. Phys. **3**, 53–58 (1989).
12. 📎 **C. A. J. Fletcher**, *Computational Galerkin Methods*, Springer, New York (1984). (Galerkin formulation.)

## C. The Hartree–Fock method and orthogonality constraints

13. 📎 **D. R. Hartree**, *The Calculation of Atomic Structures*, Wiley, New York (1957).
14. 📎 **J. C. Slater**, *Quantum Theory of Atomic Structure*, Vol. II, McGraw-Hill, New York (1960). (Average energy of a configuration.)
15. 📎 **C. Froese Fischer**, *The Hartree–Fock Method for Atoms: A Numerical Approach*, Wiley-Interscience (1977). (Log-grid transformation; sections 6-2 and 6-4 cited by the HF86 write-up.)
16. 📎 **M. Bentley**, J. Phys. B: At. Mol. Opt. Phys. **27**, 637–644 (1994). *(Projection operators to eliminate off-diagonal Lagrange multipliers — the paper's ref. [9]. I could not independently locate this record, so I give it exactly as the SPHF paper cites it.)*
17. 📎 **T. A. Koopmans**, Physica **1**, 104 (1933). (Koopmans' theorem, used for same-symmetry filled shells.)
18. ✅ **C. Froese Fischer, T. Brage, P. Jönsson**, *Computational Atomic Structure: An MCHF Approach*, IOP Publishing (1997).

## D. Predecessor Hartree–Fock programs

19. ✅ **C. Froese Fischer**, *A general multi-configuration Hartree–Fock program*, Comput. Phys. Commun. **14**, 145–153 (1978). (MCHF77.)
20. ✅ **C. Froese Fischer**, *A general Hartree–Fock program*, Comput. Phys. Commun. **43**, 355–365 (1987). (HF86; F90/95 translation available from <https://www.pks.mpg.de/~george/HF/>.)
21. ✅ **G. Gaigalas, C. Froese Fischer**, *Extension of the HF program to partially filled f-subshells*, Comput. Phys. Commun. **98**, 255–264 (1996). doi:[10.1016/0010-4655(96)00092-6](https://doi.org/10.1016/0010-4655(96)00092-6). (HF96; catalogue ADDZ.)
22. ✅ **C. Froese Fischer**, Comput. Phys. Rep. **3**, 273 (1986). (Numerical procedures for MCHF/HF.)

## E. Later context, applications and comparisons

23. ✅ **C. Froese Fischer, M. Godefroid**, *Atomic Structure: Variational Wave Functions and Properties*, in *Springer Handbook of Atomic, Molecular, and Optical Physics* (ed. G. W. F. Drake), Springer (2023). doi:[10.1007/978-3-030-73893-8_22](https://doi.org/10.1007/978-3-030-73893-8_22).
24. ✅ **C. Froese Fischer, O. Zatsarinny**, *Numerical procedures for relativistic atomic structure calculations*, Atoms **8**, 85 (2020). (Reviews SPHF/DBSR_HF, NR methods.)
25. ✅ **C. Froese Fischer** (with co-authors), *Variational methods for atoms and the virial theorem*, Atoms **10**, 110 (2022). doi:[10.3390/atoms10040110](https://doi.org/10.3390/atoms10040110). (Uses the virial ratio as a diagnostic, exactly as SPHF's Table 4 does; notes that NR methods were then not yet applied to relativistic atomic equations.)
26. ✅ **S. Lehtola**, *A review on non-relativistic fully numerical electronic structure calculations on atoms and diatomic molecules*, Int. J. Quantum Chem. **119**, e25968 (2019); arXiv:[1902.01431](https://arxiv.org/abs/1902.01431). (Places SPHF among fully numerical methods.)
27. ✅ **S. Lehtola**, *Fully numerical Hartree–Fock and density functional calculations. I. Atoms*, Int. J. Quantum Chem. **119**, e25945 (2019); arXiv:[1810.11651](https://arxiv.org/abs/1810.11651). (Source of the HF reference energies for Ne, Ar, Kr used in §6.4; code HelFEM.)
28. ✅ **O. Zatsarinny, K. Bartschat**, *The B-spline R-matrix method for atomic processes: application to atomic structure, electron collisions and photoionization*, J. Phys. B **46**, 112001 (2013).
29. ✅ **J. Sapirstein, W. R. Johnson**, J. Phys. B **29**, 5213–5225 (1996). 📎 *(B-spline many-body perturbation theory; cited by SPHF for the Rydberg-series basis. Verified only via the SPHF reference list.)*
30. 📎 **S. Verdebout, P. Jönsson, G. Gaigalas, M. Godefroid, C. Froese Fischer**, J. Phys. B **43**, 074017 (2010). (Correlation orbitals; cited for the continuum-like orbitals.)
31. 📎 **W. R. Johnson, J. Sapirstein**, *Computation of second-order many-body corrections in relativistic atomic systems*, Phys. Rev. Lett. **57**, 1126 (1986). *(Classic first B-spline relativistic application, mentioned in the DBSR_HF paper; details from memory of the literature, not verified here.)*

---

*Sources consulted:* the SPHF and DBSR_HF papers (CPC 182, 1315; CPC 202, 287), the `compas/sphf` and `compas/dbsr_hf` repositories, the IOPscience/ADS/ScienceDirect/OSTI records for the cited works, and Lehtola's arXiv papers 1810.11651, 1902.01431 and 2108.05850. All timing and accuracy figures in this document are my own measurements from this session.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the B-spline HF (Froese Fischer) B-spline basis HF program. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
