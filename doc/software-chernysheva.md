# The Chernysheva–Cherepkov–Radojević Hartree–Fock Program for Atoms: A Review

> **Scope and evidential basis.** This review covers the numerical self-consistent-field (SCF) Hartree–Fock (HF) program for atoms published by L. V. Chernysheva, N. A. Cherepkov and V. Radojević in *Computer Physics Communications* **11**, 57 (1976), its frozen-core companion (*Comput. Phys. Commun.* **18**, 87 (1979)), and the wider ATOM program system in which both were absorbed.
>
> **Read this first — what I could and could not verify.** The two original CPC papers are paywalled. I could read their metadata, keyword lists, abstracts and reference lists, but **not their full text**. The algorithmic statements below therefore fall into three tiers, and I mark them throughout:
>
> - **[Verified]** — stated in a source I actually read (keywords, abstracts, the 2022 *Atoms* review, the 2020/2021 Waide–Green–Gribakin B-spline papers, etc.).
> - **[Standard theory]** — textbook HF material that any radial-HF program of this kind must implement; I am supplying it from general knowledge, not from the CCR papers.
> - **[Inference]** — my reasoned reconstruction. Treat with appropriate caution, and check against the original listing before relying on it.
>
> Where I could not determine something (exact source-line counts, precise convergence criteria, the open-shell coupling scheme), I say so rather than guess.

---

## 1. Executive summary

The CCR program is a **fully numerical, radial, single-configuration Hartree–Fock solver for atoms**, written in Fortran, that grew out of the theory group at the Ioffe Institute (Leningrad) led by M. Ya. Amusia. It is not a stand-alone curiosity: it is the **HF engine of the ATOM system**, the family of codes that first made many-body calculations (RPAE) of atomic photoionization routine. According to a recent editorial, one of the key computational developments of ATOM was an efficient numerical solution of the HF equations, "making it a standard starting point for higher-level calculations."

Its lasting significance is less about raw HF accuracy (modern B-spline and finite-difference codes now reach much higher precision) and more about **three things**:

1. **A clean separation of ground-state SCF from "frozen-core" excited-state orbitals** — the design that lets a single code supply both occupied orbitals and a complete set of discrete and continuum excited orbitals for downstream correlation calculations.
2. **Robust, published, portable numerics** (Numerov integration, tridiagonal elimination, energy-eigenvalue adjustment) from an era of severe hardware limits.
3. **A long reuse tail.** The frozen-core HF code was later adapted by the convergent close-coupling (CCC) group to treat targets beyond helium, and the wider ATOM `hfgr` code is used as a benchmark by modern B-spline HF programs.

**Bottom line for a practitioner:** use it as a *historically important, still-instructive, numerically transparent reference implementation* and as a source of consistent HF + frozen-core orbital sets; do not choose it today for high-precision or relativistic work, or for open-shell atoms needing anything beyond spherical single-configuration treatment.

---

## 2. Provenance and authorship

### 2.1 The papers

| Item | Citation | Notes |
|---|---|---|
| Main program | L. V. Chernysheva, N. A. Cherepkov, V. Radojević, "Self-consistent field Hartree–Fock program for atoms," *Comput. Phys. Commun.* **11**(1), 57–73 (1976). DOI: 10.1016/0010-4655(76)90040-0 | 17 pages. Author affiliation for the corresponding author is listed as Novi Sad (Yugoslavia). |
| Frozen-core companion | L. V. Chernysheva, N. A. Cherepkov, V. Radojević, "Frozen core Hartree–Fock program for atomic discrete and continuous states," *Comput. Phys. Commun.* **18**(1), 87 (1979). DOI: 10.1016/0010-4655(79)90026-2 | Computes excited-state HF radial functions of a single electron in the "frozen core" field. |

**Naming note.** The 1976 paper is what the question calls "the program." Its exact CPC title is "Self-consistent field Hartree–Fock program for atoms." The 1979 paper is a distinct but tightly coupled program that *consumes* the output of the first. In practice, when people cite "the Chernysheva–Cherepkov–Radojević code" they often mean this pair together.

### 2.2 Institutional and historical context [Verified]

- The work sits inside the Ioffe Institute theory group under Miron Ya. Amusia, where, in 1966–67 and the years following, the RPAE was justified and implemented, and "the first computer codes of what later became the ATOM suite" were created.
- The Serbian collaboration is real and visible in the record: the 1979 paper carries the note "Supported in part by Serbian Research Funds, Yugoslavia," and one author's 1978/79 academic-year address was at LURE, Orsay. The wider Ioffe–Serbia link is described in the Amusia tribute editorial as one of the productive international ties the group built.
- The programs were first distributed as **Ioffe Institute preprints in Russian** before or alongside the CPC papers (see §11), which is why English-language citation trails are thin for the early 1970s history.

### 2.3 Lineage: from preprint to ATOM to ATOM-M [Verified]

The Springer chapter on the ATOM-M system lists the preprint chain. The HF-specific items in it are:

- Chernysheva & Cherepkov, *Numerical calculation of Hartree–Fock wave functions*, Ioffe preprint No. 337 (1971).
- Chernysheva, Cherepkov & Radojević, *Numerical solution of the Hartree–Fock self-consistent field equations for atoms*, Ioffe preprint No. 486 (1975) — the direct precursor of the 1976 CPC paper.
- Chernysheva, Cherepkov & Radojević, *Numerical calculation of Hartree–Fock wave functions of discrete and continuous states in the frozen core field*, Ioffe preprint No. 487 (1975) — the precursor of the 1979 CPC paper.
- Chernysheva, Cherepkov & Sheftel, *Numerical calculations of the coefficients in the Hartree–Fock equations*, Ioffe preprint No. 760 (1982) — the angular-coefficient machinery.
- Chernysheva, Kuchiev & Yakhontov, *Numerical calculations of ground-state atomic wave functions in Hartree–Fock–Dirac approximation*, Ioffe preprint No. 1015 (1986) — the relativistic extension.
- Later ES-computer ports (1987, 1991) and the frozen-core-for-LS-terms extension (Chernysheva & Gribakin, 1987).

The system was consolidated in the Amusia–Chernysheva handbook (1997) and modernised as **ATOM-M** (Amusia, Chernysheva, Semenov, 2016, Russian; Springer English edition 2021 as *Computation of Atomic and Molecular Processes*).

---

## 3. What the program computes

### 3.1 Physical problem [Verified + standard theory]

The 1976 program solves the **HF equations for a many-electron atom in the independent-particle, single-configuration approximation**. Its keyword list (as indexed by CPC) is a compact description of the whole method:

> atomic structure; energy level; wave functions; independent-particle approximation; single-particle model; Hartree–Fock; (many-)electron correlations; electron configuration; single-configuration approximation; atomic shell; *nl*-(sub)shell; self-consistent field; iteration convergency acceleration; Numerov method; chasing method; Gauss elimination method without pivoting; energy eigenvalue adjustment.

The 1979 companion "calculates the excited state HF radial wave function of single-electron in the 'frozen core' (FC) field of other electrons."

### 3.2 Outputs [Verified via the 2022 ATOM review]

Within ATOM, the HF module supplies:

- Ground-state wave functions in the HF (and, in the extended system, Hartree–Fock–Dirac) approximation.
- Excited-state functions "consistent with the functions of the ground state."
- Excited-state functions in the **continuous spectrum** at prescribed energies in a fixed field, with or without orthogonalisation to the ground-state functions.
- Excited-state functions in the **discrete spectrum** for given principal quantum numbers in a frozen-core field, again with optional orthogonalisation.
- Wave functions for exotic projectiles (μ-meson, positron) in the atomic field.

These sets are the raw input to dipole and Coulomb matrix elements and, downstream, RPAE.

---

## 4. Theory implemented

### 4.1 The HF equations for a spherical atom [Standard theory]

Writing each orbital as $\phi_{nlm\sigma}(\mathbf r)=r^{-1}P_{nl}(r)\,Y_{lm}(\hat{\mathbf r})\,\chi_\sigma$, angular and spin variables are integrated out analytically. This is confirmed as the general ATOM approach: angular integration and spin summation "are carried out analytically," leaving radial functions as solutions of the HF equations by successive approximations.

For a closed-shell atom the radial equations take the form (atomic units)

$$
-\frac{1}{2}\frac{d^2 P_{i}}{dr^2}+\left[-\frac{Z}{r}+\frac{l_i(l_i+1)}{2r^2}+V_{\mathrm{dir}}(r)\right]P_{i}(r)+\int_0^\infty U(r,r')\,P_{i}(r')\,dr' = E_i\,P_{i}(r),
$$

where $V_{\mathrm{dir}}$ is the local direct (Hartree) potential and $U(r,r')$ is the non-local exchange kernel. The integro-differential character of this equation (the exchange term) is what distinguishes HF from the simpler Hartree method and drives the numerical difficulty.

In the exchange-included form used by all such codes, the direct and exchange contributions are expanded in Slater radial integrals $R^k$ and $F^k, G^k$ with angular coefficients — precisely the "coefficients in the Hartree–Fock equations" that Chernysheva, Cherepkov and Sheftel programmed as a separate module (Ioffe preprint 760, 1982).

### 4.2 Single-configuration approximation and open shells

The keywords name the **single-configuration approximation**. For open-shell atoms, codes of this type either (a) spherically average the open shell, or (b) treat specific LS terms. A recent B-spline paper notes that the radial equation "may also be used for open-shell electronic configurations under the further approximation of spherical averaging," and that for neutral open-shell ground states the resulting difference "is relatively small since the variation … arises due to only one incomplete subshell" — this is a statement about that code, not about the CCR program.

**[Inference / unverified for CCR]:** the 1976 program appears to be a *configuration-average* (or single-term) treatment; extension to specific LS terms of the frozen-core electron appears as a later, separate development (Chernysheva & Gribakin, 1987, "for states with LS term"). I could not confirm from the 1976 full text exactly which open-shell coupling scheme it uses. **Check the original before depending on term-specific energies.**

### 4.3 Frozen-core approximation [Verified concept, standard theory for the physics]

In the frozen-core method the occupied ground-state orbitals $\{P_{nl}\}$ are held fixed after the SCF step. The additional electron (or a positron or meson) is then found by solving the same radial equation with $V_{\mathrm{dir}}$ and $U$ *fixed* by the ground-state orbitals. This yields:

- discrete excited orbitals (an electron added to the $N$-electron core), and
- continuum orbitals at chosen energies.

For positrons the exchange term is simply dropped and the sign of the interaction flipped. This is the physical basis for treating excited and scattering states in ATOM, and it underlies RPAE, where "single-particle transition amplitudes are determined in terms of wave functions."

---

## 5. Numerical methods [Verified from keywords; details are standard-theory reconstruction]

The CPC keyword list names four numerical ingredients. What each does in a radial HF solver is well established; I flag where I am reconstructing.

| Keyword | Role in a radial HF solver | Status |
|---|---|---|
| **Numerov method** | Fourth-order-accurate finite-difference integration of the second-order radial ODE of the form $P'' = f(r)\,P + g(r)$ on a mesh. | Named in keywords; role is standard. |
| **Chasing method (tridiagonal "progonka")** | Solution of the tridiagonal linear system that results from discretising the radial equation *including* the inhomogeneous exchange term; equivalent to the Thomas algorithm. | Named in keywords; role is **[Inference]**. |
| **Gauss elimination without pivoting** | Direct solution of linear systems arising in the discretised problem; without pivoting is safe when the matrix is diagonally dominant, as tridiagonal difference systems typically are. | Named in keywords; use is **[Inference]**. |
| **Energy eigenvalue adjustment** | Iterative correction of the orbital energy $E_i$ so the solution has the right number of nodes and satisfies boundary conditions (matching/shooting-type refinement). | Named in keywords; mechanism **[Inference]**. |
| **Iteration convergence acceleration** | Extrapolation of successive SCF iterates to speed or stabilise convergence. | Named in keywords. See §5.2. |

### 5.1 Grid and change of variable [Verified in spirit]

The 2022 ATOM review states that "the change of variable required in the calculation of wave functions is carried out in the calculation of all characteristics of atoms," and that for continuum quantities it is advisable to place most sampling points at low energy by using the electron momentum $p=E^{1/2}$ as the integration variable. This shows the system uses a **non-uniform, transformed radial mesh** rather than a uniform grid. I could not determine the exact transformation used in the 1976 code (a logarithmic-type mapping is typical for such solvers, but that is **[Inference]**).

### 5.2 Convergence acceleration — a well-attested detail [Verified via a citing paper]

A B-spline HF paper (Waide, Green & Gribakin) describes a mixing scheme to speed convergence:

$$
P^{(m)}_{i,\mathrm{est}}(r)=(1-\alpha)\,P^{(m)}_i(r)+\alpha\,P^{(m-1)}_i(r),
\qquad
\alpha=\frac{E_i^{(m)}-E_i^{(m-1)}}{E_i^{(m)}-2E_i^{(m-1)}+E_i^{(m-2)}},
$$

and states that this is "a scheme utilised by Amusia and Chernysheva in the hfgr code." This is an **Aitken-type extrapolation on the orbital energy**, used to weight the mixing of successive orbital iterates. Because `hfgr` is the descendant of the CCR program inside the handbook, this is strong indirect evidence that the CCR line uses an Aitken/$\Delta^2$-type accelerator — consistent with the CPC keyword "iteration convergency acceleration." **Caveat:** the citing paper attributes it to `hfgr` in the 1997 handbook, not to the 1976 CPC listing specifically. I am inferring the same scheme was already present in 1976 from the keyword and the lineage; that is plausible but **not directly verified.**

The same paper reports that this mixing was needed for Zn even after other stabilisation, illustrating the practical fragility of naive SCF iteration — an issue the older code addressed with this acceleration, and that the modern "annealing" of the electron–electron interaction addresses differently.

### 5.3 Boundary conditions and orthogonality

**[Standard theory]:** orbitals must vanish at the origin ($P(0)=0$) and decay at large $r$; SCF solutions for different $nl$ with the same $l$ must be orthogonal, either automatically (when the exchange operator is Hermitian and eigenvalues are distinct) or via Lagrange multipliers. The ATOM description confirms orthogonalisation to ground-state functions is an explicit option for excited states.

---

## 6. Program structure and the ATOM software model [Verified from the 2022 review]

The CCR programs are best understood as the seed of the ATOM architecture. The 2022 review describes ATOM as organised into **four module types**:

1. **Executive modules** — procedures without formal parameters that hold variable declarations, input/print of initial data, the algorithm and output; "each of which solves an independent physical problem."
2. **Specialised modules** — subroutines or functions with formal parameters.
3. **Service modules** — input, printing of initial data or results.
4. **Generic modules** — standard mathematical methods.

At the time of that review, the application-program layer contained **more than 50 executive modules, more than 10 service modules, more than 70 specialised modules and more than 16 generic modules**, plus a database of wave functions and input/output physical characteristics for each atom and process. Fortran is used *without* compiler-specific extensions "to facilitate the implementation of the ATOM system on other computers." Detailed printing of intermediate quantities is deliberately built in as a diagnostic — described as playing the role of "diagnostics in a natural experiment."

For the *HF-specific* code, two design consequences matter:

- The HF solver is deliberately **decoupled** from downstream physics: ground-state and excited-state orbitals are written to a database and reused by many executive programs (photoionization, GOS, scattering phases, Auger widths, etc.).
- The code was engineered for **portability across the Soviet computing environment** (BESM-6, ES computers), which the ATOM-M reference list makes explicit: BESM-ALGOL, BESM-6 operating-system manuals, and ES-computer ports appear among the citations.

**Language.** The 1976 program is Fortran. The 1997 handbook states the ATOM software "is written in FORTRAN and may be used on VAX or UNIX-based machines."

---

## 7. Capabilities and limitations

### 7.1 What it does well

- **Transparent, minimal HF core.** A small set of well-understood numerical building blocks; an excellent teaching and reference implementation.
- **Consistent orbital sets.** Ground-state SCF plus frozen-core discrete and continuum orbitals from the same framework, ready for RPAE and other many-body work.
- **Extensibility.** The same numerical skeleton was extended to Hartree–Fock–Dirac (1986), LS-term frozen-core states (1987), and non-diagonal energy parameters (1989 preprint No. 1319).
- **Exotic projectiles.** Positron and meson orbitals in the same field (dropping exchange for positrons).
- **Proven downstream utility.** Decades of photoionization, scattering, and Auger results in ATOM rest on it.

### 7.2 Limitations

| Limitation | Detail | Evidence |
|---|---|---|
| **Non-relativistic in the base program** | Relativistic (Dirac–Fock) treatment requires the later HFD preprint/extension. | Verified (separate 1986 preprint). |
| **Single-configuration** | No configuration mixing; electron correlation must be added *afterwards* (RPAE, MBPT). The CPC keywords list "single-configuration approximation" and "(many-)electron correlations" together. | Verified (keywords). |
| **Open-shell treatment is approximate** | Configuration/spherical averaging; term-dependent effects need later extensions. | Inference (see §4.2). |
| **Finite-difference on a radial mesh** | Accuracy is limited by mesh and Numerov order; modern B-spline HF benchmarks reach $\sim10^{-6}$ a.u. orbital-energy agreement with reference HF tables. | Verified for the *B-spline* code; a direct CCR-vs-reference precision figure is not something I found. |
| **Convergence fragility for hard cases** | Negative ions and some transition metals (e.g. Zn) can oscillate without acceleration or damping. | Verified for a modern code; the acceleration scheme is attributed to `hfgr`. |
| **Heavy atoms** | Non-relativistic HF is strained for heavy atoms; a modern paper notes the "limits of the nonrelativistic HF approximation" for the heaviest species (Rn). | Verified (for that context). |
| **Legacy Fortran** | Fixed-form, COMMON-block-era code; integration into modern workflows requires porting effort. | General knowledge / Inference. |

### 7.3 A note on the *extra nodes* phenomenon [Verified, physically instructive]

Modern work highlights that HF inner orbitals for heavier atoms can acquire **more nodes than $n-l-1$** because of the non-local exchange interaction. For Kr the 1s orbital shows two extra nodes, versus none for Ne. Any solver, including the CCR program's node-counting and eigenvalue-adjustment logic, has to cope with this. A code that assumes strictly $n-l-1$ nodes to identify a state can misidentify orbitals in such cases. I could not verify whether the CCR code makes that assumption; this is **a specific thing to test** if you use it for heavier atoms.

---

## 8. Validation, benchmarking, and downstream use

### 8.1 Benchmark role [Verified]

The `hfgr` code (the HF program of the Amusia–Chernysheva handbook, the direct descendant of the CCR line) was used by Waide, Green and Gribakin as a comparator for static dipole polarisabilities of noble-gas atoms. Their B-spline values were in each case "closer to the tabulated reference values," while `hfgr` results were "in good agreement" with them.

Their Table 5 (values in atomic units) illustrates the comparison:

| Atom | B-spline | `hfgr` | Reference |
|---|---|---|---|
| He | 0.997236 | 0.997167 | 1.3837675 |
| Ne | 1.974636 | 1.973492 | 2.6717 |
| Ar | 10.140367 | 10.13132 | 11.0747 |
| Kr | 15.861810 | 15.84444 | 16.7656 |

**Interpretation.** Both HF-level results sit well below the reference values, as expected: HF omits correlation, and static polarisabilities are correlation-sensitive. The point of the table is the *near-agreement between two independent HF implementations*, which supports the numerical soundness of `hfgr` — not agreement with experiment. I am reading the table as the source presents it; the reference column is the literature/handbook value quoted there.

### 8.2 Use in convergent close-coupling [Verified]

The CCC group adapted the HF codes for more complex single-valence-electron targets, where "the direct and exchange interaction with the inner electron core needs to be taken into account." The abstract of that work states the HF codes developed in the group of Miron Amusia "have been adapted." The Amusia tribute editorial adds that "the self-consistent field and frozen-core HF computer codes from the ATOM system have been adopted," with electron scattering on Li and double photoionisation of H⁻ and Li⁻ as examples. This is direct evidence that the **SCF + frozen-core pair** — precisely the two CCR papers — remain in active use, four-plus decades on.

### 8.3 Citation footprint [Verified, limited]

ScienceDirect lists **82 citing articles** for the 1976 paper (as of the page I retrieved), spanning collision theory (Bray & Stelbovics, *Adv. At. Mol. Opt. Phys.* 1995; Bray & Fursa, PRA 1996; Bray, PRA 1994), photoionisation time delays (Kheifets & Ivanov, PRL 2010; Guénot et al., PRA 2012), and Rydberg lifetimes (Theodosiou, PRA 1984). This shows the code was used as a *tool* by groups well outside the original one. I did not enumerate all 82.

---

## 9. Comparison with related HF codes

The table contrasts approach and role. Entries about other codes are from the sources indicated; entries about the CCR code follow the tiers above.

| Code / family | Method | Scope | Notes |
|---|---|---|---|
| **CCR (1976) / ATOM `hfgr`** | Finite-difference radial mesh; Numerov + tridiagonal solve; Aitken-type acceleration | Atoms; SCF + frozen-core excited/continuum | Basis of ATOM/RPAE. |
| **Froese Fischer MCHF family** | Numerical multiconfiguration HF; later spline variants | Atoms; MCHF beyond single configuration | Refs: Froese 1963; Froese Fischer 1970, 1972, 1978; MCHF atomic-structure package. CCR is contemporaneous and cited alongside these in numerical-HF reviews. |
| **BSHF (Waide–Green–Gribakin)** | B-spline basis → generalised eigenproblem; interaction "annealing" and mixing | Atoms and arbitrary central potentials; produces HF + frozen-core excited states | Modern successor in spirit; uses `hfgr` as a comparator. |
| **x2dhf (Kobus, Lehtola)** | 2-D finite-difference HF | Atoms and diatomic molecules | Different dimensionality; cited as a comparator in the fully-numerical-methods literature. |
| **Roothaan–Hartree–Fock (basis-set) atoms** | Analytic (e.g. Slater-type) basis | Atoms/molecules | Basis-set-limited; the fully numerical approach removes that limit. |

The 2019 arXiv review of non-relativistic fully numerical atomic and diatomic calculations lists the CCR 1976 paper (reference 405) in its survey of numerical HF programs, alongside Froese's numerical HF work and Froese Fischer's MCHF programs. That placement is the clearest independent statement of where the community sees this program: **a canonical member of the first generation of fully numerical atomic HF codes.**

---

## 10. Practical guidance

### 10.1 When to use it

- You want a **small, readable reference implementation** of radial HF to understand or validate your own solver.
- You need **consistent ground-state + frozen-core discrete/continuum orbitals** for many-body or scattering calculations, following the ATOM/RPAE or CCC pattern.
- You are reproducing or extending **historical ATOM results**.

### 10.2 When not to

- You need **sub-microhartree precision** or systematically improvable accuracy — prefer a spline/finite-element approach.
- You need **relativistic** results — use the HFD extension or a dedicated Dirac–Fock code.
- You need **term-resolved open-shell energies** or **multiconfiguration** wave functions — use MCHF-type tools.
- You need to run on modern parallel/HPC infrastructure without porting. (Given your background in C++ modernisation and HPC, note that the practical work here would be porting fixed-form Fortran to a modern layout; the numerical kernels themselves — a Numerov sweep and a tridiagonal solve — are small and straightforward to reimplement.)

### 10.3 Checks to run before trusting results

1. Compare total energies and orbital energies for closed-shell noble gases against tabulated HF references (for instance the Saito 2009 tables used in the B-spline benchmark).
2. Verify **node counts** for inner orbitals of heavier atoms (see §7.3).
3. Confirm **virial-theorem** satisfaction to gauge SCF convergence and mesh adequacy. *[Standard practice; not from the CCR papers.]*
4. Test **negative ions and Zn-type cases** for oscillation.
5. Check sensitivity of excited-state orbitals to the **outer mesh limit** and to **orthogonalisation** choices.

---

## 11. Related publications on the theory

The list is organised by theme. Entries marked **†** are ones I confirmed in the sources I read (via reference lists or bibliographic pages); I cite them as they appear there. I have **not** independently opened every listed item, and I have not padded the list with references I could not trace.

### 11.1 The CCR papers and their direct preprint chain

1. † L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Self-consistent field Hartree–Fock program for atoms*, Comput. Phys. Commun. **11**, 57–73 (1976). DOI: 10.1016/0010-4655(76)90040-0
2. † L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Frozen core Hartree–Fock program for atomic discrete and continuous states*, Comput. Phys. Commun. **18**, 87 (1979). DOI: 10.1016/0010-4655(79)90026-2
3. † L. V. Chernysheva, N. A. Cherepkov, *Numerical calculation of Hartree–Fock wave functions*, Ioffe Institute preprint No. 337, Leningrad (1971) (in Russian).
4. † L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Numerical solution of the Hartree–Fock self-consistent field equations for atoms*, Ioffe preprint No. 486 (1975) (in Russian).
5. † L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Numerical calculation of Hartree–Fock wave functions of discrete and continuous states in the frozen core field*, Ioffe preprint No. 487 (1975) (in Russian).
6. † L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Numerical calculations of photoionization cross sections and optical oscillator strengths of atoms in Hartree–Fock approximation and allowing for multielectron correlations in one transition*, Ioffe preprint No. 589 (1978) (in Russian).
7. † L. V. Chernysheva, N. A. Cherepkov, S. I. Sheftel, *Numerical calculations of the coefficients in the Hartree–Fock equations*, Ioffe preprint No. 760 (1982) (in Russian).
8. † L. V. Chernysheva, M. Ya. Amusia, V. F. Orlov, S. K. Semenov, N. A. Cherepkov, *Numerical calculations of electron wavefunctions in the frozen-core field with input of various off-diagonal energy parameters*, Ioffe preprint No. 1319 (1989) (in Russian).
9. † L. V. Chernysheva, S. K. Semenov, N. A. Cherepkov, *Numerical calculation of the Hartree–Fock wave functions of an electron in the field of a frozen core on an ES computer*, Ioffe preprint No. 1574 (1991) (in Russian).

### 11.2 Extensions: relativistic, LS-term, ports

10. † L. V. Chernysheva, M. Yu. Kuchiev, V. L. Yakhontov, *Numerical calculations of ground-state atomic wave functions in Hartree–Fock–Dirac approximation*, Ioffe preprint No. 1015 (1986) (in Russian).
11. † L. V. Chernysheva, G. F. Gribakin, *Numerical calculation of the wave functions of an electron in the field of a "frozen" core for states with LS term*, Ioffe Institute, Leningrad (1987) (in Russian).
12. † L. V. Chernysheva, B. V. Gul'tsev, M. M. Dovganich, *Numerical calculation of atomic wave functions in the Hartree–Fock approximation on an ES computer*, Ioffe Institute, Leningrad (1987) (in Russian).
13. † L. V. Chernysheva, *Complex of programs for automation of atomic calculations*, Ioffe Institute, Leningrad (1981) (in Russian).

### 11.3 The ATOM / ATOM-M system (books and reviews)

14. † M. Ya. Amusia, L. V. Chernysheva, *Computation of Atomic Processes: A Handbook for the ATOM Programs*, Institute of Physics Publishing, Bristol and Philadelphia (1997). DOI: 10.1201/9781003040859
15. † M. Ya. Amusia, L. V. Chernysheva, *Automatic System of Atomic Structure Investigations*, Nauka, Leningrad (1983) (in Russian).
16. † M. Ya. Amusia, L. V. Chernysheva, S. K. Semenov, *ATOM-M. Algorithms and Programs for the Study of Atomic and Molecular Processes*, Nauka, St. Petersburg (2016) (in Russian).
17. † M. Ya. Amusia, L. V. Chernysheva, *Computation of Atomic and Molecular Processes*, Springer Series on Atomic, Optical, and Plasma Physics, vol. 117, Springer, Cham (2021). ISBN 978-3-030-85142-2. Chapter 22, "Structure of the ATOM-M System" (pp. 357–405), DOI: 10.1007/978-3-030-85143-9_22
18. † L. V. Chernysheva, V. K. Ivanov, *ATOM Program System and Computational Experiment*, Atoms **10**(2), 52 (2022). DOI: 10.3390/atoms10020052
19. † M. Ya. Amusia, V. K. Ivanov, N. A. Cherepkov, L. V. Chernysheva, *Processes in Many-Electron Atoms*, Nauka, St. Petersburg (2006) (in Russian).
20. † M. Ya. Amusia, L. V. Chernysheva, V. G. Yarzhemsky, *Handbook of Theoretical Atomic Physics: Data for Photon Absorption, Electron Scattering, Vacancies Decay*, Springer, Berlin (2012).
21. † A. S. Kheifets, G. Gribakin, V. K. Ivanov, *"Atoms" Special Issue (Many-Electron and Multiphoton Atomic Processes: A Tribute to Miron Amusia)* (editorial), Atoms **11**(2), 18 (2023). DOI: 10.3390/atoms11020018

### 11.4 Foundational many-body theory the programs feed (RPAE and related)

22. † M. Ya. Amusia, N. A. Cherepkov, S. I. Sheftel, *Many-body correlations and the photoeffect*, Phys. Lett. A **24**, 541 (1967).
23. † M. Ya. Amusia, N. A. Cherepkov, L. V. Chernysheva, *Many-electron correlations in photoabsorption in the M-shell of Ar*, Phys. Lett. A **31**, 553 (1970). DOI: 10.1016/0375-9601(70)90348-8
24. † M. Ya. Amusia, N. A. Cherepkov, L. V. Chernysheva, *Cross section for the photoionization of noble-gas atoms with allowance for multielectron correlations*, Sov. Phys. JETP **33**, 90 (1971) [Zh. Eksp. Teor. Fiz. **60**, 160].
25. † M. Ya. Amusia, V. K. Ivanov, N. A. Cherepkov, L. V. Chernysheva, *Interference effects in photoionization of noble gas atoms outer s-subshells*, Phys. Lett. A **40**, 361 (1972). DOI: 10.1016/0375-9601(72)90531-2
26. † M. Ya. Amusia, N. A. Cherepkov, *Many-electron correlations in scattering processes*, Case Studies in Atomic Physics **5**, 47–179 (1975).
27. † M. Ya. Amusia, *Atomic Photoeffect*, Plenum Press, New York (1990).
28. † N. A. Cherepkov, L. V. Chernysheva, *Random phase approximation with exchange for open-shell atoms: Photoionization of Cl*, Phys. Lett. A **60**, 103 (1977) [listed as vol. 60, issue 2].
29. † G. Wendin, *Collective effects in atomic photoabsorption spectra. III. Collective resonance in the 4d¹⁰ shell in Xe*, J. Phys. B **4**, 1080 (1971).
30. † M. Ya. Amusia, E. G. Drukarev, V. G. Gorshkov, M. O. Kazachkov, *Two-electron photoionization of helium*, J. Phys. B **8**, 1248 (1975).
31. † M. Ya. Amusia, N. A. Cherepkov, L. V. Chernysheva, S. G. Shapiro, *Elastic scattering of slow positrons by helium*, J. Phys. B **9**, L531 (1976).
32. † J. Goldstone, *Proc. R. Soc. Lond. A* **239**, 267 (1957) — many-body perturbation theory (cited in the CCR paper's own reference list; volume and page as shown there).

### 11.5 HF theory background cited by the CCR paper and its milieu

33. † D. R. Hartree, *The Wave Mechanics of an Atom with a Non-Coulomb Central Field*, Proc. Camb. Phil. Soc. **24**, 89, 111, 426 (1928).
34. † V. Fock, Z. Phys. **61**, 126 (1930).
35. † J. C. Slater, Phys. Rev. **35**, 210 (1930).
36. † D. R. Hartree, W. Hartree, *Self-consistent field, with exchange, for beryllium*, Proc. R. Soc. A **150**, 9 (1935).
37. † D. R. Hartree, *The Calculation of Atomic Structures*, Wiley, New York (1957) (also cited in the CCR reference list).
38. † C. Froese Fischer, *The Hartree–Fock Method for Atoms: A Numerical Approach*, Wiley, New York (1977).
39. † C. Froese Fischer, *A multi-configuration Hartree–Fock program*, Comput. Phys. Commun. **1**, 151 (1970).
40. † C. Froese Fischer, *A multi-configuration Hartree–Fock program with improved stability*, Comput. Phys. Commun. **4**, 107 (1972).
41. † C. Froese Fischer, *Numerical solution of general Hartree–Fock equations for atoms*, J. Comput. Phys. **27**, 221 (1978).
42. † C. Froese, *Numerical solution of the Hartree–Fock equations*, Can. J. Phys. **41**, 1895 (1963).
43. † H. P. Kelly, in *Atomic Inner-Shell Processes*, ed. B. Crasemann, Academic Press, New York (1963).
44. † C. Roothaan, *Self-consistent field theory for open shells of electronic systems*, Rev. Mod. Phys. **32**, 179 (1960).
45. † S. Huzinaga, *Applicability of Roothaan's self-consistent field theory*, Phys. Rev. **120**, 866 (1960).
46. † E. R. Davidson, *Spin-restricted open-shell self-consistent-field theory*, Chem. Phys. Lett. **21**, 565 (1973).
47. † N. C. Handy, M. T. Marron, H. J. Silverstone, *Phys. Rev.* **180**, 45 (1969) — asymptotic behaviour of HF orbitals.
48. † V. A. Dzuba, V. V. Flambaum, P. G. Silvestrov, J. Phys. B **15**, L575 (1982) — extra nodes from exchange.
49. † M. G. Kozlov, V. V. Flambaum, Phys. Rev. A **87**, 042511 (2013) — physical implications of exchange-induced nodes.
50. † D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, Singapore (1988) (cited in the ATOM-M chapter).
51. † M. Abramowitz, I. A. Stegun (eds.), *Handbook of Mathematical Functions*, NBS AMS-55 (1964) (cited in the ATOM-M chapter).

### 11.6 Modern reuse, benchmarks and comparators

52. † I. Bray, X. Weber, D. V. Fursa, A. S. Kadyrov, B. I. Schneider, S. Pamidighantam, M. Cytowski, A. S. Kheifets, *Taking the Convergent Close-Coupling Method beyond Helium: The Utility of the Hartree–Fock Theory*, Atoms **10**(1), 22 (2022). DOI: 10.3390/atoms10010022
53. † D. T. Waide, D. G. Green, G. F. Gribakin, *BSHF: A program to solve the Hartree–Fock equations for arbitrary central potentials using a B-spline basis*, Comput. Phys. Commun. **250**, 107112 (2020).
54. † D. T. Waide, D. G. Green, G. F. Gribakin, *B-spline basis Hartree–Fock method for arbitrary central potentials: atoms, clusters and electron gas*, arXiv:2108.05850 (2021).
55. † S. L. Saito, *Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method*, At. Data Nucl. Data Tables **95**, 836 (2009) — HF reference data used as a benchmark.
56. † J. Kobus, S. Lehtola, *Review of the finite difference Hartree–Fock method for atoms and diatomic molecules, and its implementation in the x2dhf program* (arXiv:2408.03679).
57. † A review of non-relativistic fully numerical electronic structure calculations on atoms and diatomic molecules (arXiv:1902.01431), which lists the CCR 1976 paper in its survey of numerical HF programs.
58. † Application papers citing the 1976 program (examples surfaced on the CPC citation page): I. Bray, A. T. Stelbovics, *Adv. At. Mol. Opt. Phys.* **35** (1995); I. Bray, D. V. Fursa, PRA **54** (1996); I. Bray, PRA **49** (1994); A. S. Kheifets, I. A. Ivanov, PRL **105**, 233002 (2010); D. Guénot et al., PRA **85**, 053424 (2012); C. E. Theodosiou, PRA **30** (1984).

---

## 12. Gaps, uncertainties and how to close them

To keep this review honest, here is exactly what I could **not** establish and how you could resolve each item:

| Open question | Why it matters | How to resolve |
|---|---|---|
| Exact mesh/change-of-variable used in the 1976 code | Determines accuracy and the practical outer-radius limit | Read the 1976 paper's full text and program listing (available via institutional access to ScienceDirect or the CPC Program Library). |
| Precise open-shell coupling treatment | Determines validity for open-shell atoms | Same; also compare to the 1987 LS-term extension. |
| Whether the Aitken-type acceleration is in the 1976 listing itself or only in the later `hfgr` | Affects claims about original convergence behaviour | Inspect the 1976 listing directly. |
| Actual source size and platform in 1976 (CDC, IBM, BESM-6) | Portability history | CPC program summary / library entry. |
| Current availability of the code (CPC Program Library status, ATOM-M distribution) | Practical reuse | The 2021 Springer ATOM-M book states codes are "presented and made available in this book"; confirm the exact distribution channel with Springer/Chernysheva. |
| Quantified accuracy against modern HF references | Whether to use it for precision work | Run the noble-gas comparison in §10.3 against the Saito 2009 tables. |

I also did not access the full text of the 2022 supplementary list ("List of main publications in which the results were obtained using the ATOM program system"), which the *Atoms* article says is downloadable as supplementary material; that would be the best single source for an even more exhaustive applications bibliography.

---

## 13. Suggested further reading order

1. **Chernysheva & Ivanov, *Atoms* 10, 52 (2022)** — open access; the best single overview of the system and its physics.
2. **CCR 1976 and 1979 CPC papers** — the primary sources for the algorithm (institutional access needed).
3. **Amusia & Chernysheva, *Computation of Atomic Processes* (1997)** — the handbook version with the `hfgr` code and usage.
4. **Waide, Green & Gribakin (2020/2021)** — a modern, open-access counterpoint with explicit convergence treatment and a direct `hfgr` comparison.
5. **Froese Fischer, *The Hartree–Fock Method for Atoms* (1977)** — the standard numerical-HF companion text.
6. **Bray et al., *Atoms* 10, 22 (2022)** — a concrete recent example of the frozen-core HF machinery in use.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Chernysheva–Cherepkov–Radojević Hartree-Fock program for atoms. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
