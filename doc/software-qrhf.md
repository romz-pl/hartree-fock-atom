# Review of `qrhf`: a finite-difference Hartree–Fock program for atoms

## 0. Scope and evidence base (read this first)

This review is built from what I could verify from public sources. I want to be explicit about the limits, because `qrhf` is a small, sparsely documented research code and the usual sources of detail are missing.

**What I verified directly**

- The program's distribution page lists two source tarballs: `qrhf-1.0.tgz` (2014-09-23, ~960 KB) and `qrhf-2.0.tgz` (2018-03-13, ~1.9 MB).
- The author's homepage describes `qrhf` as a modified version of the Froese Fischer non-relativistic Hartree–Fock program for atoms, extended with several quasi-relativistic potentials, released under the GPL.
- The x2dhf README states that the HF orbitals shipped in `./hf_orbitals` were generated with "a finite-difference HF program for atoms (qrhf)".
- A 2022 *Molecular Physics* paper (Kobus and Kędziorski) uses "reference data of 12-digit accuracy provided by the GRASP2 and QRHF atomic codes" to validate a 2D finite-difference solver for the second-order Dirac equation.
- The underlying Froese Fischer HF program (Fortran 77 in CPC 1987, later a Fortran 90/95 version) uses `LOG(Z*R)` as the independent variable and `P/SQRT(R)` as the dependent variable, per its own documentation.

**What I could not verify**

- The contents of the tarballs. The distribution host is not reachable from my sandbox, and the tarball is binary so it cannot be read through page fetching. I did **not** read `qrhf`'s source, input format, build system, or README.
- Whether a dedicated, peer-reviewed paper describing `qrhf` exists. I found none. I searched for one specifically and only found papers that *use* it as a reference tool.

**Consequence.** Sections marked **[Verified]** rest on sources above. Sections marked **[Inferred]** are reasoned from the program's stated lineage and from standard theory; treat them as informed expectations to check against the source, not as facts about this code. Sections marked **[Not determined]** are things I could not establish. Where the answer depends on the source, I say what to look for.

---

## 1. Identity and provenance

### 1.1 What it is [Verified]

`qrhf` is a numerical (grid-based) Hartree–Fock solver for **atoms** in the central-field approximation. It is a derivative of Charlotte Froese Fischer's radial HF code, with the addition of quasi-relativistic potentials so that first-order relativistic effects can be folded into the self-consistent orbitals rather than added afterward.

### 1.2 Naming [Verified / Inferred]

The name pairs with the author's "quasi-relativistic Hartree–Fock" theme. The acronym expansion is my reading of the naming; the homepage itself only says "modified non-relativistic Hartree–Fock program … with various quasi-relativistic potentials included."

### 1.3 Relationship to x2dhf [Verified]

`qrhf` is **not** part of `x2dhf`. It is a separate, one-dimensional (radial) atomic code. Its role in the x2dhf ecosystem is as an independent generator of reference atomic orbitals used to seed SCF calculations. A related but distinct fact: x2dhf itself is a *two-dimensional* FD HF code that can also handle atoms, so for atoms the two programs overlap in purpose while differing in dimensionality and in the relativistic treatment.

| Aspect | `qrhf` | `x2dhf` |
|---|---|---|
| Systems | Atoms | Atoms and diatomic molecules |
| Grid dimensionality | 1D radial | 2D (prolate spheroidal) |
| Ancestry | Froese Fischer HF program | Laaksonen–Sundholm–Pyykkö diatomic code |
| Relativity | Quasi-relativistic potentials | Non-relativistic HF/DFT (plus one-particle relativistic work elsewhere) |
| Language | Fortran (per its Fortran lineage) [Inferred] | Fortran 90/95 + C |
| Documentation | Sparse; no dedicated paper found | CPC papers, arXiv review, user's guide |

### 1.4 Versions [Verified]

Two releases are hosted: 1.0 (2014) and 2.0 (2018). The size doubling (960 KB to 1.9 MB) suggests substantial additions between them, but I could not diff them, so **what changed is [Not determined]**.

---

## 2. Theoretical foundation

### 2.1 The central-field HF problem [Verified for the base program]

For a closed-shell or spherically averaged open-shell atom, each orbital is written as a radial function times a spherical harmonic:

$$
\psi_{nlm}(\mathbf r) = \frac{P_{nl}(r)}{r}\, Y_{lm}(\theta,\varphi).
$$

The radial HF equation has the form used in the Froese Fischer program's own documentation:

$$
P'' + \left(\frac{2Z}{r} - Y - \frac{l(l+1)}{r^{2}} - \epsilon\right)P = X + T,
$$

where $Y$ is the direct (Coulomb) potential, $X$ collects the exchange terms, and $T$ carries Lagrange-multiplier off-diagonal coupling between orbitals of the same symmetry that must remain orthogonal.

### 2.2 Coordinate transform and why it matters [Verified]

The Froese Fischer program integrates on a logarithmic grid, using $\ln(Zr)$ as the independent variable and $P/\sqrt r$ as the dependent variable. The purpose is standard: a logarithmic mesh resolves the sharp inner-shell structure near the nucleus while still reaching the long tail needed for diffuse orbitals, with a modest number of points. Because `qrhf` is described as a modification of that program, **it very likely inherits this transform [Inferred]**. Confirm in the source by looking for the mesh definition and the transformed equation.

### 2.3 Quasi-relativistic potentials [Verified as a feature; details Not determined]

The distinguishing feature is the inclusion of "various quasi-relativistic potentials." I could not determine *which* potentials are implemented. The standard candidates for a Breit–Pauli-style scalar treatment are:

$$
H_{\text{mv}} = -\frac{\alpha^{2}}{8}\, p^{4},
\qquad
H_{\text{D}} = \frac{\pi \alpha^{2}}{2}\, Z\, \delta(\mathbf r),
$$

the mass–velocity and one-body Darwin terms, in atomic units with $\alpha$ the fine-structure constant, plus a spin–orbit term

$$
H_{\text{so}} = \frac{\alpha^{2}}{2}\, \frac{Z}{r^{3}}\, \mathbf l \cdot \mathbf s.
$$

Whether `qrhf` uses exactly these, a Cowan-style scalar-relativistic Hartree–Fock form, or something else is **[Not determined]**. This matters, and it is the first thing to establish from the source: the mass–velocity operator is singular and is not bounded below, so it can only be used perturbatively or with regularisation, and different quasi-relativistic prescriptions give measurably different energies for heavy atoms.

### 2.4 Finite nuclear models [Verified for the author's wider work; Not determined for qrhf]

The 2022 paper that used `qrhf` for validation employs a Gaussian nuclear charge distribution, and x2dhf supports finite-nucleus models. Whether `qrhf` offers point, Gaussian, uniform-sphere, or Fermi nuclei is **[Not determined]**. For highly charged, high-$Z$ ions such as the $\mathrm{Kr}^{35+}$ and $\mathrm{Th}^{89+}$ test cases in that paper, nuclear size is a first-order effect, so a finite-nucleus option is what a reference code used that way would need. This is a strong hint but not confirmation.

### 2.5 Finite differences versus alternatives [Verified for general context]

The Lehtola HelFEM papers contrast finite-difference and finite-element approaches. In the diatomic case, HelFEM is strictly variational, whereas x2dhf energies are typically antivariational, approaching the converged value from below because of inaccuracies in the potential. That statement is made specifically about x2dhf versus HelFEM; **I would not assume it transfers unchanged to `qrhf`**, whose 1D radial equations are far better conditioned. Still, if you use `qrhf` energies as upper-bound benchmarks, check convergence behaviour yourself.

---

## 3. Numerical method

### 3.1 What is known [Inferred from lineage]

The Froese Fischer approach solves the radial equations by finite-difference integration on the transformed mesh, iterating to self-consistency. Its numerical procedures are shared with the MCHF family and are documented in *Computer Physics Reports* and in the book *The Hartree–Fock Method for Atoms* (Wiley, 1977).

### 3.2 What to check in the source [Not determined]

- Order of the difference scheme (the Numerov-type schemes in the lineage are fourth order in the mesh spacing for the second-derivative equation).
- Mesh parameters: step size $h$, number of points, and outer cutoff $r_{\max}$.
- Convergence criteria and whether they are on orbital residuals, energy, or both.
- How Lagrange multipliers and orthogonality are handled.
- Whether extended precision is available. The "12-digit accuracy" claim in the 2022 paper implies at least careful double precision, but I cannot say what the code does internally.

---

## 4. Capabilities

**Established** (from the sources above):

- Non-relativistic HF for atoms (inherited base capability).
- Quasi-relativistic corrections via a choice of potentials.
- Production of reference orbitals used to seed x2dhf calculations.
- Production of reference energies to about 12 significant digits for hydrogenic and few-electron highly charged test systems, per the 2022 validation.

**Plausible but unconfirmed:**

- Open-shell / spherically averaged configurations, as in the base program.
- Finite nuclear size.
- Ions and excited configurations.

**Almost certainly absent [Inferred]:**

- Density-functional methods (a HF-lineage code with no stated DFT support).
- Molecular calculations.
- Post-HF correlation methods.
- Multiconfiguration treatment. The base program is single-configuration HF; MCHF is a separate program in the family.
- Fully relativistic Dirac–Fock. The whole point of "quasi-relativistic" is to approximate this without solving the four-component equations.

---

## 5. Accuracy and validation

### 5.1 What the evidence supports [Verified]

The strongest evidence is indirect but meaningful. A peer-reviewed 2022 paper treats `qrhf`, alongside GRASP2, as a source of 12-digit reference data for several lowest $\sigma$, $\pi$, $\delta$ and $\varphi$ states of $\mathrm{H}$, $\mathrm{Kr}^{35+}$ and $\mathrm{Th}^{89+}$. That is a hydrogenic (one-electron) validation.

### 5.2 An important caveat about that validation

For a one-electron ion there is no electron–electron interaction, so the test exercises the radial grid, the relativistic potentials, and the nuclear model, but **not** the direct and exchange machinery, which is the hard part of a Hartree–Fock code. It is good evidence that the kinematic and nuclear parts are correct, and weaker evidence about many-electron behaviour. For many-electron accuracy, you would compare against tabulated numerical HF limits (for example the Froese Fischer or Clementi–Roetti-type tabulations) or against the HelFEM atomic results, which the Lehtola group reports for atoms. I did not find a published head-to-head comparison involving `qrhf`.

### 5.3 A verification recipe

If you adopt the program, a proportionate acceptance test is:

1. Hydrogenic ions at several $Z$: compare energies to the analytic Dirac and Schrödinger values, separating the relativistic correction.
2. Helium and beryllium: compare total HF energies to established numerical HF limits.
3. A closed-shell heavier atom such as neon or argon: compare to HelFEM or literature values.
4. Grid convergence study: halve the step and extend $r_{\max}$, and confirm the energy stabilises to the digits you need.

---

## 6. Software engineering assessment

I could not inspect the code, so this section is limited to what the distribution page and lineage allow.

| Dimension | Assessment | Basis |
|---|---|---|
| License | GPL | Author's homepage [Verified] |
| Distribution | Plain tarballs, no version control host found | Distribution page [Verified] |
| Releases | Two, 2014 and 2018; no evidence of later activity | Directory listing [Verified] |
| Maintenance | Appears dormant, though the author remains active on x2dhf | Inferred from dates |
| Documentation | No dedicated paper or manual located | Search [Verified absent] |
| Tests | Unknown | [Not determined] |
| Build system | Unknown | [Not determined] |
| Language | Fortran lineage, version not confirmed | Inferred |

**Practical implications for someone with a C++ and HPC background.** The code is a single-threaded radial solver; it is not a parallelism target, and its cost per atom is small. Its value is as a *reference implementation*, not a production engine. If you need to embed atomic HF in a larger pipeline, a modern library is easier to integrate; `qrhf` is best treated as an oracle to validate against.

---

## 7. Strengths and weaknesses

### 7.1 Strengths

- **Traceable lineage.** It descends from a heavily used, well-documented atomic HF program with a textbook behind it.
- **Independent reference.** It has already served as a benchmark source in a peer-reviewed study, and as the generator of orbitals distributed with x2dhf.
- **Relativity without a full Dirac solver.** Quasi-relativistic potentials give a cheap route to first-order relativistic effects inside the SCF loop.
- **GPL and freely downloadable.**

### 7.2 Weaknesses and risks

- **No dedicated documentation or paper.** This is the main risk: you cannot cite a methods paper for it, and behaviour must be learned from the source.
- **Validation is narrow in the public record.** The visible test is hydrogenic; many-electron accuracy is unpublished as far as I could find.
- **Ambiguous relativistic model.** Without knowing which potentials are used, results for heavy atoms cannot be interpreted against other relativistic treatments.
- **Dormant.** No sign of updates since 2018.
- **Discoverability.** It is hosted on a personal academic page rather than a repository, which complicates reproducibility and citation.

---

## 8. Comparison with alternatives

| Program | Method | Relativity | Scope | Notes |
|---|---|---|---|---|
| `qrhf` | 1D finite difference, HF | Quasi-relativistic potentials | Atoms | Reference-oriented; sparsely documented |
| Froese Fischer HF (CPC 1987; F90/95) | 1D finite difference, HF | Non-relativistic | Atoms | The parent program |
| MCHF / GRASP family | Multiconfiguration HF / Dirac–Fock | Breit–Pauli / fully relativistic | Atoms | GRASP2 was used as a companion reference alongside `qrhf` |
| HelFEM | Finite element, variational | Non-relativistic | Atoms and diatomics | Hundreds of LDA/GGA/meta-GGA functionals; well documented and papered |
| x2dhf | 2D finite difference | Non-relativistic HF/DFT | Atoms and diatomics | Documented in CPC and an arXiv review |

**Recommendation on choosing.** For non-relativistic atomic HF benchmarks with citable provenance, HelFEM or the Froese Fischer program are easier to defend in a paper. Use `qrhf` when you specifically want the quasi-relativistic potentials it provides, or when you need to reproduce the reference values cited in the 2022 study.

---

## 9. Open questions to resolve from the source

Downloading `qrhf-2.0.tgz` and reading it would settle these directly:

1. Which quasi-relativistic potentials are implemented, and are they applied self-consistently or perturbatively?
2. Which nuclear charge models are available?
3. What are the mesh, difference scheme, and precision?
4. What changed between 1.0 and 2.0?
5. Is there an input-format description, a test suite, or sample outputs?
6. Does it handle open shells and spherically averaged configurations?

---

# Publications related to the package's theory

The list is grouped by role. Entries are given as they appear in the sources I consulted. Items I could confirm from a citation in a source are marked ✔. Items I recall or infer but did **not** confirm against a source in this session are marked ⚠ and should be verified before you cite them.

## A. The parent atomic Hartree–Fock program

- ✔ C. Froese Fischer, *The Hartree–Fock Method for Atoms*, Wiley-Interscience, 1977. (Cited in the program documentation for the coordinate transform and equation forms, Sec. 6-2 and 6-4.)
- ✔ C. Froese Fischer, "A general Hartree–Fock program", *Comput. Phys. Commun.* **43**, 355–365 (1987). (The Fortran 77 program that the Fortran 90/95 version is derived from.)
- ✔ C. Froese Fischer, T. Brage, P. Jönsson, *Computational Atomic Structure: An MCHF Approach*, IOP Publishing, 1997. (Cited for usage examples and shared numerical procedures.)
- ✔ C. Froese Fischer, "A B-spline Hartree–Fock program", *Comput. Phys. Commun.* (2011). (Listed among the author's atomic HF numerical work; page range not confirmed.)
- ✔ C. Froese Fischer, W. Guo, Z. Shen, "Spline methods for multiconfiguration Hartree–Fock calculations", *Int. J. Quantum Chem.* **42**, 849–867 (1992).
- ✔ C. Froese Fischer, A. Senchuk, "Numerical procedures for relativistic atomic structure calculations", *Atoms* **8**, 85 (2020). (Cited as pp. 1–20 in the source; treat the volume/article number as needing verification.)
- ✔ C. Froese Fischer, M. Godefroid, T. Brage, P. Jönsson, G. Gaigalas, "Advanced multiconfiguration methods for complex atoms: I. Energies and wave functions", *J. Phys. B* **49**, 182004 (2016).
- ⚠ C. Froese Fischer, 1957 numerical atomic work, *Mon. Not. R. Astron. Soc.* (cited in the x2dhf review as the origin of fully numerical atomic calculations "already in the 1950s"; I did not confirm the full citation).
- ✔ D. R. Hartree, *The Calculation of Atomic Structures*, Wiley, 1957.

## B. Relativistic atomic structure references used alongside qrhf

- ✔ I. P. Grant, B. J. McKenzie, P. H. Norrington, D. F. Mayers, N. C. Pyper, *Comput. Phys. Commun.* **21**, 207–231 (1980). (The GRASP2 lineage, used as a companion reference in the 2022 validation.)
- ✔ I. P. Grant, *Relativistic Quantum Theory of Atoms and Molecules*, Springer, 2007.
- ✔ F. A. Parpia, C. Froese Fischer, I. P. Grant, "GRASP92", *Comput. Phys. Commun.* **94**, 249 (1996).
- ✔ C. Froese Fischer, G. Gaigalas, P. Jönsson, J. Bieroń, "GRASP2018", *Comput. Phys. Commun.* **237**, 184–187 (2019).
- ✔ O. Zatsarinny, C. Froese Fischer, "DBSR: a B-spline Dirac–Hartree–Fock program", *Comput. Phys. Commun.* **202**, 287–303 (2016).
- ✔ R. D. Cowan, *The Theory of Atomic Structure and Spectra*, Univ. of California Press, 1981. (The pseudo-relativistic HFR code, which includes Blume–Watson spin–orbit, mass-variation and one-body Darwin terms, is the standard quasi-relativistic comparator; the book citation is standard but was not itself retrieved in this session.)

## C. The author's use of qrhf and closely related work

- ✔ J. Kobus, A. Kędziorski, "Two-dimensional, finite-difference method of solving the Dirac equation for diatomic molecules revisited", *Mol. Phys.* **120**, e2092563 (2022). (Uses `qrhf` and GRASP2 as 12-digit reference atomic codes.)
- ✔ J. Kobus, "A finite difference Hartree–Fock program for atoms and diatomic molecules", *Comput. Phys. Commun.* **184**, 799–811 (2013). doi:10.1016/j.cpc.2012.09.033.
- ✔ J. Kobus, S. Lehtola, "Review of the finite difference Hartree–Fock method for atoms and diatomic molecules, and its implementation in the x2dhf program", *Comput. Phys. Commun.* **311**, 109576 (2025). doi:10.1016/j.cpc.2025.109576; arXiv:2408.03679.
- ✔ J. Kobus, L. Laaksonen, D. Sundholm, "A numerical Hartree–Fock program for diatomic molecules", *Comput. Phys. Commun.* **98**, 346–358 (1996). doi:10.1016/0010-4655(96)00098-7.
- ✔ J. Kobus, "Hartree–Fock limit values of multipole moments, polarizabilities and hyperpolarizabilities for atoms and diatomic molecules", *Phys. Rev. A* **91**, 022501 (2015). doi:10.1103/PhysRevA.91.022501.
- ✔ J. Kobus, "Finite-difference versus finite-element methods", *Chem. Phys. Lett.* **202**, 7–12 (1993). doi:10.1016/0009-2614(93)85342-L.
- ✔ J. Kobus, "Vectorizable algorithm for the (multicolour) successive overrelaxation method", *Comput. Phys. Commun.* **78**, 247–255 (1994). doi:10.1016/0010-4655(94)90003-5.
- ✔ J. Kobus, "Overview of finite difference Hartree–Fock method algorithm, implementation and application", *AIP Conf. Proc.* **1504**, 189–208 (2012).
- ✔ J. Kobus, "Numerical Hartree–Fock methods for diatomic molecules", in *Handbook of Molecular Physics and Quantum Chemistry*, ed. S. Wilson, Wiley, 2002.
- ✔ L. Laaksonen, P. Pyykkö, D. Sundholm, "Fully numerical Hartree–Fock methods for molecules", *Comput. Phys. Rep.* **4**, 313–344 (1986). doi:10.1016/0167-7977(86)90021-3.

## D. Finite-difference versus finite-element and modern numerical atomic HF

- ✔ S. Lehtola, "Fully numerical Hartree–Fock and density functional calculations. I. Atoms", *Int. J. Quantum Chem.* e25945 (2019). doi:10.1002/qua.25945; arXiv:1810.11651.
- ✔ S. Lehtola, "Fully numerical Hartree–Fock and density functional calculations. II. Diatomic molecules", *Int. J. Quantum Chem.* e25944 (2019). doi:10.1002/qua.25944; arXiv:1810.11653.
- ✔ S. Lehtola, "Fully numerical calculations on atoms with fractional occupations and range-separated exchange functionals", *Phys. Rev. A* **101**, 012516 (2020). doi:10.1103/PhysRevA.101.012516; arXiv:1908.02528.
- ✔ S. Lehtola, "A review on non-relativistic fully numerical electronic structure calculations on atoms and diatomic molecules", *Int. J. Quantum Chem.* (2019); arXiv:1902.01431.
- ✔ S. Lehtola, "Assessment of initial guesses for self-consistent field calculations. Superposition of atomic potentials: simple yet efficient", *J. Chem. Theory Comput.* **15**, 1593–1604 (2019). doi:10.1021/acs.jctc.8b01089; arXiv:1810.11659. (Relevant because `qrhf` orbitals seed SCF, and x2dhf now supports superposition-of-atomic-potentials starts.)
- ✔ S. Lehtola, C. Steigemann, M. J. T. Oliveira, M. A. L. Marques, "Recent developments in LIBXC", *SoftwareX* **7**, 1–5 (2018). doi:10.1016/j.softx.2017.11.002.
- ✔ H. Åström, S. Lehtola, "Systematic study of hard-wall confinement induced effects on atomic electronic structure", *J. Phys. Chem. A* **129**, 2791 (2025). doi:10.1021/acs.jpca.4c05641; arXiv:2408.11595.

## E. Suggested background not confirmed in this session

These are standard references for the quasi-relativistic terms discussed in Section 2.3. I include them because they are the natural sources for the physics, but I did not retrieve or verify them here, so treat them as pointers to look up.

- ⚠ H. A. Bethe, E. E. Salpeter, *Quantum Mechanics of One- and Two-Electron Atoms*, Springer, 1957. (Source of the Breit–Pauli reduction and its mass–velocity, Darwin, and spin–orbit terms.)
- ⚠ R. E. Moss, *Advanced Molecular Quantum Mechanics*, Chapman & Hall, 1973.
- ⚠ Literature on Pauli-type operators and the unboundedness of the mass–velocity term, including regularisation approaches (for example the zeroth-order regular approximation, ZORA, and Douglas–Kroll–Hess theory) as the modern alternatives to a Breit–Pauli-style quasi-relativistic treatment.

---

*Method note.* Verification of the primary source was limited by network access to the distribution host. To make this review authoritative, download `qrhf-2.0.tgz`, read its README and main source, and resolve the questions in Section 9; the Sections 2.3, 2.4, 3.2 and 4 claims marked *[Not determined]* would then become firm.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the qrhf Finite-difference Hartree-Fock program for atoms. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
