# Hartree-Fock Software Packages Dedicated to Single Atoms

A survey of programs whose primary purpose is Hartree-Fock (HF) or closely related self-consistent-field (SCF) calculations on **isolated atoms** (and ions). General-purpose molecular quantum chemistry codes that merely happen to run atoms are listed separately at the end for context.

**Scope note.** "Exhaustive" is a goal, not a guarantee: many atomic codes are legacy Fortran programs distributed via journal program libraries, and some are no longer maintained. Entries below were verified against published program descriptions and repositories; where a code goes beyond plain HF (multiconfiguration, relativistic, DFT), that is stated.

---

## 1. Fully numerical (grid / finite-difference / finite-element) HF for atoms

| Package | Method | Description |
|---|---|---|
| **HelFEM** | Finite element, nonrelativistic HF/DFT | A suite of programs for finite element calculations on atoms and diatomic molecules at the Hartree-Fock or density-functional levels. The atomic solver uses a basis of the form $\chi_{nlm}(r,\theta,\phi)=r^{-1}B_n(r)Y_l^m(\hat{\mathbf r})$ with finite element shape functions, allowing arbitrary accuracy. Supports hybrid functionals and electric properties, and hundreds of LDA/GGA/meta-GGA functionals via Libxc. Repository: <https://github.com/susilehtola/HelFEM> |
| **x2dhf** | Finite-difference HF/DFT | Two-dimensional finite-difference program that finds virtually exact HF and DFT solutions for diatomic molecules and atoms. Solution quality depends on grid size and arithmetic precision. Also treats one-electron systems with several model potentials. Repository: <https://github.com/x2dhf/x2dhf> |
| **qrhf** | Finite-difference HF (atoms) | Finite-difference Hartree-Fock program for atoms, used alongside HelFEM to generate reference data for x2dhf tests. |
| **dftatom** | Radial Schrödinger/Dirac solver (Kohn-Sham) | Fortran 95 solver for atomic structure that computes both Schrödinger and Dirac radial wavefunctions (LDA/DFT rather than HF). Uses outward Poisson integration to reach $10^{-8}$ Hartree accuracy for heavy atoms up to uranium. Open source with Python and C wrappers. Included here because it is widely used as a numerical benchmark and radial-solver library for atomic codes. |
| **B-spline HF** (Froese Fischer) | B-spline basis HF | Fortran 95 program that replaces the usual system of differential equations with a system of matrix generalized eigenvalue problems, one per orbital. Part of the Froese Fischer atomic-structure ecosystem. |

## 2. Classic atomic SCF programs (Slater / Gaussian basis or numerical radial)

| Package | Method | Description |
|---|---|---|
| **atmscf (Columbus version)** | Basis-set-expansion HF | Revised and extended Columbus version of the 1963 Chicago atomic SCF program. Principal present use is developing Gaussian basis sets for molecular calculations. Treats ground states of essentially all atoms within LS coupling, excited states with large angular-momentum orbitals, and can include relativistic effects via effective core potentials. Fortran 90; CPC catalogue ID ADVR. |
| **Chernysheva–Cherepkov–Radojević HF program** | Numerical HF | Self-consistent-field Hartree-Fock program for atoms (*Comput. Phys. Commun.* **11**, 57, 1976), with a later frozen-core variant for atomic single-electron discrete and continuum states. A long-standing reference implementation for atomic photoionization work. |
| **Froese Fischer numerical HF programs** | Numerical HF | Foundational numerical solution of the HF equations for atoms (*Can. J. Phys.* **41**, 1895, 1963; *J. Comput. Phys.* **27**, 221, 1978; monograph *The Hartree-Fock Method for Atoms*, 1977). Ancestor of the MCHF/ATSP line below. |
| **Roothaan-Hartree-Fock atomic programs** | Analytic (Slater-type) basis HF | Family of programs implementing Roothaan's SCF theory for open-shell atoms, based on the Roothaan/Huzinaga/Davidson open-shell formalism; historically used to produce tabulated atomic HF wavefunctions. |

## 3. Multiconfiguration Hartree-Fock (MCHF), nonrelativistic

| Package | Method | Description |
|---|---|---|
| **MCHF / ATSP2K** | MCHF + Breit-Pauli | MCHF atomic-structure package using dynamic memory allocation, sparse matrix methods, and a modern angular library. Meant for large-scale calculations in an orthogonal orbital basis for groups of LS terms of arbitrary parity. Breit-Pauli operators (spin-orbit, spin-other-orbit, spin-spin, orbit-orbit) are available. Computes transition rates of all types, isotope shifts, hyperfine constants, and g-factors. Includes the single-configuration HF as a special case. |
| **Froese Fischer MCHF (original)** | MCHF | Original multiconfiguration Hartree-Fock program (*Comput. Phys. Commun.* **1**, 151, 1970; improved-stability version **4**, 107, 1972). |
| **LSGEN** | Configuration-state-list generator | Program to generate configuration-state lists of LS-coupled basis functions, used to set up MCHF calculations. |
| **Biegler-König / Hinze numerical MCSCF** | Numerical MCSCF | Nonrelativistic numerical MCSCF program for atoms, built on Froese Fischer's earlier developments with significant modifications and extensions. |

## 4. Relativistic (Dirac-Hartree-Fock / Dirac-Fock) atomic codes

| Package | Method | Description |
|---|---|---|
| **GRASP2018** | MCDHF + RCI | General Relativistic Atomic Structure Package. Fully relativistic four-component multiconfiguration Dirac-Hartree-Fock (MCDHF) and relativistic configuration interaction (RCI), suited to medium and heavy atomic systems. Fortran 95 (Froese Fischer, Gaigalas, Jönsson, Bieroń; *Comput. Phys. Commun.* **237**, 184, 2019). |
| **GRASP2K** | MCDHF | Collection of programs for large-scale relativistic calculations based on MCDHF; a modification and extension of GRASP92. Version 3 replaced the `njgraf` recoupling module with `librang`, extended coefficients of fractional parentage to make lanthanides and actinides feasible, and added `jj2lsj` (percentage composition of the wavefunction in LSJ coupling). |
| **GRASP92** | MCDHF | Predecessor of GRASP2K by Parpia, Froese Fischer and Grant. |
| **GRASP (original)** | MCDF | Original general-purpose relativistic atomic structure program (Grant et al. 1980; Dyall et al., *Comput. Phys. Commun.* **55**, 425, 1989). |
| **GRASPG** | MCDHF with CSF generators | Extension of GRASP2018 based on configuration-state-function generators (*Comput. Phys. Commun.* **312**, 109604, 2025). |

## 5. Semi-empirical Hartree-Fock atomic structure suites

| Package | Method | Description |
|---|---|---|
| **Cowan code (RCN, RCN2, RCG, RCE)** | HF with relativistic corrections (HFR) | Widely used suite for atomic spectra. **RCN** computes single-configuration radial wavefunctions for a spherically symmetrized atom via the HF method; **RCN2** computes radial (Slater) integrals including configuration interaction; **RCG** builds and diagonalizes the Hamiltonian to obtain levels, wavelengths, and E1/M1/E2 radiative rates (also autoionization rates and plane-wave Born excitation cross-sections); **RCE** least-squares fits Slater parameters to experimental levels. Relativistic effects enter as perturbations. |
| **CATS (Cowan Atomic Structure code, Los Alamos)** | HFR | Los Alamos automation of the Cowan codes (LA-11436-M), replacing several hand operations with automated steps. Part of the broader Theoretical Atomic Physics Code Development suite. |

## 6. Atomic solvers embedded in pseudopotential / basis-generation workflows

| Package | Method | Description |
|---|---|---|
| **ATOM (SIESTA project)** | All-electron atomic solver (DFT) | Program for atomic calculations and generation of pseudopotentials, used in the SIESTA workflow. DFT-based rather than HF. Listed because its documentation discusses the self-interaction cancellation that is exact only in true HF. |
| **Hartree-Fock pseudopotential generators** | Atomic HF/Dirac-Fock inversion | Atomic-level tools that invert the HF equations per angular-momentum channel to produce effective core potentials. Example: the TNDF pseudopotentials, built from Dirac-Fock all-electron states. |

## 7. Reusable libraries used by atomic SCF programs

| Package | Description |
|---|---|
| **OpenOrbitalOptimizer** | Reusable open-source C++ library for iterative solution of coupled SCF equations, motivated by the observation that orbital optimizers in ERKALE, HelFEM, Psi4, PySCF, and OpenMolcas are all disparate despite using similar algorithms. Applicable to atomic SCF workflows. |

---

## Related but not atom-specific

These general-purpose codes can run single-atom HF calculations but are not dedicated to atoms:

| Package | Note |
|---|---|
| **ERKALE** | HF/DFT program (Lehtola) with a fully numerical atomic component that shares heritage with HelFEM. |
| **PySCF** | Python SCF package. Its `atom` initial guess performs spherically averaged atomic spin-restricted HF calculations on the fly to build a minimal basis, i.e. it embeds an atomic HF solver. |
| **aromanro/HartreeFock** | C++ HF/post-HF (MP2, CCSD(T), CIS, TDHF/RPA) with Gaussian orbitals. Single atoms are run by placing one atom into a dummy molecule. |
| **Psi4, OpenMolcas** | Molecular packages listed as using orbital optimizers of the kind OpenOrbitalOptimizer aims to replace. |

---

## Caveats

- **Version and maintenance status** were not independently verified for every legacy code; several exist mainly as journal-distributed Fortran.
- **DFT-only codes** (`dftatom`, `ATOM`) are included because they are standard companions of atomic HF work, but they do not perform Hartree-Fock themselves.
- Classic Roothaan/Slater-basis atomic programs are described at the level of the literature rather than a single canonical repository, since no single maintained distribution was identified.

## Key references

- S. Lehtola, *Fully numerical Hartree-Fock and density functional calculations. I. Atoms*, Int. J. Quantum Chem. **119**, e25945 (2019), arXiv:1810.11651.
- S. Lehtola, *A review on non-relativistic fully numerical electronic structure calculations on atoms and diatomic molecules*, Int. J. Quantum Chem. **119**, e25968 (2019), arXiv:1902.01431.
- L. V. Chernysheva, N. A. Cherepkov, V. Radojević, *Self-consistent field Hartree-Fock program for atoms*, Comput. Phys. Commun. **11**, 57 (1976).
- C. Froese Fischer, G. Gaigalas, P. Jönsson, J. Bieroń, *GRASP2018*, Comput. Phys. Commun. **237**, 184 (2019).
- R. D. Cowan code overview: *Cowan Code: 50 Years of Growing Impact on Atomic Physics*, Atoms (2019).
- *dftatom*, arXiv:1209.1752.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive list of Hartree-Fock software packages dedicated to single atoms. For each package provide the short description of the package. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
