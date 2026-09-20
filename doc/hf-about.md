# The Hartree-Fock Method for a Single Atom

## Table of Contents

1. [Introduction and Scope](#1-introduction-and-scope)
2. [The Atomic Hamiltonian](#2-the-atomic-hamiltonian)
3. [The Variational Principle](#3-the-variational-principle)
4. [The Slater Determinant Ansatz](#4-the-slater-determinant-ansatz)
5. [Energy Expectation Value](#5-energy-expectation-value)
6. [Derivation of the Hartree-Fock Equations](#6-derivation-of-the-hartree-fock-equations)
7. [Central-Field Approximation and Separation of Variables](#7-central-field-approximation-and-separation-of-variables)
8. [Angular Algebra: Reduction to Radial Equations](#8-angular-algebra-reduction-to-radial-equations)
9. [The Radial Hartree-Fock Equations](#9-the-radial-hartree-fock-equations)
10. [Closed-Shell Atoms](#10-closed-shell-atoms)
11. [Open-Shell Atoms: Restricted, Unrestricted, and Multiconfigurational Variants](#11-open-shell-atoms)
12. [The Self-Consistent Field (SCF) Procedure](#12-the-self-consistent-field-scf-procedure)
13. [Numerical Solution of the Radial Equations](#13-numerical-solution-of-the-radial-equations)
14. [Basis-Set Expansion: Roothaan-Hartree-Fock](#14-basis-set-expansion-roothaan-hartree-fock)
15. [Orbital Energies and Koopmans' Theorem](#15-orbital-energies-and-koopmans-theorem)
16. [Worked Examples: He, Be, Ne](#16-worked-examples)
17. [Correlation Energy and Limitations](#17-correlation-energy-and-limitations)
18. [Relativistic Extension: Dirac-Hartree-Fock](#18-relativistic-extension-dirac-hartree-fock)
19. [Post-Hartree-Fock Methods](#19-post-hartree-fock-methods)
20. [Practical Implementation Notes](#20-practical-implementation-notes)
21. [Summary](#21-summary)
22. [References](#22-references)

---

## 1. Introduction and Scope

The **Hartree-Fock (HF) method** is the foundational mean-field approximation to the many-electron Schrödinger equation. For a single atom with $N$ electrons and nuclear charge $Z$, HF replaces the intractable $3N$-dimensional problem with a set of coupled one-electron equations, in which each electron moves in the average field of the nucleus and the remaining $N-1$ electrons.

Key ideas:

- The many-electron wavefunction is approximated by a **single Slater determinant** (or a small linear combination of determinants adapted to symmetry).
- The orbitals are optimised via the **variational principle**.
- Electron-electron repulsion enters through the **Hartree (direct/Coulomb) potential** and the **Fock (exchange) operator**, the latter arising purely from the antisymmetry of the wavefunction.
- The resulting equations are **nonlinear** and must be solved **self-consistently**.

For atoms, the spherical symmetry of the nuclear potential enables a dramatic simplification: orbitals factorise into radial functions times spherical harmonics, and the problem reduces to coupled **one-dimensional radial integro-differential equations**.

Atomic units are used throughout: $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$, so energies are in Hartree ($E_h$) and lengths in Bohr radii ($a_0$).

---

## 2. The Atomic Hamiltonian

### 2.1 Non-relativistic, Born-Oppenheimer form

For an atom with an infinitely heavy point nucleus of charge $Z$ at the origin:

$$
\hat{H} = \sum_{i=1}^{N} \hat{h}(\mathbf{r}_i) + \sum_{i<j}^{N} \frac{1}{r_{ij}}
$$

with the one-electron operator

$$
\hat{h}(\mathbf{r}) = -\frac{1}{2}\nabla^2 - \frac{Z}{r}
$$

and $r_{ij} = |\mathbf{r}_i - \mathbf{r}_j|$.

### 2.2 Symmetries

$\hat{H}$ commutes with:

- Total orbital angular momentum $\hat{\mathbf{L}}^2$ and $\hat{L}_z$
- Total spin $\hat{\mathbf{S}}^2$ and $\hat{S}_z$
- Parity $\hat{\Pi}$
- Permutation of electron coordinates

Hence exact atomic states are labelled by the **term symbol** $^{2S+1}L_\pi$ (e.g., $^3P^o$ for carbon ground state).

### 2.3 Neglected effects

The non-relativistic Hamiltonian omits: spin-orbit coupling, Darwin and mass-velocity terms, finite nuclear size, nuclear motion (mass polarisation, finite-mass correction), QED corrections (Lamb shift), and hyperfine interactions. Some are recovered in Section 18.

---

## 3. The Variational Principle

For any normalised trial function $\Psi_\text{trial}$ in the appropriate Hilbert space:

$$
E[\Psi_\text{trial}] = \frac{\langle \Psi_\text{trial} | \hat{H} | \Psi_\text{trial} \rangle}{\langle \Psi_\text{trial} | \Psi_\text{trial} \rangle} \ge E_0
$$

where $E_0$ is the exact ground-state energy. Equality holds only if $\Psi_\text{trial}$ is the exact ground state. HF minimises $E[\Psi]$ over the manifold of single Slater determinants, yielding the **best single-determinant approximation** — and an **upper bound** to $E_0$.

---

## 4. The Slater Determinant Ansatz

### 4.1 Spin-orbitals

Each electron occupies a **spin-orbital**:

$$
\chi_a(\mathbf{x}) = \phi_a(\mathbf{r})\,\sigma_a(s), \qquad \sigma_a \in \{\alpha, \beta\}
$$

where $\mathbf{x} = (\mathbf{r}, s)$ combines spatial and spin coordinates. Orthonormality is imposed:

$$
\langle \chi_a | \chi_b \rangle = \delta_{ab}
$$

### 4.2 The determinant

$$
\Psi(\mathbf{x}_1, \dots, \mathbf{x}_N) = \frac{1}{\sqrt{N!}}
\begin{vmatrix}
\chi_1(\mathbf{x}_1) & \chi_2(\mathbf{x}_1) & \cdots & \chi_N(\mathbf{x}_1) \\
\chi_1(\mathbf{x}_2) & \chi_2(\mathbf{x}_2) & \cdots & \chi_N(\mathbf{x}_2) \\
\vdots & \vdots & \ddots & \vdots \\
\chi_1(\mathbf{x}_N) & \chi_2(\mathbf{x}_N) & \cdots & \chi_N(\mathbf{x}_N)
\end{vmatrix}
$$

This automatically satisfies the **Pauli principle**: $\Psi$ changes sign under exchange of any two electrons and vanishes if two spin-orbitals coincide.

### 4.3 Invariance

The energy is invariant under any unitary transformation among occupied spin-orbitals (the determinant changes only by a phase). This freedom permits the choice of **canonical orbitals** that diagonalise the Fock operator.

### 4.4 Spin-orbital notation for atoms

Atomic orbitals are labelled by quantum numbers $(n, l, m_l, m_s)$:

$$
\chi_{n l m_l m_s}(\mathbf{x}) = \frac{P_{nl}(r)}{r}\, Y_{l m_l}(\theta,\varphi)\, \sigma_{m_s}(s)
$$

with $P_{nl}(r) = r R_{nl}(r)$ the **reduced radial function**, normalised as

$$
\int_0^\infty P_{nl}^2(r)\,dr = 1
$$

---

## 5. Energy Expectation Value

Using the Slater-Condon rules, the energy of a single determinant is

$$
E_\text{HF} = \sum_{a=1}^{N} \langle a | \hat{h} | a \rangle + \frac{1}{2}\sum_{a=1}^{N}\sum_{b=1}^{N} \Big( \langle ab | ab \rangle - \langle ab | ba \rangle \Big)
$$

where

$$
\langle ab | cd \rangle = \iint \chi_a^\ast(\mathbf{x}_1)\chi_b^\ast(\mathbf{x}_2)\,\frac{1}{r_{12}}\,\chi_c(\mathbf{x}_1)\chi_d(\mathbf{x}_2)\,d\mathbf{x}_1\,d\mathbf{x}_2
$$

Define:

- **Coulomb integral** $J_{ab} = \langle ab | ab \rangle$
- **Exchange integral** $K_{ab} = \langle ab | ba \rangle$ (non-zero only if $a$ and $b$ have the same spin)

Then

$$
E_\text{HF} = \sum_a h_{aa} + \frac{1}{2}\sum_{a,b}\left( J_{ab} - K_{ab} \right)
$$

The $a=b$ terms cancel ($J_{aa} = K_{aa}$), so **HF is free of one-electron self-interaction**.

---

## 6. Derivation of the Hartree-Fock Equations

### 6.1 Constrained variation

Minimise $E_\text{HF}$ subject to orthonormality using Lagrange multipliers $\varepsilon_{ab}$:

$$
\mathcal{L} = E_\text{HF} - \sum_{a,b} \varepsilon_{ba}\left( \langle a | b \rangle - \delta_{ab} \right)
$$

Setting $\delta\mathcal{L} = 0$ with respect to variations $\delta\chi_a^\ast$ gives

$$
\hat{f}\,\chi_a(\mathbf{x}) = \sum_b \varepsilon_{ba}\,\chi_b(\mathbf{x})
$$

### 6.2 Canonical form

Because $\varepsilon_{ab}$ is Hermitian, a unitary transformation of the occupied orbitals diagonalises it, giving the **canonical Hartree-Fock equations**:

$$
\boxed{\hat{f}\,\chi_a(\mathbf{x}) = \varepsilon_a\,\chi_a(\mathbf{x})}
$$

### 6.3 The Fock operator

$$
\hat{f}(\mathbf{x}_1) = \hat{h}(\mathbf{r}_1) + \sum_{b=1}^{N}\left[ \hat{J}_b(\mathbf{x}_1) - \hat{K}_b(\mathbf{x}_1) \right]
$$

with

**Coulomb operator** (local):

$$
\hat{J}_b(\mathbf{x}_1)\,\chi_a(\mathbf{x}_1) = \left[\int \frac{|\chi_b(\mathbf{x}_2)|^2}{r_{12}}\,d\mathbf{x}_2\right]\chi_a(\mathbf{x}_1)
$$

**Exchange operator** (non-local):

$$
\hat{K}_b(\mathbf{x}_1)\,\chi_a(\mathbf{x}_1) = \left[\int \frac{\chi_b^\ast(\mathbf{x}_2)\,\chi_a(\mathbf{x}_2)}{r_{12}}\,d\mathbf{x}_2\right]\chi_b(\mathbf{x}_1)
$$

Note that $\hat{K}_b$ acting on $\chi_a$ produces a function proportional to $\chi_b$, not $\chi_a$ — hence "non-local".

### 6.4 Self-consistency

The Fock operator depends on its own eigenfunctions, so the equations are **nonlinear**. They are solved by iteration until the input and output orbitals agree: the **self-consistent field (SCF)** method.

---

## 7. Central-Field Approximation and Separation of Variables

### 7.1 Spherical averaging

For closed subshells the total charge density is spherically symmetric, so the effective potential is central. For open subshells, one adopts the **central-field approximation** — the Fock operator is replaced by its spherical average over the angular coordinates of the occupied shell (or, in a more rigorous treatment, the energy functional is derived from the LS-coupled term).

### 7.2 Orbital form

With a central potential $V(r)$:

$$
\phi_{nlm_l}(\mathbf{r}) = \frac{P_{nl}(r)}{r}\,Y_{lm_l}(\theta,\varphi)
$$

and the effective one-electron equation becomes

$$
\left[ -\frac{1}{2}\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} + V(r) \right] P_{nl}(r) = \varepsilon_{nl}\,P_{nl}(r)
$$

with boundary conditions

$$
P_{nl}(0) = 0, \qquad P_{nl}(r) \xrightarrow{r\to\infty} 0
$$

### 7.3 Asymptotic behaviour

- Near the origin: $P_{nl}(r) \sim r^{l+1}$
- At large $r$: $P_{nl}(r) \sim e^{-\sqrt{-2\varepsilon_{nl}}\,r}$ (bound orbitals, $\varepsilon_{nl}<0$), modulated by a power-law prefactor set by the net charge seen by the outermost electron.

---

## 8. Angular Algebra: Reduction to Radial Equations

### 8.1 Multipole expansion of $1/r_{12}$

$$
\frac{1}{r_{12}} = \sum_{k=0}^{\infty} \frac{r_<^k}{r_>^{k+1}}\,P_k(\cos\theta_{12})
= \sum_{k=0}^{\infty}\frac{r_<^k}{r_>^{k+1}}\,\frac{4\pi}{2k+1}\sum_{q=-k}^{k} Y_{kq}^\ast(\Omega_1)\,Y_{kq}(\Omega_2)
$$

with $r_< = \min(r_1,r_2)$ and $r_> = \max(r_1,r_2)$.

### 8.2 Slater radial integrals

Define

$$
R^k(ab; cd) = \int_0^\infty\!\!\int_0^\infty P_a(r_1)P_b(r_2)\,\frac{r_<^k}{r_>^{k+1}}\,P_c(r_1)P_d(r_2)\,dr_1\,dr_2
$$

Special cases:

- **Direct**: $F^k(n_al_a, n_bl_b) = R^k(ab;ab)$
- **Exchange**: $G^k(n_al_a, n_bl_b) = R^k(ab;ba)$
- For equivalent electrons: $F^k(nl,nl)$

### 8.3 Gaunt coefficients (angular integrals)

Angular integrals reduce to products of **Wigner 3-j symbols**:

$$
\langle l\,m\,|\,C^{(k)}_q\,|\,l'\,m'\rangle
= (-1)^m \sqrt{(2l+1)(2l'+1)}
\begin{pmatrix} l & k & l' \\ -m & q & m' \end{pmatrix}
\begin{pmatrix} l & k & l' \\ 0 & 0 & 0 \end{pmatrix}
$$

where $C^{(k)}_q = \sqrt{4\pi/(2k+1)}\,Y_{kq}$ is the reduced spherical harmonic. The 3-j symbol with three zeros enforces the **parity selection rule**: $l + k + l'$ must be even.

### 8.4 Direct and exchange angular coefficients

After summing over magnetic quantum numbers, the energy of a configuration takes the form

$$
E = \sum_a q_a I(a) + \sum_{a}\sum_{k>0} f_k(aa)\,F^k(a,a) + \sum_{a<b} \left( \sum_{k} f_k(ab) \, F^k(a,b) + \sum_k g_k(ab) \, G^k(a,b) \right)
$$

where $q_a$ is the occupation number, $I(a)$ the one-electron radial integral

$$
I(nl) = \int_0^\infty P_{nl}(r)\left[-\frac{1}{2}\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} - \frac{Z}{r}\right]P_{nl}(r)\,dr
$$

and $f_k$, $g_k$ are **angular coefficients** determined by the coupling scheme (LS, jj, etc.).

For a **filled subshell**, only $k=0$ direct terms survive after angular averaging, giving spherical symmetry.

---

## 9. The Radial Hartree-Fock Equations

Varying $E$ with respect to $P_{nl}$ subject to $\langle P_{nl}|P_{n'l}\rangle = \delta_{nn'}$ gives, for orbital $a \equiv nl$:

$$
\boxed{
\left[ -\frac{1}{2}\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} - \frac{Z}{r} + Y_a(r) \right]P_a(r)
- X_a(r) = \varepsilon_a\,P_a(r) + \sum_{b\neq a}\delta_{l_al_b}\,\varepsilon_{ab}\,P_b(r)
}
$$

### 9.1 Direct potential (Hartree term)

$$
Y_a(r) = \sum_b q_b \sum_k \frac{f_k(ab)}{q_a}\,Y^k(bb;r)/r
$$

where the **Hartree-Slater potential function**

$$
Y^k(bb;r) = r\int_0^\infty \frac{r_<^k}{r_>^{k+1}}\,P_b^2(r')\,dr'
= \frac{1}{r^k}\int_0^r P_b^2(r')\,r'^k\,dr' + r^{k+1}\int_r^\infty \frac{P_b^2(r')}{r'^{k+1}}\,dr'
$$

For $k=0$ this is the electrostatic potential of the spherically symmetric charge distribution of shell $b$, multiplied by $r$:

$$
Y^0(bb;r) = \frac{1}{1}\left[\int_0^r P_b^2(r')\,dr' + r\int_r^\infty \frac{P_b^2(r')}{r'}\,dr'\right]
$$

### 9.2 Exchange term

$$
X_a(r) = \sum_{b\neq a}\sum_k \frac{g_k(ab)\,q_b}{q_a}\,\frac{Y^k(ab;r)}{r}\,P_b(r)
$$

with

$$
Y^k(ab;r) = r\int_0^\infty \frac{r_<^k}{r_>^{k+1}}\,P_a(r')P_b(r')\,dr'
$$

The exchange term is an **integro-differential** contribution: the value of $X_a$ at $r$ depends on $P_a$ and $P_b$ everywhere.

### 9.3 Off-diagonal Lagrange multipliers

The terms with $\varepsilon_{ab}$ enforce orthogonality between orbitals with the same $l$ but different $n$. They are non-zero when orbitals are not automatically orthogonal by symmetry. In canonical HF for closed shells they can be eliminated; for open shells they generally remain.

### 9.4 Boundary-value problem

The equations form a set of coupled nonlinear boundary-value problems: for each occupied $(n,l)$ one solves for $P_{nl}(r)$ and eigenvalue $\varepsilon_{nl}$ with $P(0)=0$, $P(\infty)=0$, and normalisation.

---

## 10. Closed-Shell Atoms

### 10.1 Restricted HF (RHF)

For closed-shell atoms (He, Be, Ne, Mg, Ar, ...), each spatial orbital is doubly occupied with opposite spins. Spin summation gives

$$
\hat{f} = \hat{h} + \sum_{b=1}^{N/2}\left(2\hat{J}_b - \hat{K}_b\right)
$$

with spatial-orbital Coulomb and exchange operators. The total wavefunction is a pure singlet $^1S$.

### 10.2 Energy

$$
E_\text{RHF} = 2\sum_{a=1}^{N/2} h_{aa} + \sum_{a=1}^{N/2}\sum_{b=1}^{N/2}\left(2J_{ab} - K_{ab}\right)
$$

### 10.3 Spherical symmetry

Every filled $(n,l)$ subshell has a spherically symmetric density:

$$
\rho_{nl}(\mathbf{r}) = \frac{2(2l+1)}{4\pi}\frac{P_{nl}^2(r)}{r^2}
$$

The total density is $\rho(r) = \sum_{nl}\rho_{nl}(r)$, spherically symmetric, giving a central Hartree potential and central (but non-local) exchange.

### 10.4 Helium as the simplest case

For He ($1s^2$, $N=2$), the exchange operator acting on $1s$ equals the Coulomb operator, so the HF equation reduces to

$$
\left[-\frac{1}{2}\nabla^2 - \frac{2}{r} + \hat{J}_{1s}\right]\phi_{1s} = \varepsilon_{1s}\,\phi_{1s}
$$

which is the **Hartree equation** with self-interaction removed: an electron feels the other electron's charge cloud only.

---

## 11. Open-Shell Atoms

### 11.1 Unrestricted HF (UHF)

Different spatial orbitals for $\alpha$ and $\beta$ spins:

$$
\phi_a^\alpha(\mathbf{r}) \neq \phi_a^\beta(\mathbf{r})
$$

Yields lower energy than RHF at the cost of **spin contamination**: $\Psi_\text{UHF}$ is not an eigenfunction of $\hat{S}^2$. The Pople-Nesbet equations give two coupled Fock matrices:

$$
\hat{f}^\alpha = \hat{h} + \hat{J}[\rho^\alpha + \rho^\beta] - \hat{K}[\rho^\alpha], \qquad
\hat{f}^\beta = \hat{h} + \hat{J}[\rho^\alpha + \rho^\beta] - \hat{K}[\rho^\beta]
$$

Spin contamination is measured by

$$
\langle \hat{S}^2\rangle_\text{UHF} - S(S+1)
$$

### 11.2 Restricted Open-shell HF (ROHF)

Doubly occupied core orbitals and singly occupied open-shell orbitals with the **same spatial function for both spins in doubly-occupied shells**. Maintains spin purity but requires more complicated coupling operators (Roothaan's, or Guest-Saunders effective Fock operator).

### 11.3 Multi-Configuration HF (MCHF) and term-dependent HF

For atoms, the physically appropriate treatment of an open shell $nl^q$ is via the **LS term** (e.g., $^3P$ for $np^2$). The energy of a term is generally **not** expressible as a single determinant; it is a linear combination of a few determinants fixed by angular momentum coupling:

$$
\Psi(\gamma LS) = \sum_{i} c_i\,\Phi_i, \qquad c_i \text{ fixed by symmetry (Clebsch-Gordan)}
$$

The energy of a term is expressed via Slater integrals with tabulated coefficients:

| Configuration | Term | Energy expression (excerpt) |
|---|---|---|
| $1s^2$ | $^1S$ | $2I(1s) + F^0(1s,1s)$ |
| $2p^2$ | $^3P$ | $2I(2p) + F^0(2p,2p) - \tfrac{5}{25}F^2(2p,2p)$ |
| $2p^2$ | $^1D$ | $2I(2p) + F^0(2p,2p) + \tfrac{1}{25}F^2(2p,2p)$ |
| $2p^2$ | $^1S$ | $2I(2p) + F^0(2p,2p) + \tfrac{10}{25}F^2(2p,2p)$ |

The **Hund's-rule** ordering $^3P < ^1D < ^1S$ follows directly from the sign and size of the $F^2$ coefficient.

### 11.4 Average-of-configuration approximation

When one does not want to specify a term, one averages over all states of the configuration. This gives a **spherical potential** and simple $F^k$ coefficients:

$$
E_\text{av}(nl^q) = q\,I(nl) + \frac{q(q-1)}{2}\left[F^0(nl,nl) - \frac{2l+1}{4l+1}\sum_{k>0}\begin{pmatrix} l & k & l\\ 0&0&0\end{pmatrix}^2 F^k(nl,nl)\right]
$$

This is the starting point for most atomic-structure codes (Froese Fischer, Cowan).

---

## 12. The Self-Consistent Field (SCF) Procedure

### 12.1 Algorithm

1. **Initial guess**: hydrogenic orbitals with screened charge $Z_\text{eff}$, Thomas-Fermi or Slater-type orbitals, or Hartree-Slater.
2. **Build potentials**: compute $Y_a(r)$ and $X_a(r)$ from current orbitals.
3. **Solve radial equations**: for each $(n,l)$, solve the linear eigenvalue problem with fixed potentials.
4. **Orthogonalise**: enforce $\langle P_{nl}|P_{n'l}\rangle=0$ (Gram-Schmidt or via Lagrange multipliers).
5. **Mix**: $P^{(i+1)} = (1-\alpha)P^{(i)} + \alpha P^{(i+1)}_\text{new}$ to stabilise convergence.
6. **Check convergence**: $\Delta E < \epsilon_E$ and $\max\|P^{(i+1)}-P^{(i)}\|<\epsilon_P$.
7. **Repeat** from step 2 until converged.

### 12.2 Pseudocode

```
initialise P_nl(r) for all occupied (n,l)
repeat
    for each occupied (n,l):
        compute Y_nl(r), X_nl(r) from current {P}
        solve  [-1/2 d²/dr² + l(l+1)/2r² - Z/r + Y] P_new - X = ε P_new
        normalise P_new
    orthogonalise {P_new} within each l
    P <- (1-α) P + α P_new
    E <- total energy from {P}
until |ΔE| < tol and max|ΔP| < tol
```

### 12.3 Convergence acceleration

- **Damping** (simple mixing) with $\alpha \sim 0.3$–$0.7$.
- **Direct Inversion in the Iterative Subspace (DIIS)**: extrapolate using stored error vectors
  $$ \mathbf{e}_i = \mathbf{F}_i\mathbf{D}_i\mathbf{S} - \mathbf{S}\mathbf{D}_i\mathbf{F}_i $$
- **Level shifting**: add a constant to virtual orbital energies.
- **Newton-Raphson / quasi-Newton** second-order SCF for stubborn open-shell cases.

### 12.4 Convergence pathologies

- **Oscillation** between two solutions (charge sloshing) — cured by damping.
- **Convergence to saddle points** or excited-state solutions — use stability analysis of the Hessian.
- **Multiple solutions**: HF energy functional is non-convex; different initial guesses can give different self-consistent solutions.

---

## 13. Numerical Solution of the Radial Equations

### 13.1 Grid

Radial equations are solved on a **logarithmic** or **exponential** grid to resolve the steep cusp near the nucleus and the long tail:

$$
r_i = r_0\left(e^{i h} - 1\right), \qquad \text{or} \qquad t = \ln r
$$

Change of variable: $P(r) \to P(t)$, with $r = e^{t}$, gives

$$
\frac{d^2P}{dr^2} = e^{-2t}\left(\frac{d^2P}{dt^2} - \frac{dP}{dt}\right)
$$

### 13.2 Finite-difference / Numerov method

The Numerov algorithm for $y'' = -g(x)\,y$ yields fourth-order accuracy:

$$
y_{i+1}\left(1 + \frac{h^2}{12}g_{i+1}\right) = 2y_i\left(1 - \frac{5h^2}{12}g_i\right) - y_{i-1}\left(1 + \frac{h^2}{12}g_{i-1}\right)
$$

### 13.3 Shooting and matching

1. Integrate **outward** from $r\approx 0$ using $P\sim r^{l+1}$.
2. Integrate **inward** from $r_\text{max}$ using the decaying asymptotic.
3. Match at the classical turning point $r_c$.
4. Adjust $\varepsilon$ until the logarithmic derivatives agree (node counting fixes $n$: $P_{nl}$ has $n-l-1$ nodes).

### 13.4 Exchange as an inhomogeneous term

In the iterative scheme the exchange term $X_a(r)$ is treated as a known inhomogeneity, turning each radial equation into a **linear inhomogeneous boundary-value problem**, solved by, e.g., the **deferred-correction** or **Green's function** methods (Froese Fischer).

### 13.5 Poisson equation for direct potentials

The direct potential $Y^0/r$ can alternatively be found by solving Poisson's equation for the spherically symmetric density:

$$
\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{dV_H}{dr}\right) = -4\pi\rho(r)
$$

or, in reduced form with $U(r) = rV_H(r)$:

$$
\frac{d^2U}{dr^2} = -4\pi r\rho(r)
$$

with $U(0)=0$, $U(\infty)=N$.

### 13.6 Spline and B-spline methods

Modern approaches expand $P_{nl}(r)$ in **B-splines** on a non-uniform knot sequence, converting the problem into a generalised matrix eigenproblem with excellent convergence and natural handling of continuum states.

---

## 14. Basis-Set Expansion: Roothaan-Hartree-Fock

### 14.1 Expansion

Expand each radial function in a finite basis:

$$
P_{nl}(r) = \sum_{p=1}^{M_l} C_{p}^{(nl)}\,\chi_p^{(l)}(r)
$$

### 14.2 Slater-type orbitals (STOs)

$$
\chi_p^{(l)}(r) = N_p\,r^{n_p}\,e^{-\zeta_p r}, \qquad n_p \ge l+1
$$

with normalisation $N_p = (2\zeta_p)^{n_p+1/2}/\sqrt{(2n_p)!}$. STOs have the correct cusp and exponential decay, and are the traditional choice for high-precision atomic HF (Clementi-Roetti tables, Bunge et al.).

### 14.3 Gaussian-type orbitals (GTOs)

$$
\chi_p^{(l)}(r) = N_p\,r^{l+1}\,e^{-\alpha_p r^2}
$$

GTOs lack the nuclear cusp and have wrong asymptotics but make two-electron integrals analytic. Contractions $\sum_k d_k\,\chi_k$ recover STO-like behaviour.

### 14.4 Roothaan equations

Substituting into the HF equation and projecting gives the **generalised matrix eigenproblem** for each symmetry $l$:

$$
\boxed{\mathbf{F}^{(l)}\,\mathbf{C}^{(l)} = \mathbf{S}^{(l)}\,\mathbf{C}^{(l)}\,\boldsymbol{\varepsilon}^{(l)}}
$$

with

$$
S_{pq} = \int_0^\infty \chi_p\,\chi_q\,dr, \qquad
F_{pq} = h_{pq} + \sum_{rs}\left[D_{rs}\,(pq|rs) - \tfrac{1}{2}D_{rs}\,(pr|qs)\right]\quad\text{(closed shell)}
$$

and density matrix $D_{rs} = 2\sum_{a}^\text{occ}C_{ra}C_{sa}$.

### 14.5 Even-tempered and well-tempered bases

Exponents chosen by geometric progression:

$$
\zeta_p = \alpha\,\beta^{p-1}, \qquad p = 1, \dots, M
$$

with two parameters $(\alpha,\beta)$ per symmetry. **Well-tempered** bases add a mild correction to $\ln\zeta_p$. These reach the **Hartree-Fock limit** to $\mu E_h$ with $M\sim 10$–$20$ per $l$.

### 14.6 The Hartree-Fock limit

The HF limit is the energy at complete basis. Accurate atomic values (numerical HF):

| Atom | $E_\text{HF}$ ($E_h$) |
|---|---|
| H | −0.500000 |
| He | −2.861680 |
| Li | −7.432727 |
| Be | −14.573023 |
| B | −24.529061 |
| C | −37.688619 |
| N | −54.400934 |
| O | −74.809398 |
| F | −99.409349 |
| Ne | −128.547098 |

---

## 15. Orbital Energies and Koopmans' Theorem

### 15.1 Meaning of $\varepsilon_a$

Multiplying $\hat f\chi_a = \varepsilon_a\chi_a$ on the left by $\chi_a^\ast$ and integrating:

$$
\varepsilon_a = h_{aa} + \sum_{b}\left(J_{ab} - K_{ab}\right)
$$

Thus **$\sum_a\varepsilon_a \neq E_\text{HF}$**; rather

$$
E_\text{HF} = \sum_a \varepsilon_a - \frac{1}{2}\sum_{a,b}\left(J_{ab}-K_{ab}\right)
$$

(the electron-electron repulsion is double-counted in $\sum\varepsilon_a$).

### 15.2 Koopmans' theorem

If orbitals are **frozen** upon removal of electron $a$, the ionisation energy is

$$
\text{IE}_a \approx -\varepsilon_a
$$

Proof sketch: compute $E_N - E_{N-1}^{(a)}$ using the same determinant with column $a$ removed; algebra gives $-\varepsilon_a$. Errors:

- Neglect of **orbital relaxation** (makes IE too large).
- Neglect of **correlation** (usually makes IE too small, partial cancellation).

For valence electrons the two errors largely cancel; for core electrons, relaxation dominates and Koopmans overestimates by several eV.

### 15.3 Brillouin's theorem

Singly excited determinants $\Phi_a^r$ do not couple directly to the HF ground state:

$$
\langle \Phi_0 | \hat H | \Phi_a^r \rangle = 0
$$

This is a direct consequence of the stationarity condition and has the important consequence that **the leading correction to the HF energy comes from doubly excited determinants**.

### 15.4 ΔSCF

Ionisation energies and excitation energies computed as differences of separate SCF total energies (ΔSCF) recover orbital relaxation and are typically much more accurate than Koopmans values.

---

## 16. Worked Examples

### 16.1 Helium ($1s^2$)

- $E_\text{HF} = -2.861680\,E_h$
- Exact non-relativistic: $E_\text{exact} = -2.903724\,E_h$
- Correlation energy: $E_\text{corr} = -0.042044\,E_h$ (about 1.4% of total, 1.14 eV)
- Orbital energy: $\varepsilon_{1s} = -0.917956\,E_h$ (Koopmans IE: 24.98 eV vs. experiment 24.59 eV)

A single-STO variational estimate uses $\phi(r)\propto e^{-\zeta r}$; minimising gives $\zeta = Z - 5/16 = 27/16$ and $E = -2.84766\,E_h$.

### 16.2 Beryllium ($1s^22s^2$)

- $E_\text{HF} = -14.573023\,E_h$
- $E_\text{exact} = -14.667356\,E_h$
- $E_\text{corr} = -0.094333\,E_h$
- Near-degeneracy: the $2s^2 \leftrightarrow 2p^2$ mixing gives a strong **static correlation** contribution that a single determinant cannot capture; MCHF with $(2s^2 + 2p^2)$ recovers about 80% of $E_\text{corr}$.

### 16.3 Neon ($1s^22s^22p^6$)

- $E_\text{HF} = -128.547098\,E_h$
- $E_\text{exact} = -128.9376\,E_h$
- $E_\text{corr} = -0.3905\,E_h$
- Closed-shell, spherically symmetric, ideal for HF. Orbital energies: $\varepsilon_{1s}=-32.772$, $\varepsilon_{2s}=-1.9304$, $\varepsilon_{2p}=-0.8504$ $E_h$.

### 16.4 Carbon $^3P$ ($1s^22s^22p^2$)

- Open-shell; energy of the $^3P$ term via Slater integrals with $f_k$, $g_k$ tabulated.
- $E_\text{HF}(^3P) = -37.688619\,E_h$.
- Excited terms $^1D$ and $^1S$ lie higher by ~1.26 and ~2.68 eV (HF) and ~1.26 and ~2.68 eV (experiment) respectively — **Hund's rule** follows from exchange.

---

## 17. Correlation Energy and Limitations

### 17.1 Definition

$$
E_\text{corr} \equiv E_\text{exact}^\text{non-rel} - E_\text{HF}
$$

Always negative by the variational principle. Löwdin's definition.

### 17.2 Types of correlation

| Type | Origin | Example |
|---|---|---|
| **Dynamical** | Short-range Coulomb cusp, instantaneous avoidance of electrons | He, Ne |
| **Static / non-dynamical** | Near-degeneracy of configurations | Be ($2s$–$2p$), C |
| **Core-valence** | Polarisation of core by valence | Alkali atoms |

### 17.3 The Coulomb cusp

The exact wavefunction satisfies Kato's cusp condition

$$
\left.\frac{\partial\bar\Psi}{\partial r_{12}}\right|_{r_{12}=0} = \frac{1}{2}\Psi(r_{12}=0)
$$

A determinant of orbitals cannot reproduce this since it is smooth at $r_{12}=0$ (apart from same-spin nodes), so HF describes the **Coulomb hole** poorly. The **Fermi hole** (same-spin exclusion) is treated exactly.

### 17.4 Typical errors of HF

- Total energy: ~0.5–1% error, but comparable in magnitude to chemically relevant energies.
- Ionisation energies: ~5–10% (Koopmans), better with ΔSCF.
- Electron affinities: poor — anions often unbound at HF (e.g., Be$^-$ unstable), because correlation contributes a large portion of binding.
- Dipole polarisabilities: overestimated by up to 20–40% for alkali atoms.
- Excitation energies: HF cannot describe excited states with the same determinant.

### 17.5 Size-consistency and symmetry breaking

- HF is size-consistent for closed-shell fragments.
- Symmetry-broken solutions (UHF lower than RHF) indicate the RHF single-determinant is unstable — **Hartree-Fock instability**.
- Spherical symmetry of an atom's HF orbitals is enforced by hand in radial codes; relaxing it (**deformed HF**) may lower the energy but violates angular-momentum eigenstate character.

---

## 18. Relativistic Extension: Dirac-Hartree-Fock

For heavy atoms (Z ≳ 30), relativistic effects become significant (scaling as $(Z\alpha)^2$).

### 18.1 Dirac-Coulomb Hamiltonian

$$
\hat H_\text{DC} = \sum_i\left[c\,\boldsymbol\alpha\cdot\mathbf{p}_i + \beta m c^2 - \frac{Z}{r_i}\right] + \sum_{i<j}\frac{1}{r_{ij}}
$$

### 18.2 Four-component orbitals

$$
\psi_{n\kappa m}(\mathbf{r}) = \frac{1}{r}
\begin{pmatrix}
P_{n\kappa}(r)\,\Omega_{\kappa m}(\theta,\varphi)\\
i\,Q_{n\kappa}(r)\,\Omega_{-\kappa m}(\theta,\varphi)
\end{pmatrix}
$$

with **relativistic quantum number** $\kappa = \mp(j+1/2)$ for $j = l\pm 1/2$, and $P$, $Q$ the large and small radial components.

### 18.3 Radial Dirac-Fock equations

$$
\begin{aligned}
\left(\frac{d}{dr} + \frac{\kappa}{r}\right)P_a - \left[\frac{2}{\alpha}+\alpha(\varepsilon_a - V_a)\right]Q_a &= \alpha\,X^Q_a\\
\left(\frac{d}{dr} - \frac{\kappa}{r}\right)Q_a + \alpha(\varepsilon_a - V_a)P_a &= -\alpha\,X^P_a
\end{aligned}
$$

with $\alpha \approx 1/137.036$ and $V_a$ the sum of nuclear and direct potentials. Coupled first-order equations.

### 18.4 Effects captured

- **Spin-orbit splitting** of $l>0$ shells into $j=l\pm1/2$.
- **Direct relativistic contraction** of $s$ and $p_{1/2}$ orbitals; indirect expansion of $d$, $f$ (through improved screening).
- Mass-velocity and Darwin terms implicitly.

### 18.5 Breit interaction and QED

Beyond Dirac-Coulomb, add the **Breit interaction** (magnetic + retardation):

$$
\hat B_{ij} = -\frac{1}{2r_{ij}}\left[\boldsymbol\alpha_i\cdot\boldsymbol\alpha_j + \frac{(\boldsymbol\alpha_i\cdot\mathbf{r}_{ij})(\boldsymbol\alpha_j\cdot\mathbf{r}_{ij})}{r_{ij}^2}\right]
$$

plus self-energy and vacuum polarisation for high-precision heavy-atom work.

### 18.6 Continuum-dissolution / negative-energy problem

The Dirac Hamiltonian is unbounded below; variational collapse is avoided using **kinetic balance** in basis expansion or the **projected (no-pair) Hamiltonian**.

---

## 19. Post-Hartree-Fock Methods

| Method | Idea | Scaling |
|---|---|---|
| **Configuration Interaction (CI)** | Linear combination of determinants | Exponential (FCI) |
| **MCHF / MCDHF** | CI + simultaneously optimised orbitals | — |
| **Møller-Plesset PT (MP2, MP4)** | Perturbation with $\hat H_0=\sum\hat f$ | $N^5$, $N^7$ |
| **Coupled Cluster (CCSD, CCSD(T))** | $\Psi=e^{\hat T}\Phi_0$ | $N^6$, $N^7$ |
| **Many-Body Perturbation Theory (MBPT)** | Diagrammatic, all-order (RMBPT) | — |
| **Explicitly correlated (R12/F12, Hylleraas)** | Include $r_{12}$ explicitly | — |
| **DFT** | Density-functional replacement for exchange-correlation | $N^3$ |

The HF determinant serves as the **reference** for essentially all systematic correlation methods; the quality of the reference (single-reference character) determines their reliability.

### 19.1 MP2 energy for an atom

$$
E^{(2)} = \frac{1}{4}\sum_{ab}^\text{occ}\sum_{rs}^\text{virt}\frac{|\langle ab||rs\rangle|^2}{\varepsilon_a+\varepsilon_b-\varepsilon_r-\varepsilon_s}
$$

For He this recovers ~ 85% of the exact correlation energy at the HF limit.

### 19.2 Coupled-cluster

$$
\hat T = \hat T_1 + \hat T_2 + \cdots, \qquad
\hat T_2 = \frac{1}{4}\sum t_{ab}^{rs}\,a_r^\dagger a_s^\dagger a_b a_a
$$

CCSD(T) is the "gold standard" for atoms and small molecules.

---

## 20. Practical Implementation Notes

### 20.1 Existing codes

- **Froese Fischer's MCHF** package (Fortran/C++)
- **GRASP / GRASP2018** (relativistic MCDHF)
- **Cowan's RCN/RCG** (average-configuration HFR, widely used in spectroscopy)
- **atomic (ATSP2K)**, **HFS**, **OPMKS/NIST**, **DIRAC**, **ADF** (basis-set).
- Python: `pyscf` (Gaussian basis), `libxc` for DFT comparisons, `atomrad`.

### 20.2 Design checklist for a radial HF code

1. Choose grid: $r_i = r_0(e^{ih}-1)$, $h\sim 0.01$–$0.03$, $r_\text{max}\sim 40$–$60/\sqrt{\varepsilon_\text{min}}$.
2. Precompute $Y^k$ integrals via cumulative sums (Simpson or Numerov-type).
3. Precompute angular coefficients $f_k, g_k$ using 3-j/6-j symbols or a table.
4. Solve radial equations with Numerov shooting + eigenvalue root-finding (bisection / secant).
5. Apply Lagrange multipliers or Gram-Schmidt for orthogonality.
6. Use DIIS or Anderson mixing on the potentials.
7. Validate on He, Ne (closed-shell) and on $^3P$ carbon (open-shell).

### 20.3 Computational complexity

For a grid of $N_r$ points and $N_o$ occupied orbitals with maximum multipole $k_\text{max}=2l_\text{max}$:

- Direct/exchange potentials: $O(N_o^2\,k_\text{max}\,N_r)$
- Radial solves: $O(N_o\,N_r)$ per iteration
- Total: cheap; atomic HF runs in milliseconds to seconds even for heavy atoms.

### 20.4 Validation benchmarks

- Total energies compared to Clementi-Roetti, Bunge-Barrientos-Bunge (1993), or Froese Fischer numerical HF.
- Virial theorem: $2\langle T\rangle = -\langle V\rangle$ ⇒ $E = -\langle T\rangle$ at the HF limit (test on $\langle T \rangle/|E| = 1$).
- Cusp condition: $P_{nl}'(0)/P_{nl}(0)$ behaviour for $s$ orbitals: $\left.\dfrac{1}{R_{ns}}\dfrac{dR_{ns}}{dr}\right|_{r=0} = -Z$.

---

## 21. Summary

The Hartree-Fock method for a single atom:

1. Approximates the $N$-electron wavefunction by a Slater determinant (or symmetry-adapted combination) of orbitals of the form $P_{nl}(r)Y_{lm}\sigma$.
2. Minimises the energy variationally, leading to coupled radial integro-differential equations with direct (Hartree) and exchange (Fock) potentials.
3. Exploits central symmetry to reduce the problem to one dimension, using Slater integrals $F^k$, $G^k$ and angular coefficients from Racah algebra.
4. Is solved by self-consistent iteration, either on a numerical grid or in a Roothaan basis expansion.
5. Provides an upper bound to the ground-state energy, orbital energies that (via Koopmans) approximate ionisation energies, and an excellent reference for correlated methods.
6. Misses the **correlation energy** ($\sim 0.04$–$1\,E_h$ across the periodic table), and must be extended to relativistic (Dirac-HF) form for heavy elements.

---

## 22. References

1. D. R. Hartree, *The Wave Mechanics of an Atom with a Non-Coulomb Central Field*, Proc. Camb. Phil. Soc. **24**, 89, 111, 426 (1928).
2. V. Fock, *Näherungsmethode zur Lösung des quantenmechanischen Mehrkörperproblems*, Z. Phys. **61**, 126 (1930).
3. J. C. Slater, *Note on Hartree's Method*, Phys. Rev. **35**, 210 (1930); *The Theory of Complex Spectra*, Phys. Rev. **34**, 1293 (1929).
4. C. C. J. Roothaan, *New Developments in Molecular Orbital Theory*, Rev. Mod. Phys. **23**, 69 (1951).
5. D. R. Hartree, *The Calculation of Atomic Structures*, Wiley (1957).
6. C. Froese Fischer, *The Hartree-Fock Method for Atoms*, Wiley (1977).
7. E. Clementi and C. Roetti, *Roothaan-Hartree-Fock Atomic Wavefunctions*, At. Data Nucl. Data Tables **14**, 177 (1974).
8. R. D. Cowan, *The Theory of Atomic Structure and Spectra*, Univ. California Press (1981).
9. A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996).
10. I. P. Grant, *Relativistic Quantum Theory of Atoms and Molecules*, Springer (2007).
11. T. Helgaker, P. Jørgensen, J. Olsen, *Molecular Electronic-Structure Theory*, Wiley (2000).
12. T. Koopmans, *Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms*, Physica **1**, 104 (1934).
13. P.-O. Löwdin, *Quantum Theory of Many-Particle Systems*, Phys. Rev. **97**, 1474 (1955).
14. T. Kato, *On the Eigenfunctions of Many-Particle Systems in Quantum Mechanics*, Commun. Pure Appl. Math. **10**, 151 (1957).
15. C. F. Bunge, J. A. Barrientos, A. V. Bunge, *Roothaan-Hartree-Fock Ground-State Atomic Wave Functions*, At. Data Nucl. Data Tables **53**, 113 (1993).

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive description of Hartree-Fock method for the singe atom. Show the output in Markdown format. Do not copy the output of the exported files into the chat.


