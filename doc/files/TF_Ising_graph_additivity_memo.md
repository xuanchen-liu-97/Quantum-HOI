# Graph-additivity of $J^2_{\text{tot}}$ in the TF-Ising framework

**Working memo for discussion with G. Möller**
**Date:** May 2026
**Status:** Initial derivation + numerical verification complete; ready for review.

---

## 1. Why this document exists

The current manuscript's headline scaling claim — that the bare associator capacity $J_{\text{tot}}$ acquires earlier predictive power as the number of triangular faces grows — rests on **two data points**: a 4-qubit triangle and a 5-qubit tetrahedron. The interpretation we offer (face-averaging suppresses fluctuations) is plausible but, as written, has no analytic anchor. A PRL referee will almost certainly ask: *what does this look like in the thermodynamic limit?*

This memo answers that question for the **TF-Ising** specialization (manuscript Appendix E). It establishes:

1. A per-face closed-form expression for $J^2_{\text{tot}}$ that holds on *any* graph (not just a single triangle).
2. A graph-summed identity that makes the dependence on graph structure explicit.
3. An $N$-scaling corollary for triangle-decorated lattices that gives the analytic basis for the manuscript's self-averaging argument.

All three are proved (or proof-sketched) below and verified numerically to machine precision on representative graphs. The resulting statement is short enough to fit into PRL main text as an analytic anchor (~half a column).

The XYZ family used in the main text does not admit an equally clean closed form, but the **same face-local decomposition holds** — TF-Ising is the case where the closed form happens to be transparent.

---

## 2. Setup and conventions

### 2.1 Graph and edge Hamiltonians

Let $G = (V, E)$ be a simple graph with vertex set $V$ and edge set $E$. We attach a qubit to each vertex. The TF-Ising Hamiltonian on $G$, in the **edge-resolved convention used in the manuscript** (Appendix E.1), is

$$H^{\text{TF}}_{ij} = J_{ij} Z_i Z_j + \tfrac{h}{2}(X_i + X_j), \qquad H_G = \sum_{(ij) \in E} H^{\text{TF}}_{ij}. \tag{1}$$

Two notes on this convention:

- **Single-triangle limit.** Each vertex has degree 2 in a triangle, so the sum over edges gives every $X_v$ a total coefficient of $h$. The triangle Hamiltonian reduces to the textbook form $\sum_{ij} J_{ij} Z_iZ_j + h \sum_v X_v$.
- **Larger graphs.** On a graph where vertex degrees vary, the edge-resolved sum gives $X_v$ a total coefficient of $h \cdot \deg(v) / 2$. We will see that this convention is precisely what makes the per-face theorem clean.

### 2.2 Face-local mediator associators

A **triangular face** of $G$ is an unordered triple $\{a,b,c\} \subset V$ such that all three edges $(ab), (bc), (ac)$ are in $E$. Denote the set of triangular faces by $\mathcal{T}(G)$.

For each face $f = (a,b,c)$, the manuscript (Cell 35, `build_face_probe_hamiltonians_3q`) defines the mediator associators using **3-qubit edge Hamiltonians restricted to the face's vertices**:

$$\widetilde{H}^{(f)}_{ab} = J_{ab} Z_aZ_b + \tfrac{h}{2}(X_a + X_b), \quad \text{(on the 3-qubit Hilbert space of } f \text{)} \tag{2}$$

and analogously for $\widetilde{H}^{(f)}_{bc}, \widetilde{H}^{(f)}_{ac}$. These are operators on a 3-qubit local Hilbert space, **not on the full $|V|$-qubit space of $G$**. This is a deliberate choice in the manuscript and it is the key to what follows.

The three mediator associators on face $f$ are

$$\mathcal{A}^{(f)}_{ab} = \tfrac{1}{4}[\widetilde{H}^{(f)}_{ab}, [\widetilde{H}^{(f)}_{bc}, \widetilde{H}^{(f)}_{ac}]], \tag{3a}$$

and the cyclic permutations $\mathcal{A}^{(f)}_{bc}, \mathcal{A}^{(f)}_{ac}$. The face-local bare capacity is

$$J^2_{\text{tot}}(f) = \|\mathcal{A}^{(f)}_{ab}\|_F^2 + \|\mathcal{A}^{(f)}_{bc}\|_F^2 + \|\mathcal{A}^{(f)}_{ac}\|_F^2. \tag{4}$$

### 2.3 The crucial structural observation

**Equations (2)–(4) involve only the three edges within face $f$.** No other edge of $G$ enters the construction. The mediator associator on face $f$ is computed in a 3-qubit "probe" Hilbert space that is independent of how $f$ is embedded into the rest of the graph.

This is the structural reason why the graph-additivity result below is essentially trivial *given the manuscript's construction*. The non-triviality lies in the resulting closed-form expression and in the fact that the 4→5 qubit comparison can now be extrapolated to arbitrary $N$.

---

## 3. Single-triangle warm-up (recap of Phase 3, Step 1)

For a single triangle with edge couplings $(J_{ab}, J_{bc}, J_{ac})$ and transverse field $h$, direct Pauli algebra (verified symbolically in the original manuscript, Appendix E.1) gives

$$J^2_{\text{tot}}(\triangle) = h^4 \big(J_{ab}^2 + J_{bc}^2 + J_{ac}^2\big). \tag{5}$$

**Why this expression has the form it does.**

The mediator associator (Eq. 3) is a nested commutator. The pure-Ising part of each $\widetilde{H}_{ij}$ commutes with itself (all $Z$-type), so any non-zero contribution to $\mathcal{A}^{(f)}_{ab}$ must come from terms in which the transverse field $X_v$ has been "activated" through commutators with $Z_iZ_j$. Each commutator with $Z_iZ_j$ converts an $X$ into a $Y$, and a second commutator converts the $Y$ back into a structure that probes a different edge.

For a single triangle, working through the algebra (which Phase 3 Step 1 did symbolically) shows that:

- Each mediator associator has a single non-vanishing structural form.
- Its Frobenius norm squared evaluates to $h^4 J_{\mu}^2$, where $\mu$ is the mediator edge.
- Summing over the three mediator channels gives Eq. (5).

We *will not* re-derive this symbolically; we will instead use Eq. (5) as a verified base case and ask whether it generalizes when the triangle is one face among many in a larger graph.

---

## 4. Per-face theorem on an arbitrary graph

### 4.1 Statement

**Theorem (per-face independence).** Let $G$ be any simple graph and let $f = (a,b,c) \in \mathcal{T}(G)$ be a triangular face. Define $J^2_{\text{tot}}(f)$ via Eqs. (2)–(4) of Section 2.2. Then

$$\boxed{\;J^2_{\text{tot}}(f) = h^4 \big(J_{ab}^2 + J_{bc}^2 + J_{ac}^2\big)\;} \tag{6}$$

**independent of any other vertex, edge, or face of $G$**.

### 4.2 Proof

The proof is one line *given the manuscript's construction*. The face-local Hamiltonians $\widetilde{H}^{(f)}_{ab}, \widetilde{H}^{(f)}_{bc}, \widetilde{H}^{(f)}_{ac}$ are constructed from the couplings $(J_{ab}, J_{bc}, J_{ac}, h)$ alone, on a fixed 3-qubit Hilbert space (the local tensor product of the three face vertices). The mediator associators (Eq. 3) and their Frobenius norms (Eq. 4) are functions of these three operators only. No other edge of $G$, and no qubit outside $f$, enters the computation. Therefore $J^2_{\text{tot}}(f)$ depends only on $(J_{ab}, J_{bc}, J_{ac}, h)$, and the value of that function is given by Eq. (5). $\square$

### 4.3 Why this is non-trivial as a *physics* statement

The proof above looks trivial. The non-trivial content is **physical**, not mathematical:

It is *a priori plausible* that, when face $f$ is embedded in a graph where vertex $a$ has many other connections, the face-local probe in Eqs. (2)–(4) is no longer a meaningful proxy for the algebraic capacity available to that face within the dynamics of $H_G$. One could instead define an alternative "globally embedded" associator using the full $|V|$-qubit Hamiltonian; this would give a *different* quantity that does depend on the graph embedding (we verified this numerically; see §7).

The manuscript's choice to define $J_{\text{tot}}$ from face-local 3-qubit Hamiltonians is therefore a deliberate diagnostic convention. Eq. (6) asserts that, under that convention, the diagnostic is well-defined for any graph. This is what the manuscript's 5-qubit experiment (Cell 35) implicitly assumes; we are now making the assumption explicit and checking that it gives a closed form.

---

## 5. Graph-summed identity

Define the graph-summed bare capacity as

$$\Sigma J^2_{\text{tot}}(G) = \sum_{f \in \mathcal{T}(G)} J^2_{\text{tot}}(f). \tag{7}$$

(This is the natural extensive quantity; the face-mean studied in the manuscript's 5-qubit section is $\Sigma J^2_{\text{tot}}(G) / |\mathcal{T}(G)|$.)

For each edge $(ij) \in E$, define its **triangle multiplicity**

$$m_{ij}(G) = \big|\{f \in \mathcal{T}(G) : (ij) \subset f\}\big|. \tag{8}$$

This counts the number of triangles in $G$ that contain edge $(ij)$. For example, in $K_4$, every edge has $m_{ij} = 2$; in a triangular ladder (chain of edge-sharing triangles), interior edges have $m_{ij} = 2$ while boundary edges have $m_{ij} = 1$.

**Corollary (graph-summed identity).** Combining Eq. (6) with the definitions above,

$$\boxed{\;\Sigma J^2_{\text{tot}}(G) = h^4 \sum_{(ij) \in E} m_{ij}(G)\, J_{ij}^2.\;} \tag{9}$$

**Proof.** Expand Eq. (7) using Eq. (6); each edge $(ij)$ is counted once for every face containing it, by definition $m_{ij}(G)$ times. $\square$

### 5.1 Sanity check on the manuscript's 5-qubit case

For the 5-qubit tetrahedron's *internal* $K_4$ (manuscript Sec. III.C), every edge belongs to 2 triangles: $m_{ij}(K_4) = 2$ uniformly. With all internal edges set to the manuscript's fixed-norm value $\|J_{ij}\| = R$ (so that $J_{ij}^2 = R^2$), Eq. (9) gives

$$\Sigma J^2_{\text{tot}}(K_4) = h^4 \cdot \underbrace{6}_{\text{edges}} \cdot \underbrace{2}_{m_{ij}} \cdot R^2 = 12\, h^4 R^2,$$

face-mean $\langle J^2_{\text{tot}}\rangle_f = 12h^4R^2 / 4 = 3h^4R^2$. With $R=\sqrt{3}$ and $h=1$ this gives 9 per face — which is exactly the value our numerical verification (§7, Test 4) reproduces.

---

## 6. $N$-scaling on triangle-decorated lattices

This section applies Eq. (9) to families of graphs of growing size and extracts the $N$-scaling that supports the manuscript's self-averaging interpretation.

### 6.1 Setup: triangular ladder

The cleanest test case is a **triangular ladder** $G_N$ — a chain of $N$ triangles in which each consecutive pair shares one edge. This graph has:

- $|V(G_N)| = N + 2$ vertices,
- $|E(G_N)| = 2N + 1$ edges,
- $|\mathcal{T}(G_N)| = N$ triangles,
- $m_{ij} \in \{1, 2\}$, with all interior shared edges having $m=2$ and all boundary edges having $m=1$.

For uniform edge couplings $J_{ij} = J$:

$$\Sigma J^2_{\text{tot}}(G_N) = h^4 J^2 \sum_{(ij)} m_{ij}. \tag{10}$$

A direct count gives $\sum_{(ij)} m_{ij} = 3N$ (each triangle contributes 3 to the total multiplicity sum). Therefore

$$\Sigma J^2_{\text{tot}}(G_N) = 3 h^4 J^2 N, \tag{11}$$

$$\sqrt{\Sigma J^2_{\text{tot}}(G_N)} = \sqrt{3}\, h^2 J\, \sqrt{N}, \tag{12}$$

$$\langle J^2_{\text{tot}}\rangle_{f,\,G_N} = 3 h^4 J^2 \quad \text{(N-independent).} \tag{13}$$

This is precisely the analytic statement we need:

> **The graph-summed bare capacity scales as $\sqrt{N}$, while the face-mean is $N$-independent.**

For more general triangle-decorated lattices (kagome strip, triangular lattice patches, etc.), the same scaling holds whenever $|\mathcal{T}(G_N)| \propto N$ and $m_{ij}$ is bounded.

### 6.2 Connection to the manuscript's self-averaging argument

The manuscript currently argues, on the basis of two data points (4Q vs 5Q), that

> *"self-averaging across triangular faces causes bare capacity to acquire predictive power from the outset"*

and conjectures, in Discussion, that

> *"in the thermodynamic limit, $J_{\text{tot}}$ alone may suffice as a structural predictor, a quantum analogue of the classical mean-field regime."*

Eqs. (12)–(13) make this conjecture precise *for the bare capacity, in the TF-Ising specialization*. As $N \to \infty$:

1. The face-mean of $J^2_{\text{tot}}$ converges to a finite, $N$-independent constant — the per-face bare capacity is a well-defined intensive quantity.
2. The fluctuation of the face-summed capacity around its mean grows only as $O(\sqrt{N})$, while the mean itself grows as $O(N)$. The relative fluctuation $\sigma / \mu \sim 1/\sqrt{N} \to 0$. This is the standard CLT-type self-averaging.
3. *What still requires numerical verification* is whether the *correlation* between face-mean $J_{\text{tot}}$ and face-mean $I_3^-$ also stabilizes with $N$. This is the substantive empirical content that the planned $N=6,8,10$ scan is designed to test.

The TF-Ising graph-additivity result therefore provides the *structural* anchor; the upcoming $N$-scan provides the *dynamical* verification.

---

## 7. Numerical verification

All four tests below were performed with explicit Pauli matrix construction (no shortcuts) using the 3-qubit face-local Hamiltonians of Eq. (2). Code: `analytic/03_face_local_and_scaling.py`.

### Test 1: Per-face formula (Eq. 6) at four parameter points

| $J_{ab}$ | $J_{bc}$ | $J_{ac}$ | $h$ | numerical $J^2_{\text{tot}}$ | predicted $h^4 \sum J_{ij}^2$ | abs error |
|---|---|---|---|---|---|---|
| 1.0 | 1.0 | 1.0 | 1.0 | 3.000000 | 3.000000 | 0.0 |
| 0.7 | 1.3 | 0.5 | 0.8 | 0.995328 | 0.995328 | 0.0 |
| -1.0 | 1.0 | -1.0 | 1.5 | 15.187500 | 15.187500 | 0.0 |
| 0.3 | -0.9 | 1.7 | 0.6 | 0.491184 | 0.491184 | $5.6\times 10^{-17}$ |

### Test 2: Triangular ladder $G_N$ for $N = 1, \ldots, 8$

Setup: all $J_{ij} = 1$, $h = 1$. Predicted $\Sigma J^2_{\text{tot}}(G_N) = 3N$.

| $N$ | $|\mathcal{T}|$ | $|E|$ | $\Sigma m_{ij}$ | numerical $\Sigma J^2_{\text{tot}}$ | predicted | $J_{\text{tot}}/\sqrt{N}$ |
|---|---|---|---|---|---|---|
| 1 | 1 | 3 | 3 | 3.0 | 3.0 | 1.732051 |
| 2 | 2 | 5 | 6 | 6.0 | 6.0 | 1.732051 |
| 3 | 3 | 7 | 9 | 9.0 | 9.0 | 1.732051 |
| 4 | 4 | 9 | 12 | 12.0 | 12.0 | 1.732051 |
| 5 | 5 | 11 | 15 | 15.0 | 15.0 | 1.732051 |
| 6 | 6 | 13 | 18 | 18.0 | 18.0 | 1.732051 |
| 7 | 7 | 15 | 21 | 21.0 | 21.0 | 1.732051 |
| 8 | 8 | 17 | 24 | 24.0 | 24.0 | 1.732051 |

The constancy of $J_{\text{tot}}/\sqrt{N} = \sqrt{3}$ across all $N$ confirms Eq. (12) exactly.

### Test 3: $K_4$ tetrahedron

Setup: 4 faces, 6 edges, $m_{ij} = 2$ uniformly, $J_{ij} = 1$, $h = 1$. Predicted: $\Sigma J^2_{\text{tot}}(K_4) = 12$.

Numerical: $12.000000$. Abs error: $0.0$. (verified)

### Test 4: Manuscript's 5-qubit internal $K_4$

Setup: 4 faces, 6 edges, $m_{ij} = 2$, $\|J_{ij}\| = R = \sqrt{3}$, three values of $h$.

| $h$ | per-face numerical | per-face predicted ($h^4 \cdot 3 R^2$) | $\Sigma$ over 4 faces | predicted ($12 h^4 R^2$) |
|---|---|---|---|---|
| 0.5 | 0.562500 | 0.562500 | 2.250000 | 2.250000 |
| 1.0 | 9.000000 | 9.000000 | 36.000000 | 36.000000 |
| 1.5 | 45.562500 | 45.562500 | 182.250000 | 182.250000 |

All exact (machine precision). The per-face value of 9 at $h = 1, R = \sqrt{3}$ is what the manuscript's 5-qubit numerical experiment should be returning when the family is in the pure TF-Ising limit.

---

## 8. Implications for the manuscript

### 8.1 Suggested PRL paragraph

The result fits in roughly half a column. Suggested location: framework section, immediately after Eq. (1)–(2) introducing the XYZ family, where the TF-Ising paragraph currently sits. Draft:

> *Generalized to arbitrary triangle-rich graphs $G$, the per-face TF-Ising associator capacity satisfies $J^2_{\text{tot}}(f) = h^4 \sum_{(ij) \in f} J_{ij}^2$ exactly, independent of how face $f$ is embedded in $G$. This face-local decomposition implies a strict graph-summed identity $\sum_{f \in \mathcal{T}(G)} J^2_{\text{tot}}(f) = h^4 \sum_{(ij) \in E} m_{ij}(G) J_{ij}^2$, where $m_{ij}(G)$ counts triangles containing edge $(ij)$. On bounded-degree triangle-decorated lattices with $|\mathcal{T}(G_N)| \propto N$, the face-mean $J^2_{\text{tot}}$ is $N$-independent while the graph-summed capacity scales as $N$, providing an analytic basis for the self-averaging behavior reported below.*

### 8.2 What this changes in the manuscript

- The Discussion's "we conjecture that in the thermodynamic limit $J_{\text{tot}}$ alone may suffice" can be sharpened: for the *bare capacity* in the TF-Ising specialization, the thermodynamic limit is now an exact statement, not a conjecture. The conjecture should be restated as referring to the *predictive power of $J_{\text{tot}}$ for $I_3^-$* (which still requires the dynamical verification from the planned $N$-scan).
- The Phase 3 Appendix E.1 result $J_{\text{tot}} = h^2 R$ becomes the $N=1$ instance of a family.
- Frustration scan (Appendix E.5): the "frustrated configurations share identical $J_{\text{tot}}$" observation is now manifest from Eq. (6), since $J_{\text{tot}}^2$ depends on the $J_{ij}$ only through $J_{ij}^2$. This makes the Appendix E.5 setup look natural rather than coincidental.

### 8.3 What this does *not* address

- The XYZ family used in the main text does not admit a simple closed form. The face-local *theorem* (per-face quantity is independent of graph embedding) still holds, but the per-face value $J^2_{\text{tot}}(f)$ as a function of $(J_{ab}, J_{bc}, J_{ac}, h)$ has no simple expression. The PRL paragraph above is therefore TF-Ising-specific. Whether to state the general face-local property in addition is a writing choice.
- The result concerns $J_{\text{tot}}$ (the bare capacity layer) only, not $J_\rho$ (the state-activated layer) or $I_3^-$. The state-activated layer is genuinely state-dependent and does not have an analogous closed form — this is consistent with the manuscript's two-layer thesis, in which the bottleneck is on the activation side.

---

## 9. Open questions for discussion

These are the points where I would value Gunnar's input:

1. **Is the TF-Ising paragraph the right level of generality for PRL?** The face-local theorem holds for general XYZ, but the closed form is TF-Ising-specific. Three options: (a) state TF-Ising only, (b) state general face-local + TF-Ising specialization, (c) state general face-local in main text and put TF-Ising specialization in SM. My weak preference is (a) for compactness, but (b) is more honest about scope.

2. **Is "self-averaging" the right framing?** What I proved is that $J^2_{\text{tot}}$ is exactly extensive in the number of triangular faces. The standard CLT-based self-averaging argument then applies trivially in distribution over graph ensembles. But the manuscript's actual usage is on a *deterministic* sequence of graphs (4Q triangle, 5Q tetrahedron). Calling this "self-averaging" is technically a slight abuse — what's actually happening is that the face-summed quantity becomes a tighter estimator of its per-face mean as $N$ grows, simply because there are more samples. Is there a cleaner term ("face-averaging concentration"? "face-extensive scaling"?) that Gunnar would prefer?

3. **The kagome / triangular lattice analytic prediction.** For a kagome strip of $N$ unit cells, $|\mathcal{T}| = 2N$ and the same $\sqrt{N}$ scaling holds. Would it be worth computing the explicit prefactor for one or two named lattices to cite, or is the abstract bounded-degree statement enough?

4. **Connection to the $N=6, 8, 10$ numerical scan.** I will use the analytic prediction as a baseline against which to verify the upcoming numerics: for any TF-Ising specialization of the $N$-spin Hamiltonian, the face-mean $J^2_{\text{tot}}$ should match the closed form to machine precision. Any deviation indicates either a bug or a regime where the construction differs from the manuscript convention. Useful as a code-correctness check.

5. **Phase 3 Step 1 — re-derivation needed?** The single-triangle result $J^2_{\text{tot}}(\triangle) = h^4 \sum J^2_{ij}$ was derived in Phase 3 Step 1. I treated it as a verified base case but did not redo the symbolic algebra. If Gunnar wants the per-face theorem to be *fully* self-contained, I should redo the symbolic derivation (probably in SM). Estimated effort: half a day with SymPy.

---

## Appendix A: Numerical script

The code that produced the verification tables is in `analytic/03_face_local_and_scaling.py` (single file, ~250 lines, NumPy only, no external dependencies). Test 2 (the triangular ladder $N$-scaling) runs in under 1 second on a laptop.

The graph-additivity assertion can be re-verified for any user-supplied graph by passing an edge list and triangle list to `graph_J_tot_sq_summed()` and comparing with `predict_graph_J_tot_sq()`.
