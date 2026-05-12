"""
qhoi_unified.py
====================================================================

Unified, performance-optimized infrastructure for quantum
higher-order information (HOI) experiments.

System family this module covers
--------------------------------
N qubits arranged as:
  - (N-1) internal nodes: 0, 1, ..., N-2 forming K_{N-1}
  - 1 environment node: N-1 attached to the last internal node (N-2)

Resulting structure:
  - Internal edges: C(N-1, 2)
  - Environment edges: 1
  - Triangular faces: C(N-1, 3)

Examples:
  N=4 -> K_3 internal + env, 1 face   (triangle network)
  N=5 -> K_4 internal + env, 4 faces  (tetrahedral extension)
  N=6 -> K_5 internal + env, 10 faces
  N=7 -> K_6 internal + env, 20 faces
  N=8 -> K_7 internal + env, 35 faces

Optimizations vs naive implementation
-------------------------------------
Layer 1 - Pauli string cache:
    Pre-build X_i X_j, Y_i Y_j, Z_i Z_j matrices once per N.

Layer 2 - Eigendecomposition cache:
    For each (theta, phi), diagonalize H_tot once. Time evolution at any
    t becomes O(d^2) instead of repeated O(d^3) matrix exponentials.

Layer 3 - Face probe cache:
    For each (theta, phi), pre-compute 3-qubit probe Hamiltonians and
    Jordan associator matrices for every face. J_tot is state-independent.

Edge-shape families
-------------------
Default main-experiment family:
    'broad_xyz' (alias: 'mild_xyz')
    Independent continuous shape components:
        s_ij^mu ~ U(0.60, 1.40),  mu in {x, y, z}.
    This preserves fixed-norm XYZ heterogeneity without imposing a
    discrete Pauli-axis coloring constraint.

Stress-test family:
    'rainbow_anisotropic'
    Strong dominant-axis coloring designed to probe rainbow/null-proximal
    conversion failure. This is not the default and should not be mixed
    with the main finite-size ensemble.

Health check
------------
The structural health check is optional and disabled by default. Main
finite-size ensemble runs should usually use a fixed seed panel without
outcome-based seed search. If enabled, the check screens only the
Hamiltonian-level min face J_tot at a reference point; it never uses I3^-
or correlation outcomes.

Quick usage
-----------
    import numpy as np
    import pandas as pd
    import qhoi_unified as qhoi

    theta_vals = np.linspace(0.1, np.pi/2, 12)
    phi_vals   = np.linspace(0.0, np.pi/2, 12)
    t_vals     = np.linspace(0.05, 1.0, 50)

    per_t, summaries = qhoi.scan_time_grid(
        theta_vals, phi_vals, t_vals,
        N=4, base_seed=1, state_seed=456,
        family='broad_xyz', health_check=False,
    )
    df = pd.DataFrame(summaries)

CLI usage
---------
    python qhoi_unified.py --N 7 --n_states 50 \
        --shape_family broad_xyz --output_dir results_N7_broad_xyz
"""

import os
import argparse
import warnings
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np


# ====================================================================
# 1. Pauli matrices and basic utilities
# ====================================================================
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)


def kron_all(op_list):
    """Kronecker product of a list of operators or state vectors."""
    out = np.array([[1.0 + 0j]])
    for op in op_list:
        out = np.kron(out, op)
    return out


def random_single_qubit_state(rng):
    """Random Haar-distributed single-qubit pure state."""
    v = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    return v / np.linalg.norm(v)


# ====================================================================
# 2. Topology generation
# ====================================================================
@dataclass
class Topology:
    N: int
    internal_nodes: List[int]
    env_node: int
    internal_edges: List[Tuple[str, int, int]]
    env_edge: Tuple[str, int, int]
    faces: Dict[str, dict]
    distal_faces: List[str]
    adjacent_faces: List[str]


def _label(i):
    """Map node index 0,1,2,... to letter labels A,B,C,...; wrap after Z."""
    if i < 26:
        return chr(65 + i)
    return f"N{i}"


def build_topology(N):
    """Build K_{N-1} internal core + one attached environment node."""
    if N < 4:
        raise ValueError("N must be >= 4 (need at least one internal triangle).")

    n_int = N - 1
    internal_nodes = list(range(n_int))
    env_node = N - 1
    attach = n_int - 1

    internal_edges = []
    for i, j in combinations(range(n_int), 2):
        name = _label(i) + _label(j)
        internal_edges.append((name, i, j))

    env_edge = (_label(attach) + _label(env_node), attach, env_node)

    faces = {}
    for a, b, c in combinations(range(n_int), 3):
        name = _label(a) + _label(b) + _label(c)
        faces[name] = {
            "nodes_global": [a, b, c],
            "edges_local": (
                _label(a) + _label(b),
                _label(b) + _label(c),
                _label(a) + _label(c),
            ),
        }

    attach_label = _label(attach)
    distal_faces = [name for name in faces if attach_label not in name]
    adjacent_faces = [name for name in faces if attach_label in name]

    return Topology(
        N=N,
        internal_nodes=internal_nodes,
        env_node=env_node,
        internal_edges=internal_edges,
        env_edge=env_edge,
        faces=faces,
        distal_faces=distal_faces,
        adjacent_faces=adjacent_faces,
    )


# ====================================================================
# 3. Edge-shape generation families
# ====================================================================
def edge_shape_broad_xyz(i, j, base_seed=1, s_lo=0.60, s_hi=1.40):
    """Broad continuous XYZ family: independent U(s_lo, s_hi) components.

    This is the default main-experiment generator. It creates edge-to-edge
    heterogeneity without imposing a discrete Pauli-axis coloring constraint.
    """
    if i > j:
        i, j = j, i
    rng = np.random.default_rng(int(base_seed) * 100003 + int(i) * 1009 + int(j))
    return rng.uniform(s_lo, s_hi, size=3)


def edge_shape_rainbow_anisotropic(
    i,
    j,
    base_seed=1,
    sub_lo=0.30,
    sub_hi=0.65,
    dom_lo=0.90,
    dom_hi=1.10,
):
    """Rainbow-favoring anisotropic family for stress testing only.

    Each edge receives one dominant Pauli axis and two subdominant axes.
    The dominant axis is (i + j + base_seed) mod 3, which systematically
    favors rainbow K_3 triangles. Use this only as a null-proximal / failure
    mode stress test, not for the main finite-size ensemble.
    """
    if i > j:
        i, j = j, i
    rng = np.random.default_rng(int(base_seed) * 100003 + int(i) * 1009 + int(j))
    dom_axis = (int(i) + int(j) + int(base_seed)) % 3
    shape = rng.uniform(sub_lo, sub_hi, size=3)
    shape[dom_axis] = rng.uniform(dom_lo, dom_hi)
    return shape


SHAPE_FAMILIES = {
    "broad_xyz": edge_shape_broad_xyz,
    "mild_xyz": edge_shape_broad_xyz,  # backward-compatible alias
    "rainbow_anisotropic": edge_shape_rainbow_anisotropic,
}
DEFAULT_SHAPE_FAMILY = "broad_xyz"


def edge_shape(i, j, base_seed=1, family=DEFAULT_SHAPE_FAMILY):
    """Dispatch to an edge-shape generator."""
    if family not in SHAPE_FAMILIES:
        raise ValueError(
            f"Unknown shape family: {family!r}. "
            f"Available: {list(SHAPE_FAMILIES.keys())}"
        )
    return SHAPE_FAMILIES[family](i, j, base_seed=base_seed)


def normalize_edge_coeffs(theta, phi, shape, R=np.sqrt(3.0), eps=1e-12):
    """Map angular seed and shape vector to fixed-norm XYZ coefficients."""
    direction = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ],
        dtype=float,
    )
    vec = direction * shape
    norm = np.linalg.norm(vec)
    if norm < eps:
        raise ValueError(f"Degenerate edge coeffs at theta={theta}, phi={phi}.")
    return R * vec / norm


def all_edge_coeffs(
    topology,
    theta,
    phi,
    R=np.sqrt(3.0),
    gE=0.20,
    base_seed=1,
    family=DEFAULT_SHAPE_FAMILY,
):
    """Build dict {edge_name: (Jx, Jy, Jz)} for all internal + env edges."""
    coeffs = {}
    for name, i, j in topology.internal_edges:
        shape = edge_shape(i, j, base_seed=base_seed, family=family)
        coeffs[name] = normalize_edge_coeffs(theta, phi, shape, R=R)

    name_e, i_e, j_e = topology.env_edge
    shape_e = edge_shape(i_e, j_e, base_seed=base_seed, family=family)
    coeffs[name_e] = normalize_edge_coeffs(theta, phi, shape_e, R=gE * R)
    return coeffs


# ====================================================================
# 4. Full-space Pauli string cache
# ====================================================================
class PauliCache:
    """Pre-built X_iX_j, Y_iY_j, Z_iZ_j matrices for every edge."""

    def __init__(self, topology):
        self.N = topology.N
        self.cache = {}
        for name, i, j in list(topology.internal_edges) + [topology.env_edge]:
            self.cache[name] = (
                self._embed(i, j, sx, sx),
                self._embed(i, j, sy, sy),
                self._embed(i, j, sz, sz),
            )

    def _embed(self, i, j, Oi, Oj):
        ops = [I2] * self.N
        ops[i] = Oi
        ops[j] = Oj
        return kron_all(ops)

    def assemble_H(self, coeffs):
        d = 2 ** self.N
        H = np.zeros((d, d), dtype=complex)
        for name, (XX, YY, ZZ) in self.cache.items():
            Jx, Jy, Jz = coeffs[name]
            H += Jx * XX + Jy * YY + Jz * ZZ
        return H


# ====================================================================
# 5. Eigendecomposition-based time evolution
# ====================================================================
class TimeEvolution:
    """Cache eigendecomposition of H; compute psi(t) cheaply for many t."""

    def __init__(self, H, psi0):
        self.E, self.V = np.linalg.eigh(H)
        self.c0 = self.V.conj().T @ psi0

    def psi_at(self, t):
        phases = np.exp(-1j * self.E * t)
        return self.V @ (phases * self.c0)


# ====================================================================
# 6. Partial trace and face probe cache
# ====================================================================
def pure_state_partial_trace(psi, keep, dims):
    """Reduced density matrix from pure state via reshape + matmul."""
    n = len(dims)
    keep = list(keep)
    trace_out = [i for i in range(n) if i not in keep]
    perm = keep + trace_out
    psi_r = psi.reshape(dims).transpose(perm)
    d_keep = int(np.prod([dims[i] for i in keep]))
    d_trace = int(np.prod([dims[i] for i in trace_out]))
    M = psi_r.reshape(d_keep, d_trace)
    return M @ M.conj().T


def partial_trace_dm(rho, keep, dims):
    """Partial trace of a mixed density matrix."""
    n = len(dims)
    keep = list(keep)
    trace_out = [i for i in range(n) if i not in keep]
    rho_t = rho.reshape(dims + dims)
    perm = keep + trace_out + [i + n for i in keep] + [i + n for i in trace_out]
    rho_t = rho_t.transpose(perm)
    d_keep = int(np.prod([dims[i] for i in keep]))
    d_trace = int(np.prod([dims[i] for i in trace_out]))
    rho_t = rho_t.reshape(d_keep, d_trace, d_keep, d_trace)
    return np.trace(rho_t, axis1=1, axis2=3)


def jordan_associator(X, Y, Z):
    """(X o Y) o Z - X o (Y o Z), with X o Y=(XY+YX)/2."""
    XY = 0.5 * (X @ Y + Y @ X)
    YZ = 0.5 * (Y @ Z + Z @ Y)
    return 0.5 * (XY @ Z + Z @ XY) - 0.5 * (X @ YZ + YZ @ X)


def _embed_3q(i, j, J):
    """Build Jx X_iX_j + Jy Y_iY_j + Jz Z_iZ_j on a 3-qubit space."""
    Jx, Jy, Jz = J
    out = np.zeros((8, 8), dtype=complex)
    for coeff, P in [(Jx, sx), (Jy, sy), (Jz, sz)]:
        ops = [I2, I2, I2]
        ops[i] = P
        ops[j] = P
        out += coeff * kron_all(ops)
    return out


class FaceProbeCache:
    """Cache the three Jordan associators and J_tot for each face."""

    def __init__(self, coeffs, topology):
        self.cache = {}
        for face_name, info in topology.faces.items():
            e01, e12, e02 = info["edges_local"]

            H_e01 = _embed_3q(0, 1, coeffs[e01])
            H_e12 = _embed_3q(1, 2, coeffs[e12])
            H_e02 = _embed_3q(0, 2, coeffs[e02])

            A_e01 = jordan_associator(H_e12, H_e01, H_e02)
            A_e12 = jordan_associator(H_e01, H_e12, H_e02)
            A_e02 = jordan_associator(H_e01, H_e02, H_e12)

            j_e01 = np.linalg.norm(A_e01, "fro")
            j_e12 = np.linalg.norm(A_e12, "fro")
            j_e02 = np.linalg.norm(A_e02, "fro")
            J_tot = float(np.sqrt(j_e01**2 + j_e12**2 + j_e02**2))

            self.cache[face_name] = {
                "A_e01": A_e01,
                "A_e12": A_e12,
                "A_e02": A_e02,
                "j_e01": float(j_e01),
                "j_e12": float(j_e12),
                "j_e02": float(j_e02),
                "J_tot": J_tot,
            }


# ====================================================================
# 7. Information observables
# ====================================================================
def von_neumann_entropy(rho, base=2, eps=1e-12):
    eigs = np.linalg.eigvalsh(rho)
    eigs = eigs[eigs > eps]
    return float(-np.sum(eigs * np.log(eigs)) / np.log(base))


def I3_of_3q(rho):
    rho_A = partial_trace_dm(rho, [0], [2, 2, 2])
    rho_B = partial_trace_dm(rho, [1], [2, 2, 2])
    rho_C = partial_trace_dm(rho, [2], [2, 2, 2])
    rho_AB = partial_trace_dm(rho, [0, 1], [2, 2, 2])
    rho_AC = partial_trace_dm(rho, [0, 2], [2, 2, 2])
    rho_BC = partial_trace_dm(rho, [1, 2], [2, 2, 2])
    return (
        von_neumann_entropy(rho_A)
        + von_neumann_entropy(rho_B)
        + von_neumann_entropy(rho_C)
        - von_neumann_entropy(rho_AB)
        - von_neumann_entropy(rho_AC)
        - von_neumann_entropy(rho_BC)
        + von_neumann_entropy(rho)
    )


def mutual_info_2q(rho2):
    rho_a = partial_trace_dm(rho2, [0], [2, 2])
    rho_b = partial_trace_dm(rho2, [1], [2, 2])
    return von_neumann_entropy(rho_a) + von_neumann_entropy(rho_b) - von_neumann_entropy(rho2)


def negativity(rho, dims=(2, 2), pt_subsys=0):
    """Negativity = sum abs(negative eigenvalues) of partial transpose."""
    n = len(dims)
    rho_t = rho.reshape(list(dims) + list(dims))
    perm = list(range(2 * n))
    perm[pt_subsys], perm[pt_subsys + n] = perm[pt_subsys + n], perm[pt_subsys]
    rho_t = rho_t.transpose(perm).reshape(int(np.prod(dims)), int(np.prod(dims)))
    eigs = np.linalg.eigvalsh(rho_t)
    return float(np.sum(np.abs(eigs[eigs < 0])))


def pair_negativity_sum_3q(rho3):
    rAB = partial_trace_dm(rho3, [0, 1], [2, 2, 2])
    rAC = partial_trace_dm(rho3, [0, 2], [2, 2, 2])
    rBC = partial_trace_dm(rho3, [1, 2], [2, 2, 2])
    return negativity(rAB) + negativity(rAC) + negativity(rBC)


def pair_mi_sum_3q(rho3):
    rAB = partial_trace_dm(rho3, [0, 1], [2, 2, 2])
    rAC = partial_trace_dm(rho3, [0, 2], [2, 2, 2])
    rBC = partial_trace_dm(rho3, [1, 2], [2, 2, 2])
    return mutual_info_2q(rAB) + mutual_info_2q(rAC) + mutual_info_2q(rBC)


def linear_entropy(rho):
    return 1.0 - float(np.real(np.trace(rho @ rho)))


def safe_corr(x, y, eps=1e-12):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or np.std(x) < eps or np.std(y) < eps:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


# ====================================================================
# 8. Initial state generation
# ====================================================================
def make_random_product_state(N, seed):
    """Random product state: (N-1) random internal qubits tensor |0>_env."""
    rng = np.random.default_rng(int(seed))
    factors = [random_single_qubit_state(rng) for _ in range(N - 1)]
    factors.append(ket0)
    return kron_all(factors).reshape(-1)


# ====================================================================
# 9. Per-face evaluation
# ====================================================================
def evaluate_face_at_t(psi_t, face_info, probe, dims):
    keep = face_info["nodes_global"]
    rho_face = pure_state_partial_trace(psi_t, keep, dims)

    I3 = I3_of_3q(rho_face)
    I3m = max(0.0, -I3)

    m01 = complex(np.trace(rho_face @ probe["A_e01"]))
    m12 = complex(np.trace(rho_face @ probe["A_e12"]))
    m02 = complex(np.trace(rho_face @ probe["A_e02"]))
    J_rho = float(np.sqrt(abs(m01) ** 2 + abs(m12) ** 2 + abs(m02) ** 2))

    J_tot = probe["J_tot"]
    eta = J_rho / J_tot if J_tot > 1e-12 else float("nan")

    return {
        "I3": float(I3),
        "I3_minus": float(I3m),
        "J_rho": J_rho,
        "J_tot": J_tot,
        "eta": float(eta),
        "SL": linear_entropy(rho_face),
        "P_pair": pair_mi_sum_3q(rho_face),
        "N_pair_sum": pair_negativity_sum_3q(rho_face),
    }


# ====================================================================
# 10. Optional structural health check
# ====================================================================
HEALTH_REFERENCE_THETA = np.pi / 4
HEALTH_REFERENCE_PHI = np.pi / 4
DEFAULT_HEALTH_THRESHOLD = 18.0


def family_min_face_jtot(
    topology,
    base_seed,
    R=np.sqrt(3.0),
    gE=0.20,
    theta_ref=HEALTH_REFERENCE_THETA,
    phi_ref=HEALTH_REFERENCE_PHI,
    family=DEFAULT_SHAPE_FAMILY,
):
    """Min face J_tot at a reference point; optional structural diagnostic."""
    try:
        coeffs = all_edge_coeffs(
            topology,
            theta_ref,
            phi_ref,
            R=R,
            gE=gE,
            base_seed=base_seed,
            family=family,
        )
    except ValueError:
        return 0.0
    probes = FaceProbeCache(coeffs, topology)
    jtots = [probes.cache[fn]["J_tot"] for fn in topology.faces]
    return float(min(jtots)) if jtots else float("inf")


def find_healthy_base_seed(
    topology,
    base_seed_start=1,
    threshold=DEFAULT_HEALTH_THRESHOLD,
    max_attempts=100,
    R=np.sqrt(3.0),
    gE=0.20,
    family=DEFAULT_SHAPE_FAMILY,
    verbose=False,
):
    """Increment base_seed until the optional structural health threshold passes."""
    for attempt in range(max_attempts):
        seed = int(base_seed_start) + attempt
        jt = family_min_face_jtot(topology, seed, R=R, gE=gE, family=family)
        if verbose and attempt > 0:
            print(f"      attempt {attempt}: seed={seed}, min J_tot={jt:.3f}")
        if jt >= threshold:
            return seed, jt
    raise RuntimeError(
        f"No healthy family found after {max_attempts} attempts starting from "
        f"base_seed={base_seed_start} (family={family!r}). "
        f"Consider lowering threshold (current={threshold})."
    )


def _resolve_family_seed(
    topology,
    base_seed,
    family,
    R,
    gE,
    health_check,
    health_threshold,
    health_max_attempts,
    verbose=False,
):
    """Internal helper for seed handling and metadata."""
    if family == "rainbow_anisotropic" and health_check:
        warnings.warn(
            "health_check=True was requested for rainbow_anisotropic. Disabling "
            "health_check so the null-proximal stress-test family is not filtered.",
            RuntimeWarning,
        )
        health_check = False

    if health_check:
        effective_seed, min_jtot = find_healthy_base_seed(
            topology,
            base_seed_start=base_seed,
            threshold=health_threshold,
            max_attempts=health_max_attempts,
            R=R,
            gE=gE,
            family=family,
            verbose=verbose,
        )
        if verbose and effective_seed != base_seed:
            print(
                f"    Health check: base_seed {base_seed} -> {effective_seed} "
                f"(min face J_tot={min_jtot:.3f}, threshold={health_threshold}, "
                f"family={family!r})"
            )
    else:
        effective_seed = int(base_seed)
        min_jtot = family_min_face_jtot(
            topology, effective_seed, R=R, gE=gE, family=family
        )

    return int(effective_seed), float(min_jtot), bool(health_check)


# ====================================================================
# 11. Main scan: one initial state
# ====================================================================
def scan_time_grid(
    theta_vals,
    phi_vals,
    t_vals,
    N=6,
    R=np.sqrt(3.0),
    gE=0.20,
    base_seed=1,
    state_seed=456,
    psi0=None,
    verbose=False,
    health_threshold=DEFAULT_HEALTH_THRESHOLD,
    health_max_attempts=100,
    health_check=False,
    family=DEFAULT_SHAPE_FAMILY,
):
    """Single-state scan over (theta, phi) grid and time array.

    Notes
    -----
    state_seed is recorded in every summary row as metadata. It is used to
    construct psi0 only when psi0 is None. If a caller passes psi0 directly
    (e.g., from scan_ensemble), state_seed should be set to a value that
    identifies that psi0; otherwise the recorded metadata will not match
    the actual initial state.
    """
    topology = build_topology(N)
    effective_seed, min_jtot, health_check = _resolve_family_seed(
        topology,
        base_seed=base_seed,
        family=family,
        R=R,
        gE=gE,
        health_check=health_check,
        health_threshold=health_threshold,
        health_max_attempts=health_max_attempts,
        verbose=verbose,
    )

    pauli = PauliCache(topology)
    dims = [2] * N

    if psi0 is None:
        psi0 = make_random_product_state(N, state_seed)

    per_t_grid_data = {float(t): [] for t in t_vals}
    n_total = len(theta_vals) * len(phi_vals)
    counter = 0

    for th in theta_vals:
        for ph in phi_vals:
            counter += 1
            if verbose and (counter % max(1, n_total // 10) == 0):
                print(f"    [{counter}/{n_total}] theta={th:.3f}, phi={ph:.3f}")

            try:
                coeffs = all_edge_coeffs(
                    topology,
                    th,
                    ph,
                    R=R,
                    gE=gE,
                    base_seed=effective_seed,
                    family=family,
                )
            except ValueError:
                continue

            H_tot = pauli.assemble_H(coeffs)
            evo = TimeEvolution(H_tot, psi0)
            probes = FaceProbeCache(coeffs, topology)

            for t in t_vals:
                psi_t = evo.psi_at(t)
                face_data = {
                    fn: evaluate_face_at_t(psi_t, fi, probes.cache[fn], dims)
                    for fn, fi in topology.faces.items()
                }

                I3m_v = np.array([face_data[f]["I3_minus"] for f in face_data])
                Jr_v = np.array([face_data[f]["J_rho"] for f in face_data])
                Jt_v = np.array([face_data[f]["J_tot"] for f in face_data])
                eta_v = np.array([face_data[f]["eta"] for f in face_data])
                SL_v = np.array([face_data[f]["SL"] for f in face_data])
                Pp_v = np.array([face_data[f]["P_pair"] for f in face_data])
                Np_v = np.array([face_data[f]["N_pair_sum"] for f in face_data])

                I3m_dist = np.array(
                    [face_data[f]["I3_minus"] for f in topology.distal_faces]
                )
                I3m_adj = np.array(
                    [face_data[f]["I3_minus"] for f in topology.adjacent_faces]
                )

                env_keep = [topology.env_edge[1], topology.env_edge[2]]
                rho_env = pure_state_partial_trace(psi_t, env_keep, dims)
                N_env = negativity(rho_env)
                I_env = mutual_info_2q(rho_env)

                row = {
                    "theta": float(th),
                    "phi": float(ph),
                    "t": float(t),
                    "N": int(N),
                    "base_seed_requested": int(base_seed),
                    "effective_base_seed": int(effective_seed),
                    "state_seed": int(state_seed),
                    "family": str(family),
                    "family_min_jtot": float(min_jtot),
                    "health_check": bool(health_check),
                    "I3_minus_facemean": float(np.mean(I3m_v)),
                    "I3_minus_facestd": float(np.std(I3m_v)),
                    "J_rho_facemean": float(np.mean(Jr_v)),
                    "J_rho_facestd": float(np.std(Jr_v)),
                    "J_tot_facemean": float(np.mean(Jt_v)),
                    "J_tot_facestd": float(np.std(Jt_v)),
                    "eta_facemean": float(np.nanmean(eta_v)),
                    "SL_facemean": float(np.mean(SL_v)),
                    "P_pair_facemean": float(np.mean(Pp_v)),
                    "N_pair_sum_facemean": float(np.mean(Np_v)),
                    "I3_minus_distal_mean": float(np.mean(I3m_dist))
                    if I3m_dist.size
                    else float("nan"),
                    "I3_minus_adj_mean": float(np.mean(I3m_adj))
                    if I3m_adj.size
                    else float("nan"),
                    "N_env": float(N_env),
                    "I_env": float(I_env),
                    "face_data": face_data,
                }
                per_t_grid_data[float(t)].append(row)

    time_summaries = []
    for t in t_vals:
        rows = per_t_grid_data[float(t)]
        if not rows:
            continue
        target = np.array([r["I3_minus_facemean"] for r in rows])
        Jr = np.array([r["J_rho_facemean"] for r in rows])
        Jt = np.array([r["J_tot_facemean"] for r in rows])
        SL = np.array([r["SL_facemean"] for r in rows])
        Pp = np.array([r["P_pair_facemean"] for r in rows])
        Np = np.array([r["N_pair_sum_facemean"] for r in rows])
        Ne = np.array([r["N_env"] for r in rows])

        time_summaries.append(
            {
                "t": float(t),
                "N": int(N),
                "n_points": len(rows),
                "base_seed_requested": int(base_seed),
                "effective_base_seed": int(effective_seed),
                "state_seed": int(state_seed),
                "family": str(family),
                "family_min_jtot": float(min_jtot),
                "health_check": bool(health_check),
                "corr_J_rho": safe_corr(Jr, target),
                "corr_J_tot": safe_corr(Jt, target),
                "corr_SL": safe_corr(SL, target),
                "corr_P_pair": safe_corr(Pp, target),
                "corr_N_pair": safe_corr(Np, target),
                "corr_N_env": safe_corr(Ne, target),
                "target_max": float(np.nanmax(target)),
                "target_mean": float(np.nanmean(target)),
                "J_tot_min_grid": float(np.nanmin(Jt)),
                "J_tot_mean_grid": float(np.nanmean(Jt)),
            }
        )

    return per_t_grid_data, time_summaries


# ====================================================================
# 12. Ensemble scan: many initial states
# ====================================================================
def scan_ensemble(
    theta_vals,
    phi_vals,
    t_vals,
    N=6,
    n_states=50,
    base_seed=1,
    state_seed_start=1000,
    R=np.sqrt(3.0),
    gE=0.20,
    output_dir=None,
    verbose=False,
    health_threshold=DEFAULT_HEALTH_THRESHOLD,
    health_max_attempts=100,
    health_check=False,
    family=DEFAULT_SHAPE_FAMILY,
):
    """Run an initial-state ensemble for a shared Hamiltonian shape family."""
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    topology = build_topology(N)
    effective_seed, min_jtot, health_check = _resolve_family_seed(
        topology,
        base_seed=base_seed,
        family=family,
        R=R,
        gE=gE,
        health_check=health_check,
        health_threshold=health_threshold,
        health_max_attempts=health_max_attempts,
        verbose=verbose,
    )

    if verbose:
        print(
            f"Ensemble family={family!r}, requested seed={base_seed}, "
            f"effective seed={effective_seed}, min face J_tot={min_jtot:.3f}, "
            f"health_check={health_check}"
        )

    all_summaries = []
    safe_family = str(family).replace("/", "_")

    for k in range(n_states):
        seed_k = state_seed_start + k
        psi0_k = make_random_product_state(N, seed_k)
        if verbose:
            print(f"  state {k + 1}/{n_states} (state_seed={seed_k})")

        _, summaries = scan_time_grid(
            theta_vals,
            phi_vals,
            t_vals,
            N=N,
            R=R,
            gE=gE,
            base_seed=effective_seed,
            state_seed=seed_k,
            psi0=psi0_k,
            verbose=False,
            health_check=False,
            family=family,
        )

        for s in summaries:
            s["base_seed_requested"] = int(base_seed)
            s["effective_base_seed"] = int(effective_seed)
            s["state_seed"] = int(seed_k)
            s["state_idx"] = int(k)
            s["family"] = str(family)
            s["family_min_jtot"] = float(min_jtot)
            s["health_check"] = bool(health_check)
        all_summaries.extend(summaries)

        if output_dir is not None:
            try:
                import pandas as pd

                filename = (
                    f"N{N}_{safe_family}_seed{effective_seed:04d}_state{k:04d}.csv"
                )
                pd.DataFrame(summaries).to_csv(
                    os.path.join(output_dir, filename), index=False
                )
            except ImportError:
                pass

    return all_summaries


# ====================================================================
# 13. Smoke test
# ====================================================================
def smoke_test():
    print("=== qhoi_unified smoke test ===")
    print()

    for N in [4, 5, 6]:
        topo = build_topology(N)
        print(f"N={N}: internal nodes={topo.internal_nodes}, env={topo.env_node}")
        print(f"      internal edges: {len(topo.internal_edges)}, faces: {len(topo.faces)}")
        print(f"      distal: {topo.distal_faces}, adjacent: {topo.adjacent_faces}")
        print()

    per_t, _ = scan_time_grid(
        [0.5], [0.5], [0.25],
        N=4, base_seed=1, state_seed=42,
        family="broad_xyz", health_check=False,
    )
    row = per_t[0.25][0]
    print("--- N=4 single point, broad_xyz ---")
    print(f"  I3_minus (ABC): {row['face_data']['ABC']['I3_minus']:.6f}")
    print(f"  J_rho    (ABC): {row['face_data']['ABC']['J_rho']:.6f}")
    print(f"  J_tot    (ABC): {row['face_data']['ABC']['J_tot']:.6f}")
    print(f"  N_env:          {row['N_env']:.6f}")
    print(f"  I_env:          {row['I_env']:.6f}")

    per_t6, _ = scan_time_grid(
        [0.5], [0.5], [0.25],
        N=6, base_seed=1, state_seed=42,
        family="broad_xyz", health_check=False,
    )
    row6 = per_t6[0.25][0]
    print()
    print("--- N=6 single point, broad_xyz ---")
    print(f"  Faces evaluated: {len(row6['face_data'])}")
    print(f"  Face-mean I3_minus: {row6['I3_minus_facemean']:.6f} ± {row6['I3_minus_facestd']:.6f}")
    print(f"  Face-mean J_rho:    {row6['J_rho_facemean']:.6f}")
    print(f"  Face-mean J_tot:    {row6['J_tot_facemean']:.6f}")
    print(f"  Distal-face I3:     {row6['I3_minus_distal_mean']:.6f}")
    print(f"  Adjacent-face I3:   {row6['I3_minus_adj_mean']:.6f}")

    print()
    print("Smoke test passed.")


# ====================================================================
# 14. CLI driver
# ====================================================================
def main():
    parser = argparse.ArgumentParser(
        description="qhoi_unified: ensemble runs of N-qubit HOI scans."
    )
    parser.add_argument("--smoke", action="store_true", help="Run smoke test.")
    parser.add_argument("--N", type=int, default=6, help="Total number of qubits (>=4).")
    parser.add_argument("--n_states", type=int, default=50, help="Initial-state ensemble size.")
    parser.add_argument("--gE", type=float, default=0.20)
    parser.add_argument("--R", type=float, default=float(np.sqrt(3.0)))
    parser.add_argument("--n_grid", type=int, default=12, help="Grid size in each angle.")
    parser.add_argument("--base_seed", type=int, default=1, help="Hamiltonian shape seed.")
    parser.add_argument("--state_seed_start", type=int, default=1000)
    parser.add_argument("--output_dir", type=str, default=None)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--health_threshold",
        type=float,
        default=DEFAULT_HEALTH_THRESHOLD,
        help=f"Minimum face J_tot at reference point. Default={DEFAULT_HEALTH_THRESHOLD}.",
    )
    parser.add_argument("--health_max_attempts", type=int, default=100)
    parser.add_argument(
        "--health_check",
        action="store_true",
        help="Enable optional structural health check / seed search. Disabled by default.",
    )
    parser.add_argument(
        "--no_health_check",
        action="store_true",
        help="Backward-compatible flag. Health check is already disabled unless --health_check is supplied.",
    )
    parser.add_argument(
        "--shape_family",
        type=str,
        default=DEFAULT_SHAPE_FAMILY,
        choices=list(SHAPE_FAMILIES.keys()),
        help=f"Edge shape family. Default={DEFAULT_SHAPE_FAMILY!r}.",
    )

    args = parser.parse_args()

    if args.smoke:
        smoke_test()
        return

    if args.output_dir is None:
        args.output_dir = f"results_N{args.N}_{args.shape_family}"

    theta_vals = np.linspace(0.1, np.pi / 2, args.n_grid)
    phi_vals = np.linspace(0.0, np.pi / 2, args.n_grid)
    t_vals = np.sort(
        np.unique(
            np.concatenate(
                [
                    np.linspace(0.05, 0.15, 5),
                    np.linspace(0.15, 0.25, 9),
                    np.linspace(0.25, 0.45, 11),
                    np.linspace(0.45, 0.70, 6),
                    np.linspace(0.70, 0.90, 9),
                    np.linspace(0.90, 1.00, 3),
                ]
            )
        )
    )

    health_enabled = bool(args.health_check and not args.no_health_check)
    if args.shape_family == "rainbow_anisotropic" and health_enabled:
        warnings.warn(
            "--health_check requested for rainbow_anisotropic. It will be disabled "
            "inside scan_ensemble to preserve the stress-test family.",
            RuntimeWarning,
        )

    print("=== qhoi_unified ensemble run ===")
    print(f"N = {args.N}  (internal K_{args.N - 1}, {len(list(combinations(range(args.N - 1), 3)))} faces)")
    print(f"n_states  = {args.n_states}")
    print(f"grid      = {args.n_grid} x {args.n_grid}")
    print(f"|t_vals|  = {len(t_vals)}")
    print(f"family    = {args.shape_family!r}")
    print(
        f"health_threshold = {args.health_threshold}"
        + (" (ENABLED)" if health_enabled else " (DISABLED)")
    )
    print(f"output    = {args.output_dir}")
    print()

    scan_ensemble(
        theta_vals,
        phi_vals,
        t_vals,
        N=args.N,
        n_states=args.n_states,
        base_seed=args.base_seed,
        state_seed_start=args.state_seed_start,
        R=args.R,
        gE=args.gE,
        output_dir=args.output_dir,
        verbose=args.verbose,
        health_threshold=args.health_threshold,
        health_max_attempts=args.health_max_attempts,
        health_check=health_enabled,
        family=args.shape_family,
    )

    print("Done.")


if __name__ == "__main__":
    main()
