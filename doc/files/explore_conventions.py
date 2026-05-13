"""
Step 2: Graph-additivity test on multi-triangle graphs.

Two transverse-field conventions:

  (A) Edge-resolved (manuscript Appendix E):
        H = sum_{(ij) in E} [ J_ij Z_i Z_j + (h/2)(X_i + X_j) ]
      This puts coefficient h*deg(v)/2 in front of each X_v.

  (B) Global / textbook:
        H = sum_{(ij) in E} J_ij Z_i Z_j + h sum_{v in V} X_v
      Each X_v has fixed coefficient h, independent of the graph.

On a single triangle every vertex has degree 2, so (A) and (B) coincide
with effective field h. On K_4 (tetrahedron, deg=3 everywhere) they do
not coincide: (A) effectively uses 3h/2.

For each convention we test the conjecture:
  J_tot^2(G) = h_eff^4 * sum_{triangle f in G} sum_{(ij) in f} J_ij^2

where h_eff is the per-vertex X coefficient under the chosen convention.
For convention (A), h_eff varies per face (depending on incident vertex
degrees). For convention (B), h_eff = h uniformly.

This script reports:
- J_tot^2 numerically computed
- Predicted J_tot^2 under both conventions
- Whether one of them matches face-additively
"""

import numpy as np
from itertools import combinations

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for op in ops:
        out = np.kron(out, op)
    return out


def embed(n, sites_ops):
    ops = [I2] * n
    for site, op in sites_ops.items():
        ops[site] = op
    return kron_all(ops)


def commutator(A, B):
    return A @ B - B @ A


def edge_zz(n, i, j, J_ij):
    """J_ij Z_i Z_j on n-qubit space."""
    return J_ij * embed(n, {i: sz, j: sz})


def vertex_x(n, v, h):
    """h X_v on n-qubit space."""
    return h * embed(n, {v: sx})


def build_H_edge_resolved(n_qubits, edges_with_J, h):
    """Convention (A): H = sum_ij [J_ij ZZ + (h/2)(X_i + X_j)]."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for (i, j), J_ij in edges_with_J.items():
        H += edge_zz(n_qubits, i, j, J_ij)
        H += vertex_x(n_qubits, i, h / 2.0)
        H += vertex_x(n_qubits, j, h / 2.0)
    return H


def build_H_global(n_qubits, edges_with_J, h):
    """Convention (B): H = sum_ij J_ij ZZ + h sum_v X_v."""
    H = np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    for (i, j), J_ij in edges_with_J.items():
        H += edge_zz(n_qubits, i, j, J_ij)
    for v in range(n_qubits):
        H += vertex_x(n_qubits, v, h)
    return H


def edge_H_in_global_convention(n_qubits, i, j, J_ij, h, deg_i, deg_j):
    """Edge Hamiltonian as it would 'split' from convention (B), per the
    manuscript's mediator definition. Each X_v contributes h/deg(v) per
    incident edge. Used for face-resolved associator computation.
    """
    H = edge_zz(n_qubits, i, j, J_ij)
    H += vertex_x(n_qubits, i, h / deg_i)
    H += vertex_x(n_qubits, j, h / deg_j)
    return H


def edge_H_in_edge_resolved_convention(n_qubits, i, j, J_ij, h):
    """Edge Hamiltonian under convention (A): each edge carries (h/2)(X_i+X_j)."""
    H = edge_zz(n_qubits, i, j, J_ij)
    H += vertex_x(n_qubits, i, h / 2.0)
    H += vertex_x(n_qubits, j, h / 2.0)
    return H


def mediator_associator(H_mu, H_a, H_b):
    return 0.25 * commutator(H_mu, commutator(H_a, H_b))


def frob_norm_sq(M):
    return float(np.real(np.trace(M.conj().T @ M)))


def find_triangles(edge_set, n_qubits):
    """Return list of triangles (sorted tuples of 3 vertices)."""
    triangles = []
    vertices = list(range(n_qubits))
    for tri in combinations(vertices, 3):
        a, b, c = tri
        if (
            tuple(sorted((a, b))) in edge_set
            and tuple(sorted((b, c))) in edge_set
            and tuple(sorted((a, c))) in edge_set
        ):
            triangles.append(tri)
    return triangles


def degree_map(edges, n_qubits):
    deg = {v: 0 for v in range(n_qubits)}
    for (i, j) in edges:
        deg[i] += 1
        deg[j] += 1
    return deg


def compute_J_tot_sq(n_qubits, edges_with_J, h, convention):
    """Sum ||A_mu||_F^2 over all (face, mediator-edge) pairs."""
    edges_sorted = {tuple(sorted(e)): J for e, J in edges_with_J.items()}
    edge_set = set(edges_sorted.keys())
    triangles = find_triangles(edge_set, n_qubits)
    deg = degree_map(edges_sorted.keys(), n_qubits)

    total = 0.0
    per_face = {}
    for tri in triangles:
        a, b, c = tri
        e_ab = tuple(sorted((a, b)))
        e_bc = tuple(sorted((b, c)))
        e_ac = tuple(sorted((a, c)))

        if convention == "A":
            H_ab = edge_H_in_edge_resolved_convention(n_qubits, a, b, edges_sorted[e_ab], h)
            H_bc = edge_H_in_edge_resolved_convention(n_qubits, b, c, edges_sorted[e_bc], h)
            H_ac = edge_H_in_edge_resolved_convention(n_qubits, a, c, edges_sorted[e_ac], h)
        elif convention == "B":
            H_ab = edge_H_in_global_convention(n_qubits, a, b, edges_sorted[e_ab], h, deg[a], deg[b])
            H_bc = edge_H_in_global_convention(n_qubits, b, c, edges_sorted[e_bc], h, deg[b], deg[c])
            H_ac = edge_H_in_global_convention(n_qubits, a, c, edges_sorted[e_ac], h, deg[a], deg[c])
        else:
            raise ValueError(convention)

        A_ab = mediator_associator(H_ab, H_bc, H_ac)
        A_bc = mediator_associator(H_bc, H_ab, H_ac)
        A_ac = mediator_associator(H_ac, H_ab, H_bc)

        face_J2 = frob_norm_sq(A_ab) + frob_norm_sq(A_bc) + frob_norm_sq(A_ac)
        per_face[tri] = face_J2
        total += face_J2

    return total, per_face, triangles, deg


def predict_face_additive(edges_with_J, h, triangles, convention, deg):
    """Predicted J_tot^2 if face-additive with appropriate per-face h_eff."""
    edges_sorted = {tuple(sorted(e)): J for e, J in edges_with_J.items()}
    total = 0.0
    per_face = {}
    for tri in triangles:
        a, b, c = tri
        e_ab = tuple(sorted((a, b)))
        e_bc = tuple(sorted((b, c)))
        e_ac = tuple(sorted((a, c)))
        J2_sum = (
            edges_sorted[e_ab] ** 2
            + edges_sorted[e_bc] ** 2
            + edges_sorted[e_ac] ** 2
        )
        if convention == "A":
            # In edge-resolved convention, mediator face uses (h/2)+(h/2)=h
            # per vertex *within the face only*, since each face's three
            # edges contribute (h/2) twice per vertex. So h_eff = h.
            h_eff = h
        elif convention == "B":
            # In global convention, the X_v term is split as h/deg(v) per
            # incident edge. Within a triangle face, each vertex receives
            # contributions from 2 of its incident edges (the two in the face).
            # So effective field at vertex v on face f is h * 2/deg(v).
            # NOTE: this means face-additivity will NOT hold cleanly in (B)
            # if vertices on the face have different degrees. We compute
            # the prediction anyway to see how badly it fails.
            #
            # As a simple guess, use the geometric mean across face vertices.
            h_eff_a = h * 2 / deg[a]
            h_eff_b = h * 2 / deg[b]
            h_eff_c = h * 2 / deg[c]
            h_eff = (h_eff_a * h_eff_b * h_eff_c) ** (1 / 3)
        per_face[tri] = h_eff**4 * J2_sum
        total += h_eff**4 * J2_sum
    return total, per_face


def run_graph_test(name, n_qubits, edges_with_J, h):
    print(f"\n{'='*70}")
    print(f"Graph: {name}")
    print(f"  n_qubits = {n_qubits}, edges = {list(edges_with_J.keys())}")
    print(f"  h = {h}")
    print(f"{'='*70}")

    for convention in ["A", "B"]:
        print(f"\n--- Convention {convention} "
              f"({'edge-resolved' if convention == 'A' else 'global'}) ---")
        J2_num, per_face_num, triangles, deg = compute_J_tot_sq(
            n_qubits, edges_with_J, h, convention
        )
        J2_pred, per_face_pred = predict_face_additive(
            edges_with_J, h, triangles, convention, deg
        )
        print(f"  Found {len(triangles)} triangle face(s): {triangles}")
        print(f"  Vertex degrees: {deg}")
        for tri in triangles:
            print(f"    face {tri}: numerical={per_face_num[tri]:.6f}, "
                  f"predicted (face-additive)={per_face_pred[tri]:.6f}, "
                  f"ratio={per_face_num[tri]/per_face_pred[tri] if per_face_pred[tri]>1e-15 else float('nan'):.4f}")
        print(f"  TOTAL J_tot^2 numerical:  {J2_num:.6f}")
        print(f"  TOTAL J_tot^2 predicted:  {J2_pred:.6f}")
        print(f"  abs diff: {abs(J2_num - J2_pred):.4e}")
        print(f"  Face-additive holds: {abs(J2_num - J2_pred) < 1e-10}")


if __name__ == "__main__":
    # ---- Test 1: single triangle (sanity, both conventions equal) ----
    run_graph_test(
        "Single triangle K_3",
        n_qubits=3,
        edges_with_J={(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0},
        h=1.0,
    )

    # ---- Test 2: two triangles sharing an edge ("diamond" K_4 minus one edge) ----
    # Vertices 0,1,2,3. Triangles: (0,1,2) and (0,1,3). Shared edge: (0,1).
    run_graph_test(
        "Diamond / two triangles sharing edge (0,1)",
        n_qubits=4,
        edges_with_J={
            (0, 1): 1.0,
            (0, 2): 1.0,
            (1, 2): 1.0,
            (0, 3): 1.0,
            (1, 3): 1.0,
        },
        h=1.0,
    )

    # ---- Test 3: K_4 (tetrahedron, 4 triangle faces, deg 3) ----
    run_graph_test(
        "K_4 tetrahedron (4 triangle faces)",
        n_qubits=4,
        edges_with_J={
            (0, 1): 1.0,
            (0, 2): 1.0,
            (0, 3): 1.0,
            (1, 2): 1.0,
            (1, 3): 1.0,
            (2, 3): 1.0,
        },
        h=1.0,
    )

    # ---- Test 4: bowtie (two triangles sharing a vertex) ----
    # Vertices 0,1,2 form one triangle; 0,3,4 form another. Shared: vertex 0.
    run_graph_test(
        "Bowtie / two triangles sharing vertex 0",
        n_qubits=5,
        edges_with_J={
            (0, 1): 1.0,
            (1, 2): 1.0,
            (0, 2): 1.0,
            (0, 3): 1.0,
            (3, 4): 1.0,
            (0, 4): 1.0,
        },
        h=1.0,
    )

    # ---- Test 5: K_4 with asymmetric J's, h=0.7 ----
    run_graph_test(
        "K_4 with asymmetric couplings",
        n_qubits=4,
        edges_with_J={
            (0, 1): 0.5,
            (0, 2): 1.2,
            (0, 3): 0.3,
            (1, 2): 0.9,
            (1, 3): 1.7,
            (2, 3): 0.4,
        },
        h=0.7,
    )
