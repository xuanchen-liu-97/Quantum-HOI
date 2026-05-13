"""
Step 1: Numerical verification of J^2_tot = h^4 * sum_ij J_ij^2
on a single triangle ABC under the TF-Ising Hamiltonian.

This re-establishes the Phase3 Step 1 result with explicit
construction, as a baseline before generalizing to arbitrary graphs.

H_ij^TF = J_ij Z_i Z_j + (h/2)(X_i + X_j)         (edge-resolved form)
H_face  = sum_{(ij) in face} H_ij^TF
        = J_AB Z_A Z_B + J_BC Z_B Z_C + J_AC Z_A Z_C + h(X_A + X_B + X_C)

Mediator associator (3-qubit form):
A_mu = (1/4) [H_mu, [H_other1, H_other2]]
J_tot^2 = ||A_AB||_F^2 + ||A_BC||_F^2 + ||A_AC||_F^2

Predicted: J_tot^2 = h^4 (J_AB^2 + J_BC^2 + J_AC^2)
"""

import numpy as np

# Pauli matrices
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
    """sites_ops = {site_index: pauli_op}; identities elsewhere."""
    ops = [I2] * n
    for site, op in sites_ops.items():
        ops[site] = op
    return kron_all(ops)


def commutator(A, B):
    return A @ B - B @ A


def edge_TF_hamiltonian(n, i, j, J_ij, h):
    """H_ij^TF = J_ij Z_i Z_j + (h/2)(X_i + X_j) on n-qubit space."""
    H = J_ij * embed(n, {i: sz, j: sz})
    H += (h / 2.0) * embed(n, {i: sx})
    H += (h / 2.0) * embed(n, {j: sx})
    return H


def mediator_associator(H_mu, H_a, H_b):
    """A_mu = (1/4) [H_mu, [H_a, H_b]] (nested commutator form)."""
    return 0.25 * commutator(H_mu, commutator(H_a, H_b))


def frob_norm_sq(M):
    return float(np.real(np.trace(M.conj().T @ M)))


def run_single_triangle(J_AB, J_BC, J_AC, h, verbose=True):
    """Build H on the 3-qubit triangle, compute the three mediator
    associators, and return J_tot^2 along with the analytic prediction."""
    n = 3
    A, B, C = 0, 1, 2

    # Edge Hamiltonians; note we use the *full* sum-of-edges convention,
    # which gives total H_face with each X_i appearing twice (once per
    # incident edge), matching the manuscript convention H = sum_ij H_ij.
    H_AB = edge_TF_hamiltonian(n, A, B, J_AB, h)
    H_BC = edge_TF_hamiltonian(n, B, C, J_BC, h)
    H_AC = edge_TF_hamiltonian(n, A, C, J_AC, h)

    # Mediator associators on this face
    A_AB = mediator_associator(H_AB, H_BC, H_AC)
    A_BC = mediator_associator(H_BC, H_AB, H_AC)
    A_AC = mediator_associator(H_AC, H_AB, H_BC)

    j_AB_sq = frob_norm_sq(A_AB)
    j_BC_sq = frob_norm_sq(A_BC)
    j_AC_sq = frob_norm_sq(A_AC)
    J_tot_sq = j_AB_sq + j_BC_sq + j_AC_sq

    # Analytic prediction
    J_tot_sq_predicted = h**4 * (J_AB**2 + J_BC**2 + J_AC**2)

    if verbose:
        print(f"  J_AB={J_AB}, J_BC={J_BC}, J_AC={J_AC}, h={h}")
        print(f"    j_AB^2 = {j_AB_sq:.10f}")
        print(f"    j_BC^2 = {j_BC_sq:.10f}")
        print(f"    j_AC^2 = {j_AC_sq:.10f}")
        print(f"    J_tot^2 (numerical)  = {J_tot_sq:.10f}")
        print(f"    J_tot^2 (predicted)  = {J_tot_sq_predicted:.10f}")
        print(f"    abs error            = {abs(J_tot_sq - J_tot_sq_predicted):.2e}")
    return J_tot_sq, J_tot_sq_predicted


if __name__ == "__main__":
    print("=" * 70)
    print("Step 1: TF-Ising single-triangle verification of J^2_tot formula")
    print("=" * 70)

    test_cases = [
        # Baseline: equal couplings
        (1.0, 1.0, 1.0, 1.0),
        # Pure Ising (h=0): should be exactly zero
        (1.0, 1.0, 1.0, 0.0),
        # Asymmetric edges
        (0.7, 1.3, 0.5, 0.8),
        # Matches manuscript convention sqrt(3) total norm
        (np.sqrt(1.0), np.sqrt(1.0), np.sqrt(1.0), 1.5),
        # Negative couplings (FM-like)
        (-1.0, -1.0, -1.0, 1.0),
        # Mixed signs (frustrated): should not affect J_tot^2
        (1.0, 1.0, -1.0, 1.0),
        # Random-ish
        (0.3, -0.9, 1.7, 0.6),
    ]

    all_passed = True
    max_err = 0.0
    for case in test_cases:
        print()
        num, pred = run_single_triangle(*case)
        err = abs(num - pred)
        max_err = max(max_err, err)
        if err > 1e-10:
            all_passed = False
            print("    *** MISMATCH ***")

    print()
    print("=" * 70)
    print(f"Maximum absolute error across {len(test_cases)} cases: {max_err:.2e}")
    print(f"All passed (tol 1e-10): {all_passed}")
    print("=" * 70)
