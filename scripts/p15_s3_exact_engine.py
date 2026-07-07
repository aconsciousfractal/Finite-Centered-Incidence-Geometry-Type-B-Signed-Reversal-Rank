# -*- coding: utf-8 -*-
"""
P15-S3 exact replay engine.

This script is the deterministic, stdlib-only consolidation of the P15 discovery
scripts needed for the S3 gate.  It rebuilds the signed group B_n, the layer
BtildeI_n(k,l), the reversal rho, the sign-aware V_ref channel M_ref, and the
sign-blind pair channel M_pair.  All ranks and projector checks use exact
integer/Fraction arithmetic.  No numpy, no floats, no randomized checks.

The S3 scope is finite replay and exact local certification.  It does not prove
the H_m recurrence certificate, odd-n pair positivity in general, Type-D scout
claims, or the C2xC2 square theorem; those remain later gates.

ASCII output only.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


SCHEMA = "p15.s3.exact_engine.v1"


# ---------------------------------------------------------------------------
# Basic exact linear algebra


def f(x: int | Fraction) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x, 1)


def zero_matrix(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


def identity_matrix(n: int) -> list[list[Fraction]]:
    return [[Fraction(1 if i == j else 0, 1) for j in range(n)] for i in range(n)]


def rank_exact(matrix: list[list[int | Fraction]]) -> int:
    """Reduced-row-echelon rank over Q."""
    if not matrix:
        return 0
    mat = [[f(x) for x in row] for row in matrix]
    rows = len(mat)
    cols = len(mat[0])
    r = 0
    for c in range(cols):
        pivot = None
        for i in range(r, rows):
            if mat[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        pv = mat[r][c]
        mat[r] = [x / pv for x in mat[r]]
        for i in range(rows):
            if i != r and mat[i][c] != 0:
                q = mat[i][c]
                mat[i] = [mat[i][j] - q * mat[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def matmul(
    a: list[list[int | Fraction]],
    b: list[list[int | Fraction]],
) -> list[list[Fraction]]:
    rows = len(a)
    mid = len(b)
    cols = len(b[0]) if b else 0
    out: list[list[Fraction]] = [[Fraction(0, 1) for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = f(a[i][k])
            if aik == 0:
                continue
            for j in range(cols):
                out[i][j] += aik * f(b[k][j])
    return out


def standard_projected_matrix(matrix: list[list[int]]) -> list[list[Fraction]]:
    """Return P M P for P=I-J/n; rank equals rank on 1-perp."""
    n = len(matrix)
    p = [[Fraction(1 if i == j else 0, 1) - Fraction(1, n) for j in range(n)] for i in range(n)]
    return matmul(matmul(p, matrix), p)


def matrix_equal(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> bool:
    if len(a) != len(b):
        return False
    for row_a, row_b in zip(a, b):
        if len(row_a) != len(row_b):
            return False
        for xa, xb in zip(row_a, row_b):
            if f(xa) != f(xb):
                return False
    return True


# ---------------------------------------------------------------------------
# B_n and the BtildeI layers


SignedPerm = tuple[tuple[int, ...], tuple[int, ...]]


def reversal(n: int) -> tuple[int, ...]:
    return tuple(n - 1 - i for i in range(n))


def signed_perms(n: int):
    for pi in itertools.permutations(range(n)):
        for eps in itertools.product((1, -1), repeat=n):
            yield pi, eps


def positive_hit_counts(pi: tuple[int, ...], eps: tuple[int, ...], rho: tuple[int, ...]) -> tuple[int, int]:
    n = len(pi)
    k = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
    l = sum(1 for i in range(n) if pi[i] == rho[i] and eps[i] == 1)
    return k, l


def fixed_count(pi: tuple[int, ...]) -> int:
    return sum(1 for i, j in enumerate(pi) if i == j)


def neg_fixed_count(pi: tuple[int, ...], eps: tuple[int, ...]) -> int:
    return sum(1 for i in range(len(pi)) if pi[i] == i and eps[i] == -1)


def build_layers(n: int) -> dict[tuple[int, int], list[SignedPerm]]:
    rho = reversal(n)
    layers: dict[tuple[int, int], list[SignedPerm]] = defaultdict(list)
    for pi, eps in signed_perms(n):
        layers[positive_hit_counts(pi, eps, rho)].append((pi, eps))
    return dict(layers)


def channel_matrices(layer: list[SignedPerm], n: int) -> tuple[list[list[int]], list[list[int]]]:
    """Return M_ref and M_pair as row=value, col=position matrices."""
    m_ref = zero_matrix(n, n)
    m_pair = zero_matrix(n, n)
    for pi, eps in layer:
        for i in range(n):
            j = pi[i]
            m_ref[j][i] += eps[i]
            m_pair[j][i] += 1
    return m_ref, m_pair


def projector_plus_matrix(n: int) -> list[list[Fraction]]:
    rho = reversal(n)
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        out[i][i] += Fraction(1, 2)
        out[rho[i]][i] += Fraction(1, 2)
    return out


def projector_plus_pairs_matrix(n: int) -> list[list[Fraction]]:
    """Odd-n projector on non-center reversal pairs, zero on the center."""
    rho = reversal(n)
    center = (n - 1) // 2
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        if i == center:
            continue
        out[i][i] += Fraction(1, 2)
        out[rho[i]][i] += Fraction(1, 2)
    return out


def scalar_times_matrix(scalar: int | Fraction, matrix: list[list[int | Fraction]]) -> list[list[Fraction]]:
    return [[f(scalar) * f(x) for x in row] for row in matrix]


def add_matrices(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> list[list[Fraction]]:
    return [[f(x) + f(y) for x, y in zip(row_a, row_b)] for row_a, row_b in zip(a, b)]


def center_idempotent(n: int) -> list[list[Fraction]]:
    c = (n - 1) // 2
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    out[c][c] = Fraction(1, 1)
    return out


# ---------------------------------------------------------------------------
# Exact closed forms from the excerpt


def poly_pow_2var(base: dict[tuple[int, int], int], power: int) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {(0, 0): 1}
    for _ in range(power):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (i, j), coeff in out.items():
            for (a, b), d in base.items():
                nxt[(i + a, j + b)] += coeff * d
        out = dict(nxt)
    return out


def h_even_table(n: int) -> dict[tuple[int, int], int]:
    """H(f,r) for even n via the two-alphabet rook GF."""
    assert n % 2 == 0
    base = {(2, 0): 1, (1, 0): 2, (0, 2): 1, (0, 1): 2, (0, 0): 1}
    n_table = poly_pow_2var(base, n // 2)
    h: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b), count in n_table.items():
        if count == 0:
            continue
        fac = factorial(n - a - b) * count
        for ff in range(a + 1):
            cf = comb(a, ff) * ((-1) ** (a - ff))
            for rr in range(b + 1):
                cr = comb(b, rr) * ((-1) ** (b - rr))
                h[(ff, rr)] += fac * cf * cr
    return dict(h)


def size_even(n: int, k: int, l: int, h_table: dict[tuple[int, int], int] | None = None) -> int:
    if h_table is None:
        h_table = h_even_table(n)
    total = 0
    for (ff, rr), h in h_table.items():
        total += h * comb(ff, k) * comb(rr, l) * (2 ** (n - ff - rr))
    return total


def a_even_closed_form(n: int, k: int, h_table: dict[tuple[int, int], int] | None = None) -> int:
    """Excerpt 9.7 Step 4 formula for the non-center M_ref scalar, even n."""
    if h_table is None:
        h_table = h_even_table(n)
    numerator = 0
    for (ff, rr), h in h_table.items():
        if ff < 1:
            continue
        cf = (comb(ff - 1, k - 1) if k >= 1 else 0) - comb(ff - 1, k)
        numerator += ff * h * cf * comb(rr, k) * (2 ** (n - ff - rr))
    assert numerator % n == 0
    return numerator // n


def d_even_diag(n: int, p: int, h_table: dict[tuple[int, int], int] | None = None) -> int:
    if p < 0:
        return 0
    return size_even(n, p, p, h_table)


# ---------------------------------------------------------------------------
# Permanent by Ryser and the domination replay


def permanent_ryser(matrix: list[list[int]]) -> int:
    n = len(matrix)
    if n == 0:
        return 1
    total = 0
    cols = range(n)
    for mask in range(1, 1 << n):
        bits = mask.bit_count()
        prod = 1
        for i in range(n):
            row_sum = 0
            for j in cols:
                if mask & (1 << j):
                    row_sum += matrix[i][j]
            prod *= row_sum
            if prod == 0:
                break
        total += ((-1) ** (n - bits)) * prod
    return total


def minor_matrix(matrix: list[list[int]], row: int, col: int) -> list[list[int]]:
    return [
        [matrix[i][j] for j in range(len(matrix)) if j != col]
        for i in range(len(matrix))
        if i != row
    ]


def residual_matrix_for_pattern(n: int, p_set: frozenset[int], q_set: frozenset[int]):
    """Residual weighted bijection matrix after fixing positive id/rev hit sets."""
    rho = reversal(n)
    center = (n - 1) // 2 if n % 2 else None
    if len(p_set) != len(q_set):
        return None
    if center is not None and ((center in p_set) != (center in q_set)):
        return None
    for i in p_set & q_set:
        if rho[i] != i:
            return None

    forced: dict[int, int] = {}
    for i in p_set:
        forced[i] = i
    for i in q_set:
        if i in forced and forced[i] != rho[i]:
            return None
        forced[i] = rho[i]
    images = list(forced.values())
    if len(set(images)) != len(images):
        return None

    domain = [i for i in range(n) if i not in forced]
    codomain = [j for j in range(n) if j not in set(images)]
    if len(domain) != len(codomain):
        return None

    matrix: list[list[int]] = []
    d_cells: list[tuple[int, int]] = []
    for row, i in enumerate(domain):
        mat_row = []
        for col, j in enumerate(codomain):
            weight = 1 if (j == i or j == rho[i]) else 2
            mat_row.append(weight)
            if j == i:
                d_cells.append((row, col))
        matrix.append(mat_row)
    return matrix, d_cells


def check_pattern_permanent(n: int, p_set: frozenset[int], q_set: frozenset[int]):
    residual = residual_matrix_for_pattern(n, p_set, q_set)
    if residual is None:
        return None
    matrix, d_cells = residual
    denom = permanent_ryser(matrix)
    numer = sum(permanent_ryser(minor_matrix(matrix, row, col)) for row, col in d_cells)
    if denom == 0:
        return None
    return numer, denom


# ---------------------------------------------------------------------------
# Signed permutation multiplication for centrality checks


def mul_signed(a: SignedPerm, b: SignedPerm) -> SignedPerm:
    pa, ea = a
    pb, eb = b
    n = len(pa)
    p = tuple(pa[pb[i]] for i in range(n))
    e = tuple(ea[pb[i]] * eb[i] for i in range(n))
    return p, e


def inv_signed(a: SignedPerm) -> SignedPerm:
    pa, ea = a
    n = len(pa)
    p = [0] * n
    e = [0] * n
    for i in range(n):
        p[pa[i]] = i
        e[pa[i]] = ea[i]
    return tuple(p), tuple(e)


def conjugate(g: SignedPerm, h: SignedPerm) -> SignedPerm:
    return mul_signed(mul_signed(h, g), inv_signed(h))


# ---------------------------------------------------------------------------
# R1 enumeration replay


def positive_fixed_marginal(n: int) -> list[int]:
    out = [0] * (n + 1)
    for pi in itertools.permutations(range(n)):
        ff = fixed_count(pi)
        for k in range(ff + 1):
            out[k] += comb(ff, k) * (2 ** (n - ff))
    return out


def dplus(n: int) -> int:
    return positive_fixed_marginal(n)[0]


def run_r1() -> dict:
    d_seq = [dplus(n) for n in range(0, 9)]
    expected = [1, 1, 5, 29, 233, 2329, 27949, 391285, 6260561]
    recurrence_ok = all(d_seq[n] == 2 * n * d_seq[n - 1] + ((-1) ** n) for n in range(1, len(d_seq)))
    binomial_ok = True
    marginals = {}
    for n in range(1, 8):
        marginal = positive_fixed_marginal(n)
        marginals[str(n)] = marginal
        for k, value in enumerate(marginal):
            if value != comb(n, k) * d_seq[n - k]:
                binomial_ok = False
    return {
        "id": "R1_A000354_replay",
        "status": "PASS" if d_seq == expected and recurrence_ok and binomial_ok else "FAIL",
        "range": "n=0..8 for Dplus, n=1..7 for marginals",
        "dplus_sequence": d_seq,
        "expected_A000354_prefix": expected,
        "recurrence_ok": recurrence_ok,
        "binomial_distribution_ok": binomial_ok,
        "marginals_n_1_to_7": marginals,
        "boundary": "Enumeration is classical: A000354 positive-fixed marginal and Type-A two-diagonal lift; not A007016 for the signed count.",
    }


# ---------------------------------------------------------------------------
# R2 orbital algebra replay


def centralizer_rho(n: int) -> list[tuple[int, ...]]:
    rho = reversal(n)
    return [
        tau
        for tau in itertools.permutations(range(n))
        if all(tau[rho[i]] == rho[tau[i]] for i in range(n))
    ]


def pair_orbits(n: int, group: list[tuple[int, ...]]) -> dict[tuple[int, int], int]:
    orbit: dict[tuple[int, int], int] = {}
    oid = 0
    for i in range(n):
        for j in range(n):
            if (i, j) in orbit:
                continue
            stack = [(i, j)]
            seen = set()
            while stack:
                a, b = stack.pop()
                if (a, b) in seen:
                    continue
                seen.add((a, b))
                for tau in group:
                    stack.append((tau[a], tau[b]))
            for item in seen:
                orbit[item] = oid
            oid += 1
    return orbit


def values_on_orbits(matrix: list[list[int]], orbit: dict[tuple[int, int], int]):
    values: dict[int, int] = {}
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            oid = orbit[(i, j)]
            value = matrix[j][i]
            if oid in values and values[oid] != value:
                return False, values
            values[oid] = value
    return True, dict(sorted(values.items()))


def orbital_matrices(n: int, orbit: dict[tuple[int, int], int]) -> list[list[list[int]]]:
    count = max(orbit.values()) + 1
    mats = []
    for oid in range(count):
        mats.append([[1 if orbit[(i, j)] == oid else 0 for i in range(n)] for j in range(n)])
    return mats


def is_constant_on_orbits(matrix: list[list[int]], orbit: dict[tuple[int, int], int]) -> bool:
    return values_on_orbits(matrix, orbit)[0]


def orbital_closure_and_commutativity(n: int, orbit: dict[tuple[int, int], int]) -> tuple[bool, bool]:
    mats = orbital_matrices(n, orbit)
    closed = True
    commutative = True
    for a in range(len(mats)):
        for b in range(len(mats)):
            ab = matmul(mats[a], mats[b])
            ba = matmul(mats[b], mats[a])
            if not is_constant_on_orbits([[int(x) for x in row] for row in ab], orbit):
                closed = False
            if not matrix_equal(ab, ba):
                commutative = False
    return closed, commutative


def exact_noncentral_witness(n: int, k: int) -> dict:
    layers = build_layers(n)
    x = layers[(k, k)]
    x_set = set(x)
    group = list(signed_perms(n))
    for h in group:
        for g in x:
            cg = conjugate(g, h)
            if cg not in x_set:
                return {
                    "checked": True,
                    "is_central": False,
                    "n": n,
                    "k": k,
                    "witness_g": [list(g[0]), list(g[1])],
                    "witness_h": [list(h[0]), list(h[1])],
                    "conjugate_kl": list(positive_hit_counts(cg[0], cg[1], reversal(n))),
                }
    return {"checked": True, "is_central": True, "n": n, "k": k}


def run_r2() -> dict:
    cases = [(4, 1), (4, 2), (5, 1), (5, 2), (6, 1)]
    rows = []
    all_ok = True
    noncentral = exact_noncentral_witness(4, 1)
    if noncentral.get("is_central"):
        all_ok = False
    for n, k in cases:
        layers = build_layers(n)
        x = layers[(k, k)]
        m_ref, m_pair = channel_matrices(x, n)
        group = centralizer_rho(n)
        orbit = pair_orbits(n, group)
        ref_ok, ref_values = values_on_orbits(m_ref, orbit)
        pair_ok, _ = values_on_orbits(m_pair, orbit)
        closed, commutative = orbital_closure_and_commutativity(n, orbit)
        rank_ref = rank_exact(m_ref)
        rank_pair_std = rank_exact(standard_projected_matrix(m_pair))
        dim_plus = (n + 1) // 2
        case_ok = (
            ref_ok
            and pair_ok
            and closed
            and rank_ref == dim_plus
            and rank_pair_std == dim_plus - 1
        )
        if not case_ok:
            all_ok = False
        rows.append(
            {
                "n": n,
                "k": k,
                "size": len(x),
                "num_orbits": max(orbit.values()) + 1,
                "M_ref_in_orbital_algebra": ref_ok,
                "M_pair_in_orbital_algebra": pair_ok,
                "orbital_algebra_closed": closed,
                "orbital_algebra_commutative": commutative,
                "M_ref_orbit_values": ref_values,
                "rank_M_ref": rank_ref,
                "rank_M_pair_std": rank_pair_std,
                "expected_ref_rank": dim_plus,
                "expected_pair_std_rank": dim_plus - 1,
            }
        )
    return {
        "id": "R2_matching_scheme_replay",
        "status": "PASS" if all_ok else "FAIL",
        "range": "n=4,5,6 selected diagonal cases",
        "noncentrality_exact_witness": noncentral,
        "rows": rows,
        "boundary": "Scheme gives the orbital algebra and idempotent multiplicity; P15 contribution is collapse plus positivity.",
    }


# ---------------------------------------------------------------------------
# Channel/collapse/scalar checks


def actual_lambda_pair_even(m_pair: list[list[int]], n: int) -> int:
    rho = reversal(n)
    i0 = 0
    ri0 = rho[i0]
    generic = [j for j in range(n) if j not in (i0, ri0)]
    cp = m_pair[generic[0]][i0] if generic else 0
    return m_pair[i0][i0] + m_pair[ri0][i0] - 2 * cp


def channel_case(n: int, k: int, layers: dict[tuple[int, int], list[SignedPerm]]) -> dict:
    x = layers[(k, k)]
    m_ref, m_pair = channel_matrices(x, n)
    rank_ref = rank_exact(m_ref)
    rank_pair_std = rank_exact(standard_projected_matrix(m_pair))
    expected_ref = (n + 1) // 2
    expected_pair = expected_ref - 1
    snf = sum(neg_fixed_count(pi, eps) for pi, eps in x)
    size = len(x)
    rho = reversal(n)
    center = (n - 1) // 2 if n % 2 else None
    i0 = 0 if center != 0 else 1
    a_actual = m_ref[i0][i0]
    result = {
        "n": n,
        "k": k,
        "size": size,
        "sum_neg_fixed": snf,
        "a_actual": a_actual,
        "rank_M_ref": rank_ref,
        "rank_M_pair_std": rank_pair_std,
        "expected_rank_M_ref": expected_ref,
        "expected_rank_M_pair_std": expected_pair,
        "rank_ok": rank_ref == expected_ref and rank_pair_std == expected_pair,
    }
    if n % 2 == 0:
        expected = scalar_times_matrix(2 * a_actual, projector_plus_matrix(n))
        h = h_even_table(n)
        a_closed = a_even_closed_form(n, k, h)
        layer_formula_size = size_even(n, k, k, h)
        trace_rhs = k * size - (k + 1) * len(layers.get((k + 1, k), []))
        lambda_formula = Fraction(2 * ((k - 1) * size + snf), n - 2)
        lambda_actual = actual_lambda_pair_even(m_pair, n)
        result.update(
            {
                "parity": "even",
                "collapse_ok": matrix_equal(m_ref, expected),
                "a_closed_form": a_closed,
                "a_closed_form_ok": a_actual == a_closed,
                "size_closed_form_ok": size == layer_formula_size,
                "trace_flip_ok": n * a_actual == trace_rhs,
                "lambda_pair_actual": str(lambda_actual),
                "lambda_pair_formula": str(lambda_formula),
                "lambda_pair_formula_ok": Fraction(lambda_actual, 1) == lambda_formula,
            }
        )
    else:
        assert center is not None
        a_c = m_ref[center][center]
        expected = add_matrices(
            scalar_times_matrix(2 * a_actual, projector_plus_pairs_matrix(n)),
            scalar_times_matrix(a_c, center_idempotent(n)),
        )
        h_even = h_even_table(n - 1)
        a_c_closed = d_even_diag(n - 1, k - 1, h_even) - d_even_diag(n - 1, k, h_even)
        trace_rhs = k * size - snf
        result.update(
            {
                "parity": "odd",
                "a_c_actual": a_c,
                "collapse_ok": matrix_equal(m_ref, expected),
                "a_c_closed_form": a_c_closed,
                "a_c_closed_form_ok": a_c == a_c_closed,
                "odd_trace_ok": (n - 1) * a_actual + a_c == trace_rhs,
            }
        )
    return result


def run_channels() -> dict:
    cases_by_n = {
        4: [1, 2],
        5: [1, 2],
        6: [1, 2, 3],
        7: [1, 2, 3],
    }
    rows = []
    all_ok = True
    for n, ks in cases_by_n.items():
        layers = build_layers(n)
        for k in ks:
            if (k, k) not in layers or not layers[(k, k)]:
                continue
            row = channel_case(n, k, layers)
            rows.append(row)
            row_ok = row["rank_ok"] and row["collapse_ok"]
            if row["parity"] == "even":
                row_ok = row_ok and row["a_closed_form_ok"] and row["size_closed_form_ok"] and row["trace_flip_ok"] and row["lambda_pair_formula_ok"]
            else:
                row_ok = row_ok and row["a_c_closed_form_ok"] and row["odd_trace_ok"]
            if not row_ok:
                all_ok = False
    return {
        "id": "ambient_channel_exact_replay",
        "status": "PASS" if all_ok else "FAIL",
        "range": "n=4..7 diagonal nonempty k>=1 cases; even closed forms checked for n=4,6",
        "rows": rows,
    }


def run_closed_form_sweep() -> dict:
    rows = []
    all_ok = True
    for n in range(4, 13, 2):
        h = h_even_table(n)
        values = []
        for k in range(1, n // 2 + 1):
            size = size_even(n, k, k, h)
            if size == 0:
                continue
            a_value = a_even_closed_form(n, k, h)
            values.append({"k": k, "size": size, "a": a_value, "a_positive": a_value > 0})
            if a_value <= 0:
                all_ok = False
        rows.append({"n": n, "values": values})
    return {
        "id": "even_closed_form_sweep",
        "status": "PASS" if all_ok else "FAIL",
        "range": "even n=4..12, nonempty k>=1",
        "rows": rows,
    }


def run_permanent_domination() -> dict:
    rows = []
    all_ok = True
    for n in (4, 5, 6):
        for k in range(1, (n + 1) // 2 + 1):
            checked = 0
            feasible = 0
            max_ratio = Fraction(0, 1)
            max_pattern = None
            for p_tuple in itertools.combinations(range(n), k):
                p_set = frozenset(p_tuple)
                for q_tuple in itertools.combinations(range(n), k):
                    q_set = frozenset(q_tuple)
                    out = check_pattern_permanent(n, p_set, q_set)
                    if out is None:
                        continue
                    feasible += 1
                    numer, denom = out
                    checked += 1
                    ratio = Fraction(numer, denom)
                    if ratio > max_ratio:
                        max_ratio = ratio
                        max_pattern = {"P": sorted(p_set), "Q": sorted(q_set)}
                    if numer > denom:
                        all_ok = False
            if checked:
                rows.append(
                    {
                        "n": n,
                        "k": k,
                        "feasible_patterns": feasible,
                        "checked_patterns": checked,
                        "max_E_neg_fixed": str(max_ratio),
                        "max_pattern": max_pattern,
                        "all_patterns_le_1": max_ratio <= 1,
                    }
                )
    return {
        "id": "permanent_domination_ryser_replay",
        "status": "PASS" if all_ok else "FAIL",
        "range": "all feasible positive-hit patterns for n=4,5,6 and k>=1 up to ceil(n/2)",
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# Report assembly


def make_report() -> dict:
    checks = [
        run_r1(),
        run_r2(),
        run_channels(),
        run_closed_form_sweep(),
        run_permanent_domination(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S3 engine and replay plan",
        "engine": "scripts/p15_s3_exact_engine.py",
        "arithmetic": "exact integer and Fraction arithmetic; no floats; no numpy",
        "status": "PASS" if not failed else "FAIL",
        "failed_checks": failed,
        "checks": checks,
        "non_goals": [
            "No H_m recurrence certificate; deferred to P15-S8.",
            "No general odd-n pair-standard positivity proof; deferred to P15-S5/S7.",
            "No Type-D theorem promotion; Dtilde remains scout-only.",
            "No C2xC2 square closeout; deferred to P15-S6.",
            "No public claim promotion.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S3 exact engine")
    print("schema:", report["schema"])
    print("arithmetic:", report["arithmetic"])
    print("status:", report["status"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"], "-", check["range"])
        if check["id"] == "ambient_channel_exact_replay":
            for row in check["rows"]:
                extra = "a=%s" % row["a_actual"]
                if row["parity"] == "odd":
                    extra += " ac=%s" % row["a_c_actual"]
                if row["parity"] == "even":
                    extra += " lambda_pair=%s" % row["lambda_pair_actual"]
                print(
                    "  n={n} k={k} size={size} ranks=({rank_M_ref},{rank_M_pair_std}) "
                    "expected=({expected_rank_M_ref},{expected_rank_M_pair_std}) "
                    "collapse={collapse_ok} {extra}".format(extra=extra, **row)
                )
        elif check["id"] == "permanent_domination_ryser_replay":
            for row in check["rows"]:
                print(
                    "  n={n} k={k} patterns={checked_patterns} maxE={max_E_neg_fixed} le1={all_patterns_le_1}".format(
                        **row
                    )
                )
        elif check["id"] == "R2_matching_scheme_replay":
            print("  noncentral witness is_central:", check["noncentrality_exact_witness"]["is_central"])
            for row in check["rows"]:
                print(
                    "  n={n} k={k} orbits={num_orbits} ranks=({rank_M_ref},{rank_M_pair_std}) "
                    "closed={orbital_algebra_closed} comm={orbital_algebra_commutative}".format(**row)
                )
        elif check["id"] == "R1_A000354_replay":
            print("  Dplus:", check["dplus_sequence"])
            print("  recurrence_ok:", check["recurrence_ok"], "binomial_distribution_ok:", check["binomial_distribution_ok"])
        elif check["id"] == "even_closed_form_sweep":
            for row in check["rows"]:
                vals = ", ".join("k=%d:a=%d" % (v["k"], v["a"]) for v in row["values"])
                print("  n=%d %s" % (row["n"], vals))


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S3 exact replay engine")
    parser.add_argument("--write-json", help="write deterministic JSON report to this path")
    args = parser.parse_args()
    report = make_report()
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
