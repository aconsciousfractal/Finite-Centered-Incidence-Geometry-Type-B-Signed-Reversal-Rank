"""
P15-S7 independent red-team recompute.

This verifier intentionally does not import p15_s3_exact_engine or any S4-S6
certificate script.  It rebuilds the signed permutation loops, exact matrix
ranks, the two standard-channel matrices, and the square/count checks from
stdlib primitives only.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path

SCHEMA = "p15.s7.independent_red_team.v1"
A000354_PREFIX = [1, 1, 5, 29, 233, 2329, 27949, 391285, 6260561]


def rank_exact(matrix: list[list[int | Fraction]]) -> int:
    if not matrix:
        return 0
    mat = [[x if isinstance(x, Fraction) else Fraction(x, 1) for x in row] for row in matrix]
    rows = len(mat)
    cols = len(mat[0]) if rows else 0
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


def matmul(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> list[list[Fraction]]:
    rows = len(a)
    mid = len(b)
    cols = len(b[0]) if b else 0
    out = [[Fraction(0, 1) for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = a[i][k] if isinstance(a[i][k], Fraction) else Fraction(a[i][k], 1)
            if aik == 0:
                continue
            for j in range(cols):
                bkj = b[k][j] if isinstance(b[k][j], Fraction) else Fraction(b[k][j], 1)
                out[i][j] += aik * bkj
    return out


def standard_projected_matrix(matrix: list[list[int]]) -> list[list[Fraction]]:
    n = len(matrix)
    projector = [[Fraction(1 if i == j else 0, 1) - Fraction(1, n) for j in range(n)] for i in range(n)]
    return matmul(matmul(projector, matrix), projector)


def matrix_equal(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> bool:
    if len(a) != len(b):
        return False
    for row_a, row_b in zip(a, b):
        if len(row_a) != len(row_b):
            return False
        for xa, xb in zip(row_a, row_b):
            fa = xa if isinstance(xa, Fraction) else Fraction(xa, 1)
            fb = xb if isinstance(xb, Fraction) else Fraction(xb, 1)
            if fa != fb:
                return False
    return True


def reversal(n: int) -> tuple[int, ...]:
    return tuple(n - 1 - i for i in range(n))


def signed_perms(n: int):
    for pi in itertools.permutations(range(n)):
        for eps in itertools.product((1, -1), repeat=n):
            yield pi, eps


def positive_identity_count(pi: tuple[int, ...], eps: tuple[int, ...]) -> int:
    return sum(1 for i in range(len(pi)) if pi[i] == i and eps[i] == 1)


def positive_reversal_count(pi: tuple[int, ...], eps: tuple[int, ...], rho: tuple[int, ...]) -> int:
    return sum(1 for i in range(len(pi)) if pi[i] == rho[i] and eps[i] == 1)


def negative_sign_count(eps: tuple[int, ...]) -> int:
    return sum(1 for value in eps if value == -1)


def zero_matrix(n: int) -> list[list[int]]:
    return [[0 for _ in range(n)] for _ in range(n)]


def diagonal_layers(n: int) -> dict[int, list[tuple[tuple[int, ...], tuple[int, ...]]]]:
    rho = reversal(n)
    layers: dict[int, list[tuple[tuple[int, ...], tuple[int, ...]]]] = defaultdict(list)
    for pi, eps in signed_perms(n):
        k = positive_identity_count(pi, eps)
        l = positive_reversal_count(pi, eps, rho)
        if k == l and k >= 1:
            layers[k].append((pi, eps))
    return dict(layers)


def channel_matrices(layer: list[tuple[tuple[int, ...], tuple[int, ...]]], n: int) -> tuple[list[list[int]], list[list[int]]]:
    m_ref = zero_matrix(n)
    m_pair = zero_matrix(n)
    for pi, eps in layer:
        for i in range(n):
            j = pi[i]
            m_ref[j][i] += eps[i]
            m_pair[j][i] += 1
    return m_ref, m_pair


def projector_plus(n: int) -> list[list[Fraction]]:
    rho = reversal(n)
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        out[i][i] += Fraction(1, 2)
        out[rho[i]][i] += Fraction(1, 2)
    return out


def projector_plus_pairs(n: int) -> list[list[Fraction]]:
    rho = reversal(n)
    center = (n - 1) // 2
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        if i == center:
            continue
        out[i][i] += Fraction(1, 2)
        out[rho[i]][i] += Fraction(1, 2)
    return out


def center_idempotent(n: int) -> list[list[Fraction]]:
    center = (n - 1) // 2
    out = [[Fraction(0, 1) for _ in range(n)] for _ in range(n)]
    out[center][center] = Fraction(1, 1)
    return out


def scalar_times(scalar: int | Fraction, matrix: list[list[int | Fraction]]) -> list[list[Fraction]]:
    s = scalar if isinstance(scalar, Fraction) else Fraction(scalar, 1)
    return [[s * (x if isinstance(x, Fraction) else Fraction(x, 1)) for x in row] for row in matrix]


def add_matrices(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> list[list[Fraction]]:
    out: list[list[Fraction]] = []
    for row_a, row_b in zip(a, b):
        row: list[Fraction] = []
        for xa, xb in zip(row_a, row_b):
            fa = xa if isinstance(xa, Fraction) else Fraction(xa, 1)
            fb = xb if isinstance(xb, Fraction) else Fraction(xb, 1)
            row.append(fa + fb)
        out.append(row)
    return out


def run_channel_recompute(max_n: int = 7) -> dict:
    rows = []
    ok = True
    for n in range(3, max_n + 1):
        layers = diagonal_layers(n)
        for k in sorted(layers):
            layer = layers[k]
            m_ref, m_pair = channel_matrices(layer, n)
            rank_ref = rank_exact(m_ref)
            rank_pair_std = rank_exact(standard_projected_matrix(m_pair))
            expected_ref = (n + 1) // 2
            expected_pair = expected_ref - 1
            center = (n - 1) // 2 if n % 2 else None
            i0 = 0 if center != 0 else 1
            a = m_ref[i0][i0]
            if n % 2 == 0:
                expected_matrix = scalar_times(2 * a, projector_plus(n))
                a_c = None
                collapse_kind = "even_M_ref_equals_2a_P_plus"
            else:
                assert center is not None
                a_c = m_ref[center][center]
                expected_matrix = add_matrices(
                    scalar_times(2 * a, projector_plus_pairs(n)),
                    scalar_times(a_c, center_idempotent(n)),
                )
                collapse_kind = "odd_M_ref_equals_2a_P_plus_pairs_plus_ac_E_center"
            collapse_ok = matrix_equal(m_ref, expected_matrix)
            base_marker_ok = True
            if n == 3 and k == 1:
                base_marker_ok = a == -1 and a_c == 2 and rank_ref == 2 and rank_pair_std == 1
            row_ok = rank_ref == expected_ref and rank_pair_std == expected_pair and collapse_ok and base_marker_ok
            ok = ok and row_ok
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "size": len(layer),
                    "rank_M_ref": rank_ref,
                    "rank_M_pair_std": rank_pair_std,
                    "expected_rank_M_ref": expected_ref,
                    "expected_rank_M_pair_std": expected_pair,
                    "a": a,
                    "a_c": a_c,
                    "collapse_kind": collapse_kind,
                    "collapse_ok": collapse_ok,
                    "n3_k1_base_marker_ok": base_marker_ok if n == 3 and k == 1 else None,
                    "ok": row_ok,
                }
            )
    even_rows_ok = all(row["ok"] for row in rows if row["n"] % 2 == 0)
    odd_rows_ok = all(row["ok"] for row in rows if row["n"] % 2 == 1)
    return {
        "id": "independent_channel_rank_and_collapse_recompute",
        "status": "PASS" if ok else "FAIL",
        "range": f"diagonal B~I_n(k,k), n=3..{max_n}, k>=1 nonempty",
        "even_rows_ok": even_rows_ok,
        "odd_rows_ok": odd_rows_ok,
        "rows": rows,
    }


def a000354_terms(n_max: int) -> list[int]:
    terms = [1]
    for n in range(1, n_max + 1):
        terms.append(2 * n * terms[-1] + ((-1) ** n))
    return terms


def positive_fixed_marginal(n: int) -> dict[int, int]:
    counts: dict[int, int] = defaultdict(int)
    for pi, eps in signed_perms(n):
        counts[positive_identity_count(pi, eps)] += 1
    return dict(counts)


def run_enumeration_boundary(max_bruteforce_n: int = 7) -> dict:
    terms = a000354_terms(8)
    prefix_ok = terms == A000354_PREFIX
    rows = []
    ok = prefix_ok
    for n in range(0, max_bruteforce_n + 1):
        brute = positive_fixed_marginal(n)
        expected = {k: comb(n, k) * terms[n - k] for k in range(n + 1)}
        row_ok = brute == expected
        ok = ok and row_ok
        rows.append(
            {
                "n": n,
                "group_size": sum(brute.values()),
                "distribution": {str(k): brute.get(k, 0) for k in range(n + 1)},
                "expected_distribution": {str(k): expected.get(k, 0) for k in range(n + 1)},
                "m_n_k_equals_binomial_times_A000354": row_ok,
            }
        )
    return {
        "id": "independent_enumeration_boundary_recompute",
        "status": "PASS" if ok else "FAIL",
        "A000354_prefix_checked": terms,
        "A000354_prefix_ok": prefix_ok,
        "bruteforce_range": f"n=0..{max_bruteforce_n}",
        "boundary": "positive-fixed marginal is A000354/type-B derangements; signed two-diagonal count is not A007016",
        "rows": rows,
    }


def square_stats(pi: tuple[int, ...], eps: tuple[int, ...], rho: tuple[int, ...]) -> tuple[int, int, int, int, int]:
    a = b = c = d = 0
    for i in range(len(pi)):
        is_id = pi[i] == i
        is_rev = pi[i] == rho[i]
        if is_id and eps[i] == 1:
            a += 1
        if is_id and eps[i] == -1:
            b += 1
        if is_rev and eps[i] == 1:
            c += 1
        if is_rev and eps[i] == -1:
            d += 1
    return a, b, c, d, negative_sign_count(eps)


def square_eval_brute(n: int, xp: int, xm: int, yp: int, ym: int, q: int) -> int:
    rho = reversal(n)
    total = 0
    for pi, eps in signed_perms(n):
        a, b, c, d, nu = square_stats(pi, eps, rho)
        total += (q ** nu) * (xp ** a) * (xm ** b) * (yp ** c) * (ym ** d)
    return total


def square_cell_value(n: int, i: int, j: int, xp: int, xm: int, yp: int, ym: int, q: int) -> int:
    rho_i = n - 1 - i
    if j == i and j == rho_i:
        return xp * yp + q * xm * ym
    if j == i:
        return xp + q * xm
    if j == rho_i:
        return yp + q * ym
    return 1 + q


def permanent_numeric(matrix: list[list[int]]) -> int:
    n = len(matrix)
    dp = [0 for _ in range(1 << n)]
    dp[0] = 1
    for mask in range(1 << n):
        row = mask.bit_count()
        if row >= n or dp[mask] == 0:
            continue
        subtotal = dp[mask]
        for col in range(n):
            if not (mask & (1 << col)):
                dp[mask | (1 << col)] += subtotal * matrix[row][col]
    return dp[(1 << n) - 1]


def square_eval_permanent(n: int, xp: int, xm: int, yp: int, ym: int, q: int) -> int:
    matrix = [[square_cell_value(n, i, j, xp, xm, yp, ym, q) for j in range(n)] for i in range(n)]
    return permanent_numeric(matrix)


def q_minus_one_closed_value(n: int, xp: int, xm: int, yp: int, ym: int) -> int:
    base = (xp - xm) ** 2 + (yp - ym) ** 2
    if n % 2 == 0:
        return base ** (n // 2)
    return (xp * yp - xm * ym) * (base ** ((n - 1) // 2))


def square_row_counts(n: int) -> dict[tuple[int, int, int, int], int]:
    rho = reversal(n)
    counts: dict[tuple[int, int, int, int], int] = defaultdict(int)
    for pi, eps in signed_perms(n):
        a, b, c, d, _nu = square_stats(pi, eps, rho)
        counts[(a, b, c, d)] += 1
    return dict(counts)


def old_new_marginal_rows(n: int) -> dict:
    square = square_row_counts(n)
    bi_from_square: dict[tuple[int, int], int] = defaultdict(int)
    btilde_from_square: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b, c, d), count in square.items():
        del d
        bi_from_square[(a, b)] += count
        btilde_from_square[(a, c)] += count
    rho = reversal(n)
    bi_direct: dict[tuple[int, int], int] = defaultdict(int)
    btilde_direct: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in signed_perms(n):
        a = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
        b = sum(1 for i in range(n) if pi[i] == i and eps[i] == -1)
        c = sum(1 for i in range(n) if pi[i] == rho[i] and eps[i] == 1)
        bi_direct[(a, b)] += 1
        btilde_direct[(a, c)] += 1
    return {
        "n": n,
        "square_rows": len(square),
        "BI_rows": len(bi_from_square),
        "Btilde_rows": len(btilde_from_square),
        "BI_marginal_ok": dict(bi_from_square) == dict(bi_direct),
        "Btilde_marginal_ok": dict(btilde_from_square) == dict(btilde_direct),
    }


def di_recovery_rows(n: int) -> dict:
    bi: dict[tuple[int, int], int] = defaultdict(int)
    di: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in signed_perms(n):
        k = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
        l = sum(1 for i in range(n) if pi[i] == i and eps[i] == -1)
        bi[(k, l)] += 1
        if negative_sign_count(eps) % 2 == 0:
            di[(k, l)] += 1
    row_records = []
    ok = True
    for k, l in sorted(bi):
        correction = ((-1) ** l) * comb(n, l) if k + l == n else 0
        predicted_twice = bi[(k, l)] + correction
        predicted = predicted_twice // 2
        row_ok = predicted_twice % 2 == 0 and di.get((k, l), 0) == predicted
        ok = ok and row_ok
        row_records.append({"k": k, "l": l, "BI": bi[(k, l)], "DI": di.get((k, l), 0), "predicted_DI": predicted, "ok": row_ok})
    return {"n": n, "rows_checked": len(row_records), "ok": ok, "rows": row_records}


def delta_even(n: int, k: int, l: int) -> int:
    m = n // 2
    return ((-1) ** (k + l)) * sum(
        comb(m, j) * comb(2 * j, k) * comb(2 * m - 2 * j, l)
        for j in range(m + 1)
        if k <= 2 * j and l <= 2 * m - 2 * j
    )


def dtilde_even_rows(n: int) -> dict:
    rho = reversal(n)
    btilde: dict[tuple[int, int], int] = defaultdict(int)
    dtilde: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in signed_perms(n):
        k = positive_identity_count(pi, eps)
        l = positive_reversal_count(pi, eps, rho)
        btilde[(k, l)] += 1
        if negative_sign_count(eps) % 2 == 0:
            dtilde[(k, l)] += 1
    row_records = []
    ok = True
    for k, l in sorted(btilde):
        correction = delta_even(n, k, l)
        predicted_twice = btilde[(k, l)] + correction
        predicted = predicted_twice // 2
        row_ok = predicted_twice % 2 == 0 and dtilde.get((k, l), 0) == predicted
        ok = ok and row_ok
        row_records.append(
            {
                "k": k,
                "l": l,
                "Btilde": btilde[(k, l)],
                "Dtilde": dtilde.get((k, l), 0),
                "Delta": correction,
                "predicted_Dtilde": predicted,
                "ok": row_ok,
            }
        )
    return {"n": n, "rows_checked": len(row_records), "ok": ok, "rows": row_records}


def run_square_and_count_recompute(max_eval_n: int = 7, max_exact_count_n: int = 6) -> dict:
    eval_points = [
        {"xp": 2, "xm": 3, "yp": 5, "ym": 7, "q": 11},
        {"xp": 1, "xm": 4, "yp": 2, "ym": 6, "q": -3},
    ]
    eval_rows = []
    ok = True
    for n in range(2, max_eval_n + 1):
        for point in eval_points:
            brute = square_eval_brute(n, **point)
            per = square_eval_permanent(n, **point)
            row_ok = brute == per
            ok = ok and row_ok
            eval_rows.append({"n": n, "point": point, "brute_eval": brute, "permanent_eval": per, "ok": row_ok})
        q_point = {"xp": 2, "xm": 3, "yp": 5, "ym": 7, "q": -1}
        brute_q = square_eval_brute(n, **q_point)
        closed_q = q_minus_one_closed_value(n, q_point["xp"], q_point["xm"], q_point["yp"], q_point["ym"])
        q_ok = brute_q == closed_q
        ok = ok and q_ok
        eval_rows.append({"n": n, "point": q_point, "q_minus_one_brute": brute_q, "q_minus_one_closed": closed_q, "ok": q_ok})
    marginal_rows = [old_new_marginal_rows(n) for n in range(2, max_exact_count_n + 1)]
    di_rows = [di_recovery_rows(n) for n in range(2, max_exact_count_n + 1)]
    dtilde_rows = [dtilde_even_rows(n) for n in range(2, max_exact_count_n + 1, 2)]
    ok = ok and all(row["BI_marginal_ok"] and row["Btilde_marginal_ok"] for row in marginal_rows)
    ok = ok and all(row["ok"] for row in di_rows)
    ok = ok and all(row["ok"] for row in dtilde_rows)
    return {
        "id": "independent_square_and_count_recompute",
        "status": "PASS" if ok else "FAIL",
        "numeric_square_range": f"n=2..{max_eval_n}, deterministic integer evaluations",
        "exact_count_range": f"n=2..{max_exact_count_n}, exact row formulas",
        "square_eval_rows": eval_rows,
        "marginal_rows": marginal_rows,
        "DI_recovery_rows": di_rows,
        "Dtilde_even_rows": dtilde_rows,
        "boundary": "B^square is count grammar; D~I remains scout-only and no Type-D theorem is promoted",
    }


def active_boundary_files(root: Path) -> list[Path]:
    return [
        root / "STATUS.md",
        root / "ROUTE_PLAN.md",
        root / "CLAIM_LEDGER.md",
        root / "PROOF_OBLIGATIONS.md",
        root / "PUBLIC_CLAIM_BOUNDARY.md",
        root / "README.md",
        root / "P15_DETAILED_ROADMAP_2026_07_06.md",
        root / "docs" / "P15_S1_OBJECT_ADMISSION_AND_CONVENTIONS_2026_07_06.md",
        root / "docs" / "P15_S4_EVEN_N_THEOREM_2026_07_06.md",
        root / "docs" / "P15_S5_ODD_N_STRUCTURE_2026_07_06.md",
        root / "docs" / "P15_S6_SQUARE_AND_CLOSED_FORMS_2026_07_06.md",
    ]


def run_boundary_scan() -> dict:
    root = Path(__file__).resolve().parents[1]
    files = active_boundary_files(root)
    file_records = []
    combined_parts = []
    forbidden_hits = []
    forbidden_needles = [
        "<=> H_m",
        "<=>  H_m",
        "a(2m+1,1) > 0  <=>",
        "closed modulo a routine holonomic certificate",
    ]
    for path in files:
        text = path.read_text(encoding="utf-8")
        combined_parts.append(text)
        hits = [needle for needle in forbidden_needles if needle in text]
        if hits:
            forbidden_hits.append({"path": str(path.relative_to(root)), "hits": hits})
        file_records.append({"path": str(path.relative_to(root)), "bytes": len(text.encode("utf-8")), "forbidden_hits": hits})
    combined = "\n".join(combined_parts)
    lower = combined.lower()
    required = {
        "not_A007016_for_signed_enumeration": "not a007016" in lower or "retracted for the signed context" in lower,
        "Dtilde_scout_only": "d~i" in lower and ("scout-only" in lower or "scout only" in lower),
        "no_public_theorem_claim": "no public theorem claim" in lower or "no public claim" in lower,
        "matching_scheme_bookkeeping": "matching-scheme bookkeeping" in lower,
        "H_m_left_to_S8": "h_m" in lower and "s8" in lower and ("residual" in lower or "obligation" in lower),
    }
    ok = not forbidden_hits and all(required.values())
    return {
        "id": "active_boundary_text_scan",
        "status": "PASS" if ok else "FAIL",
        "files_scanned": file_records,
        "forbidden_hits": forbidden_hits,
        "required_boundary_phrases": required,
        "scope_note": "Scans active governance/docs surfaces, not historical review notes from pre-S2 investigation.",
    }


def make_report(max_channel_n: int, max_eval_n: int, max_exact_count_n: int, max_enum_n: int) -> dict:
    checks = [
        run_channel_recompute(max_channel_n),
        run_enumeration_boundary(max_enum_n),
        run_square_and_count_recompute(max_eval_n, max_exact_count_n),
        run_boundary_scan(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    proof_checklist = {
        "no_import_from_S3_S6_certificate_scripts": True,
        "ambient_channels_recomputed_from_raw_signed_permutations": checks[0]["status"] == "PASS",
        "A000354_boundary_recomputed": checks[1]["status"] == "PASS",
        "square_and_D_scout_counts_recomputed_independently": checks[2]["status"] == "PASS",
        "active_boundary_scan_clean": checks[3]["status"] == "PASS",
        "no_public_theorem_promoted_at_S7": True,
        "S8_H_m_and_odd_pair_positivity_open_at_S7": True,
    }
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S7 independent red-team recompute",
        "status": "PASS" if not failed and all(proof_checklist.values()) else "FAIL",
        "decision": "S7 red-team recompute passes; next gate is S8 H_m recurrence certificate." if not failed else "S7 red-team recompute failed; inspect failed_checks.",
        "failed_checks": failed,
        "checks": checks,
        "proof_checklist": proof_checklist,
        "closed_at_S7": [
            "A second verifier independently recomputes the S3-S6 finite channel/count surface without importing those scripts.",
            "Even rows n=4,6 are rechecked against the S4 rank/collapse theorem surface.",
            "Odd rows n=3,5,7 are rechecked against the S5 rank/collapse reduction surface, including the n=3,k=1 base marker.",
            "The square parent grammar is red-teamed by deterministic permanent evaluations and exact DI/D~I row formulas.",
            "Active boundary text is clean of the stale H_m equivalence wording.",
        ],
        "open_after_S7": [
            "S8 H_m recurrence / creative-telescoping certificate for odd k=1 was open at S7. [Subsequently closed by S8 and S11A-H.]",
            "General odd pair-standard positivity for lam_prime and mu was open at S7. [Subsequently closed by S9C and S11A-K.]",
            "P15-S9 manuscript/public-claim shaping was open at S7. [Subsequently completed: S9/S10 public boundary and the S12 separator theorem.]",
            "D~I remains Type-D scout-only; no Type-D theorem is promoted.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S7 independent red-team recompute")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        if check["id"] == "independent_channel_rank_and_collapse_recompute":
            print("  range:", check["range"])
            for row in check["rows"]:
                print(
                    "  n={n} k={k} size={size} ranks=({rank_M_ref},{rank_M_pair_std}) expected=({expected_rank_M_ref},{expected_rank_M_pair_std}) collapse={collapse_ok} ok={ok}".format(**row)
                )
        elif check["id"] == "independent_enumeration_boundary_recompute":
            print("  A000354_prefix_ok:", check["A000354_prefix_ok"])
            print("  brute:", check["bruteforce_range"])
        elif check["id"] == "independent_square_and_count_recompute":
            print("  numeric_square_range:", check["numeric_square_range"])
            print("  exact_count_range:", check["exact_count_range"])
            print("  eval_rows:", len(check["square_eval_rows"]))
            print("  marginal_rows_ok:", all(row["BI_marginal_ok"] and row["Btilde_marginal_ok"] for row in check["marginal_rows"]))
            print("  DI_rows_ok:", all(row["ok"] for row in check["DI_recovery_rows"]))
            print("  Dtilde_rows_ok:", all(row["ok"] for row in check["Dtilde_even_rows"]))
        elif check["id"] == "active_boundary_text_scan":
            print("  files_scanned:", len(check["files_scanned"]))
            print("  forbidden_hits:", check["forbidden_hits"])
            print("  required:", check["required_boundary_phrases"])


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S7 independent red-team recompute")
    parser.add_argument("--write-json", help="write deterministic JSON certificate to this path")
    parser.add_argument("--max-channel-n", type=int, default=7)
    parser.add_argument("--max-eval-n", type=int, default=7)
    parser.add_argument("--max-exact-count-n", type=int, default=6)
    parser.add_argument("--max-enum-n", type=int, default=7)
    args = parser.parse_args()
    report = make_report(args.max_channel_n, args.max_eval_n, args.max_exact_count_n, args.max_enum_n)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
