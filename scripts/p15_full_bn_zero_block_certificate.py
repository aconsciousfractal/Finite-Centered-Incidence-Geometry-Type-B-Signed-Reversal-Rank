#!/usr/bin/env python3
"""
Exact finite certificate for one native full-B_n zero block in P15.

Target:
    M_{(1.1 | 3.2.1)}(B~I_8(3,3)) = 0.

This is not a theorem about all native B_n irreps.  It is a finite witness
that the naive full-fingerprint R-even saturation conjecture fails at n=8.

Method:
    Use the same little-group model as the modular full-fingerprint scout.
    Here lambda_plus=(1,1), so S^{lambda_plus} is the sign representation
    of S_2.  The S^{(3,2,1)} seminormal matrices are computed exactly over Q;
    all 720 permutation matrices have denominators dividing 96.  The layer
    operator is summed after multiplying those S_6 matrices by 96.  A zero
    scaled integer matrix proves the rational matrix is zero.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from collections import defaultdict
from dataclasses import dataclass, asdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

Perm = Tuple[int, ...]
Mask = int
Pair = Tuple[int, int]


def gcd_lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


@lru_cache(None)
def standard_tableaux(shape: Tuple[int, ...]) -> Tuple[Tuple[Tuple[int, int], ...], ...]:
    n = sum(shape)
    if n == 0:
        return ((),)
    corners = []
    for r, rowlen in enumerate(shape):
        if rowlen > 0 and (r == len(shape) - 1 or shape[r + 1] < rowlen):
            corners.append((r, rowlen - 1))
    out = []
    for cell in corners:
        r, _ = cell
        new_shape = list(shape)
        new_shape[r] -= 1
        while new_shape and new_shape[-1] == 0:
            new_shape.pop()
        for tab in standard_tableaux(tuple(new_shape)):
            out.append(tab + (cell,))
    return tuple(out)


def identity_fraction_matrix(d: int) -> List[List[Fraction]]:
    return [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]


def fraction_matmul(A: List[List[Fraction]], B: List[List[Fraction]]) -> List[List[Fraction]]:
    n = len(A)
    k = len(B)
    m = len(B[0])
    Z = [[Fraction(0) for _ in range(m)] for __ in range(n)]
    for i in range(n):
        Zi = Z[i]
        for t in range(k):
            a = A[i][t]
            if not a:
                continue
            Bt = B[t]
            for j in range(m):
                if Bt[j]:
                    Zi[j] += a * Bt[j]
    return Z


def perm_word_right(perm: Perm) -> List[int]:
    arr = list(range(len(perm)))
    word: List[int] = []
    for pos, val in enumerate(perm):
        loc = arr.index(val)
        while loc > pos:
            arr[loc - 1], arr[loc] = arr[loc], arr[loc - 1]
            word.append(loc - 1)
            loc -= 1
    if tuple(arr) != tuple(perm):
        raise AssertionError((arr, perm))
    return word


class SpechtSeminormalExact:
    def __init__(self, shape: Tuple[int, ...]):
        self.shape = tuple(shape)
        self.n = sum(shape)
        self.tabs = list(standard_tableaux(self.shape))
        self.dim = len(self.tabs)
        self.index = {t: i for i, t in enumerate(self.tabs)}
        self.gens = self._build_adjacent_generators()
        self.perm_cache: Dict[Perm, List[List[Fraction]]] = {}

    def _build_adjacent_generators(self) -> List[List[List[Fraction]]]:
        out = []
        for i in range(self.n - 1):
            M = [[Fraction(0) for _ in range(self.dim)] for __ in range(self.dim)]
            for col, T in enumerate(self.tabs):
                r1, c1 = T[i]
                r2, c2 = T[i + 1]
                axial = (c2 - r2) - (c1 - r1)
                if axial == 1:
                    M[col][col] += 1
                elif axial == -1:
                    M[col][col] -= 1
                else:
                    invr = Fraction(1, axial)
                    L = list(T)
                    L[i], L[i + 1] = L[i + 1], L[i]
                    rowp = self.index[tuple(L)]
                    M[col][col] += invr
                    M[rowp][col] += 1 + invr
            out.append(M)
        return out

    def mat_perm(self, perm: Perm) -> List[List[Fraction]]:
        perm = tuple(perm)
        cached = self.perm_cache.get(perm)
        if cached is not None:
            return cached
        M = identity_fraction_matrix(self.dim)
        for s in perm_word_right(perm):
            M = fraction_matmul(M, self.gens[s])
        self.perm_cache[perm] = M
        return M

    def relation_checks(self) -> Dict[str, bool]:
        I = identity_fraction_matrix(self.dim)
        involution = True
        braid = True
        commute = True
        for i, G in enumerate(self.gens):
            involution = involution and (fraction_matmul(G, G) == I)
            if i + 1 < len(self.gens):
                left = fraction_matmul(fraction_matmul(G, self.gens[i + 1]), G)
                right = fraction_matmul(fraction_matmul(self.gens[i + 1], G), self.gens[i + 1])
                braid = braid and (left == right)
            for j in range(i + 2, len(self.gens)):
                left = fraction_matmul(G, self.gens[j])
                right = fraction_matmul(self.gens[j], G)
                commute = commute and (left == right)
        return {"s_i_squared": involution, "braid": braid, "distant_commute": commute}


def scaled_permutation_matrices(shape: Tuple[int, ...]) -> Tuple[int, Dict[Perm, np.ndarray], Dict[str, bool]]:
    S = SpechtSeminormalExact(shape)
    checks = S.relation_checks()
    perms = list(itertools.permutations(range(sum(shape))))
    denom = 1
    mats_q: Dict[Perm, List[List[Fraction]]] = {}
    for perm in perms:
        M = S.mat_perm(perm)
        mats_q[perm] = M
        for row in M:
            for x in row:
                denom = gcd_lcm(denom, x.denominator)
    mats_z: Dict[Perm, np.ndarray] = {}
    for perm, M in mats_q.items():
        Z = np.zeros((S.dim, S.dim), dtype=np.int64)
        for i, row in enumerate(M):
            for j, x in enumerate(row):
                y = x * denom
                if y.denominator != 1:
                    raise AssertionError((perm, i, j, x, denom))
                Z[i, j] = int(y.numerator)
        mats_z[perm] = Z
    return denom, mats_z, checks


@lru_cache(None)
def positions_of_mask(n: int, mask: Mask) -> Tuple[int, ...]:
    return tuple(i for i in range(n) if (mask >> i) & 1)


def masks_of_size(n: int, size: int) -> List[Mask]:
    out = []
    for comb in itertools.combinations(range(n), size):
        mask = 0
        for i in comb:
            mask |= 1 << i
        out.append(mask)
    return out


def rho_perm(n: int) -> Perm:
    return tuple(n - 1 - i for i in range(n))


def restricted_perm_data(pi: Perm, mask: Mask, n: int) -> Tuple[Mask, Perm]:
    dom = positions_of_mask(n, mask)
    imgs = tuple(pi[i] for i in dom)
    image_mask = 0
    for v in imgs:
        image_mask |= 1 << v
    image_pos = positions_of_mask(n, image_mask)
    rank = {v: j for j, v in enumerate(image_pos)}
    sigma = tuple(rank[v] for v in imgs)
    return image_mask, sigma


def sign_of_perm(perm: Perm) -> int:
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            inv += int(perm[i] > perm[j])
    return -1 if inv % 2 else 1


def sign_poly_target_coeff(pi: Perm, char_mask: Mask, target: Pair) -> int:
    n = len(pi)
    rho = rho_perm(n)
    poly: Dict[Pair, int] = {(0, 0): 1}
    k_target, l_target = target
    for i in range(n):
        alpha = 1 if pi[i] == i else 0
        beta = 1 if pi[i] == rho[i] else 0
        minus_weight = -1 if ((char_mask >> i) & 1) else 1
        new: Dict[Pair, int] = defaultdict(int)
        for (a, c), val in poly.items():
            ap = a + alpha
            cp = c + beta
            if ap <= k_target and cp <= l_target:
                new[(ap, cp)] += val
            new[(a, c)] += val * minus_weight
        poly = dict(new)
    return poly.get(target, 0)


def layer_size(n: int, target: Pair) -> int:
    return sum(sign_poly_target_coeff(pi, 0, target) for pi in itertools.permutations(range(n)))


@dataclass(frozen=True)
class CertificateSummary:
    status: str
    claim: str
    n: int
    k: int
    l: int
    lambda_plus: str
    lambda_minus: str
    dimension: int
    subset_count: int
    specht_plus_dim: int
    specht_minus_dim: int
    specht_minus_common_denominator: int
    layer_size: int
    permutations_scanned: int
    source_subsets_per_permutation: int
    active_terms: int
    sum_abs_coefficients: int
    max_abs_scaled_entry: int
    nonzero_scaled_entries: int
    relation_checks_pass: bool
    relation_checks: Dict[str, bool]
    runtime_seconds: float
    boundary: str


def write_markdown(summary: CertificateSummary, path: Path) -> None:
    s = asdict(summary)
    lines = [
        "# P15 Full Bn Zero Block Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Method",
        "",
        "The certificate uses the same little-group model as the modular full-fingerprint scout.",
        "For `lambda_plus=(1,1)`, the plus Specht module is the one-dimensional sign representation of `S_2`.",
        "For `lambda_minus=(3,2,1)`, all Young seminormal `S_6` permutation matrices are computed exactly over `Q`.",
        f"Their common denominator is `{summary.specht_minus_common_denominator}`; the layer operator is summed after clearing this denominator.",
        "If the cleared integer matrix is zero, the original rational Fourier block is zero.",
        "",
        "## Exact Audit",
        "",
        "| field | value |",
        "|---|---:|",
    ]
    for key in [
        "dimension",
        "subset_count",
        "specht_minus_dim",
        "layer_size",
        "permutations_scanned",
        "source_subsets_per_permutation",
        "active_terms",
        "sum_abs_coefficients",
        "max_abs_scaled_entry",
        "nonzero_scaled_entries",
    ]:
        lines.append(f"| `{key}` | {s[key]} |")
    lines += [
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "This is a finite exact witness against the naive full native `B_n` saturation conjecture. It does not promote a full-fingerprint theorem and does not change the scoped two-standard-channel P15 theorem.",
        "",
        "Prossima task: use the family certificate for the full `lambda_plus=(1,1), lambda_minus |- 6` zero family, then seek a conceptual explanation or return to final PDF readthrough.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run_certificate(out_json: Path | None, out_md: Path | None) -> CertificateSummary:
    t0 = time.time()
    n = 8
    target = (3, 3)
    lam_plus = (1, 1)
    lam_minus = (3, 2, 1)
    a = sum(lam_plus)
    b = sum(lam_minus)
    if a + b != n:
        raise AssertionError("bad bipartition")

    denom, minus_mats, checks = scaled_permutation_matrices(lam_minus)
    relation_checks_pass = all(checks.values())
    if not relation_checks_pass:
        raise AssertionError(checks)

    subset_masks = masks_of_size(n, a)
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full_mask = (1 << n) - 1
    block_dim = len(next(iter(minus_mats.values())))
    dim = len(subset_masks) * block_dim
    agg = np.zeros((len(subset_masks), len(subset_masks), block_dim, block_dim), dtype=np.int64)

    active_terms = 0
    sum_abs_coefficients = 0
    perms_scanned = 0
    restrict_cache: Dict[Tuple[Perm, Mask], Tuple[Mask, Perm]] = {}

    for pi in itertools.permutations(range(n)):
        perms_scanned += 1
        pi = tuple(pi)
        for src_si, A_mask in enumerate(subset_masks):
            B_mask = full_mask ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            active_terms += 1
            sum_abs_coefficients += abs(coeff)
            key_a = (pi, A_mask)
            data_a = restrict_cache.get(key_a)
            if data_a is None:
                data_a = restricted_perm_data(pi, A_mask, n)
                restrict_cache[key_a] = data_a
            Ap_mask, sigmaA = data_a
            key_b = (pi, B_mask)
            data_b = restrict_cache.get(key_b)
            if data_b is None:
                data_b = restricted_perm_data(pi, B_mask, n)
                restrict_cache[key_b] = data_b
            _Bp_mask, sigmaB = data_b
            dst_si = subset_index[Ap_mask]
            signA = sign_of_perm(sigmaA)
            agg[dst_si, src_si] += coeff * signA * minus_mats[sigmaB]

    max_abs_scaled_entry = int(np.max(np.abs(agg))) if agg.size else 0
    nonzero_scaled_entries = int(np.count_nonzero(agg))
    status = "PASS" if max_abs_scaled_entry == 0 and nonzero_scaled_entries == 0 else "FAIL"
    summary = CertificateSummary(
        status=status,
        claim="M_{(1.1 | 3.2.1)}(B~I_8(3,3)) = 0 in the rational little-group model",
        n=n,
        k=target[0],
        l=target[1],
        lambda_plus="1.1",
        lambda_minus="3.2.1",
        dimension=dim,
        subset_count=len(subset_masks),
        specht_plus_dim=1,
        specht_minus_dim=block_dim,
        specht_minus_common_denominator=denom,
        layer_size=layer_size(n, target),
        permutations_scanned=perms_scanned,
        source_subsets_per_permutation=len(subset_masks),
        active_terms=active_terms,
        sum_abs_coefficients=sum_abs_coefficients,
        max_abs_scaled_entry=max_abs_scaled_entry,
        nonzero_scaled_entries=nonzero_scaled_entries,
        relation_checks_pass=relation_checks_pass,
        relation_checks=checks,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite exact full-fingerprint defect witness only; not a public full-B_n fingerprint theorem.",
    )
    if out_json is not None:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md is not None:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = run_certificate(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
