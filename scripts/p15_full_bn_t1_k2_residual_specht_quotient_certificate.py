#!/usr/bin/env python3
"""
Residual Specht-quotient certificate for T_1^(8,2) in the P15 n=8 atlas.

This closes the three residual nonzero lambda_plus=1 extra-defect mechanisms
left after the extended quotient/sign-twist package:

    (4,1,1,1) by oriented 3-subsets modulo oriented-pair boundaries;
    (3,3,1) by the antisymmetric two-triples-plus-point module modulo the
              central projector onto the lower Specht types;
    (3,2,2) by conjugate/sign-twist rank transfer.

This is a finite n=8,k=2 certificate.  It is not an all-n theorem and not a
full native B_n fingerprint classification.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from p15_full_bn_prespecht_fixed_plus_scout import (
    block_matrix_mod,
    build_prespecht_terms,
    mod_rank,
    needed_specht_matrices,
)
from p15_full_bn_zero_block_certificate import SpechtSeminormalExact

PRIME = 1000003
ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class OrientedTripleProfile:
    target_shape: str
    natural_dim: int
    lower_rank: int
    natural_operator_rank: int
    image_mod_lower_rank: int
    image_intersection_lower_rank: int
    source_lower_image_rank: int
    source_lower_image_outside_lower_rank: int
    specht_dim: int
    specht_block_rank: int
    pass_conditions: Dict[str, bool]


@dataclass(frozen=True)
class AntiTripleProjectorProfile:
    target_shape: str
    conjugate_partner: str
    fiber_dim: int
    projector_lower_shapes: List[str]
    projector_rank_per_fiber: int
    projector_idempotent_residual_rank: int
    full_lower_rank: int
    natural_operator_rank: int
    image_mod_lower_rank: int
    image_intersection_lower_rank: int
    source_lower_image_rank: int
    source_lower_image_outside_lower_rank: int
    specht_dim: int
    specht_block_rank: int
    conjugate_specht_dim: int
    conjugate_specht_block_rank: int
    pass_conditions: Dict[str, bool]


@dataclass(frozen=True)
class ResidualSummary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    prime: int
    subset_count: int
    group_keys: int
    oriented_triple: OrientedTripleProfile
    anti_triple_projector: AntiTripleProjectorProfile
    runtime_seconds: float
    boundary: str
    next_task: str


def part_str(part: Tuple[int, ...]) -> str:
    return ".".join(str(x) for x in part) if part else "empty"


def complement_points(n: int, plus_mask: int) -> List[int]:
    return [x for x in range(n) if not ((plus_mask >> x) & 1)]


def canonical_sorted_with_sign(values: List[int]) -> Tuple[Tuple[int, ...], int]:
    inv = 0
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            inv += int(values[i] > values[j])
    return tuple(sorted(values)), (-1 if inv % 2 else 1)


def canonical_pair_of_triples(A, B) -> Tuple[Tuple[int, ...], Tuple[int, ...], int]:
    A = tuple(sorted(A))
    B = tuple(sorted(B))
    if A < B:
        return A, B, 1
    return B, A, -1


def specht_rank(shape: Tuple[int, ...], subset_count: int, grouped, p: int) -> Tuple[int, int]:
    sigmas = {sigma for (_dst, _src, sigma) in grouped}
    mats, dim = needed_specht_matrices(shape, sigmas, p)
    block = block_matrix_mod(grouped, subset_count, mats, dim, p)
    return dim, mod_rank(block, p)


def oriented_triple_basis(n: int, plus_masks: List[int]):
    basis: List[Tuple[int, Tuple[int, int, int]]] = []
    index: Dict[Tuple[int, Tuple[int, int, int]], int] = {}
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for triple in itertools.combinations(comp, 3):
            key = (plus_i, tuple(triple))
            index[key] = len(basis)
            basis.append(key)
    return basis, index


def oriented_triple_operator(n: int, plus_masks: List[int], grouped, p: int) -> np.ndarray:
    basis, index = oriented_triple_basis(n, plus_masks)
    M = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for (dst, src, sigma), coeff in grouped.items():
        src_comp = complement_points(n, plus_masks[src])
        dst_comp = complement_points(n, plus_masks[dst])
        for triple in itertools.combinations(src_comp, 3):
            ranks = [src_comp.index(x) for x in triple]
            image, sign = canonical_sorted_with_sign([dst_comp[sigma[t]] for t in ranks])
            M[index[(dst, image)], index[(src, tuple(triple))]] = (
                M[index[(dst, image)], index[(src, tuple(triple))]] + coeff * sign
            ) % p
    return M


def oriented_triple_lower(n: int, plus_masks: List[int], p: int) -> np.ndarray:
    basis, index = oriented_triple_basis(n, plus_masks)
    columns: List[np.ndarray] = []
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for a, b in itertools.combinations(comp, 2):
            v = np.zeros((len(basis),), dtype=np.int64)
            for z in comp:
                if z == a or z == b:
                    continue
                triple, sign = canonical_sorted_with_sign([a, b, z])
                v[index[(plus_i, triple)]] = (v[index[(plus_i, triple)]] + sign) % p
            columns.append(v)
    return np.stack(columns, axis=1)


def oriented_triple_profile(n: int, plus_masks: List[int], grouped, p: int) -> OrientedTripleProfile:
    M = oriented_triple_operator(n, plus_masks, grouped, p)
    B = oriented_triple_lower(n, plus_masks, p)
    rank_M = mod_rank(M, p)
    rank_B = mod_rank(B, p)
    concat_rank = mod_rank(np.hstack([B, M]), p)
    MB = (M @ B) % p
    source_lower_rank = mod_rank(MB, p)
    source_lower_outside = mod_rank(np.hstack([B, MB]), p) - rank_B
    dim, sprank = specht_rank((4, 1, 1, 1), len(plus_masks), grouped, p)
    image_mod_lower = concat_rank - rank_B
    intersection = rank_M + rank_B - concat_rank
    checks = {
        "natural_rank_is_120": rank_M == 120,
        "lower_rank_is_120": rank_B == 120,
        "image_mod_lower_is_72": image_mod_lower == 72,
        "source_lower_image_rank_is_48": source_lower_rank == 48,
        "source_lower_maps_inside_lower": source_lower_outside == 0,
        "specht_rank_matches_quotient_rank": sprank == image_mod_lower,
    }
    return OrientedTripleProfile(
        target_shape="4.1.1.1",
        natural_dim=M.shape[0],
        lower_rank=rank_B,
        natural_operator_rank=rank_M,
        image_mod_lower_rank=image_mod_lower,
        image_intersection_lower_rank=intersection,
        source_lower_image_rank=source_lower_rank,
        source_lower_image_outside_lower_rank=source_lower_outside,
        specht_dim=dim,
        specht_block_rank=sprank,
        pass_conditions=checks,
    )


def anti331_fiber_basis():
    n = 7
    basis: List[Tuple[Tuple[int, ...], Tuple[int, ...], int]] = []
    index: Dict[Tuple[Tuple[int, ...], Tuple[int, ...], int], int] = {}
    for z in range(n):
        rem = [x for x in range(n) if x != z]
        for A in itertools.combinations(rem, 3):
            B = tuple(x for x in rem if x not in A)
            C, D, _sign = canonical_pair_of_triples(A, B)
            key = (C, D, z)
            if key not in index:
                index[key] = len(basis)
                basis.append(key)
    return basis, index


def cycle_type(perm: Tuple[int, ...]) -> Tuple[int, ...]:
    seen = [False] * len(perm)
    out: List[int] = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        j = i
        length = 0
        while not seen[j]:
            seen[j] = True
            length += 1
            j = perm[j]
        out.append(length)
    return tuple(sorted(out, reverse=True))


def partitions(n: int, max_part: int | None = None):
    if n == 0:
        yield ()
        return
    if max_part is None or max_part > n:
        max_part = n
    for first in range(max_part, 0, -1):
        rest_n = n - first
        for rest in partitions(rest_n, min(first, rest_n) if rest_n else 0):
            yield (first,) + rest


def perm_from_cycle_type(ct: Tuple[int, ...]) -> Tuple[int, ...]:
    out = list(range(sum(ct)))
    start = 0
    for length in ct:
        cyc = list(range(start, start + length))
        start += length
        for i, x in enumerate(cyc):
            out[x] = cyc[(i + 1) % length]
    return tuple(out)


def anti331_lower_projector(p: int) -> Tuple[np.ndarray, int]:
    basis, index = anti331_fiber_basis()
    lower_shapes = [(6, 1), (5, 2), (5, 1, 1), (4, 3)]
    coeff_by_class: Dict[Tuple[int, ...], int] = {}
    for ct in partitions(7):
        rep = perm_from_cycle_type(ct)
        coeff = 0
        for shape in lower_shapes:
            S = SpechtSeminormalExact(shape)
            mat = S.mat_perm(rep)
            chi = sum(mat[i][i] for i in range(S.dim))
            coeff += S.dim * int(chi)
        coeff_by_class[ct] = (coeff % p) * pow(math.factorial(7), p - 2, p) % p

    E = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for perm in itertools.permutations(range(7)):
        coeff = coeff_by_class[cycle_type(tuple(perm))]
        if coeff == 0:
            continue
        for col, (A, B, z) in enumerate(basis):
            image_A = tuple(sorted(perm[x] for x in A))
            image_B = tuple(sorted(perm[x] for x in B))
            image_z = perm[z]
            C, D, sign = canonical_pair_of_triples(image_A, image_B)
            row = index[(C, D, image_z)]
            E[row, col] = (E[row, col] + coeff * sign) % p
    idempotent_residual_rank = mod_rank((E @ E - E) % p, p)
    return E, idempotent_residual_rank


def anti331_basis(n: int, plus_masks: List[int]):
    basis: List[Tuple[int, Tuple[int, ...], Tuple[int, ...], int]] = []
    index: Dict[Tuple[int, Tuple[int, ...], Tuple[int, ...], int], int] = {}
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for z in comp:
            rem = [x for x in comp if x != z]
            for A in itertools.combinations(rem, 3):
                B = tuple(x for x in rem if x not in A)
                C, D, _sign = canonical_pair_of_triples(A, B)
                key = (plus_i, C, D, z)
                if key not in index:
                    index[key] = len(basis)
                    basis.append(key)
    return basis, index


def anti331_operator(n: int, plus_masks: List[int], grouped, p: int) -> np.ndarray:
    basis, index = anti331_basis(n, plus_masks)
    M = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for (dst, src, sigma), coeff in grouped.items():
        src_comp = complement_points(n, plus_masks[src])
        dst_comp = complement_points(n, plus_masks[dst])
        for _plus_i, A, B, z in [row for row in basis if row[0] == src]:
            image_A = tuple(sorted(dst_comp[sigma[src_comp.index(x)]] for x in A))
            image_B = tuple(sorted(dst_comp[sigma[src_comp.index(x)]] for x in B))
            image_z = dst_comp[sigma[src_comp.index(z)]]
            C, D, sign = canonical_pair_of_triples(image_A, image_B)
            M[index[(dst, C, D, image_z)], index[(src, A, B, z)]] = (
                M[index[(dst, C, D, image_z)], index[(src, A, B, z)]] + coeff * sign
            ) % p
    return M


def block_diagonal_projector(E: np.ndarray, blocks: int) -> np.ndarray:
    dim = E.shape[0]
    out = np.zeros((blocks * dim, blocks * dim), dtype=np.int64)
    for b in range(blocks):
        out[b * dim:(b + 1) * dim, b * dim:(b + 1) * dim] = E
    return out


def anti331_profile(n: int, plus_masks: List[int], grouped, p: int) -> AntiTripleProjectorProfile:
    E, idem_residual = anti331_lower_projector(p)
    M = anti331_operator(n, plus_masks, grouped, p)
    B = block_diagonal_projector(E, len(plus_masks))
    rank_M = mod_rank(M, p)
    rank_B = mod_rank(B, p)
    concat_rank = mod_rank(np.hstack([B, M]), p)
    MB = (M @ B) % p
    source_lower_rank = mod_rank(MB, p)
    source_lower_outside = mod_rank(np.hstack([B, MB]), p) - rank_B
    dim, sprank = specht_rank((3, 3, 1), len(plus_masks), grouped, p)
    cdim, crank = specht_rank((3, 2, 2), len(plus_masks), grouped, p)
    image_mod_lower = concat_rank - rank_B
    intersection = rank_M + rank_B - concat_rank
    checks = {
        "projector_rank_per_fiber_is_49": mod_rank(E, p) == 49,
        "projector_is_idempotent": idem_residual == 0,
        "full_lower_rank_is_392": rank_B == 392,
        "natural_rank_is_192": rank_M == 192,
        "image_mod_lower_is_60": image_mod_lower == 60,
        "source_lower_image_rank_is_132": source_lower_rank == 132,
        "source_lower_maps_inside_lower": source_lower_outside == 0,
        "specht_331_rank_matches_quotient_rank": sprank == image_mod_lower,
        "conjugate_322_rank_matches_331": crank == sprank,
    }
    return AntiTripleProjectorProfile(
        target_shape="3.3.1",
        conjugate_partner="3.2.2",
        fiber_dim=E.shape[0],
        projector_lower_shapes=["6.1", "5.2", "5.1.1", "4.3"],
        projector_rank_per_fiber=mod_rank(E, p),
        projector_idempotent_residual_rank=idem_residual,
        full_lower_rank=rank_B,
        natural_operator_rank=rank_M,
        image_mod_lower_rank=image_mod_lower,
        image_intersection_lower_rank=intersection,
        source_lower_image_rank=source_lower_rank,
        source_lower_image_outside_lower_rank=source_lower_outside,
        specht_dim=dim,
        specht_block_rank=sprank,
        conjugate_specht_dim=cdim,
        conjugate_specht_block_rank=crank,
        pass_conditions=checks,
    )


def certify(out_json: Path | None, out_md: Path | None) -> ResidualSummary:
    t0 = time.time()
    n = 8
    k = 2
    plus_masks, grouped, _perms, _raw_active, _raw_abs, _zero_removed = build_prespecht_terms(n, k, 1)
    oriented = oriented_triple_profile(n, plus_masks, grouped, PRIME)
    anti = anti331_profile(n, plus_masks, grouped, PRIME)
    status = "PASS" if all(oriented.pass_conditions.values()) and all(anti.pass_conditions.values()) else "FAIL"
    summary = ResidualSummary(
        status=status,
        claim=(
            "The three residual nonzero lambda_plus=1 extra-defect mechanisms of T_1^(8,2) are now accounted for: "
            "(4,1,1,1) is the top quotient of oriented 3-subsets modulo oriented-pair boundaries; "
            "(3,3,1) is the top quotient of the antisymmetric two-triples-plus-point module modulo the central "
            "projector onto lower Specht types; and (3,2,2) follows by sign-twist conjugacy."
        ),
        n=n,
        k=k,
        lambda_plus="1",
        complement_size=7,
        prime=PRIME,
        subset_count=len(plus_masks),
        group_keys=len(grouped),
        oriented_triple=oriented,
        anti_triple_projector=anti,
        runtime_seconds=round(time.time() - t0, 3),
        boundary=(
            "Finite n=8,k=2 residual mechanism certificate only. This completes the lambda_plus=1 fixed-plus "
            "family accounting at n=8, but it does not prove an all-n theorem and does not classify other fixed-plus families."
        ),
        next_task="Keep this package appendix-scoped; optional rational audits and any remaining families are separate finite checks, not a classification theorem.",
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: ResidualSummary, path: Path) -> None:
    o = summary.oriented_triple
    a = summary.anti_triple_projector
    lines = [
        "# P15 T1 k=2 Residual Specht-Quotient Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Residual Mechanisms",
        "",
        "| residual shape | mechanism | natural rank | lower rank | quotient rank | Specht rank |",
        "|---|---|---:|---:|---:|---:|",
        f"| `{o.target_shape}` | oriented 3-subsets modulo oriented-pair boundaries | {o.natural_operator_rank} | {o.lower_rank} | {o.image_mod_lower_rank} | {o.specht_block_rank} |",
        f"| `{a.target_shape}` | antisymmetric two-triples-plus-point modulo lower central projector | {a.natural_operator_rank} | {a.full_lower_rank} | {a.image_mod_lower_rank} | {a.specht_block_rank} |",
        f"| `{a.conjugate_partner}` | sign-twist conjugate of `{a.target_shape}` | - | - | {a.image_mod_lower_rank} | {a.conjugate_specht_block_rank} |",
        "",
        "## Details",
        "",
        "For `(4,1,1,1)`, the oriented 3-subset module has dimension 280 across the eight fibers. The oriented-pair boundary subspace has rank 120, the operator maps it into itself with image rank 48, and the quotient image has rank 72, matching the direct Specht block.",
        "",
        "For `(3,3,1)`, use the antisymmetric module of two complementary triples plus a singleton. One fiber has dimension 70 and decomposes into lower Specht types plus the top `(3,3,1)`. The central lower projector onto `(6,1)`, `(5,2)`, `(5,1,1)`, and `(4,3)` has rank 49 per fiber and is idempotent. Across eight fibers this gives lower rank 392; the quotient image has rank 60, matching the direct Specht block `(3,3,1)`. The conjugate partner `(3,2,2)` has the same rank by the previously certified sign-twist row/column equivalence.",
        "",
        "## Checks",
        "",
        "### `(4,1,1,1)`",
        "",
        "```json",
        json.dumps(o.pass_conditions, indent=2),
        "```",
        "",
        "### `(3,3,1)` / `(3,2,2)`",
        "",
        "```json",
        json.dumps(a.pass_conditions, indent=2),
        "```",
        "",
        "## Sources",
        "",
        "```text",
        "scripts/p15_full_bn_t1_k2_residual_specht_quotient_certificate.py",
        "certified/P15_FULL_BN_T1_K2_RESIDUAL_SPECHT_QUOTIENT_CERTIFICATE.json",
        "certified/P15_FULL_BN_T1_K2_RESIDUAL_SPECHT_QUOTIENT_CERTIFICATE.md",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: " + summary.next_task,
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = certify(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
