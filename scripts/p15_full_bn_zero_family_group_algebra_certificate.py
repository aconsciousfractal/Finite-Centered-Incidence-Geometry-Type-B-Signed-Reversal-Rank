#!/usr/bin/env python3
"""
Group-algebra certificate for the P15 n=8 zero family.

It proves the stronger finite statement behind

    M_{(1.1 | nu)}(B~I_8(3,3)) = 0 for every nu |- 6.

Instead of checking each Specht shape, it works before applying any Specht
representation.  In the little-group model, lambda_plus=(1,1) contributes the
sign of the induced permutation on the two plus-color coordinates.  After the
signed epsilon-sum, every surviving complement S_6 basis term appears exactly
with two plus-sector bijections, with equal unsigned coefficient and opposite
sign.  Hence each matrix block is zero in Q[S_6], and all Specht blocks vanish.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from p15_full_bn_zero_block_certificate import (
    Mask,
    Perm,
    masks_of_size,
    positions_of_mask,
    restricted_perm_data,
    sign_of_perm,
    sign_poly_target_coeff,
)


@dataclass(frozen=True)
class GroupAlgebraSummary:
    status: str
    claim: str
    n: int
    k: int
    l: int
    lambda_plus: str
    complement_size: int
    subset_count: int
    permutations_scanned: int
    active_terms: int
    group_algebra_keys: int
    nonzero_group_algebra_keys: int
    cancelling_pairs: int
    bad_pairs: int
    signed_coefficient_patterns: Dict[str, int]
    source_destination_type_patterns: Dict[str, int]
    all_terms_abs_coeff: int
    total_abs_before_cancellation: int
    total_abs_after_cancellation: int
    boundary: str
    runtime_seconds: float


def rho(n: int, i: int) -> int:
    return n - 1 - i


def subset_type(n: int, mask: Mask) -> str:
    pts = positions_of_mask(n, mask)
    if len(pts) != 2:
        return f"size_{len(pts)}"
    return "rho_pair" if rho(n, pts[0]) == pts[1] else "split_pairs"


def build_detailed_terms(n: int, target: Tuple[int, int], plus_size: int):
    subset_masks = masks_of_size(n, plus_size)
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full = (1 << n) - 1
    terms = defaultdict(list)
    permutations_scanned = 0
    for pi in itertools.permutations(range(n)):
        permutations_scanned += 1
        pi = tuple(pi)
        for src_si, A_mask in enumerate(subset_masks):
            B_mask = full ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            Ap_mask, sigmaA = restricted_perm_data(pi, A_mask, n)
            _Bp_mask, sigmaB = restricted_perm_data(pi, B_mask, n)
            dst_si = subset_index[Ap_mask]
            signed_coeff = coeff * sign_of_perm(sigmaA)
            terms[(src_si, dst_si, sigmaB)].append(
                {
                    "signed_coeff": signed_coeff,
                    "unsigned_coeff": coeff,
                    "sigmaA_sign": sign_of_perm(sigmaA),
                    "A_mask": A_mask,
                    "Ap_mask": Ap_mask,
                    "sigmaA": sigmaA,
                    "pi_on_A": tuple(pi[i] for i in positions_of_mask(n, A_mask)),
                    "pi": pi,
                }
            )
    return subset_masks, terms, permutations_scanned


def certify(out_json: Path | None, out_md: Path | None) -> GroupAlgebraSummary:
    t0 = time.time()
    n = 8
    target = (3, 3)
    plus_size = 2
    subset_masks, grouped, permutations_scanned = build_detailed_terms(n, target, plus_size)

    active_terms = sum(len(v) for v in grouped.values())
    nonzero_keys = 0
    cancelling_pairs = 0
    bad_pairs = 0
    signed_patterns: Counter = Counter()
    type_patterns: Counter = Counter()
    total_abs_before = 0
    total_abs_after = 0
    abs_coeffs = set()

    for (src_si, dst_si, sigmaB), vals in grouped.items():
        signed_sum = sum(v["signed_coeff"] for v in vals)
        total_abs_after += abs(signed_sum)
        if signed_sum != 0:
            nonzero_keys += 1
        signed_patterns[tuple(sorted(v["signed_coeff"] for v in vals))] += 1
        total_abs_before += sum(abs(v["signed_coeff"]) for v in vals)
        abs_coeffs.update(abs(v["unsigned_coeff"]) for v in vals)

        if len(vals) == 2 and vals[0]["signed_coeff"] + vals[1]["signed_coeff"] == 0:
            v0, v1 = vals
            A_mask = v0["A_mask"]
            A_pts = positions_of_mask(n, A_mask)
            comp_pts = [i for i in range(n) if not ((A_mask >> i) & 1)]
            pi0 = v0["pi"]
            pi1 = v1["pi"]
            same_complement = all(pi0[i] == pi1[i] for i in comp_pts)
            swapped_plus = pi0[A_pts[0]] == pi1[A_pts[1]] and pi0[A_pts[1]] == pi1[A_pts[0]]
            same_unsigned = v0["unsigned_coeff"] == v1["unsigned_coeff"]
            opposite_plus_sign = v0["sigmaA_sign"] == -v1["sigmaA_sign"]
            if same_complement and swapped_plus and same_unsigned and opposite_plus_sign:
                cancelling_pairs += 1
            else:
                bad_pairs += 1
        else:
            bad_pairs += 1

        src_type = subset_type(n, subset_masks[src_si])
        dst_type = subset_type(n, subset_masks[dst_si])
        unsigned_tuple = tuple(sorted(v["unsigned_coeff"] for v in vals))
        sign_tuple = tuple(sorted(v["sigmaA_sign"] for v in vals))
        type_patterns[(src_type, dst_type, unsigned_tuple, sign_tuple)] += 1

    status = "PASS" if nonzero_keys == 0 and bad_pairs == 0 and abs_coeffs == {4} else "FAIL"
    summary = GroupAlgebraSummary(
        status=status,
        claim="The lambda_plus=(1,1) block of B~I_8(3,3) is zero in every complement Q[S_6] block before applying Specht modules.",
        n=n,
        k=target[0],
        l=target[1],
        lambda_plus="1.1",
        complement_size=6,
        subset_count=len(subset_masks),
        permutations_scanned=permutations_scanned,
        active_terms=active_terms,
        group_algebra_keys=len(grouped),
        nonzero_group_algebra_keys=nonzero_keys,
        cancelling_pairs=cancelling_pairs,
        bad_pairs=bad_pairs,
        signed_coefficient_patterns={str(k): v for k, v in sorted(signed_patterns.items(), key=lambda kv: str(kv[0]))},
        source_destination_type_patterns={str(k): v for k, v in sorted(type_patterns.items(), key=lambda kv: str(kv[0]))},
        all_terms_abs_coeff=4 if abs_coeffs == {4} else -1,
        total_abs_before_cancellation=total_abs_before,
        total_abs_after_cancellation=total_abs_after,
        boundary="Finite group-algebra cancellation certificate only; no full native B_n fingerprint theorem or classification theorem.",
        runtime_seconds=round(time.time() - t0, 3),
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: GroupAlgebraSummary, path: Path) -> None:
    lines = [
        "# P15 Full Bn Zero Family Group-Algebra Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Interpretation",
        "",
        "For `lambda_plus=(1,1)`, the plus-sector Specht module is the sign representation of `S_2`.",
        "After summing signs over `B~I_8(3,3)`, every surviving complement `S_6` basis term occurs exactly twice: the two terms have the same complement map, the two plus images swapped, equal unsigned coefficient `4`, and opposite plus-sector sign.",
        "Therefore the alternating sum is already zero in the relevant block of `Q[S_6]`. Applying any Specht module `S^nu`, `nu |- 6`, preserves zero.",
        "",
        "## Audit",
        "",
        "| field | value |",
        "|---|---:|",
        f"| `subset_count` | {summary.subset_count} |",
        f"| `permutations_scanned` | {summary.permutations_scanned} |",
        f"| `active_terms` | {summary.active_terms} |",
        f"| `group_algebra_keys` | {summary.group_algebra_keys} |",
        f"| `nonzero_group_algebra_keys` | {summary.nonzero_group_algebra_keys} |",
        f"| `cancelling_pairs` | {summary.cancelling_pairs} |",
        f"| `bad_pairs` | {summary.bad_pairs} |",
        f"| `all_terms_abs_coeff` | {summary.all_terms_abs_coeff} |",
        f"| `total_abs_before_cancellation` | {summary.total_abs_before_cancellation} |",
        f"| `total_abs_after_cancellation` | {summary.total_abs_after_cancellation} |",
        "",
        "## Pattern",
        "",
        "```text",
        f"signed_coefficient_patterns = {summary.signed_coefficient_patterns}",
        f"source_destination_type_patterns = {summary.source_destination_type_patterns}",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: turn this finite cancellation into a short human-readable lemma, or return to final PDF/freeze routing.",
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