#!/usr/bin/env python3
"""
Exact finite certificate for the remaining n=8,k=3 zero block

    M_{(2.2.2.1.1 | empty)}(B~I_8(3,3)) = 0.

This is the twelfth ZERO row in the n=8,k=3 full native B_n atlas.  It is
separate from the lambda_plus=(1,1) zero family.  Since lambda_minus is empty,
the block is an S_8 Specht block with scalar sign-count weights.
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

from p15_full_bn_zero_block_certificate import (
    Perm,
    SpechtSeminormalExact,
    gcd_lcm,
    layer_size,
    sign_poly_target_coeff,
)


@dataclass(frozen=True)
class PurePlusSummary:
    status: str
    claim: str
    n: int
    k: int
    l: int
    lambda_plus: str
    lambda_minus: str
    specht_dim: int
    nonzero_underlying_permutations: int
    layer_size: int
    weight_patterns: Dict[str, int]
    common_denominator_seen: int
    max_abs_numerator_seen: int
    nonzero_entries: int
    relation_checks_pass: bool
    relation_checks: Dict[str, bool]
    runtime_seconds: float
    boundary: str


def rho(n: int, i: int) -> int:
    return n - 1 - i


def nonzero_weight_rows(n: int, target: Tuple[int, int]) -> List[Tuple[Perm, int, int, int]]:
    rows: List[Tuple[Perm, int, int, int]] = []
    k, l = target
    for pi0 in itertools.permutations(range(n)):
        pi = tuple(pi0)
        I = sum(1 for i, v in enumerate(pi) if v == i)
        R = sum(1 for i, v in enumerate(pi) if v == rho(n, i))
        if I >= k and R >= l:
            # lambda_minus is empty, so the sign sum is just the number of sign choices.
            w = sign_poly_target_coeff(pi, 0, target)
            if w:
                rows.append((pi, I, R, w))
    return rows


def certify(out_json: Path | None, out_md: Path | None) -> PurePlusSummary:
    t0 = time.time()
    n = 8
    target = (3, 3)
    shape = (2, 2, 2, 1, 1)
    rows = nonzero_weight_rows(n, target)
    S = SpechtSeminormalExact(shape)
    checks = S.relation_checks()
    relation_checks_pass = all(checks.values())
    if not relation_checks_pass:
        raise AssertionError(checks)

    agg = [[0 for _ in range(S.dim)] for __ in range(S.dim)]
    denom_seen = 1
    max_num_seen = 0
    weight_patterns: Dict[str, int] = {}

    for pi, I, R, w in rows:
        weight_patterns[str((I, R, w))] = weight_patterns.get(str((I, R, w)), 0) + 1
        M = S.mat_perm(pi)
        for a in range(S.dim):
            for b in range(S.dim):
                x = M[a][b]
                if x:
                    denom_seen = gcd_lcm(denom_seen, x.denominator)
                    max_num_seen = max(max_num_seen, abs(x.numerator))
                    agg[a][b] += w * x

    nonzero_entries = sum(1 for row in agg for x in row if x != 0)
    status = "PASS" if nonzero_entries == 0 and sum(w for _pi, _I, _R, w in rows) == layer_size(n, target) else "FAIL"
    summary = PurePlusSummary(
        status=status,
        claim="M_{(2.2.2.1.1 | empty)}(B~I_8(3,3)) = 0 in the rational S_8 Specht model",
        n=n,
        k=target[0],
        l=target[1],
        lambda_plus="2.2.2.1.1",
        lambda_minus="empty",
        specht_dim=S.dim,
        nonzero_underlying_permutations=len(rows),
        layer_size=layer_size(n, target),
        weight_patterns=dict(sorted(weight_patterns.items())),
        common_denominator_seen=denom_seen,
        max_abs_numerator_seen=max_num_seen,
        nonzero_entries=nonzero_entries,
        relation_checks_pass=relation_checks_pass,
        relation_checks=checks,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite exact pure-plus zero-block witness only; no full native B_n fingerprint theorem or classification theorem.",
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: PurePlusSummary, path: Path) -> None:
    lines = [
        "# P15 Full Bn Pure-Plus Zero Block Certificate",
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
        "This is the remaining `ZERO` row in the n=8, k=3 full native `B_n` atlas after the `(1.1|nu)`, `nu |- 6`, zero family is removed.",
        "Because `lambda_minus` is empty, signs carry no nontrivial character. The coefficient of an underlying permutation `pi in S_8` is exactly the number of sign choices producing three positive identity hits and three positive reversal hits.",
        "Only 102 underlying permutations have nonzero coefficient. Summing their exact Young-seminormal `S_8` matrices on shape `(2,2,2,1,1)`, with these integer weights, gives the zero rational matrix.",
        "",
        "## Exact Audit",
        "",
        "| field | value |",
        "|---|---:|",
        f"| `specht_dim` | {summary.specht_dim} |",
        f"| `nonzero_underlying_permutations` | {summary.nonzero_underlying_permutations} |",
        f"| `layer_size` | {summary.layer_size} |",
        f"| `common_denominator_seen` | {summary.common_denominator_seen} |",
        f"| `max_abs_numerator_seen` | {summary.max_abs_numerator_seen} |",
        f"| `nonzero_entries` | {summary.nonzero_entries} |",
        "",
        "## Weight Patterns",
        "",
        "```text",
        f"{summary.weight_patterns}",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: record the two zero mechanisms together in the n=8 defect note, then decide whether to search for a conceptual proof of the pure-plus zero or keep it as a finite certified witness.",
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
