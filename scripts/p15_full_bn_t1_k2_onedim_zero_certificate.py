#!/usr/bin/env python3
"""
Exact one-dimensional quotient certificate for T_1^(8,2).

The pre-Specht scout shows that the fixed-plus family lambda_plus=1 at
B~I_8(2,2) has two zero rows: complement shapes (7) and (1^7).  This script
certifies the stronger cellwise reason.

After compressing to T_1^(8,2) in Mat_8(Q[S_7]), only cells j=i or j=rho(i)
survive after applying either one-dimensional complement character.  For every
surviving cell before character evaluation there are exactly four complement
permutation terms, with coefficients 15,-5,-5,-5 and a common complement parity.
Therefore both the trivial and sign character sums vanish cellwise.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from p15_full_bn_prespecht_fixed_plus_scout import build_prespecht_terms, mod_rank, perm_sign


@dataclass(frozen=True)
class T1OneDimSummary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    subset_count: int
    group_keys: int
    active_cells: int
    diagonal_cells: int
    reversal_cells: int
    bad_cells: int
    coefficient_patterns: Dict[str, int]
    parity_patterns: Dict[str, int]
    trivial_nonzero_entries: int
    sign_nonzero_entries: int
    trivial_rank_mod_p: int
    sign_rank_mod_p: int
    runtime_seconds: float
    boundary: str


def rho(n: int, i: int) -> int:
    return n - 1 - i


def matrix_for_character(grouped, subset_count: int, char) -> List[List[int]]:
    M = [[0 for _ in range(subset_count)] for __ in range(subset_count)]
    for (dst, src, sigma), coeff in grouped.items():
        M[dst][src] += coeff * char(sigma)
    return M


def rank_mod_list(M: List[List[int]], p: int) -> int:
    import numpy as np
    return mod_rank(np.array(M, dtype=np.int64), p)


def certify(out_json: Path | None, out_md: Path | None) -> T1OneDimSummary:
    t0 = time.time()
    n = 8
    k = 2
    subset_masks, grouped, _perms, _raw_active, _raw_abs, _zero_removed = build_prespecht_terms(n, k, 1)
    cells = defaultdict(list)
    for (dst, src, sigma), coeff in grouped.items():
        cells[(dst, src)].append((sigma, coeff, perm_sign(sigma)))

    coeff_patterns: Counter = Counter()
    parity_patterns: Counter = Counter()
    bad_cells = 0
    diagonal_cells = 0
    reversal_cells = 0
    for (dst, src), vals in cells.items():
        coeff_tuple = tuple(sorted(coeff for _sigma, coeff, _parity in vals))
        parity_tuple = tuple(sorted(parity for _sigma, _coeff, parity in vals))
        coeff_patterns[str(coeff_tuple)] += 1
        parity_patterns[str(parity_tuple)] += 1
        if dst == src:
            diagonal_cells += 1
        elif dst == rho(n, src):
            reversal_cells += 1
        else:
            bad_cells += 1
        if len(vals) != 4 or coeff_tuple != (-5, -5, -5, 15) or len(set(p for _s, _c, p in vals)) != 1:
            bad_cells += 1

    M_triv = matrix_for_character(grouped, len(subset_masks), lambda sigma: 1)
    M_sign = matrix_for_character(grouped, len(subset_masks), perm_sign)
    triv_nonzero = sum(1 for row in M_triv for x in row if x)
    sign_nonzero = sum(1 for row in M_sign for x in row if x)
    status = "PASS" if bad_cells == 0 and triv_nonzero == 0 and sign_nonzero == 0 and len(cells) == 16 else "FAIL"
    summary = T1OneDimSummary(
        status=status,
        claim="The two one-dimensional complement quotients of T_1^(8,2), namely lambda_minus=(7) and (1^7), vanish cellwise: each active identity/reversal cell has coefficients 15,-5,-5,-5 with common complement parity, so both trivial and sign character sums are zero.",
        n=n,
        k=k,
        lambda_plus="1",
        complement_size=7,
        subset_count=len(subset_masks),
        group_keys=len(grouped),
        active_cells=len(cells),
        diagonal_cells=diagonal_cells,
        reversal_cells=reversal_cells,
        bad_cells=bad_cells,
        coefficient_patterns={k: v for k, v in sorted(coeff_patterns.items())},
        parity_patterns={k: v for k, v in sorted(parity_patterns.items())},
        trivial_nonzero_entries=triv_nonzero,
        sign_nonzero_entries=sign_nonzero,
        trivial_rank_mod_p=rank_mod_list(M_triv, 1000003),
        sign_rank_mod_p=rank_mod_list(M_sign, 1000003),
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite n=8,k=2 fixed-plus one-dimensional quotient certificate only; no full native B_n fingerprint theorem or classification theorem.",
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: T1OneDimSummary, path: Path) -> None:
    lines = [
        "# P15 T1 k=2 One-Dimensional Zero Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Mechanism",
        "",
        "Compress `B~I_8(2,2)` at fixed plus shape `lambda_plus=1` to",
        "",
        "```text",
        "T_1^(8,2) in Mat_8(Q[S_7]).",
        "```",
        "",
        "The two zero rows in the n=8 defect CSV are the one-dimensional complement characters `(7)` and `(1^7)`. The certificate shows that their vanishing is cellwise, not a numerical rank accident.",
        "",
        "For every active source/destination cell, the destination is either the same point or its reversal. The four complement terms in that cell have coefficient pattern `(-5,-5,-5,15)`. Their complement permutation parity is constant inside the cell. Therefore both character sums vanish:",
        "",
        "```text",
        "15 - 5 - 5 - 5 = 0,",
        "sgn * (15 - 5 - 5 - 5) = 0.",
        "```",
        "",
        "## Audit",
        "",
        "| field | value |",
        "|---|---:|",
        f"| `group_keys` | {summary.group_keys} |",
        f"| `active_cells` | {summary.active_cells} |",
        f"| `diagonal_cells` | {summary.diagonal_cells} |",
        f"| `reversal_cells` | {summary.reversal_cells} |",
        f"| `bad_cells` | {summary.bad_cells} |",
        f"| `trivial_nonzero_entries` | {summary.trivial_nonzero_entries} |",
        f"| `sign_nonzero_entries` | {summary.sign_nonzero_entries} |",
        f"| `trivial_rank_mod_p` | {summary.trivial_rank_mod_p} |",
        f"| `sign_rank_mod_p` | {summary.sign_rank_mod_p} |",
        "",
        "Coefficient patterns:",
        "",
        "```text",
        f"{summary.coefficient_patterns}",
        "```",
        "",
        "Parity patterns:",
        "",
        "```text",
        f"{summary.parity_patterns}",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: keep this package appendix-scoped; optional rational audits and any remaining families are separate finite checks, not a classification theorem.",
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
