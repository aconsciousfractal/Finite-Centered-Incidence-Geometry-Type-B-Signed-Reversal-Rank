#!/usr/bin/env python3
"""
Sign-twist/conjugate-rank certificate for T_1^(8,2).

For Specht modules of S_7, S^{lambda'} is S^lambda tensor sign.  On the
compressed fixed-plus operator T_1^(8,2), the complement sign character factors
cellwise as a row sign times a column sign on the eight plus fibers. Therefore
sign-twisting the complement representation is rank-preserving for every
complement partition.

This is a finite n=8,k=2 certificate. It does not explain the self-conjugate or
middle residual ranks; it only proves the conjugate partner ranks follow from
already explained rows.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from p15_full_bn_prespecht_fixed_plus_scout import build_prespecht_terms, perm_sign

ROOT = Path(__file__).resolve().parents[1]
PRESPECHT_JSON = ROOT / "certified" / "P15_FULL_BN_PRESPECHT_FIXED_PLUS_SCOUT.json"


@dataclass(frozen=True)
class ConjugateRankPair:
    shape: str
    conjugate: str
    rank: int | None
    conjugate_rank: int | None
    rank_match: bool | None
    reading: str


@dataclass(frozen=True)
class SignTwistSummary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    group_keys: int
    active_cells: int
    parity_constant_cells: int
    bad_parity_cells: int
    row_signs: Dict[str, int]
    column_signs: Dict[str, int]
    factorization_bad_terms: int
    diagonal_cells: int
    reversal_cells: int
    conjugate_rank_pairs: List[ConjugateRankPair]
    runtime_seconds: float
    boundary: str
    next_task: str


def part_str(part: Tuple[int, ...]) -> str:
    return ".".join(str(x) for x in part) if part else "empty"


def parse_part(s: str) -> Tuple[int, ...]:
    if s == "empty":
        return ()
    return tuple(int(x) for x in s.split("."))


def conjugate(part: Tuple[int, ...]) -> Tuple[int, ...]:
    if not part:
        return ()
    return tuple(sum(1 for x in part if x >= j) for j in range(1, max(part) + 1))


def load_t1_ranks() -> Dict[str, int]:
    if not PRESPECHT_JSON.exists():
        return {}
    data = json.loads(PRESPECHT_JSON.read_text(encoding="utf-8"))
    for case in data.get("cases", []):
        if case.get("label") == "T_1^(8,2)":
            return {row["lambda_minus"]: int(row["rank_mod_p"]) for row in case.get("rows", [])}
    return {}


def make_signs(n: int) -> Tuple[Dict[int, int], Dict[int, int]]:
    col: Dict[int, int] = {}
    row: Dict[int, int] = {}
    for i in range(n):
        r = n - 1 - i
        col[i] = 1 if i < r else -1
    for i in range(n):
        row[i] = -col[i]
    return row, col


def certify(out_json: Path | None, out_md: Path | None) -> SignTwistSummary:
    t0 = time.time()
    n = 8
    k = 2
    _plus_masks, grouped, _perms, _raw_active, _raw_abs, _zero_removed = build_prespecht_terms(n, k, 1)
    cells: Dict[Tuple[int, int], set[int]] = defaultdict(set)
    for (dst, src, sigma), _coeff in grouped.items():
        cells[(dst, src)].add(perm_sign(sigma))

    row_sign, col_sign = make_signs(n)
    bad_terms = 0
    for (dst, src, sigma), _coeff in grouped.items():
        if perm_sign(sigma) != row_sign[dst] * col_sign[src]:
            bad_terms += 1

    bad_cells = sum(1 for signs in cells.values() if len(signs) != 1)
    diagonal = sum(1 for (dst, src) in cells if dst == src)
    reversal = sum(1 for (dst, src) in cells if dst == n - 1 - src)

    ranks = load_t1_ranks()
    selected = [
        "6.1",
        "5.2",
        "5.1.1",
        "4.3",
        "4.2.1",
        "4.1.1.1",
        "3.3.1",
    ]
    pairs: List[ConjugateRankPair] = []
    for shape in selected:
        conj = part_str(conjugate(parse_part(shape)))
        rank = ranks.get(shape)
        crank = ranks.get(conj)
        if shape == conj:
            reading = "self-conjugate; sign twist gives no new partner mechanism"
        elif shape in {"3.3.1"}:
            reading = "middle pair; sign twist proves rank equality, and the later residual certificate explains the rank"
        else:
            reading = "partner rank follows from the quotient-explained dominant-side shape"
        pairs.append(ConjugateRankPair(
            shape=shape,
            conjugate=conj,
            rank=rank,
            conjugate_rank=crank,
            rank_match=(rank == crank) if rank is not None and crank is not None else None,
            reading=reading,
        ))

    status = "PASS" if bad_terms == 0 and bad_cells == 0 and diagonal == 8 and reversal == 8 and all(p.rank_match is not False for p in pairs) else "FAIL"
    summary = SignTwistSummary(
        status=status,
        claim=(
            "For T_1^(8,2), the complement sign character factors as row_sign(dst)*column_sign(src) "
            "on every nonzero compressed group-algebra term. Hence sign-twisting the complement "
            "Specht representation, equivalently conjugating the complement partition, preserves rank."
        ),
        n=n,
        k=k,
        lambda_plus="1",
        group_keys=len(grouped),
        active_cells=len(cells),
        parity_constant_cells=len(cells) - bad_cells,
        bad_parity_cells=bad_cells,
        row_signs={str(i): row_sign[i] for i in range(n)},
        column_signs={str(i): col_sign[i] for i in range(n)},
        factorization_bad_terms=bad_terms,
        diagonal_cells=diagonal,
        reversal_cells=reversal,
        conjugate_rank_pairs=pairs,
        runtime_seconds=round(time.time() - t0, 3),
        boundary=(
            "Finite n=8,k=2 fixed-plus sign-twist certificate only. It transfers ranks to conjugate "
            "complement shapes, but it is not a standalone full extra-defect classification. The residual "
            "middle/self-conjugate mechanisms are handled by the later residual Specht quotient certificate."
        ),
        next_task=(
            "Attack the next compact one-row fixed-plus family T_4^(8,3) (lambda_plus=4, n=8,k=3), "
            "and keep this package appendix-scoped; rational T_1/T_2/T_4 audits are separate optional checks."
        ),
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: SignTwistSummary, path: Path) -> None:
    lines = [
        "# P15 T1 k=2 Sign-Twist Conjugate Symmetry Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Cell Factorization",
        "",
        "The compressed operator has exactly 16 active source/destination cells: 8 diagonal cells and 8 reversal cells. In every active cell all complement permutations have the same parity.",
        "",
        "The parity factors as:",
        "",
        "```text",
        "sgn(sigma_B) = row_sign[dst] * column_sign[src]",
        "```",
        "",
        "Audit:",
        "",
        "| field | value |",
        "|---|---:|",
        f"| group keys | {summary.group_keys} |",
        f"| active cells | {summary.active_cells} |",
        f"| parity-constant cells | {summary.parity_constant_cells} |",
        f"| bad parity cells | {summary.bad_parity_cells} |",
        f"| factorization bad terms | {summary.factorization_bad_terms} |",
        f"| diagonal cells | {summary.diagonal_cells} |",
        f"| reversal cells | {summary.reversal_cells} |",
        "",
        "## Consequence",
        "",
        "Since `S^(lambda') = S^lambda tensor sgn`, the block for `lambda'` is obtained from the block for `lambda` by invertible row and column sign matrices. Therefore conjugate complement shapes have equal rank in this fixed-plus family.",
        "",
        "| shape | conjugate | rank | conjugate rank | reading |",
        "|---|---|---:|---:|---|",
    ]
    for p in summary.conjugate_rank_pairs:
        lines.append(f"| `{p.shape}` | `{p.conjugate}` | {p.rank} | {p.conjugate_rank} | {p.reading} |")
    lines += [
        "",
        "## Sources",
        "",
        "```text",
        "scripts/p15_full_bn_t1_k2_sign_twist_conjugate_certificate.py",
        "artifacts/certified/P15_FULL_BN_T1_K2_SIGN_TWIST_CONJUGATE_CERTIFICATE.json",
        "artifacts/certified/P15_FULL_BN_T1_K2_SIGN_TWIST_CONJUGATE_CERTIFICATE.md",
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
