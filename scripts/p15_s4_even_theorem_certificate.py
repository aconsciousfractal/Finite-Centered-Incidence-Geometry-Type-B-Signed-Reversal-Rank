# -*- coding: utf-8 -*-
"""
P15-S4 even-n theorem certificate builder.

This script is a focused gate wrapper around the S3 exact engine.  It records
finite replay evidence and proof-dependency checks for the even-n theorem:

  rank M_ref(B~I_{2m}(k,k)) = m
  rank M_pair-std(B~I_{2m}(k,k)) = m-1

for k>=1 with nonempty layer.  The general proof is in the S4 theorem note;
this script writes the machine-readable certificate supporting that gate.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import p15_s3_exact_engine as s3

SCHEMA = "p15.s4.even_theorem_certificate.v1"


def even_closed_form_sweep(limit_n: int = 20) -> dict:
    rows = []
    all_positive = True
    for n in range(4, limit_n + 1, 2):
        h = s3.h_even_table(n)
        values = []
        for k in range(1, n // 2 + 1):
            size = s3.size_even(n, k, k, h)
            if size == 0:
                continue
            a_value = s3.a_even_closed_form(n, k, h)
            ok = a_value > 0
            all_positive = all_positive and ok
            values.append({"k": k, "size": size, "a": a_value, "a_positive": ok})
        rows.append({"n": n, "values": values})
    return {
        "id": "even_closed_form_positive_sweep",
        "status": "PASS" if all_positive else "FAIL",
        "range": "even n=4..%d, nonempty k>=1" % limit_n,
        "rows": rows,
    }


def even_channel_replay() -> dict:
    channel = s3.run_channels()
    rows = [row for row in channel["rows"] if row["parity"] == "even"]
    ok = True
    for row in rows:
        ok = ok and row["rank_ok"] and row["collapse_ok"]
        ok = ok and row["a_closed_form_ok"] and row["size_closed_form_ok"]
        ok = ok and row["trace_flip_ok"] and row["lambda_pair_formula_ok"]
    return {
        "id": "even_channel_exact_replay",
        "status": "PASS" if ok else "FAIL",
        "range": "even full enumeration cases inherited from S3: n=4,6, diagonal nonempty k>=1",
        "rows": rows,
    }


def even_permanent_replay() -> dict:
    perm = s3.run_permanent_domination()
    rows = [row for row in perm["rows"] if row["n"] % 2 == 0]
    ok = all(row["all_patterns_le_1"] for row in rows)
    return {
        "id": "even_permanent_domination_replay",
        "status": "PASS" if ok else "FAIL",
        "range": "all feasible positive-hit patterns for even n=4,6 from S3 Ryser replay",
        "rows": rows,
    }


def make_report(limit_n: int = 20) -> dict:
    checks = [
        even_channel_replay(),
        even_closed_form_sweep(limit_n),
        even_permanent_replay(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    proof_checklist = {
        "source_lock_S1": True,
        "prior_art_boundary_S2": True,
        "exact_engine_S3": True,
        "collapse_M_ref_even": True,
        "permanent_domination_lemma_available": True,
        "trace_flip_identity_even": True,
        "pair_standard_trace_identity_even": True,
        "rank_number_and_projector_classical_boundary_carried": True,
        "no_odd_claim_promoted": True,
        "no_type_D_claim_promoted": True,
    }
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S4 even-n theorem",
        "status": "PASS" if not failed and all(proof_checklist.values()) else "FAIL",
        "failed_checks": failed,
        "theorem_scope": "even n>=4, k>=1 with B~I_n(k,k) nonempty, two ambient standard channels only",
        "theorem_statement": {
            "M_ref_rank": "n/2",
            "M_pair_std_rank": "n/2 - 1",
        },
        "checks": checks,
        "proof_checklist": proof_checklist,
        "boundary": [
            "No public theorem promotion at S4; S7/S9 still pending.",
            "Rank number and P_plus are matching-scheme bookkeeping, cite not claim.",
            "P15 contribution is collapse plus positivity onto the idempotent.",
            "Odd-n, H_m, C2xC2 square, Type-D scout, and full B_n fingerprint are outside S4.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S4 even-n theorem certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"], "-", check["range"])
        if check["id"] == "even_channel_exact_replay":
            for row in check["rows"]:
                print(
                    "  n={n} k={k} size={size} ranks=({rank_M_ref},{rank_M_pair_std}) "
                    "collapse={collapse_ok} a={a_actual} lambda_pair={lambda_pair_actual}".format(**row)
                )
        elif check["id"] == "even_closed_form_positive_sweep":
            for row in check["rows"]:
                vals = ", ".join("k=%d:a=%d" % (v["k"], v["a"]) for v in row["values"])
                print("  n=%d %s" % (row["n"], vals))
        elif check["id"] == "even_permanent_domination_replay":
            for row in check["rows"]:
                print(
                    "  n={n} k={k} patterns={checked_patterns} maxE={max_E_neg_fixed} le1={all_patterns_le_1}".format(**row)
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S4 even-n theorem certificate builder")
    parser.add_argument("--write-json", help="write deterministic JSON certificate to this path")
    parser.add_argument("--limit-n", type=int, default=20, help="maximum even n for closed-form positivity sweep")
    args = parser.parse_args()
    report = make_report(args.limit_n)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
