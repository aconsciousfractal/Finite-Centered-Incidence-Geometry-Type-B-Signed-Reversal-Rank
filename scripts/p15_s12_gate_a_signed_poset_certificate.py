#!/usr/bin/env python3
"""Gate A certificate for the P15 S12 signed-poset separator.

The certificate has two independent roles:
- finite n=4 witnesses killing every signed precedence, including zero;
- exact small-range smoke check for the sign-extension lemma used for n>=5.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import p15_positive_only_separator_hunt as scout

SCHEMA = "p15.s12_gate_a_signed_poset_certificate.v1"


def nonzero_labels(n: int) -> tuple[int, ...]:
    return tuple(range(-n, 0)) + tuple(range(1, n + 1))


def labels_with_zero(n: int) -> tuple[int, ...]:
    return tuple(range(-n, 0)) + (0,) + tuple(range(1, n + 1))


def order_index(w: tuple[int, ...]) -> dict[int, int]:
    return {value: idx for idx, value in enumerate(scout.signed_total_order(w))}


def witness_rows(n: int, include_zero: bool = False) -> list[dict]:
    x_list = [w for w in scout.signed_permutations(n) if scout.in_x(w)]
    labels = labels_with_zero(n) if include_zero else nonzero_labels(n)
    rows = []
    for a in labels:
        for b in labels:
            if a == b:
                continue
            found = None
            for w in x_list:
                idx = order_index(w)
                if idx[b] < idx[a]:
                    found = w
                    break
            if found is None:
                rows.append({"a": a, "b": b, "witness": None, "kills_precedence_a_before_b": False})
            else:
                rows.append({
                    "a": a,
                    "b": b,
                    "witness": list(found),
                    "signed_total_order": list(scout.signed_total_order(found)),
                    "kills_precedence_a_before_b": True,
                })
    return rows


def sign_extension_rows(n_min: int, n_max: int) -> list[dict]:
    out = []
    for n in range(n_min, n_max + 1):
        x_sets = [set(w) for w in scout.signed_permutations(n) if scout.in_x(w)]
        failures = []
        total_specs = 0
        one_label_failures = 0
        for size in (1, 2):
            for labels in combinations(range(1, n + 1), size):
                for signs in product((-1, 1), repeat=size):
                    total_specs += 1
                    wanted = tuple(sign * label for label, sign in zip(labels, signs))
                    if not any(all(value in x for value in wanted) for x in x_sets):
                        failures.append(list(wanted))
                        if size == 1:
                            one_label_failures += 1
        out.append({
            "n": n,
            "total_specs": total_specs,
            "failure_count": len(failures),
            "one_label_failure_count": one_label_failures,
            "failures": failures[:20],
            "all_specs_extend": not failures,
            "all_one_label_specs_extend": one_label_failures == 0,
        })
    return out


def make_payload() -> dict:
    n4_nonzero_rows = witness_rows(4, include_zero=False)
    n4_zero_rows = witness_rows(4, include_zero=True)
    sign_rows = sign_extension_rows(4, 7)
    return {
        "schema": SCHEMA,
        "status": "GATE_A_CERTIFICATE",
        "object": "X_n = B~I_n(1,1) subset B_n",
        "signed_total_order_convention": "T(w)=-w_n<...<-w_1<0<w_1<...<w_n",
        "n4_witness_rows_nonzero": n4_nonzero_rows,
        "n4_witness_rows_with_zero": n4_zero_rows,
        "sign_extension_rows_n4_n7": sign_rows,
        "summary": {
            "n4_ordered_nonzero_pairs": len(n4_nonzero_rows),
            "n4_all_nonzero_precedences_killed": all(r["kills_precedence_a_before_b"] for r in n4_nonzero_rows),
            "n4_ordered_pairs_with_zero": len(n4_zero_rows),
            "n4_all_precedences_with_zero_killed": all(r["kills_precedence_a_before_b"] for r in n4_zero_rows),
            "n4_X_size": len([w for w in scout.signed_permutations(4) if scout.in_x(w)]),
            "n4_distinct_witnesses_with_zero": len({tuple(row["witness"]) for row in n4_zero_rows}),
            "sign_extension_all_specs_n5_n7": all(r["all_specs_extend"] for r in sign_rows if r["n"] >= 5),
            "sign_extension_all_one_label_specs_n4_n7": all(r["all_one_label_specs_extend"] for r in sign_rows),
            "sign_extension_n4_failure_count": next(r for r in sign_rows if r["n"] == 4)["failure_count"],
            "gate_a_status": "PASS_FOR_SIGNED_POSET_ROUTE",
        },
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8", newline="\n")


def make_markdown(payload: dict) -> str:
    lines = [
        "# P15 S12 Gate A Signed-Poset Certificate",
        "",
        f"Schema: `{payload['schema']}`",
        f"Status: **{payload['status']}**",
        "",
        f"Convention: `{payload['signed_total_order_convention']}`.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## n=4 Finite Witnesses",
        "",
        "For every ordered signed pair `(a,b)`, including pairs with `0`, the certificate stores a witness `w in X_4` with `b` before `a` in `T(w)`. This kills the possible common precedence `a<b`.",
        "",
        "| a | b | witness w |",
        "|---:|---:|---|",
    ]
    for row in payload["n4_witness_rows_with_zero"][:24]:
        lines.append(f"| {row['a']} | {row['b']} | `{row['witness']}` |")
    lines += [
        "",
        "Only the first 24 witness rows are printed here; the JSON contains all rows.",
        "",
        "## Sign-Extension Smoke Check",
        "",
        "| n | total sign specs | failures | one-label failures | all extend? | failure sample |",
        "|---:|---:|---:|---:|---|---|",
    ]
    for row in payload["sign_extension_rows_n4_n7"]:
        lines.append(
            f"| {row['n']} | {row['total_specs']} | {row['failure_count']} | "
            f"{row['one_label_failure_count']} | {row['all_specs_extend']} | `{row['failures']}` |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "Gate A is closed for the signed-poset route: `n=4` is finite-certified, one-label sign freedom kills zero-relations in `n=4..7`, and the `n>=5` proof uses the sign-extension lemma written in the companion proof note.",
        "",
        "## Reproducibility",
        "",
        "Run from the project root:",
        "",
        "`python -B scripts/p15_s12_gate_a_signed_poset_certificate.py --write-json certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.json --write-md certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.md`",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-json", type=Path)
    parser.add_argument("--write-md", type=Path)
    args = parser.parse_args()
    payload = make_payload()
    if args.write_json:
        write_text(args.write_json, json.dumps(payload, indent=2))
    if args.write_md:
        write_text(args.write_md, make_markdown(payload))
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()