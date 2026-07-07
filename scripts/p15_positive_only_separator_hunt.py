#!/usr/bin/env python3
"""Positive-only separator hunt for P15.

Exploratory exact checks for X_n = B~I_n(1,1):
- Reiner-style signed-poset common-precedence criterion on signed total orders;
- necessary translate-tiler divisibility |X_n| | |B_n|;
- permanent/DP count for |X_n|.

This is a research scout, not a promoted P15 theorem certificate.
"""
from __future__ import annotations

import argparse
import json
from itertools import permutations, product
from math import factorial, gcd
from pathlib import Path
from typing import Iterable

SCHEMA = "p15.positive_only_separator_hunt.v1"


def signed_permutations(n: int) -> Iterable[tuple[int, ...]]:
    for perm in permutations(range(1, n + 1)):
        for signs in product((-1, 1), repeat=n):
            yield tuple(signs[i] * perm[i] for i in range(n))


def pos_fixed_count(w: tuple[int, ...]) -> int:
    return sum(1 for i, v in enumerate(w, start=1) if v == i)


def pos_reversal_count(w: tuple[int, ...]) -> int:
    n = len(w)
    return sum(1 for i, v in enumerate(w, start=1) if v == n + 1 - i)


def in_x(w: tuple[int, ...]) -> bool:
    return pos_fixed_count(w) == 1 and pos_reversal_count(w) == 1


def signed_total_order(w: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-v for v in reversed(w)) + (0,) + w


def precedence_pairs(order: tuple[int, ...]) -> set[tuple[int, int]]:
    return {(a, b) for i, a in enumerate(order) for b in order[i + 1 :]}


def transitive_closure(pairs: set[tuple[int, int]], elements: tuple[int, ...]) -> set[tuple[int, int]]:
    reach = {a: set() for a in elements}
    for a, b in pairs:
        reach[a].add(b)
    changed = True
    while changed:
        changed = False
        for a in elements:
            add: set[int] = set()
            for b in tuple(reach[a]):
                add.update(reach[b])
            before = len(reach[a])
            reach[a].update(add)
            changed = changed or len(reach[a]) != before
    return {(a, b) for a in elements for b in reach[a] if a != b}


def satisfies_relation(w: tuple[int, ...], relation: set[tuple[int, int]]) -> bool:
    pos = {a: i for i, a in enumerate(signed_total_order(w))}
    return all(pos[a] < pos[b] for a, b in relation)


def common_precedence_analysis(n: int) -> dict:
    x_list: list[tuple[int, ...]] = []
    common: set[tuple[int, int]] | None = None
    for w in signed_permutations(n):
        if not in_x(w):
            continue
        x_list.append(w)
        pairs = precedence_pairs(signed_total_order(w))
        common = pairs if common is None else common & pairs
    if common is None:
        common = set()
    elements = tuple(range(-n, 0)) + (0,) + tuple(range(1, n + 1))
    closure = transitive_closure(common, elements)
    group_size = (2**n) * factorial(n)
    if not closure:
        lb_count = group_size
        exact_poset = len(x_list) == group_size
    else:
        x_set = set(x_list)
        lb_count = 0
        exact_poset = True
        for w in signed_permutations(n):
            ok = satisfies_relation(w, closure)
            if ok:
                lb_count += 1
            if ok != (w in x_set):
                exact_poset = False
    nonzero_common = {(a, b) for a, b in common if a != 0 and b != 0}
    with_zero_common = common - nonzero_common
    return {
        "n": n,
        "X_size_enumerated": len(x_list),
        "group_size": group_size,
        "common_precedence_size": len(common),
        "common_precedence_nonzero_size": len(nonzero_common),
        "common_precedence_with_zero_size": len(with_zero_common),
        "closure_size": len(closure),
        "LB_PX_size": lb_count,
        "exact_signed_poset": exact_poset,
        "common_precedence_sample": sorted(list(common), key=lambda p: (p[0], p[1]))[:12],
    }


def count_x_dp(n: int) -> int:
    full = (1 << n) - 1
    dp: dict[int, list[list[int]]] = {0: [[1, 0], [0, 0]]}
    for mask in range(1 << n):
        if mask not in dp:
            continue
        i = mask.bit_count()
        if i >= n:
            continue
        table = dp[mask]
        for j in range(n):
            bit = 1 << j
            if mask & bit:
                continue
            rev = n - 1 - i
            if j == i and j == rev:
                terms = [(0, 0, 1), (1, 1, 1)]
            elif j == i:
                terms = [(0, 0, 1), (1, 0, 1)]
            elif j == rev:
                terms = [(0, 0, 1), (0, 1, 1)]
            else:
                terms = [(0, 0, 2)]
            newmask = mask | bit
            out = dp.setdefault(newmask, [[0, 0], [0, 0]])
            for a in range(2):
                for b in range(2):
                    val = table[a][b]
                    if not val:
                        continue
                    for da, db, weight in terms:
                        na, nb = a + da, b + db
                        if na <= 1 and nb <= 1:
                            out[na][nb] += val * weight
    return dp[full][1][1]


def divisibility_rows(n_min: int, n_max: int) -> list[dict]:
    rows = []
    for n in range(n_min, n_max + 1):
        x_size = count_x_dp(n)
        group_size = (2**n) * factorial(n)
        rows.append({
            "n": n,
            "X_size": x_size,
            "B_n_size": group_size,
            "B_n_mod_X": group_size % x_size if x_size else None,
            "gcd_B_n_X": gcd(group_size, x_size),
            "divides_B_n": bool(x_size and group_size % x_size == 0),
            "tile_count_if_divides": group_size // x_size if x_size and group_size % x_size == 0 else None,
        })
    return rows

def compose_signed(g: tuple[int, ...], a: tuple[int, ...]) -> tuple[int, ...]:
    out = []
    for x in a:
        y = g[abs(x) - 1]
        out.append(y if x > 0 else -y)
    return tuple(out)


def exact_left_translate_tiling(n: int) -> dict:
    group = list(signed_permutations(n))
    group_set = set(group)
    x_set = frozenset(w for w in group if in_x(w))
    if not x_set or len(group) % len(x_set):
        return {
            "n": n,
            "group_size": len(group),
            "X_size": len(x_set),
            "target_tile_count": None,
            "distinct_left_translates": None,
            "left_translate_tiling_exists": False,
            "chosen_left_translate_reps": [],
            "reason": "empty_or_cardinality_obstructed",
        }
    translates = []
    translate_reps = []
    seen = set()
    for g in group:
        tile = frozenset(compose_signed(g, x) for x in x_set)
        if tile not in seen:
            seen.add(tile)
            translates.append(tile)
            translate_reps.append(g)
    target = len(group) // len(x_set)
    by_element = {w: [] for w in group}
    for idx, tile in enumerate(translates):
        for w in tile:
            by_element[w].append(idx)
    chosen: list[int] = []

    def search(covered: set[tuple[int, ...]]) -> bool:
        if len(chosen) == target:
            return len(covered) == len(group_set)
        if len(covered) + (target - len(chosen)) * len(x_set) < len(group_set):
            return False
        best = None
        for w in group_set - covered:
            candidates = [idx for idx in by_element[w] if idx not in chosen and translates[idx].isdisjoint(covered)]
            if best is None or len(candidates) < len(best):
                best = candidates
                if not best:
                    return False
        if best is None:
            return True
        for idx in best:
            chosen.append(idx)
            if search(covered | set(translates[idx])):
                return True
            chosen.pop()
        return False

    exists = search(set())
    return {
        "n": n,
        "group_size": len(group),
        "X_size": len(x_set),
        "target_tile_count": target,
        "distinct_left_translates": len(translates),
        "left_translate_tiling_exists": exists,
        "chosen_left_translate_reps": [list(translate_reps[idx]) for idx in chosen] if exists else [],
        "reason": "exact_cover_search" if exists else "exact_cover_search_no_cover",
    }


def make_payload() -> dict:
    poset_rows = [common_precedence_analysis(n) for n in range(3, 8)]
    div_rows = divisibility_rows(3, 12)
    tiler_exact_rows = [exact_left_translate_tiling(3)]
    enum_match = []
    for r in poset_rows:
        n = r["n"]
        dp_count = count_x_dp(n)
        enum_match.append({
            "n": n,
            "enumerated": r["X_size_enumerated"],
            "dp": dp_count,
            "match": r["X_size_enumerated"] == dp_count,
        })
    return {
        "schema": SCHEMA,
        "object": "X_n = B~I_n(1,1) subset B_n",
        "status": "RESEARCH_SCOUT_NOT_PROMOTED",
        "signed_total_order_convention": "w -> -w_n < ... < -w_1 < 0 < w_1 < ... < w_n",
        "count_formula": "|X_n| = [xy] per(W_n(x,y)), W_ij=1+x if j=i only, 1+y if j=rho(i) only, 1+xy at the odd center, and 2 otherwise.",
        "poset_rows_n3_n7": poset_rows,
        "divisibility_rows_n3_n12": div_rows,
        "tiler_exact_rows": tiler_exact_rows,
        "enumeration_vs_dp_checks": enum_match,
        "summary": {
            "all_enumeration_dp_checks_match": all(r["match"] for r in enum_match),
            "common_precedence_empty_n4_n7": all(r["common_precedence_size"] == 0 for r in poset_rows if r["n"] >= 4),
            "nonzero_common_precedence_empty_n4_n7": all(r["common_precedence_nonzero_size"] == 0 for r in poset_rows if r["n"] >= 4),
            "all_poset_exact_false_n3_n7": all(not r["exact_signed_poset"] for r in poset_rows),
            "divisibility_fails_n4_n12": all(not r["divides_B_n"] for r in div_rows if r["n"] >= 4),
            "divisibility_fails_n3_n12": all(not r["divides_B_n"] for r in div_rows),
            "n3_divisibility_exception": next(r for r in div_rows if r["n"] == 3)["divides_B_n"],
            "n3_exact_translate_tiling_exists": tiler_exact_rows[0]["left_translate_tiling_exists"],
        },
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8", newline="\n")


def make_markdown(payload: dict) -> str:
    lines = [
        "# P15 Positive-Only Separator Hunt",
        "",
        f"Schema: `{payload['schema']}`",
        "Status: **RESEARCH-SCOUT**, not a promoted P15 theorem.",
        "",
        "## Convention",
        "",
        f"Signed total order: `{payload['signed_total_order_convention']}`.",
        "",
        "## Signed-Poset Common-Precedence Table",
        "",
        "| n | |X_n| | |B_n| mod |X_n| | common precedence | nonzero common | |L_B(P_X)| | exact signed-poset? |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    div_by_n = {r["n"]: r for r in payload["divisibility_rows_n3_n12"]}
    for r in payload["poset_rows_n3_n7"]:
        d = div_by_n[r["n"]]
        lines.append(
            f"| {r['n']} | {r['X_size_enumerated']} | {d['B_n_mod_X']} | "
            f"{r['common_precedence_size']} | {r['common_precedence_nonzero_size']} | "
            f"{r['LB_PX_size']} | {r['exact_signed_poset']} |"
        )
    lines += [
        "",
        "## Exact Tiler Exception",
        "",
        "| n | |X_n| | |B_n| | target tile count | distinct left translates | tiling exists? |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for r in payload["tiler_exact_rows"]:
        lines.append(
            f"| {r['n']} | {r['X_size']} | {r['group_size']} | {r['target_tile_count']} | "
            f"{r['distinct_left_translates']} | {r['left_translate_tiling_exists']} |"
        )
    lines += [
        "",
        "## Divisibility Table",
        "",
        "| n | |X_n| | |B_n| | |B_n| mod |X_n| | gcd | divides? |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for r in payload["divisibility_rows_n3_n12"]:
        lines.append(
            f"| {r['n']} | {r['X_size']} | {r['B_n_size']} | {r['B_n_mod_X']} | "
            f"{r['gcd_B_n_X']} | {r['divides_B_n']} |"
        )
    lines += [
        "",
        "## Count Formula",
        "",
        payload["count_formula"],
        "The script evaluates this permanent coefficient by a row-by-row bitmask DP truncated to bidegree `(1,1)`.",
        "",
        "## Summary",
        "",
    ]
    for k, v in payload["summary"].items():
        lines.append(f"- `{k}`: `{v}`")
    lines += [
        "",
        "## Interpretation",
        "",
        "The n=4..7 common-precedence data support the P11-style non-signed-poset route: the common relation is empty, so the induced signed poset is the antichain and its extension set is all of `B_n`, while `X_n` is a proper subset.",
        "",
        "The divisibility obstruction is strong from `n=4` onward: `|X_n|` does not divide `|B_n|` for every `4<=n<=12`, so exact translate tiling is cardinality-obstructed throughout that range. The row `n=3` is exceptional: divisibility holds, and this script records representatives for an actual exact left-translate tiling.",
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
    print(make_markdown(payload))


if __name__ == "__main__":
    main()
