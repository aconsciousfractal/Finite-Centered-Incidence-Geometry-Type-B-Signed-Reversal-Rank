# P15-S11 Bounded Raw Rank Smoke Check Certificate

Date: 2026-07-07
Status: PASS

This verifier does not import S3 or S7. It rebuilds signed permutations, diagonal layers, channel matrices, standard projection, and exact ranks over Q for a bounded set of rows.

## Rank Rows

```text
n=3 k=1 size=2 ranks=(2,1) expected=(2,1) ok=True
n=4 k=1 size=40 ranks=(2,1) expected=(2,1) ok=True
n=4 k=2 size=2 ranks=(2,1) expected=(2,1) ok=True
n=5 k=1 size=364 ranks=(3,2) expected=(3,2) ok=True
n=5 k=2 size=42 ranks=(3,2) expected=(3,2) ok=True
n=5 k=3 size=2 ranks=(3,2) expected=(3,2) ok=True
n=6 k=1 size=4080 ranks=(3,2) expected=(3,2) ok=True
n=6 k=2 size=132 ranks=(3,2) expected=(3,2) ok=True
```

## Nonempty Boundary Replay

```text
n=2 k=1 size=0 expected_nonempty=False observed_nonempty=False ok=True
n=2 k=2 size=0 expected_nonempty=False observed_nonempty=False ok=True
n=3 k=1 size=2 expected_nonempty=True observed_nonempty=True ok=True
n=3 k=2 size=0 expected_nonempty=False observed_nonempty=False ok=True
n=4 k=1 size=40 expected_nonempty=True observed_nonempty=True ok=True
n=4 k=2 size=2 expected_nonempty=True observed_nonempty=True ok=True
n=4 k=3 size=0 expected_nonempty=False observed_nonempty=False ok=True
n=5 k=1 size=364 expected_nonempty=True observed_nonempty=True ok=True
n=5 k=2 size=42 expected_nonempty=True observed_nonempty=True ok=True
n=5 k=3 size=2 expected_nonempty=True observed_nonempty=True ok=True
n=6 k=1 size=4080 expected_nonempty=True observed_nonempty=True ok=True
n=6 k=2 size=132 expected_nonempty=True observed_nonempty=True ok=True
n=6 k=3 size=0 expected_nonempty=False observed_nonempty=False ok=True
n=6 k=4 size=0 expected_nonempty=False observed_nonempty=False ok=True
```

## Boundary

This is a compact smoke check only. The theorem still relies on the S4/S5/S8/S9/S9C proof chain and the S11 manuscript lemmas.
