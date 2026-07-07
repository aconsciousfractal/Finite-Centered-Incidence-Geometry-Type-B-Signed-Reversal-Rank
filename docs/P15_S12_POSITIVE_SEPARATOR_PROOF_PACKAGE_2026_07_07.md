# P15 S12 Positive-Only Separator Proof Package

Status: **S12-PROOF-PACKAGE / GATE A CLOSED / GATE B CLOSED / RED-TEAM FIXES APPLIED**.

Object, in the S12 one-indexed word convention:

`B_n = {w=(w_1,...,w_n): {|w_1|,...,|w_n|}={1,...,n}}`,

`rho(i)=n+1-i`,

`X_n = B~I_n(1,1) = {w in B_n : #{i : w_i=i}=1 and #{i : w_i=rho(i)}=1}`.

This is the same object as the zero-indexed manuscript convention after replacing `i=1,...,n` by `i=0,...,n-1` and `rho(i)=n+1-i` by `rho(i)=n-1-i`.

The S12 target is a P11-style add-on theorem:

> For `n>=4`, the positive-only signed-reversal layer `X_n` is neither a Reiner-style signed-poset linear-extension set nor an exact left-translate tiler of `B_n`.

The `n=3` row is explicitly exceptional for tiling: `X_3` has an exact left-translate tiling with `24` tiles, recorded in the Gate B certificate.

## 1. Signed-Poset Separator

For a signed permutation `w=(w_1,...,w_n)`, use the Reiner-style total order

`T(w) = -w_n < ... < -w_1 < 0 < w_1 < ... < w_n`.

Let `Sigma_n={-n,...,-1,0,1,...,n}`. For `X subset B_n`, define the common signed-precedence relation

`C(X) = intersection_{w in X} {(a,b) in Sigma_n^2 : a precedes b in T(w)}`.

If `X=L_B(P)` for a Reiner-style signed poset `P`, every strict relation of `P`, and every relation forced by transitive closure, must occur in `C(X)`. Therefore, if `C(X)=empty`, then `P` is the antichain and `L_B(P)=B_n`. Since `X_n` is a proper subset of `B_n`, this proves the signed-poset obstruction.

### Lemma S12.1, Witness Form

For every `n>=4` and every ordered pair `a,b in Sigma_n` with `a != b`, there exists `w in X_n` such that `b` precedes `a` in `T(w)`.

This says that no signed-precedence relation, including a relation with `0`, is common to all elements of `X_n`.

### Proof Route

For `n>=5`, the key is the sign-extension lemma:

**Sign-extension lemma.** Any prescribed signs on at most two absolute labels can be extended to an element of `X_n`.

Given ordered signed labels `a,b`:

- if `a,b` are nonzero and `a != -b`, prescribe label `|a|` with the sign of `a` and label `|b|` with the sign opposite to `b`; then `b` lies before `0` and `a` lies after `0` in `T(w)`;
- if `a=-b`, prescribe label `|a|` with the sign of `a`; then `b=-a` lies before `0` and `a` lies after `0`;
- if one of `a,b` is `0`, prescribe the one nonzero label with the sign that puts it on the required side of `0`.

The sign-extension construction chooses the required positive fixed/reversal hits away from negatively prescribed labels. If the odd center is used as an auxiliary hit, it must be the intended center hit; auxiliary noncenter hits are chosen away from the center to avoid creating both a fixed and a reversal hit. This center exclusion is the red-team fix applied to the Gate A proof note.

The case `n=4` is finite. The exact witness-coverage certificate kills all `72` ordered pairs in `Sigma_4`, including zero-relations.

### Current Evidence

`P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.md` records:

| item | value |
|---|---:|
| `n4_ordered_nonzero_pairs` | 56 |
| `n4_ordered_pairs_with_zero` | 72 |
| `n4_all_precedences_with_zero_killed` | True |
| `n4_X_size` | 40 |

Verdict for this half: **closed** by `P15_S12_GATE_A_SIGNED_POSET_PROOF_2026_07_07.md` and `P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.json`.

## 2. Translate-Tiler Divisibility Obstruction

An exact left-translate tiling of `B_n` by `X_n` is a partition of `B_n` by sets `gX_n`. Such a tiling requires `|X_n|` to divide `|B_n|=2^n n!`.

The case `n=3` is a real exception: `|X_3|=2` divides `|B_3|=48`, and the Gate B certificate records an exact left-translate tiling with `24` representative translates. Therefore every non-tiler theorem claim starts at `n>=4`.

## 3. Exact Rook Count

Let `N_n=|X_n|`. The rook polynomial for selecting compatible positive fixed/reversal events is

`R_{2m}(A,B,Z)=(1+2AZ+2BZ+A^2Z^2+B^2Z^2)^m`,

`R_{2m+1}(A,B,Z)=R_{2m}(A,B,Z)(1+AZ+BZ+ABZ)`.

Here `A` marks selected fixed-positive events, `B` marks selected reversal-positive events, and `Z` marks distinct forced cells. The center term `ABZ` is the odd case where one positive center cell contributes both events.

By inclusion-exclusion for exactly one fixed-positive and one reversal-positive event,

`N_n = sum_{p,q>=1,u} (-1)^(p+q) p q [A^p B^q Z^u] R_n(A,B,Z) 2^(n-u)(n-u)!`.

The minicertificate verifies this rook count against the independent DP count for `n=3..12`.

## 4. Differentiated Generators

For quotient estimates, divide by `G_n=2^n n!`. Let

`S_n(Z) = (A d/dA)(B d/dB) R_n(A,B,Z) |_{A=B=-1}`.

Then

`N_n/G_n = sum_u [Z^u] S_n(Z) / (2^u (n)_u)`.

The differentiated forms are compact.

For `n=2m`,

`S_{2m}(Z)=4m(m-1) Z^2(1-Z)^2(1-4Z+2Z^2)^(m-2)`.

For `n=2m+1`,

`S_{2m+1}(Z)=Z(1-4Z+2Z^2)^m + 4m(m-1)Z^2(1-Z)^3(1-4Z+2Z^2)^(m-2)`.

These formulas are the starting point for the all-`n` quotient-window proof.

## 5. Quotient Windows

The quotient windows are:

| range | target window |
|---|---|
| `n=4` | `9 < G_n/N_n < 10` |
| odd `n>=5` | `10 < G_n/N_n < 11` |
| even `6<=n<=20` | `11 < G_n/N_n < 12` |
| even `n>=22` | `10 < G_n/N_n < 11` |

This proves non-divisibility whenever the windows are established, because the quotient lies strictly between consecutive integers.

The old heuristic `11 < G_n/N_n < 12` for all even `n` is false. At `n=22`, the exact floor is already `10`.

## 6. Infinite Tail Proof

For even `n=2m>=22`, substitute `Z=-x` and factor

`1-4Z+2Z^2 = (1+(2+sqrt(2))x)(1+(2-sqrt(2))x)`.

The coefficients become positive elementary symmetric functions after the sign change. Newton ultra-log-concavity gives decreasing normalized terms, so alternating partial sums bound the quotient. Symbolic truncation proves `L_9>1/11` and `U_6<1/10` for all `m>=11`.

For odd `n=2m+1>=23`, the proof writes

`S_{2m+1}(-x)=x(1+4x+2x^2)^(m-2) Q_m(x)`.

The sign-pattern proof treats the top support separately: for degrees `2m-1,2m,2m+1`, positivity comes from the positive `q_2,q_3,q_4` terms, not from the lower two-term bound. The normalized monotonicity is then reduced to positivity of

`H_m(x)=2((2m+1)C_m(x)-xC_m'(x))-C_m(x)/x`,

and the Gate B proof note gives the coefficient proof using Newton/log-concavity and the explicit polynomial `P_m(x)`.

Verdict for this half: **closed** by `P15_S12_GATE_B_TILER_PROOF_2026_07_07.md` and `P15_S12_GATE_B_TILER_CERTIFICATE.json`.

## 7. Gates

Gate A, signed-poset: **closed**. The sign-extension proof is written with the odd-center auxiliary exclusion, and `n=4` is handled by the finite witness certificate including zero-relations.

Gate B, tiler: **closed**. Finite checks, the `n=3` exact tiling witness, rook formula, symbolic truncation inequalities, and tail monotonicity proof are recorded in `P15_S12_GATE_B_TILER_PROOF_2026_07_07.md` and `P15_S12_GATE_B_TILER_CERTIFICATE.json`.

S12 is now ready for a final post-fix red-team scan as a P11-style theorem block for `n>=4`, with `n=3` explicitly recorded as the tiling exception.
