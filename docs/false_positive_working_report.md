# False-Positive Working Report

This is a working note based on the audited Parquet corpus in
`data/run_corpus/v1/`. It is intentionally organized so that we can decide
 later which parts belong in the main writeup and which parts are better kept
 as appendix or internal notes.

Full machine-readable data is in:
`data/run_corpus/v1/analysis/false_positive_analysis.json`

## Scope

- Corpus: canonical unique-coverage run corpus up to `25T`
- Audited windows: exact run lengths `6..14`
- Total audited windows: `72,646,980`
- Total false positives: `338,322`
- Overall false-positive rate: `0.4657%`

## Likely Keep

### 1. False positives are overwhelmingly a small-prime phenomenon

Overall failing-prime counts:

| Prime | Count | Share |
|---|---:|---:|
| `3` | `149,606` | `44.22%` |
| `5` | `134,466` | `39.75%` |
| `2` | `52,327` | `15.47%` |
| `7` | `1,922` | `0.57%` |
| `11` | `1` | `<0.01%` |

So `{2,3,5}` explain `99.43%` of all false positives.

### 2. Most false positives are near-misses

Deficit means `demand - supply` at the first failing prime.

| Deficit | Count | Share |
|---|---:|---:|
| `1` | `308,119` | `91.07%` |
| `2` | `24,513` | `7.25%` |
| `3` | `5,093` | `1.51%` |
| `4+` | `597` | `0.18%` |

This strongly supports the idea that most false positives are "one-slack"
 failures rather than large structural collapses.

### 3. Longer runs are proportionally more fragile

| Run length | Audited | False positives | Rate |
|---|---:|---:|---:|
| `6` | `64,416,047` | `263,229` | `0.4086%` |
| `7` | `7,310,113` | `61,644` | `0.8433%` |
| `8` | `818,645` | `11,340` | `1.3852%` |
| `9` | `91,158` | `1,746` | `1.9154%` |
| `10` | `9,877` | `309` | `3.1285%` |
| `11` | `1,035` | `47` | `4.5411%` |
| `12` | `93` | `5` | `5.3763%` |
| `13` | `11` | `2` | `18.1818%` |
| `14` | `1` | `0` | `0%` |

Absolute counts are dominated by short runs, but conditional fragility rises
 with run length.

### 4. Nearly all false positives live in the short-run regime

- Lengths `6` and `7` account for `324,873 / 338,322 = 96.03%`
- Lengths `6`, `7`, and `8` account for `336,213 / 338,322 = 99.38%`
- Lengths `9..14` account for only `2,109 / 338,322 = 0.62%`

This matters because it means most of the false-positive story is already
 visible in the short-run regime.

## Strong Exploratory Signal

### 5. The same obstruction often persists across adjacent lengths

- Unique false-positive `n` values: `301,182`
- Repeated false-positive `n` values: `33,064`

Most common repeated length patterns:

| Length pattern | Count |
|---|---:|
| `[6,7]` | `25,893` |
| `[7,8]` | `3,028` |
| `[6,7,8]` | `2,873` |
| `[8,9]` | `427` |
| `[7,8,9]` | `337` |
| `[6,7,8,9]` | `328` |

This is evidence that false positives are not isolated accidents. A single
 arithmetic obstruction often survives when the window is extended by one term.

### 6. The `p=3` failures show a sharp residue-class pattern mod `9`

Overall `p=3` false-positive hotspots:

| Residue mod `9` | Audited | `p=3` false positives | Rate | Lift vs baseline |
|---|---:|---:|---:|---:|
| `0` | `8,726,227` | `40,201` | `0.4607%` | `2.24x` |
| `3` | `8,738,052` | `40,154` | `0.4595%` | `2.23x` |
| `4` | `8,596,679` | `32,710` | `0.3805%` | `1.85x` |
| `1` | `8,604,937` | `32,706` | `0.3801%` | `1.85x` |
| `6` | `7,151,977` | `3,436` | `0.0480%` | `0.23x` |
| `7` | `6,852,241` | `355` | `0.0052%` | `0.03x` |
| `2` | `8,585,031` | `44` | `0.0005%` | `0.00x` |
| `5` | `8,569,308` | `0` | `0%` | `0x` |
| `8` | `6,822,528` | `0` | `0%` | `0x` |

The striking fact is not just that some residues are elevated. It is that some
 residues are essentially absent. That looks much more like a real modular law
 than random noise.

### 7. The `p=5` failures are heavily concentrated in residue classes mod `25`

Overall `p=5` false-positive hotspots:

| Residue mod `25` | Audited | `p=5` false positives | Rate | Lift vs baseline |
|---|---:|---:|---:|---:|
| `0` | `4,193,939` | `57,665` | `1.3750%` | `7.43x` |
| `5` | `4,190,931` | `57,120` | `1.3629%` | `7.36x` |
| `6` | `2,725,215` | `6,002` | `0.2202%` | `1.19x` |
| `1` | `3,423,012` | `6,176` | `0.1804%` | `0.97x` |
| `10` | `3,203,915` | `5,144` | `0.1606%` | `0.87x` |

The main effect is very concentrated: residues `0` and `5` mod `25` dominate
 the `p=5` obstruction.

## Possible Interpretation

The computational picture is unusually coherent:

- False positives are almost always caused by very small primes.
- They are usually short by exactly one valuation unit.
- The same `n` frequently fails at adjacent lengths.
- For `p=3` and `p=5`, the failures cluster in very specific residue classes.

Taken together, this looks much more like a small-prime boundary law than a
 diffuse or chaotic phenomenon.

## Probably Appendix / Internal Only

- Full per-length residue tables
- Per-trillion bucket counts
- Rare-prime outliers such as the single `p=11` case
- Long lists of repeated examples

These are useful for internal checking, but probably too detailed for the main
 narrative unless we later discover a specific theorem-level use for them.
