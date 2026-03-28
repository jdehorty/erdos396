# Validation Certificates

This directory contains the public computational certificate bundle for the
completed `validate` runs for `k=8` through `k=13`.

The full copied raw run directories are archived under `artifacts/validate_runs/`.
This directory is the smaller reviewer-facing subset.

For each certified `k`, the bundle includes:

- the two final `validate_report_*.json` files copied from the completed runs
- one run manifest per validated range, with UTC start time, commit, binary
  SHA-256, exact command line, and the host actually used for that partition
- `SHA256SUMS` for every file in this directory
- one `validate_k*_partition_check.txt` per certified `k`, produced by
  `certificates/scripts/check_validate_report.py --check-partition`
- `validate_k9_witness_1019547844_check.txt`, an independent stdlib-only
  verification of the corrected `k=9` witness found by `validate`

The trust model is the one documented in `docs/trust.md`: these are
computational certificates backed by checked coverage invariants, explicit build
metadata, and independent report verification.

## Certified Results

### `k=8`

Certified minimum witness:

- `339,949,252`

Validated partitions:

- `[0, 169,974,626)`
- `[169,974,626, 339,949,252)`

Binary hashes by partition:

- `[0, 169,974,626)`: `f9723953d896128bf7bcd836dbb897604a2c624d80c75ce492fd751fc2d216b2`
- `[169,974,626, 339,949,252)`: `65980f908f6527106f99e39163b6883dc37dd5ea6a51a30e0d53058367cb58c5`

Embedded report build metadata:

- git commit: `f0dfb4169c5fb726f29753612ca9bd0fea98abad`
- git dirty: `true`
- rustc: `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 169,974,626)`: `checked=169974617`, `screen_pass=8`,
  `false_alarms=8`, `witnesses=0`, `rate_per_sec=46801905.480290964`
- `[169,974,626, 339,949,252)`: `checked=169974626`, `screen_pass=1`,
  `false_alarms=1`, `witnesses=0`, `rate_per_sec=39333781.52707239`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

### `k=9`

Certified minimum witness:

- `1,019,547,844`

Validated partitions:

- `[0, 8,804,882,496)`
- `[8,804,882,496, 17,609,764,993)`

Binary hashes by partition:

- `[0, 8,804,882,496)`: `f9723953d896128bf7bcd836dbb897604a2c624d80c75ce492fd751fc2d216b2`
- `[8,804,882,496, 17,609,764,993)`: `65980f908f6527106f99e39163b6883dc37dd5ea6a51a30e0d53058367cb58c5`

Embedded report build metadata:

- git commit: `f0dfb4169c5fb726f29753612ca9bd0fea98abad`
- git dirty: `true`
- rustc: `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 8,804,882,496)`: `checked=8804882486`, `screen_pass=10`,
  `false_alarms=9`, `witnesses=[1019547844]`, `rate_per_sec=45917486.15994614`
- `[8,804,882,496, 17,609,764,993)`: `checked=8804882497`,
  `screen_pass=0`, `false_alarms=0`, `witnesses=0`,
  `rate_per_sec=32881116.335832387`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

Correction note:

- the March 16, 2026 `validate` pass disproved the earlier public claim that
  `17,609,764,993` was the smallest `k=9` witness
- the first partition report records the corrected witness directly in its
  `witnesses` list
- `validate_k9_witness_1019547844_check.txt` independently rechecks the corrected
  witness with `certificates/scripts/verify_witness.py`

### `k=10`

Certified minimum witness:

- `17,609,764,994`

Validated partitions:

- `[0, 8,804,882,497)`
- `[8,804,882,497, 17,609,764,994)`

Binary hashes by partition:

- `[0, 8,804,882,497)`: `f9723953d896128bf7bcd836dbb897604a2c624d80c75ce492fd751fc2d216b2`
- `[8,804,882,497, 17,609,764,994)`: `65980f908f6527106f99e39163b6883dc37dd5ea6a51a30e0d53058367cb58c5`

Embedded report build metadata:

- git commit: `f0dfb4169c5fb726f29753612ca9bd0fea98abad`
- git dirty: `true`
- rustc: `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 8,804,882,497)`: `checked=8804882486`, `screen_pass=12`,
  `false_alarms=12`, `witnesses=0`, `rate_per_sec=45599604.88285182`
- `[8,804,882,497, 17,609,764,994)`: `checked=8804882497`,
  `screen_pass=0`, `false_alarms=0`, `witnesses=0`,
  `rate_per_sec=34377613.67847653`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

### `k=11`

Certified minimum witness:

- `1,070,858,041,585`

Validated partitions:

- `[0, 610,796,592,505)`
- `[610,796,592,505, 1,070,858,041,585)`

Binary hashes by partition:

- `[0, 610,796,592,505)`: `2fbebf402631f3be5ee8f7d645f3fc6e2eeb2aee0122ac3462d299fcbb3cba32`
- `[610,796,592,505, 1,070,858,041,585)`: `1cf2aa9103307ca5e22c4473110e23542ba4e30c3d8c9ef1f3459f7b0fde36c2`

Embedded report build metadata:

- `[0, 610,796,592,505)`: git commit `f0dfb4169c5fb726f29753612ca9bd0fea98abad`, git dirty `true`, rustc `rustc 1.87.0 (17067e9ac 2025-05-09)`
- `[610,796,592,505, 1,070,858,041,585)`: git commit `f0dfb4169c5fb726f29753612ca9bd0fea98abad`, git dirty `false`, rustc `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 610,796,592,505)`: `checked=610796592493`, `screen_pass=11`,
  `false_alarms=11`, `witnesses=0`, `rate_per_sec=40499234.25234593`
- `[610,796,592,505, 1,070,858,041,585)`: `checked=460061449080`,
  `screen_pass=2`, `false_alarms=2`, `witnesses=0`,
  `rate_per_sec=30054332.64434419`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

Provenance note:

- the server-b partition manifest records `git_dirty: true` at launch, while the
  embedded build metadata in the JSON report records `git_dirty: false`; the
  bundle preserves both records verbatim
- the commit hash, rustc version, and binary SHA-256 provenance are preserved
  in the copied artifacts for independent inspection

### `k=12`

Certified minimum witness:

- `5,048,891,644,646`

Validated partitions:

- `[0, 2,880,387,737,783)`
- `[2,880,387,737,783, 5,048,891,644,646)`

Binary hashes by partition:

- `[0, 2,880,387,737,783)`: `2c892f9ac59087ac1d62a40176b709ef29b5948974c990e6f874578d1db7ebec`
- `[2,880,387,737,783, 5,048,891,644,646)`: `e2563bef0ea904f15d95a622d5fcb4b00ce0254c5ce809ba87e7b230c5fc19b1`

Embedded report build metadata:

- git commit: `f0dfb4169c5fb726f29753612ca9bd0fea98abad`
- git dirty: `false`
- rustc: `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 2,880,387,737,783)`: `checked=2880387737770`, `screen_pass=17`,
  `false_alarms=17`, `witnesses=0`, `rate_per_sec=38810293.49235762`
- `[2,880,387,737,783, 5,048,891,644,646)`: `checked=2168503906863`,
  `screen_pass=0`, `false_alarms=0`, `witnesses=0`,
  `rate_per_sec=29232513.871942077`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

### `k=13`

Certified minimum witness:

- `18,253,129,921,842`

Validated partitions:

- `[0, 9,684,496,742,351)`
- `[9,684,496,742,351, 18,253,129,921,842)`

Binary hashes by partition:

- `[0, 9,684,496,742,351)`: `2c892f9ac59087ac1d62a40176b709ef29b5948974c990e6f874578d1db7ebec`
- `[9,684,496,742,351, 18,253,129,921,842)`: `e2563bef0ea904f15d95a622d5fcb4b00ce0254c5ce809ba87e7b230c5fc19b1`

Embedded report build metadata:

- git commit: `f0dfb4169c5fb726f29753612ca9bd0fea98abad`
- git dirty: `false`
- rustc: `rustc 1.87.0 (17067e9ac 2025-05-09)`

Report summary:

- `[0, 9,684,496,742,351)`: `checked=9684496742337`, `screen_pass=15`,
  `false_alarms=15`, `witnesses=0`, `rate_per_sec=37945728.17781568`
- `[9,684,496,742,351, 18,253,129,921,842)`: `checked=8568633179491`,
  `screen_pass=0`, `false_alarms=0`, `witnesses=0`,
  `rate_per_sec=28567494.392981406`

Independent checker result:

- both reports: `OK`
- partition check: `OK`

## How To Recheck

Run:

```bash
sha256sum -c certificates/validate/SHA256SUMS

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k8_0_169974626.json \
  certificates/validate/validate_report_k8_169974626_339949252.json \
  --check-partition

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k9_0_8804882496.json \
  certificates/validate/validate_report_k9_8804882496_17609764993.json \
  --check-partition

python3 certificates/scripts/verify_witness.py -k 9 -n 1019547844

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k10_0_8804882497.json \
  certificates/validate/validate_report_k10_8804882497_17609764994.json \
  --check-partition

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k11_0_610796592505.json \
  certificates/validate/validate_report_k11_610796592505_1070858041585.json \
  --check-partition

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k12_0_2880387737783.json \
  certificates/validate/validate_report_k12_2880387737783_5048891644646.json \
  --check-partition

python3 certificates/scripts/check_validate_report.py \
  certificates/validate/validate_report_k13_0_9684496742351.json \
  certificates/validate/validate_report_k13_9684496742351_18253129921842.json \
  --check-partition
```

If these checks pass, then under the trust model in `docs/trust.md`, the public
record here certifies minimality for `k=8` through `k=13`. In particular, it
certifies that no `k=8` witness exists below `339,949,252`, that
`1,019,547,844` is a valid `k=9` witness and no smaller `k=9` witness exists,
and that no `k=10`, `k=11`, `k=12`, or `k=13` witness exists below the listed
minimum for that `k`.
