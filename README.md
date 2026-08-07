# Reproducibility package v1.1.0

This package supports the manuscript **Auditing dissolved-oxygen monitoring
adequacy under sparse record uncertainty: low-oxygen duration and cumulative
deficit**.

## Release status

Version 1.1.0 is publicly released at
<https://github.com/lzwmy80-droid/do-monitoring-adequacy-audit/releases/tag/v1.1.0>.
Submission documents are not included. No persistent DOI has been assigned to
this GitHub release.

## Evidence included

- the frozen exact endpoint runtime and 63 unit/adversarial tests;
- the eight-station MDDNR/MARACOOS analysis and corrected preprocessing audit;
- all nine metadata-eligible Maryland USGS dissolved-oxygen series, including
  official NWIS raw JSON, lossless normalized CSV, source URLs, qualifiers,
  identifiers, and SHA-256 hashes;
- the frozen USGS cross-network transfer stress test (36 station-blocks, 1,080
  `k=0` phase cases, and 1,080 fixed-`k=1` cases);
- revised block-first Table 2, external Table 3, and Supplementary Tables S1,
  S4, and S5;
- immutable, timing-scrubbed reference outputs for automated comparison.

The UTF-8 directory names are retained because the frozen scripts derive paths
from this layout. Python and Git support these names on Windows, macOS, and
Linux.

## Quick verification

Python 3.12 is recommended.

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.venv\Scripts\python -m pip install -r requirements-lock.txt
.venv\Scripts\python reproduce.py
```

macOS/Linux:

```bash
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python reproduce.py
```

The quick route runs 63 tests, validates package hashes, verifies the 17 frozen
MDDNR source objects and 18 frozen USGS source objects, independently recomputes
the MDDNR calibration slope, and rechecks the frozen USGS denominators,
containment, phase geometry, and `k=0`-within-`k=1` nesting.

## Full reproduction

```bash
python reproduce.py --full
```

The full route recomputes the MDDNR analyses and assets, recomputes the USGS
cross-network stress test from the included frozen source snapshot, regenerates
the revision-O tables, and compares numerical outputs against immutable
references. Timing-only fields are excluded from the release outputs and from
scientific equality checks. Allow about 25--35 minutes on a modern desktop; the
USGS portion accounted for about 21 minutes on the development machine.

No network access is used by either route. The official-source adapter
`download_usgs_external_do.py` is included for provenance and an independently
refreshed snapshot, but it is deliberately not called during frozen
reproduction because a new access timestamp would define a new snapshot.

## Scientific boundaries

The empirical results are conditional monitoring stress tests. The 15-minute
record is an operational reference, not latent continuous truth; the transferred
slope is a declared scenario label, not a physical guarantee; and an
endpoint-specific retained/omitted certificate is an optimization witness, not
a sensor-error label. The USGS result is a cross-network transportability stress
test within Maryland, not external validation or population inference.

## Licence boundaries

The MIT License in `LICENSE` applies only to original author-owned software in
this repository, copyright 2026 Ziwen Luo. It does not relicense source data or
third-party software. See `DATA_LICENSE.md` for the MDDNR/MARACOOS CC0-1.0 and
USGS U.S. federal public-domain boundaries, and `THIRD_PARTY_NOTICES.md` for
software dependencies.
