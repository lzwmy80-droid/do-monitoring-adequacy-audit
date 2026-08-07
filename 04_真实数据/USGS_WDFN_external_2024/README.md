# USGS Water Data for the Nation external-network cohort

Status: frozen local research inputs for the DIR-01 cross-network transfer
stress test. Accessed 2 August 2026 (local date; acquisition timestamp is also
recorded in UTC in `download_manifest.json`).

## Authority and citation

- Data authority: U.S. Geological Survey (USGS), Water Data for the Nation /
  National Water Information System (NWIS).
- Parameter: dissolved oxygen, USGS code `00300`, unit `mg/l`.
- Metadata discovery: official USGS Water Data OGC time-series endpoint.
- Historical values: official NWIS Instantaneous Values service.
- Recommended citation: U.S. Geological Survey (2026), *U.S. Geological Survey
  National Water Information System database*, accessed 2 August 2026,
  <https://doi.org/10.5066/F7P55KJN>.
- USGS citation guidance: <https://waterdata.usgs.gov/citation/>.
- USGS public-domain guidance:
  <https://www.usgs.gov/faqs/are-usgs-reportspublications-copyrighted>.

USGS source records are not relicensed as author code. Any future MIT licence
in this project applies only to code for which the author owns the copyright.

## Frozen cohort and processing boundary

The cohort contains every one of the nine Maryland primary instantaneous DO
series identified by the pre-computation metadata rule in
`USGS_EXTERNAL_REPLICATION_FREEZE_NOTE_20260802.md`. Stations were not removed
or substituted because of DO values, slope compatibility, endpoint width, or
the presence or absence of low oxygen.

Raw JSON responses are retained byte-for-byte under `raw/`. The `normalized/`
CSV files are lossless adapters that preserve the source timestamp, UTC
timestamp, value, qualifier, method identifier, frozen OGC time-series ID and
NWIS series identity. `download_manifest.json` records exact queries and
SHA-256 hashes for both representations.

Only records with legacy NWIS qualifier exactly `A` enter the analysis. USGS
defines `A` as approved for publication after processing and review. This
status does not make the observations a latent continuous truth. Native 5-min
series are deterministically subselected at exact UTC quarter hours; no value
is averaged or interpolated.

## Reproduction

From the project root:

```text
python 03_代码与测试/current/download_usgs_external_do.py
python 03_代码与测试/current/run_usgs_external_replication.py
python 03_代码与测试/current/verify_usgs_external_replication.py
```

The downloader reuses existing raw responses by default, preserving the frozen
snapshot. Removing or replacing the raw snapshot is a new acquisition and
must not be represented as byte-identical reproduction of this cohort.
