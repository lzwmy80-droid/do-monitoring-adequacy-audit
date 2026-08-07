# USGS external-network replication freeze note

Freeze time: 2026-08-02, before endpoint-width computation.  
Status: `FROZEN_EXTERNAL_TRANSFER_STRESS_TEST`.

## Purpose

This analysis is a cross-network transfer stress test, not an environmental
validation study. It asks whether the frozen Maryland DNR/MARACOOS scenario
labels and the cadence-width pattern can be reproduced on an independent
federal monitoring system without recalibration or result-based station
selection.

## Authority, access, and reuse

- Source: U.S. Geological Survey Water Data for the Nation / National Water
  Information System (NWIS).
- Official API landing page:
  <https://api.waterdata.usgs.gov/ogcapi/v0/?f=html>.
- The frozen time-series identities came from the OGC metadata endpoint.
  Historical observations are retrieved from the official NWIS Instantaneous
  Values service at <https://waterservices.usgs.gov/nwis/iv/> because the OGC
  service rate-limited anonymous bulk access during acquisition. Both services
  expose USGS Water Data for the Nation records; the adapter preserves the OGC
  time-series identity alongside the NWIS site/parameter identity.
- USGS states that USGS-authored data and information are generally U.S.
  public-domain works and asks users to credit the agency. The dataset is cited
  using the USGS-recommended persistent DOI
  <https://doi.org/10.5066/F7P55KJN>.
- Continuous observations are high-frequency automated-sensor data; parameter
  code `00300` denotes dissolved oxygen concentration in mg/L.
- Only observations whose legacy NWIS qualifier is exactly `A` enter the
  analysis. USGS defines `A` as approved for publication after processing and
  review. Other or provisional observations, if present, remain documented in
  the raw response but do not enter an endpoint calculation.
- NOAA NERRS/SWMP was considered because its estuarine design and three-stage
  QA/QC are attractive. It was not downloaded because its delivery workflow
  requires personal contact information and acceptance of programme-specific
  citation/acknowledgement conditions. No identity was invented and no access
  condition was bypassed.

## Metadata-only inclusion rule

Query the USGS time-series metadata endpoint for all records satisfying:

1. `parameter_code=00300`;
2. `computation_identifier=Instantaneous`;
3. `primary=Primary`;
4. `state_name=Maryland`;
5. unit `mg/l`;
6. period of record begins no later than `2024-06-01T00:00:00Z` and ends no
   earlier than `2024-09-01T00:00:00Z`.

All matching series are retained in lexicographic order by
`monitoring_location_id` and `time_series_id`. No station is replaced because
its DO values, compatibility, or endpoint widths are inconvenient. Discovery
responses used to check cadence and approval may carry DO values, but values
were not used to define the eligible set. The frozen query yields nine series:

| Monitoring location | Time-series ID |
|---|---|
| USGS-01491000 | 79480fea74964d81bf3822a6a66ee1c7 |
| USGS-01579550 | abaf8aaaf7bd495483198ad33d277024 |
| USGS-01594441 | 979a9c0c3a4842e8bd56bd05b649ab28 |
| USGS-01638500 | 7329a5d59c164090aee93904299d4a7a |
| USGS-01643580 | 4d259ae612cb4d82bd36c7ccd1051434 |
| USGS-01646500 | bb18f3558934465194f8da1d5f1c3df3 |
| USGS-01649190 | 42075f9554bd4ee0b1909d52959482ae |
| USGS-01649500 | 80262e6c200045eb8755f9f1780b799b |
| USGS-01650800 | d4a0194d99aa4663b9c935003d7ed3f0 |

## Frozen observation window and preprocessing

- The NWIS retrieval envelope uses local dates `2024-05-31/2024-09-02` so that
  conversion across daylight-saving offsets cannot truncate either UTC
  boundary. The frozen analysis interval is the half-open UTC interval
  `[2024-06-01T00:00:00Z, 2024-09-01T00:00:00Z)`.
- Preserve every returned raw row and source URL.
- Require agency code `USGS`, parameter code `00300`, the NWIS Instantaneous
  Values statistic code `00000`, unit `mg/l`, unique frozen time-series
  identity, and finite numeric value. The OGC metadata describes the same
  series using `computation_identifier=Instantaneous`.
- Retain only observations with qualifier exactly `A` for analysis.
- Convert timestamps to UTC and sort stably; duplicate UTC timestamps are an
  analysis stop, not silently averaged.
- Build the common 15-min operational reference by retaining only observations
  exactly on UTC quarter hours (`minute` in `00,15,30,45`, `second=0`). For
  native 5-min stations this is deterministic sub-selection, not interpolation
  or averaging.
- Apply the existing seven-day rule unchanged: 673 points, 167.5-168.5 h span,
  positive adjacent time differences no greater than 18 min, earliest
  row-disjoint qualifying windows, maximum four windows per station.
- A selected series that supplies no qualifying block remains in the reported
  denominator. No alternative period, interpolation, or replacement station is
  permitted.

## Frozen scientific parameters

- Demonstration threshold: `H=5 mg/L`.
- Physical lower state bound: `B=0 mg/L`.
- Path label: Maryland DNR calibration q99.9,
  `L=13.203959999999936 mg/L/h`.
- `L` is transferred unchanged. USGS values do not recalibrate, repair, or
  validate it as a physical bound.
- Cadences: `0.5, 1, 2, 4 h`, with every index phase.
- Main external analysis: exact-inlier `k=0`, with and without the nonnegative
  state floor as in the frozen Table 2 pipeline.
- Paired sensitivity: fixed `k=1`, only on the identical cases that are feasible
  at `k=0` and whose complete 15-min reference block is compatible with the
  frozen `L`.

## Required reporting and stop rules

- Report all nine selected series, all qualifying blocks, and the full
  compatible/incompatible denominator.
- Report cadence summaries only conditional on compatible feasible cases;
  overlapping cadence phases are not independent replicates.
- Verify deterministic containment of the retained 15-min operational
  reference, `k=0` nesting within `k=1`, physical width ranges, endpoint witness
  budgets, and source/core-runtime hashes.
- If no external block is compatible, report failure of label transfer; do not
  increase `L`.
- If containment or nesting fails, stop manuscript editing and audit the
  adapter/implementation.
- Do not call this an external validation, population sample, confidence
  analysis, ecological-burden estimate, or proof that `H` is locally regulatory.
- Do not add a third network, new threshold, new model, or Figure 4 to search
  for a more attractive result.
