# Supplementary Table S1. Frozen data inventory and protocol

Data were downloaded on 2026-07-16 from the MARACOOS ERDDAP delivery of the
Maryland DNR Continuous Monitoring Program. The frozen download manifest has
SHA-256 `c06cb797bed4541135b0726ccb07ea04df1eb51929f60645cfce60d530733d2d`.
The source is public and the catalog record reports a CC0-1.0 license.

The station-selection rule was fixed before endpoint analysis: sort eligible
MDDNR datasets by dataset identifier, exclude explicit paired Bottom/Surface
titles, search warm seasons from 2025 backward to 2017, retain the most recent
June–August record with at least 5,000 rows, and stop after eight stations.
DO values were not used for station selection. Four station-years were assigned
to slope-label calibration and four to the frozen holdout analysis.

| Role | Dataset ID | Year | Raw rows | Valid DO rows | Valid time span (UTC) |
|---|---|---:|---:|---:|---|
| Calibration | mddnr_Arundel_on_the_Bay | 2019 | 6,677 | 6,676 | 2019-06-20 16:30 to 2019-09-01 00:00 |
| Calibration | mddnr_Budds_Landing | 2018 | 10,275 | 10,275 | 2018-06-01 00:00 to 2018-09-01 00:00 |
| Calibration | mddnr_Dares_Beach | 2017 | 8,017 | 7,025 | 2017-06-01 00:00 to 2017-08-23 15:00 |
| Calibration | mddnr_Harris_Creek_Downstream | 2018 | 6,526 | 6,526 | 2018-06-01 00:00 to 2018-09-01 00:00 |
| Holdout | mddnr_Bishopville_Prong | 2018 | 8,833 | 8,833 | 2018-06-01 00:00 to 2018-09-01 00:00 |
| Holdout | mddnr_Camp_Tockwogh | 2017 | 8,832 | 8,832 | 2017-06-01 00:00:14 to 2017-08-31 23:45:14 |
| Holdout | mddnr_Greys_Creek | 2018 | 8,833 | 8,265 | 2018-06-01 00:00 to 2018-09-01 00:00 |
| Holdout | mddnr_Harris_Creek_Upstream | 2018 | 8,779 | 7,636 | 2018-06-01 00:00 to 2018-09-01 00:00 |

## Preprocessing and frozen parameters

- Variable: `mass_concentration_of_oxygen_in_sea_water`, expressed in mg/L.
- Published QA/QC values were used; missing/non-numeric DO values were not
  imputed.
- Stream identity was determined without examining DO values, using the frozen
  sample-depth proximity rule in the manifest and pilot metadata.
- The calibration-only slope labels used 28,968 eligible adjacent slopes:
  q99.9 \(L=13.20396\) mg/L/h and observed maximum \(L=25.52\) mg/L/h.
  These are sensitivity labels, not validated physical path bounds.
- \(H=5\) mg/L was fixed as a policy-relevant demonstration threshold. A
  5 mg/L criterion appears in Maryland rules, but designated-use criteria vary;
  \(H\) is not presented as a universal ecological threshold for every site
  and time.
- The physical analysis fixes \(x(t)\ge0\) mg/L.
- Figure 2 uses the earliest complete 24-h, 15-min block at each holdout
  station, seven frozen \(L\) values, and station-specific \(k_{\min}(L)\).
- Table 2 is a separate \(k=0\) sampling baseline. It uses up to four earliest
  row-disjoint seven-day blocks per holdout station, all phases at 0.5, 1, 2,
  and 4 h cadence, and the calibration q99.9 slope label.
- Of 480 pooled station–window–cadence–phase cases, 468 were algorithm-feasible
  and 240 belonged to the reference-compatible feasible subset. Those 240
  cases arose from 8/16 blocks and 2/4 holdout stations and are not independent
  windows.

## Source links

- MARACOOS ERDDAP search:
  <https://erddap.maracoos.org/erddap/search/index.json?page=1&itemsPerPage=1000&searchFor=mddnr>
- Maryland DNR continuous-monitoring station table:
  <https://www.eyesonthebay.dnr.maryland.gov/eyesonthebay/ConMonStationTable.cfm>
- Maryland designated-use DO criteria:
  <https://dsd.maryland.gov/regulations/Pages/26.08.02.03-3.aspx>
