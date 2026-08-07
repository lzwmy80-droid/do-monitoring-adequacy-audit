# Table S4. USGS external-network source and preprocessing audit

The cohort was frozen from metadata before endpoint computation. All raw responses and adapted files are hash-locked. Only the legacy NWIS qualifier exactly `A` (approved for publication after processing and review) entered the frozen UTC interval. Native 5-min records were subselected at exact UTC quarter hours without averaging or interpolation. The NWIS instantaneous series statistic code is `00000`.

| USGS site | Native gap (min) | Raw response rows | Frozen-interval rows | Approved-A rows | Quarter-hour rows | 7-day blocks | DO range (mg/L) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `USGS-01491000` | 15 | 9,113 | 8,825 | 8,825 | 8,825 | 4 | 4.40–9.90 |
| `USGS-01579550` | 5 | 26,973 | 26,349 | 26,349 | 8,783 | 4 | 4.60–12.00 |
| `USGS-01594441` | 5 | 27,329 | 26,465 | 26,465 | 8,820 | 4 | 5.00–8.70 |
| `USGS-01638500` | 15 | 9,104 | 8,821 | 8,821 | 8,821 | 4 | 4.50–12.30 |
| `USGS-01643580` | 15 | 9,080 | 8,796 | 8,796 | 8,796 | 4 | 4.90–12.70 |
| `USGS-01646500` | 15 | 9,099 | 8,811 | 8,811 | 8,811 | 4 | 4.20–16.30 |
| `USGS-01649190` | 5 | 27,259 | 26,395 | 26,395 | 8,795 | 4 | 6.70–9.60 |
| `USGS-01649500` | 5 | 27,318 | 26,454 | 26,454 | 8,817 | 4 | 2.20–15.20 |
| `USGS-01650800` | 5 | 27,218 | 26,355 | 26,355 | 8,778 | 4 | 5.90–9.80 |

## Exact URLs and checksums

- `USGS-01491000` — CHOPTANK RIVER NEAR GREENSBORO, MD; time-series ID `79480fea74964d81bf3822a6a66ee1c7`; statistic `00000`: raw SHA-256 `605d286662a556aacd5d64f70f80f0acc1afd12465ce1205bc88d59497a08443`; normalized SHA-256 `1338ee9e3b241654140b47d2e71f9225d672f587a7bd529fce10f8ee7ec248a6`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01491000&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01579550` — SUSQUEHANNA RIVER NEAR DARLINGTON, MD; time-series ID `abaf8aaaf7bd495483198ad33d277024`; statistic `00000`: raw SHA-256 `a51bd667dceca5ca762da609790d8b28ae9f286552fa462e90a8da7ec033f99b`; normalized SHA-256 `7369c39c86c8964fdeb63702488067a102b590ded7401617dcbf1f78dd48dce4`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01579550&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01594441` — Patuxent R at Gov Bridge Natural Area nr Bowie MD; time-series ID `979a9c0c3a4842e8bd56bd05b649ab28`; statistic `00000`: raw SHA-256 `5e56e834801e596c8c5e7e7849ca58f61960f18310c10064512f7a471209f0ff`; normalized SHA-256 `3da92e49f980bd01bdcbe4eda77ad5ef694c79bd006f21bb2de995a9446b802d`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01594441&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01638500` — POTOMAC RIVER AT POINT OF ROCKS, MD; time-series ID `7329a5d59c164090aee93904299d4a7a`; statistic `00000`: raw SHA-256 `5370de8689d1e44f943caaee53e660ded93425e878239d32bcc1d784d47923bc`; normalized SHA-256 `57759408e58aef2df124e2f99520c686fb8936154ccd1baf6e62b2511f2dc444`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01638500&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01643580` — MONOCACY RIVER NEAR DICKERSON, MD; time-series ID `4d259ae612cb4d82bd36c7ccd1051434`; statistic `00000`: raw SHA-256 `af3cbbe8dc4bdf6ee30e6a2e6551360c4246e3c5cffbf05e99d974f496b2ef05`; normalized SHA-256 `0e130483674a81d3bef4db258cf9173ab11607d979aef25cedc2a36eab4296b1`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01643580&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01646500` — POTOMAC RIVER NEAR WASH, DC LITTLE FALLS PUMP STA; time-series ID `bb18f3558934465194f8da1d5f1c3df3`; statistic `00000`: raw SHA-256 `d38041e4997c67382c1d99d621c63889771c068b9c41055168bee50a10f407f9`; normalized SHA-256 `9aabc1f272c2bcaea4f648579f56b8dee659e0491cb968c597252059e227d290`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01646500&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01649190` — PAINT BRANCH NEAR COLLEGE PARK, MD; time-series ID `42075f9554bd4ee0b1909d52959482ae`; statistic `00000`: raw SHA-256 `a2f69437fda6e9a5f2baf45b0f5fdf07e17936b66a3602c90d0ba9634a93cea0`; normalized SHA-256 `6e329662f78cbc3773154188722052e2bc7b8a0ef9abd45a7eac6d932508e2b4`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01649190&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01649500` — NORTHEAST BRANCH ANACOSTIA RIVER AT RIVERDALE, MD; time-series ID `80262e6c200045eb8755f9f1780b799b`; statistic `00000`: raw SHA-256 `3636788bde8c8dcd3a3d0ea983e267fac855de52e570a6ba6fd9556b63db7460`; normalized SHA-256 `289657456208bedf056052e6252f033ace749c8de78df5826e378e5128679132`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01649500&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.
- `USGS-01650800` — SLIGO CREEK NEAR TAKOMA PARK, MD; time-series ID `d4a0194d99aa4663b9c935003d7ed3f0`; statistic `00000`: raw SHA-256 `04e4e304152f3f1873455a438e57e69a43355f22dde0eb74fd627930bdf2ac9f`; normalized SHA-256 `2860e1e70ab28d3cab5e8d88f7d6f7712edcc32518909936256f7676c56cc0e5`; query <https://waterservices.usgs.gov/nwis/iv/?format=json&sites=01650800&parameterCd=00300&startDT=2024-05-31&endDT=2024-09-02&siteStatus=all>.

Recommended source citation: U.S. Geological Survey (2026), U.S. Geological Survey National Water Information System database, accessed 2 August 2026, <https://doi.org/10.5066/F7P55KJN>.
