# Data sources and legal boundaries

## MDDNR/MARACOOS

The eight station-year files under
`04_真实数据/MARACOOS_MDDNR_discovery_v3` were obtained from the public
MARACOOS ERDDAP delivery of the Maryland Department of Natural Resources
Continuous Monitoring Program. The associated Data.gov catalog records identify
the data as public and CC0-1.0. Exact station catalog pages, frozen query URLs,
and hashes are recorded in the manifest and Supplementary Table S1.

- Example catalog record: https://catalog.data.gov/dataset/mddnr-station-bishopville-prong
- MARACOOS ERDDAP: https://erddap.maracoos.org/erddap/

## USGS Water Data for the Nation / NWIS

The nine series under `04_真实数据/USGS_WDFN_external_2024` were retrieved
from the official USGS NWIS Instantaneous Values service after identities were
frozen from official Water Data for the Nation metadata. Raw responses and
lossless normalized adapters are both retained. Only observations with legacy
qualifier exactly `A` (approved for publication after processing and review)
entered analysis.

- Dataset citation DOI: https://doi.org/10.5066/F7P55KJN
- USGS citation guidance: https://waterdata.usgs.gov/citation/
- USGS public-domain guidance: https://www.usgs.gov/faqs/are-usgs-reportspublications-copyrighted
- USGS data-licensing guidance: https://www.usgs.gov/data-management/data-licensing

The author-owned code licence does not change either source-data status.
Redistribution does not imply endorsement by MDDNR, MARACOOS, NOAA, Data.gov,
USGS, or the U.S. Government.
