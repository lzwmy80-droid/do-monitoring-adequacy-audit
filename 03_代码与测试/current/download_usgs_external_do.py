# -*- coding: utf-8 -*-
"""Download the pre-frozen USGS dissolved-oxygen replication inputs.

The eligible time-series identities were frozen from the USGS Water Data OGC
time-series metadata endpoint before endpoint computation.  Historical values
are retrieved from the official USGS NWIS Instantaneous Values web service,
which exposes the same published Water Data for the Nation observations and
the publication qualifier.  Raw JSON is retained verbatim; a lossless tabular
adapter and a hash manifest are generated for analysis and audit.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUTPUT = PROJECT / "04_真实数据" / "USGS_WDFN_external_2024"
RAW = OUTPUT / "raw"
NORMALIZED = OUTPUT / "normalized"

# Values were not used to choose or order these series.  The identities are
# the complete frozen Maryland metadata result documented in the freeze note.
SERIES = (
    ("01491000", "79480fea74964d81bf3822a6a66ee1c7"),
    ("01579550", "abaf8aaaf7bd495483198ad33d277024"),
    ("01594441", "979a9c0c3a4842e8bd56bd05b649ab28"),
    ("01638500", "7329a5d59c164090aee93904299d4a7a"),
    ("01643580", "4d259ae612cb4d82bd36c7ccd1051434"),
    ("01646500", "bb18f3558934465194f8da1d5f1c3df3"),
    ("01649190", "42075f9554bd4ee0b1909d52959482ae"),
    ("01649500", "80262e6c200045eb8755f9f1780b799b"),
    ("01650800", "d4a0194d99aa4663b9c935003d7ed3f0"),
)

START_DATE = "2024-05-31"
END_DATE = "2024-09-02"
ANALYSIS_START_UTC = "2024-06-01T00:00:00Z"
ANALYSIS_END_UTC = "2024-09-01T00:00:00Z"
PARAMETER_CODE = "00300"
EXPECTED_UNIT = "mg/l"
EXPECTED_STATISTIC = "00000"
BASE_URL = "https://waterservices.usgs.gov/nwis/iv/"
OGC_METADATA_URL = (
    "https://api.waterdata.usgs.gov/ogcapi/v0/collections/"
    "time-series-metadata/items?f=json&parameter_code=00300&"
    "computation_identifier=Instantaneous&primary=Primary&"
    "state_name=Maryland&limit=1000"
)
CITATION_URL = "https://waterdata.usgs.gov/citation/"
PUBLIC_DOMAIN_URL = (
    "https://www.usgs.gov/faqs/are-usgs-reportspublications-copyrighted"
)
USER_AGENT = (
    "DIR01-USGS-external-replication/1.0 "
    "(scholarly reproducibility; contact supplied upon publication)"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, attempts: int = 5) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504}:
                raise
            if attempt == attempts - 1:
                raise
            retry_after = error.headers.get("Retry-After")
            pause = float(retry_after) if retry_after else 2.0 ** (attempt + 1)
            time.sleep(min(pause, 30.0))
        except urllib.error.URLError:
            if attempt == attempts - 1:
                raise
            time.sleep(min(2.0 ** (attempt + 1), 30.0))
    raise RuntimeError("unreachable retry state")


def _query_url(site: str) -> str:
    query = urllib.parse.urlencode(
        {
            "format": "json",
            "sites": site,
            "parameterCd": PARAMETER_CODE,
            "startDT": START_DATE,
            "endDT": END_DATE,
            "siteStatus": "all",
        }
    )
    return f"{BASE_URL}?{query}"


def _series_rows(
    payload: dict[str, object], site: str, time_series_id: str
) -> tuple[list[dict[str, str]], dict[str, object]]:
    candidates = payload.get("value", {}).get("timeSeries", [])
    matching: list[dict[str, object]] = []
    for candidate in candidates:
        source = candidate.get("sourceInfo", {})
        site_codes = source.get("siteCode", [])
        agency_site_pairs = {
            (str(code.get("agencyCode", "")), str(code.get("value", "")))
            for code in site_codes
        }
        variable = candidate.get("variable", {})
        pcodes = {
            str(code.get("value", ""))
            for code in variable.get("variableCode", [])
        }
        unit = str(variable.get("unit", {}).get("unitCode", "")).lower()
        options = variable.get("options", {}).get("option", [])
        statistics = {
            str(option.get("optionCode", ""))
            for option in options
            if str(option.get("name", "")).lower() == "statistic"
        }
        if (
            ("USGS", site) in agency_site_pairs
            and PARAMETER_CODE in pcodes
            and unit == EXPECTED_UNIT
            and EXPECTED_STATISTIC in statistics
        ):
            matching.append(candidate)
    if len(matching) != 1:
        raise ValueError(
            f"{site}: expected exactly one USGS 00300 mg/l instantaneous "
            f"series, found {len(matching)}"
        )

    candidate = matching[0]
    source = candidate["sourceInfo"]
    variable = candidate["variable"]
    value_groups = candidate.get("values", [])
    if not value_groups:
        raise ValueError(f"{site}: the matched time series has no value group")
    qualifier_definitions = {
        str(item.get("qualifierCode", "")): str(
            item.get("qualifierDescription", "")
        )
        for group in value_groups
        for item in group.get("qualifier", [])
    }
    all_method_ids = sorted(
        {
            str(item.get("methodID", ""))
            for group in value_groups
            for item in group.get("method", [])
        }
    )
    method_ids = ";".join(all_method_ids)
    rows: list[dict[str, str]] = []
    for group in value_groups:
        group_method_ids = ";".join(
            sorted(
                {str(item.get("methodID", "")) for item in group.get("method", [])}
            )
        )
        for item in group.get("value", []):
            qualifiers = [str(value) for value in item.get("qualifiers", [])]
            timestamp_source = str(item.get("dateTime", ""))
            parsed = datetime.fromisoformat(
                timestamp_source.replace("Z", "+00:00")
            )
            timestamp_utc = parsed.astimezone(timezone.utc).isoformat().replace(
                "+00:00", "Z"
            )
            rows.append(
                {
                    "monitoring_location_id": f"USGS-{site}",
                    "site_number": site,
                    "time_series_id": time_series_id,
                    "site_name": str(source.get("siteName", "")),
                    "agency_code": "USGS",
                    "parameter_code": PARAMETER_CODE,
                    "variable_name": str(variable.get("variableName", "")),
                    "unit": EXPECTED_UNIT,
                    "statistic_code": EXPECTED_STATISTIC,
                    "method_ids": group_method_ids,
                    "timestamp_source": timestamp_source,
                    "time_utc": timestamp_utc,
                    "value_mgL": str(item.get("value", "")),
                    "qualifiers": ";".join(qualifiers),
                    "qualifier_descriptions": ";".join(
                        qualifier_definitions.get(code, "")
                        for code in qualifiers
                    ),
                }
            )
    rows.sort(key=lambda row: row["time_utc"])
    metadata = {
        "monitoring_location_id": f"USGS-{site}",
        "site_number": site,
        "time_series_id": time_series_id,
        "site_name": str(source.get("siteName", "")),
        "latitude": source.get("geoLocation", {})
        .get("geogLocation", {})
        .get("latitude"),
        "longitude": source.get("geoLocation", {})
        .get("geogLocation", {})
        .get("longitude"),
        "huc_code": next(
            (
                item.get("value")
                for item in source.get("siteProperty", [])
                if item.get("name") == "hucCd"
            ),
            None,
        ),
        "series_name": str(candidate.get("name", "")),
        "variable_name": str(variable.get("variableName", "")),
        "unit": EXPECTED_UNIT,
        "statistic_code": EXPECTED_STATISTIC,
        "method_ids": method_ids,
        "value_group_count": len(value_groups),
        "qualifier_definitions": qualifier_definitions,
        "row_count": len(rows),
        "first_time_utc": rows[0]["time_utc"] if rows else None,
        "last_time_utc": rows[-1]["time_utc"] if rows else None,
    }
    return rows, metadata


def run() -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    NORMALIZED.mkdir(parents=True, exist_ok=True)
    accessed = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    fieldnames = [
        "monitoring_location_id",
        "site_number",
        "time_series_id",
        "site_name",
        "agency_code",
        "parameter_code",
        "variable_name",
        "unit",
        "statistic_code",
        "method_ids",
        "timestamp_source",
        "time_utc",
        "value_mgL",
        "qualifiers",
        "qualifier_descriptions",
    ]
    records: list[dict[str, object]] = []
    for index, (site, time_series_id) in enumerate(SERIES):
        url = _query_url(site)
        raw_path = RAW / f"USGS-{site}_00300_20240531_20240902.json"
        if raw_path.exists():
            raw_bytes = raw_path.read_bytes()
        else:
            raw_bytes = _download(url)
            raw_path.write_bytes(raw_bytes)
        payload = json.loads(raw_bytes.decode("utf-8"))
        rows, metadata = _series_rows(payload, site, time_series_id)
        csv_path = NORMALIZED / f"USGS-{site}_00300.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        records.append(
            {
                **metadata,
                "query_url": url,
                "raw_file": raw_path.relative_to(OUTPUT).as_posix(),
                "raw_sha256": _sha256(raw_path),
                "normalized_file": csv_path.relative_to(OUTPUT).as_posix(),
                "normalized_sha256": _sha256(csv_path),
            }
        )
        if index + 1 < len(SERIES):
            time.sleep(1.0)

    manifest = {
        "status": "downloaded_and_adapted",
        "accessed_utc": accessed,
        "source_authority": "U.S. Geological Survey",
        "source_database": "USGS Water Data for the Nation / NWIS",
        "metadata_selection_endpoint": OGC_METADATA_URL,
        "value_service": BASE_URL,
        "citation_url": CITATION_URL,
        "public_domain_guidance_url": PUBLIC_DOMAIN_URL,
        "recommended_dataset_citation": (
            "U.S. Geological Survey, 2026, U.S. Geological Survey National "
            "Water Information System database, accessed August 2, 2026, "
            "at https://doi.org/10.5066/F7P55KJN."
        ),
        "parameter_code": PARAMETER_CODE,
        "unit": EXPECTED_UNIT,
        "legacy_iv_statistic_code": EXPECTED_STATISTIC,
        "retrieval_envelope_local_dates": [START_DATE, END_DATE],
        "frozen_analysis_interval_utc_half_open": [
            ANALYSIS_START_UTC,
            ANALYSIS_END_UTC,
        ],
        "selection_note": (
            "All nine metadata-eligible series were retained before endpoint "
            "computation; value direction did not affect selection."
        ),
        "series": records,
        "download_script": "03_代码与测试/current/download_usgs_external_do.py",
        "download_script_sha256": _sha256(Path(__file__).resolve()),
    }
    manifest_path = OUTPUT / "download_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest_path


if __name__ == "__main__":
    print(run())
