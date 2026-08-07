from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


EXPERIMENTS = Path(__file__).resolve().parents[1] / "experiments"
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

from future_partial_identification_maracoos_pilot import (  # noqa: E402
    cross_group_collision_diagnostics,
)


class CrossGroupCollisionDiagnosticsTests(unittest.TestCase):
    def test_strict_half_cadence_boundary(self) -> None:
        present = pd.Series(pd.to_datetime(["2020-01-01T00:00:00Z"], utc=True))
        inside = pd.Series(pd.to_datetime(["2020-01-01T00:07:29Z"], utc=True))
        boundary = pd.Series(pd.to_datetime(["2020-01-01T00:07:30Z"], utc=True))

        self.assertEqual(cross_group_collision_diagnostics(present, inside), (1, 1, 1.0))
        self.assertEqual(
            cross_group_collision_diagnostics(present, boundary),
            (0, 1, 0.0),
        )

    def test_multiday_gap_is_not_misread_as_subminute_collision(self) -> None:
        present = pd.Series(pd.to_datetime(["2020-01-01T00:00:00Z"], utc=True))
        missing = pd.Series(pd.to_datetime(["2020-01-06T00:00:00Z"], utc=True))

        self.assertEqual(
            cross_group_collision_diagnostics(present, missing),
            (0, 1, 0.0),
        )


if __name__ == "__main__":
    unittest.main()

