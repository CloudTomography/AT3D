#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities to build extended AT3D grid-data CSV files.

This module is intentionally independent from the simulation script so that
MATLAB/1DRT inversion outputs can be converted into a richer CSV schema.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd


# Core columns already used by at3d.util.load_from_csv(..., density='cv').
BASE_COLUMNS = ["x", "y", "z", "cv", "reff", "veff", "lat", "lon"]

# Optional surface columns (2D fields repeated over z is allowed).
SURFACE_COLUMNS = [
    "surface_model",  # lambertian | wave_fresnel | diner | ocean_unpolarized | rpv_unpolarized
    "surface_albedo",
    "surface_wind_speed",
    "surface_pigmentation",
    "surface_refractive_index_real",
    "surface_refractive_index_imag",
    "surface_parameter_1",
    "surface_parameter_2",
    "surface_parameter_3",
]


@dataclass
class ExtendedGridOptions:
    """Controls which optional fields are materialized."""

    mode_count: int = 2
    include_surface: bool = True
    include_per_grid_refractive_index: bool = True


def _mode_columns(mode_count: int) -> List[str]:
    cols: List[str] = []
    for mode in range(1, mode_count + 1):
        cols.extend(
            [
                f"mode{mode}_fraction",
                f"mode{mode}_reff",
                f"mode{mode}_veff",
                f"mode{mode}_mr",
                f"mode{mode}_mi",
            ]
        )
    return cols


def _ensure_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = [name for name in required if name not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def validate_extended_grid_data(df: pd.DataFrame, options: ExtendedGridOptions) -> None:
    """Validate an extended grid-data table before writing.

    Notes
    -----
    - `mode*_fraction` should sum to 1 for rows with finite values.
    - `surface_model` is optional; when absent, caller should fallback to
      `at3d.surface.lambertian(...)` in the simulation stage.
    """

    required = list(BASE_COLUMNS)
    required.extend(_mode_columns(options.mode_count))
    _ensure_columns(df, required)

    if options.include_surface:
        # Surface columns are optional by design. Validate only if provided.
        provided_surface_cols = [c for c in SURFACE_COLUMNS if c in df.columns]
        if "surface_model" in provided_surface_cols:
            valid_models = {
                "lambertian",
                "wave_fresnel",
                "diner",
                "ocean_unpolarized",
                "rpv_unpolarized",
            }
            bad = sorted(set(df["surface_model"].dropna().astype(str)) - valid_models)
            if bad:
                raise ValueError(f"Unsupported surface_model values: {bad}")

    fractions = [f"mode{m}_fraction" for m in range(1, options.mode_count + 1)]
    frac_sum = df[fractions].sum(axis=1)
    finite = np.isfinite(frac_sum)
    if np.any(np.abs(frac_sum[finite] - 1.0) > 1e-3):
        raise ValueError("mode fractions must sum to 1 (±1e-3) for finite rows")


def write_extended_grid_csv(
    output_path: str | Path,
    df: pd.DataFrame,
    nx: int,
    ny: int,
    nz: int,
    dx_km: float,
    dy_km: float,
    z_levels_km: Sequence[float],
    options: Optional[ExtendedGridOptions] = None,
    extra_metadata: Optional[Mapping[str, str]] = None,
) -> Path:
    """Write an extended AT3D-compatible grid CSV.

    The first four header lines keep backward compatibility with
    `at3d.util.load_from_csv`.
    """

    options = options or ExtendedGridOptions()
    validate_extended_grid_data(df, options)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    z_line = ",".join(f"{float(v):.6f}" for v in z_levels_km)
    if len(z_levels_km) != nz:
        raise ValueError("len(z_levels_km) must equal nz")

    with out.open("w", encoding="utf-8") as f:
        f.write("# extended AT3D grid data\n")
        if extra_metadata:
            for k, v in extra_metadata.items():
                f.write(f"# {k}: {v}\n")
        f.write(f"{nx},{ny},{nz}\n")
        f.write(f"{dx_km:.6f},{dy_km:.6f}\n")
        f.write(f"{z_line}\n")
        df.to_csv(f, index=False)

    return out


def build_extended_grid_dataframe(base_columns: Mapping[str, np.ndarray], options: Optional[ExtendedGridOptions] = None) -> pd.DataFrame:
    """Build a normalized DataFrame from dict-like arrays.

    This helper is intended as the Python side counterpart of a MATLAB 1DRT
    inversion exporter.
    """

    options = options or ExtendedGridOptions()
    records = {k: np.asarray(v).reshape(-1) for k, v in base_columns.items()}
    df = pd.DataFrame(records)
    validate_extended_grid_data(df, options)
    return df
