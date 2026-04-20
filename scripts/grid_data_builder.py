#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities to build AT3D grid-data CSV files.

Design goals
------------
1) Keep AT3D CSV header compatibility.
2) Support both retrieval-1D style and LES-style sources.
3) Allow optional multi-mode/surface/per-grid-refractive-index columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xarray as xr


BASE_COLUMNS = ["x", "y", "z", "cv", "reff", "veff", "lat", "lon"]
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
    mode_count: int = 2
    include_surface: bool = True
    include_per_grid_refractive_index: bool = True


@dataclass
class GridGeometry:
    nx: int
    ny: int
    nz: int
    dx_km: float
    dy_km: float
    z_levels_km: Sequence[float]


def _mode_columns(mode_count: int) -> List[str]:
    cols: List[str] = []
    for mode in range(1, mode_count + 1):
        cols.extend([
            f"mode{mode}_fraction",
            f"mode{mode}_reff",
            f"mode{mode}_veff",
            f"mode{mode}_mr",
            f"mode{mode}_mi",
        ])
    return cols


def _ensure_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = [name for name in required if name not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _broadcast_2d_to_3d(arr2d: np.ndarray, nz: int) -> np.ndarray:
    return np.repeat(arr2d[:, :, None], nz, axis=2)


def _coarse_mean_2d(arr: np.ndarray, cy: int, cx: int) -> np.ndarray:
    ny, nx = arr.shape
    ny_c = ny // cy
    nx_c = nx // cx
    arr_t = arr[: ny_c * cy, : nx_c * cx]
    return arr_t.reshape(ny_c, cy, nx_c, cx).mean(axis=(1, 3))


def _coarse_mean_3d(arr: np.ndarray, cz: int, cy: int, cx: int) -> np.ndarray:
    nz, ny, nx = arr.shape
    nz_c = nz // cz
    ny_c = ny // cy
    nx_c = nx // cx
    arr_t = arr[: nz_c * cz, : ny_c * cy, : nx_c * cx]
    return arr_t.reshape(nz_c, cz, ny_c, cy, nx_c, cx).mean(axis=(1, 3, 5))


def _flatten_xyz(arr3d_zyx: np.ndarray) -> np.ndarray:
    # Input order: (z, y, x); output follows x/y/z index columns.
    return np.transpose(arr3d_zyx, (2, 1, 0)).reshape(-1)


def generate_latlon(lat0: float, lon0: float, dx_m: float, dy_m: float, ny: int, nx: int) -> Tuple[np.ndarray, np.ndarray]:
    """Generate regular-grid lat/lon from origin + spacing (LES fallback helper)."""
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(np.deg2rad(lat0))
    lat_1d = lat0 + (np.arange(ny) * dy_m / 1000.0) / km_per_deg_lat
    lon_1d = lon0 + (np.arange(nx) * dx_m / 1000.0) / km_per_deg_lon
    return np.meshgrid(lat_1d, lon_1d, indexing='ij')


def validate_extended_grid_data(df: pd.DataFrame, options: ExtendedGridOptions) -> None:
    required = list(BASE_COLUMNS)
    required.extend(_mode_columns(options.mode_count))
    _ensure_columns(df, required)

    if options.include_surface:
        if "surface_model" in df.columns:
            valid_models = {
                "lambertian", "wave_fresnel", "diner", "ocean_unpolarized", "rpv_unpolarized"
            }
            bad = sorted(set(df["surface_model"].dropna().astype(str)) - valid_models)
            if bad:
                raise ValueError(f"Unsupported surface_model values: {bad}")

    fractions = [f"mode{m}_fraction" for m in range(1, options.mode_count + 1)]
    frac_sum = df[fractions].sum(axis=1)
    finite = np.isfinite(frac_sum)
    if np.any(np.abs(frac_sum[finite] - 1.0) > 1e-3):
        raise ValueError("mode fractions must sum to 1 (±1e-3) for finite rows")


def build_extended_grid_dataframe(base_columns: Mapping[str, np.ndarray], options: Optional[ExtendedGridOptions] = None) -> pd.DataFrame:
    options = options or ExtendedGridOptions()
    records = {k: np.asarray(v).reshape(-1) for k, v in base_columns.items()}
    df = pd.DataFrame(records)
    validate_extended_grid_data(df, options)
    return df


def write_extended_grid_csv(
    output_path: str | Path,
    df: pd.DataFrame,
    geometry: GridGeometry,
    options: Optional[ExtendedGridOptions] = None,
    extra_metadata: Optional[Mapping[str, str]] = None,
) -> Path:
    options = options or ExtendedGridOptions()
    validate_extended_grid_data(df, options)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if len(geometry.z_levels_km) != geometry.nz:
        raise ValueError("len(z_levels_km) must equal nz")
    z_line = ",".join(f"{float(v):.6f}" for v in geometry.z_levels_km)

    with out.open("w", encoding="utf-8") as f:
        f.write("# extended AT3D grid data\n")
        if extra_metadata:
            for k, v in extra_metadata.items():
                f.write(f"# {k}: {v}\n")
        f.write(f"{geometry.nx},{geometry.ny},{geometry.nz}\n")
        f.write(f"{geometry.dx_km:.6f},{geometry.dy_km:.6f}\n")
        f.write(f"{z_line}\n")
        df.to_csv(f, index=False)

    return out


def build_from_retrieval_1d_netcdf(
    nc_path: str | Path,
    z_levels_km: Sequence[float],
    dx_km: float,
    dy_km: float,
    mode_count: int = 2,
    lat_var: str = "lat",
    lon_var: str = "lon",
    cv_var: str = "Cv_total",
    fine_fraction_var: str = "FineModeFraction",
    reff_vars: Optional[Sequence[str]] = None,
    veff_vars: Optional[Sequence[str]] = None,
    mr_vars: Optional[Sequence[str]] = None,
    mi_vars: Optional[Sequence[str]] = None,
    default_veff: float = 0.10,
    fallback_lat0: Optional[float] = None,
    fallback_lon0: Optional[float] = None,
) -> Tuple[pd.DataFrame, GridGeometry, ExtendedGridOptions]:
    """Convert a retrieval-1D netCDF (2D fields) into extended 3D grid dataframe.

    By default this uses common variable names from your 1DRT export. Missing fields
    are filled with conservative defaults so the pipeline can still run.
    """
    ds = xr.open_dataset(nc_path)

    cv2d = np.asarray(ds[cv_var].values, dtype=float)
    ny, nx = cv2d.shape
    nz = len(z_levels_km)

    if reff_vars is None:
        reff_vars = ["R_g_1", "R_g_2"][:mode_count]
    if veff_vars is None:
        veff_vars = ["Sigma_g_1", "Sigma_g_2"][:mode_count]
    if mr_vars is None:
        mr_vars = ["Refractive_Index_Real_Mode_1", "Refractive_Index_Real_Mode_2"][:mode_count]
    if mi_vars is None:
        mi_vars = ["Refractive_Index_Imaginary_Mode_1", "Refractive_Index_Imaginary_Mode_2"][:mode_count]

    # lat/lon: prefer file, else regular-grid fallback
    if lat_var in ds and lon_var in ds:
        lat2d = np.asarray(ds[lat_var].values, dtype=float)
        lon2d = np.asarray(ds[lon_var].values, dtype=float)
    elif fallback_lat0 is not None and fallback_lon0 is not None:
        lat2d, lon2d = generate_latlon(fallback_lat0, fallback_lon0, dx_km * 1000.0, dy_km * 1000.0, ny, nx)
    else:
        lat2d = np.zeros((ny, nx), dtype=float)
        lon2d = np.zeros((ny, nx), dtype=float)

    cv3d = _broadcast_2d_to_3d(cv2d, nz)
    lat3d = _broadcast_2d_to_3d(lat2d, nz)
    lon3d = _broadcast_2d_to_3d(lon2d, nz)

    mode_fraction_2d: List[np.ndarray] = []
    if mode_count == 1:
        mode_fraction_2d = [np.ones_like(cv2d)]
    else:
        fine = np.asarray(ds[fine_fraction_var].values, dtype=float) if fine_fraction_var in ds else np.full_like(cv2d, 0.5)
        mode_fraction_2d.append(fine)
        if mode_count >= 2:
            mode_fraction_2d.append(1.0 - fine)
        for _ in range(3, mode_count + 1):
            mode_fraction_2d.append(np.zeros_like(cv2d))

    base: Dict[str, np.ndarray] = {}
    X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="xy")
    base["x"] = X.reshape(-1)
    base["y"] = Y.reshape(-1)
    base["z"] = Z.reshape(-1)
    base["cv"] = _flatten_xyz(cv3d)
    # Backward-compatible scalar reff/veff columns use mode1.
    mode1_reff2d = np.asarray(ds[reff_vars[0]].values, dtype=float) if reff_vars and reff_vars[0] in ds else np.full_like(cv2d, 0.2)
    mode1_veff2d = np.asarray(ds[veff_vars[0]].values, dtype=float) if veff_vars and veff_vars[0] in ds else np.full_like(cv2d, default_veff)
    base["reff"] = _flatten_xyz(_broadcast_2d_to_3d(mode1_reff2d, nz))
    base["veff"] = _flatten_xyz(_broadcast_2d_to_3d(mode1_veff2d, nz))
    base["lat"] = _flatten_xyz(lat3d)
    base["lon"] = _flatten_xyz(lon3d)

    for i_mode in range(mode_count):
        m = i_mode + 1
        frac2d = mode_fraction_2d[i_mode]
        reff2d = np.asarray(ds[reff_vars[i_mode]].values, dtype=float) if i_mode < len(reff_vars) and reff_vars[i_mode] in ds else mode1_reff2d
        veff2d = np.asarray(ds[veff_vars[i_mode]].values, dtype=float) if i_mode < len(veff_vars) and veff_vars[i_mode] in ds else np.full_like(cv2d, default_veff)
        mr2d = np.asarray(ds[mr_vars[i_mode]].values, dtype=float) if i_mode < len(mr_vars) and mr_vars[i_mode] in ds else np.full_like(cv2d, np.nan)
        mi2d = np.asarray(ds[mi_vars[i_mode]].values, dtype=float) if i_mode < len(mi_vars) and mi_vars[i_mode] in ds else np.full_like(cv2d, np.nan)

        base[f"mode{m}_fraction"] = _flatten_xyz(_broadcast_2d_to_3d(frac2d, nz))
        base[f"mode{m}_reff"] = _flatten_xyz(_broadcast_2d_to_3d(reff2d, nz))
        base[f"mode{m}_veff"] = _flatten_xyz(_broadcast_2d_to_3d(veff2d, nz))
        base[f"mode{m}_mr"] = _flatten_xyz(_broadcast_2d_to_3d(mr2d, nz))
        base[f"mode{m}_mi"] = _flatten_xyz(_broadcast_2d_to_3d(mi2d, nz))

    options = ExtendedGridOptions(mode_count=mode_count, include_surface=True, include_per_grid_refractive_index=True)
    df = build_extended_grid_dataframe(base, options)
    geom = GridGeometry(nx=nx, ny=ny, nz=nz, dx_km=dx_km, dy_km=dy_km, z_levels_km=list(z_levels_km))
    return df, geom, options


def build_from_les_arrays(
    qvapor_zyx: np.ndarray,
    z_levels_m: np.ndarray,
    lat_2d: np.ndarray,
    lon_2d: np.ndarray,
    dx_m: float,
    dy_m: float,
    coarse: Tuple[int, int, int] = (2, 4, 4),
    rho_air_kgm3: float = 1.0,
    default_reff: float = 12.0,
    default_veff: float = 0.10,
) -> Tuple[pd.DataFrame, GridGeometry, ExtendedGridOptions]:
    """Build extended dataframe from LES-like arrays (compatible with your WRF workflow)."""
    cz, cy, cx = coarse
    qv_c = _coarse_mean_3d(qvapor_zyx, cz, cy, cx)
    lat_c = _coarse_mean_2d(lat_2d, cy, cx)
    lon_c = _coarse_mean_2d(lon_2d, cy, cx)
    z_c = np.asarray(z_levels_m[: (len(z_levels_m) // cz) * cz], dtype=float).reshape(-1, cz).mean(axis=1)

    lwc_gm3 = qv_c * rho_air_kgm3 * 1000.0
    nz, ny, nx = lwc_gm3.shape

    X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="xy")
    base: Dict[str, np.ndarray] = {
        "x": X.reshape(-1),
        "y": Y.reshape(-1),
        "z": Z.reshape(-1),
        "cv": _flatten_xyz(lwc_gm3),
        "reff": _flatten_xyz(np.full_like(lwc_gm3, default_reff, dtype=float)),
        "veff": _flatten_xyz(np.full_like(lwc_gm3, default_veff, dtype=float)),
        "lat": _flatten_xyz(_broadcast_2d_to_3d(lat_c, nz)),
        "lon": _flatten_xyz(_broadcast_2d_to_3d(lon_c, nz)),
        "mode1_fraction": np.ones(nx * ny * nz),
        "mode1_reff": _flatten_xyz(np.full_like(lwc_gm3, default_reff, dtype=float)),
        "mode1_veff": _flatten_xyz(np.full_like(lwc_gm3, default_veff, dtype=float)),
        "mode1_mr": np.full(nx * ny * nz, np.nan),
        "mode1_mi": np.full(nx * ny * nz, np.nan),
    }
    options = ExtendedGridOptions(mode_count=1, include_surface=True, include_per_grid_refractive_index=False)
    df = build_extended_grid_dataframe(base, options)
    geom = GridGeometry(
        nx=nx,
        ny=ny,
        nz=nz,
        dx_km=(dx_m * cx) / 1000.0,
        dy_km=(dy_m * cy) / 1000.0,
        z_levels_km=(z_c / 1000.0).tolist(),
    )
    return df, geom, options



def parse_z_levels(z: str) -> List[float]:
    """Accept "0.01:0.5:20" or "0.01,0.51,1.01"."""
    if ":" in z:
        a, b, c = [float(v) for v in z.split(":")]
        return list(np.arange(a, c + 1e-12, b))
    return [float(v) for v in z.split(",") if v.strip()]


def run_retrieval_case(
    input_nc: str,
    output_csv: str,
    dx_km: float,
    dy_km: float,
    z_levels_km: Sequence[float],
    mode_count: int = 2,
) -> Path:
    """Spyder-friendly wrapper: one function call to build CSV from retrieval nc."""
    df, geom, options = build_from_retrieval_1d_netcdf(
        nc_path=input_nc,
        z_levels_km=z_levels_km,
        dx_km=dx_km,
        dy_km=dy_km,
        mode_count=mode_count,
    )
    return write_extended_grid_csv(output_csv, df, geom, options)


if __name__ == "__main__":
    # =========================
    # Spyder editable run block
    # =========================
    # 直接修改下面变量后，在 Spyder 中 Run File 即可。
    input_nc = "../data/retrieval_1d/2019_0806_1839_N_Pxl25_3_3.nc"
    output_csv = "../data/synthetic_cloud_fields/jpl_les/retrieval_2019_0806_1839_extended.csv"
    dx_km = 0.16
    dy_km = 0.16
    z_levels_km = parse_z_levels("0.01:0.5:20")
    mode_count = 2

    out = run_retrieval_case(
        input_nc=input_nc,
        output_csv=output_csv,
        dx_km=dx_km,
        dy_km=dy_km,
        z_levels_km=z_levels_km,
        mode_count=mode_count,
    )
    print(f"✅ wrote: {out}")
