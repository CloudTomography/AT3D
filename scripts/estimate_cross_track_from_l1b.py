#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimate cross-track configuration fields for config_v5b.yaml from AirMSPI L1B + MetNav.
"""
from __future__ import annotations

import argparse
import re
import glob
import sys
import io
from pathlib import Path
from typing import Dict, List, Tuple, Any

import h5py
import numpy as np
import pandas as pd
import xarray as xr


TIME_RE = re.compile(r"_(\d{8})_(\d{6})Z")


def _parse_base_time_from_name(path: Path) -> pd.Timestamp:
    m = TIME_RE.search(path.name)
    if not m:
        raise ValueError(f"Cannot parse base time from filename: {path.name}")
    return pd.to_datetime(f"{m.group(1)} {m.group(2)}", format="%Y%m%d %H%M%S", utc=True)


def _read_h5_array(f: h5py.File, dataset_path: str) -> np.ndarray:
    if dataset_path not in f:
        raise KeyError(f"Dataset not found: {dataset_path}")
    arr = np.asarray(f[dataset_path][...], dtype=float)
    fill = f[dataset_path].attrs.get("_FillValue", None)
    if fill is not None:
        if isinstance(fill, np.ndarray):
            fill = float(np.ravel(fill)[0])
        else:
            fill = float(fill)
        arr = np.where(arr == fill, np.nan, arr)
    return arr


def _find_first_existing(f: h5py.File, candidates: List[str]) -> str:
    for c in candidates:
        if c in f:
            return c
    raise KeyError(f"None of datasets exist: {candidates}")


def _scan_time_window(l1b_path: Path, time_dataset: str) -> Tuple[pd.Timestamp, pd.Timestamp, Tuple[float, float], Tuple[float, float]]:
    def _infer_geo_path(base_path: str, field_name: str) -> str:
        return base_path.replace("Time_in_seconds_from_epoch", field_name)

    with h5py.File(l1b_path, "r") as f:
        t = _read_h5_array(f, time_dataset)
        lat = _read_h5_array(f, _infer_geo_path(time_dataset, "Latitude"))
        lon = _read_h5_array(f, _infer_geo_path(time_dataset, "Longitude"))
    valid = np.isfinite(t)
    if not np.any(valid):
        raise ValueError(f"No valid time values in {l1b_path}")
    tmin = float(np.nanmin(t[valid]))
    tmax = float(np.nanmax(t[valid]))
    base_time = _parse_base_time_from_name(l1b_path)
    # Some AirMSPI products store seconds relative to scene base time in filename,
    # while others may already be epoch seconds. Heuristic:
    # - small values (e.g., 0~86400) => relative seconds
    # - very large values (>=1e8) => epoch seconds
    if np.nanmax(np.abs([tmin, tmax])) < 1e8:
        t_start_utc, t_end_utc = (
            base_time + pd.to_timedelta(tmin, unit="s"),
            base_time + pd.to_timedelta(tmax, unit="s"),
        )
    else:
        t_start_utc, t_end_utc = pd.to_datetime(tmin, unit="s", utc=True), pd.to_datetime(tmax, unit="s", utc=True)

    geo_valid = (
        np.isfinite(lat) & np.isfinite(lon) & np.isfinite(t)
        & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
        & (lat != 0) & (lon != 0)
    )
    if not np.any(geo_valid):
        raise ValueError(f"No valid geo pixels for boundary midpoint in {l1b_path}")
    t_geo = np.where(geo_valid, t, np.nan)

    # Boundary near scan start/end: use lowest/highest 1% time among valid pixels.
    q_lo = np.nanpercentile(t_geo, 1.0)
    q_hi = np.nanpercentile(t_geo, 99.0)
    start_mask = geo_valid & (t_geo <= q_lo)
    end_mask = geo_valid & (t_geo >= q_hi)
    if not np.any(start_mask):
        start_mask = geo_valid & np.isclose(t_geo, tmin, atol=1e-6)
    if not np.any(end_mask):
        end_mask = geo_valid & np.isclose(t_geo, tmax, atol=1e-6)

    start_latlon = (float(np.nanmedian(lat[start_mask])), float(np.nanmedian(lon[start_mask])))
    end_latlon = (float(np.nanmedian(lat[end_mask])), float(np.nanmedian(lon[end_mask])))
    return t_start_utc, t_end_utc, start_latlon, end_latlon


def _latlon_to_xy_km(lat: np.ndarray, lon: np.ndarray, lat0: float, lon0: float) -> Tuple[np.ndarray, np.ndarray]:
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(np.deg2rad(lat0))
    x_north = (lat - lat0) * km_per_deg_lat
    y_east = (lon - lon0) * km_per_deg_lon
    return x_north, y_east


def _load_retrieval_sw_corner_latlon(retrieval_nc: Path) -> Tuple[float, float]:
    ds = xr.open_dataset(retrieval_nc)
    lat_name = "latitude" if "latitude" in ds else ("lat" if "lat" in ds else None)
    lon_name = "longitude" if "longitude" in ds else ("lon" if "lon" in ds else None)
    if lat_name is None or lon_name is None:
        raise ValueError("retrieval nc must contain latitude/longitude (or lat/lon)")
    lat = np.asarray(ds[lat_name].values, dtype=float)
    lon = np.asarray(ds[lon_name].values, dtype=float)
    # TODO: remove `(lat != 0) & (lon != 0)` once upstream retrieval nc data no longer
    # contains zero-valued invalid geolocation placeholders.
    mask = (
        np.isfinite(lat)
        & np.isfinite(lon)
        & (np.abs(lat) <= 90)
        & (np.abs(lon) <= 180)
        & (lat != 0)
        & (lon != 0)
    )
    if not np.any(mask):
        raise ValueError("No valid lat/lon in retrieval nc")
    return float(np.nanmin(lat[mask])), float(np.nanmin(lon[mask]))


def _load_metnav_table(path: Path, reference_date_utc: pd.Timestamp | None = None) -> pd.DataFrame:
    def _read_ict(p: Path) -> pd.DataFrame:
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        if not lines:
            raise ValueError(f"Empty ICT file: {p}")
        # Prefer explicit header row detection (works better for FIREX-AQ ICT variants).
        header_idx = None
        for i, ln in enumerate(lines):
            ls = ln.strip()
            if not ls:
                continue
            if ls.lower().startswith("time_start,") or ls.lower().startswith("time,"):
                header_idx = i
        if header_idx is not None:
            payload = "".join(lines[header_idx:])
            df_try = pd.read_csv(io.StringIO(payload))
        else:
            first = lines[0].strip().split(",")[0].strip()
            n_header = int(float(first))
            # ICARTT convention: variable names are usually on header line (n_header-1).
            df_try = pd.read_csv(p, skiprows=max(n_header - 1, 0))
            if df_try.shape[1] <= 1:
                df_try = pd.read_csv(p, skiprows=max(n_header - 1, 0), delim_whitespace=True)
        df_try.columns = [str(c).strip() for c in df_try.columns]
        df_try = df_try.replace([-9999, -8888, -7777], np.nan)
        return df_try

    if path.suffix.lower() in {".nc", ".nc4"}:
        ds = xr.open_dataset(path)
        df = ds.to_dataframe().reset_index()
    elif path.suffix.lower() == ".ict":
        df = _read_ict(path)
    else:
        try:
            df = pd.read_csv(path)
        except pd.errors.ParserError:
            # Some MetNav text products are ICARTT-like.
            df = _read_ict(path)
    cols = {c.lower(): c for c in df.columns}
    time_col = (
        cols.get("time_in_seconds_from_epoch")
        or cols.get("time")
        or cols.get("timestamp")
        or cols.get("time_start")
        or cols.get("utc")
    )
    if time_col is None:
        raise ValueError("MetNav file needs time column (Time_in_seconds_from_epoch / time / timestamp)")
    if np.issubdtype(df[time_col].dtype, np.number):
        t_num = df[time_col].astype(float)
        # MetNav daily products often store Time_Start as seconds from 00:00:00 (UTC)
        # instead of epoch seconds.
        if str(time_col).lower() in {"time_start", "utc"}:
            if reference_date_utc is None:
                raise ValueError(
                    "Numeric Time_Start/UTC detected in MetNav, but reference_date_utc is None."
                )
            day0 = pd.Timestamp(reference_date_utc).tz_convert("UTC").normalize()
            t = day0 + pd.to_timedelta(t_num, unit="s")
        else:
            t = pd.to_datetime(t_num, unit="s", utc=True)
    else:
        t = pd.to_datetime(df[time_col], utc=True)
    df = df.copy()
    df["time_utc"] = t
    return df.sort_values("time_utc")


def _interp_metnav(df: pd.DataFrame, t_utc: pd.Timestamp, col_candidates: List[str]) -> float:
    cols = {c.lower(): c for c in df.columns}
    src = None
    for c in col_candidates:
        if c.lower() in cols:
            src = cols[c.lower()]
            break
    if src is None:
        raise ValueError(f"Missing MetNav column, tried: {col_candidates}")
    x = df["time_utc"].astype("int64").values.astype(float) / 1e9
    y = pd.to_numeric(df[src], errors="coerce").values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    return float(np.interp(t_utc.value / 1e9, x[m], y[m]))


def estimate_cross_track(l1b_files: List[Path], metnav_path: Path, retrieval_nc: Path, time_dataset: str) -> List[Dict[str, Any]]:
    first_base_time = _parse_base_time_from_name(sorted(l1b_files)[0])
    met = _load_metnav_table(metnav_path, reference_date_utc=first_base_time)
    lat0, lon0 = _load_retrieval_sw_corner_latlon(retrieval_nc)
    per_view: List[Dict[str, Any]] = []
    for idx, p in enumerate(sorted(l1b_files), start=1):
        t_start, t_end, start_ground_latlon, end_ground_latlon = _scan_time_window(p, time_dataset)

        lat1 = _interp_metnav(met, t_start, ["Latitude"])
        lon1 = _interp_metnav(met, t_start, ["Longitude"])
        alt1_m = _interp_metnav(met, t_start, ["Altitude", "MSL_GPS_Altitude", "HAE_GPS_Altitude"])

        lat2 = _interp_metnav(met, t_end, ["Latitude"])
        lon2 = _interp_metnav(met, t_end, ["Longitude"])
        alt2_m = _interp_metnav(met, t_end, ["Altitude", "MSL_GPS_Altitude", "HAE_GPS_Altitude"])

        x1, y1 = _latlon_to_xy_km(np.array([lat1]), np.array([lon1]), lat0, lon0)
        x2, y2 = _latlon_to_xy_km(np.array([lat2]), np.array([lon2]), lat0, lon0)

        # Pitch from L1B scan boundary ground midpoint and aircraft position (not from MetNav pitch variable).
        sg_lat, sg_lon = start_ground_latlon
        eg_lat, eg_lon = end_ground_latlon
        dnx1, dey1 = _latlon_to_xy_km(np.array([lat1]), np.array([lon1]), sg_lat, sg_lon)
        dnx2, dey2 = _latlon_to_xy_km(np.array([lat2]), np.array([lon2]), eg_lat, eg_lon)
        hdist1_km = float(np.hypot(dnx1[0], dey1[0]))
        hdist2_km = float(np.hypot(dnx2[0], dey2[0]))
        alt1_km = float(alt1_m / 1000.0)
        alt2_km = float(alt2_m / 1000.0)
        pitch1 = float(np.degrees(np.arctan2(max(alt1_km, 0.0), max(hdist1_km, 1e-6))))
        pitch2 = float(np.degrees(np.arctan2(max(alt2_km, 0.0), max(hdist2_km, 1e-6))))

        per_view.append({
            "view_index": idx,
            "l1b_file": str(p),
            "cross_track_x1": float(x1[0]),
            "cross_track_y1": float(y1[0]),
            "cross_track_z1": float(alt1_m / 1000.0),
            "cross_track_x2": float(x2[0]),
            "cross_track_y2": float(y2[0]),
            "cross_track_z2": float(alt2_m / 1000.0),
            "cross_track_pitch_start_deg": float(pitch1),
            "cross_track_pitch_end_deg": float(pitch2),
            "scan_start_utc": t_start.isoformat(),
            "scan_end_utc": t_end.isoformat(),
            "reference_sw_lat": lat0,
            "reference_sw_lon": lon0,
        })
    return per_view


def main(argv: List[str] | None = None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--l1b_glob", default=None, help="e.g. 'data/l1b/Case_21_*_V006.hdf'")
    ap.add_argument("--metnav", default=None, help="MetNav file (.csv/.nc)")
    ap.add_argument("--retrieval_nc", default=None, help="retrieval nc for SW corner reference")
    ap.add_argument("--time_dataset", default="/HDFEOS/GRIDS/935nm_band/Data_Fields/Time_in_seconds_from_epoch")
    args = ap.parse_args(argv)

    # Spyder / notebook friendly defaults (avoid SystemExit when no CLI args are passed).
    if not args.l1b_glob and not args.metnav and not args.retrieval_nc:
        print("⚠️ No CLI arguments provided. Please pass --l1b_glob/--metnav/--retrieval_nc, or edit defaults in this file.")
        print("Example:")
        print("  python scripts/estimate_cross_track_from_l1b.py \\")
        print("    --l1b_glob 'data/l1b/Case_21_*.hdf' \\")
        print("    --metnav 'data/metnav/MetNav_AircraftInSitu_ER2_Data_1.csv' \\")
        print("    --retrieval_nc 'data/retrieval_1d/2019_0806_1839_N_Pxl25_3_3.nc'")
        return
    missing = [name for name in ["l1b_glob", "metnav", "retrieval_nc"] if getattr(args, name) in (None, "")]
    if missing:
        raise ValueError(f"Missing required arguments: {missing}")

    l1b_files = [Path(p) for p in sorted(glob.glob(args.l1b_glob))]
    if not l1b_files:
        raise FileNotFoundError(f"No files matched --l1b_glob: {args.l1b_glob}")

    est_list = estimate_cross_track(
        l1b_files=l1b_files,
        metnav_path=Path(args.metnav),
        retrieval_nc=Path(args.retrieval_nc),
        time_dataset=args.time_dataset,
    )
    for est in est_list:
        print(f"\n# view_index: {est['view_index']}")
        print(f"# l1b_file: {est['l1b_file']}")
        print("trajectory:")
        for k in [
            "cross_track_x1", "cross_track_y1", "cross_track_z1",
            "cross_track_x2", "cross_track_y2", "cross_track_z2",
            "cross_track_pitch_start_deg", "cross_track_pitch_end_deg",
        ]:
            print(f"  {k}: {est[k]:.6f}")
        print(f"# scan_start_utc: {est['scan_start_utc']}")
        print(f"# scan_end_utc:   {est['scan_end_utc']}")
        print(f"# reference_sw_latlon: ({est['reference_sw_lat']:.6f}, {est['reference_sw_lon']:.6f})")


if __name__ == "__main__":
    # =========================
    # Spyder editable run block
    # =========================
    # 在 Spyder 里直接 Run File 时可修改下面路径，无需手输 CLI。
    l1b_glob = "K:/AirMSPI_FIREX-AQ_Terrain-projected_Georegistered_Radiance_Data_6-20241212_232427/Case_21_AirMSPI_ER2_GRP_TERRAIN_20190806_*.hdf"
    metnav = "../data/FIREXAQ_MetNav_AircraftInSitu_ER2_Data_1-20260427_163432/FIREXAQ-METNAV_ER2_20190806_R0.ict"
    retrieval_nc = "../results/ncFile/test143/Case_21_2019_0806_1839/2019_0806_1839_N_Pxl25_3_3.nc"
    time_dataset = "/HDFEOS/GRIDS/935nm_band/Data Fields/Time_in_seconds_from_epoch"

    if len(sys.argv) > 1:
        main()
    else:
        if not l1b_glob or not metnav or not retrieval_nc:
            print("⚠️ 请在脚本底部 Spyder run block 中填写 l1b_glob / metnav / retrieval_nc。")
        else:
            main([
                "--l1b_glob", l1b_glob,
                "--metnav", metnav,
                "--retrieval_nc", retrieval_nc,
                "--time_dataset", time_dataset,
            ])
