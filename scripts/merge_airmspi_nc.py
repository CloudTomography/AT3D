#!/usr/bin/env python3
"""Merge multi-view, multi-band AirMSPI results into one NetCDF.

Output variables and dimensions are aligned with retrieval-style products:
- Full resolution (dim_x, dim_y): datalon, datalat, dataElevation, dataLand_water_mask
- Downsampled (dim_x_downsampling, dim_y_downsampling): lon, lat, elevation, Land_water_mask
- Multi-view/multi-band (dim_x_downsampling, dim_y_downsampling, dim_view, dim_band):
  I, Q, U, DoLP, ErrI, ErrQ, ErrU, ErrDoLP, theta0, thetav, faipfai0
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr


def _crop_to_factor(a: np.ndarray, factor: int) -> np.ndarray:
    ny, nx = a.shape
    ny2 = (ny // factor) * factor
    nx2 = (nx // factor) * factor
    if ny2 == 0 or nx2 == 0:
        raise ValueError(f"shape {a.shape} is too small for factor={factor}")
    return a[:ny2, :nx2]


def _downsample_mean2d(a: np.ndarray, factor: int) -> np.ndarray:
    c = _crop_to_factor(a, factor)
    ny, nx = c.shape
    return np.nanmean(c.reshape(ny // factor, factor, nx // factor, factor), axis=(1, 3))


def _choose_first(ds: xr.Dataset, names: Iterable[str]) -> str:
    for n in names:
        if n in ds:
            return n
    raise KeyError(f"None of the candidate variables exist: {list(names)}")


def _to_numpy2d(ds: xr.Dataset, varname: str) -> np.ndarray:
    arr = ds[varname].values
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"{varname} must be 2D after squeeze, got {arr.shape}")
    return np.asarray(arr, dtype=np.float64)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True, help="Input NetCDF files (per-view/per-band or mixed)")
    p.add_argument("--output", required=True, help="Output merged NetCDF file")
    p.add_argument("--factor", type=int, default=25, help="Downsampling factor (default: 25)")
    args = p.parse_args()

    input_files = [Path(i) for i in args.inputs]
    if not input_files:
        raise ValueError("No input files provided")

    opened = [xr.open_dataset(f) for f in input_files]

    lat_name = _choose_first(opened[0], ["datalat", "latitude", "lat"])
    lon_name = _choose_first(opened[0], ["datalon", "longitude", "lon"])
    elev_name = _choose_first(opened[0], ["dataElevation", "elevation"])
    lwm_name = _choose_first(opened[0], ["dataLand_water_mask", "Land_water_mask", "land_water_mask"])

    datalat = _to_numpy2d(opened[0], lat_name)
    datalon = _to_numpy2d(opened[0], lon_name)
    data_elev = _to_numpy2d(opened[0], elev_name)
    data_lwm = _to_numpy2d(opened[0], lwm_name)

    lat_ds = _downsample_mean2d(datalat, args.factor)
    lon_ds = _downsample_mean2d(datalon, args.factor)
    elev_ds = _downsample_mean2d(data_elev, args.factor)
    lwm_ds = _downsample_mean2d(data_lwm, args.factor)

    stack = []
    bands_ref = None
    for ds in opened:
        i_name = _choose_first(ds, ["I"])
        q_name = _choose_first(ds, ["Q"])
        u_name = _choose_first(ds, ["U"])
        dolp_name = _choose_first(ds, ["DoLP"])
        erri_name = _choose_first(ds, ["ErrI"])
        errq_name = _choose_first(ds, ["ErrQ"])
        erru_name = _choose_first(ds, ["ErrU"])
        errdolp_name = _choose_first(ds, ["ErrDoLP"])
        theta0_name = _choose_first(ds, ["theta0"])
        thetav_name = _choose_first(ds, ["thetav"])
        faipfai0_name = _choose_first(ds, ["faipfai0"])

        band_name = "band" if "band" in ds.dims else ("dim_band" if "dim_band" in ds.dims else None)
        if band_name is None:
            # fallback: infer single band
            bands = np.array([np.nan], dtype=np.float64)
        else:
            bands = np.asarray(ds[band_name].values, dtype=np.float64)

        if bands_ref is None:
            bands_ref = bands
        elif len(bands_ref) != len(bands):
            raise ValueError("All inputs must have the same band size")

        stack.append(
            dict(
                I=np.asarray(ds[i_name].values, dtype=np.float64),
                Q=np.asarray(ds[q_name].values, dtype=np.float64),
                U=np.asarray(ds[u_name].values, dtype=np.float64),
                DoLP=np.asarray(ds[dolp_name].values, dtype=np.float64),
                ErrI=np.asarray(ds[erri_name].values, dtype=np.float64),
                ErrQ=np.asarray(ds[errq_name].values, dtype=np.float64),
                ErrU=np.asarray(ds[erru_name].values, dtype=np.float64),
                ErrDoLP=np.asarray(ds[errdolp_name].values, dtype=np.float64),
                theta0=np.asarray(ds[theta0_name].values, dtype=np.float64),
                thetav=np.asarray(ds[thetav_name].values, dtype=np.float64),
                faipfai0=np.asarray(ds[faipfai0_name].values, dtype=np.float64),
            )
        )

    # Expected shape per file: (dim_x_downsampling, dim_y_downsampling, dim_band)
    keys = list(stack[0].keys())
    merged = {k: np.stack([s[k] for s in stack], axis=2) for k in keys}  # -> x, y, view, band

    nx_ds, ny_ds = lat_ds.shape
    nview = len(stack)
    nband = merged["I"].shape[-1]

    out = xr.Dataset(
        data_vars={
            "DoLP": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["DoLP"]),
            "datalon": (("dim_x", "dim_y"), datalon),
            "lon": (("dim_x_downsampling", "dim_y_downsampling"), lon_ds),
            "Height_AirMSPI": (("height_dim", "height_dim2"), np.array([[20.0]], dtype=np.float64)),
            "thetav": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["thetav"]),
            "Land_water_mask": (("dim_x_downsampling", "dim_y_downsampling"), lwm_ds),
            "elevation": (("dim_x_downsampling", "dim_y_downsampling"), elev_ds),
            "datalat": (("dim_x", "dim_y"), datalat),
            "Band_AirMSPI": (("dim_band", "band_scalar"), np.asarray(bands_ref).reshape(-1, 1)),
            "Q": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["Q"]),
            "dataLand_water_mask": (("dim_x", "dim_y"), data_lwm),
            "lat": (("dim_x_downsampling", "dim_y_downsampling"), lat_ds),
            "ErrI": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["ErrI"]),
            "U": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["U"]),
            "ErrU": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["ErrU"]),
            "dataElevation": (("dim_x", "dim_y"), data_elev),
            "faipfai0": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["faipfai0"]),
            "theta0": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["theta0"]),
            "ErrQ": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["ErrQ"]),
            "I": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["I"]),
            "ErrDoLP": (("dim_x_downsampling", "dim_y_downsampling", "dim_view", "dim_band"), merged["ErrDoLP"]),
        },
        coords={
            "dim_x": np.arange(datalat.shape[0]),
            "dim_y": np.arange(datalat.shape[1]),
            "dim_x_downsampling": np.arange(nx_ds),
            "dim_y_downsampling": np.arange(ny_ds),
            "dim_view": np.arange(nview),
            "dim_band": np.arange(nband),
            "height_dim": [0],
            "height_dim2": [0],
            "band_scalar": [0],
        },
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(args.output)
    for ds in opened:
        ds.close()
    print(f"Saved merged file: {args.output}")
    print(f"dims: dim_x={datalat.shape[0]}, dim_y={datalat.shape[1]}, "
          f"dim_x_downsampling={nx_ds}, dim_y_downsampling={ny_ds}, dim_view={nview}, dim_band={nband}")


if __name__ == "__main__":
    main()
