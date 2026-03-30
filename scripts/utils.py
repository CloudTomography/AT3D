#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:30:03 2025

@author: benting
"""
import os
import sys
import time
import yaml
import argparse
import numpy as np
import xarray as xr
from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

# External libs
import at3d
from scipy.ndimage import map_coordinates
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
from at3dclass import (
    OutputConfig, SensorConfig, BandsConfig, DownsampleConfig, SaveVersionsConfig,
    PlotConfig, GroundGridConfig, ComputeConfig, GroundCropConfig, SceneConfig,
    CameraConfig, SolarConfig, SolverConfig, AerosolConfig
)

def load_config(cfg_path: str):
    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f)

    out = OutputConfig(**cfg["output"])
    trajectory_cfg = cfg["sensor"].get("trajectory", {})
    sen = SensorConfig(
        type=cfg["sensor"]["type"],
        altitude_km=float(cfg["sensor"]["altitude_km"]),
        fov_deg=float(cfg["sensor"]["fov_deg"]),
        resolution_km=float(cfg["sensor"]["resolution_km"]),
        views_names=list(cfg["sensor"]["views"]["names"]),
        views_zenith_deg=[float(v) for v in cfg["sensor"]["views"]["zenith_deg"]],
        views_azimuth_deg=[float(v) for v in cfg["sensor"]["views"]["azimuth_deg"]],
        trajectory_mode=str(trajectory_cfg.get("mode", "auto")),
        fallback_heading_deg=float(trajectory_cfg.get("fallback_heading_deg", 0.0)),
        manual_flight_azimuth_deg=(
            None
            if trajectory_cfg.get("manual_flight_azimuth_deg", None) is None
            else float(trajectory_cfg.get("manual_flight_azimuth_deg"))
        ),
        camera_relative_roll_deg=float(trajectory_cfg.get("camera_relative_roll_deg", 0.0)),
        camera_align_with_flight_heading=bool(trajectory_cfg.get("camera_align_with_flight_heading", False)),
        apply_flight_azimuth_offset_to_vaa=bool(trajectory_cfg.get("apply_flight_azimuth_offset_to_vaa", False)),
        camera_image_transpose=bool(trajectory_cfg.get("camera_image_transpose", False)),
        camera_image_flip_lr=bool(trajectory_cfg.get("camera_image_flip_lr", False)),
    )
    bnd = BandsConfig(
        wavelength_nm=[int(w) for w in cfg["bands"]["wavelength_nm"]],
        is_polarized=[bool(v) for v in cfg["bands"]["is_polarized"]],
    )
    dsm = DownsampleConfig(**cfg["downsample"])
    svr = SaveVersionsConfig(versions=list(cfg["save_versions"]))
    plt_cfg = PlotConfig(**cfg["plot"])
    grd = GroundGridConfig(**cfg["ground_grid"])
    comp = ComputeConfig(**cfg.get("compute", {}))
    gcr = GroundCropConfig(
        x_range=list(cfg.get("ground_crop", {}).get("x_range", [0.0, 10.0])),
        y_range=list(cfg.get("ground_crop", {}).get("y_range", [0.0, 6.0]))
    )

    scene_cfg = SceneConfig(**cfg.get("scene", {
        "input_path": "../data/synthetic_cloud_fields/jpl_les/firex_case21_test143_multiview-single-layer_coarse_remove_singlepixeledge.txt",
        "lookat_center_km": [6.0, 1.0, 0.0],
    }))
    cam_cfg = CameraConfig(**cfg.get("camera", {"x_resolution": 2000, "y_resolution": 1200}))
    solar_cfg = SolarConfig(**cfg.get("solar", {"sza_deg": 35.0, "saa_deg": 325.0}))
    solver_cfg = SolverConfig(**cfg.get("solver", {}))
    aerosol_cfg = AerosolConfig(**cfg.get("aerosol", {}))

    return (cfg, out, sen, bnd, dsm, svr, plt_cfg, grd, comp, gcr,
            scene_cfg, cam_cfg, solar_cfg, solver_cfg, aerosol_cfg)


def ensure_dirs(root: str, versions: List[str]):
    os.makedirs(root, exist_ok=True)
    subdirs = {}
    for v in versions:
        p = os.path.join(root, v)
        os.makedirs(p, exist_ok=True)
        subdirs[v] = p
    return subdirs


def crop_to_multiple(a: np.ndarray, factor: int) -> np.ndarray:
    ny, nx = a.shape
    new_ny = (ny // factor) * factor
    new_nx = (nx // factor) * factor
    if new_ny <= 0 or new_nx <= 0:
        raise ValueError(f"Image too small for factor {factor}: shape={a.shape}")
    return a[:new_ny, :new_nx]


def downsample_block(a, factor: int, method: str = "mean") -> np.ndarray:
    # robust: accept xarray or numpy
    if not isinstance(a, np.ndarray):
        a = np.asarray(a)
    a2 = crop_to_multiple(a, factor)
    if method == "mean":
        return block_reduce(a2, block_size=(factor, factor), func=np.nanmean)
    else:
        raise ValueError(f"Unsupported downsample method: {method}")


def compute_column_aod(medium): 
    AOD = {} 
    for key, dataset in medium.items(): 
        if "extinction" in dataset: 
            ext = dataset['extinction'].values # (x,y,z) 
            z = dataset['z'].values # (z) 
            AOD[key] = np.trapz(ext, z, axis=2) 
            ext_min = np.nanmin(ext); 
            ext_max = np.nanmax(ext) 
            print("ext range (3D):", ext_min, ext_max) 
        # total AOD 
        AOD['total'] = sum(v for v in AOD.values()) 
        return AOD

def compute_column_ssa(medium):
    SSA = {}

    for key, dataset in medium.items():
        if "extinction" in dataset:
            ext = dataset['extinction'].values      # (x,y,z)
            ssa = dataset['ssalb'].values           # (x,y,z)
            z = dataset['z'].values                 # (z)

            num = np.trapz(ext * ssa, z, axis=2)
            den = np.trapz(ext, z, axis=2)
            SSA[key] = np.where(den > 0, num / den, np.nan)

    # total SSA (weighted)
    if "cloud" in SSA and "rayleigh" in SSA:
        # weighting by extinction contribution
        ext_total = sum(medium[k]['extinction'].values for k in SSA.keys())
        z = medium['cloud']['z'].values

        num_tot = 0
        for k in SSA:
            num_tot += np.trapz(
                medium[k]['extinction'].values * medium[k]['ssalb'].values,
                z, axis=2
            )
        den_tot = np.trapz(ext_total, z, axis=2)

        SSA["total"] = np.where(den_tot > 0, num_tot / den_tot, np.nan)

    return SSA

def compute_cloud_column_ssa(medium, tau_min=1e-4):
    """
    Compute column-integrated SSA for cloud ONLY.

    Parameters
    ----------
    medium : dict
        medium['cloud'] must contain 'extinction', 'ssalb', 'z'
    tau_min : float
        Minimum column optical depth threshold

    Returns
    -------
    SSA_cloud : 2D ndarray (x, y)
    tau_cloud : 2D ndarray (x, y)
    """

    aerosol = medium['aerosol']

    ext = aerosol['extinction'].values      # (x, y, z)
    ssa = aerosol['ssalb'].values           # (x, y, z)
    z   = aerosol['z'].values               # (z)

    # --- column optical depth ---
    tau = np.trapz(ext, z, axis=2)
    ssa_min = np.nanmin(ssa); ssa_max = np.nanmax(ssa)

    print("ssalb range (3D):", ssa_min, ssa_max)
    print("tau range (2D):", np.nanmin(tau), np.nanmax(tau))
    
    SSA_col = np.trapz(ext*ssa, z, axis=2) / np.maximum(tau, 1e-12)
    bad = (tau > 1e-3) & (SSA_col < 0.9)
    print("bad pixels:", np.sum(bad), "of", bad.size)

    SSA_aerosol = np.full_like(tau, np.nan)

    mask = tau > tau_min
    SSA_aerosol[mask] = (
        np.trapz(ext * ssa, z, axis=2)[mask] /
        tau[mask]
    )

    return SSA_aerosol, tau

def plot_field(field, title, save_path=None, vmin=0.95, vmax=1.0):
    plt.figure(figsize=(4, 3))
    im = plt.imshow(
        field.T,
        origin="lower",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax
    )
    plt.colorbar(im)
    plt.title(title)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved: {save_path}")

    plt.close()
    
def build_multiband_xarray(AOD_all, SSA_all):

    wavelengths = sorted(AOD_all.keys())  # e.g. [355, 380, 445, ...]

    # 读取空间维度 (假设所有波段一样)
    sample = next(iter(AOD_all.values()))
    ny, nx = sample["total"].shape

    # 创建一个 dataset
    ds = xr.Dataset(
        {
            "AOD_total": (["wavelength", "y", "x"], 
                np.stack([AOD_all[w]["total"] for w in wavelengths])),

            "AOD_aerosol": (["wavelength", "y", "x"],
                np.stack([AOD_all[w]["aerosol"] for w in wavelengths])),

            # "AOD_rayleigh": (["wavelength", "y", "x"],
            #     np.stack([AOD_all[w]["rayleigh"] for w in wavelengths])),

            # "SSA_total": (["wavelength", "y", "x"],
            #     np.stack([SSA_all[w]["total"] for w in wavelengths])),

            "SSA_aerosol": (["wavelength", "y", "x"],
                np.stack([SSA_all[w]["aerosol"] for w in wavelengths])),

            # "SSA_rayleigh": (["wavelength", "y", "x"],
            #     np.stack([SSA_all[w]["rayleigh"] for w in wavelengths])),
        },

        coords={
            "wavelength": wavelengths,
            "x": np.arange(nx),
            "y": np.arange(ny),
        }
    )

    return ds
