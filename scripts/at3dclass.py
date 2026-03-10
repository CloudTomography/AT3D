#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:33:43 2025

@author: benting
"""
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


@dataclass
class OutputConfig:
    root_dir: str
    save_png: bool
    save_nc: bool
    plot_mode: str  # "skip" | "if_missing" | "overwrite"
    png_dpi: int = 300

@dataclass
class SensorConfig:
    type: str  # "perspective_projection" or "orthographic_projection"
    altitude_km: float
    fov_deg: float
    resolution_km: float
    views_names: List[str]
    views_zenith_deg: List[float]
    views_azimuth_deg: List[float]

@dataclass
class BandsConfig:
    wavelength_nm: List[int]
    is_polarized: List[bool]

@dataclass
class DownsampleConfig:
    factor: int
    method: str  # "mean"

@dataclass
class SaveVersionsConfig:
    versions: List[str]

@dataclass
class PlotConfig:
    enable_3d_geometry: bool
    enable_ground_image: bool
    colormap: str

@dataclass
class GroundGridConfig:
    x_min: float
    x_max: float
    nx: int
    y_min: float
    y_max: float
    ny: int

@dataclass
class ComputeConfig:
    n_jobs: Optional[int] = None

@dataclass
class GroundCropConfig:
    x_range: List[float]
    y_range: List[float]

@dataclass
class SceneConfig:
    input_path: str
    lookat_center_km: List[float]

@dataclass
class CameraConfig:
    x_resolution: int
    y_resolution: int

@dataclass
class SolarConfig:
    sza_deg: float
    saa_deg: float

@dataclass
class SolverConfig:
    x_boundary_condition: str = "open"
    y_boundary_condition: str = "open"
    num_mu_bins: int = 48
    num_phi_bins: int = 96
    split_accuracy: float = 0.01
    deltam: bool = True

@dataclass
class AerosolConfig:
    refractive_index_real: float = 1.4
    refractive_index_imag: float = 0.001
    particle_density: float = 1.6
