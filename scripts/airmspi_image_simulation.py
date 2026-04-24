#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SimulatingRadiances_v5b.py
----------------------------------
Per-wavelength AT3D simulation pipeline with:
- Original camera-view PNGs via imshow()
- Ground-projected, cropped terrain PNGs using user's crop/plot functions
- Configurable ground crop ranges via YAML (ground_crop.x_range / y_range)
- Same NPZ / NC outputs as v3c, plus expanded metadata (angles, lat/lon, elevation, mask)
- NEW: Read lat/lon columns from wrf_to_shdom_coarse.txt and propagate to outputs

Created on Thu Nov 20 13:19:59 2025

Author: Benting Chen
"""

import os
import sys
import time
import tempfile
import yaml
import argparse
import json
import hashlib
import subprocess
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple, Optional
from datetime import datetime

# External libs
import at3d
from scipy.ndimage import map_coordinates
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import utils
from at3dclass import (
    OutputConfig, SensorConfig, BandsConfig, DownsampleConfig, SaveVersionsConfig,
    PlotConfig, GroundGridConfig, ComputeConfig, GroundCropConfig, SceneConfig,
    CameraConfig, SolarConfig, SolverConfig, AerosolConfig
)
import cartopy.crs as ccrs
from scipy.ndimage import distance_transform_edt
from scipy.interpolate import RegularGridInterpolator

_MIE_TABLE_CACHE: Dict[Tuple[float, float, float], Any] = {}
# =============================
#%% Section 1: Config dataclasses
# =============================

# @dataclass
# class OutputConfig:
#     root_dir: str
#     save_png: bool
#     save_nc: bool
#     plot_mode: str  # "skip" | "if_missing" | "overwrite"
#     png_dpi: int = 300

# @dataclass
# class SensorConfig:
#     type: str  # "perspective_projection" or "orthographic_projection"
#     altitude_km: float
#     fov_deg: float
#     resolution_km: float
#     views_names: List[str]
#     views_zenith_deg: List[float]
#     views_azimuth_deg: List[float]

# @dataclass
# class BandsConfig:
#     wavelength_nm: List[int]
#     is_polarized: List[bool]

# @dataclass
# class DownsampleConfig:
#     factor: int
#     method: str  # "mean"

# @dataclass
# class SaveVersionsConfig:
#     versions: List[str]

# @dataclass
# class PlotConfig:
#     enable_3d_geometry: bool
#     enable_ground_image: bool
#     colormap: str

# @dataclass
# class GroundGridConfig:
#     x_min: float
#     x_max: float
#     nx: int
#     y_min: float
#     y_max: float
#     ny: int

# @dataclass
# class ComputeConfig:
#     n_jobs: Optional[int] = None

# @dataclass
# class GroundCropConfig:
#     x_range: List[float]
#     y_range: List[float]


# =============================
#%% Section 2: Utilities
# =============================


# =============================
#%% Section 3: Geometry / Registration
# =============================

def calculate_up_vector(position, look_at_point, world_up=np.array([0, 0, 1])):
    forward_vec = look_at_point - position
    forward_vec = forward_vec / np.linalg.norm(forward_vec)
    if np.abs(np.dot(forward_vec, world_up)) > 0.9999:
        right_vec = np.array([1, 0, 0])
        up_vec = np.cross(right_vec, forward_vec)
        if np.linalg.norm(up_vec) < 1e-6:
            right_vec = np.array([0, -1, 0])
            up_vec = np.cross(right_vec, forward_vec)
    else:
        right_vec = np.cross(forward_vec, world_up)
        right_vec = right_vec / np.linalg.norm(right_vec)
        up_vec = np.cross(right_vec, forward_vec)
    return up_vec

def _rotate_vector_about_axis(v: np.ndarray, axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rodrigues rotation."""
    ang = np.deg2rad(float(angle_deg))
    k = np.asarray(axis, dtype=float)
    kn = np.linalg.norm(k)
    if kn < 1e-12 or abs(ang) < 1e-12:
        return np.asarray(v, dtype=float)
    k = k / kn
    vv = np.asarray(v, dtype=float)
    return (vv * np.cos(ang) +
            np.cross(k, vv) * np.sin(ang) +
            k * np.dot(k, vv) * (1.0 - np.cos(ang)))

def calculate_sensor_trajectory(sensor_zenith_list,
                                sensor_azimuth_list,
                                look_at_point=np.array([0.0, 0.0, 0.0]),
                                sensor_altitude=20.0,
                                trajectory_mode="adjacent_views",
                                fallback_heading_deg=0.0,
                                manual_flight_azimuth_deg=None,
                                camera_relative_roll_deg=0.0,
                                camera_align_with_flight_heading=False):
    """
    Build sensor positions and up-vectors.

    trajectory_mode:
      - "adjacent_views": infer heading from neighboring view positions.
      - "sensor_azimuth": use sensor_azimuth_list as per-view heading.
    """
    positions = []
    up_vectors = []
    for zenith, azimuth in zip(sensor_zenith_list, sensor_azimuth_list):
        zen = np.deg2rad(zenith)
        azi = np.deg2rad(azimuth)
        v = np.array([
            np.sin(zen) * np.cos(azi),
            np.sin(zen) * np.sin(azi),
            -np.cos(zen)
        ])
        R = sensor_altitude / np.cos(zen)
        pos = look_at_point - R * v
        pos[2] = sensor_altitude
        positions.append(pos)
    positions = np.array(positions)

    n_view = len(positions)
    if n_view == 0:
        raise ValueError("No sensor views provided.")

    mode = str(trajectory_mode).lower()
    if mode == "manual_azimuth":
        # backward-compatible alias
        mode = "sensor_azimuth"

    if n_view < 2:
        # 单视角：直接使用该视角 azimuth（用户要求）
        heading = np.deg2rad(float(sensor_azimuth_list[0] if len(sensor_azimuth_list) > 0 else fallback_heading_deg))
        traj_dirs = np.array([[np.cos(heading), np.sin(heading), 0.0]], dtype=float)
    elif mode in {"adjacent_views", "auto", "neighbor_views"}:
        traj = np.gradient(positions, axis=0)
        traj_dirs = []
        for d in traj:
            d = np.asarray(d, dtype=float)
            d[2] = 0.0
            nrm = np.linalg.norm(d)
            if nrm < 1e-9:
                heading = np.deg2rad(float(fallback_heading_deg))
                d = np.array([np.cos(heading), np.sin(heading), 0.0], dtype=float)
            else:
                d = d / nrm
            traj_dirs.append(d)
        traj_dirs = np.asarray(traj_dirs)
    elif mode in {"sensor_azimuth", "sensor_azimuth_list"}:
        traj_dirs = []
        for az in sensor_azimuth_list:
            heading = np.deg2rad(float(az))
            traj_dirs.append(np.array([np.cos(heading), np.sin(heading), 0.0], dtype=float))
        traj_dirs = np.asarray(traj_dirs)
    else:
        raise ValueError("Unsupported trajectory_mode. Use 'adjacent_views' or 'sensor_azimuth'.")

    for pos, traj_dir in zip(positions, traj_dirs):
        look_dir = look_at_point - pos
        look_dir = look_dir / np.linalg.norm(look_dir)
        if camera_align_with_flight_heading:
            up_vec = np.cross(traj_dir, look_dir)
            up_norm = np.linalg.norm(up_vec)
            if up_norm < 1e-9:
                up_vec = calculate_up_vector(pos, look_at_point)
            else:
                up_vec = up_vec / up_norm
        else:
            # World-aligned camera frame: ensures no implicit heading-driven image rotation.
            up_vec = calculate_up_vector(pos, look_at_point)
        # Camera-vs-aircraft relative angle (roll around look direction).
        if abs(float(camera_relative_roll_deg)) > 1e-12:
            up_vec = _rotate_vector_about_axis(up_vec, look_dir, camera_relative_roll_deg)
            up_vec = up_vec / max(np.linalg.norm(up_vec), 1e-12)
        up_vectors.append(up_vec)

    return positions, up_vectors


def _rotation_matrix_x(angle_deg):
    a = np.deg2rad(angle_deg)
    return np.array([
        [1, 0, 0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a),  np.cos(a)]
    ])


def _rotation_matrix_y(angle_deg):
    a = np.deg2rad(angle_deg)
    return np.array([
        [ np.cos(a), 0, np.sin(a)],
        [0,          1, 0],
        [-np.sin(a), 0, np.cos(a)]
    ])


def _rotation_matrix_z(angle_deg):
    # heading: cw from +y → convert to math (ccw from +x)
    a = np.deg2rad(90 - angle_deg)
    return np.array([
        [np.cos(a), -np.sin(a), 0],
        [np.sin(a),  np.cos(a), 0],
        [0,          0,         1]
    ])


def calculate_sensor_trajectory_from_aircraft(
        look_at_point=np.array([0.0, 0.0, 0.0]),
        sensor_altitude=20.0,

        heading_angle_deg=0.0,
        pitch_angle_deg=0.0,
        roll_angle_deg=0.0,

        camera_pitch_relative_deg=0.0,   # 沿机身方向前后摆（关键参数）
        camera_roll_relative_deg=0.0,    # 相机自身roll

        n_views=1
):
    """
    Aircraft-based sensor geometry.

    Assumed body frame:
      x = forward
      y = right
      z = up

    Camera nominal nadir-looking:
      look = -z

    Image 'up' direction in body frame is chosen as +x
    (i.e. top of image points toward aircraft forward direction when no extra roll)
    """

    positions = []
    up_vectors = []

    # --- 1. camera basis in body frame ---
    cam_look_body = np.array([0.0, 0.0, -1.0])
    cam_up_body = np.array([1.0, 0.0, 0.0])   # image top = aircraft forward

    # --- 2. camera relative rotation ---
    # first tilt camera relative to aircraft
    R_cam_pitch = _rotation_matrix_y(camera_pitch_relative_deg)
    cam_look_body = R_cam_pitch @ cam_look_body
    cam_up_body = R_cam_pitch @ cam_up_body

    # then camera self-roll around LOS
    if abs(camera_roll_relative_deg) > 1e-12:
        cam_up_body = _rotate_vector_about_axis(
            cam_up_body, cam_look_body, camera_roll_relative_deg
        )

    # --- 3. aircraft attitude rotation ---
    R_roll = _rotation_matrix_x(roll_angle_deg)
    R_pitch = _rotation_matrix_y(pitch_angle_deg)
    R_heading = _rotation_matrix_z(heading_angle_deg)
    R_aircraft = R_heading @ R_pitch @ R_roll

    # --- 4. transform to world ---
    look_dir_world = R_aircraft @ cam_look_body
    up_vec_world = R_aircraft @ cam_up_body

    look_dir_world = look_dir_world / np.linalg.norm(look_dir_world)

    # make up vector exactly perpendicular to look direction
    up_vec_world = up_vec_world - np.dot(up_vec_world, look_dir_world) * look_dir_world
    up_vec_world = up_vec_world / np.linalg.norm(up_vec_world)

    # --- 5. sensor position ---
    # line from sensor to look_at_point follows look_dir_world
    # sensor = look_at_point - s * look_dir_world
    if abs(look_dir_world[2]) < 1e-12:
        raise ValueError("Look direction is horizontal; cannot intersect z=sensor_altitude plane.")

    s = (look_at_point[2] - sensor_altitude) / look_dir_world[2]
    pos = look_at_point - s * look_dir_world

    # ensure exact altitude
    pos[2] = sensor_altitude

    # --- 6. replicate for n_views ---
    for _ in range(n_views):
        positions.append(pos.copy())
        up_vectors.append(up_vec_world.copy())

    return np.array(positions), np.array(up_vectors)


def calculate_sensor_trajectory_cross_track(
        n_views,
        x1, y1, z1,
        x2, y2, z2,
        spacing,
        scan1_deg,
        scan2_deg,
        delscan_deg):
    """
    SHDOM-like cross track geometry sampler.
    A sequence of camera positions is generated from (x1,y1,z1) to (x2,y2,z2),
    and each scan applies a cross-track rotation around the along-track axis.
    If `n_views` is None all generated (scan_position, scan_angle) samples are returned.
    """
    if n_views is not None and n_views <= 0:
        raise ValueError("n_views must be > 0 for cross track mode.")
    if abs(float(delscan_deg)) < 1e-12:
        raise ValueError("cross_track_delscan_deg cannot be 0.")

    p1 = np.array([x1, y1, z1], dtype=float)
    p2 = np.array([x2, y2, z2], dtype=float)
    track_vec = p2 - p1
    track_len = np.linalg.norm(track_vec)
    if track_len < 1e-12:
        raise ValueError("cross track start/end positions are identical.")
    along_track = track_vec / track_len

    if spacing <= 0.0:
        scan_positions = np.array([p1], dtype=float)
    else:
        n_scans = int(np.floor(track_len / spacing)) + 1
        dists = np.arange(n_scans, dtype=float) * spacing
        dists = np.clip(dists, 0.0, track_len)
        scan_positions = p1[None, :] + dists[:, None] * along_track[None, :]
        if np.linalg.norm(scan_positions[-1] - p2) > 1e-8:
            scan_positions = np.vstack([scan_positions, p2[None, :]])

    if np.sign(scan2_deg - scan1_deg) == np.sign(delscan_deg):
        scan_angles = np.arange(scan1_deg, scan2_deg + 0.5 * delscan_deg, delscan_deg, dtype=float)
    else:
        scan_angles = np.array([scan1_deg], dtype=float)

    samples = [(pos, ang) for pos in scan_positions for ang in scan_angles]
    if len(samples) == 0:
        raise ValueError("No valid cross track samples were generated.")
    if n_views is None:
        selected_samples = samples
    else:
        if len(samples) < n_views:
            reps = int(np.ceil(n_views / len(samples)))
            samples = samples * reps
        selected_samples = samples[:n_views]

    positions = []
    lookat_vectors = []
    up_vectors = []
    base_look = np.array([0.0, 0.0, -1.0], dtype=float)
    world_up = np.array([0.0, 0.0, 1.0], dtype=float)
    for pos, ang in selected_samples:
        # SHDOM cross-track convention:
        # scan angle sign is defined in scanner coordinates; in this NEU setup we
        # need the opposite rotation sign to keep left/right VAA ordering consistent.
        look_dir = _rotate_vector_about_axis(base_look, along_track, -float(ang))
        look_dir = look_dir / np.linalg.norm(look_dir)
        lookat = np.asarray(pos, dtype=float) + look_dir

        right = np.cross(look_dir, along_track)
        if np.linalg.norm(right) < 1e-12:
            right = np.cross(look_dir, world_up)
        up_vec = np.cross(right, look_dir)
        if np.linalg.norm(up_vec) < 1e-12:
            up_vec = calculate_up_vector(np.asarray(pos, dtype=float), lookat)
        else:
            up_vec = up_vec / np.linalg.norm(up_vec)

        positions.append(np.asarray(pos, dtype=float))
        lookat_vectors.append(np.asarray(lookat, dtype=float))
        up_vectors.append(np.asarray(up_vec, dtype=float))
    return np.asarray(positions), np.asarray(lookat_vectors), np.asarray(up_vectors)


def cross_track_scan_projection(
        wavelength,
        stokes,
        x1, y1, z1,
        x2, y2, z2,
        spacing,
        scan1_deg,
        scan2_deg,
        delscan_deg,
        pitch_start_deg=None,
        pitch_end_deg=None,
        pitch_list_deg=None):
    """
    Build a single cross-track scan sensor without perspective projection.
    Pixels are organized as [scan_index, cross_track_angle_index].
    """
    if abs(float(delscan_deg)) < 1e-12:
        raise ValueError("cross_track_delscan_deg cannot be 0.")
    p1 = np.array([x1, y1, z1], dtype=float)
    p2 = np.array([x2, y2, z2], dtype=float)
    track_vec = p2 - p1
    track_len = np.linalg.norm(track_vec)
    if track_len < 1e-12:
        raise ValueError("cross track start/end positions are identical.")
    along_track = track_vec / track_len

    if spacing <= 0.0:
        scan_positions = np.array([p1], dtype=float)
    else:
        n_scans = int(np.floor(track_len / spacing)) + 1
        dists = np.arange(n_scans, dtype=float) * spacing
        dists = np.clip(dists, 0.0, track_len)
        scan_positions = p1[None, :] + dists[:, None] * along_track[None, :]
        if np.linalg.norm(scan_positions[-1] - p2) > 1e-8:
            scan_positions = np.vstack([scan_positions, p2[None, :]])

    if np.sign(scan2_deg - scan1_deg) == np.sign(delscan_deg):
        scan_angles = np.arange(scan1_deg, scan2_deg + 0.5 * delscan_deg, delscan_deg, dtype=float)
    else:
        scan_angles = np.array([scan1_deg], dtype=float)
    if scan_angles.size == 0:
        raise ValueError("No cross-track scan angles generated.")

    base_look = np.array([0.0, 0.0, -1.0], dtype=float)
    world_up = np.array([0.0, 0.0, 1.0], dtype=float)
    scan_pitch_deg = _resolve_cross_track_scan_pitch_angles(
        n_scans=scan_positions.shape[0],
        pitch_start_deg=pitch_start_deg,
        pitch_end_deg=pitch_end_deg,
        pitch_list_deg=pitch_list_deg,
    )
    x, y, z, mu, phi = [], [], [], [], []
    for iscan, pos in enumerate(scan_positions):
        pitch_deg = float(scan_pitch_deg[iscan])
        for ang in scan_angles:
            # Keep sign convention consistent with calculate_sensor_trajectory_cross_track.
            look_dir = _rotate_vector_about_axis(base_look, along_track, -float(ang))
            right = np.cross(look_dir, along_track)
            if np.linalg.norm(right) < 1e-12:
                right = np.cross(look_dir, world_up)
            right = right / max(np.linalg.norm(right), 1e-12)
            if abs(pitch_deg) > 1e-12:
                look_dir = _rotate_vector_about_axis(look_dir, right, pitch_deg)
            look_dir = look_dir / max(np.linalg.norm(look_dir), 1e-12)  # sensor->scene
            v_out = -look_dir  # scene->sensor
            x.append(pos[0]); y.append(pos[1]); z.append(pos[2])
            mu.append(v_out[2])
            phi.append((np.arctan2(v_out[1], v_out[0]) + 2*np.pi) % (2*np.pi))

    sensor = shdom_cross_track_sensor_wrapper(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        z=np.asarray(z, dtype=float),
        mu=np.asarray(mu, dtype=float),
        phi=np.asarray(phi, dtype=float),
        stokes=stokes,
        wavelength=wavelength,
        fill_ray_variables=True,
        image_shape=(int(scan_angles.size), int(scan_positions.shape[0]))
    )
    sensor.attrs.update({
        'projection': 'CrossTrackScan',
        'x_resolution': int(scan_angles.size),
        'y_resolution': int(scan_positions.shape[0]),
        'cross_track_scan_angles_deg': scan_angles.astype(float).tolist(),
        'cross_track_scan_positions': scan_positions.astype(float).tolist(),
        'cross_track_scan_pitch_deg': scan_pitch_deg.astype(float).tolist(),
    })
    return sensor, scan_positions, scan_angles, scan_pitch_deg


def shdom_cross_track_sensor_wrapper(
        x, y, z, mu, phi, stokes, wavelength, fill_ray_variables=True, image_shape=None):
    """
    Wrapper-style sensor builder that directly constructs a SHDOM/AT3D sensor dataset.
    """
    sensor = at3d.sensor.make_sensor_dataset(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        z=np.asarray(z, dtype=float),
        mu=np.asarray(mu, dtype=float),
        phi=np.asarray(phi, dtype=float),
        stokes=stokes,
        wavelength=wavelength,
        fill_ray_variables=bool(fill_ray_variables)
    )
    if image_shape is not None:
        nx, ny = int(image_shape[0]), int(image_shape[1])
        sensor['image_shape'] = xr.DataArray(
            [nx, ny],
            coords={'image_dims': ['nx', 'ny']},
            dims='image_dims'
        )
    return sensor


def _resolve_cross_track_scan_pitch_angles(n_scans, pitch_start_deg=None, pitch_end_deg=None, pitch_list_deg=None):
    """
    Resolve along-track per-scan pitch angles (deg).
    Priority:
      1) pitch_list_deg (manual list, length must equal n_scans),
      2) pitch_start_deg + pitch_end_deg (linear interpolation),
      3) zeros.
    """
    if pitch_list_deg is not None:
        pitch_arr = np.asarray(pitch_list_deg, dtype=float).reshape(-1)
        if pitch_arr.size != int(n_scans):
            raise ValueError(
                f"cross_track_pitch_list_deg length {pitch_arr.size} != n_scans {int(n_scans)}."
            )
        return pitch_arr
    if (pitch_start_deg is None) ^ (pitch_end_deg is None):
        raise ValueError("cross_track_pitch_start_deg and cross_track_pitch_end_deg must be both set.")
    if pitch_start_deg is not None and pitch_end_deg is not None:
        return np.linspace(float(pitch_start_deg), float(pitch_end_deg), int(n_scans), dtype=float)
    return np.zeros(int(n_scans), dtype=float)
def project_to_ground_lookat(img, pos, look_at_point, up_vector, fov_diag_deg, Xg, Yg):
    ny, nx = img.shape
    img_ground = np.full(Xg.shape, np.nan, dtype=np.float32)
    fov_diag_rad = np.deg2rad(fov_diag_deg)
    z_cam = np.array(look_at_point) - np.array(pos); z_cam /= np.linalg.norm(z_cam)
    x_cam = np.cross(np.array(up_vector), z_cam)
    if np.linalg.norm(x_cam) < 1e-9:
        world_x = np.array([1.0, 0.0, 0.0]); x_cam = np.cross(z_cam, world_x)
        if np.linalg.norm(x_cam) < 1e-9:
            world_y = np.array([0.0, 1.0, 0.0]); x_cam = np.cross(world_y, z_cam)
    x_cam /= np.linalg.norm(x_cam)
    y_cam = np.cross(z_cam, x_cam)
    rotation_matrix = np.array([x_cam, y_cam, z_cam])
    P_ground = np.stack([Xg, Yg, np.zeros_like(Xg)], axis=-1)
    pos_arr = np.array(pos).reshape(1, 1, 3)
    vec_world = P_ground - pos_arr
    norm = np.linalg.norm(vec_world, axis=-1, keepdims=True); norm[norm == 0] = 1e-9
    vec_world_norm = vec_world / norm
    vec_cam = np.einsum('ij,hwj->hwi', rotation_matrix, vec_world_norm)
    vx, vy, vz = vec_cam[..., 0], vec_cam[..., 1], vec_cam[..., 2]
    cos_half_fov = np.cos(fov_diag_rad / 2.0)
    valid_mask = vz > cos_half_fov
    u = np.full_like(vx, np.nan); v = np.full_like(vy, np.nan)
    u[valid_mask] = vx[valid_mask] / vz[valid_mask]
    v[valid_mask] = vy[valid_mask] / vz[valid_mask]
    aspect_ratio = nx / ny
    tan_half_fov_diag = np.tan(fov_diag_rad / 2.0)
    tan_half_fov_h = tan_half_fov_diag * aspect_ratio / np.sqrt(aspect_ratio**2 + 1)
    tan_half_fov_v = np.tan(fov_diag_rad / 2.0) / np.sqrt(aspect_ratio**2 + 1)
    px = (u / tan_half_fov_h) * (nx / 2.0) + (nx / 2.0 - 0.5)
    py = (-v / tan_half_fov_v) * (ny / 2.0) + (ny / 2.0 - 0.5)
    valid_indices = np.where(valid_mask)
    coords_to_sample = np.vstack((py[valid_indices], px[valid_indices]))
    sampled_values = map_coordinates(img, coords_to_sample, order=1, mode='constant', cval=np.nan)
    img_ground[valid_indices] = sampled_values
    return img_ground

def reproject_to_ground(sensor, ground_z=0.0):
    """
    Reproject camera image pixels to ground coordinates (assuming a flat ground plane z = ground_z).

    Parameters
    ----------
    sensor : xr.Dataset
        The dataset output from `perspective_projection()`.
    ground_z : float
        The z (Up) coordinate of the ground plane [km]. Default is 0.

    Returns
    -------
    ground_coords : xr.Dataset
        Dataset containing the ground-projected (x, y, z=ground_z) coordinates
        corresponding to each image pixel.
    """

    # ---- 从 sensor 中恢复几何信息 ----
    nx = int(sensor.attrs['x_resolution'])
    ny = int(sensor.attrs['y_resolution'])
    v_out_map = compute_vout_map_from_sensor(sensor)
    rays_world = -v_out_map.reshape(-1, 3).T  # sensor -> scene

    cam_x = np.asarray(sensor.cam_x.data, dtype=float).reshape(-1)
    cam_y = np.asarray(sensor.cam_y.data, dtype=float).reshape(-1)
    cam_z = np.asarray(sensor.cam_z.data, dtype=float).reshape(-1)

    # ---- 计算与地面 (z=ground_z) 的交点 ----
    dir_z = rays_world[2, :]
    t = (ground_z - cam_z) / dir_z  # 每个射线到地面的比例因子

    # 仅保留“向下看且在相机前方”的有效交点：
    # - dir_z < 0: 射线朝向地面
    # - t > 0: 交点在相机前方
    valid = (dir_z < 0) & (t > 0)

    x_ground = np.full_like(t, np.nan, dtype=float)
    y_ground = np.full_like(t, np.nan, dtype=float)
    z_ground = np.full_like(t, np.nan, dtype=float)

    x_ground[valid] = cam_x[valid] + t[valid] * rays_world[0, valid]
    y_ground[valid] = cam_y[valid] + t[valid] * rays_world[1, valid]
    z_ground[valid] = ground_z

    # ---- 整理为 DataArray ----
    x_ground = x_ground.reshape(ny, nx)
    y_ground = y_ground.reshape(ny, nx)
    z_ground = z_ground.reshape(ny, nx)

    ground_coords = xr.Dataset({
        "x_ground": (("y", "x"), x_ground),
        "y_ground": (("y", "x"), y_ground),
        "z_ground": (("y", "x"), z_ground)
    })

    ground_coords.attrs = {
        "description": "Ground-projected coordinates from camera pixels",
        "ground_z": ground_z,
        "camera_position": [float(np.nanmean(cam_x)), float(np.nanmean(cam_y)), float(np.nanmean(cam_z))]
    }

    return ground_coords

def compute_vout_map_from_sensor(sensor_ds):
    """
    Returns v_out_map (ny, nx, 3): propagation direction from scene -> sensor
    in world coordinates (x North, y East, z Up).
    """
    nx = int(sensor_ds.attrs['x_resolution'])
    ny = int(sensor_ds.attrs['y_resolution'])
    mu = np.asarray(sensor_ds.cam_mu.data, dtype=float).reshape(ny, nx)
    phi = np.asarray(sensor_ds.cam_phi.data, dtype=float).reshape(ny, nx)
    sin_theta = np.sqrt(np.clip(1.0 - mu**2, 0.0, 1.0))
    v_out = np.stack([
        sin_theta * np.cos(phi),
        sin_theta * np.sin(phi),
        mu
    ], axis=-1)
    return v_out


def _get_flight_azimuth_offset_deg_from_context(context_cfg: Any) -> float:
    """
    Flight heading angle used to rotate VAA into flight-referenced convention.
    Always applied; default heading is 0° (North).
    """
    try:
        if isinstance(context_cfg, dict):
            return float(context_cfg.get("heading_angle_deg", 0.0))
        return float(getattr(context_cfg, "aircraft_heading_deg", 0.0))
    except Exception:
        return 0.0


def _get_image_orientation_from_sensor_cfg(sen_cfg: Any) -> Tuple[bool, bool]:
    transpose = bool(getattr(sen_cfg, "camera_image_transpose", False))
    flip_lr = bool(getattr(sen_cfg, "camera_image_flip_lr", False))
    return transpose, flip_lr


def _apply_image_orientation(arr: np.ndarray, transpose: bool, flip_lr: bool) -> np.ndarray:
    out = np.asarray(arr)
    if transpose:
        if out.ndim == 2:
            out = out.T
        elif out.ndim == 3:
            out = np.transpose(out, (1, 0, 2))
    if flip_lr:
        out = np.fliplr(out)
    return out


def _to_aircraft_eye_view(arr: np.ndarray, flip_vertical: bool = True) -> np.ndarray:
    """
    Convert raw camera image to a WYSIWYG view (as seen from aircraft cockpit).
    By default this is a 180° image rotation; for cross-track mode users may
    request keeping the vertical direction unchanged.
    """
    out = np.asarray(arr)
    if out.ndim < 2:
        return out
    out = np.fliplr(out)
    if flip_vertical:
        out = np.flipud(out)
    return out


def _ensure_2d_shape(arr: np.ndarray, target_shape: Tuple[int, int]) -> np.ndarray:
    """Ensure 2D array matches target shape; allow implicit transpose if needed."""
    out = np.asarray(arr)
    if out.shape == target_shape:
        return out
    if out.ndim == 2 and out.T.shape == target_shape:
        return out.T
    raise ValueError(f"Shape mismatch: got {out.shape}, expected {target_shape}")


def _compute_angle_maps_from_sensor(
    sensor_ds,
    solar_azimuth_deg: float,
    solar_zenith_deg: float,
    heading_angle_deg: float = 0.0,
    apply_heading_offset: bool = False,
    transpose: bool = False,
    flip_lr: bool = False,
):
    """
    Compute VZA/VAA/RAA/Scattering-angle maps from a sensor dataset using one
    unified geometry path (shared by simulation and precheck/replot code).
    """
    v_out_map = compute_vout_map_from_sensor(sensor_ds)
    v_out_map = _apply_image_orientation(v_out_map, transpose, flip_lr)

    vz = np.clip(v_out_map[..., 2], -1.0, 1.0)
    vza_map = np.degrees(np.arccos(vz))

    vaa_map = (np.degrees(np.arctan2(v_out_map[..., 1], v_out_map[..., 0])) + 360.0) % 360.0
    # Camera-image convention: enforce 0° from image center toward "up" (not down).
    # vaa_map = (360 - vaa_map) % 360
    # vaa_map = (vaa_map + 180.0) % 360.0
    if apply_heading_offset:
        vaa_map = ((vaa_map - float(heading_angle_deg) + 360.0) % 360.0)

    saa = (float(solar_azimuth_deg) + 360.0) % 360.0
    sza = float(solar_zenith_deg)
    raa_map = ((vaa_map - saa))

    mu0 = np.cos(np.radians(sza))
    mu = np.cos(np.radians(vza_map))
    cos_sca = -mu0 * mu + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu**2) * np.cos(np.radians(raa_map))
    sca_angle = np.degrees(np.arccos(np.clip(cos_sca, -1.0, 1.0)))
    return vza_map, vaa_map, raa_map, sca_angle

def assign_latlon_from_grid(xg, yg, wrf_x, wrf_y, xlats, xlons):
    """
    将相机重投影地面坐标 (xg, yg) 映射为对应的 (lat, lon)。

    Parameters:
        xg, yg : (ny, nx)
            相机每个像素的地面投影坐标（单位：km）。
        wrf_x, wrf_y : 1D or 2D
            WRF 网格的 x, y 坐标（km），通常为 1D。
        xlats, xlons : (ny_wrf, nx_wrf)
            每个 WRF 网格点的纬度、经度。

    Returns:
        lat_img, lon_img : (ny, nx)
            每个相机像元的 (lat, lon)
    """

    # 如果是 2D 网格点，转换为 1D
    if wrf_x.ndim == 2:
        wrf_x = wrf_x[0, :]
    if wrf_y.ndim == 2:
        wrf_y = wrf_y[:, 0]

    # Normalize lat/lon array orientation for interpolator axes (x, y).
    # Common case is (ny, nx); convert to (nx, ny) by transpose.
    nx = wrf_x.size
    ny = wrf_y.size
    if xlats.shape == (ny, nx):
        lat_field = xlats.T
        lon_field = xlons.T
    elif xlats.shape == (nx, ny):
        lat_field = xlats
        lon_field = xlons
    else:
        raise ValueError(
            f"lat/lon shape {xlats.shape} does not match wrf grid "
            f"(ny,nx)=({ny},{nx}) or (nx,ny)=({nx},{ny})."
        )

    lat_interp = RegularGridInterpolator(
        (wrf_x, wrf_y),
        lat_field,
        bounds_error=False,
        fill_value=np.nan
    )

    lon_interp = RegularGridInterpolator(
        (wrf_x, wrf_y),
        lon_field,
        bounds_error=False,
        fill_value=np.nan
    )
    
    

    # 展开 pixel 坐标
    pts = np.stack([xg.ravel(), yg.ravel()], axis=-1)
    
    
    # plt.scatter(xg, yg, s=1)
    # xxx,yyy = np.meshgrid(wrf_x, wrf_y)
    # plt.scatter(xxx, yyy, s=1,c='red')
    # plt.xlabel("xg")
    # plt.ylabel("yg")
    # plt.gca().set_aspect("equal")
    # plt.show()
    
    

    # 执行插值
    lat_img = lat_interp(pts).reshape(xg.shape)
    
    lon_img = lon_interp(pts).reshape(xg.shape)
    
    # plt.imshow(lat_img)
    # plt.show()
    
    # plt.imshow(lon_img)
    # plt.show()

    return lat_img, lon_img
#%% Section 3.1: User-provided helpers for ground cropping & plotting
# =============================

def centers_to_edges_2d(xc, yc):
    ny, nx = xc.shape
    xv = 0.25*(xc[:-1,:-1] + xc[1:,:-1] + xc[:-1,1:] + xc[1:,1:])
    yv = 0.25*(yc[:-1,:-1] + yc[1:,:-1] + yc[:-1,1:] + yc[1:,1:])
    def _extrap1d(a, axis):
        a1 = np.take(a, indices=0, axis=axis)
        a2 = np.take(a, indices=1, axis=axis)
        a_end = np.take(a, indices=-1, axis=axis)
        a_end2 = np.take(a, indices=-2, axis=axis)
        a0 = 2*a1 - a2
        a_last = 2*a_end - a_end2
        return np.concatenate([np.expand_dims(a0, axis=axis),
                               a,
                               np.expand_dims(a_last, axis=axis)], axis=axis)
    xv = _extrap1d(xv, axis=0)
    yv = _extrap1d(yv, axis=0)
    xv = _extrap1d(xv, axis=1)
    yv = _extrap1d(yv, axis=1)
    return xv, yv

def crop_by_world_box(data, xg, yg, x_range, y_range):
    x1, x2 = sorted(x_range)
    y1, y2 = sorted(y_range)
    mask = (xg >= x1) & (xg <= x2) & (yg >= y1) & (yg <= y2)
    rows = np.where(np.any(mask, axis=1))[0]
    cols = np.where(np.any(mask, axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        raise ValueError("No data within crop range. Check x_range/y_range.")
    row_slice = slice(rows.min(), rows.max() + 1)
    col_slice = slice(cols.min(), cols.max() + 1)
    if hasattr(data, "where") and hasattr(data, "isel"):
        data_c = data.where(mask)
        xg_c = np.where(mask, xg, np.nan)
        yg_c = np.where(mask, yg, np.nan)
        if hasattr(data, "dims"):
            dims = list(data.dims)
            if ("imgdim1" in dims) and ("imgdim0" in dims):
                data_c = data_c.isel(imgdim1=row_slice, imgdim0=col_slice)
            else:
                indexer = {dims[-2]: row_slice, dims[-1]: col_slice}
                data_c = data_c.isel(**indexer)
        else:
            data_c = data[row_slice, col_slice]
        xg_c = xg[row_slice, col_slice]
        yg_c = yg[row_slice, col_slice]
        return data_c, xg_c, yg_c
    data_c = data[row_slice, col_slice]
    xg_c = xg[row_slice, col_slice]
    yg_c = yg[row_slice, col_slice]
    return data_c, xg_c, yg_c

def plot_on_ground(data, xg, yg, title='', cmap='viridis', vmin=None, vmax=None, save_path=None, show=False):
    xv, yv = centers_to_edges_2d(xg, yg)
    fig, ax = plt.subplots(figsize=(7, 6))
    pm = ax.pcolormesh(xv, yv, data, shading='flat', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x_ground [km] (North +)')
    ax.set_ylabel('y_ground [km] (East +)')
    ax.set_title(title)
    fig.colorbar(pm, ax=ax, label='value')
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✅ 图像已保存: {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)
        
def plot_image(data, xg, yg,
               cmap='viridis', vmin=None, vmax=None, 
               title='', xlabel = '',ylable = '',
               save_path=None, show=False):
    xc, yc = np.meshgrid(xg, yg)
    fig, ax = plt.subplots(figsize=(7, 6))
    pm = ax.pcolormesh(xc, yc, data, shading='flat', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylable)
    ax.set_title(title)
    fig.colorbar(pm, ax=ax, label='value')
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✅ 图像已保存: {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


PLOT_REQUIRED_KEYS = ("I", "Q", "U", "VZA", "VAA", "RAA", "Scattering_Angle")


def _build_level_npz_from_original(target_npz_path: str, overwrite: bool = False) -> Optional[str]:
    """
    若目标层级 NPZ 缺失（或指定覆盖），则自动从 sibling original/*.npz 重建。
    支持层级: registered / downsampled / downsampled_registered。
    """
    if os.path.exists(target_npz_path) and not overwrite:
        return target_npz_path

    target_level = os.path.basename(os.path.dirname(target_npz_path))
    if target_level not in {"registered", "downsampled", "downsampled_registered"}:
        return None

    out_root = os.path.dirname(os.path.dirname(target_npz_path))
    base = os.path.basename(target_npz_path)
    original_npz = os.path.join(out_root, "original", base)
    if not os.path.exists(original_npz):
        return None

    arr = np.load(original_npz, allow_pickle=True)
    try:
        files = set(arr.files)
        required = {"I", "Q", "U", "VZA", "VAA", "RAA", "Scattering_Angle", "x", "y"}

        def _from_sensor_dict_payload():
            if "sensor_dict" not in files:
                return None
            try:
                sensor_dict_obj = arr["sensor_dict"].item()
                context_cfg = arr["context"].item() if "context" in files else {}
                sensor_keys = list(sensor_dict_obj.keys())
                if len(sensor_keys) == 0:
                    return None
                view_hint = _infer_view_hint_from_filename(target_npz_path)
                selected_key = sensor_keys[0]
                if view_hint:
                    for k in sensor_keys:
                        if str(k).lower().startswith(str(view_hint).lower()):
                            selected_key = k
                            break

                transpose, flip_lr = _get_image_orientation_from_sensor_cfg(sen_cfg)
                sim = sensor_dict_obj.get_images(selected_key)[0]
                I0 = _apply_image_orientation(sim.I, transpose, flip_lr)
                Q0 = _apply_image_orientation(sim.Q, transpose, flip_lr) if hasattr(sim, "Q") else np.zeros_like(I0)
                U0 = _apply_image_orientation(sim.U, transpose, flip_lr) if hasattr(sim, "U") else np.zeros_like(I0)
                sensor_ds = sensor_dict_obj[selected_key]["sensor_list"][0]
                v_out_map = _apply_image_orientation(compute_vout_map_from_sensor(sensor_ds), transpose, flip_lr)
                ground = reproject_to_ground(sensor_ds, ground_z=0.0)
                x0 = _apply_image_orientation(ground.x_ground.values, transpose, flip_lr)
                y0 = _apply_image_orientation(ground.y_ground.values, transpose, flip_lr)
                vz = np.clip(v_out_map[..., 2], -1.0, 1.0)
                vza0 = np.degrees(np.arccos(vz))
                vaa0 = (np.degrees(np.arctan2(v_out_map[..., 1], v_out_map[..., 0])) + 360.0) % 360.0
                vaa0 = (vaa0 + 180.0) % 360.0
                vaa0 = ((vaa0 - _get_flight_azimuth_offset_deg_from_context(context_cfg) + 360.0) % 360.0)
                saa = (context_cfg.get("solar_azimuth", 0.0) + 360.0) % 360.0 if isinstance(context_cfg, dict) else 0.0
                sza = context_cfg.get("theta_0", np.nan) if isinstance(context_cfg, dict) else np.nan
                raa0 = ((vaa0 - saa + 180.0) % 360.0) - 180.0
                mu0 = np.cos(np.radians(sza))
                mu = np.cos(np.radians(vza0))
                cos_sca = -mu0 * mu + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu**2) * np.cos(np.radians(raa0))
                sca0 = np.degrees(np.arccos(np.clip(cos_sca, -1.0, 1.0)))
                target_shape = x0.shape
                I0 = _ensure_2d_shape(I0, target_shape)
                Q0 = _ensure_2d_shape(Q0, target_shape)
                U0 = _ensure_2d_shape(U0, target_shape)
                vza0 = _ensure_2d_shape(vza0, target_shape)
                vaa0 = _ensure_2d_shape(vaa0, target_shape)
                raa0 = _ensure_2d_shape(raa0, target_shape)
                sca0 = _ensure_2d_shape(sca0, target_shape)

                return dict(
                    I=I0, Q=Q0, U=U0,
                    VZA=vza0, VAA=vaa0, RAA=raa0, Scattering_Angle=sca0,
                    x=x0, y=y0,
                    theta0=np.full_like(I0, sza, dtype=np.float32),
                    thetav=vza0,
                    faipfai0=raa0,
                    lat=np.full_like(I0, np.nan, dtype=np.float32),
                    lon=np.full_like(I0, np.nan, dtype=np.float32),
                    elevation=np.zeros_like(I0, dtype=np.float32),
                    Land_water_mask=np.ones_like(I0, dtype=np.float32),
                    Height_AirMSPI=20000,
                )
            except Exception:
                return None

        if required.issubset(files):
            I = arr["I"]; Q = arr["Q"]; U = arr["U"]
            VZA = arr["VZA"]; VAA = arr["VAA"]; RAA = arr["RAA"]; SCA = arr["Scattering_Angle"]
            xg = arr["x"]; yg = arr["y"]
            theta0 = arr["theta0"] if "theta0" in files else np.full_like(I, np.nan, dtype=np.float32)
            thetav = arr["thetav"] if "thetav" in files else VZA
            faipfai0 = arr["faipfai0"] if "faipfai0" in files else RAA
            lat = arr["lat"] if "lat" in files else np.full_like(I, np.nan, dtype=np.float32)
            lon = arr["lon"] if "lon" in files else np.full_like(I, np.nan, dtype=np.float32)
            elevation = arr["elevation"] if "elevation" in files else np.zeros_like(I, dtype=np.float32)
            land_water_mask = arr["Land_water_mask"] if "Land_water_mask" in files else np.ones_like(I, dtype=np.float32)
            h_airmspi = arr["Height_AirMSPI"] if "Height_AirMSPI" in files else 20000
        else:
            payload0 = _from_sensor_dict_payload()
            if payload0 is None:
                return None
            I = payload0["I"]; Q = payload0["Q"]; U = payload0["U"]
            VZA = payload0["VZA"]; VAA = payload0["VAA"]; RAA = payload0["RAA"]; SCA = payload0["Scattering_Angle"]
            xg = payload0["x"]; yg = payload0["y"]
            theta0 = payload0["theta0"]; thetav = payload0["thetav"]; faipfai0 = payload0["faipfai0"]
            lat = payload0["lat"]; lon = payload0["lon"]; elevation = payload0["elevation"]
            land_water_mask = payload0["Land_water_mask"]; h_airmspi = payload0["Height_AirMSPI"]

        factor = 2
        method = "mean"
        if "dsm" in files:
            try:
                dsm_cfg = arr["dsm"].item()
                factor = int(getattr(dsm_cfg, "factor", factor))
                method = str(getattr(dsm_cfg, "method", method))
            except Exception:
                pass

        x_range = (float(np.nanmin(xg)), float(np.nanmax(xg)))
        y_range = (float(np.nanmin(yg)), float(np.nanmax(yg)))
        if "context" in files:
            try:
                ctx = arr["context"].item()
                if isinstance(ctx, dict):
                    xr = ctx.get("cloud_x_range")
                    yr = ctx.get("cloud_y_range")
                    if xr is not None and len(xr) == 2:
                        x_range = (max(float(min(xr)), x_range[0]), min(float(max(xr)), x_range[1]))
                    if yr is not None and len(yr) == 2:
                        y_range = (max(float(min(yr)), y_range[0]), min(float(max(yr)), y_range[1]))
            except Exception:
                pass

        def _crop_all():
            out = {}
            out["I"], out["x"], out["y"] = crop_by_world_box(I, xg, yg, x_range, y_range)
            out["Q"], _, _ = crop_by_world_box(Q, xg, yg, x_range, y_range)
            out["U"], _, _ = crop_by_world_box(U, xg, yg, x_range, y_range)
            out["VZA"], _, _ = crop_by_world_box(VZA, xg, yg, x_range, y_range)
            out["VAA"], _, _ = crop_by_world_box(VAA, xg, yg, x_range, y_range)
            out["RAA"], _, _ = crop_by_world_box(RAA, xg, yg, x_range, y_range)
            out["Scattering_Angle"], _, _ = crop_by_world_box(SCA, xg, yg, x_range, y_range)
            out["theta0"], _, _ = crop_by_world_box(theta0, xg, yg, x_range, y_range)
            out["thetav"], _, _ = crop_by_world_box(thetav, xg, yg, x_range, y_range)
            out["faipfai0"], _, _ = crop_by_world_box(faipfai0, xg, yg, x_range, y_range)
            out["lat"], _, _ = crop_by_world_box(lat, xg, yg, x_range, y_range)
            out["lon"], _, _ = crop_by_world_box(lon, xg, yg, x_range, y_range)
            out["elevation"], _, _ = crop_by_world_box(elevation, xg, yg, x_range, y_range)
            out["Land_water_mask"], _, _ = crop_by_world_box(land_water_mask, xg, yg, x_range, y_range)
            out["DoLP"] = np.sqrt(out["Q"] ** 2 + out["U"] ** 2) / np.maximum(out["I"], 1e-12)
            return out

        def _downsample_all(d):
            out = {}
            for k, v in d.items():
                if k == "Height_AirMSPI":
                    out[k] = v
                else:
                    out[k] = utils.downsample_block(v, factor, method)
            out["DoLP"] = np.sqrt(out["Q"] ** 2 + out["U"] ** 2) / np.maximum(out["I"], 1e-12)
            return out

        base_payload = dict(
            I=I, Q=Q, U=U,
            VZA=VZA, VAA=VAA, RAA=RAA, Scattering_Angle=SCA,
            x=xg, y=yg,
            theta0=theta0, thetav=thetav, faipfai0=faipfai0,
            lat=lat, lon=lon, elevation=elevation, Land_water_mask=land_water_mask,
            Height_AirMSPI=h_airmspi,
            DoLP=np.sqrt(Q ** 2 + U ** 2) / np.maximum(I, 1e-12),
        )

        if target_level == "registered":
            payload = _crop_all()
            payload["Height_AirMSPI"] = h_airmspi
        elif target_level == "downsampled":
            payload = _downsample_all(base_payload)
            payload["Height_AirMSPI"] = h_airmspi
        else:
            payload = _downsample_all(_crop_all())
            payload["Height_AirMSPI"] = h_airmspi

        os.makedirs(os.path.dirname(target_npz_path), exist_ok=True)
        np.savez_compressed(target_npz_path, **payload)
        return target_npz_path
    finally:
        arr.close()


def _infer_view_hint_from_filename(path: str) -> Optional[str]:
    stem = os.path.splitext(os.path.basename(path))[0]
    parts = stem.split("_", 1)
    if len(parts) == 2 and parts[1].strip():
        return parts[1].strip()
    return None


def inspect_simulation_npz(npz_path: str) -> Dict[str, Any]:
    arr = np.load(npz_path, allow_pickle=True)
    try:
        keys = set(arr.files)
        required = set(PLOT_REQUIRED_KEYS)
        info = {
            "path": npz_path,
            "keys": sorted(keys),
            "is_plot_ready": required.issubset(keys),
            "missing_plot_keys": sorted(required - keys),
            "has_sensor_dict": "sensor_dict" in keys,
            "replot_source": "direct_fields" if required.issubset(keys) else ("sensor_dict" if "sensor_dict" in keys else "insufficient"),
        }
        if "sensor_dict" in keys:
            try:
                sdict = arr["sensor_dict"].item()
                info["sensor_views"] = [str(k) for k in sdict.keys()]
            except Exception:
                info["sensor_views"] = []
        return info
    finally:
        arr.close()


def plot_simulation_results(result_path, output_dir=None, option="panel", show=False,
                            rebuild_if_missing=True, overwrite_npz=False):
    """
    读取 simulation 结果（.npz 或 .nc）并重绘图像。

    option:
      - "panel" (or "option1"): I/Q/U 一张图；角度一张图
      - "separate" (or "option2"): 每个变量单独绘图
    """
    option = str(option).lower()
    alias = {"option1": "panel", "option2": "separate"}
    option = alias.get(option, option)
    if option not in {"panel", "separate"}:
        raise ValueError("option must be 'panel'/'option1' or 'separate'/'option2'")

    requested_level = os.path.basename(os.path.dirname(result_path)).lower()
    if output_dir is None:
        output_dir = os.path.dirname(result_path)
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(result_path)[1].lower()
    cloud_box = None  # (x_range, y_range) from input txt/cloud grid if available
    if ext == ".npz":
        if rebuild_if_missing and (overwrite_npz or (not os.path.exists(result_path))):
            _build_level_npz_from_original(result_path, overwrite=overwrite_npz)

        def _try_load_from_sensor_dict(npz, view_hint=None):
            """Fallback: recover selected-view I/Q/U + geometry from metadata npz(sensor_dict)."""
            if "sensor_dict" not in npz.files:
                return None
            try:
                sensor_dict_obj = npz["sensor_dict"].item()
                sen_cfg = npz["sen"].item() if "sen" in npz.files else None
                context_cfg = npz["context"].item() if "context" in npz.files else {}

                sensor_keys = list(sensor_dict_obj.keys())
                if len(sensor_keys) == 0:
                    return None
                selected_key = sensor_keys[0]
                if view_hint:
                    for k in sensor_keys:
                        if str(k).lower() == str(view_hint).lower():
                            selected_key = k
                            break

                transpose, flip_lr = _get_image_orientation_from_sensor_cfg(sen_cfg)
                sim = sensor_dict_obj.get_images(selected_key)[0]
                I0 = _apply_image_orientation(sim.I, transpose, flip_lr)
                Q0 = _apply_image_orientation(sim.Q, transpose, flip_lr) if hasattr(sim, "Q") else np.zeros_like(I0)
                U0 = _apply_image_orientation(sim.U, transpose, flip_lr) if hasattr(sim, "U") else np.zeros_like(I0)

                sensor_ds = sensor_dict_obj[selected_key]["sensor_list"][0]
                v_out_map = _apply_image_orientation(compute_vout_map_from_sensor(sensor_ds), transpose, flip_lr)
                vz = np.clip(v_out_map[..., 2], -1.0, 1.0)
                vza0 = np.degrees(np.arccos(vz))
                vaa0 = (np.degrees(np.arctan2(v_out_map[..., 1], v_out_map[..., 0])) + 360.0) % 360.0
                vaa0 = (vaa0 + 180.0) % 360.0
                vaa0 = ((vaa0 - _get_flight_azimuth_offset_deg_from_context(context_cfg) + 360.0) % 360.0)

                saa = (context_cfg.get("solar_azimuth", 0.0) + 360.0) % 360.0
                sza = context_cfg.get("theta_0", np.nan)
                raa0 = ((vaa0 - saa + 180.0) % 360.0) - 180.0

                mu0 = np.cos(np.radians(sza))
                mu = np.cos(np.radians(vza0))
                cos_sca = -mu0 * mu + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu**2) * np.cos(np.radians(raa0))
                sca0 = np.degrees(np.arccos(np.clip(cos_sca, -1.0, 1.0)))
                ground = reproject_to_ground(sensor_ds, ground_z=0.0)
                x0 = _apply_image_orientation(ground.x_ground.values, transpose, flip_lr)
                y0 = _apply_image_orientation(ground.y_ground.values, transpose, flip_lr)
                target_shape = x0.shape
                I0 = _ensure_2d_shape(I0, target_shape)
                Q0 = _ensure_2d_shape(Q0, target_shape)
                U0 = _ensure_2d_shape(U0, target_shape)
                vza0 = _ensure_2d_shape(vza0, target_shape)
                vaa0 = _ensure_2d_shape(vaa0, target_shape)
                raa0 = _ensure_2d_shape(raa0, target_shape)
                sca0 = _ensure_2d_shape(sca0, target_shape)

                if sen_cfg is not None and hasattr(sen_cfg, "views_names") and len(sen_cfg.views_names) > 0:
                    view_name = str(sen_cfg.views_names[0])
                else:
                    view_name = str(selected_key)
                return I0, Q0, U0, vza0, vaa0, raa0, sca0, x0, y0, view_name
            except Exception:
                return None

        def _extract_plot_fields(npz):
            files = set(npz.files)
            if set(PLOT_REQUIRED_KEYS).issubset(files):
                return {
                    "I": npz["I"], "Q": npz["Q"], "U": npz["U"],
                    "VZA": npz["VZA"], "VAA": npz["VAA"], "RAA": npz["RAA"],
                    "Scattering_Angle": npz["Scattering_Angle"],
                    "x": npz["x"] if "x" in files else None,
                    "y": npz["y"] if "y" in files else None,
                }
            if {"I", "Q", "U"}.issubset(files):
                vza = npz["VZA"] if "VZA" in files else npz["thetav"]
                vaa = npz["VAA"] if "VAA" in files else np.full_like(vza, np.nan, dtype=np.float32)
                raa = npz["RAA"] if "RAA" in files else npz["faipfai0"]
                if "Scattering_Angle" in files:
                    sca = npz["Scattering_Angle"]
                elif "sca_angle" in files:
                    sca = npz["sca_angle"]
                else:
                    theta0 = npz["theta0"]
                    mu0 = np.cos(np.radians(theta0))
                    mu = np.cos(np.radians(vza))
                    cos_sca = -mu0 * mu + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu**2) * np.cos(np.radians(raa))
                    sca = np.degrees(np.arccos(np.clip(cos_sca, -1.0, 1.0)))
                return {
                    "I": npz["I"], "Q": npz["Q"], "U": npz["U"],
                    "VZA": vza, "VAA": vaa, "RAA": raa, "Scattering_Angle": sca,
                    "x": npz["x"] if "x" in files else None,
                    "y": npz["y"] if "y" in files else None,
                }
            return None

        def _open_npz_with_iqu(npz_path, allow_sensor_fallback=False):
            npz = np.load(npz_path, allow_pickle=True)
            payload = _extract_plot_fields(npz)
            if payload is not None:
                npz.close()
                return payload, npz_path
            if allow_sensor_fallback:
                fallback = _try_load_from_sensor_dict(npz, view_hint=_infer_view_hint_from_filename(npz_path))
                if fallback is not None:
                    npz.close()
                    return fallback, npz_path
            npz.close()
            return None, None


        def _infer_is_cross_track(npz_path):
            try:
                npz = np.load(npz_path, allow_pickle=True)
                mode = None
                if "context" in npz.files:
                    ctx = npz["context"].item()
                    if isinstance(ctx, dict):
                        mode = str(ctx.get("trajectory_mode", "")).lower()
                if (not mode) and ("sen" in npz.files):
                    sen_obj = npz["sen"].item()
                    mode = str(getattr(sen_obj, "trajectory_mode", "")).lower()
                npz.close()
                return mode == "cross_track"
            except Exception:
                return False

        def _try_load_cloud_box(npz_path):
            """Load cloud_x/y ranges (input txt box) from context if present."""
            try:
                npz = np.load(npz_path, allow_pickle=True)
                files = set(npz.files)
                if "context" not in files:
                    npz.close()
                    return None
                ctx = npz["context"].item()
                npz.close()
                if not isinstance(ctx, dict):
                    return None
                xr = ctx.get("cloud_x_range")
                yr = ctx.get("cloud_y_range")
                if xr is None or yr is None or len(xr) != 2 or len(yr) != 2:
                    return None
                return (tuple(sorted((float(xr[0]), float(xr[1])))),
                        tuple(sorted((float(yr[0]), float(yr[1])))))
            except Exception:
                return None

        arr, used_path = _open_npz_with_iqu(result_path, allow_sensor_fallback=False)
        cloud_box = _try_load_cloud_box(result_path)
        is_cross_track_camera = _infer_is_cross_track(result_path)
        if arr is None and requested_level == "original":
            # 对 original 路径优先尝试其自身 metadata 回退，避免误跳转到 registered 层级。
            arr, used_path = _open_npz_with_iqu(result_path, allow_sensor_fallback=True)
        if arr is None and requested_level != "original":
            # 支持传入 original 元数据 npz：自动查找同名前缀的有效结果文件
            base = os.path.basename(result_path)
            cur_dir = os.path.dirname(result_path)
            parent = os.path.dirname(cur_dir)
            search_dirs = ["downsampled_registered", "registered", "downsampled", "original"]
            for d in search_dirs:
                cand = os.path.join(parent, d, base)
                if not os.path.exists(cand):
                    continue
                arr, used_path = _open_npz_with_iqu(cand, allow_sensor_fallback=False)
                if cloud_box is None:
                    cloud_box = _try_load_cloud_box(cand)
                if arr is not None:
                    break
        if cloud_box is None:
            # 从 sibling original 获取 context（多数派生层级不含 context）
            base = os.path.basename(result_path)
            cur_dir = os.path.dirname(result_path)
            parent = os.path.dirname(cur_dir)
            original_cand = os.path.join(parent, "original", base)
            if os.path.exists(original_cand):
                cloud_box = _try_load_cloud_box(original_cand)
        if arr is None:
            # 最后兜底：仍允许直接从 metadata-only NPZ(sensor_dict) 读取，
            # 但这通常是相机投影图，不一定等同于 registered/downsampled_registered。
            arr, used_path = _open_npz_with_iqu(result_path, allow_sensor_fallback=True)

        if used_path is not None:
            is_cross_track_camera = _infer_is_cross_track(used_path)

        if arr is None:
            arr_dbg = np.load(result_path, allow_pickle=True)
            keys = list(arr_dbg.files)
            arr_dbg.close()
            raise ValueError(
                "NPZ does not contain I/Q/U arrays. "
                f"Current keys: {keys}. "
                "Please pass a data NPZ (e.g., downsampled_registered/*.npz)."
            )

        x_plot = None
        y_plot = None
        if isinstance(arr, tuple):
            I, Q, U, vza, vaa, raa, sca, x_plot, y_plot, _ = arr
        else:
            I = arr["I"]
            Q = arr["Q"]
            U = arr["U"]
            vza = arr["VZA"]
            vaa = arr["VAA"]
            raa = arr["RAA"]
            sca = arr["Scattering_Angle"]
            x_plot = arr.get("x")
            y_plot = arr.get("y")
        actual_level = os.path.basename(os.path.dirname(used_path or result_path)).lower()
        force_image_plot = (actual_level == "original")

        if (
            cloud_box is not None and
            isinstance(x_plot, np.ndarray) and isinstance(y_plot, np.ndarray) and
            x_plot.shape == I.shape and y_plot.shape == I.shape and
            (not force_image_plot)
        ):
            x_range, y_range = cloud_box
            try:
                x_base = x_plot
                y_base = y_plot
                I, x_plot, y_plot = crop_by_world_box(I, x_base, y_base, x_range, y_range)
                Q, _, _ = crop_by_world_box(Q, x_base, y_base, x_range, y_range)
                U, _, _ = crop_by_world_box(U, x_base, y_base, x_range, y_range)
                vza, _, _ = crop_by_world_box(vza, x_base, y_base, x_range, y_range)
                vaa, _, _ = crop_by_world_box(vaa, x_base, y_base, x_range, y_range)
                raa, _, _ = crop_by_world_box(raa, x_base, y_base, x_range, y_range)
                sca, _, _ = crop_by_world_box(sca, x_base, y_base, x_range, y_range)
            except Exception:
                pass
    elif ext == ".nc":
        ds = xr.open_dataset(result_path)
        view_idx = 0
        level = "downsampled_registered"
        I = ds[f"I_{level}"].isel(view=view_idx).values
        Q = ds[f"Q_{level}"].isel(view=view_idx).values
        U = ds[f"U_{level}"].isel(view=view_idx).values
        vza = ds[f"VZA_{level}"].isel(view=view_idx).values if f"VZA_{level}" in ds else ds[f"thetav_{level}"].isel(view=view_idx).values
        vaa = ds[f"VAA_{level}"].isel(view=view_idx).values if f"VAA_{level}" in ds else np.full_like(vza, np.nan, dtype=np.float32)
        raa = ds[f"RAA_{level}"].isel(view=view_idx).values if f"RAA_{level}" in ds else ds[f"faipfai0_{level}"].isel(view=view_idx).values
        if f"Scattering_Angle_{level}" in ds:
            sca = ds[f"Scattering_Angle_{level}"].isel(view=view_idx).values
        else:
            theta0 = ds[f"theta0_{level}"].isel(view=view_idx).values
            mu0 = np.cos(np.radians(theta0))
            mu = np.cos(np.radians(vza))
            cos_sca = -mu0 * mu + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu**2) * np.cos(np.radians(raa))
            sca = np.degrees(np.arccos(np.clip(cos_sca, -1.0, 1.0)))
        ds.close()
    else:
        raise ValueError("Unsupported result file. Use .npz or .nc")

    def _symmetric_limits_about_zero(data):
        max_abs = np.nanmax(np.abs(data))
        if (not np.isfinite(max_abs)) or max_abs <= 0:
            max_abs = 1.0
        return -float(max_abs), float(max_abs)

    def _plot_field(ax, data, title, cmap, vmin=None, vmax=None):
        use_ground = (
            isinstance(x_plot, np.ndarray) and isinstance(y_plot, np.ndarray) and
            x_plot.shape == data.shape and y_plot.shape == data.shape and
            (not force_image_plot)
        )
        if use_ground:
            xv, yv = centers_to_edges_2d(x_plot, y_plot)
            im = ax.pcolormesh(yv, xv, data, shading="flat", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel("x_ground [km]")
            ax.set_ylabel("y_ground [km]")
            if cloud_box is not None:
                x_range, y_range = cloud_box
                ax.set_xlim(*x_range)
                ax.set_ylim(*y_range)
        else:
            data_show = _to_aircraft_eye_view(data, flip_vertical=(not is_cross_track_camera))
            im = ax.imshow(data_show, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        return im

    if option == "panel":
        fig1, ax1 = plt.subplots(2, 2, figsize=(10, 8))
        iqud_maps = [I, Q, U, np.sqrt(Q**2 + U**2) / np.maximum(I, 1e-12)]
        iqud_names = ["I", "Q", "U", "DoLP"]
        iqud_cmaps = ["viridis", "RdBu_r", "RdBu_r", "viridis"]
        for ax, data, name, cmap in zip(ax1.flatten(), iqud_maps, iqud_names, iqud_cmaps):
            if name in {"Q", "U"}:
                vmin, vmax = _symmetric_limits_about_zero(data)
            else:
                vmin = vmax = None
            im = _plot_field(ax, data, name, cmap, vmin=vmin, vmax=vmax)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout()
        out1 = os.path.join(output_dir, "IQU_panel.png")
        plt.savefig(out1, dpi=300, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig1)

        fig2, ax2 = plt.subplots(2, 2, figsize=(10, 8))
        angle_maps = [vza, vaa, raa, sca]
        angle_titles = ["VZA", "VAA", "RAA", "Scattering Angle"]
        angle_cmaps = ["viridis", "viridis", "RdBu_r", "viridis"]
        for ax, data, title, cmap in zip(ax2.flatten(), angle_maps, angle_titles, angle_cmaps):
            im = _plot_field(ax, data, title, cmap)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout()
        out2 = os.path.join(output_dir, "Angles_panel.png")
        plt.savefig(out2, dpi=300, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig2)
    else:
        single_map = {
            "I": (I, "viridis"),
            "Q": (Q, "RdBu_r"),
            "U": (U, "RdBu_r"),
            "DoLP": (np.sqrt(Q**2 + U**2) / np.maximum(I, 1e-12), "viridis"),
            "VZA": (vza, "viridis"),
            "VAA": (vaa, "viridis"),
            "RAA": (raa, "RdBu_r"),
            "Scattering_Angle": (sca, "viridis"),
        }
        for name, (data, cmap) in single_map.items():
            fig, ax = plt.subplots(figsize=(6, 5))
            if name in {"Q", "U"}:
                vmin, vmax = _symmetric_limits_about_zero(data)
            else:
                vmin = vmax = None
            im = _plot_field(ax, data, name, cmap, vmin=vmin, vmax=vmax)
            plt.colorbar(im, ax=ax)
            plt.tight_layout()
            outp = os.path.join(output_dir, f"{name}.png")
            plt.savefig(outp, dpi=300, bbox_inches="tight")
            if show:
                plt.show()
            else:
                plt.close(fig)


def replot_simulation_results_with_config(result_path: str, cfg_path: str = "config_v5b.yaml",
                                          output_dir: Optional[str] = None, show: bool = False,
                                          rebuild_if_missing: bool = True, overwrite_npz: bool = False):
    """
    Read plotting layout from config.plot.replot_layout and call plot_simulation_results.
    """
    (cfg, out_cfg, sen_cfg, bnd_cfg, dsm_cfg, svr_cfg, plt_cfg, grd_cfg, comp_cfg, crop_cfg,
     scene_cfg, cam_cfg, solar_cfg, solver_cfg, aerosol_cfg) = utils.load_config(cfg_path)
    return plot_simulation_results(
        result_path=result_path,
        output_dir=output_dir,
        option=getattr(plt_cfg, "replot_layout", "panel"),
        show=show,
        rebuild_if_missing=rebuild_if_missing,
        overwrite_npz=overwrite_npz,
    )


def plot_all_wavelength_results(root_dir: str,
                                cfg_path: str = "config_v5b.yaml",
                                level: str = "registered",
                                view_name: Optional[str] = None,
                                show: bool = False,
                                rebuild_if_missing: bool = True,
                                overwrite_npz: bool = False) -> List[str]:
    """
    Replot every wavelength configured in YAML from existing NPZ products.

    Returns
    -------
    list[str]
        List of NPZ file paths that were replotted successfully.
    """
    (cfg, out_cfg, sen_cfg, bnd_cfg, dsm_cfg, svr_cfg, plt_cfg, grd_cfg, comp_cfg, crop_cfg,
     scene_cfg, cam_cfg, solar_cfg, solver_cfg, aerosol_cfg) = utils.load_config(cfg_path)

    level = str(level).lower()
    root = Path(root_dir)
    if not root.exists():
        raise FileNotFoundError(f"root_dir not found: {root_dir}")
    npz_dir = root / level
    if not npz_dir.exists():
        raise FileNotFoundError(f"Level folder not found: {npz_dir}")

    replotted: List[str] = []
    for w in bnd_cfg.wavelength_nm:
        if view_name:
            cand = npz_dir / f"{int(w)}nm_{view_name}.npz"
            candidates = [cand] if cand.exists() else []
        else:
            candidates = sorted(npz_dir.glob(f"{int(w)}nm_*.npz"))

        if not candidates:
            print(f"⚠️ No NPZ found for wavelength {w} nm in {npz_dir}")
            continue

        for npz_path in candidates:
            outp = npz_path.parent / f"plots_{npz_path.stem}"
            plot_simulation_results(
                result_path=str(npz_path),
                output_dir=str(outp),
                option=getattr(plt_cfg, "replot_layout", "panel"),
                show=show,
                rebuild_if_missing=rebuild_if_missing,
                overwrite_npz=overwrite_npz,
            )
            replotted.append(str(npz_path))
    return replotted

# =============================
#%% Section 4: Scene/Sensors (single band) with lat/lon ingestion
# =============================

def make_ground_grid(cfg: GroundGridConfig) -> Tuple[np.ndarray, np.ndarray]:
    x = np.linspace(cfg.x_min, cfg.x_max, cfg.nx)
    y = np.linspace(cfg.y_min, cfg.y_max, cfg.ny)
    Xg, Yg = np.meshgrid(x, y)
    return Xg, Yg

def _read_lat_lon_from_txt(cloud_scatterer: xr.Dataset):
    """
    从已加载的 cloud_scatterer（例如 wrf_to_shdom_coarse.txt 转成的 xarray.Dataset）
    中读取经纬度信息（假定包含 x, y, z, lwc, reff, lat, lon 字段）。

    返回:
        lat_2d, lon_2d : ndarray, shape (ny, nx)
    """
    # 检查关键字段是否存在
    required_fields = ["x", "y", "lat", "lon"]
    for f in required_fields:
        if f not in cloud_scatterer:
            raise ValueError(f"cloud_scatterer 缺少必要字段 '{f}'，无法提取 lat/lon。")

    lat_da = cloud_scatterer["lat"]
    lon_da = cloud_scatterer["lon"]

    # 优先处理 xarray 维度信息：若 lat/lon 含 x,y 以及额外维度(如 z)，
    # 则沿额外维度取首层，再转成 (y, x)。
    if {"x", "y"}.issubset(set(lat_da.dims)) and {"x", "y"}.issubset(set(lon_da.dims)):
        for d in list(lat_da.dims):
            if d not in {"x", "y"}:
                lat_da = lat_da.isel({d: 0})
        for d in list(lon_da.dims):
            if d not in {"x", "y"}:
                lon_da = lon_da.isel({d: 0})
        lat_2d = lat_da.transpose("y", "x").values
        lon_2d = lon_da.transpose("y", "x").values
        return np.asarray(lat_2d), np.asarray(lon_2d)

    # 回退：按平铺点表重建 2D（支持 lat/lon 是 1D 或更高维 flatten 的情况）
    x_coord = np.asarray(cloud_scatterer["x"]).ravel()
    y_coord = np.asarray(cloud_scatterer["y"]).ravel()
    lat_col = np.asarray(lat_da).ravel()
    lon_col = np.asarray(lon_da).ravel()

    # 若 lat/lon 是 x-y-z 展平，x/y 可能是坐标轴长度而非同长度点表；
    # 这里广播出同长度的点坐标。
    if lat_col.size != x_coord.size or lon_col.size != y_coord.size:
        X, Y = np.meshgrid(np.asarray(cloud_scatterer["x"]).ravel(),
                           np.asarray(cloud_scatterer["y"]).ravel(),
                           indexing="xy")
        x_pts = X.ravel()
        y_pts = Y.ravel()
        rep = max(1, int(lat_col.size // max(1, x_pts.size)))
        x_col = np.tile(x_pts, rep)[:lat_col.size]
        y_col = np.tile(y_pts, rep)[:lat_col.size]
    else:
        x_col = x_coord
        y_col = y_coord

    xs = np.unique(x_col)
    ys = np.unique(y_col)
    nx, ny = len(xs), len(ys)
    lat_2d = np.full((ny, nx), np.nan, dtype=np.float32)
    lon_2d = np.full((ny, nx), np.nan, dtype=np.float32)

    # 对每个 (x,y) 仅取首个有效点（常见于存在 z 层时）。
    order = np.lexsort((x_col, y_col))
    x_sorted = x_col[order]
    y_sorted = y_col[order]
    lat_sorted = lat_col[order]
    lon_sorted = lon_col[order]

    seen = set()
    for xx, yy, la, lo in zip(x_sorted, y_sorted, lat_sorted, lon_sorted):
        key = (float(xx), float(yy))
        if key in seen:
            continue
        seen.add(key)
        ix = np.searchsorted(xs, xx)
        iy = np.searchsorted(ys, yy)
        if 0 <= ix < nx and 0 <= iy < ny:
            lat_2d[iy, ix] = la
            lon_2d[iy, ix] = lo

    return lat_2d, lon_2d

def build_scene_and_sensors_single_band(sen: SensorConfig,
                                        wavelength_nm: int,
                                        is_polarized: bool,
                                        out_path:str,
                                        scene_cfg: SceneConfig,
                                        cam_cfg: CameraConfig,
                                        solar_cfg: SolarConfig,
                                        solver_cfg: SolverConfig,
                                        aerosol_cfg: AerosolConfig):
    t_stage = {}
    t0_all = time.perf_counter()
    input_path = scene_cfg.input_path
    t0 = time.perf_counter()

    def _load_csv_numeric_only(csv_path: str):
        try:
            return at3d.util.load_from_csv(csv_path, density=None, origin=(0.0, 0.0))
        except ValueError as e:
            if "could not convert string to float" not in str(e):
                raise
            with open(csv_path, "r", encoding="utf-8") as f:
                lines = f.readlines()
            if len(lines) < 5:
                raise
            # Keep AT3D CSV preamble exactly as expected by at3d.util.load_from_csv:
            # line1 comment, line2 nx/ny/nz, line3 dx/dy, line4 z-levels.
            header_lines = lines[:4]
            df = pd.read_csv(csv_path, skiprows=4)
            df_num = df.apply(pd.to_numeric, errors="coerce")
            # drop columns that are effectively non-numeric (e.g. surface_model='diner')
            bad_cols = [c for c in df.columns if df_num[c].isna().all() and not df[c].isna().all()]
            if bad_cols:
                print(f"⚠️ Dropping non-numeric CSV columns for at3d.load_from_csv: {bad_cols}")
            keep_cols = [c for c in df.columns if c not in bad_cols]
            with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, encoding="utf-8") as tf:
                tmp_path = tf.name
                tf.writelines(header_lines)
                df_num[keep_cols].to_csv(tf, index=False)
            try:
                return at3d.util.load_from_csv(tmp_path, density=None, origin=(0.0, 0.0))
            finally:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

    # Compatible with both retrieval-style CSV (cv) and LES-style CSV (lwc/density).
    cloud_scatterer = _load_csv_numeric_only(input_path)
    density_candidates = ["cv", "lwc", "density"]
    density_name = next((name for name in density_candidates if name in cloud_scatterer.data_vars), None)
    if density_name is None:
        raise ValueError(f"No density variable found in {input_path}. Tried {density_candidates}")
    if density_name != "density":
        cloud_scatterer = cloud_scatterer.rename_vars({density_name: "density"})
        cloud_scatterer.attrs["density_name"] = density_name

    mode_selection = str(getattr(aerosol_cfg, "mode_selection", "both")).lower()
    if mode_selection not in {"both", "mode1", "mode2"}:
        raise ValueError(f"aerosol.mode_selection must be one of ['both','mode1','mode2'], got: {mode_selection}")

    def _mode_selected(mode_id: str) -> bool:
        if mode_selection == "both":
            return True
        return mode_id == mode_selection

    # Build scalar reff/veff for AT3D from extended mode columns when available.
    mode_reff_vars = sorted([v for v in cloud_scatterer.data_vars if v.startswith('mode') and v.endswith('_reff')])
    mode_veff_vars = sorted([v for v in cloud_scatterer.data_vars if v.startswith('mode') and v.endswith('_veff')])

    if ("reff" not in cloud_scatterer) and mode_reff_vars:
        num = np.zeros(cloud_scatterer.density.shape, dtype=float)
        den = np.zeros(cloud_scatterer.density.shape, dtype=float)
        for v in mode_reff_vars:
            mode_id = v.split('_')[0]  # mode1
            if not _mode_selected(mode_id):
                continue
            frac_name = f"{mode_id}_fraction"
            frac = np.asarray(cloud_scatterer[frac_name].data, dtype=float) if frac_name in cloud_scatterer else np.ones_like(num)
            val = np.asarray(cloud_scatterer[v].data, dtype=float)
            mask = np.isfinite(frac) & np.isfinite(val) & (frac > 0)
            num[mask] += frac[mask] * val[mask]
            den[mask] += frac[mask]
        reff_eff = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        cloud_scatterer["reff"] = (['x', 'y', 'z'], reff_eff)

    if ("veff" not in cloud_scatterer) and mode_veff_vars:
        num = np.zeros(cloud_scatterer.density.shape, dtype=float)
        den = np.zeros(cloud_scatterer.density.shape, dtype=float)
        for v in mode_veff_vars:
            mode_id = v.split('_')[0]
            if not _mode_selected(mode_id):
                continue
            frac_name = f"{mode_id}_fraction"
            frac = np.asarray(cloud_scatterer[frac_name].data, dtype=float) if frac_name in cloud_scatterer else np.ones_like(num)
            val = np.asarray(cloud_scatterer[v].data, dtype=float)
            mask = np.isfinite(frac) & np.isfinite(val) & (frac > 0)
            num[mask] += frac[mask] * val[mask]
            den[mask] += frac[mask]
        veff_eff = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        cloud_scatterer["veff"] = (['x', 'y', 'z'], veff_eff)

    # Provide conservative defaults for legacy LES files.
    if "reff" not in cloud_scatterer:
        cloud_scatterer["reff"] = (['x', 'y', 'z'], np.full(cloud_scatterer.density.shape, 12.0, dtype=float))
    if "veff" not in cloud_scatterer:
        cloud_scatterer["veff"] = (['x', 'y', 'z'], np.full(cloud_scatterer.density.shape, 0.10, dtype=float))

    cloud_scatterer["density"] *= 1

    # Guard for extended CSV optional fields (e.g., mode2_mr/mode2_mi) that may be all-NaN.
    # at3d.grid.resample_onto_grid asserts each variable is not entirely NaN.
    for _name in list(cloud_scatterer.data_vars):
        _arr = np.asarray(cloud_scatterer[_name].data)
        if np.issubdtype(_arr.dtype, np.number) and np.all(np.isnan(_arr)):
            cloud_scatterer[_name] = (cloud_scatterer[_name].dims, np.zeros_like(_arr, dtype=float))
    t_stage["load_cloud_csv"] = time.perf_counter() - t0
    latlon_source = "text"
    try:
        lat_2d, lon_2d = _read_lat_lon_from_txt(cloud_scatterer)
    except Exception:
        latlon_source = "zeros_fallback_after_read_error"
        xs = cloud_scatterer.x.data.size
        ys = cloud_scatterer.y.data.size
        lat_2d = np.zeros((ys, xs), dtype=np.float32)
        lon_2d = np.zeros((ys, xs), dtype=np.float32)
    try:
        cloud_scatterer = cloud_scatterer.assign(lat=(('y','x'), lat_2d),
                                                 lon=(('y','x'), lon_2d))
    except Exception:
        pass
    
    
    # Build cloud-grid geolocation from CSV grid itself.
    # Priority: lat/lon read from file; fallback to regular-grid synthesis.
    cloud_x = np.asarray(cloud_scatterer.x.data, dtype=float)
    cloud_y = np.asarray(cloud_scatterer.y.data, dtype=float)
    ny_cloud = cloud_y.size
    nx_cloud = cloud_x.size

    latlon_is_valid = (
        lat_2d.shape == (ny_cloud, nx_cloud)
        and lon_2d.shape == (ny_cloud, nx_cloud)
        and np.isfinite(lat_2d).any()
        and np.isfinite(lon_2d).any()
        and not (np.nanmax(np.abs(lat_2d)) == 0 and np.nanmax(np.abs(lon_2d)) == 0)
    )
    if latlon_is_valid:
        xlats, xlons = lat_2d, lon_2d
    else:
        latlon_source = "synthetic_regular_grid_fallback"
        lat0 = 35.0
        lon0 = -112.0
        km_per_deg_lat = 111.32
        km_per_deg_lon = 111.32 * np.cos(np.deg2rad(lat0))
        lat_1d = lat0 + (cloud_y / km_per_deg_lat)
        lon_1d = lon0 + (cloud_x / km_per_deg_lon)
        xlons, xlats = np.meshgrid(lon_1d, lat_1d, indexing='xy')
    print(
        f"🧭 lat/lon source: {latlon_source}; "
        f"lat range=({np.nanmin(xlats):.6f}, {np.nanmax(xlats):.6f}), "
        f"lon range=({np.nanmin(xlons):.6f}, {np.nanmax(xlons):.6f})"
    )
    
    # k = 15  # vertical level

    # LWC_k = cloud_scatterer["density"].isel(z=k).values    # shape (nx, ny)
    # # lat_k = lat_2d[..., k]     # 保证 shape 相同
    # # lon_k = lon_2d[..., k]
    
    
    # def fill_nan_nearest(arr):
    #     arr2 = arr.copy()
    #     mask = np.isnan(arr2)
    #     if np.any(mask):
    #         idx = distance_transform_edt(mask,
    #                                      return_distances=False,
    #                                      return_indices=True)
    #         arr2[mask] = arr2[tuple(idx[i] for i in range(arr2.ndim))][mask]
    #     return arr2
    
    # # lon_k = fill_nan_nearest(lon_k)
    # # lat_k = fill_nan_nearest(lat_k)
    # LWC_k = fill_nan_nearest(LWC_k)
    

    
    # plt.figure(figsize=(10, 8))
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.set_title(f"Cloud field (density) at z={k}", fontsize=15)
    # ax.coastlines()
    
    # # 绘制云场
    # pcm = ax.pcolormesh(xlons, xlats, LWC_k,
    #                     cmap="viridis",
    #                     shading="auto",
    #                     transform=ccrs.PlateCarree())
    
    # plt.colorbar(pcm, ax=ax, label="Cloud Density (kg/m3)")
    # plt.show()

    atmosphere = xr.open_dataset('../data/ancillary/AFGL_summer_mid_lat.nc')
    reduced_atmosphere = atmosphere.sel({'z': atmosphere.coords['z'].data[atmosphere.coords['z'].data <= 4.0]})
    _ = at3d.grid.combine_z_coordinates([reduced_atmosphere, cloud_scatterer])
    rte_grid = at3d.grid.make_grid(cloud_scatterer.x.diff('x')[0], cloud_scatterer.x.data.size,
                                   cloud_scatterer.y.diff('y')[0], cloud_scatterer.y.data.size,
                                   np.append(0, cloud_scatterer.z.data))

    # Only resample density in z; keep geophysical/microphysical descriptor fields
    # (lat/lon, mode fractions, mode reff/veff, refractive indices, etc.) vertically invariant.
    cloud_scatterer_on_rte_grid = at3d.grid.resample_onto_grid(rte_grid, cloud_scatterer[['density']])

    reserved_grid_vars = {'density', 'delx', 'dely', 'nx', 'ny', 'nz'}
    static_like_vars = [
        name for name in cloud_scatterer.data_vars
        if name not in reserved_grid_vars
    ]
    nz_target = cloud_scatterer_on_rte_grid.sizes['z']
    for name in static_like_vars:
        da = cloud_scatterer[name]
        arr = np.asarray(da.data)
        if 'z' in da.dims:
            # Keep original vertical-invariant intent: use first z slice and broadcast.
            z_axis = da.dims.index('z')
            arr2d = np.take(arr, indices=0, axis=z_axis)
            arr3d = np.repeat(arr2d[..., None], nz_target, axis=2)
            cloud_scatterer_on_rte_grid[name] = (('x', 'y', 'z'), arr3d)
        else:
            # 2D field, broadcast over z.
            if da.dims == ('x', 'y'):
                arr3d = np.repeat(arr[..., None], nz_target, axis=2)
                cloud_scatterer_on_rte_grid[name] = (('x', 'y', 'z'), arr3d)
            elif da.dims == ('y', 'x'):
                arr_xy = np.transpose(arr, (1, 0))
                arr3d = np.repeat(arr_xy[..., None], nz_target, axis=2)
                cloud_scatterer_on_rte_grid[name] = (('x', 'y', 'z'), arr3d)
            else:
                # Generic fallback for non-standard dims.
                if da.ndim == 0:
                    arr3d = np.full((cloud_scatterer_on_rte_grid.sizes['x'],
                                     cloud_scatterer_on_rte_grid.sizes['y'],
                                     nz_target), float(arr))
                    cloud_scatterer_on_rte_grid[name] = (('x', 'y', 'z'), arr3d)
                elif {'x', 'y'}.issubset(set(da.dims)):
                    arr_xy = np.asarray(da.transpose('x', 'y').data)
                    arr3d = np.repeat(arr_xy[..., None], nz_target, axis=2)
                    cloud_scatterer_on_rte_grid[name] = (('x', 'y', 'z'), arr3d)
                else:
                    # Skip unsupported auxiliary variable dimensions.
                    print(f"⚠️ Skip static var '{name}' with dims={da.dims} (unsupported for z-broadcast)")

    # Ensure grid spacing variables remain scalar for at3d.checks.check_grid.
    cloud_scatterer_on_rte_grid['delx'] = xr.DataArray(float(np.asarray(rte_grid.delx).reshape(-1)[0]))
    cloud_scatterer_on_rte_grid['dely'] = xr.DataArray(float(np.asarray(rte_grid.dely).reshape(-1)[0]))

    # Optional single-mode selection: scale total density by selected mode fraction.
    if mode_selection in {"mode1", "mode2"}:
        frac_name = f"{mode_selection}_fraction"
        if frac_name in cloud_scatterer_on_rte_grid:
            frac = np.asarray(cloud_scatterer_on_rte_grid[frac_name].data, dtype=float)
            frac = np.clip(np.nan_to_num(frac, nan=0.0, posinf=0.0, neginf=0.0), 0.0, 1.0)
            dens = np.asarray(cloud_scatterer_on_rte_grid["density"].data, dtype=float)
            cloud_scatterer_on_rte_grid["density"] = (
                cloud_scatterer_on_rte_grid["density"].dims,
                dens * frac
            )
        else:
            print(f"⚠️ mode_selection={mode_selection} but '{frac_name}' not found; density unchanged.")

    # Recompute scalar reff/veff from selected mode(s) on RTE grid when mode fields exist.
    mode_reff_vars_rte_all = sorted([v for v in cloud_scatterer_on_rte_grid.data_vars if v.startswith('mode') and v.endswith('_reff')])
    mode_veff_vars_rte_all = sorted([v for v in cloud_scatterer_on_rte_grid.data_vars if v.startswith('mode') and v.endswith('_veff')])
    if mode_reff_vars_rte_all:
        num = np.zeros(cloud_scatterer_on_rte_grid["density"].shape, dtype=float)
        den = np.zeros(cloud_scatterer_on_rte_grid["density"].shape, dtype=float)
        for v in mode_reff_vars_rte_all:
            mode_id = v.split('_')[0]
            if not _mode_selected(mode_id):
                continue
            frac_name = f"{mode_id}_fraction"
            frac = np.asarray(cloud_scatterer_on_rte_grid[frac_name].data, dtype=float) if frac_name in cloud_scatterer_on_rte_grid else np.ones_like(num)
            arr = np.asarray(cloud_scatterer_on_rte_grid[v].data, dtype=float)
            m = np.isfinite(arr) & np.isfinite(frac) & (frac > 0)
            num[m] += frac[m] * arr[m]
            den[m] += frac[m]
        reff_eff = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        cloud_scatterer_on_rte_grid["reff"] = (('x', 'y', 'z'), reff_eff)
    if mode_veff_vars_rte_all:
        num = np.zeros(cloud_scatterer_on_rte_grid["density"].shape, dtype=float)
        den = np.zeros(cloud_scatterer_on_rte_grid["density"].shape, dtype=float)
        for v in mode_veff_vars_rte_all:
            mode_id = v.split('_')[0]
            if not _mode_selected(mode_id):
                continue
            frac_name = f"{mode_id}_fraction"
            frac = np.asarray(cloud_scatterer_on_rte_grid[frac_name].data, dtype=float) if frac_name in cloud_scatterer_on_rte_grid else np.ones_like(num)
            arr = np.asarray(cloud_scatterer_on_rte_grid[v].data, dtype=float)
            m = np.isfinite(arr) & np.isfinite(frac) & (frac > 0)
            num[m] += frac[m] * arr[m]
            den[m] += frac[m]
        veff_eff = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        cloud_scatterer_on_rte_grid["veff"] = (('x', 'y', 'z'), veff_eff)

    size_distribution_function = at3d.size_distribution.gamma
    size_distribution_function = at3d.size_distribution.lognormal
    # cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.reff.dims,
    #                                        np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=0.1))
    sensor_dict = at3d.containers.SensorsDict()
    center = np.array(scene_cfg.lookat_center_km, dtype=float)
    t0 = time.perf_counter()
    # AT3D/SHDOM world coordinates follow x=North, y=East, z=Up (NEU).
    # This axis ordering is left-handed (x × y = -z), so avoid accidental x/y swaps.
    center_NEU = center.copy()
    mode_lc = str(sen.trajectory_mode).lower()
    view_names = list(sen.views_names)
    if mode_lc == "cross_track":
        required = [
            sen.cross_track_x1, sen.cross_track_y1, sen.cross_track_z1,
            sen.cross_track_x2, sen.cross_track_y2, sen.cross_track_z2
        ]
        if any(v is None for v in required):
            raise ValueError(
                "cross_track mode requires trajectory.cross_track_x1/y1/z1 and x2/y2/z2 in config."
            )
        stokes = ['I', 'Q', 'U'] if is_polarized else ['I']
        base_name = str(view_names[0]) if len(view_names) > 0 else "cross_track"
        cross_sensor, scan_positions, scan_angles, scan_pitch_deg = cross_track_scan_projection(
            wavelength=wavelength_nm/1000,
            stokes=stokes,
            x1=sen.cross_track_x1, y1=sen.cross_track_y1, z1=sen.cross_track_z1,
            x2=sen.cross_track_x2, y2=sen.cross_track_y2, z2=sen.cross_track_z2,
            spacing=sen.cross_track_spacing,
            scan1_deg=sen.cross_track_scan1_deg,
            scan2_deg=sen.cross_track_scan2_deg,
            delscan_deg=sen.cross_track_delscan_deg,
            pitch_start_deg=sen.cross_track_pitch_start_deg,
            pitch_end_deg=sen.cross_track_pitch_end_deg,
            pitch_list_deg=sen.cross_track_pitch_list_deg,
        )
        key = f"{base_name}_{int(wavelength_nm)}nm"
        sensor_dict.add_sensor(key, cross_sensor)
        view_names = [base_name]
        sen.views_names = view_names
        position_vectors = np.array([scan_positions[0]], dtype=float)
        lookat_vectors = np.full((1, 3), np.nan, dtype=float)
        up_vectors = np.full((1, 3), np.nan, dtype=float)
    else:
        position_vectors, up_vectors = calculate_sensor_trajectory(
            sensor_zenith_list=sen.views_zenith_deg,
            sensor_azimuth_list=sen.views_azimuth_deg,
            look_at_point=center,
            sensor_altitude=sen.altitude_km,
            trajectory_mode=sen.trajectory_mode,
            fallback_heading_deg=sen.fallback_heading_deg,
            manual_flight_azimuth_deg=sen.manual_flight_azimuth_deg,
            camera_relative_roll_deg=sen.camera_relative_roll_deg,
            camera_align_with_flight_heading=sen.camera_align_with_flight_heading,
        )
        position_vectors, up_vectors = calculate_sensor_trajectory_from_aircraft(
            look_at_point=center_NEU,
            sensor_altitude=sen.altitude_km,
            heading_angle_deg=sen.aircraft_heading_deg,
            pitch_angle_deg=sen.aircraft_pitch_deg,
            roll_angle_deg=sen.aircraft_roll_deg,
            camera_pitch_relative_deg=sen.camera_pitch_relative_deg,
            camera_roll_relative_deg=sen.camera_roll_relative_deg,
            n_views=1
            )
        lookat_vectors = [center for _ in view_names]
        scan_pitch_deg = None
    if mode_lc != "cross_track":
        for name, pos, look, up in zip(view_names, position_vectors, lookat_vectors, up_vectors):
            stokes = ['I', 'Q', 'U'] if is_polarized else ['I']
            if sen.type == "perspective_projection":
                sensor = at3d.sensor.perspective_projection(
                    wavelength=wavelength_nm/1000,
                    fov=sen.fov_deg,
                    x_resolution=int(cam_cfg.x_resolution),
                    y_resolution=int(cam_cfg.y_resolution),
                    position_vector=pos,
                    lookat_vector=look,
                    up_vector=up,
                    stokes=stokes
                )
            else:
                raise NotImplementedError(
                    f"sensor.type={sen.type} is not implemented in v5b; "
                    "currently only perspective_projection is supported."
                )
            key = f"{name}_{int(wavelength_nm)}nm"
            sensor_dict.add_sensor(key, sensor)
    t_stage["build_sensors"] = time.perf_counter() - t0
    wavelengths = sensor_dict.get_unique_solvers()
    mie_mono_tables = OrderedDict()
    t0 = time.perf_counter()
    for wavelength in wavelengths:
        cache_key = (
            float(wavelength),
            float(aerosol_cfg.refractive_index_real),
            float(aerosol_cfg.refractive_index_imag),
        )
        if cache_key in _MIE_TABLE_CACHE:
            mie_mono_tables[wavelength] = _MIE_TABLE_CACHE[cache_key]
        else:
            mie_mono_tables[wavelength] = at3d.mie.get_mono_table(
                particle_type='Aerosol',
                wavelength_band=(wavelength, wavelength),
                max_integration_radius=65.0,
                minimum_effective_radius=0.1,
                refractive_index=aerosol_cfg.refractive_index_real - aerosol_cfg.refractive_index_imag*1j,
                relative_dir='../mie_tables',
                verbose=False
            )
            _MIE_TABLE_CACHE[cache_key] = mie_mono_tables[wavelength]
    t_stage["mie_table"] = time.perf_counter() - t0
        
    rho_p_gcm3 = 1.6                 # 你选定的颗粒密度
    ho_p_kgm3 = rho_p_gcm3 * 1000.0 # kg/m^3    
    particle_density_gm3 = rho_p_gcm3 * 1e6  # g/cm^3 → g/m^3
    
    # === sanitize microphysics (avoid NaN/Inf/outliers causing interpolation range errors) ===
    density_data = np.asarray(cloud_scatterer_on_rte_grid.density.data, dtype=float)
    if float(aerosol_cfg.density_floor) > 0:
        density_data[density_data < float(aerosol_cfg.density_floor)] = 0.0
        cloud_scatterer_on_rte_grid['density'] = (cloud_scatterer_on_rte_grid.density.dims, density_data)
    clear_air = (~np.isfinite(density_data)) | (density_data <= 0)

    reff_data = np.asarray(cloud_scatterer_on_rte_grid.reff.data, dtype=float)
    veff_data = np.asarray(cloud_scatterer_on_rte_grid.veff.data, dtype=float)

    bad_reff = (~np.isfinite(reff_data)) | (reff_data <= 0)
    bad_veff = (~np.isfinite(veff_data)) | (veff_data <= 0)
    if np.any(bad_reff):
        reff_data[bad_reff] = float(aerosol_cfg.reff_default)
    if np.any(bad_veff):
        veff_data[bad_veff] = float(aerosol_cfg.veff_default)

    # Avoid interpolation artifacts in clear-air cells.
    reff_data[clear_air] = float(aerosol_cfg.reff_default)
    veff_data[clear_air] = float(aerosol_cfg.veff_default)

    # Conservative clipping for aerosol retrieval products.
    reff_data = np.clip(reff_data, float(aerosol_cfg.reff_clip_min), float(aerosol_cfg.reff_clip_max))
    veff_data = np.clip(veff_data, float(aerosol_cfg.veff_clip_min), float(aerosol_cfg.veff_clip_max))

    cloud_scatterer_on_rte_grid['reff'] = (cloud_scatterer_on_rte_grid.reff.dims, reff_data)
    cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.veff.dims, veff_data)
    print(
        f"🧪 microphysics clip range: reff[{aerosol_cfg.reff_clip_min}, {aerosol_cfg.reff_clip_max}], "
        f"veff[{aerosol_cfg.veff_clip_min}, {aerosol_cfg.veff_clip_max}]"
    )

    # === 从 scatterer 中动态获取范围（finite-only, multi-mode aware）===
    reff_candidates = [reff_data]
    veff_candidates = [veff_data]
    mode_reff_vars_rte = sorted([v for v in cloud_scatterer_on_rte_grid.data_vars if v.startswith('mode') and v.endswith('_reff')])
    mode_veff_vars_rte = sorted([v for v in cloud_scatterer_on_rte_grid.data_vars if v.startswith('mode') and v.endswith('_veff')])

    for v in mode_reff_vars_rte:
        mode_id = v.split('_')[0]
        if not _mode_selected(mode_id):
            continue
        frac_name = f"{mode_id}_fraction"
        frac = np.asarray(cloud_scatterer_on_rte_grid[frac_name].data, dtype=float) if frac_name in cloud_scatterer_on_rte_grid else np.ones_like(reff_data)
        arr = np.asarray(cloud_scatterer_on_rte_grid[v].data, dtype=float)
        m = np.isfinite(arr) & np.isfinite(frac) & (frac > 0)
        if np.any(m):
            reff_candidates.append(arr[m])
    for v in mode_veff_vars_rte:
        mode_id = v.split('_')[0]
        if not _mode_selected(mode_id):
            continue
        frac_name = f"{mode_id}_fraction"
        frac = np.asarray(cloud_scatterer_on_rte_grid[frac_name].data, dtype=float) if frac_name in cloud_scatterer_on_rte_grid else np.ones_like(veff_data)
        arr = np.asarray(cloud_scatterer_on_rte_grid[v].data, dtype=float)
        m = np.isfinite(arr) & np.isfinite(frac) & (frac > 0)
        if np.any(m):
            veff_candidates.append(arr[m])

    reff_all = np.concatenate([np.ravel(a) for a in reff_candidates])
    veff_all = np.concatenate([np.ravel(a) for a in veff_candidates])
    reff_finite = reff_all[np.isfinite(reff_all)]
    veff_finite = veff_all[np.isfinite(veff_all)]
    reff_min = float(np.min(reff_finite))
    reff_max = float(np.max(reff_finite))
    veff_min = float(np.min(veff_finite))
    veff_max = float(np.max(veff_finite))

    # === 留 margin（非常重要）===
    margin_reff = 0.1   # μm
    margin_veff = 0.01

    reff_grid = np.linspace(
        max(0.05, reff_min - margin_reff),
        reff_max + margin_reff,
        40
    )

    veff_grid = np.linspace(
        max(0.01, veff_min - margin_veff),
        veff_max + margin_veff,
        10    # ⚠️ 如果 veff 几乎不变，直接用 1
    )

    # Ensure microphysics coordinates lie inside interpolation grid (numerical safety).
    cloud_scatterer_on_rte_grid['reff'] = (
        cloud_scatterer_on_rte_grid.reff.dims,
        np.clip(cloud_scatterer_on_rte_grid.reff.data, reff_grid.min(), reff_grid.max())
    )
    cloud_scatterer_on_rte_grid['veff'] = (
        cloud_scatterer_on_rte_grid.veff.dims,
        np.clip(cloud_scatterer_on_rte_grid.veff.data, veff_grid.min(), veff_grid.max())
    )


    optical_property_generator = at3d.medium.OpticalPropertyGenerator(
        'aerosol', 
        mie_mono_tables,
        size_distribution_function,
        reff=reff_grid,
        veff=veff_grid,
        density_normalization='density',  # ← 关键,
        particle_density=aerosol_cfg.particle_density
    )
    t0 = time.perf_counter()
    optical_properties = optical_property_generator(cloud_scatterer_on_rte_grid)
    t_stage["optical_property_generator"] = time.perf_counter() - t0
    
    rho = cloud_scatterer_on_rte_grid.density.data  # g/m^3
    z = rte_grid.z.data                              # km (你已经确认)
    dz_m = np.diff(z) * 1000.0                       # m
    
    Mcol_py = (rho[..., :-1] * dz_m).sum(axis=-1)    # g/m^2
    print("Mcol_py max:", float(Mcol_py.max()))
    rayleigh_scattering = at3d.rayleigh.to_grid(wavelengths, atmosphere, rte_grid)
    theta_0 = 180.0 - float(solar_cfg.sza_deg)
    solarmu = np.cos(np.deg2rad(theta_0))
    # TODO: 
    solar_azimuth = float(solar_cfg.saa_deg)
    solvers_dict = at3d.containers.SolversDict()
    config = at3d.configuration.get_config()
    # config["x_boundary_condition"]='periodic'
    # config["y_boundary_condition"]='periodic'
    
    # config["ip_flag"] = 3          # ⭐ 关键：Independent Pixel (1D)
    config["x_boundary_condition"] = solver_cfg.x_boundary_condition
    config["y_boundary_condition"] = solver_cfg.y_boundary_condition
    config["num_mu_bins"] = int(solver_cfg.num_mu_bins)
    config["num_phi_bins"] = int(solver_cfg.num_phi_bins)
    config["split_accuracy"] = float(solver_cfg.split_accuracy)
    config["deltam"] = bool(solver_cfg.deltam)
    if solver_cfg.adapt_grid_factor is not None:
        config["adapt_grid_factor"] = float(solver_cfg.adapt_grid_factor)
    if solver_cfg.cell_to_point_ratio is not None:
        config["cell_to_point_ratio"] = float(solver_cfg.cell_to_point_ratio)
    if solver_cfg.max_total_mb is not None:
        config["max_total_mb"] = float(solver_cfg.max_total_mb)
    
    
    w = wavelength_nm/1000
    medium = {'aerosol': optical_properties[w], 'rayleigh': rayleigh_scattering[w]}
    
    AOD = utils.compute_column_aod(medium)
    # SSA = utils.compute_column_ssa(medium)
    SSA = {}
    SSA['aerosol'], tau_cloud = utils.compute_cloud_column_ssa(medium)


    utils.plot_field(AOD['aerosol'], 
                     f"AOD total @ {w*1000} nm", 
                     os.path.join(out_path, f"AOD_total_{w*1000}nm.png"),
                     vmin=0,
                     vmax=4)
    # utils.plot_field(SSA['total'], f"SSA total @ {w} nm", os.path.join(out_path, f"SSA_total_{w}nm.png"))
    
    utils.plot_field(
        SSA['aerosol'],
        f"Aerosol SSA @ {w*1000} nm",
        os.path.join(out_path, f"SSA_aerosol_{w*1000}nm.png"),
        vmin=0.90,
        vmax=0.96
    )
    
    # plt.figure(figsize=(4, 3))
    # im = plt.imshow(
    #     SSA['aerosol'].T,
    #     origin="lower",
    #     cmap="viridis",
    #     vmin=0.74,
    #     vmax=1
    # )
    # plt.colorbar(im, label=f"Aerosol SSA @ {w} nm")
    # plt.title(f"Aerosol SSA @ {w} nm")
    # plt.tight_layout()

    # plt.savefig(os.path.join(out_path, f"SSA_aerosol_{w}nm.png"), dpi=300)
    # plt.close()
    
    # utils.plot_field(AOD['cloud'], f"AOD cloud @ {w} nm", os.path.join(out_path, f"AOD_cloud_{w}nm.png"))
    # utils.plot_field(AOD['rayleigh'], f"AOD Rayleigh @ {w} nm", os.path.join(out_path, f"AOD_rayleigh_{w}nm.png"))
    
    tab = mie_mono_tables[w]
    print("wavelength_center in table =", tab.attrs["wavelength_center"])
    print("wavelength_band in table   =", tab.attrs["wavelength_band"])
    
    np.savez_compressed(f"medium_{w}nm.npz",
                        medium=medium)
    
    num_stokes = 3 if is_polarized else 1

    def _surface_param_2d(name: str, default_value: float) -> np.ndarray:
        if name in cloud_scatterer_on_rte_grid:
            arr = np.asarray(cloud_scatterer_on_rte_grid[name].data, dtype=float)
            if arr.ndim == 3:
                arr2d = arr[:, :, 0]
            elif arr.ndim == 2:
                arr2d = arr
            else:
                arr2d = np.full((rte_grid.x.size, rte_grid.y.size), float(default_value), dtype=float)
        else:
            arr2d = np.full((rte_grid.x.size, rte_grid.y.size), float(default_value), dtype=float)
        arr2d = np.nan_to_num(arr2d, nan=float(default_value), posinf=float(default_value), neginf=float(default_value))
        return arr2d

    enable_brdf = bool(getattr(scene_cfg, "enable_brdf", False))
    enable_bpdf = bool(getattr(scene_cfg, "enable_bpdf", False))
    brdf_model = str(getattr(scene_cfg, "brdf_model", "diner")).lower()
    delx_km = float(np.asarray(rte_grid.delx).reshape(-1)[0])
    dely_km = float(np.asarray(rte_grid.dely).reshape(-1)[0])

    if enable_brdf and brdf_model == "diner":
        A = _surface_param_2d("a0_surf", float(getattr(scene_cfg, "brdf_default_a", 0.0)))
        K = _surface_param_2d("k0_surf", float(getattr(scene_cfg, "brdf_default_k", 1.0)))
        B = _surface_param_2d("b0_surf", float(getattr(scene_cfg, "brdf_default_b", 0.0)))
        if enable_bpdf:
            E = _surface_param_2d("e0_surface", float(getattr(scene_cfg, "bpdf_default_e", 0.0)))
            wind_vv = _surface_param_2d("wind_vv", float(getattr(scene_cfg, "bpdf_default_wind_vv", 5.0)))
            wind_wd = _surface_param_2d("wind_wd", float(getattr(scene_cfg, "bpdf_default_wind_wd", 0.0)))
            if np.any(np.abs(wind_wd) > 1e-9):
                print("⚠️ wind_wd is currently not used by at3d.surface.diner and will be ignored.")
            zeta = E
            sigma = np.sqrt(np.maximum(0.003 + 0.00512 * np.maximum(wind_vv, 0.0), 1e-6) / 2.0)
        else:
            zeta = np.zeros_like(A)
            sigma = np.zeros_like(A)
        surface_ds = at3d.surface.diner(A=A, K=K, B=B, ZETA=zeta, SIGMA=sigma, delx=delx_km, dely=dely_km)
    else:
        if enable_brdf and brdf_model != "diner":
            print(f"⚠️ Unsupported brdf_model='{brdf_model}', fallback to lambertian.")
        lambert_alb = float(getattr(scene_cfg, "lambertian_albedo", 0.0))
        surface_ds = at3d.surface.lambertian(lambert_alb)

    t0 = time.perf_counter()
    solvers_dict.add_solver(
        w,
        at3d.solver.RTE(
            numerical_params=config,
            surface=surface_ds,
            source=at3d.source.solar(w, solarmu, solar_azimuth),
            medium=medium,
            num_stokes=num_stokes
        )
    )
    t_stage["build_rte_solver"] = time.perf_counter() - t0
    context = dict(center=center,
                   view_names=view_names,
                   position_vectors=position_vectors,
                   lookat_vectors=lookat_vectors,
                   up_vectors=up_vectors,
                   theta_0=theta_0,
                   solar_azimuth=solar_azimuth,
                   trajectory_mode=sen.trajectory_mode,
                   manual_flight_azimuth_deg=sen.manual_flight_azimuth_deg,
                   fallback_heading_deg=sen.fallback_heading_deg,
                   camera_relative_roll_deg=sen.camera_relative_roll_deg,
                   camera_align_with_flight_heading=sen.camera_align_with_flight_heading,
                   camera_image_transpose=sen.camera_image_transpose,
                   camera_image_flip_lr=sen.camera_image_flip_lr,
                   heading_angle_deg=sen.aircraft_heading_deg,
                   cross_track_nbytes=sen.cross_track_nbytes,
                   cross_track_scale=sen.cross_track_scale,
                   cross_track_x1=sen.cross_track_x1,
                   cross_track_y1=sen.cross_track_y1,
                   cross_track_z1=sen.cross_track_z1,
                   cross_track_x2=sen.cross_track_x2,
                   cross_track_y2=sen.cross_track_y2,
                   cross_track_z2=sen.cross_track_z2,
                   cross_track_spacing=sen.cross_track_spacing,
                   cross_track_scan1_deg=sen.cross_track_scan1_deg,
                   cross_track_scan2_deg=sen.cross_track_scan2_deg,
                   cross_track_delscan_deg=sen.cross_track_delscan_deg,
                   cross_track_pitch_start_deg=sen.cross_track_pitch_start_deg,
                   cross_track_pitch_end_deg=sen.cross_track_pitch_end_deg,
                   cross_track_pitch_list_deg=sen.cross_track_pitch_list_deg,
                   cross_track_scan_positions=(scan_positions.tolist() if mode_lc == "cross_track" else None),
                   cross_track_scan_angles_deg=(scan_angles.tolist() if mode_lc == "cross_track" else None),
                   cross_track_scan_pitch_deg=(scan_pitch_deg.tolist() if mode_lc == "cross_track" else None),
                   lat=lat_2d,
                   lon=lon_2d,
                   latlon_source=latlon_source,
                   cloud_x=cloud_x,
                   cloud_y=cloud_y,
                   cloud_x_range=(float(np.nanmin(cloud_x)), float(np.nanmax(cloud_x))),
                   cloud_y_range=(float(np.nanmin(cloud_y)), float(np.nanmax(cloud_y))))
    t_stage["total_scene_build"] = time.perf_counter() - t0_all
    print("⏱️ Scene build timing (s): " + ", ".join([f"{k}={v:.2f}" for k, v in t_stage.items()]))
    return sensor_dict, solvers_dict, context, medium, AOD, SSA, xlats, xlons


# =============================
#%% Section 5: Build versions for ONE band (with camera & terrain PNGs + metadata)
# =============================

def build_versions_single_band(sensor_dict,
                               xlats, xlons,
                               wavelength_nm: int,
                               is_polarized: bool,
                               sen: SensorConfig,
                               grd: GroundGridConfig,
                               dsm: DownsampleConfig,
                               plot_cfg: PlotConfig,
                               out_cfg: OutputConfig,
                               out_subdirs: Dict[str, str],
                               crop_cfg: GroundCropConfig,
                               context: Dict):
    Xg, Yg = make_ground_grid(grd)
    
    
    
    
    first_key = list(sensor_dict.keys())[0]
    first_image = sensor_dict.get_images(first_key)[0]
    cam_ny, cam_nx = first_image.I.T.shape
    V = len(sen.views_names)
    I_orig = np.full((V, cam_ny, cam_nx), np.nan, dtype=np.float32)
    Q_orig = np.full_like(I_orig, np.nan); U_orig = np.full_like(I_orig, np.nan)
    DoLP_orig = np.full_like(I_orig, np.nan)
    I_reg  = np.full((V, grd.ny, grd.nx), np.nan, dtype=np.float32)
    Q_reg  = np.full_like(I_reg,  np.nan); U_reg  = np.full_like(I_reg,  np.nan)
    DoLP_reg = np.full_like(I_reg, np.nan)
    ny_ds = (cam_ny // dsm.factor)
    nx_ds = (cam_nx // dsm.factor)
    I_ds = np.full((V, ny_ds, nx_ds), np.nan, dtype=np.float32)
    Q_ds = np.full_like(I_ds, np.nan); U_ds = np.full_like(I_ds, np.nan)
    DoLP_ds = np.full_like(I_ds, np.nan)
    ny_gds = None
    nx_gds = None
    I_reg_ds = None
    Q_reg_ds = None
    U_reg_ds = None
    DoLP_reg_ds = None
    VZA_reg_ds = None
    VAA_reg_ds = None
    RAA_reg_ds = None
    SCA_reg_ds = None
    thetav_o = np.zeros((V, cam_ny, cam_nx), dtype=np.float32)
    theta0_o = np.zeros_like(thetav_o); faipfai0_o = np.zeros_like(thetav_o)
    thetav_r = np.zeros((V, grd.ny, grd.nx), dtype=np.float32)
    theta0_r = np.zeros_like(thetav_r); faipfai0_r = np.zeros_like(thetav_r)
    # Camera lat/lon fill (mean of context field)
    lat_cam = np.full((cam_ny, cam_nx), np.nan, dtype=np.float32)
    lon_cam = np.full((cam_ny, cam_nx), np.nan, dtype=np.float32)

    # Downsampled metadata
    lat_cam_ds = utils.downsample_block(lat_cam, dsm.factor, "mean")
    lon_cam_ds = utils.downsample_block(lon_cam, dsm.factor, "mean")
    # lat_grd_ds = downsample_block(lat_grd, dsm.factor, "mean")
    # lon_grd_ds = downsample_block(lon_grd, dsm.factor, "mean")
    
    
    
    def _unit(v, eps=1e-12):
        n = np.linalg.norm(v)
        return v / max(n, eps)

    def _fit_to_shape(a: np.ndarray, target_shape: Tuple[int, int], fill_value=np.nan):
        """Crop/pad 2D array to target_shape (top-left aligned)."""
        ty, tx = target_shape
        out = np.full((ty, tx), fill_value, dtype=a.dtype if np.issubdtype(a.dtype, np.floating) else np.float32)
        sy = min(ty, a.shape[0])
        sx = min(tx, a.shape[1])
        out[:sy, :sx] = a[:sy, :sx]
        return out

    def _stokes_rotate_QU(Q, U, chi):
        """Rotate (Q,U) by angle chi (radians) using standard Stokes rotation."""
        c2 = np.cos(2.0 * chi)
        s2 = np.sin(2.0 * chi)
        Qp = Q * c2 + U * s2
        Up = -Q * s2 + U * c2
        return Qp, Up
    
    def rotate_to_scattering_plane(Q, U, pos, lookat, up_vec, theta0_deg, solar_azimuth_deg):
        """
        Rotate Q/U from the current sensor reference (defined by up_vec)
        to the scattering-plane reference.
    
        Coordinate convention assumed consistent with your code:
          x: North (+), y: East (+), z: Up (+)
          azimuth angles measured in x-y plane from +x toward +y (math convention)
        """
    
        # --- outgoing/view direction v (direction of propagation toward sensor) ---
        # Use from lookat point -> sensor position (same geometry used to build sensor)
        v = _unit(np.array(pos) - np.array(lookat))
    
        # --- incoming solar direction s (direction of propagation from sun -> scene) ---
        th = np.deg2rad(theta0_deg)
        ph = np.deg2rad(solar_azimuth_deg)
        s = _unit(np.array([np.sin(th) * np.cos(ph),
                            np.sin(th) * np.sin(ph),
                            -np.cos(th)]))
    
        # --- current reference axis e_cam: up projected onto plane ⟂ v ---
        up = _unit(np.array(up_vec))
        e_cam = up - np.dot(up, v) * v
        e_cam = _unit(e_cam)
    
        # --- scattering-plane reference axis e_scat: in scattering plane, ⟂ v ---
        n = np.cross(s, v)
        n_norm = np.linalg.norm(n)
        if n_norm < 1e-10:
            # sun and view are (nearly) colinear: scattering plane ill-defined
            return Q, U, 0.0
        n = n / n_norm
        e_scat = _unit(np.cross(n, v))
    
        # --- signed rotation angle chi taking e_cam -> e_scat about axis v ---
        cos_chi = np.clip(np.dot(e_cam, e_scat), -1.0, 1.0)
        sin_chi = np.dot(v, np.cross(e_cam, e_scat))
        chi = np.arctan2(sin_chi, cos_chi)
    
        Qs, Us = _stokes_rotate_QU(Q, U, chi)
        return Qs, Us, chi
    
    def _unit_vec(a, eps=1e-12):
        n = np.linalg.norm(a, axis=-1, keepdims=True)
        return a / np.maximum(n, eps)

    def rotate_qu_pixelwise_to_scattering_plane(Q, U, v_out_map, up_vec_world, theta0_deg_shdom, solar_azimuth_deg):
        """
        Q,U: (ny,nx)
        v_out_map: (ny,nx,3) scene->sensor unit vectors
        up_vec_world: (3,) the SAME up_vector used to build the sensor (world coords)
        theta0_deg_shdom: SHDOM theta relative to +z (e.g., 145 deg), so solar goes downward
        """
        ny, nx = Q.shape
    
        # solar incoming direction s (sun -> scene), world coords
        th = np.deg2rad(theta0_deg_shdom)
        ph = np.deg2rad(solar_azimuth_deg)
        s = np.array([np.sin(th)*np.cos(ph),
                      np.sin(th)*np.sin(ph),
                      -np.cos(th)], dtype=float)
        s = s / np.linalg.norm(s)
    
        v = _unit_vec(v_out_map)
    
        # e_cam: project up onto plane ⟂ v  (pixelwise)
        up = np.array(up_vec_world, dtype=float)
        up = up / np.linalg.norm(up)
        dot_up_v = (v[...,0]*up[0] + v[...,1]*up[1] + v[...,2]*up[2])[...,None]
        e_cam = up[None,None,:] - dot_up_v * v
        e_cam = _unit_vec(e_cam)
    
        # scattering-plane axis e_scat
        n = np.cross(s[None,None,:], v)          # (ny,nx,3)
        n_norm = np.linalg.norm(n, axis=-1)
        valid = n_norm > 1e-10
        n = _unit_vec(n)
    
        e_scat = np.cross(n, v)
        e_scat = _unit_vec(e_scat)
    
        # signed chi: e_cam -> e_scat about axis v
        cos_chi = np.clip(np.sum(e_cam * e_scat, axis=-1), -1.0, 1.0)
        sin_chi = np.sum(v * np.cross(e_cam, e_scat), axis=-1)
        chi = np.arctan2(sin_chi, cos_chi)
    
        # rotate
        c2 = np.cos(2.0*chi)
        s2 = np.sin(2.0*chi)
        Qs = Q*c2 + U*s2
        Us = -Q*s2 + U*c2
    
        # where scattering plane ill-defined, keep original
        Qs = np.where(valid, Qs, Q)
        Us = np.where(valid, Us, U)
    
        return Qs, Us, chi
    
    for iv, (name, pos, look, up) in enumerate(zip(
            sen.views_names, context["position_vectors"], context["lookat_vectors"], context["up_vectors"])):
        key = f"{name}_{int(wavelength_nm)}nm"
        sim = sensor_dict.get_images(key)[0]
        sensor = sensor_dict[key]
        sensor_ds = sensor['sensor_list'][0]
        ground = reproject_to_ground(sensor_ds, ground_z=0.0)
        xg = ground.x_ground.values
        yg = ground.y_ground.values
        # Use cloud grid coordinates from CSV/scatterer (km), not hard-coded extents.
        wrf_x = np.asarray(context.get("cloud_x"), dtype=float)
        wrf_y = np.asarray(context.get("cloud_y"), dtype=float)
        lat_img, lon_img = assign_latlon_from_grid(xg, yg, wrf_x, wrf_y, xlats, xlons)


        
        
        transpose, flip_lr = _get_image_orientation_from_sensor_cfg(sen)
        I = _apply_image_orientation(sim.I, transpose, flip_lr)
        Q = _apply_image_orientation(sim.Q, transpose, flip_lr) if is_polarized and hasattr(sim, "Q") else np.zeros_like(I)
        U = _apply_image_orientation(sim.U, transpose, flip_lr) if is_polarized and hasattr(sim, "U") else np.zeros_like(I)
        xg = _apply_image_orientation(xg, transpose, flip_lr)
        yg = _apply_image_orientation(yg, transpose, flip_lr)
        lat_img = _apply_image_orientation(lat_img, transpose, flip_lr)
        lon_img = _apply_image_orientation(lon_img, transpose, flip_lr)
        
        
        theta0 = context.get("theta_0")
        phi0   = context.get("solar_azimuth", 0.0)
        # Q, U, chi = rotate_to_scattering_plane(
        #     Q, U,
        #     pos=pos,
        #     lookat=look,
        #     up_vec=up,
        #     theta0_deg=theta0,
        #     solar_azimuth_deg=phi0
        # )
        
        v_out_map = compute_vout_map_from_sensor(sensor_ds)
        
        v_out_map = _apply_image_orientation(v_out_map, transpose, flip_lr)
        target_shape = v_out_map.shape[:2]
        I = _ensure_2d_shape(I, target_shape)
        Q = _ensure_2d_shape(Q, target_shape)
        U = _ensure_2d_shape(U, target_shape)
        xg = _ensure_2d_shape(xg, target_shape)
        yg = _ensure_2d_shape(yg, target_shape)
        lat_img = _ensure_2d_shape(lat_img, target_shape)
        lon_img = _ensure_2d_shape(lon_img, target_shape)
        
        vza_map, vaa_map, raa_map, sca_angle = _compute_angle_maps_from_sensor(
            sensor_ds=sensor_ds,
            solar_azimuth_deg=context.get("solar_azimuth", 0.0),
            solar_zenith_deg=context.get("theta_0", np.nan),
            heading_angle_deg=_get_flight_azimuth_offset_deg_from_context(context),
            apply_heading_offset=bool(getattr(sen, "apply_flight_azimuth_offset_to_vaa", False)),
            transpose=transpose,
            flip_lr=flip_lr,
        )

        is_cross_track_mode = str(sensor_ds.attrs.get("projection", "")).lower() == "crosstrackscan"

        if out_cfg.save_png and (out_cfg.plot_mode != "skip"):
            angle_products = {
                "VZA": (vza_map, plot_cfg.colormap),
                "VAA": (vaa_map, plot_cfg.colormap),
                "RAA": (raa_map, "RdBu_r"),
                "Scattering_Angle": (sca_angle, plot_cfg.colormap),
            }
            for angle_name, (angle_map, angle_cmap) in angle_products.items():
                angle_png = os.path.join(
                    out_subdirs["original"],
                    f"{angle_name}_{int(wavelength_nm)}nm_{name}_camera.png"
                )
                if out_cfg.plot_mode == "overwrite" or not os.path.exists(angle_png):
                    plot_image(
                        _to_aircraft_eye_view(angle_map, flip_vertical=not is_cross_track_mode),
                        np.arange(1, angle_map.shape[1] + 2),
                        np.arange(1, angle_map.shape[0] + 2),
                        title=f"{name} {angle_name} {int(wavelength_nm)} nm",
                        cmap=angle_cmap,
                        save_path=angle_png,
                        show=False
                    )

        # Q, U, chi = rotate_qu_pixelwise_to_scattering_plane(
        #     Q, U,
        #     v_out_map=v_out_map,
        #     up_vec_world=up,               # 这里用你构造 sensor 时的那个 up_vector
        #     theta0_deg_shdom=theta0,
        #     solar_azimuth_deg=phi0
        # )
        
        DoLP = np.sqrt(Q**2 + U**2) / np.maximum(I, 1e-12)
        I_orig[iv] = I; Q_orig[iv] = Q; U_orig[iv] = U; DoLP_orig[iv] = DoLP
        
        
        mu0 = np.abs(np.cos(np.deg2rad(theta0)))
        
        I_brf = I * np.pi / mu0
        Q_brf = Q * np.pi / mu0
        U_brf = U * np.pi / mu0
        if out_cfg.save_png and (out_cfg.plot_mode != "skip"):
            cam_png_I = os.path.join(out_subdirs["original"], f"I_{int(wavelength_nm)}nm_{name}_camera.png")
            cam_png_Q = os.path.join(out_subdirs["original"], f"Q_{int(wavelength_nm)}nm_{name}_camera.png")
            cam_png_U = os.path.join(out_subdirs["original"], f"U_{int(wavelength_nm)}nm_{name}_camera.png")
            cam_png_DoLP = os.path.join(out_subdirs["original"], f"DoLP_{int(wavelength_nm)}nm_{name}_camera.png")
            if out_cfg.plot_mode == "overwrite" or not os.path.exists(cam_png_I):

                plot_image(_to_aircraft_eye_view(I_brf, flip_vertical=not is_cross_track_mode), np.arange(1,I.shape[1]+2), np.arange(1,I.shape[0]+2), 
                           title=f"{name} I {int(wavelength_nm)} nm",
                           cmap=plot_cfg.colormap, 
                           save_path=cam_png_I, show=False)
     
                        
                if is_polarized:
                    q_lim = np.nanmax(np.abs(Q_brf))
                    u_lim = np.nanmax(np.abs(U_brf))
                    if (not np.isfinite(q_lim)) or q_lim <= 0:
                        q_lim = 1.0
                    if (not np.isfinite(u_lim)) or u_lim <= 0:
                        u_lim = 1.0
                    plot_image(_to_aircraft_eye_view(Q_brf, flip_vertical=not is_cross_track_mode), np.arange(1,I.shape[1]+2), np.arange(1,I.shape[0]+2),
                               title=f"{name} Q {int(wavelength_nm)} nm",
                               cmap="RdBu_r",
                               vmin=-q_lim, vmax=q_lim,
                               save_path=cam_png_Q, show=False)
                    plot_image(_to_aircraft_eye_view(U_brf, flip_vertical=not is_cross_track_mode), np.arange(1,I.shape[1]+2), np.arange(1,I.shape[0]+2),
                               title=f"{name} U {int(wavelength_nm)} nm",
                               cmap="RdBu_r",
                               vmin=-u_lim, vmax=u_lim,
                               save_path=cam_png_U, show=False)
                    plot_image(_to_aircraft_eye_view(DoLP, flip_vertical=not is_cross_track_mode), np.arange(1,I.shape[1]+2), np.arange(1,I.shape[0]+2),
                               title=f"{name} DoLP {int(wavelength_nm)} nm",
                               cmap=plot_cfg.colormap,
                               save_path=cam_png_DoLP, show=False)
                
                
                # fig, ax = plt.subplots(figsize=(7, 6))
                # ax.imshow(I, origin='lower', cmap=plot_cfg.colormap)
                # fig.colorbar(pm, ax=ax, label='value')
                # ax.set_title(f"{name} {int(wavelength_nm)} nm (Original View)")
                # plt.tight_layout()
                # plt.savefig(cam_png, dpi=out_cfg.png_dpi, bbox_inches='tight')
                # plt.close(fig)
        # I_g = project_to_ground_lookat(I, pos, look, up, sen.fov_deg, Xg, Yg)
        # Q_g = project_to_ground_lookat(Q, pos, look, up, sen.fov_deg, Xg, Yg) if is_polarized else np.zeros_like(I_g)
        # U_g = project_to_ground_lookat(U, pos, look, up, sen.fov_deg, Xg, Yg) if is_polarized else np.zeros_like(I_g)
        # DoLP_g = np.sqrt(Q_g**2 + U_g**2) / np.maximum(I_g, 1e-12)
        # I_reg[iv] = I_g; Q_reg[iv] = Q_g; U_reg[iv] = U_g; DoLP_reg[iv] = DoLP_g
        
        I_brf_d   = utils.downsample_block(I_brf, dsm.factor, dsm.method)
        Q_brf_d   = utils.downsample_block(Q_brf, dsm.factor, dsm.method)
        U_brf_d   = utils.downsample_block(U_brf, dsm.factor, dsm.method)
        vza_d     = utils.downsample_block(vza_map, dsm.factor, dsm.method)
        vaa_d     = utils.downsample_block(vaa_map, dsm.factor, dsm.method)
        raa_d     = utils.downsample_block(raa_map, dsm.factor, dsm.method)
        sca_d     = utils.downsample_block(sca_angle, dsm.factor, dsm.method)
        DoLP_brf_d = np.sqrt(Q_brf_d**2 + U_brf_d**2) / np.maximum(I_brf_d, 1e-12)
        I_ds[iv] = I_brf_d; Q_ds[iv] = Q_brf_d; U_ds[iv] = U_brf_d; DoLP_ds[iv] = DoLP_brf_d
        
        
        # Ground crop window should match original CSV cloud x/y coverage.
        x_range = tuple(context.get("cloud_x_range", (grd.x_min, grd.x_max)))
        y_range = tuple(context.get("cloud_y_range", (grd.y_min, grd.y_max)))

        # Intersect with actual projected coordinates to avoid empty windows.
        x_min, x_max = np.nanmin(xg), np.nanmax(xg)
        y_min, y_max = np.nanmin(yg), np.nanmax(yg)
        x_range = (max(min(x_range), x_min), min(max(x_range), x_max))
        y_range = (max(min(y_range), y_min), min(max(y_range), y_max))

        I_brf_g, xg_g, yg_g, = crop_by_world_box(I_brf, xg, yg, x_range, y_range)
        
        Q_brf_g, _, _, = crop_by_world_box(Q_brf, xg, yg, x_range, y_range)
        U_brf_g, _, _, = crop_by_world_box(U_brf, xg, yg, x_range, y_range)
        vza_g, _, _, = crop_by_world_box(vza_map, xg, yg, x_range, y_range)
        vaa_g, _, _, = crop_by_world_box(vaa_map, xg, yg, x_range, y_range)
        raa_g, _, _, = crop_by_world_box(raa_map, xg, yg, x_range, y_range)
        sca_g, _, _, = crop_by_world_box(sca_angle, xg, yg, x_range, y_range)
        DoLP_brf_g = np.sqrt(Q_brf_g**2 + U_brf_g**2) / np.maximum(I_brf_g, 1e-12)
        lat_img_g, _, _, = crop_by_world_box(lat_img, xg, yg, x_range, y_range)
        lon_img_g, _, _, = crop_by_world_box(lon_img, xg, yg, x_range, y_range)
        
        
        
        if out_cfg.save_png and (out_cfg.plot_mode != "skip"):
            terr_I_png = os.path.join(out_subdirs["registered"], f"I_{int(wavelength_nm)}nm_{name}_terrain.png")
            terr_Q_png = os.path.join(out_subdirs["registered"], f"Q_{int(wavelength_nm)}nm_{name}_terrain.png")
            terr_U_png = os.path.join(out_subdirs["registered"], f"U_{int(wavelength_nm)}nm_{name}_terrain.png")
            if out_cfg.plot_mode == "overwrite" or not os.path.exists(terr_I_png):
                try:
                    I_g_c, Xg_c, Yg_c = crop_by_world_box(I_brf_g, xg_g, yg_g, crop_cfg.x_range, crop_cfg.y_range)
                except Exception:
                    I_g_c, Xg_c, Yg_c = I_brf_g, xg_g, yg_g
                plot_on_ground(I_g_c, Xg_c, Yg_c,
                               title=f"{name} I {int(wavelength_nm)} nm",
                               cmap=plot_cfg.colormap, save_path=terr_I_png, show=False)
                if is_polarized:
                    Q_g_c, _, _ = crop_by_world_box(Q_brf_g, xg_g, yg_g, crop_cfg.x_range, crop_cfg.y_range)
                    U_g_c, _, _ = crop_by_world_box(U_brf_g, xg_g, yg_g, crop_cfg.x_range, crop_cfg.y_range)
                    qg_lim = np.nanmax(np.abs(Q_g_c))
                    ug_lim = np.nanmax(np.abs(U_g_c))
                    if (not np.isfinite(qg_lim)) or qg_lim <= 0:
                        qg_lim = 1.0
                    if (not np.isfinite(ug_lim)) or ug_lim <= 0:
                        ug_lim = 1.0
                    plot_on_ground(Q_g_c, Xg_c, Yg_c,
                                   title=f"{name} Q {int(wavelength_nm)} nm",
                                   cmap="RdBu_r", vmin=-qg_lim, vmax=qg_lim,
                                   save_path=terr_Q_png, show=False)
                    plot_on_ground(U_g_c, Xg_c, Yg_c,
                                   title=f"{name} U {int(wavelength_nm)} nm",
                                   cmap="RdBu_r", vmin=-ug_lim, vmax=ug_lim,
                                   save_path=terr_U_png, show=False)

                angle_products_reg = {
                    "VZA": (vza_g, plot_cfg.colormap),
                    "VAA": (vaa_g, plot_cfg.colormap),
                    "RAA": (raa_g, "RdBu_r"),
                    "Scattering_Angle": (sca_g, plot_cfg.colormap),
                }
                for angle_name, (angle_map_reg, angle_cmap_reg) in angle_products_reg.items():
                    angle_reg_png = os.path.join(
                        out_subdirs["registered"],
                        f"{angle_name}_{int(wavelength_nm)}nm_{name}_terrain.png"
                    )
                    try:
                        angle_reg_c, Xg_angle_c, Yg_angle_c = crop_by_world_box(
                            angle_map_reg, xg_g, yg_g, crop_cfg.x_range, crop_cfg.y_range
                        )
                    except Exception:
                        angle_reg_c, Xg_angle_c, Yg_angle_c = angle_map_reg, xg_g, yg_g
                    plot_on_ground(
                        angle_reg_c, Xg_angle_c, Yg_angle_c,
                        title=f"{name} {angle_name} {int(wavelength_nm)} nm",
                        cmap=angle_cmap_reg,
                        save_path=angle_reg_png,
                        show=False
                    )
        
        I_gd = utils.downsample_block(I_brf_g, dsm.factor, dsm.method)
        Q_gd = utils.downsample_block(Q_brf_g, dsm.factor, dsm.method)
        U_gd = utils.downsample_block(U_brf_g, dsm.factor, dsm.method)
        DoLP_gd = np.sqrt(Q_gd**2 + U_gd**2) / np.maximum(I_gd, 1e-12)
        lat_img_gd = utils.downsample_block(lat_img_g, dsm.factor, dsm.method)
        lon_img_gd = utils.downsample_block(lon_img_g, dsm.factor, dsm.method)
        vza_gd = utils.downsample_block(vza_g, dsm.factor, dsm.method)
        vaa_gd = utils.downsample_block(vaa_g, dsm.factor, dsm.method)
        raa_gd = utils.downsample_block(raa_g, dsm.factor, dsm.method)
        sca_gd = utils.downsample_block(sca_g, dsm.factor, dsm.method)
        if I_reg_ds is None:
            ny_gds, nx_gds = I_gd.shape
            I_reg_ds = np.full((V, ny_gds, nx_gds), np.nan, dtype=np.float32)
            Q_reg_ds = np.full_like(I_reg_ds, np.nan); U_reg_ds = np.full_like(I_reg_ds, np.nan)
            DoLP_reg_ds = np.full_like(I_reg_ds, np.nan)
            VZA_reg_ds = np.full_like(I_reg_ds, np.nan)
            VAA_reg_ds = np.full_like(I_reg_ds, np.nan)
            RAA_reg_ds = np.full_like(I_reg_ds, np.nan)
            SCA_reg_ds = np.full_like(I_reg_ds, np.nan)
        if I_gd.shape != (ny_gds, nx_gds):
            I_gd = _fit_to_shape(I_gd, (ny_gds, nx_gds))
            Q_gd = _fit_to_shape(Q_gd, (ny_gds, nx_gds))
            U_gd = _fit_to_shape(U_gd, (ny_gds, nx_gds))
            DoLP_gd = _fit_to_shape(DoLP_gd, (ny_gds, nx_gds))
            vza_gd = _fit_to_shape(vza_gd, (ny_gds, nx_gds))
            vaa_gd = _fit_to_shape(vaa_gd, (ny_gds, nx_gds))
            raa_gd = _fit_to_shape(raa_gd, (ny_gds, nx_gds))
            sca_gd = _fit_to_shape(sca_gd, (ny_gds, nx_gds))
        I_reg_ds[iv] = I_gd
        Q_reg_ds[iv] = Q_gd
        U_reg_ds[iv] = U_gd
        DoLP_reg_ds[iv] = DoLP_gd
        VZA_reg_ds[iv] = vza_gd
        VAA_reg_ds[iv] = vaa_gd
        RAA_reg_ds[iv] = raa_gd
        SCA_reg_ds[iv] = sca_gd
        
        xg_gd = utils.downsample_block(xg_g, dsm.factor, dsm.method)
        
        yg_gd = utils.downsample_block(yg_g, dsm.factor, dsm.method)
        
        # Q_gd   = utils.downsample_block(Q_g, dsm.factor, dsm.method)
        # U_gd   = utils.downsample_block(U_g, dsm.factor, dsm.method)
        # DoLP_gd = np.sqrt(Q_gd**2 + U_gd**2) / np.maximum(I_gd, 1e-12)
        # I_reg_ds[iv] = I_gd; Q_reg_ds[iv] = Q_gd; U_reg_ds[iv] = U_gd; DoLP_reg_ds[iv] = DoLP_gd
        thetav_o[iv][:] = sen.views_zenith_deg[iv]
        theta0_o[iv][:] = context.get("theta_0", 180 - 35)
        faipfai0_o[iv][:] = sen.views_azimuth_deg[iv] - (context.get("solar_azimuth", 325.0 - 360))
        thetav_r[iv][:] = sen.views_zenith_deg[iv]
        theta0_r[iv][:] = context.get("theta_0", 180 - 35)
        faipfai0_r[iv][:] = sen.views_azimuth_deg[iv] - (context.get("solar_azimuth", 325.0 - 360))
        elevation_o = np.full_like(I_brf, 0, dtype=np.float32)
        Land_water_mask_o = np.full_like(I_brf, 1, dtype=np.float32)
        elevation_r = np.full_like(I_brf_g, 0, dtype=np.float32)
        Land_water_mask_r = np.full_like(I_brf_g, 1, dtype=np.float32)
        elevation_ds = np.full_like(I_brf_d, 0, dtype=np.float32)
        Land_water_mask_ds = np.full_like(I_brf_d, 1, dtype=np.float32)
        elevation_gds = np.full_like(I_gd, 0, dtype=np.float32)
        Land_water_mask_gds = np.full_like(I_gd, 1, dtype=np.float32)
        prefix = f"{int(wavelength_nm)}nm_{name}"
        
        # # === 静态地理网格 (degree per km 简化近似) ===
        # # 左下角为 (35N, -112E)
        # # 1°纬度约 111 km；经度按纬度 cos(φ) 缩放
        # lat0, lon0 = 35.0, -112.0
        # km_per_deg_lat = 111.0
        # km_per_deg_lon = 111.0 * np.cos(np.deg2rad(lat0))
        
        # # 实际覆盖范围 (x方向 15 km, y方向 6 km)
        # x_len_km = grd.x_max - grd.x_min
        # y_len_km = grd.y_max - grd.y_min
        
        # dlat = y_len_km / km_per_deg_lat   # ≈ 0.054°
        # dlon = x_len_km / km_per_deg_lon   # ≈ 0.136°
        
        # # 构建地理坐标网格，与 downsampled_registered 图像尺寸相同
        # lat_grd_ds = np.linspace(lat0, lat0 + dlat, ny_gds).reshape(-1, 1) * np.ones((1, nx_gds))
        # lon_grd_ds = np.ones((ny_gds, 1)) * np.linspace(lon0, lon0 + dlon, nx_gds)
        
        np.savez_compressed(
            os.path.join(out_subdirs["original"], f"{prefix}.npz"),
            I=I_brf, Q=Q_brf, U=U_brf, DoLP=DoLP,
            VZA=vza_map, VAA=vaa_map, RAA=raa_map, Scattering_Angle=sca_angle,
            x=xg, y=yg,
            theta0=theta0_o[iv], thetav=thetav_o[iv], faipfai0=faipfai0_o[iv],
            lat=lat_img, lon=lon_img, elevation=elevation_o,
            Height_AirMSPI=20000, Land_water_mask=Land_water_mask_o,
            sensor_dict=sensor_dict, wavelength_nm=wavelength_nm, is_polarized=is_polarized, sen=sen,
            grd=grd, dsm=dsm, plot_cfg=plot_cfg, out_cfg=out_cfg, out_subdirs=out_subdirs,
            crop_cfg=crop_cfg, context=context
        )
        np.savez_compressed(
            os.path.join(out_subdirs["registered"], f"{prefix}.npz"),
            I=I_brf_g, Q=Q_brf_g, U=U_brf_g, DoLP=DoLP_brf_g,
            VZA=vza_g, VAA=vaa_g, RAA=raa_g, Scattering_Angle=sca_g,
            x=xg_g, y=yg_g,
            theta0=theta0_r[iv][:I_brf_g.shape[0], :I_brf_g.shape[1]],
            thetav=thetav_r[iv][:I_brf_g.shape[0], :I_brf_g.shape[1]],
            faipfai0=faipfai0_r[iv][:I_brf_g.shape[0], :I_brf_g.shape[1]],
            lat=lat_img_g, lon=lon_img_g, elevation=elevation_r[:I_brf_g.shape[0], :I_brf_g.shape[1]],
            Height_AirMSPI=20000, Land_water_mask=Land_water_mask_r[:I_brf_g.shape[0], :I_brf_g.shape[1]]
        )
        np.savez_compressed(
            os.path.join(out_subdirs["downsampled"], f"{prefix}.npz"),
            I=I_brf_d, Q=Q_brf_d, U=U_brf_d, DoLP=DoLP_brf_d,
            VZA=vza_d, VAA=vaa_d, RAA=raa_d, Scattering_Angle=sca_d,
            theta0=utils.downsample_block(theta0_o[iv], dsm.factor),
            thetav=utils.downsample_block(thetav_o[iv], dsm.factor),
            faipfai0=utils.downsample_block(faipfai0_o[iv], dsm.factor),
            lat=lat_cam_ds, lon=lon_cam_ds, elevation=elevation_ds,
            Height_AirMSPI=20000, Land_water_mask=Land_water_mask_ds
        )
        np.savez_compressed(os.path.join(out_subdirs["downsampled_registered"], f"{prefix}.npz"),
                            I=I_gd, Q=Q_gd, U=U_gd, DoLP=DoLP_gd,
                            VZA=vza_gd, VAA=vaa_gd, RAA=raa_gd, Scattering_Angle=sca_gd,
                            x=xg_gd,y=yg_gd,
                            theta0=utils.downsample_block(theta0_r[iv], dsm.factor),
                            thetav=utils.downsample_block(thetav_r[iv], dsm.factor),
                            faipfai0=utils.downsample_block(faipfai0_r[iv], dsm.factor),
                            lat=lat_img_gd, lon=lon_img_gd, elevation=elevation_gds,
                            Height_AirMSPI=20000, Land_water_mask=Land_water_mask_gds)
    if I_reg_ds is None:
        ny_gds = max(1, grd.ny // dsm.factor)
        nx_gds = max(1, grd.nx // dsm.factor)
        I_reg_ds = np.full((V, ny_gds, nx_gds), np.nan, dtype=np.float32)
        Q_reg_ds = np.full_like(I_reg_ds, np.nan); U_reg_ds = np.full_like(I_reg_ds, np.nan)
        DoLP_reg_ds = np.full_like(I_reg_ds, np.nan)
        VZA_reg_ds = np.full_like(I_reg_ds, np.nan)
        VAA_reg_ds = np.full_like(I_reg_ds, np.nan)
        RAA_reg_ds = np.full_like(I_reg_ds, np.nan)
        SCA_reg_ds = np.full_like(I_reg_ds, np.nan)

    ds = xr.Dataset(
        data_vars=dict(
            I_original=(["view", "y", "x"], I_orig),
            Q_original=(["view", "y", "x"], Q_orig),
            U_original=(["view", "y", "x"], U_orig),
            DoLP_original=(["view", "y", "x"], DoLP_orig),
            thetav_original=(["view", "y", "x"], thetav_o),
            theta0_original=(["view", "y", "x"], theta0_o),
            faipfai0_original=(["view", "y", "x"], faipfai0_o),
            I_registered=(["view", "y_g", "x_g"], I_reg),
            Q_registered=(["view", "y_g", "x_g"], Q_reg),
            U_registered=(["view", "y_g", "x_g"], U_reg),
            DoLP_registered=(["view", "y_g", "x_g"], DoLP_reg),
            thetav_registered=(["view", "y_g", "x_g"], thetav_r),
            theta0_registered=(["view", "y_g", "x_g"], theta0_r),
            faipfai0_registered=(["view", "y_g", "x_g"], faipfai0_r),
            I_downsampled=(["view", "y_ds", "x_ds"], I_ds),
            Q_downsampled=(["view", "y_ds", "x_ds"], Q_ds),
            U_downsampled=(["view", "y_ds", "x_ds"], U_ds),
            DoLP_downsampled=(["view", "y_ds", "x_ds"], DoLP_ds),
            I_downsampled_registered=(["view", "y_gds", "x_gds"], I_reg_ds),
            Q_downsampled_registered=(["view", "y_gds", "x_gds"], Q_reg_ds),
            U_downsampled_registered=(["view", "y_gds", "x_gds"], U_reg_ds),
            DoLP_downsampled_registered=(["view", "y_gds", "x_gds"], DoLP_reg_ds),
            VZA_downsampled_registered=(["view", "y_gds", "x_gds"], VZA_reg_ds),
            VAA_downsampled_registered=(["view", "y_gds", "x_gds"], VAA_reg_ds),
            RAA_downsampled_registered=(["view", "y_gds", "x_gds"], RAA_reg_ds),
            Scattering_Angle_downsampled_registered=(["view", "y_gds", "x_gds"], SCA_reg_ds),
        ),
        coords=dict(
            view=("view", np.array(sen.views_names, dtype=str)),
            y=np.arange(cam_ny), x=np.arange(cam_nx),
            y_g=np.linspace(grd.y_min, grd.y_max, grd.ny),
            x_g=np.linspace(grd.x_min, grd.x_max, grd.nx),
            y_ds=np.arange(ny_ds), x_ds=np.arange(nx_ds),
            y_gds=np.arange(ny_gds), x_gds=np.arange(nx_gds),
        ),
        attrs=dict(description=f"AT3D simulated AirMSPI dataset for {int(wavelength_nm)} nm")
    )
    return ds


def _record_experiment_snapshot(cfg_path: str, cfg: Dict[str, Any], out_root: str) -> None:
    """
    Save per-run config snapshot and append a lightweight experiment log entry.
    """
    try:
        os.makedirs(out_root, exist_ok=True)
        ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")

        cfg_text = ""
        if cfg_path and os.path.exists(cfg_path):
            with open(cfg_path, "r", encoding="utf-8") as f:
                cfg_text = f.read()
        else:
            cfg_text = yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True)

        cfg_hash = hashlib.md5(cfg_text.encode("utf-8")).hexdigest()
        snapshot_path = os.path.join(out_root, f"config_snapshot_{ts}.yaml")
        with open(snapshot_path, "w", encoding="utf-8") as f:
            f.write(cfg_text)

        git_commit = "unknown"
        try:
            git_commit = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=os.path.dirname(__file__),
                text=True
            ).strip()
        except Exception:
            pass

        sensor_cfg = cfg.get("sensor", {}) if isinstance(cfg, dict) else {}
        traj_cfg = sensor_cfg.get("trajectory", {}) if isinstance(sensor_cfg, dict) else {}

        record = {
            "timestamp_utc": ts,
            "root_dir": out_root,
            "cfg_path": cfg_path,
            "cfg_hash_md5": cfg_hash,
            "git_commit": git_commit,
            "trajectory_mode": traj_cfg.get("mode"),
            "manual_flight_azimuth_deg": traj_cfg.get("manual_flight_azimuth_deg"),
            "fallback_heading_deg": traj_cfg.get("fallback_heading_deg"),
            "views_zenith_deg": sensor_cfg.get("views", {}).get("zenith_deg") if isinstance(sensor_cfg.get("views", {}), dict) else None,
            "views_azimuth_deg": sensor_cfg.get("views", {}).get("azimuth_deg") if isinstance(sensor_cfg.get("views", {}), dict) else None,
            "snapshot_file": os.path.basename(snapshot_path),
        }

        log_path = os.path.join(out_root, "experiment_log.jsonl")
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")
    except Exception as e:
        print(f"⚠️ Failed to record experiment snapshot: {e}")


# =============================
#%% Section 6: Main (per-wavelength)
# =============================

def main(cfg_path: str = "config_v5b.yaml", only_band: Optional[int] = None,
         n_jobs: Optional[int] = None, overwrite: bool = False, clean_after_band: bool = True,
         precheck_only: bool = False):
    (cfg, out_cfg, sen_cfg, bnd_cfg, dsm_cfg, svr_cfg, plt_cfg, grd_cfg, comp_cfg, crop_cfg,
     scene_cfg, cam_cfg, solar_cfg, solver_cfg, aerosol_cfg) = utils.load_config(cfg_path)
    _record_experiment_snapshot(cfg_path, cfg, out_cfg.root_dir)
    subdirs = utils.ensure_dirs(out_cfg.root_dir, svr_cfg.versions)
    if n_jobs is None:
        n_jobs = comp_cfg.n_jobs or max(1, os.cpu_count() // 2)
    print(f"💻 Using n_jobs = {n_jobs}")
    bands = list(zip(bnd_cfg.wavelength_nm, bnd_cfg.is_polarized))
    if only_band is not None:
        bands = [(w, pol) for (w, pol) in bands if int(w) == int(only_band)]
        if not bands:
            raise ValueError(f"Requested band {only_band} not found in config.")
            
    spectral_AOD_all = {}
    spectral_SSA_all = {}

    def _save_precheck_panels(sensor_dict, context, AOD, SSA, wavelength_nm):
        preview_dir = os.path.join(out_cfg.root_dir, "precheck")
        os.makedirs(preview_dir, exist_ok=True)

        def _aod_grid_to_surface_via_camera(sensor_ds):
            """
            1) Interpolate grid AOD to camera pixels (camera image domain).
            2) Reproject that camera image back to surface using the same per-pixel ground intersections.
            This guarantees pixel-by-pixel correspondence with angle maps.
            """
            cloud_x = np.asarray(context.get("cloud_x"), dtype=float).ravel()
            cloud_y = np.asarray(context.get("cloud_y"), dtype=float).ravel()
            aod_grid = np.asarray(AOD.get("aerosol"), dtype=float)
            if aod_grid.ndim != 2 or cloud_x.size < 2 or cloud_y.size < 2:
                return None, None, None

            # Normalize AOD grid shape to (ny, nx) for RegularGridInterpolator((y, x), values)
            if aod_grid.shape == (cloud_x.size, cloud_y.size):
                aod_grid = aod_grid.T
            elif aod_grid.shape != (cloud_y.size, cloud_x.size):
                return None, None, None

            ground = reproject_to_ground(sensor_ds, ground_z=0.0)
            transpose, flip_lr = _get_image_orientation_from_sensor_cfg(sen_cfg)
            xg = _apply_image_orientation(ground.x_ground.values, transpose, flip_lr)
            yg = _apply_image_orientation(ground.y_ground.values, transpose, flip_lr)

            interp = RegularGridInterpolator(
                (cloud_y, cloud_x),
                aod_grid,
                bounds_error=False,
                fill_value=np.nan
            )

            # Step-1: grid -> camera pixels
            pts = np.stack([yg.ravel(), xg.ravel()], axis=-1)
            aod_camera = interp(pts).reshape(xg.shape)
            # Step-2: camera image -> surface (same pixel rays / same xg, yg)
            aod_surface = aod_camera.copy()
            return xg, yg, aod_surface

        for view_name in sen_cfg.views_names:
            key = f"{view_name}_{int(wavelength_nm)}nm"
            if key not in sensor_dict:
                continue
            sensor_ds = sensor_dict[key]["sensor_list"][0]
            ground = reproject_to_ground(sensor_ds, ground_z=0.0)
            transpose, flip_lr = _get_image_orientation_from_sensor_cfg(sen_cfg)
            xg = _apply_image_orientation(ground.x_ground.values, transpose, flip_lr)
            yg = _apply_image_orientation(ground.y_ground.values, transpose, flip_lr)
            vza_map, vaa_map, raa_map, sca_angle = _compute_angle_maps_from_sensor(
                sensor_ds=sensor_ds,
                solar_azimuth_deg=context.get("solar_azimuth", 0.0),
                solar_zenith_deg=context.get("theta_0", np.nan),
                heading_angle_deg=_get_flight_azimuth_offset_deg_from_context(sen_cfg),
                apply_heading_offset=bool(getattr(sen_cfg, "apply_flight_azimuth_offset_to_vaa", False)),
                transpose=transpose,
                flip_lr=flip_lr,
            )
            _, _, aod_surface = _aod_grid_to_surface_via_camera(sensor_ds)

            panel_maps = [vza_map, vaa_map, aod_surface, raa_map, sca_angle]
            titles = ["VZA", "VAA", "AOD (grid→camera→surface)", "RAA", "Scattering Angle"]
            cmaps = ["viridis", "viridis", "viridis", "viridis", "viridis"]

            # Figure-1: camera image coordinates (pixel domain)
            fig_cam, axes_cam = plt.subplots(2, 3, figsize=(15, 8))
            for ax, data, title, cmap in zip(axes_cam.ravel(), panel_maps, titles, cmaps):
                im = ax.imshow(np.asarray(data), origin="lower", cmap=cmap)
                ax.set_title(f"{title} (camera)")
                ax.set_xlabel("x_pixel")
                ax.set_ylabel("y_pixel")
                fig_cam.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            axes_cam.ravel()[-1].axis("off")
            plt.tight_layout()
            plt.savefig(os.path.join(preview_dir, f"angles_panel_camera_{int(wavelength_nm)}nm_{view_name}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(fig_cam)

            # Figure-2: ground-projected coordinates
            fig, axes = plt.subplots(2, 3, figsize=(15, 8))
            xv, yv = centers_to_edges_2d(xg, yg)
            xlim = (np.nanmin(xg), np.nanmax(xg))
            ylim = (np.nanmin(yg), np.nanmax(yg))
            for ax, data, title, cmap in zip(axes.ravel(), panel_maps, titles, cmaps):
                im = ax.pcolormesh(yv, xv, data, shading="flat", cmap=cmap)
                ax.set_title(f"{title} (ground)")
                ax.set_xlabel("y_ground [km]")
                ax.set_ylabel("x_ground [km]")
                ax.set_aspect("equal", adjustable="box")
                ax.set_xlim(*ylim)
                ax.set_ylim(*xlim)
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            # 2x3 panel has one extra slot
            axes.ravel()[-1].axis("off")
            plt.tight_layout()
            plt.savefig(os.path.join(preview_dir, f"angles_panel_ground_{int(wavelength_nm)}nm_{view_name}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        fields = [AOD.get("aerosol"), SSA.get("aerosol")]
        titles = [f"AOD aerosol @ {int(wavelength_nm)} nm", f"SSA aerosol @ {int(wavelength_nm)} nm"]
        for ax, data, title in zip(axes, fields, titles):
            im = ax.imshow(np.asarray(data), origin="lower", cmap="viridis")
            ax.set_title(title)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig(os.path.join(preview_dir, f"aod_ssa_panel_{int(wavelength_nm)}nm.png"),
                    dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"🧪 Saved precheck maps to: {preview_dir}")

    for i, (w, pol) in enumerate(bands, start=1):
        print(f"\n==============================")
        print(f"🚀 Processing wavelength {w} nm  ({i}/{len(bands)})")
        print(f"==============================")
        start_t = time.time()
        nc_path = os.path.join(out_cfg.root_dir, f"AirMSPI_{int(w)}nm.nc")
        if out_cfg.save_nc and os.path.exists(nc_path) and not overwrite:
            print(f"⏭️  Skip existing {nc_path} (use --overwrite to regenerate)")
            continue
        (sensor_dict, solvers_dict, 
         context, medium, 
         AOD, SSA, 
         xlats, xlons) = build_scene_and_sensors_single_band(
            sen_cfg, w, pol, out_cfg.root_dir, scene_cfg, cam_cfg, solar_cfg, solver_cfg, aerosol_cfg)
        
        spectral_AOD_all[w] = AOD     # AOD 是一个 dict: {'cloud': 2-D, 'rayleigh': 2-D, 'total': 2-D}
        spectral_SSA_all[w] = SSA     # SSA 是一个 dict
        _save_precheck_panels(sensor_dict, context, AOD, SSA, w)
        if precheck_only:
            print(f"⏭️ precheck_only=True, skip RTE solve for {int(w)} nm.")
            continue

        sensor_dict.get_measurements(solvers_dict, n_jobs=n_jobs, verbose=True)
        ds = build_versions_single_band(sensor_dict, xlats, xlons, w, pol, sen_cfg, grd_cfg, dsm_cfg,
                                        PlotConfig(True, True, plt_cfg.colormap), out_cfg,
                                        subdirs, crop_cfg, context)
        if out_cfg.save_nc:
            encoding = {v: {"zlib": True, "complevel": 2} for v in ds.data_vars}
            ds.to_netcdf(nc_path, encoding=encoding)
            print(f"✅ Saved NetCDF: {nc_path}")
            
            
            
        if clean_after_band:
            del sensor_dict, solvers_dict, ds
        dur = time.time() - start_t
        print(f"⏱️  {int(w)} nm done in {dur/60:.2f} min")
    print("\n✅ All requested wavelengths processed.")
    if out_cfg.save_nc:
        if len(spectral_AOD_all) == 0:
            print("⚠️ Skip spectral_AOD_SSA.nc: no AOD/SSA data collected (all bands skipped or precheck_only mode).")
        else:
            ds = utils.build_multiband_xarray(spectral_AOD_all, spectral_SSA_all)
            out_path = os.path.join(out_cfg.root_dir, "spectral_AOD_SSA.nc")
            ds.to_netcdf(out_path)
            print(f"📦 Saved spectral AOD/SSA NetCDF to {out_path}")


__all__ = [
    "main",
    "build_scene_and_sensors_single_band",
    "build_versions_single_band",
    "calculate_up_vector",
    "calculate_sensor_trajectory",
    "project_to_ground_lookat",
    "make_ground_grid",
    "centers_to_edges_2d",
    "crop_by_world_box",
    "plot_on_ground",
    "inspect_simulation_npz",
    "plot_simulation_results",
    "replot_simulation_results_with_config",
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Per-wavelength AT3D simulation pipeline (v5b)")
    parser.add_argument("--config", type=str, default="config_v5b.yaml", help="Path to config file")
    parser.add_argument("--band", type=int, default=None, help="Run only this wavelength (nm)")
    parser.add_argument("--n_jobs", type=int, default=None, help="Number of parallel jobs (overrides config)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing NetCDF files")
    parser.add_argument("--no_clean", action="store_true", help="Do not cleanup per-band data (debug)")
    parser.add_argument("--precheck_only", action="store_true",
                        help="Only build geometry/medium and save precheck angle+AOD/SSA maps; skip RTE solve")
    args = parser.parse_args()
    main(cfg_path=args.config, only_band=args.band, n_jobs=args.n_jobs,
         overwrite=args.overwrite, clean_after_band=not args.no_clean,
         precheck_only=args.precheck_only)
