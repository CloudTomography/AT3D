# Scripts 
Scripts are divided into three categories: generate, render, optimize.
To learn more about each individual script command line flags use
```sh
python script.py --help
```
&nbsp;

## Generate
Generation scripts create atmospheric components to be used for rendering or inversion.

Generate Mie scattering tables for AirMSPI wavelengths
```sh
python scripts/generate_mie_tables.py  \
    --start_reff 4.0 --end_reff 25.0 --num_reff 50 --start_veff 0.01 --end_veff 0.2 --num veff 50 \
    --radius_cutoff 65.0 --wavelength 0.355 0.38 0.445 0.47 0.555 0.66 0.865 0.935
```

&nbsp;

## Render
Rendering scripts solve the forward problem of simulating images from atmospheric optical properties
A render script takes a generator as input (see: shdom/generate.py for more info) and outputs radiance images of the medium using SHDOM rte solver. The resulting measurements and ground-truth atmospheric parameters are saved for subsequent use by an optimize script.

Render a single voxel atmosphere at 9 view angles, 20m resolution, 672nm
```sh
python scripts/render_radiance_toa.py experiments/single_voxel/monochromatic 0.672\
        --generator SingleVoxel --domain_size 1.0 --x_res 0.02 --y_res 0.02 --nx 10 --ny 10 --nz 10  \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --extinction 5.0 --reff 10.0 --add_rayleigh
```

Render an LES cloud field (rico) at 9 view angles, 10m resolution, 672nm, with a rayleigh scattering atmosphere and parallelization
```sh
python scripts/render_radiance_toa.py experiments/rico32x36x25/monochromatic 0.672 \
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt \
        --x_res 0.01 --y_res 0.01 --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40
```

Render a single voxel atmosphere at 9 view angles, 5 spectral bands 20m resolution
```sh
python scripts/render_radiance_toa.py experiments/single_voxel/polychromatic 0.935 0.865  0.672 0.55 0.445\
        --generator SingleVoxel --nx 10 --ny 10 --nz 10 --domain_size 1.0 --x_res 0.02 --y_res 0.02 \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --solar_zenith 165 --solar_azimuth 90 --lwc 0.135 --reff 12.3 --veff 0.1 --add_rayleigh --n_jobs 40
```        
        
Render an LES cloud field (rico) at 9 view angles, 3 spectral bands, 10m resolution, with a rayleigh scattering atmosphere
```sh
python scripts/render_radiance_toa.py experiments/rico32x36x25/polychromatic 0.672 0.55 0.445 \
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt --x_res 0.01 --y_res 0.01 \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40
```
Render a Stochastically generated 'Blob' cloud at 9 view angles, 3 spectral bands, 10m resolution, with a rayleigh scattering atmosphere
```sh
python scripts/render_radiance_toa.py experiments/stochastic_blob/polychromatic 0.672 0.55 0.445 \
        --generator StochasticCloud --x_res 0.01 --y_res 0.01  \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40 --cloud_geometry Blob --domain_size 1.0 --nx 25 --ny 25 --nz 25 --beta -1.6667 \
        --lwc_mean 0.15 --lwc_snr 0.5 --reff_mean 10.0 --reff_snr 0.1 --solar_zenith 135 --aspect_ratio 1.0
```

Render a Stochastically generated 'Cuboid' cloud at 9 view angles, 3 spectral bands, 10m resolution, with a rayleigh scattering atmosphere. Vertical:horizontal aspect ratio is 2
```sh
python scripts/render_radiance_toa.py experiments/stochastic_cuboid/polychromatic 0.672 0.55 0.445 \
        --generator StochasticCloud --x_res 0.01 --y_res 0.01  \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40 --cloud_geometry Cuboid --domain_size 1.0 --nx 25 --ny 25 --nz 25 --beta -1.6667 \
        --lwc_mean 0.15 --lwc_snr 0.5 --reff_mean 10.0 --reff_snr 0.1 --solar_zenith 135 --aspect_ratio 2.0
```

&nbsp;

## Optimize

Optimization scripts solve the inverse problem of recovering medium properties from radiance measurements. 
The inputs are measurements and the outputs are estimated medium properties. Below are some example usage of the optimize scripts.

#### Estimate Extinction

---
Estimate extinction with the ground truth phase function, grid and cloud mask for precomputed LES measurements
```sh
python scripts/optimize_extinction_lbfgs.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase  --use_forward_mask \
        --init Homogeneous --extinction 0.01 --log log_name --n_jobs 40 
```  
        
Estimate extinction with the ground truth phase function and grid for precomputed LES measurements
```sh 
python scripts/optimize_extinction_lbfgs.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase   \
        --init Homogeneous --extinction 0.01 --radiance_threshold 0.055 --log log_name --n_jobs 40 
```

### Estimate Micro-physical Properties

---
Estimate lwc with ground-truth effective radius and effective variance for precomputed single voxel measurements
```sh
python scripts/optimize_microphysics_lbfgs.py \
        --input_dir experiments/single_voxel/polychromatic --add_rayleigh \
        --use_forward_grid --use_forward_mask --use_forward_veff --use_forward_reff \
        --init Homogeneous --lwc 0.01 --n_jobs 40
```
 
Estimate effective radius with ground-truth lwc and effective variance for precomputed single voxel  measurements
```sh
python scripts/optimize_microphysics_lbfgs.py \
        --input_dir experiments/single_voxel/polychromatic --add_rayleigh \
        --use_forward_grid --use_forward_mask --use_forward_lwc --use_forward_veff \
        --init Homogeneous --reff 10.0 --n_jobs 40 
```
    
Estimate effective variance with ground-truth lwc and effective radius for precomputed single voxel  measurements
```sh
python scripts/optimize_microphysics_lbfgs.py \
        --input_dir experiments/single_voxel/polychromatic --add_rayleigh \
        --use_forward_grid --use_forward_mask --use_forward_lwc --use_forward_reff \
        --init Homogeneous --veff 0.1 --n_jobs 40 
```