# Scripts
--- 
Scripts are divided into three categories: generate, render, optimize.
To learn more about each indivudual script command line flags use
```sh
python script.py --help
```

## Generate
Generation scripts create atmospheric components to be used for rendering or inversion.

**Generate Mie scattering tables for AirMSPI wavelengths**
```sh
python scripts/generate_mie_tables.py  \
    --start_reff 4.0 --end_reff 25.0 --num_reff 50 --start_veff 0.01 --end_veff 0.2 --num veff 50 \
    --radius_cutoff 65.0 --wavelength 0.355 0.38 0.445 0.47 0.555 0.66 0.865 0.935
```


## Render
Rendering scripts solve the forward problem of simulating images from atmospheric optical properties
A render script takes a generator as input (see: shdom/generate.py for more info) and outputs radiance images of the medium using SHDOM rte solver. The resulting measurements and ground-truth atmospheric parameters are saved for subsequent use by an optimize script.

**Render a single voxel atmosphere at 9 view angles, 20m resolution, 672nm**
```sh
python scripts/render_monochromatic_radiance_toa.py experiments/single_voxel/monochromatic \
        --generator SingleVoxel --lwc 0.1 --reff 10.0 --domain_size 1.0 \
        --x_res 0.02 --y_res 0.02 --nx 10 --ny 10 --nz 10  \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --mie_table_path mie_tables/polydisperse/Water_672nm.scat --add_rayleigh
```

**Render an LES cloud field (rico) at 9 view angles, 10m resolution, 672nm, with a rayleigh scattering atmosphere**
```sh
python scripts/render_monochromatic_radiance_toa.py experiments/rico32x36x25/monochromatic \
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt \
        --x_res 0.01 --y_res 0.01 --mie_table_path mie_tables/polydisperse/Water_672nm.scat \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40 
```

**Render a single voxel atmosphere at 9 view angles, 5 spectral bands 20m resolution**
```sh
python scripts/render_polychromatic_radiance_toa.py experiments/single_voxel/polychromatic \
        --generator SingleVoxel --nx 10 --ny 10 --nz 10 --domain_size 1.0 --x_res 0.02 --y_res 0.02 \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --wavelength 0.935 0.865  0.672 0.55 0.445  --solar_flux 0.86 1.008 1.553 1.853 1.868 \
        --solar_zenith 165 --solar_azimuth 90 --lwc 0.135 --reff 12.3 --add_rayleigh --n_jobs 40
```        
        
**Render an LES cloud field (rico) at 9 view angles, 3 spectral bands, 10m resolution, with a rayleigh scattering atmosphere**
```sh
python scripts/render_polychromatic_radiance_toa.py experiments/rico32x36x25/polychromatic \
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt \
        --x_res 0.01 --y_res 0.01 --wavelength 0.672 0.55 0.445  --solar_flux 1.553 1.853 1.868 \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40
```

&nbsp;

## Optimize

Optimization scripts solve the inverse problem of recovering medium properties from radiance measurements. 
The inputs are measurements and the outputs are estimated medium properties. Below are some example usage of the optimize scripts.

**Optimize extinction with the ground truth phase function, grid and cloud mask for precomputed measurements**
```sh
python scripts/optimize_extinction.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase  --use_forward_mask \
        --init Homogeneous --extinction 0.01 --log log_name --n_jobs 40 
```  
        
**Optimize extinction with the ground truth phase function and grid for precomputed measurements**
```sh 
python scripts/optimize_extinction.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase   \
        --init Homogeneous --extinction 0.01 --radiance_threshold 0.055 --log log_name --n_jobs 40 
```


**Optimize a single voxel for unknown microphysical properties**
```sh
python scripts/render_polychromatic_radiance_toa.py experiments/single_voxel/polychromatic \
        --generator SingleVoxel  --nx 10 --ny 10 --nz 10   --domain_size 1.0 \
        --x_res 0.1 --y_res 0.1 --wavelength 0.672 0.55 0.445 --solar_flux 1.553 1.853 1.868 \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --lwc 0.135 --reff 12.3 --add_rayleigh --n_jobs 40
 ```      
        
**Optimize microphysics with ground-truth effective variance and grid for precomputed measurements**
```sh
python scripts/optimize_microphysics.py \
        --input_dir experiments/rico32x36x25/polychromatic --add_rayleigh\
        --use_forward_grid  --use_forward_veff  --use_forward_mask \
        --init Homogeneous --lwc 0.001  --reff 15.0 --log log_name --n_jobs 40 
```
    