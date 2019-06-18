# Scripts

Scripts are divided into two categories

Render
--- 
Rendering scripts solve the forward problem of simulating images from atmospheric optical properties.
A rendering script takes a generator (see: shdom/generate.py for more info) as input and outputs radiance images of the medium using SHDOM.
The resulting measurements and ground-truth are saved for subsequent use by an optimize script.

Examples:
1. Render a single voxel atmosphere at 9 view angles, 100m resolution, 672nm
python scripts/render_monochromatic_radiance_toa.py experiments/single_voxel/monochromatic \
        --generator SingleVoxel --lwc 0.1 --reff 10.0 --domain_size 1.0 \
        --x_res 0.1 --y_res 0.1 --nx 10 --ny 10 --nz 10  \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --mie_table_path mie_tables/polydisperse/Water_672nm.scat --add_rayleigh

2. Render an LES cloud field (rico) at 9 view angles, 10m resolution, 672nm, with a rayleigh scattering atmosphere
python scripts/render_monochromatic_radiance_toa.py experiments/rico32x36x25/monochromatic \
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt \
        --x_res 0.01 --y_res 0.01 --mie_table_path mie_tables/polydisperse/Water_672nm.scat \
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40 

2. Render an LES cloud field (rico) at 9 view angles, 3 spectral bands, 10m resolution, with a rayleigh scattering atmosphere
python scripts/render_polychromatic_radiance_toa.py experiments/rico32x36x25/polychromatic\
        --generator LesFile --path synthetic_cloud_fields/jpl_les/rico32x36x25.txt \
        --x_res 0.01 --y_res 0.01 --wavelength 0.672 0.55 0.445  --solar_flux 1.553 1.853 1.868\
        --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 \
        --add_rayleigh --n_jobs 40


Optimize
---
Optimization scripts solve the inverse problem of recovering medium properties from radiance measurements. 
The inputs are measurements and the outputs are estimated medium parameters.

Examples:
1. Optimize extinction with the ground truth phase function, grid and cloud mask for precomputed measurements. Keep a tensorboard log.
python scripts/optimize_extinction.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase  --use_forward_mask \
        --init Homogeneous --extinction 0.01 --log log_name
        
2. Optimize extinction with the ground truth phase function and grid for precomputed measurements. Keep a tensorboard log.
python scripts/optimize_extinction.py \
        --input_dir experiments/rico32x36x25/monochromatic --add_rayleigh \
        --use_forward_grid  --use_forward_albedo --use_forward_phase   \
        --init Homogeneous --extinction 0.01 --radiance_threshold 0.055 --log log_name

3. Optimize microphysics with ground-truth effective variance and grid for precomputed measurements. Keep a tensorboard log.
python scripts/optimize_extinction.py \
        --input_dir experiments/rico32x36x25/polychromatic --add_rayleigh\
        --use_forward_grid  --use_forward_veff  --use_forward_mask \
        --init Homogeneous --lwc 0.001  --reff 15.0 --log log_name