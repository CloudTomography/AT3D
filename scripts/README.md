# Scripts

Scripts are divided into two categories

--- 
1. Render
Rendering scripts solve the forward problem of simulating images from atmospheric optical properties.
A rendering script takes a generator (see: shdom/generate.py for more info) as input and outputs radiance images of the medium using SHDOM.
The resulting measurements and ground-truth are saved for subsequent use by an optimize script.

Examples:
1. Render an LES cloud field (rico) at 9 view angles, 10m resolution, 672nm, with a rayleigh scattering atmosphere
python scripts/render_radiance_toa.py rico32x36x25 --generator LesFile --path synthetic_cloud_fields/jpl_les/32x36x25.txt --x_res 0.01 --y_res 0.01  --azimuth 90 90 90 90 0 -90 -90 -90 -90 --zenith 70.5 60 45.6 26.1 0.0 26.1 45.6 60 70.5 --n_jobs 30 --mie_table_path mie_tables/polydisperse/Water_672nm.scat

---

2. Optimize
Optimization scripts solve the inverse problem of recovering medium properties from radiance measurements. 
The inputs are measurements and the outputs are estimated medium parameters.

Examples:
1. Optimize extinction with the ground truth phase function and grid for precomputed measurements. Keep a tensorboard log.
python scripts/optimize_extinction.py --input_dir  experiments/rico32x36x25/  --mie_table_path mie_tables/polydisperse/Water_672nm.scat --add_rayleigh --wavelength 0.672 --use_forward_grid  --use_forward_albedo --use_forward_phase  --init Homogeneous --extinction 0.01 --log radiance_mask_treshold_0.0585 
