import pyshdom
import time
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path

try:
    from roipoly import RoiPoly # https://github.com/jdoepfert/roipoly.py
except ModuleNotFoundError as e:
    print(e)

directory = 'data/airmspi/'
pattern = 'AirMSPI_ER2_GRP_ELLIPSOID_20130206_*_NorthPacificOcean-32N123W_{view_angle}_F01_V005.hdf'
view_angles = ['661A', '478A', '291A', '000N', '291F', '478F', '589F', '661F']
datetime = time.strftime("%d-%b-%Y-%H:%M:%S")

for view_angle in view_angles:
    filepath = list(Path(directory).glob(pattern.format(view_angle=view_angle)))[0]
    dataset = pyshdom.util.load_airmspi_data(filepath, bands='660nm').expand_dims(view_angle=[view_angle])

    # Select region of interest according to 660nm image
    radiance = dataset.I.to_dataframe().dropna(how='all').to_xarray().sortby('XDim').sortby('YDim').I.squeeze()
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.imshow(radiance)
    roi = RoiPoly(fig=fig, ax=ax, color='r', close_fig=True)
    mask = xr.DataArray(roi.get_mask(radiance.values), coords=radiance.coords)


    # Crop ROI in all bands and save to netcdf with two groups.
    outpath = filepath.with_name('[{}]Cropped_'.format(datetime) + filepath.name).with_suffix('.nc')
    print('Saving data to file: {}'.format(outpath))
    rad_ds, pol_ds = pyshdom.util.load_airmspi_data(filepath)
    radiance = radiance.reset_coords(drop=True)
    rad_ds = rad_ds.expand_dims(view_angle=[view_angle]).sel(radiance.coords).sortby('XDim').sortby('YDim').where(mask, drop=True)
    pol_ds = pol_ds.expand_dims(view_angle=[view_angle]).sel(radiance.coords).sortby('XDim').sortby('YDim').where(mask, drop=True)
    rad_ds.to_netcdf(outpath, group='radiance')
    pol_ds.to_netcdf(outpath, group='polarization', mode='a')

