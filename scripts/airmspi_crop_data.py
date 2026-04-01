import at3d
import time
import argparse
import xarray as xr
from pathlib import Path
from datetime import datetime
import re

directory = 'data/airmspi/'
pattern = 'AirMSPI_ER2_GRP_ELLIPSOID_20130206_*_NorthPacificOcean-32N123W_{view_angle}_F01_V005.hdf'

current_time = time.strftime("%d-%b-%Y-%H:%M:%S")

parser = argparse.ArgumentParser()
parser.add_argument('--view_angles',
                    nargs='+',
                    type=str,
                    default=['661A', '589A', '478A', '291A', '000N', '291F', '478F', '589F', '661F'],
                    help='(default value: %(default)s)  angles to crop.')
args = parser.parse_args()

for view_angle in args.view_angles:
    filepath = list(Path(directory).glob(pattern.format(view_angle=view_angle)))[0]
    dataset = at3d.preprocessing.load_airmspi_data(filepath, bands='660nm')

    # Interactively select a region of interest according to 660nm image
    image = dataset.I.to_dataframe().dropna(how='all').to_xarray().sortby('XDim').sortby('YDim').I.squeeze()
    roi = at3d.preprocessing.get_roi(image)
    mask = xr.DataArray(roi.get_mask(image.values), coords=image.coords)


    # Crop ROI in all bands and save to netcdf with two groups.
    outpath = filepath.with_name('{}_Cropped_'.format(current_time) + filepath.name).with_suffix('.nc')
    print('Saving data to file: {}'.format(outpath))
    rad_ds, pol_ds = at3d.preprocessing.load_airmspi_data(filepath)
    image = image.reset_coords(drop=True)
    rad_ds = rad_ds.expand_dims(view_angle=[view_angle]).sel(image.coords).sortby('XDim').sortby('YDim').where(mask, drop=True)
    pol_ds = pol_ds.expand_dims(view_angle=[view_angle]).sel(image.coords).sortby('XDim').sortby('YDim').where(mask, drop=True)

    # Use regular expressions to extract date and time from the path and add it to dataset
    time = xr.DataArray(datetime.strptime(re.findall(r"\d{8}_\d{6}Z", str(filepath))[0][:-1], "%Y%m%d_%H%M%S"),
                        coords={'view_angle': [view_angle]}, dims=['view_angle'], name='time')
    rad_ds = xr.merge((rad_ds, time))
    pol_ds = xr.merge((pol_ds, time))

    rad_ds.to_netcdf(outpath, group='radiance')
    pol_ds.to_netcdf(outpath, group='polarization', mode='a')
