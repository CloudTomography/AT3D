"""
Utility functions for at3d. These are not critical to its operation.
"""
import typing

import numpy as np
from collections import OrderedDict
import netCDF4 as nc
import xarray as xr
import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess
import os
import warnings
from IPython.display import display
from ipywidgets import fixed, interactive
from pathlib import Path

import at3d.core
import at3d.solver
import at3d.grid

def github_version():
    """
    Get github version.

    Returns
    -------
    github_version: str,
        The current GitHub version.

    Raises
    ------
    warning if there are uncomitted changes in at3d.
    """
    github_dir = Path(os.environ['at3d_DIR']).joinpath('at3d')
    uncomitted_changes = subprocess.check_output(["git", "diff", "--name-only", github_dir]).strip().decode('UTF-8')
    if uncomitted_changes:
        warnings.warn('There are uncomitted changes in the {}/at3d directories: {}'.format(github_dir, uncomitted_changes))
    github_version = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode('UTF-8')
    github_version = github_version + ' + uncomitted changes' if uncomitted_changes else github_version
    return github_version

def slider_select_file(dir, filetype=None):
    """
    Slider for interactive selection of a file

    Parameters
    ----------
    dir: str,
        The directory from which to choose.
    filetype:  str, optional,
        Filetypes to display. If filetype is None then all files in the directory are displayed.

    Returns
    -------
    file: interactive,
        An interactive slider utility.
    """

    def select_path(i, paths):
        print(paths[i])
        return paths[i]

    filetype = '*' if filetype is None else '*.' + filetype
    paths = [str(path) for path in Path(dir).rglob('{}'.format(filetype))]
    file = interactive(select_path, i=(0, len(paths)-1), paths=fixed(paths));
    display(file)
    return file

def get_phase_function(legcoef, angles, phase_elements='All'):
    """Calculates phase function from legendre tables.

    If a multi-dimensional table is passed then phase functions for all
    microphysical dimensions (including individual radii if available)
    will be sampled at the specified `angles`.

    Parameters
    ----------
    legcoef : xr.DataArray
        Contains the legendre/Wigner coefficients for the phase function.
        should be produced by mie.get_mie_mono or  mie.get_poly_table or
        in the same format.
    angles : array_like of floats
        scattering angles to sample the phase function at in degrees.
        should be a 1D array_like
    phase_elements : str, list/tuple of strings
        valid values are from P11, P22, P33, P44, P12, P34, or All.

    Returns
    -------
    phase_array : xr.DataArray
        Contains the phase function at the sampled `angles` for each of the provided
        legendre/Wigner series, for the specified `phase_elements`

    Raises
    ------
    ValueError
        If `phase_elements` is not composed of valid strings
    TypeError
        If `phase_elements` is not of type ``str``, ``tuple`` or ``list``.

    See Also
    --------
    mie.get_poly_table
    mie.get_mono_table

    Example
    -------
    >>> legcoef = mie_mono_table.legendre[:,:,50:55]
    #select 5 radii from the mie_mono_table.

    >>> phase_array = get_phase_function(legcoef,
                                         np.linspace(0.0,180.0,361),
                                         phase_elements='All')
    >>> phase_array
    <xarray.DataArray 'phase_function' (phase_elements: 6, scattering_angle: 361, radius: 5)>
    array([[[ 5.06317383e-03,  6.17436739e-03,  7.50618055e-03,
              9.09753796e-03,  1.09932469e-02],
            [ 5.06292842e-03,  6.17406424e-03,  7.50580709e-03,
              9.09707788e-03,  1.09926835e-02],
            [ 5.06219361e-03,  6.17315574e-03,  7.50468718e-03,
              9.09570046e-03,  1.09909941e-02],
            ...,
            [ 2.58962740e-03,  3.01506580e-03,  3.48956930e-03,
              4.01409063e-03,  4.58831480e-03],
            [ 2.58983485e-03,  3.01530003e-03,  3.48983146e-03,
              4.01438121e-03,  4.58863331e-03],
            [ 2.58990400e-03,  3.01537826e-03,  3.48991877e-03,
              4.01447807e-03,  4.58873948e-03]],

           [[ 5.06317383e-03,  6.17436739e-03,  7.50618055e-03,
              9.09753796e-03,  1.09932479e-02],
            [ 5.06292889e-03,  6.17406424e-03,  7.50580709e-03,
              9.09707882e-03,  1.09926844e-02],
            [ 5.06219361e-03,  6.17315574e-03,  7.50468718e-03,
              9.09570139e-03,  1.09909950e-02],
    ...
            [-4.16639125e-07, -4.87358477e-07, -5.66944379e-07,
             -6.55818667e-07, -7.54250607e-07],
            [-1.04164208e-07, -1.21844508e-07, -1.41741438e-07,
             -1.63960422e-07, -1.88568734e-07],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00]],

           [[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00],
            [ 4.72140493e-10,  6.52067178e-10,  8.93253749e-10,
              1.21347665e-09,  1.63433822e-09],
            [ 1.88840676e-09,  2.60805422e-09,  3.57272101e-09,
              4.85350737e-09,  6.53681553e-09],
            ...,
            [ 1.69794856e-09,  2.33928010e-09,  3.19921267e-09,
              4.34327951e-09,  5.85340043e-09],
            [ 4.24516838e-10,  5.84860882e-10,  7.99858901e-10,
              1.08589548e-09,  1.46345203e-09],
            [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              0.00000000e+00,  0.00000000e+00]]], dtype=float32)
    Coordinates:
      * phase_elements    (phase_elements) <U3 'P11' 'P22' 'P33' 'P44' 'P12' 'P34'
      * scattering_angle  (scattering_angle) float64 0.0 0.5 1.0 ... 179.5 180.0
      * radius            (radius) float32 0.1186 0.1224 0.1263 0.1303 0.1343
    """
    pelem_dict = {'P11': 1, 'P22': 2, 'P33': 3, 'P44': 4, 'P12': 5, 'P34': 6}
    if phase_elements == 'All':
        phase_elements = list(pelem_dict.keys())
    elif isinstance(phase_elements, (typing.List, typing.Tuple)):
        for element in phase_elements:
            if element not in pelem_dict:
                raise ValueError("Invalid value for phase_elements '{}' "
                                 "Valid values are '{}'".format(element, pelem_dict.keys()))
    elif phase_elements in pelem_dict:
        phase_elements = [phase_elements]
    else:
        raise TypeError("phase_elements argument should be either 'All' or a list/tuple of strings"
                        "from {}".format(pelem_dict.keys()))

    coord_sizes = {name:legcoef[name].size for name in legcoef.coords
                   if name not in ('stokes_index', 'legendre_index', 'table_index')}

    coord_arrays = [np.arange(size) for size in coord_sizes.values()]
    coord_indices = np.meshgrid(*coord_arrays, indexing='ij')
    flattened_coord_indices = [coord.ravel() for coord in coord_indices]

    phase_functions_full = []

    loop_max = 1
    if flattened_coord_indices:
        loop_max = flattened_coord_indices[0].shape[0]

    for i in range(loop_max):
        index = tuple([slice(0, legcoef.stokes_index.size),
                       slice(0, legcoef.legendre_index.size)] +
                      [coord[i] for coord, size in zip(flattened_coord_indices, coord_sizes.values())
                       if size > 1])
        single_legcoef = legcoef.data[index]

        phase_functions = []
        for phase_element in phase_elements:
            pelem = pelem_dict[phase_element]

            phase = at3d.core.transform_leg_to_phase(
                maxleg=legcoef.legendre_index.size - 1,
                nphasepol=6,
                legcoef=single_legcoef,
                pelem=pelem,
                nleg=legcoef.legendre_index.size - 1,
                nangle=len(angles),
                angle=angles
            )
            phase_functions.append(phase)
        phase_functions_full.append(np.stack(phase_functions, axis=0))

    small_coord_sizes = {name:size for name, size in coord_sizes.items() if size > 1}
    coords = {
        'phase_elements': np.array(phase_elements),
        'scattering_angle': angles
        }
    for name in coord_sizes:
        coords[name] = legcoef.coords[name]

    phase_functions_full = np.stack(phase_functions_full, axis=-1).reshape(
        [len(phase_elements), len(angles)] + list(small_coord_sizes.values())
    )

    phase_array = xr.DataArray(
        name='phase_function',
        dims=['phase_elements', 'scattering_angle'] + list(small_coord_sizes.keys()),
        data=phase_functions_full,
        coords=coords
    )
    return phase_array

def plot_cell_grid(solver, y, visualize=None):
    """
    Visualize the adaptive grid.

    Plots all cell lines in the X-Z plane within the nearest base grid cell
    to the specified y value.

    Parameters
    ----------
    solver : at3d.solver.RTE
        The solver object.
    y : float
        The Y-plane of the grid to visualize.
    visualize : str
        If None then only the adaptive grid lines are visualized.
        If 'extinct' then the cell averaged extinction is visualized.
        If 'adaptcrit' then the maximum adaptive grid splitting criterion
        for each cell is visualized.

    Notes
    -----
    The adaptive grid cells are plotted on top of each other without any transparency
    so the most recently formed ones (which should be highest resolution) are the
    ones that are shown.
    """
    diff = solver._ygrid-y
    yinds = np.where(np.abs(diff) == np.abs(diff).min())[0]
    if diff[yinds] < 0.0:
        y_low = yinds
        y_high = yinds+1
    else:
        y_low = yinds-1
        y_high = yinds
    fig, ax = plt.subplots(figsize=(8, 8))
    if visualize is not None:
        cellextinct, adaptcrit = at3d.core.output_cell_split(
            gridptr=solver._gridptr,
            gridpos=solver._gridpos,
            nstokes=solver._nstokes,
            total_ext=solver._total_ext,
            shptr=solver._shptr,
            source=solver._source,
            ncells=solver._ncells
       )

        if visualize == 'extinct':
            minv = cellextinct.min()
            maxv = cellextinct.max()*1.1
            label = 'Cell averaged Extinction'
        elif visualize == 'adaptcrit':
            if not solver.check_solved(verbose=False):
                raise at3d.exceptions.SHDOMError(
                    "at3d.solver.RTE object has to be solved to visualize "
                    "adaptive splitting."
                )
            adapt = np.max(adaptcrit, axis=0)
            minv = adapt.min()
            maxv = solver._splitacc
            label = 'Maximum cell splitting criterion.'

        cmap = plt.cm.gray_r
        norm = plt.Normalize(minv, maxv)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label=label)

    for IC in range(solver._ncells):
        X1 = solver._gridpos[0, solver._gridptr[0, IC]-1]
        Z1 = solver._gridpos[2, solver._gridptr[0, IC]-1]
        X2 = solver._gridpos[0, solver._gridptr[1, IC]-1]
        Z2 = solver._gridpos[2, solver._gridptr[1, IC]-1]
        X3 = solver._gridpos[0, solver._gridptr[5, IC]-1]
        Z3 = solver._gridpos[2, solver._gridptr[5, IC]-1]
        X4 = solver._gridpos[0, solver._gridptr[4, IC]-1]
        Z4 = solver._gridpos[2, solver._gridptr[4, IC]-1]
        y = solver._gridpos[1, solver._gridptr[:, IC]-1]
        if (y.min() >= solver._ygrid[y_low]-1e-6) & (y.max() <= solver._ygrid[y_high]+1e-6):
            if visualize == 'extinct':
                facecolor = cmap(norm(cellextinct[IC]))
            elif visualize == 'adaptcrit':
                facecolor = cmap(norm(adapt[IC]))
            else:
                facecolor = 'none'
            rect = patches.Rectangle((X1, Z1), X2-X1, Z4-Z1, linewidth=0.5,
                                     edgecolor='black', facecolor=facecolor)
            ax.add_patch(rect)
            #print(colors.to_rgba('black', vis[IC]))
    zdomain = solver._zgrid[-1] - solver._zgrid[0]
    xdomain = solver._xgrid[-1] - solver._xgrid[0]
    ax.set_xlim(-0.1*xdomain, 1.1*xdomain)
    ax.set_ylim(solver._zgrid[0]-0.1*zdomain, solver._zgrid[-1]+0.1*zdomain)
    ax.set_ylabel('Z (km)')
    ax.set_xlabel('X (km)')
    ax.set_title('Adaptive grid lines in the y=[{:2.3f}, {:2.3f}] plane'.format(
        solver._ygrid[y_low][0], solver._ygrid[y_high][0])
        )
    plt.show()


def planck_function(temperature, wavelength, c=2.99792458e8, h=6.62606876e-34, k=1.3806503e-23):
    """
    temperature
        units, Kelvin
    wavelength
        units, micrometers
    radiance
        units, Watts/m^2/micrometer/steradian (SHDOM units)
    """
    wavelength = wavelength*1e-6
    radiance = 2*h*c**2/ wavelength**5 *1.0/(np.exp((h*c)/(wavelength*k*temperature)) - 1.0)*1e-6
    return radiance

def cell_average_comparison(reference, other, variable_name):
    """
    calculates average values of 'variable name' in the cells
    of reference's grid for both reference and other (other is on a different grid.)
    """
    ref_vol, ref_val, other_vol, other_val = at3d.core.cell_average(
        xgrid1=reference.x.data,
        ygrid1=reference.y.data,
        zgrid1=reference.z.data,
        xgrid2=other.x.data,
        ygrid2=other.y.data,
        zgrid2=other.z.data,
        values1=reference[variable_name].data,
        values2=other[variable_name].data
    )
    cell_average_ref = np.zeros(ref_vol.shape)
    cell_average_ref[np.where(ref_vol > 0.0)] = ref_val[np.where(ref_vol > 0.0)] / \
        ref_vol[np.where(ref_vol > 0.0)]
    cell_average_other = np.zeros(other_vol.shape)
    cell_average_other[np.where(other_vol > 0.0)] = other_val[np.where(other_vol > 0.0)] \
        /other_vol[np.where(other_vol > 0.0)]
    return cell_average_ref, cell_average_other

def load_2parameter_lwc_file(file_name, density='lwc'):
    """
    TODO
    Function that loads a scatterer from the '2 parameter lwc file' format used by
    SHDOM and i3rc monte carlo model.
    """
    header = pd.read_csv(file_name, nrows=4)
    nx, ny, nz = np.fromstring(header['2 parameter LWC file'][0], sep=' ').astype(np.int)
    dx, dy = np.fromstring(header['2 parameter LWC file'][1], sep=' ').astype(np.float)
    z = np.fromstring(header['2 parameter LWC file'][2], sep=' ').astype(np.float)
    temperature = np.fromstring(header['2 parameter LWC file'][3], sep=' ').astype(np.float)
    dset = at3d.grid.make_grid(dx, nx, dy, ny, z)

    data = np.genfromtxt(file_name, skip_header=5)

    lwc = np.zeros((nx, ny, nz))*np.nan
    reff = np.zeros((nx, ny, nz))*np.nan

    i, j, k = data[:, 0].astype(np.int)-1, data[:, 1].astype(np.int)-1, data[:, 2].astype(np.int)-1
    lwc[i, j, k] = data[:, 3]
    reff[i, j, k] = data[:, 4]

    dset['density'] = xr.DataArray(
        data=lwc,
        dims=['x', 'y', 'z']
    )

    dset['reff'] = xr.DataArray(
        data=reff,
        dims=['x', 'y', 'z']
    )

    dset['temperature'] = xr.DataArray(
        data=temperature,
        dims=['z']
    )

    dset.attrs['density_name'] = density
    dset.attrs['file_name'] = file_name

    return dset

def to_2parameter_lwc_file(file_name, cloud_scatterer, atmosphere=None, fill_temperature=280.0):
    """
    TODO
    Write lwc & reff to the '2 parameter lwc' file format used by i3rc MonteCarlo model and SHDOM.
    atmosphere should contain the temperature. It is interpolated to the specified z grid.
    If no atmosphere is included then a fill_temperature is used (Temperature is required
    in the file).
    """

    nx, ny, nz = cloud_scatterer.density.shape
    dx, dy = (cloud_scatterer.x[1]-cloud_scatterer.x[0]).data, (cloud_scatterer.y[1] - cloud_scatterer.y[0]).data
    z = cloud_scatterer.z.data

    if atmosphere is not None:
        temperature = atmosphere.interp({'z': cloud_scatterer.z}).temperature.data
    else:
        temperature = np.ones(z.shape)*fill_temperature

    i, j, k = np.meshgrid(np.arange(1, nx+1), np.arange(1, ny+1), np.arange(1, nz+1), indexing='ij')

    lwc = cloud_scatterer.density.data.ravel()
    reff = cloud_scatterer.reff.data.ravel()

    z_string = ''
    for z_value in z:
        if z_value == z[-1]:
            z_string += '{}'.format(z_value)
        else:
            z_string += '{} '.format(z_value)

    t_string = ''
    for index, temp_value in enumerate(temperature):
        if index == len(temperature) - 1:
            t_string += '{:5.2f}'.format(temp_value)
        else:
            t_string += '{:5.2f} '.format(temp_value)

    with open(file_name, "w") as f:
        f.write('2 parameter LWC file\n')
        f.write(' {} {} {}\n'.format(nx, ny, nz))
        f.write('{} {}\n'.format(dx, dy))
        f.write('{}\n'.format(z_string))
        f.write('{}\n'.format(t_string))
        for x, y, z, l, r in zip(i.ravel(), j.ravel(), k.ravel(), lwc.ravel(), reff.ravel()):
            f.write('{} {} {} {:5.4f} {:3.2f}\n'.format(x, y, z, l, r))

def load_from_csv(path, density=None, origin=(0.0,0.0)):
    """
    TODO
    """
    df = pd.read_csv(path, comment='#', skiprows=4, index_col=['x', 'y', 'z'])
    nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int, delimiter=',')
    dx, dy = np.genfromtxt(path, max_rows=1, dtype=float, skip_header=2, delimiter=',')
    z = xr.DataArray(np.genfromtxt(path, max_rows=1, dtype=float, skip_header=3, delimiter=','), coords=[range(nz)], dims=['z'])

    dset = at3d.grid.make_grid(dx, nx, dy, ny, z)
    i, j, k = zip(*df.index)

    for name in df.columns:
        #initialize with np.nans so that empty data is np.nan
        variable_data = np.zeros((dset.sizes['x'], dset.sizes['y'], dset.sizes['z']))*np.nan
        variable_data[i, j, k] = df[name]
        dset[name] = (['x', 'y', 'z'], variable_data)

    if density is not None:
        assert density in dset.data_vars, \
        "density variable: '{}' must be in the file".format(density)

        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset

def load_from_netcdf(path, density=None):
    """
        TODO
    A shallow wrapper around open_dataset that sets the density_name.
    """
    dset = xr.open_dataset(path)

    if density is not None:
        if density not in dset.data_vars:
            raise ValueError("density variable: '{}' must be in the file".format(density))
        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset


def load_forward_model(file_name):
    """
    Load from netCDF the solvers and sensors that were saved using
    `save_forward_model` below.
    """
    dataset = nc.Dataset(file_name)

    groups = dataset.groups
    sensors = groups['sensors'].groups
    solvers = groups['solvers'].groups
    sensor_dict = at3d.containers.SensorsDict()
    solver_dict = at3d.containers.SolversDict()

    for key,sensor in sensors.items():
        sensor_list = []
        for i, image in sensor.groups.items():
            sensor_dataset = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'sensors/'+str(key)+'/'+str(i)]))
            sensor_dataset['stokes'] = (['stokes_index'],sensor_dataset['stokes'].data.astype(np.bool))
            sensor_dataset['use_subpixel_rays'] = sensor_dataset['use_subpixel_rays'].data.astype(np.bool)
            sensor_dict.add_sensor(key, sensor_dataset)
        #     sensor_list.append(xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
        #         'sensors/'+str(key)+'/'+str(i)])))
        # sensor_dict[key] = {'sensor_list':sensor_list}

    for key, solver in solvers.items():
        numerical_params = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/numerical_parameters']))
        numerical_params['deltam'] = numerical_params.deltam.data.astype(np.bool)
        numerical_params['high_order_radiance'] = numerical_params.high_order_radiance.data.astype(np.bool)
        num_stokes = numerical_params.num_stokes.data
        surface = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/surface']))
        source = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/source']))

        mediums = OrderedDict()
        for name, med in solver['medium'].groups.items():
            mediums[name] = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/medium/'+str(name)]))

        if 'atmosphere' in solver.groups.keys():
            atmosphere = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'solvers/'+str(key)+'/atmosphere']))
        else:
            atmosphere=None

        numerical_params['transcut'] = 1e-5
        numerical_params['angle_set'] = 2

        solver_dict.add_solver(float(key), at3d.solver.RTE(numerical_params=numerical_params,
                                            medium=mediums,
                                           source=source,
                                           surface=surface,
                                            num_stokes=num_stokes,
                                            name=None
                                           )
                                           )
        rte_grid = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/grid']))

    return sensor_dict, solver_dict, rte_grid

def save_forward_model(file_name, sensors, solvers):
    """
    Save sensors and solvers together to netcdf.

    This is useful for saving rendered measurements and the specifications that
    were used to make them. The SHDOM solution itself is not currently incorporated
    here. See solver.RTE.save_solution for how to save that.

    Parameters
    ----------
    file_name : str
        The path to the file that will be created and saved to.
        If it exists an integer is appended to the file name until
        it can be safely created.
    sensors : at3d.containers.SensorsDict
        The container with all of the sensors information.
    solvers : at3d.containers.SolversDict
        The container with all of the solver.RTE objects.
    """

    if not isinstance(sensors, at3d.containers.SensorsDict):
        raise TypeError(
            "`sensors` should be an instance of '{}'".format(at3d.containers.SensorsDict)
        )
    if not isinstance(solvers, at3d.containers.SolversDict):
        raise TypeError(
            "`sensors` should be an instance of '{}'".format(at3d.containers.SolversDict)
        )

    # make a safe file name. We always want one that works
    # as we might have spent a lot of time to make the input data.
    counter = 1
    initial_file_name = file_name
    while os.path.exists(file_name):
        if counter > 1:
            file_name = file_name[:-5] + '_{}'.format(counter) + '.nc'
        else:
            file_name = file_name[:-3] + '_{}'.format(counter) + '.nc'
        counter += 1

    if counter > 1:

        warnings.warn("file_name '{}' already exists, your file is now at '{}'".format(
            initial_file_name,file_name), category=RuntimeWarning
             )

    # intialize the save file.
    save_file = nc.Dataset(file_name, 'w')
    save_file.close()

    # do the mediums first in their own loop because otherwise the whole thing crashes?!?
    for key, solver in solvers.items():
        for name, med in solver.medium.items():
            med.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'medium/'+str(name))

    for i, (key,sensor) in enumerate(sensors.items()):
        for j, image in enumerate(sensor['sensor_list']):
            if (i==0) & (j==0):
                image.to_netcdf(file_name, 'a', group = 'sensors/'+str(key)+'/'+str(j))

            else:
                image.to_netcdf(file_name, 'a', group = 'sensors/'+str(key)+'/'+str(j))


    for i, (key, solver) in enumerate(solvers.items()):
        numerical_params = solver.numerical_params
        numerical_params['num_stokes'] = solver._nstokes
        solver.numerical_params.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'numerical_parameters')
        solver.surface.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'surface')
        solver.source.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'source')
        solver._grid.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'grid')
        if solver.atmosphere is not None:
            solver.atmosphere.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'atmosphere')

def adjoint_linear_interpolation(reference_coords, data_array):
    """
    Performs the adjoint of linear interpolation for 3D regular
    grids.

    Parameters
    ----------
    reference_coords : xr.Dataset, xr.DataArray
        A dataset or dataarray which has the coordinates which
        you want to adjointly interpolate data to. Should have 3D
        coordinates labeled 'x', 'y', 'z'.
    data_array : xr.DataArray
        data on the new coordinates to be interpolated back.
        Should have 3D coordinates labeled 'x', 'y', 'z'.

    Returns
    -------
    interpolated : xr.DataArray
        The interpolated data on the reference_coords
    """
    if not isinstance(reference_coords, (xr.Dataset, xr.DataArray)):
        raise TypeError(
            "`reference_coords` should be an xr.Dataset "
            "or xr.DataArray not of type '{}'".format(type(reference_coords))
        )
    if not isinstance(data_array, xr.DataArray):
        raise TypeError(
            "`data_array` should be an "
            "xr.DataArray not of type '{}'".format(type(data_array))
        )
    for dim_name in ('x', 'y', 'z'):
        if dim_name not in reference_coords.coords:
            raise KeyError(
                "Expected coordinate name '{}' in `reference_coords`".format(dim_name)
            )
        if dim_name not in data_array.coords:
            raise KeyError(
                "Expected coordinate name '{}' in `data_array`".format(dim_name)
            )
        if data_array[dim_name].min() < reference_coords[dim_name].min() or \
            data_array[dim_name].max() > reference_coords[dim_name].max():
            raise ValueError(
                "`data_array` '{}' coordinate is out of range of `reference_coords`.".format(
                    dim_name
                )
            )

    interpolated_array = at3d.core.adjoint_linear_interpolation(
        xgrid=reference_coords.x.data,
        ygrid=reference_coords.y.data,
        zgrid=reference_coords.z.data,
        inxgrid=data_array.x.data,
        inygrid=data_array.y.data,
        inzgrid=data_array.z.data,
        field=data_array.data
    )
    interpolated = xr.DataArray(
        data=interpolated_array,
        dims=['x','y','z'],
        coords={
            'x': reference_coords.x.data,
            'y': reference_coords.y.data,
            'z': reference_coords.z.data
        }
    )
    return interpolated

def save_transforms_to_file(unknown_scatterers, file_name):
    """
    Save the names and inputs of the transforms to netcdf for later use.

    This only works for objects which follow the convention that the
    inputs with 'name' are attributes as '_name'.

    """

    groups = []
    for scatterer_name, unknown_scatterer in unknown_scatterers.items():
        for variable_name, (coordinate_transform, state_to_grid_transform) in unknown_scatterer.variables.items():
            groups.append(scatterer_name + '/' + variable_name)

    group_dset = xr.Dataset(
        data_vars={
            'groups': (['ngroups'],np.array(groups).astype(str))
        }
    )
    group_dset.to_netcdf(file_name)

    for scatterer_name, unknown_scatterer in unknown_scatterers.items():
        for variable_name, (coordinate_transform, state_to_grid_transform) in unknown_scatterer.variables.items():
            group_name = scatterer_name + '/' + variable_name
            data_vars ={
                'state_to_grid_transform': type(state_to_grid_transform).__name__,
                'coordinate_transform': type(coordinate_transform).__name__,
            }
            for transform in (coordinate_transform, state_to_grid_transform):
                input_names = transform.__init__.__code__.co_varnames[1:]
                max_dim_name = 0
                for name in input_names:
                    data = getattr(transform, '_'+name)
                    data = np.array(data)
                    data_vars[name] = (['dim'+str(i+max_dim_name) for i in range(data.ndim)], data)
                    max_dim_name += data.ndim

            dset = xr.Dataset(data_vars=data_vars)
            dset.to_netcdf(file_name, group=group_name, mode='a')

def load_transforms_from_file(file_name):
    """
    Loads the transforms saved to netcdf.

    Only works for transforms defined in at3d.transforms as that module
    is searched for the correct transform object name to instantiate the new
    transform.
    """

    state_to_grid_transforms = {}
    coordinate_transforms = {}

    with xr.load_dataset(file_name) as out:
        for group in out.groups:
            group_name = str(group.data.astype(str))
            with xr.load_dataset(file_name, group=group_name) as nc:

                state_to_grid_transform = getattr(at3d.transforms, str(nc.state_to_grid_transform.data))
                input_names = state_to_grid_transform.__init__.__code__.co_varnames[1:]
                args = {name: nc[name].data for name in input_names}
                state_to_grid_transform = state_to_grid_transform(**args)
                state_to_grid_transforms[group_name] = state_to_grid_transform

                coordinate_transform = getattr(at3d.transforms, str(nc.coordinate_transform.data))
                input_names = coordinate_transform.__init__.__code__.co_varnames[1:]
                args = {name: nc[name].data for name in input_names}
                coordinate_transform = coordinate_transform(**args)
                coordinate_transforms[group_name] = coordinate_transform

    return state_to_grid_transforms, coordinate_transforms
