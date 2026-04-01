"""
This module contains the space carving routines and geometric projection
that can be used to initialize cloud masking.
"""
import warnings
import xarray as xr
import numpy as np

import at3d.core
import at3d.checks
import at3d.solver
import at3d.util
import at3d.grid

class SpaceCarver:
    """
    Performs a space carving operation to mask points in a volume based on masked
    projections of the volume.

    This can be used to estimate the cloudy volume from multi-angle imager cloud
    masks.
    Uses the SHDOM grid structure and ray tracing to perform the space carving.
    A grid point is flagged by a ray's mask if the ray passes through a cell to
    which the gridpoint is adjacent. This is also the criteria for the ray having a
    non-zero contribution to the gradient at that point.
    This object can also perform line integrals through an arbitrary field defined
    on the SHDOM grid (linear interpolation kernel) as well as the adjoint operation
    (backprojection of ray weights to the grid points). These functions can be used
    to perform matrix free linear tomography.

    Parameters
    ----------
    grid : xr.Dataset
        A valid SHDOM grid object containing the x, y, z coordinates. See grid.py
        for more details.
    bcflag : int
        Sets the horizontal boundary conditions for the ray tracing. The default
        is 3 for open boundary conditions, 0 for periodic in both horizontal
        directions.
    """
    def __init__(self, grid, bcflag=3):
        at3d.checks.check_grid(grid)
        self._grid = grid
        if np.any([var in self._grid.data_vars for var in ('nx', 'ny', 'nz')]):
            warnings.warn(
                "SHDOM grid resolution parameters (nx/ny/nz) were detected in "
                " in `grid`. SpaceCarver does not support grids other than the property grid "
                " so these are ignored."
            )
        #different values are untested: changing these could result in failure of ray_tracing.
        #biggest issue is that NPTS differs from property grid points even on base grid
        #when periodic boundaries are used. So don't do that.
        #But 2D/1D mode SHOULD work (untested) though have to be careful about
        #where rays enter the domain.
        self._setup_grid(self._grid, ipflag=0, bcflag=bcflag)

    def _setup_grid(self, grid, ipflag, bcflag):
        """
        Set the base grid and define the SHDOM grid structures
        for ray tracing.

        Parameters
        ----------
        grid : xr.Dataset
            Dataset containing x, y, z coordinates. must be a valid
            SHDOM grid.
        ipflag : int
            This flag sets the dimensionality of the radiative transfer 1D/2D/3D etc.
            0=3D. see default_config.json and shdom.txt for definition.
        bcflag : int
            This flag sets the horizontal boundary conditions as either periodic
            or open. 3=open in both dimensions.
            see default_config.json and shdom.txt for definition.
        """

        def ibits(val, bit, ret_val):
            if val & 2 ** bit:
                return ret_val
            return 0

        self._ipflag = ipflag
        self._bcflag = bcflag

        self._pa = at3d.solver.ShdomPropertyArrays()

        # Set shdom property array
        self._pa.npx = grid.dims['x']
        self._pa.npy = grid.dims['y']
        self._pa.npz = grid.dims['z']
        self._pa.xstart = grid.x[0]
        self._pa.ystart = grid.y[0]

        self._pa.delx = grid.delx
        self._pa.dely = grid.dely
        self._pa.zlevels = grid.z

        # Initialize shdom internal grid sizes to property array grid
        self._nx = self._pa.npx
        self._ny = self._pa.npy
        self._nz = self._pa.npz
        self._maxpg = grid.dims['x'] * grid.dims['y'] * grid.dims['z']
        self._maxnz = grid.dims['z']

        # Set the full domain grid sizes (variables end in t)
        self._nxt, self._nyt, self._npxt, self._npyt = self._nx, self._ny, self._pa.npx, self._pa.npy
        self._nx, self._ny, self._nz = max(1, self._nx), max(1, self._ny), max(2, self._nz)

        # gridtype = 'P': Z levels taken from property file
        self._gridtype = 'P'

        # Set up base grid point actual size (NX1xNY1xNZ)
        self._nx1, self._ny1 = self._nx + 1, self._ny + 1
        if self._bcflag & 5 or ibits(self._ipflag, 0, 1):
            self._nx1 -= 1
        if self._bcflag & 7 or ibits(self._ipflag, 1, 1):
            self._ny1 -= 1
        self._nbpts = self._nx1 * self._ny1 * self._nz

        # Calculate the number of base grid cells depending on the BCFLAG
        self._nbcells = (self._nz - 1) * (self._nx + ibits(self._bcflag, 0, 1) - ibits(self._bcflag, 2, 1)) * \
                (self._ny + ibits(self._bcflag, 1, 1) - ibits(self._bcflag, 3, 1))

        self._xgrid, self._ygrid, self._zgrid = at3d.core.new_grids(
            bcflag=self._bcflag,
            gridtype=self._gridtype,
            npx=self._pa.npx,
            npy=self._pa.npy,
            nx=self._nx,
            ny=self._ny,
            nz=self._nz,
            xstart=self._pa.xstart,
            ystart=self._pa.ystart,
            delxp=self._pa.delx,
            delyp=self._pa.dely,
            zlevels=self._pa.zlevels
        )

        self._npts, self._ncells, self._gridpos, self._gridptr, self._neighptr, \
        self._treeptr, self._cellflags = at3d.core.init_cell_structure(
            maxig=self._nbpts,
            maxic=1.1*self._nbpts,
            bcflag=self._bcflag,
            ipflag=self._ipflag,
            nx=self._nx,
            ny=self._ny,
            nz=self._nz,
            nx1=self._nx1,
            ny1=self._ny1,
            xgrid=self._xgrid,
            ygrid=self._ygrid,
            zgrid=self._zgrid
        )
        self._nbcells = self._ncells


    def carve(self, sensor_masks, agreement=None, linear_mode=True):
        """
        Performs a space carving operation using the binary masks from the
        sensors.

        Parameters
        ----------
        sensor_masks : List/Tuple of xr.Dataset or at3d.containers.SensorsDict
            A group of datasets containing geometry of the rays and also their
            binary masks. See sensor.py for more details.
        agreement : Tuple/List
            Sets the logic for defining when a gridpoint is masks based on
            the ray masks from each of the elements of `sensor_masks`.
            If None, then the individual gridpoint masks are returned to be
            combined by the user.
            Otherwise, the first element of agreement is the threshold fraction of masked
            to total rays that intersect each grid point. If greater, then the
            gridpoint is set as cloudy for the contribution from that SENSOR
            (elements in `sensor_masks`).
            The second element sets the fraction of SENSORS that have to be masked
            at each gridpoint for that gridpoint to be masked.
        linear_mode : bool
            A flag which decides whether a linear or nearest neighbor interpolation
            kernel is used for the back propagation of weights. This is the
            adjoint of the SpaceCarver.project method below.

        Returns
        -------
        space_carved : xr.Dataset
            Contains the masked counts and total counts of ray intersections with
            grid points as well as the average weight, for each element of `sensor_masks`.
            Additionally, if `agreement` is supplied a 'mask' variable is present based
            on thresholding and averaging of the volumetric masks from each of the
            sensors.
        """
        if isinstance(sensor_masks, xr.Dataset):
            sensor_list = [sensor_masks]
        elif isinstance(sensor_masks, type([])):
            sensor_list = sensor_masks
        elif isinstance(sensor_masks, at3d.containers.SensorsDict):
            sensor_list = []
            for instrument in sensor_masks:
                sensor_list.extend(sensor_masks[instrument]['sensor_list'])

        counts = np.zeros((len(sensor_list), 2, self._nx, self._ny, self._nz), dtype=int)
        adjoint_weights = np.zeros((len(sensor_list), self._nx, self._ny, self._nz))
        for i, sensor in enumerate(sensor_list):

            camx = sensor['ray_x'].data
            camy = sensor['ray_y'].data
            camz = sensor['ray_z'].data
            cammu = sensor['ray_mu'].data
            camphi = sensor['ray_phi'].data
            flags = sensor['cloud_mask'].data
            if 'distance_limits' in sensor.data_vars:
                distance_limits = sensor['distance_limits'].data
            else:
                distance_limits = np.zeros(flags.shape)
            if 'weights' in sensor.data_vars:
                weights = sensor['weights'].data
            else:
                weights = np.ones(flags.shape)

            assert camx.ndim == camy.ndim == camz.ndim == cammu.ndim == camphi.ndim == 1
            total_pix = camphi.size

            carved_volume, count = at3d.core.space_carve(
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                npts=self._npts,
                ncells=self._ncells,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                gridpos=self._gridpos,
                camx=camx,
                camy=camy,
                camz=camz,
                cammu=cammu,
                camphi=camphi,
                npix=total_pix,
                flags=flags,
                weights=weights,
                linear=linear_mode,
                distance_limits=distance_limits
            )
            grid_counts = count.reshape(2, self._nx1, self._ny1, self._nz)
            grid_carved_volume = carved_volume.reshape(self._nx1, self._ny1, self._nz)
            if self._nx1 != self._nx:
                counts[i] += grid_counts[:, :-1, :-1]
                adjoint_weights[i] += grid_carved_volume[:-1, :-1]
            else:
                counts[i] += grid_counts
                adjoint_weights[i] += grid_carved_volume

        space_carved = xr.Dataset(
                        data_vars={
                            'cloudy_counts': (['nsensors', 'x', 'y', 'z'], counts[:, 0]),
                            'total_counts': (['nsensors', 'x', 'y', 'z'], np.sum(counts[:, :2], axis=1)),
                            'weights':(['nsensors', 'x', 'y', 'z'], adjoint_weights),
                        },
                        coords={'x':self._grid.x,
                                'y':self._grid.y,
                                'z':self._grid.z}
                        )

        if agreement is not None:
            mask = ((space_carved.cloudy_counts/space_carved.total_counts > agreement[0]).astype(
                int).mean('nsensors') >= agreement[1]).astype(int)
            masked_weights = (space_carved.weights/space_carved.total_counts).mean('nsensors')
            space_carved['mask'] = mask
            space_carved['masked_weights'] = masked_weights

        return space_carved

    def project(self, weights, sensors):
        """Performs line integrations assuming a linear interpolation kernel.

        The supplied `weights` are first interpolated onto the grid used by the
        SpaceCarver then line integrations are performed for each ray in the `sensors`.
        The `sensors` are modified in place to include the result of the line integrations.

        Parameters
        ----------
        weights : xr.Dataset
            A dataset containing a 'density' variable which will be integrated assuming
            linear interpolation between points.
        sensors : List/Tuple or at3d.containers.SensorsDict or xr.Dataset
            Contains the ray geometry for defining the line integrations. See
            sensor.py for more details.
        """
        if isinstance(sensors, xr.Dataset):
            sensor_list = [sensors]
        elif isinstance(sensors, type([])):
            sensor_list = sensors
        elif isinstance(sensors, at3d.containers.SensorsDict):
            sensor_list = []
            for instrument in sensors:
                sensor_list.extend(sensors[instrument]['sensor_list'])
        weights = weights.copy(deep=True)
        weights = weights.drop_vars(
            [name for name in weights.variables if name not in ('x', 'y', 'z', 'density')]
            )

        resampled_weights = at3d.grid.resample_onto_grid(self._grid, weights).density.data.ravel()
        for sensor in sensor_list:
            camx = sensor['ray_x'].data.astype(np.float64)
            camy = sensor['ray_y'].data.astype(np.float64)
            camz = sensor['ray_z'].data.astype(np.float64)
            cammu = sensor['ray_mu'].data.astype(np.float64)
            camphi = sensor['ray_phi'].data.astype(np.float64)
            ray_weights = sensor['ray_weight'].data
            pixel_indices = sensor['pixel_index'].data
            total_pix = sensor.npixels.size
            nrays = sensor.nrays.size

            path = at3d.core.project(
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                npts=self._npts,
                ncells=self._ncells,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                gridpos=self._gridpos,
                camx=camx,
                camy=camy,
                camz=camz,
                cammu=cammu,
                camphi=camphi,
                npix=total_pix,
                weights=resampled_weights,
                ray_weights=ray_weights,
                pixel_indices=pixel_indices,
                nrays=nrays
            )
            sensor['weights'] = ('npixels', path)

    def shadow_mask(self, volume_mask, sensors, solar_mu, solar_azimuth):
        """
        TODO
        """
        if isinstance(sensors, xr.Dataset):
            sensor_list = [sensors]
        elif isinstance(sensors, type([])):
            sensor_list = sensors
        elif isinstance(sensors, at3d.containers.SensorsDict):
            sensor_list = []
            for instrument in sensors:
                sensor_list.extend(sensors[instrument]['sensor_list'])

        weights = np.zeros(volume_mask.shape)
        weights[np.where(volume_mask.data > 0.0)] = 1.0
        weights_data = xr.Dataset(
            data_vars={
                'density': (['x','y','z'],weights)
            },
            coords={
            'x': self._grid.x,
            'y': self._grid.y,
            'z': self._grid.z
            }
        )
        grid_sensor = at3d.sensor.make_sensor_dataset(wavelength=0.672,
            x=self._gridpos[0,:], y=self._gridpos[1,:],
            z=self._gridpos[2,:],
            mu=np.array([solar_mu]*self._npts),phi=np.array([solar_azimuth]*self._npts),
                                            stokes=['I'],
                                                   fill_ray_variables=True)

        self.project(weights_data, grid_sensor)
        sundistance = grid_sensor.weights.data
        sundistance[np.where(weights_data.density.data.ravel() == 0.0)] = 0.0
        for sensor in sensor_list:
            camx = sensor['ray_x'].data.astype(np.float64)
            camy = sensor['ray_y'].data.astype(np.float64)
            camz = sensor['ray_z'].data.astype(np.float64)
            cammu = sensor['ray_mu'].data.astype(np.float64)
            camphi = sensor['ray_phi'].data.astype(np.float64)
            ray_weights = sensor['ray_weight'].data
            pixel_indices = sensor['pixel_index'].data
            total_pix = sensor.npixels.size
            nrays = sensor.nrays.size

            paths = at3d.core.get_shadow(
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                npts=self._npts,
                ncells=self._ncells,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                gridpos=self._gridpos,
                camx=camx,
                camy=camy,
                camz=camz,
                cammu=cammu,
                camphi=camphi,
                npix=total_pix,
                weights=sundistance.ravel(),
                ray_weights=ray_weights,
                pixel_indices=pixel_indices,
                nrays=nrays
            )
            sensor['sun_distance'] = (['npixels'], paths)

        sundistance_grid = xr.Dataset(
            data_vars={
                'sun_distance': (
                    ['x', 'y', 'z'],sundistance.reshape((self._nx1, self._ny1, self._nz)))
            },
            coords={
                'x': self._grid.x,
                'y': self._grid.y,
                'z': self._grid.z
            }
        )
        return sundistance_grid
