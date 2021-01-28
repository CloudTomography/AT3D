"""
This module contains the space carving routines and geometric projection
that can be used to initialize cloud masking.
"""
import xarray as xr
import numpy as np

import pyshdom.core
import pyshdom.checks
import pyshdom.solver
import pyshdom.util
import pyshdom.grid

class SpaceCarver:
    """
    An object that sets up an SHDOM grid structure through which
    rays can be traced to perform space carving and ray integration
    operations.
    """
    def __init__(self, grid):

        pyshdom.checks.check_grid(grid)
        self._grid = grid
        #different values are untested: changing these could result in failure of ray_tracing.
        self._setup_grid(self._grid, ipflag=0, bcflag=3)

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

        self._pa = pyshdom.solver.ShdomPropertyArrays()

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
        self._nbpts = self._nx * self._ny * self._nz

        # Calculate the number of base grid cells depending on the BCFLAG
        self._nbcells = (self._nz - 1) * (self._nx + ibits(self._bcflag, 0, 1) - ibits(self._bcflag, 2, 1)) * \
                (self._ny + ibits(self._bcflag, 1, 1) - ibits(self._bcflag, 3, 1))

        self._xgrid, self._ygrid, self._zgrid = pyshdom.core.new_grids(
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
        self._treeptr, self._cellflags = pyshdom.core.init_cell_structure(
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


    def carve(self, sensor_masks, agreement=None, linear_mode=False):
        """
        Performs a space carving operation on
        """
        if isinstance(sensor_masks, xr.Dataset):
            sensor_list = [sensor_masks]
        elif isinstance(sensor_masks, type([])):
            sensor_list = sensor_masks
        elif isinstance(sensor_masks, pyshdom.containers.SensorsDict):
            sensor_list = []
            for instrument in sensor_masks:
                sensor_list.extend(sensor_masks[instrument]['sensor_list'])

        volume = np.zeros((len(sensor_list), 3, self._nx, self._ny, self._nz))
        for i, sensor in enumerate(sensor_list):

            camx = sensor['ray_x'].data
            camy = sensor['ray_y'].data
            camz = sensor['ray_z'].data
            cammu = sensor['ray_mu'].data
            camphi = sensor['ray_phi'].data
            flags = sensor['cloud_mask'].data
            weights = sensor['weights'].data

            assert camx.ndim == camy.ndim == camz.ndim == cammu.ndim == camphi.ndim == 1
            total_pix = camphi.size

            carved_volume = pyshdom.core.space_carve(
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
                linear=linear_mode
            )
            volume[i] += carved_volume.reshape(3, self._nx, self._ny, self._nz)

        space_carved = xr.Dataset(
                        data_vars={
                            'cloudy_counts': (['nsensors', 'x', 'y', 'z'], volume[:, 0]),
                            'total_counts': (['nsensors', 'x', 'y', 'z'], np.sum(volume[:, :2], axis=1)),
                            'weights':(['nsensors', 'x', 'y', 'z'], volume[:, 2]),
                        },
                        coords={'x':self._grid.x,
                                'y':self._grid.y,
                                'z':self._grid.z}
                        )

        if agreement is not None:
            mask = ((space_carved.cloudy_counts/space_carved.total_counts > agreement[0]).astype(
                    np.int).mean('nsensors') >= agreement[1]).astype(np.int)
            masked_weights = (space_carved.weights/space_carved.total_counts).mean('nsensors')
            space_carved['mask'] = mask
            space_carved['masked_weights'] = masked_weights

        return space_carved

    def project(self, weights, sensors, return_matrix=False):
        """
        TODO
        """
        if isinstance(sensors, xr.Dataset):
            sensor_list = [sensors]
        elif isinstance(sensors, type([])):
            sensor_list = sensors
        elif isinstance(sensors, pyshdom.containers.SensorsDict):
            sensor_list = []
            for instrument in sensors:
                sensor_list.extend(sensors[instrument]['sensor_list'])
        weights = weights.copy(deep=True)
        weights = weights.drop_vars([name for name in weights.variables
                           if name not in ('x', 'y', 'z', 'density')])

        resampled_weights = pyshdom.grid.resample_onto_grid(self._grid, weights).density.data.ravel()
        for sensor in sensor_list:
            camx = sensor['ray_x'].data
            camy = sensor['ray_y'].data
            camz = sensor['ray_z'].data
            cammu = sensor['ray_mu'].data
            camphi = sensor['ray_phi'].data
            ray_weights = sensor['ray_weight'].data
            pixel_indices = sensor['pixel_index'].data
            total_pix = sensor.npixels.size
            nrays = sensor.nrays.size

            if return_matrix:
                raise NotImplementedError
            else:
                matrix_size = 1#(self._nx+self._ny+self._nz)
                matrix_ptrs = np.zeros((matrix_size, total_pix), dtype=np.int)+1
                matrix = np.zeros((matrix_size, total_pix))

            path, matrix, matrix_ptrs = pyshdom.core.project(
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
                matrix=matrix,
                matrix_ptrs=matrix_ptrs,
                return_matrix=False,
                weights=resampled_weights,
                ray_weights=ray_weights,
                pixel_indices=pixel_indices,
                nrays=nrays
            )
            sensor['integrated_weights'] = ('npixels', path)
