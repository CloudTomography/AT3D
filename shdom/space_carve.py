import xarray as xr
import numpy as np
import shdom

class SpaceCarver(object):

    @shdom.checks.dataset_checks(grid=shdom.checks.check_grid)
    def __init__(self, grid):

        self._grid = grid
        self._setup_grid(self._grid, ipflag=0,bcflag=3) #changing these could result in failure of ray_tracing.

    def _setup_grid(self, grid, ipflag,bcflag):
        """
        Set the base grid and related grid structures for SHDOM.

        Parameters
        ----------
        grid: xarray
            The grid containing x, y, z coordinates
        """

        def ibits(val, bit, ret_val):
            if val & 2 ** bit:
                return ret_val
            return 0

        self._ipflag = ipflag
        self._bcflag = bcflag

        self._pa = shdom.solver.ShdomPropertyArrays()

        # Set shdom property array
        self._pa.npx = grid.dims['x']
        self._pa.npy = grid.dims['y']
        self._pa.npz = grid.dims['z']
        self._pa.xstart = grid.x[0]
        self._pa.ystart = grid.y[0]

        # TODO: check constant dx and dy
        self._pa.delx = grid.x[1] - grid.x[0]
        self._pa.dely = grid.y[1] - grid.y[0]
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
        self._nbcells = (self._nz - 1) * (self._nx + ibits(self._bcflag, 0, 1) - \
                                          ibits(self._bcflag, 2, 1)) * (
                                    self._ny + ibits(self._bcflag, 1, 1) - ibits(self._bcflag, 3, 1))

        self._xgrid, self._ygrid, self._zgrid = shdom.core.new_grids(
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
        self._treeptr, self._cellflags = shdom.core.init_cell_structure(
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


    def carve(self,sensor_masks, agreement=None, linear_mode=False):
        #volume = np.zeros((self._nx, self._ny, self._nz), dtype=np.int)

        if isinstance(sensor_masks,xr.Dataset):
            sensor_masks = [sensor_masks]
        volume = np.zeros((len(sensor_masks), 3, self._nx, self._ny, self._nz))
        for i,  sensor in enumerate(sensor_masks):

            camx = sensor['ray_x'].data
            camy = sensor['ray_y'].data
            camz = sensor['ray_z'].data
            cammu = sensor['ray_mu'].data
            camphi = sensor['ray_phi'].data
            flags = sensor['cloud_mask'].data
            weights = sensor['weights'].data

            assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
            total_pix = camphi.size

            carved_volume = shdom.core.space_carve(
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
            volume[i] += carved_volume.reshape(3,self._nx, self._ny, self._nz)

        space_carved = xr.Dataset(
                        data_vars = {
                            'cloudy_counts': (['nsensors','x','y','z'], volume[:,0]),
                            'total_counts': (['nsensors','x','y','z'], np.sum(volume[:,:2],axis=1)),
                            'weights':(['nsensors','x','y','z'], volume[:,2]),
                        },
                        coords={'x':self._grid.x,
                                'y':self._grid.y,
                                'z':self._grid.z}
                        )

        if agreement is not None:
            mask = ((space_carved.cloudy_counts/space_carved.total_counts > agreement[0]).astype(np.int).mean('nsensors') >= agreement[1]).astype(np.int)
            masked_weights = (space_carved.weights/space_carved.total_counts).mean('nsensors')
            space_carved['mask'] = mask
            space_carved['masked_weights'] = masked_weights

        return space_carved
