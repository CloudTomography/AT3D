import sys
import xarray as xr
import numpy as np
from shdom import core, util

class RTE(object):
    """
    Radiative Trasnfer solver object.
    This object contains the interface to SHDOM internal structures and methods.

    Parameters
    ----------
    medium: list or xr.dataset
        A list or dataset containing the optical properties of the different scatter types within the medium
    numerical_params: xr.dataset
        A dataset containing the numerical parameters requiered for the RTE solution. These can be loaded
        from a config file (see ancillary_data/config.cfg).
    num_stokes: int, default=1
        The number of stokes for which to solve the RTE can be 1, 3, or 4.
        num_stokes=1 means unpolarized.
        num_stokes=3 means linear polarization.
        num_stokes=4 means full polarization.
    name: str, optional
       The name for the solver. Will be used when printing solution iteration messages.
       If non specified a default <type> <wavelength> is given, where <type> is Radiance for num_stokes=1 and
       Polarized for num_stokes>1

    Notes
    -----
    k-distribution not supported.
    """

    def __init__(self, numerical_params, medium, source, surface, num_stokes=1, name=None):

        # Check number of stokes and setup type of solver
        assert num_stokes in [1, 3, 4], 'num_stokes should be {1, 3, 4}'
        self._type = 'Radiance' if num_stokes == 1 else 'Polarization'
        self._nstokes = num_stokes
        self._nstleg = 1 if num_stokes == 1 else 6

        self.medium = [medium] if not isinstance(medium, list) else medium

        # Check that all optical scatterers have the same wavelength
        wavelengths = [scatterer.attrs['wavelength center'] for scatterer in self.medium if \
                       scatterer.attrs['wavelength center'] is not None]
        assert len(wavelengths) > 0, 'At least one scatterer has to have a wavelength defined'
        assert np.allclose(wavelengths[0], wavelengths), 'scatterers have different wavelengths {}'.format(wavelengths)
        self.wavelength = wavelengths[0]

        # For legacy purposes
        self._wavelen = util.float_round(self.wavelength)

        # Setup a name for the solver
        self._name = '{} {:1.3f} micron'.format(self._type, self.wavelength) if name is None else name

        # TODO: Start mpi (if available).
        self._masterproc = core.start_mpi()

        # Link to the properties array module.
        self._pa = ShdomPropertyArrays()

        # Setup numerical parameters
        self._setup_numerical_params(numerical_params)
        self.numerical_params = numerical_params

        # Setup source
        self._setup_source(source)
        self.source = source

        # Setup surface
        self._setup_surface(surface)
        self.surface = surface

        # k-distribution not supported yet
        self._kdist = False
        self._ng = 1
        self._delg = np.ones(1, order='F')
        self._pa.nzckd = 0
        self._baseout = False
        self._npart = 1

        # No iterations have taken place
        self._iters = 0

        # TODO: check that all scatterers are on the same grid
        self._grid = xr.Dataset({'x': self.medium[0].coords['x'],
                                 'y': self.medium[0].coords['y'],
                                 'z': self.medium[0].coords['z']})
        self._setup_grid(self._grid)

        # Temperature is used for thermal radiation. Not supported yet.
        self._pa.tempp = np.zeros(shape=(self._nbpts,), dtype=np.float32)

    def _setup_source(self, source):
        """
        TODO
        """
        # Source parameters
        # TODO: check that srctpye is supported
        self._srctype = source.srctype.data
        self._solarflux = source.solarflux.data
        self._solarmu = source.solarmu.data
        self._solaraz = source.solaraz.data
        self._skyrad = source.skyrad.data
        self._units = source.units.data
        self._waveno = source.wavenumber.data

    def _setup_surface(self, surface):
        """
        TODO
        """
        # Surface parameters
        # TODO: check that surface parameters are implementerd
        self._maxsfcpars = surface.maxsfcpars.data
        self._sfctype = str(surface.sfctype.data)
        self._gndalbedo = surface.gndalbedo.data
        self._gndtemp = surface.gndtemp.data
        self._nxsfc = surface.nxsfc.data
        self._nysfc = surface.nysfc.data
        self._delxsfc = surface.delxsfc.data
        self._delysfc = surface.delysfc.data
        self._nsfcpar = surface.nsfcpar.data
        self._sfcparms = surface.sfcparms.data
        self._sfcgridparms = surface.sfcgridparms.data

    def _setup_numerical_params(self, numerical_params):
        """
        Set the numerical parameters of the SHDOM forward solver.

        Parameters
        ----------
        numerical_params : NumericalParameters
            An object which encapsulate numerical parameters such as number of azimuthal and zenith bins.
        """
        # TODO: check input parameters (e.g. nmu should be even and grater than 2 nphi should be grater than 1)
        self._nmu = max(2, 2 * int((numerical_params.num_mu_bins.data + 1) / 2))
        self._nphi = max(1, numerical_params.num_phi_bins.data)
        self._deltam = numerical_params.deltam.data
        self._max_total_mb = numerical_params.max_total_mb.data
        self._adapt_grid_factor = float(numerical_params.adapt_grid_factor.data)
        self._accelflag = numerical_params.acceleration_flag.data
        self._num_sh_term_factor = numerical_params.num_sh_term_factor.data
        self._cell_to_point_ratio = numerical_params.cell_to_point_ratio.data
        self._splitacc = numerical_params.split_accuracy.data
        self._shacc = numerical_params.spherical_harmonics_accuracy.data
        self._solacc = numerical_params.solution_accuracy.data
        self._highorderrad = numerical_params.high_order_radiance.data

        # Setup horizontal boundary conditions
        # TODO: check that BC is strictly open or periodic
        self._bcflag = 0
        if numerical_params.x_boundary_condition == 'open':
            self._bcflag += 1
        if numerical_params.y_boundary_condition == 'open':
            self._bcflag += 2

        # Independent pixel not supported yet
        self._ipflag = numerical_params.ip_flag.data

    def _setup_grid(self, grid):
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

        self._xgrid, self._ygrid, self._zgrid = core.new_grids(
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

        self._setup_memory()

        self._npts, self._ncells, self._gridpos, self._gridptr, self._neighptr, \
        self._treeptr, self._cellflags = core.init_cell_structure(
            maxig=self._maxig,
            maxic=self._maxic,
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

    def _setup_memory(self):
        """A utility function to initialize internal memory parameters."""

        # Make ml and mm from nmu and nphi
        # ML is the maximum meridional mode, MM is the maximum azimuthal mode,
        # and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for
        # sines and cosines).
        # nphi0max: The maximum number of azimuth angles actually used;
        # for NCS=1 (cosine modes only) NPHI0=INT((NPHI+2)/2),
        # otherwise NPHI0=NPHI.
        self._ncs = 2
        self._ml = self._nmu - 1
        self._mm = max(0, int(self._nphi / 2) - 1)
        self._nlm = (2 * self._mm + 1) * (self._ml + 1) - self._mm * (self._mm + 1)
        if self._ncs == 1:
            self._nphi0max = int((self._nphi + 2) / 2)
        elif self._ncs == 2:
            self._nphi0max = self._nphi
        else:
            raise AttributeError
        self._memword = self._nmu * (2 + 2 * self._nphi + 2 * self._nlm + 2 * 33 * 32)

        # Guess maximum number of grid points, cells, SH vector size needed
        # but don't let MAX_TOTAL_MB be exceeded
        if self._max_total_mb * 1024 ** 2 > 1.75 * sys.maxsize:
            self._max_total_mb > 1.75 * sys.maxsize / 1024 ** 2
            print('MAX_TOTAL_MB reduced to fit memory model: ', self._max_total_mb)

        if self._splitacc < 0.0:
            self._adapt_grid_factor = 1.0

        self._big_arrays = 3
        if not self._accelflag:
            self._big_arrays = 2

        wantmem = self._adapt_grid_factor * self._nbpts * (
                28 + 16.5 * self._cell_to_point_ratio + \
                self._nphi0max * self._nstokes + self._num_sh_term_factor * self._nstokes * self._nlm * self._big_arrays)

        REDUCE = min(1.0, ((self._max_total_mb * 1024 ** 2) / 4 - self._memword) / wantmem)
        self._adapt_grid_factor *= REDUCE
        assert self._adapt_grid_factor > 1.0, 'max_total_mb memory limit exceeded with just base grid.'

        if REDUCE < 1.0:
            print('adapt_grid_factor reduced to ', self._adapt_grid_factor)

        wantmem *= REDUCE
        assert wantmem < sys.maxsize, 'number of words of memory exceeds max integer size: %d' % wantmem
        self._maxig = int(self._adapt_grid_factor * self._nbpts)
        self._maxic = int(self._cell_to_point_ratio * self._maxig)
        self._maxiv = int(self._num_sh_term_factor * self._nlm * self._maxig)
        self._maxido = self._maxig * self._nphi0max

        assert 4.0 * (
                    self._maxiv + self._maxig) * self._nstokes < sys.maxsize, 'size of big sh arrays (maxiv) probably exceeds max integer number of bytes: %d' % self._maxiv
        assert 4.0 * 8.0 * self._maxic <= sys.maxsize, 'size of gridptr array (8*maxic) probably exceeds max integer number of bytes: %d' % 8 * self._maxic

        self._maxnbc = int(self._maxig * 3 / self._nz)
        self._maxsfcpars = 4
        if self._sfctype.endswith('L'):
            self._maxbcrad = 2 * self._maxnbc
        else:
            self._maxbcrad = int((2 + self._nmu * self._nphi0max / 2) * self._maxnbc)

    def _init_solution(self):
        """
        TODO: improve this
        Initilize the solution (I, J fields) from the direct transmission and a simple layered model.
        Many arrays are initialized including the dirflux (the direct solar transmission).
        """

        # Iterate over all particle types and aggragate scattering table
        self._pa.extinctp = np.zeros(shape=[self._nbpts, len(self.medium)], dtype=np.float32)
        self._pa.albedop = np.zeros(shape=[self._nbpts, len(self.medium)], dtype=np.float32)
        self._pa.iphasep = np.zeros(shape=[self._nbpts, len(self.medium)], dtype=np.int32)
        for i, scatterer in enumerate(self.medium):
            self._pa.extinctp[:, i] = scatterer.extinction.data.ravel()
            self._pa.albedop[:, i] = scatterer.ssalb.data.ravel()
            self._pa.iphasep[:, i] = scatterer.table_index.data.ravel() + self._pa.iphasep.max()

        #TODO 0 indices may appear in the clear sky for the first optical scatterer.
        #set them to 1. This should be revisited after all checks are developed for the
        #workflow that creates table_index.
        self._pa.iphasep[np.where(self._pa.iphasep == 0)] = 1
        # Concatenate all scatterer tables into one table
        legendre_table = xr.concat([scatterer.legcoef for scatterer in self.medium], dim='table_index')

        self._pa.numphase = legendre_table.sizes['table_index']

        # Determine the number of legendre coefficient for a given angular resolution
        #Note this is corrected to self._ml rather than self._mm based on original SHDOM and
        #what PREPARE_PROP assumes.
        self._nleg = self._ml + 1 if self._deltam else self._ml

        self._nleg = max(legendre_table.sizes['legendre_index'] - 1, self._nleg)
        self._nscatangle = max(36, min(721, 2 * self._nleg))

        # Check if legendre table needs padding

        if self._nleg > legendre_table.sizes['legendre_index']:
            legendre_table = legendre_table.pad({'legendre_index': (0, 1 + self._nleg - legendre_table.sizes['legendre_index'])},
                               constant_values=0.0)

        # Check if scalar or vector RTE
        if self._nstokes == 1:
            self._pa.legenp = legendre_table.isel(stokes_index=0).data.ravel(order='F').astype(np.float32)
        else:
            self._pa.legenp = legendre_table.data.ravel(order='F').astype(np.float32)

        self._maxasym = np.max(legendre_table.isel(stokes_index=0, legendre_index=1).data / 3.0)
        self._maxpgl = self._maxpg * legendre_table.sizes['legendre_index']

        if legendre_table.sizes['table_index'] > 0:
            self._maxigl = legendre_table.sizes['table_index'] * (legendre_table.sizes['legendre_index'] + 1)
        else:
            self._maxigl = self._maxig * (legendre_table.sizes['legendre_index'] + 1)

        self._npart = len(self.medium)
        self._nstphase = min(self._nstleg, 2)

        # Transfer property arrays into internal grid structures
        self._temp, self._planck, self._extinct, self._albedo, self._legen, self._iphase, \
        self._total_ext, self._extmin, self._scatmin, self._albmax = core.transfer_pa_to_grid(
            nstleg=self._nstleg,
            npart=self._npart,
            extinctp=self._pa.extinctp,
            albedop=self._pa.albedop,
            iphasep=self._pa.iphasep,
            delx=self._pa.delx,
            dely=self._pa.dely,
            xstart=self._pa.xstart,
            ystart=self._pa.ystart,
            zlevels=self._pa.zlevels,
            tempp=self._pa.tempp,
            legenp=self._pa.legenp,
            nzckd=self._pa.nzckd,
            zckd=self._pa.zckd,
            gasabs=self._pa.gasabs,
            ml=self._ml,
            mm=self._mm,
            numphase=self._pa.numphase,
            deltam=self._deltam,
            units=self._units,
            waveno=self._waveno,
            wavelen=self._wavelen,
            gridpos=self._gridpos,
            nleg=self._nleg,
            maxig=self._maxig,
            npx=self._pa.npx,
            npy=self._pa.npy,
            npz=self._pa.npz,
            srctype=self._srctype,
            npts=self._npts)

        # Restart solution criteria
        self._oldnpts = 0
        self._solcrit = 1.0
        self._nang, self._nphi0, self._mu, self._phi, self._wtdo, \
        self._sfcgridparms, self._ntoppts, self._nbotpts, self._bcptr, self._rshptr, self._shptr, self._oshptr, \
        self._source, self._delsource, self._radiance, self._fluxes, self._dirflux, self._ylmsun, self._uniformzlev, \
        self._pa.extdirp, self._bcrad, self._cx, self._cy, self._cz, self._cxinv, self._cyinv, self._czinv, \
        self._ipdirect, self._di, self._dj, self._dk, self._epss, self._epsz, self._xdomain, self._ydomain, \
        self._delxd, self._delyd, self._deljdot, self._deljold, self._deljnew, self._jnorm, self._fftflag, \
        self._cmu1, self._cmu2, self._wtmu, self._cphi1, self._cphi2, self._wphisave, self._work, self._work1, \
        self._work2, self._uniform_sfc_brdf, self._sfc_brdf_do = core.init_solution(
            nstleg=self._nstleg,
            nstokes=self._nstokes,
            nx=self._nx,
            ny=self._ny,
            nz=self._nz,
            nx1=self._nx1,
            ny1=self._ny1,
            ml=self._ml,
            mm=self._mm,
            maxiv=self._maxiv,
            maxic=self._maxic,
            maxig=self._maxig,
            maxido=self._maxido,
            shacc=self._shacc,
            ncs=self._ncs,
            nlm=self._nlm,
            nmu=self._nmu,
            nphi=self._nphi,
            nleg=self._nleg,
            numphase=self._pa.numphase,
            bcflag=self._bcflag,
            ipflag=self._ipflag,
            deltam=self._deltam,
            srctype=self._srctype,
            solarflux=self._solarflux,
            solarmu=self._solarmu,
            solaraz=self._solaraz,
            skyrad=self._skyrad,
            sfctype=self._sfctype,
            gndtemp=self._gndtemp,
            gndalbedo=self._gndalbedo,
            nxsfc=self._nxsfc,
            nysfc=self._nysfc,
            delxsfc=self._delxsfc,
            delysfc=self._delysfc,
            nsfcpar=self._nsfcpar,
            sfcparms=self._sfcparms,
            units=self._units,
            waveno=self._waveno,
            wavelen=self.wavelength,
            accelflag=self._accelflag,
            xgrid=self._xgrid,
            ygrid=self._ygrid,
            zgrid=self._zgrid,
            temp=self._temp,
            planck=self._planck,
            extinct=self._extinct,
            albedo=self._albedo,
            legen=self._legen,
            iphase=self._iphase,
            maxnbc=self._maxnbc,
            maxbcrad=self._maxbcrad,
            npts=self._npts,
            gridpos=self._gridpos,
            ncells=self._ncells,
            gridptr=self._gridptr,
            treeptr=self._treeptr,
            npx=self._pa.npx,
            npy=self._pa.npy,
            npz=self._pa.npz,
            delx=self._pa.delx,
            dely=self._pa.dely,
            xstart=self._pa.xstart,
            ystart=self._pa.ystart,
            zlevels=self._pa.zlevels,
            tempp=self._pa.tempp,
            extinctp=self._pa.extinctp,
            nbpts=self._nbpts,
            albedop=self._pa.albedop,
            legenp=self._pa.legenp,
            iphasep=self._pa.iphasep,
            nzckd=self._pa.nzckd,
            zckd=self._pa.zckd,
            gasabs=self._pa.gasabs,
            npart=self._npart,
            total_ext=self._total_ext,
            nbcells=self._nbcells,
            maxsfcpars=self._maxsfcpars,
            nphi0max=self._nphi0max
        )

    def _precompute_phase(self):
        """
        Precompute angular scattering for the entire legendre table.
        Preform a negativity check. (negcheck=True).
        """
        self._phasetab = core.precompute_phase_check(
            negcheck=True,
            nscatangle=self._nscatangle,
            numphase=self._pa.numphase,
            nstphase=self._nstphase,
            nstokes=self._nstokes,
            nstleg=self._nstleg,
            nleg=self._nleg,
            ml=self._ml,
            nlm=self._nlm,
            legen=self._legen,
            deltam=self._deltam
        )

    def _make_direct(self):
        """
        Compute the direct transmission from the solar beam.
        """
        self._dirflux, self._pa.extdirp, self._cx, self._cy, self._cz, self._cxinv, self._cyinv, self._czinv, \
        self._di, self._dj, self._dk, self._ipdirect, self._delxd, self._delyd, self._xdomain, self._ydomain, \
        self._epss, self._epsz, self._uniformzlev = \
            core.make_direct(
                nstleg=self._nstleg,
                dirflux=self._dirflux,
                extdirp=self._pa.extdirp,
                npts=self._npts,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                ml=self._ml,
                nleg=self._nleg,
                solarflux=self._solarflux,
                solarmu=self._solarmu,
                solaraz=self._solaraz,
                gridpos=self._gridpos,
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,
                numphase=self._pa.numphase,
                delx=self._pa.delx,
                dely=self._pa.dely,
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,
                tempp=self._pa.tempp,
                extinctp=self._pa.extinctp,
                albedop=self._pa.albedop,
                legenp=self._pa.legenp,
                iphasep=self._pa.iphasep,
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs
            )

    def solve(self, maxiter, init_solution=True, verbose=True):
        """
        Main solver routine. This routine is comprised of two parts:
          1. Initialization, optional
          2. Solution iterations

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        init_solution: boolean, default=True
            If False the solution is initialized according to the existing radiance and source function saved within the RteSolver object (previously computed)
            If True or no prior solution (I,J fields) exists then an initialization is preformed (part 1.).
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        # Part 1: Initialize solution (from a layered model)
        if init_solution:
            self._init_solution()

        # Part 2: Solution itertaions
        self._sfcgridparms, self._solcrit, \
            self._iters, self._temp, self._planck, self._extinct, self._albedo, self._legen, self._iphase, \
            self._ntoppts, self._nbotpts, self._bcptr, self._bcrad, self._npts, self._gridpos, self._ncells, self._gridptr, \
            self._neighptr, self._treeptr, self._cellflags, self._rshptr, self._shptr, self._oshptr, self._source, \
            self._delsource, self._radiance, self._fluxes, self._dirflux, self._uniformzlev, self._pa.extdirp, \
            self._oldnpts, self._total_ext, self._deljdot, self._deljold, self._deljnew, self._jnorm, \
            self._work, self._work1, self._work2  = core.solution_iterations(
                uniform_sfc_brdf=self._uniform_sfc_brdf,
                sfc_brdf_do=self._sfc_brdf_do,
                work=self._work,
                work1=self._work1,
                work2=self._work2,
                bcrad=self._bcrad,
                fluxes=self._fluxes,
                nang=self._nang,
                nphi0=self._nphi0,
                maxnbc=self._maxnbc,
                ntoppts=self._ntoppts,
                nbotpts=self._nbotpts,
                uniformzlev=self._uniformzlev,
                extmin=self._extmin,
                scatmin=self._scatmin,
                cx=self._cx,
                cy=self._cy,
                cz=self._cz,
                cxinv=self._cxinv,
                cyinv=self._cyinv,
                czinv=self._czinv,
                ipdirect=self._ipdirect,
                di=self._di,
                dj=self._dj,
                dk=self._dk,
                nphi0max=self._nphi0max,
                epss=self._epss,
                epsz=self._epsz,
                xdomain=self._xdomain,
                ydomain=self._ydomain,
                delxd=self._delxd,
                delyd=self._delyd,
                albmax=self._albmax,
                deljdot=self._deljdot,
                deljold=self._deljold,
                deljnew=self._deljnew,
                jnorm=self._jnorm,
                fftflag=self._fftflag,
                cmu1=self._cmu1,
                cmu2=self._cmu2,
                wtmu=self._wtmu,
                cphi1=self._cphi1,
                cphi2=self._cphi2,
                wphisave=self._wphisave,
                nbpts=self._nbpts,
                npart=self._npart,
                extinct=self._extinct,
                albedo=self._albedo,
                legen=self._legen,
                total_ext=self._total_ext,
                extinctp=self._pa.extinctp,
                albedop=self._pa.albedop,
                planck=self._planck,
                iphase=self._iphase,
                iphasep=self._pa.iphasep,
                nstokes=self._nstokes,
                nstleg=self._nstleg,
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,
                delx=self._pa.delx,
                dely=self._pa.dely,
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,
                tempp=self._pa.tempp,
                legenp=self._pa.legenp,
                extdirp=self._pa.extdirp,
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs,
                solcrit=self._solcrit,
                nx=self._nx,
                ny=self._ny,
                nx1=self._nx1,
                ny1=self._ny1,
                nz=self._nz,
                ml=self._ml,
                mm=self._mm,
                ncs=self._ncs,
                nlm=self._nlm,
                nmu=self._nmu,
                nphi=self._nphi,
                numphase=self._pa.numphase,
                mu=self._mu,
                phi=self._phi,
                wtdo=self._wtdo,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                srctype=self._srctype,
                highorderrad=self._highorderrad,
                solarflux=self._solarflux,
                solarmu=self._solarmu,
                solaraz=self._solaraz,
                skyrad=self._skyrad,
                sfctype=self._sfctype,
                gndtemp=self._gndtemp,
                gndalbedo=self._gndalbedo,
                nxsfc=self._nxsfc,
                nysfc=self._nysfc,
                delxsfc=self._delxsfc,
                delysfc=self._delysfc,
                nsfcpar=self._nsfcpar,
                sfcparms=self._sfcparms,
                sfcgridparms=self._sfcgridparms,
                units=self._units,
                waveno=self._waveno,
                wavelen=self._wavelen,
                accelflag=self._accelflag,
                solacc=self._solacc,
                maxiter=maxiter,
                splitacc=self._splitacc,
                shacc=self._shacc,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                temp=self._temp,
                maxbcrad=self._maxbcrad,
                bcptr=self._bcptr,
                npts=self._npts,
                gridpos=self._gridpos,
                ncells=self._ncells,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags,
                rshptr=self._rshptr,
                shptr=self._shptr,
                oshptr=self._oshptr,
                source=self._source,
                delsource=self._delsource,
                radiance=self._radiance,
                dirflux=self._dirflux,
                nleg=self._nleg,
                maxiv=self._maxiv,
                maxic=self._maxic,
                maxig=self._maxig,
                maxido=self._maxido,
                verbose=verbose,
                oldnpts=self._oldnpts,
                ylmsun=self._ylmsun,
                runname=self._name
            )

    def integrate_to_sensor(self, sensor, n_jobs=1):
        """
        TODO
        Integrates the source function along rays with geometry specified in sensor.
        """

        camx = sensor['cam_x'].data
        camy = sensor['cam_y'].data
        camz = sensor['cam_z'].data
        cammu = sensor['cam_mu'].data
        camphi = sensor['cam_phi'].data
        #TODO
        #Some checks on the dimensions: this kind of thing.
        assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
        total_pix = sensor.sizes['total_pixels']

        phase_check = hasattr(self, '_phasetab')
        if not phase_check:
            self._precompute_phase()

        #TODO Raise warning if observables are not included.

        #TODO figure out parallelization.

        output = core.render(
            nstphase=self._nstphase,
            ylmsun=self._ylmsun,
            phasetab=self._phasetab,
            nscatangle=self._nscatangle,
            ncs=self._ncs,
            nstokes=self._nstokes,
            nstleg=self._nstleg,
            camx=camx,
            camy=camy,
            camz=camz,
            cammu=cammu,
            camphi=camphi,
            npix=total_pix,
            nx=self._nx,
            ny=self._ny,
            nz=self._nz,
            bcflag=self._bcflag,
            ipflag=self._ipflag,
            npts=self._npts,
            ncells=self._ncells,
            ml=self._ml,
            mm=self._mm,
            nlm=self._nlm,
            numphase=self._pa.numphase,
            nmu=self._nmu,
            nphi0max=self._nphi0max,
            nphi0=self._nphi0,
            maxnbc=self._maxnbc,
            ntoppts=self._ntoppts,
            nbotpts=self._nbotpts,
            nsfcpar=self._nsfcpar,
            gridptr=self._gridptr,
            neighptr=self._neighptr,
            treeptr=self._treeptr,
            shptr=self._shptr,
            bcptr=self._bcptr,
            cellflags=self._cellflags,
            iphase=self._iphase[:self._npts],
            deltam=self._deltam,
            solarmu=self._solarmu,
            solaraz=self._solaraz,
            gndtemp=self._gndtemp,
            gndalbedo=self._gndalbedo,
            skyrad=self._skyrad,
            waveno=self._waveno,
            wavelen=self._wavelen,
            mu=self._mu,
            phi=self._phi.reshape(self._nmu, -1),
            wtdo=self._wtdo.reshape(self._nmu, -1),
            xgrid=self._xgrid,
            ygrid=self._ygrid,
            zgrid=self._zgrid,
            gridpos=self._gridpos,
            sfcgridparms=self._sfcgridparms,
            bcrad=self._bcrad,
            extinct=self._extinct[:self._npts],
            albedo=self._albedo[:self._npts],
            legen=self._legen,
            dirflux=self._dirflux,
            fluxes=self._fluxes,
            source=self._source,
            srctype=self._srctype,
            sfctype=self._sfctype,
            units=self._units,
            total_ext=self._total_ext[:self._npts],
            npart=self._npart)

        #merge output across parallelization pixel split before doing this.
        sensor['I'] = xr.DataArray(
            data=output[0],
            dims='total_pixels',
            attrs={
                'long_name': 'Radiance'
            }
        )
        if self._nstokes > 1:
            sensor['Q']= xr.DataArray(
                data=output[1],
                dims='total_pixels',
                attrs={
                    'long_name': 'Stokes Parameter for Linear Polarization (Q)'
                }
            )
            sensor['U']= xr.DataArray(
                data=output[2],
                dims='total_pixels',
                attrs={
                    'long_name': 'Stokes Parameter for Linear Polarization (U)'
                }
            )
        if self._nstokes == 4:
            sensor['V']= xr.DataArray(
                data=output[3],
                dims='total_pixels',
                attrs={
                    'long_name': 'Stokes Paramaeter for Circular Polarization (V)'
                }
            )

        return output, sensor

    @property
    def num_iterations(self):
        return self._iters


class ShdomPropertyArrays(object):
    """
    Shdom property array module.
    Contain the parameters that are in the original SHDOM_PROPERTY_ARRAYS module
    in shdom90.f90 .

    Parameters
    ----------
    npx: Number of x grid points
    npy: Number of y grid points
    npz: Number of z grid points
    numphase: Number of phase function enteries in the table
    delx: Delta-x spacing
    dely: Delta-y spacing
    xstart: Starting position of x coordinates
    ystart: Starting position of y coordinates
    zlevels: Altitude grid points
    tempp: Temperatures at altitude grid points
    extinctp: Extinction grid
    albedop: Single scattering albedo grid
    legenp: legendre phase function table
    extdirp: optical depth integrated along sun rays
    iphasep: pointer to phase function table
    nzckd: Number of correlated-k distribution points
    zckd: correlated-k distribution grid
    gasabs: Gas absorption
    """
    def __init__(self):
        self.npx = None
        self.npy = None
        self.npz = None
        self.numphase = None
        self.delx = None
        self.dely = None
        self.xstart = None
        self.ystart = None
        self.zlevels = None
        self.tempp = None
        self.extinctp = None
        self.albedop = None
        self.legenp = None
        self.extdirp = None
        self.iphasep = None
        self.nzckd = None
        self.zckd = None
        self.gasabs = None
