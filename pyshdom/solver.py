"""
This module defines an `RTE` object which sets up an SHDOM solution to the
RTE equation through calls to fortran subroutines. See src/polarized/shdomsub1.f etc.
xr.Datasets prepared from sensor, source, surface, medium/mie/size_distribution are
used here.
containers.py contains a `SolversDict` object which stores multiple solver.RTE
objects and can be used for parallelization of solving the RTE etc.
"""

import sys
import warnings
import copy
import typing
from collections import OrderedDict
import psutil

import xarray as xr
import numpy as np

import pyshdom.core
import pyshdom.util
import pyshdom.checks


class ShdomPropertyArrays(object):
    """
    Shdom property array module.
    Contains the parameters that were in the original SHDOM_PROPERTY_ARRAYS module
    in shdom90.f90.

    This contains arrays that have been formatted for input to SHDOM routines
    but have not yet been fully preprocessed through delta scaling etc. In the
    original SHDOM, the arrays here are formed by PROPGEN and associated routines.
    In our case, the data is simply rearranged from python inputs as those are
    already defined on 'SHDOM-like' grids (see grid.py).

    Parameters
    ----------
    npx: int
        Number of x grid points
    npy: int
        Number of y grid points
    npz: int
        Number of z grid points
    numphase: int
        Number of phase function enteries in the table
    delx: float32
        Delta-x spacing
    dely: float32
        Delta-y spacing
    xstart: float32
        Starting position of x coordinates
    ystart: float32
        Starting position of y coordinates
    zlevels: np.ndarray(shape=(npz,), dtype=np.float32))
        Altitude grid points
    tempp: np.ndarray
        Temperatures at grid points
    extinctp: np.ndarray
        Extinction on grid
    albedop: np.ndarray
        Single scattering albedo grid
    legenp: np.ndarray
        legendre phase function table
    extdirp: np.ndarray
        Delta-M scaled extinction on grid, used for calculating direct beam.
    iphasep: np.ndarray, dtype=np.int32
        pointer to phase function table
    nzckd: int
        Number of correlated-k distribution points
    zckd: np.ndarray
        correlated-k distribution z levels.
    gasabs: np.ndarray, shape=(npz,)
        Gas absorption extinction on vertical levels.
    nlegp: int
        Number of legendre moments in full phase functio expansions (legenp)
    max_num_micro: int
        The maximum number of phase functions to mix per grid point across
        all scattering species.
    phasewtp: np.ndarray
        interpolation weights for mixtures of phase functions.
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
        self.nlegp = None
        self.max_num_micro = None
        self.phasewtp = None

class RTE:
    """
    Radiative Transfer solver object.
    This object contains the interface to SHDOM's internal structures and methods.

    Rather than reading inputs from files, the xr.Datasets are directly passed to
    instantiate this class.
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
    def __init__(self, numerical_params, medium, source, surface, num_stokes=1, name=None,
                 atmosphere=None):

        # Check number of stokes and setup type of solver
        if num_stokes not in (1, 3, 4):
            raise ValueError("num_stokes should be in (1, 3, 4) not '{}'".format(num_stokes))
        self._type = 'Radiance' if num_stokes == 1 else 'Polarization'
        self._nstokes = num_stokes
        self._nstleg = 1 if num_stokes == 1 else 6

        # If the second character is 'O' then 'max scattering' interpolation
        # of phase functions occurs. If the second character is 'N' then
        # linear mixing of phase functions occurs.
        self._interpmethod = 'ON'
        # phasemax is only used for linear mixing of phase functions
        # and is the threshold for the weight for neglecting the
        # contributions of other phase functions.
        self._phasemax = 0.999
        self._adjflag = False # Used for testing the adjoint source.

        self._correctinterpolate = True
        # Flag for switching between the 'Radiance' and 'Visualization'
        # methods for calculating radiances in original SHDOM.
        # This is added for verification against original SHDOM.
        # The methods only differ in calculations when line integrations
        # pass through a region with varying grid resolution.
        # `True` corresponds to the 'Visualization' method which utilizes
        # a more accurate interpolation for the source-extinction
        # product than the 'Radiance' method.

        self.source = self._setup_source(source)
        self.medium, self._grid = self._setup_medium(medium)

        # Setup a name for the solver
        self._name = '{} {:1.3f} micron'.format(self._type, self.wavelength) if name is None else name
        #Start mpi (if available). This is a dummy routine. MPI is not currently supported.
        self._masterproc = pyshdom.core.start_mpi()

        # Link to the properties array module.
        self._pa = ShdomPropertyArrays()

        # k-distribution not supported yet
        self._kdist = False
        self._ng = 1
        self._delg = np.ones(1, order='F')
        self._pa.nzckd = 0
        self._baseout = False
        self._npart = 1

        # No iterations have taken place
        self._iters = 0

        self.numerical_params = self._setup_numerical_params(numerical_params)
        self._setup_grid(self._grid)
        self.surface = self._setup_surface(surface)
        # atmosphere includes temperature for thermal radiation
        self.atmosphere = self._setup_atmosphere(atmosphere)

        #these warnings can only be done after both surface and numerical params.
        flux = 0.0
        if self._srctype in ('T', 'B'):
            #Surface flux is typically the warmest and therefore largest flux
            #so it is used as guidance for the splitting accuracy.
            flux = np.pi*pyshdom.util.planck_function(self._gndtemp, self.wavelength)
        if self._srctype != 'T':
            flux += self._solarflux
        if (self._splitacc > 0.0) & ((self._splitacc < 0.001*flux) | (self._splitacc > 0.1*flux)):
            warnings.warn("Splitting accuracy parameter is not a fraction but scales with fluxes."
                          " splitacc/nominal flux = {}".format(self._splitacc/flux))
        if (self._shacc > 0.0) & (self._shacc > 0.03*flux):
            warnings.warn("spherical_harmonics_accuracy is not a fraction but scales with fluxes."
                          " SHACC/nominal flux = {}".format(self._shacc/flux))

        #this is called at initialization so that warnings about the optical
        #thickness across cells in the medium can be called to warn a user
        #before they try to run RTE.solve().
        self._prepare_optical_properties()

        #here is where an xr.Dataset containing a solved solution will go
        #if it is loaded (self.load_solution.) in preparation to be read into
        #memory for the self.solve method. (in self._init_solution.)
        self._restore_data = None
        self._setup_grid_flag = True

        #set the cached spherical_harmonics and net flux divergence to None
        self._netfluxdiv = None
        self._shterms = None

        #Initialize solution criterion here so that we can use it to check
        #if RTE is 'solved'.
        self._solcrit = None

        #initialize these attributes that are only filled with the
        #true memory-related numerical parameters after an SHDOM solution.
        self._maxmb_out = None
        self._adapt_grid_factor_out = None
        self._shterm_fac_out = None
        self._cell_point_out = None


    @property
    def final_maxmb(self):
        return self._maxmb_out

    @property
    def final_adapt_grid_factor(self):
        return self._adapt_grid_factor_out

    @property
    def final_shterm_factor(self):
        return self._shterm_fac_out

    @property
    def final_cell_point_ratio(self):
        return self._cell_point_out

    @property
    def solution_accuracy(self):
        return self._solacc

    def set_solution_accuracy(self, val):
        """
        Update the solution accuracy to allow a more accurate solution to
        be iterated towards without reinitializing the solver.
        """
        self._solacc = val
        self.numerical_params['solacc'] = val
        #important to update both consistently.

    def solve(self, maxiter, init_solution=True, setup_grid=True, verbose=True,
              solve=True):
        """
        Main solver routine. This routine is comprised of two parts:
          1. Initialization, optional
          2. Solution iterations

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution to SHDOM.
        setup_grid: boolean
            If True then a new grid is initialized. If False then the grid
            (including adaptive grid points)
        init_solution: boolean, default=True
            If False then a solution is initialized. This is overwritten
            to True if no existing Radiance/Source function fields exist.
            The solution initialization will use a Radiance/Source field provided
            by RTE.load_solution or use a 1D two-stream model.
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        if not isinstance(verbose, np.bool):
            raise TypeError("`verbose` should be a boolean.")
        if not isinstance(init_solution, np.bool):
            raise TypeError("`init_solution` should be a boolean.")
        if not isinstance(setup_grid, np.bool):
            raise TypeError("`setup_grid` should be a boolean.")

        # Part 1: Initialize solution (from a 1D layered model)
        # or from a loaded solution if available.
        if np.any([not hasattr(self, attr) for attr in ('_radiance', '_source', '_fluxes')]) \
        and (not init_solution):
            warnings.warn(
                "RTE object does not have initialized Radiance/Source fields as such, the"
                " `init_solution` flag has been overwritten and a solution will be initialized.")
            init_solution = True
        if init_solution:
            self._init_solution(
                setup_grid=setup_grid,
            )
        if not self.check_solved(verbose=False):
            if maxiter <= self._iters:
                warnings.warn(
                    "The solver is not converged to the specified accuracy but maxiter "
                    "has already been exceeded. Please increase `maxiter`."
                )

        #set the cached spherical_harmonics and net flux divergence to None
        #as a new solution has been formed. This is also done in ._init_solution
        #AND here as either could be done without the other and either way the
        #cached values are no longer representative.
        self._netfluxdiv = None
        self._shterms = None

        # Part 2: Solution itertaions
        # This is the time consuming part, equivalent to SOLVE_RTE in SHDOM.
        # All of these arrays are initialized in _init_solution or _setup_grid.
        # And are modified in-place to reflect the solved RTE.
        self._sfcgridparms, self._solcrit, self._iters, self._temp, self._planck, \
        self._extinct, self._albedo, self._legen, self._iphase, self._ntoppts, \
        self._nbotpts, self._bcptr, self._bcrad, self._npts, self._gridpos, \
        self._ncells, self._gridptr, self._neighptr, self._treeptr, self._cellflags, \
        self._rshptr, self._shptr, self._oshptr, self._source, self._delsource, \
        self._radiance, self._fluxes, self._dirflux, self._uniformzlev, \
        self._pa.extdirp, self._oldnpts, self._total_ext, self._deljdot, \
        self._deljold, self._deljnew, self._jnorm, self._work, self._work1, \
        self._work2, ierr, errmsg, self._phaseinterpwt, self._cpu_time \
         = pyshdom.core.solution_iterations(
            verbose=verbose,
            solve=solve,
            maxnmicro=self._pa.max_num_micro,
            phasewtp=self._pa.phasewtp,
            nlegp=self._pa.nlegp,
            phasemax=self._phasemax,
            phaseinterpwt=self._phaseinterpwt,
            interpmethod=self._interpmethod,
            iterfixsh=self._iterfixsh,
            iter=self._iters,
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
            wavelen=self.wavelength,
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
            oldnpts=self._oldnpts,
            ylmsun=self._ylmsun,
            runname=self._name
        )
        pyshdom.checks.check_errcode(ierr, errmsg)
        nsh = self._shptr[self._npts]
        self._maxmb_out = 4*(self._nmu*(2+2*self._nphi + 2*self._nlm+2*33*32) \
                        + 4.5*self._maxpg + self._maxpgl + self._nstleg*self._pa.numphase*(self._nleg + 1) \
                        + 16.5*self._ncells + self._npts*(28+self._nphi0max*self._nstokes) \
                        + self._nstokes*nsh*self._big_arrays)/(1024**2)
        self._adapt_grid_factor_out = self._npts/self._nbpts
        self._shterm_fac_out = nsh/(self._nlm*self._npts)
        self._cell_point_out = self._ncells/self._npts
        # if verbose:
        #     warnings.warn("Actual MAX_TOTAL_MB: {:.2f}".format(self._maxmb_out))
        #     warnings.warn("Actual adapt_grid_factor: {:.4f}".format(self._adapt_grid_factor_out))
        #     warnings.warn("Actual cell_point_ratio: {:.4f}".format(self._cell_point_out))

    def integrate_to_sensor(self, sensor, single_scatter=False):
        """Calculates the StokesVector at specified geometry using an RTE solution.

        Integrates the source function along rays with positions and
        directions specified in sensor. This is the method SHDOM uses to calculate
        Radiances. Each 'ray' StokesVector is as close to an idealized (delta function)
        sampling of the radiance field as the discretization used by SHDOM allows.
        As such the ray values are not area averaged.

        Parameters
        ----------
        sensor : xr.Dataset
            A valid pyshdom sensor dataset (see sensor.py) that contains AT LEAST
            the ray geometries required to perform the Source function integration.

        Returns
        -------
        sensor : xr.Dataset
            The same sensor that was input but modified in-place by the addition
            of the simulated Stokes Vector.
        """
        if not isinstance(sensor, xr.Dataset):
            raise TypeError("`sensor` should be an xr.Dataset not "
                            "of type '{}''".format(type(sensor)))
        pyshdom.checks.check_hasdim(sensor, ray_mu='nrays', ray_phi='nrays',
                                    ray_x='nrays', ray_y='nrays', ray_z='nrays',
                                    stokes='stokes_index')

        if 'nimage' in sensor.stokes.dims:
            stokes_averaged = sensor.stokes.any('nimage')
        else:
            stokes_averaged = sensor.stokes
        if stokes_averaged.sum('stokes_index') > self._nstokes:
            raise ValueError("'{}' Stokes components are required by sensor but RTE "
                             "only has nstokes={}".format(stokes_averaged.data,
                                                          self._nstokes)
                            )

        camx = sensor['ray_x'].data
        camy = sensor['ray_y'].data
        camz = sensor['ray_z'].data
        cammu = sensor['ray_mu'].data
        camphi = sensor['ray_phi'].data
        total_pix = sensor.sizes['nrays']

        self.check_solved()
        self._precompute_phase()

        output, ierr, errmsg = pyshdom.core.render(
            correctinterpolate=self._correctinterpolate,
            singlescatter=single_scatter,
            transcut=self._transcut,
            tautol=self._tautol,
            maxnmicro=self._pa.max_num_micro,
            interpmethod=self._interpmethod,
            phaseinterpwt=self._phaseinterpwt[:,:self._npts,:],
            phasemax=self._phasemax,
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
            iphase=self._iphase[:,:self._npts],
            deltam=self._deltam,
            solarmu=self._solarmu,
            solaraz=self._solaraz,
            gndtemp=self._gndtemp,
            gndalbedo=self._gndalbedo,
            skyrad=self._skyrad,
            waveno=self._waveno,
            wavelen=self.wavelength,
            mu=self._mu,
            phi=self._phi.reshape(self._nmu, -1),
            wtdo=self._wtdo.reshape(self._nmu, -1),
            xgrid=self._xgrid,
            ygrid=self._ygrid,
            zgrid=self._zgrid,
            gridpos=self._gridpos,
            sfcgridparms=self._sfcgridparms,
            bcrad=copy.deepcopy(self._bcrad), #deep copied as it is modified in place
                                              #which is otherwise bad for parallelization.
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
        pyshdom.checks.check_errcode(ierr, errmsg)
        sensor['I'] = xr.DataArray(
            data=output[0],
            dims='nrays',
            attrs={
                'long_name': 'Radiance'
            }
        )
        if self._nstokes > 1:
            sensor['Q'] = xr.DataArray(
                data=output[1],
                dims='nrays',
                attrs={
                    'long_name': 'Stokes Parameter for Linear Polarization (Q)'
                }
            )
            sensor['U'] = xr.DataArray(
                data=output[2],
                dims='nrays',
                attrs={
                    'long_name': 'Stokes Parameter for Linear Polarization (U)'
                }
            )
        if self._nstokes == 4:
            sensor['V'] = xr.DataArray(
                data=output[3],
                dims='nrays',
                attrs={
                    'long_name': 'Stokes Parameter for Circular Polarization (V)'
                }
            )
        return sensor

    def optical_path(self, sensor, deltam_scaled_path=False):
        """Calculates the optical paths along specified rays by integrating
        the extinction field.

        Parameters
        ----------
        sensor : xr.Dataset
            A valid pyshdom sensor dataset (see sensor.py) that contains AT LEAST
            the ray geometries required to define the line integration of extinction.
        deltam_scaled_path : bool
            If True then the optical path is calculated for the delta-M scaled
            extinction, if False then the total extinction of the field is used.

        Returns
        -------
        sensor : xr.Dataset
            The same `sensor` as input but modified in place to include the
            calculated optical paths.
        """

        if not isinstance(sensor, xr.Dataset):
            raise TypeError("`sensor` should be an xr.Dataset "
                            "not of type '{}''".format(type(sensor)))
        pyshdom.checks.check_hasdim(sensor, ray_mu='nrays', ray_phi='nrays',
                                    ray_x='nrays', ray_y='nrays', ray_z='nrays')

        camx = sensor['ray_x'].data
        camy = sensor['ray_y'].data
        camz = sensor['ray_z'].data
        cammu = sensor['ray_mu'].data
        camphi = sensor['ray_phi'].data
        total_pix = sensor.sizes['nrays']

        optical_path = pyshdom.core.optical_depth(
            maxnmicro=self._pa.max_num_micro,
            interpmethod=self._interpmethod,
            phasemax=self._phasemax,
            phaseinterpwt=self._phaseinterpwt[:, :self._npts],
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
            extinct=self._extinct[:self._npts],
            albedo=self._albedo[:self._npts],
            iphase=self._iphase[:, :self._npts],
            legen=self._legen,
            npart=self._npart,
            nstleg=self._nstleg,
            deltam=self._deltam,
            deltampath=deltam_scaled_path,
            nleg=self._nleg,
            ml=self._ml
            )
        if deltam_scaled_path:
            sensor['optical_path_deltam'] = (['nrays'], optical_path)
        else:
            sensor['optical_path'] = (['nrays'], optical_path)
        return sensor

    def transmission_integral(self, sensor, field):

        if not isinstance(sensor, xr.Dataset):
            raise TypeError("`sensor` should be an xr.Dataset "
                            " not of type '{}''".format(type(sensor)))
        pyshdom.checks.check_hasdim(sensor, ray_mu='nrays', ray_phi='nrays',
                                    ray_x='nrays', ray_y='nrays', ray_z='nrays')

        camx = sensor['ray_x'].data
        camy = sensor['ray_y'].data
        camz = sensor['ray_z'].data
        cammu = sensor['ray_mu'].data
        camphi = sensor['ray_phi'].data
        total_pix = sensor.sizes['nrays']

        if self.check_solved(verbose=False):
            raise pyshdom.exceptions.SHDOMError(
                "This function can only be run before RTE.solve()"
                )

        transmission_integral = pyshdom.core.transmission_integral(
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
            total_ext=self._total_ext[:self._npts],
            field=field,
            transcut=self._transcut,
            tautol=self._tautol
            )
        sensor['transmission_integral'] = (['nrays'], transmission_integral)
        return sensor



    def min_optical_path(self, sensor, deltam_scaled_path=False, do_all=False):

        if not isinstance(sensor, xr.Dataset):
            raise TypeError("`sensor` should be an xr.Dataset "
                            " not of type '{}''".format(type(sensor)))
        pyshdom.checks.check_hasdim(sensor, ray_mu='nrays', ray_phi='nrays',
                                    ray_x='nrays', ray_y='nrays', ray_z='nrays')

        camx = sensor['ray_x'].data
        camy = sensor['ray_y'].data
        camz = sensor['ray_z'].data
        cammu = sensor['ray_mu'].data
        camphi = sensor['ray_phi'].data
        total_pix = sensor.sizes['nrays']

        if self.check_solved(verbose=False):
            raise pyshdom.exceptions.SHDOMError(
                "This function can only be run before RTE.solve()"
                )
        if do_all:
            #optical_path = np.zeros((npixels, self._npts))
            paths_size = total_pix
        else:
            paths_size = 1
        optical_path = pyshdom.core.min_optical_depth(
            maxnmicro=self._pa.max_num_micro,
            interpmethod=self._interpmethod,
            phasemax=self._phasemax,
            phaseinterpwt=self._phaseinterpwt[:, :self._npts],
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
            extinct=self._extinct[:self._npts],
            albedo=self._albedo[:self._npts],
            iphase=self._iphase[:, :self._npts],
            legen=self._legen,
            npart=self._npart,
            nstleg=self._nstleg,
            deltam=self._deltam,
            deltampath=deltam_scaled_path,
            paths_size=paths_size,
            nleg=self._nleg,
            ml=self._ml
            )
        name = 'min_optical_path'
        if deltam_scaled_path:
            name = name + '_deltam'
        if do_all:
            optical_path_dataset = xr.Dataset(
                data_vars={
                    name:(['x', 'y', 'z', 'npixels'],
                          optical_path[:self._nbpts, :].reshape(
                              self._nx, self._ny, self._nz, -1)
                         )},
                coords={'x': self._grid.x,
                        'y': self._grid.y,
                        'z': self._grid.z,
                       },
            )
        else:
            optical_path_dataset = xr.Dataset(
                data_vars={
                    name:(['x', 'y', 'z'],
                          optical_path[:self._nbpts, 0].reshape(
                              self._nx, self._ny, self._nz)
                         )},
                coords={'x': self._grid.x,
                        'y': self._grid.y,
                        'z': self._grid.z,
                        },
            )
        return optical_path_dataset

    def check_solved(self, verbose=True):
        """
        A simple check on whether the solver solution has actually converged.

        Useful for determining if outputs should be truested. May be used just to
        print a warning or in some contexts the flag may be used to raise an Exception,
        for example.
        """
        #note that evaluation is sequential. If the first criterion is True,
        #the subsequent in an or statement are not even evaluated so this won't
        #cause an error.
        flag = True
        if self._solcrit is None or (self._solcrit >= self._solacc):
            flag = False
            if verbose:
                warnings.warn(
                    "Solution has not converged to the "
                    "specified accuracy log10(SOLCRIT): {:3f} > log10(SOLUTION_ACCURACY): {:3f}. "
                    "Calculated quantities will not be accurate.".format(
                        np.log10(self._solcrit), np.log10(self._solacc)
                        )
                    )

        return flag


    @property
    def spherical_harmonics(self):
        """
        Calculates the mean intensity and 3 spatial fluxes (Fx, Fy, Fz) of the
        radiation field.

        This is equivalent SH_OUT from SHDOM in a dataset but only on the tidy base grid,
        so the same as the netCDF output. If output on the adaptive grid points
        access the self.shterms property as all values are stored there in
        array form.
        The root-mean-square of the higher order radiance terms is also
        returned if the `high_order_radiance` was set to True in the
        numerical_parameters used to initialize the `RTE` object.

        Returns
        -------
        sh_out_dataset : xr.Dataset
            Dataset containing mean intensity and 3 spatial fluxes (Fx, Fy, Fz)
            on SHDOM's base grid.

        Notes
        -----
        This is basically a wrapper for COMPUTE_SH in src/shdomsub2.f
        """

        self.check_solved()
        if self._highorderrad:
            nshout = 5
        else:
            nshout = 4
        if self._shterms is None:
            shterms = pyshdom.core.compute_sh(nshout=nshout,
                                              nstokes=self._nstokes,
                                              npts=self._npts,
                                              srctype=self._srctype,
                                              solarmu=self._solarmu,
                                              solaraz=self._solaraz,
                                              dirflux=self._dirflux[:self._npts],
                                              rshptr=self._rshptr[:self._npts+1],
                                              ml=self._ml,
                                              mm=self._mm,
                                              radiance=self._radiance,
                                             )
            self._shterms = shterms

        if len(self._xgrid) == self._nx1:
            xcoord = self._xgrid
        elif len(self._xgrid)-1 == self._nx1:
            xcoord = self._xgrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )
        if len(self._ygrid) == self._ny1:
            ycoord = self._ygrid
        elif len(self._ygrid)-1 == self._ny1:
            ycoord = self._ygrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )
        sh_out_dataset = xr.Dataset(
            data_vars={
                'mean_intensity': (['x', 'y', 'z'], self._shterms[0, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                'Fx': (['x', 'y', 'z'], self._shterms[1, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                'Fy': (['x', 'y', 'z'], self._shterms[2, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                'Fz': (['x', 'y', 'z'], self._shterms[3, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                },

            coords={'x': xcoord,
                    'y': ycoord,
                    'z': self._zgrid,
                   },
            attrs={
                'long_names':{'Fx': 'Net Flux in x direction',
                              'Fy': 'Net Flux in y direction',
                              'Fz': 'Net Flux in z direction'},
                }
        )
        if self._highorderrad & (self._shterms.shape[0] == 5):
            sh_out_dataset['rms_higher_rad'] = (['x', 'y', 'z'],
                                                self._shterms[-1, :self._nbpts].reshape(
                                                    self._nx1, self._ny1, self._nz)/ \
                                                    (np.sqrt(np.pi*4.0*self._nlm)))
        return sh_out_dataset

    @property
    def adaptive_shterms(self):
        """
        Access the mean intensity and 3 spatial fluxes (Fx, Fy, Fz) at the
        valid adaptive grid points.
        """
        if self._shterms is None:
            self.spherical_harmonics
        return self._shterms[:, :self._npts]

    @property
    def adaptive_fluxes(self):
        """
        Access the hemispherical fluxes at all valid grid points.
        """
        return self._fluxes[:, :self._npts]

    @property
    def fluxes(self):
        """Calculates the hemispherical fluxes from SHDOM but only on the tidy base grid.

        To access the hemispherical fluxes at the adaptive grid points use the
        self.fluxes property.
        """
        self.check_solved()

        if len(self._xgrid) == self._nx1:
            xcoord = self._xgrid
        elif len(self._xgrid)-1 == self._nx1:
            xcoord = self._xgrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )
        if len(self._ygrid) == self._ny1:
            ycoord = self._ygrid
        elif len(self._ygrid)-1 == self._ny1:
            ycoord = self._ygrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )

        fluxes = xr.Dataset(
            data_vars={
                'flux_down': (['x', 'y', 'z'], self._fluxes[0, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                'flux_up': (['x', 'y', 'z'], self._fluxes[1, :self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                'flux_direct': (['x', 'y', 'z'], self._dirflux[:self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz)),
                },
            coords={'x': xcoord,
                    'y': ycoord,
                    'z': self._zgrid,
                   },
            attrs={
                'long_names': {'flux_down': 'Downwelling Hemispherical Flux',
                               'flux_up': 'Upwelling Hemispherical Flux',
                               'flux_direct': 'Downwelling Direct Solar Beam Flux (On Horizontal Surface)'},
                'units': 'Same units as input flux.'
            }
            )
        return fluxes

    @property
    def adaptive_net_flux_div(self):
        """
        Access the net flux divergence on all the used adaptive grid points.
        """
        if self._netfluxdiv is None:
            self.net_flux_divergence
        return self._netfluxdiv[:self._npts]

    @property
    def net_flux_divergence(self):
        """
        Computes the net flux divergence from SHDOM in a dataset but only on the tidy base grid.

        If you want the net flux divergence on the adaptive grid points then access the
        self.adaptive_net_flux_div property.
        """
        self.check_solved()
        if self._netfluxdiv is None:
            netfluxdiv = pyshdom.core.compute_netfluxdiv(
                nstokes=self._nstokes,
                npts=self._npts,
                rshptr=self._rshptr[:self._npts+1],
                srctype=self._srctype,
                solarmu=self._solarmu,
                extinct=self._extinct[:self._npts],
                albedo=self._albedo[:self._npts],
                planck=self._planck[:self._npts],
                dirflux=self._dirflux[:self._npts],
                radiance=self._radiance,
                npart=self._npart
                )
            self._netfluxdiv = netfluxdiv

        if len(self._xgrid) == self._nx1:
            xcoord = self._xgrid
        elif len(self._xgrid)-1 == self._nx1:
            xcoord = self._xgrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )
        if len(self._ygrid) == self._ny1:
            ycoord = self._ygrid
        elif len(self._ygrid)-1 == self._ny1:
            ycoord = self._ygrid[:-1]
        else:
            raise pyshdom.exceptions.SHDOMError(
                "Inconsistent sizes of RTE grid and property grid. "
                "There has been a mistake in interpretation."
                )

        netfluxdiv_dataset = xr.Dataset(
            data_vars={
                'net_flux_div':(['x', 'y', 'z'], self._netfluxdiv[:self._nbpts].reshape(
                    self._nx1, self._ny1, self._nz))
                },
            coords={'x': xcoord,
                    'y': ycoord,
                    'z': self._zgrid,
                   },
            attrs={
                'long_names': {'net_flux_div': 'Net Flux Divergence'},
                'units': ''
            }
            )
        return netfluxdiv_dataset

    def calculate_direct_beam_derivative(self):
        """
        Calculate the geometry of the direct beam at each point and solver.
        Solver is modified in-place.
        If the solver does not have any solar source then empty arrays
        are added so that the signature of the gradient call doesn't need
        to change for each source type.
        """
        if self._srctype != 'T':
            #calculate the solar direct beam on the base grid
            #which ensures the solver has the required information to
            #calculate the derivative.
            self._make_direct()

            direct_derivative_path, direct_derivative_ptr, ierr, errmsg = \
                pyshdom.core.make_direct_derivative(
                    npts=self._npts,
                    bcflag=self._bcflag,
                    gridpos=self._gridpos,
                    npx=self._pa.npx,
                    npy=self._pa.npy,
                    npz=self._pa.npz,
                    delx=self._pa.delx,
                    dely=self._pa.dely,
                    xstart=self._pa.xstart,
                    ystart=self._pa.ystart,
                    zlevels=self._pa.zlevels,
                    ipdirect=self._ipdirect,
                    di=self._di,
                    dj=self._dj,
                    dk=self._dk,
                    epss=self._epss,
                    epsz=self._epsz,
                    xdomain=self._xdomain,
                    ydomain=self._ydomain,
                    cx=self._cx,
                    cy=self._cy,
                    cz=self._cz,
                    cxinv=self._cxinv,
                    cyinv=self._cyinv,
                    czinv=self._czinv,
                    uniformzlev=self._uniformzlev,
                    delxd=self._delxd,
                    delyd=self._delyd
                )
            pyshdom.checks.check_errcode(ierr, errmsg)
        else:
            direct_derivative_ptr = np.zeros(
                (8*(self._pa.npx + self._pa.npy + self._pa.npz), self._npts),
                dtype=np.int32,
                order='F'
            )
            direct_derivative_path = np.zeros(
                (8*(self._pa.npx + self._pa.npy + self._pa.npz), self._npts),
                dtype=np.float32,
                order='F'
            )
        self._direct_derivative_ptr = direct_derivative_ptr
        self._direct_derivative_path = direct_derivative_path

    def calculate_microphysical_partial_derivatives(self, derivative_information):
        """
        Calculate the derivatives of optical properties with respect to the unknowns
        (microphysical or optical).

        Uses an interpolation method and table_data supplied from a
        pyshdom.containers.UnknownScatterers object. Does the delta-M scaling of the
        partial derivatives if appropriate.

        Parameters
        ----------

        """
        self._precompute_phase()
        #solver_derivative_table = #all_derivative_tables[key]
        num_derivatives = sum([len(scatterer_derivative_information.values()) for
                               name, scatterer_derivative_information in derivative_information.items()
                              ])
        self._num_derivatives = np.array(num_derivatives, dtype=np.int32)
        unknown_scatterer_indices = []

        dext = np.zeros(shape=[self._maxpg, num_derivatives], dtype=np.float32)
        dalb = np.zeros(shape=[self._maxpg, num_derivatives], dtype=np.float32)

        # find maximum number of phase pointers across all species.
        num_micros = []
        max_legendre = []
        unknown_scatterer_indices = []
        table_phase_derivative_flag = []

        i = 0
        for scatterer_name, scatterer_derivative_information in derivative_information.items():
            scatterer_index = np.where(scatterer_name == np.array(list(self.medium)))[0][0]
            for variable_derivative in scatterer_derivative_information.values():
                num_micros.append(variable_derivative.num_micro.size)
                max_legendre.append(variable_derivative.legendre_index.size)
                table_phase_flag = 0
                if variable_derivative['derivative_method'] == 'exact':
                    table_phase_flag = 1
                table_phase_derivative_flag.append(table_phase_flag)
                unknown_scatterer_indices.append(scatterer_index+1)
            i += 1

        deriv_max_num_micro = max(num_micros)
        max_legendre = max(max_legendre)
        self._unknown_scatterer_indices = np.array(unknown_scatterer_indices).astype(np.int32)
        self._table_phase_derivative_flag = np.array(table_phase_derivative_flag).astype(np.int32)

        diphase = np.zeros(
            shape=[deriv_max_num_micro, self._maxpg, num_derivatives],
            dtype=np.int32
        )

        dphasewt = np.zeros(
            shape=[deriv_max_num_micro, self._maxpg, num_derivatives],
            dtype=np.float32
        )

        dleg_table = []
        i = 0
        for scat_name, scatterer_derivative_information in derivative_information.items():
            phase_size_sum = 0
            for name, scatterer in self.medium.items():
                if name == scat_name:
                    break
                phase_size_sum += scatterer.sizes['table_index']

            for variable_derivative in scatterer_derivative_information.values():
                assert variable_derivative['derivative_method'] in ('exact', 'table'), 'Bad `derivative_method`'
                dext[:, i] = variable_derivative.extinction.data.ravel()
                dalb[:, i] = variable_derivative.ssalb.data.ravel()
                dphasewt[..., i] = np.pad(
                    variable_derivative.phase_weights.data.reshape((-1,self._maxpg)),
                    ((0,deriv_max_num_micro - variable_derivative.sizes['num_micro']), (0,0)),
                    mode='constant' # default pads with zeros.
                )

                dleg_max = phase_size_sum
                if variable_derivative['derivative_method'] == 'exact':
                    # in this case phase pointer points to dleg table.
                    dleg_max = sum([table.sizes['table_index'] for table in dleg_table])
                    dleg_table.append(variable_derivative.legcoef.pad(
                        {'legendre_index': (0, max_legendre - variable_derivative.legcoef.sizes['legendre_index'])},
                        constant_values=0.0
                    ))
                diphase[..., i] = np.pad(
                    variable_derivative.table_index.data.reshape((-1,self._maxpg)) + dleg_max,
                    ((0,deriv_max_num_micro - variable_derivative.sizes['num_micro']), (0,0)),
                    mode='constant' # default pads with zeros.
                )
                i += 1

        #COPIED FROM solver.RTE
        #In regions which are not covered by any optical scatterer they have an iphasep of 0.
        #In original SHDOM these would be pointed to the rayleigh phase function (which is always included
        #in the legendre table even if there is no rayleigh extinction.)
        #Here, instead we set them to whatever the first phase function is.
        #An arbitrary valid choice can be made as the contribution from these grid points is zero.
        diphase[np.where(diphase == 0)] = 1

        # Concatenate all legendre tables into one table
        legendre_table = xr.concat(dleg_table, dim='table_index')
        if self._pa.nlegp + 1 > legendre_table.sizes['legendre_index']:
            legendre_table = legendre_table.pad(
                {'legendre_index':
                 (0, 1 + self._nleg - legendre_table.sizes['legendre_index'])
                }, constant_values=0.0
            )
        self._dnumphase = legendre_table.sizes['table_index']
        dleg = legendre_table.data

        scaling_factor = np.atleast_3d(np.array([2.0*i+1.0 for i in range(0, self._pa.nlegp+1)]))
        dleg[0, 0, :] = 0.0
        dleg = dleg[:self._nstleg] / scaling_factor

        self._dext = dext
        self._dalb = dalb
        self._diphasep = diphase
        self._dphasewtp = dphasewt
        # temperature derivatives are not yet supported in the python interface.
        self._dtemp = np.zeros((dext.shape))

        # compute lut of phase function derivatives evaluated at scattering angles.
        self._dphasetab, ierr, errmsg = pyshdom.core.precompute_phase_check_grad(
            negcheck=False,
            nstphase=self._nstphase,
            nstleg=self._nstleg,
            nscatangle=self._nscatangle,
            nstokes=self._nstokes,
            dnumphase=self._dnumphase,
            ml=self._ml,
            nlm=self._nlm,
            nleg=self._pa.nlegp,
            dleg=dleg,
            deltam=self._deltam
        )
        pyshdom.checks.check_errcode(ierr, errmsg)
        # now that we have dphasetab we can truncate dleg
        # to the RTE accuracy.
        self._dleg = dleg[:, :self._nleg+1]
        # make property grid to RTE grid pointers and interpolation weights.
        self._optinterpwt, self._interpptr, ierr, errmsg, \
        self._dalbm, self._dextm, self._dfj = \
        pyshdom.core.prepare_deriv_interps(
            gridpos=self._gridpos[:, :self._npts],
            npx=self._pa.npx,
            npy=self._pa.npy,
            npz=self._pa.npz,
            npts=self._npts,
            maxpg=self._maxpg,
            delx=self._pa.delx,
            dely=self._pa.dely,
            xstart=self._pa.xstart,
            ystart=self._pa.ystart,
            zlevels=self._pa.zlevels,
            legen=self._legen,
            numphase=self._pa.numphase,
            nstleg=self._nstleg,
            dleg=self._dleg,
            dnumphase=self._dnumphase,
            phasewtp=self._pa.phasewtp,
            dphasewtp=self._dphasewtp,
            iphasep=self._pa.iphasep,
            diphasep=self._diphasep,
            nleg=self._nleg,
            maxnmicro=self._pa.max_num_micro,
            albedop=self._pa.albedop,
            extinctp=self._pa.extinctp,
            npart=self._npart,
            dalb=self._dalb,
            dext=self._dext,
            numder=self._num_derivatives,
            partder=self._unknown_scatterer_indices,
            doexact=self._table_phase_derivative_flag,
            ml=self._ml,
            deltam=self._deltam,
            albedo=self._albedo[:self._npts],
            iphase=self._iphase[:,:self._npts],
            phaseinterpwt=self._phaseinterpwt[:,:self._npts],
            phasemax=self._phasemax,
            interpmethod=self._interpmethod
        )
        pyshdom.checks.check_errcode(ierr, errmsg)


    def load_solution(self, input_dataset, load_radiance=True):
        """
        Loads the grid and radiance/source spherical harmonics of another
        (e.g. solved) solver.RTE object.

        These properties must have been exported to an xr.Dataset using RTE.save_solution
        and be consistent with the initialization of this RTE object.
        This overwrites the current SHDOM grid related attributes but the radiance
        related attributes are not overwritten until self._init_solution has been called.

        Parameters
        ----------
        input_dataset : xr.Dataset
            The dataset containing the SHDOM grid and radiance/source and pointer
            arrays.
        load_radiance : bool
            If False then only the grid is loaded. If True then both are loaded.

        Raises
        ------
        ValueError
            If `input_dataset` contains a base grid that is inconsistent with self or
            the spherical harmonics and number of Stokes components are inconsistent with self.

        Notes
        -----
        The information used here does not provide a bit perfect restart file for an
        SHDOM solution. A 'restarted' SHDOM solution initialized using saved spherical harmonics
        will differ with a 'single execution' SHDOM solution due to the loss of information
        about the convergence of the SOURCE and sequence acceleration if it was used.
        However, the 'restarted' SHDOM solution SHOULD agree with 'single execution' to
        within the accuracy of the technique.
        """
        if input_dataset['nx'] != self._nx or input_dataset['ny'] != self._ny \
            or input_dataset['nz'] != self._nz:
            raise ValueError(
                "Incompatible grid sizes in loaded solution. self: (nx={}, ny={}, nz={}). "
                "!= loaded: (nx={}, ny={}, nz={}).".format(
                    self._nx, self._ny, self._nz,
                    input_dataset['nx'].data, input_dataset['ny'].data,
                    input_dataset['nz'].data)
                )
        if np.any(input_dataset.xgrid.data != self._xgrid) or np.any(input_dataset.ygrid.data != self._ygrid) \
            or np.any(input_dataset.zgrid.data != self._zgrid):
            raise ValueError(
                "Incompatible grids (xgrid, ygrid, zgrid) between loaded solution and self."
            )

        #read grid into memory here and overwrite normal grid from
        #initialization, no need to worry about this as this is small.
        if self._gridpos.shape[1] < input_dataset.npts.data:
            raise pyshdom.exceptions.SHDOMError(
                "Cannot load solution as loaded grid has more points ({}) "
                "than supported by the current RTE's memory parameters ({})".format(
                    input_dataset.npts.data, self._gridpos.shape[1]
                )
            )
        if self._gridptr.shape[1] < input_dataset.ncells.data:
            raise pyshdom.exceptions.SHDOMError(
                "Cannot load solution as loaded grid has more cells ({}) "
                "than supported by the current RTE's memory parameters ({})".format(
                    input_dataset.ncells.data, self._gridptr.shape[1]
                )
            )

        self._npts, self._ncells, self._nbcells = input_dataset.npts.data, \
        input_dataset.ncells.data, input_dataset.nbcells.data
        self._gridpos[:, :self._npts] = input_dataset.gridpos.data
        self._gridptr[:, :self._ncells] = input_dataset.gridptr.data
        self._treeptr[:, :self._ncells] = input_dataset.treeptr.data
        self._cellflags[:self._ncells] = input_dataset.cellflags.data
        self._neighptr[:, :self._ncells] = input_dataset.neighptr.data

        self._setup_grid_flag = False
        # interpolate the optical properties onto the grid points
        # (including the adaptive ones.)
        self._prepare_optical_properties()
        if load_radiance:
            if input_dataset['nstokes'].data != self._nstokes:
                raise ValueError(
                    "Incompatible number of Stokes components"
                    " in loaded solution. self: nstokes={} != loaded: nstokes={}".format(
                        self._nstokes, input_dataset['nstokes'].data)
                    )
            if input_dataset['ml'] > self._ml or input_dataset['mm'] > self._mm:
                raise ValueError(
                    "Incompatible SphericalHarmonics in loaded solution."
                    " self: (ml={}, mm={}) < loaded: (ml={}, mm={}).".format(
                        self._ml, self._mm, input_dataset['ml'].data, input_dataset['mm'].data
                    )
                )
            for name in ('radiance', 'source', 'fluxes', 'shptr', 'rshptr'):
                if name not in input_dataset.data_vars:
                    raise KeyError(
                        "`input_dataset` does not appear to contain variable '{}' "
                        "which is required.")
            #only store the data here. Without calling input_dataset['radiance'].source.data etc
            #the data won't be read into memory (assuming open_dataset was used to open the file).
            #We just store it here and it will be utilized when RTE.solve is called
            #(through self._init_solution)
            self._restore_data = input_dataset

    def save_solution(self, save_radiances=True):
        """Saves the adaptive grid and radiance/source arrays to an xr.Dataset.

        This can be used to save the solutions to expensive radiative transfer solutions
        and they can then be loaded with self.load_solution and additional
        radiances or linearizations can be evaluated with little extra expense.
        The downside is that they take up a bunch of memory!

        Parameters
        ----------
        save_radiances : bool
            If True then both the adaptive grid and the spherical harmonic expansion
            of the radiance and source fields is returned in `output_dataset`.
            If False then only the adaptive grid structure is returned.

        Returns
        -------
        output_dataset : xr.Dataset
            The dataset containing the adaptive grid and possibly the spherical
            harmonics which describe a solved solution.

        Notes
        -----
        The information returned here does not provide a bit perfect restart file for an
        SHDOM solution. A 'restarted' SHDOM solution initialized using saved spherical harmonics
        will differ with a 'single execution' SHDOM solution due to the loss of information
        about the convergence of the SOURCE and sequence acceleration if it was used.
        However, the 'restarted' SHDOM solution SHOULD agree with 'single execution' to
        within the accuracy of the technique.
        """

        output_dataset = xr.Dataset()
        #specific information to confirm that the solver's grids and
        #spherical harmonics are compatible.
        #e.g. RESTORE_STATE/SAVE_STATE in src/shdomsub3.f
        output_dataset['nstokes'] = self._nstokes
        output_dataset['nx'] = self._nx
        output_dataset['ny'] = self._ny
        output_dataset['nz'] = self._nz
        output_dataset['ml'] = self._ml
        output_dataset['mm'] = self._mm
        output_dataset['nlm'] = self._nlm

        output_dataset['npts'] = self._npts
        output_dataset['ncells'] = self._ncells
        output_dataset['nbcells'] = self._nbcells
        output_dataset['xgrid'] = (['nx1'], self._xgrid)
        output_dataset['ygrid'] = (['ny1'], self._ygrid)
        output_dataset['zgrid'] = (['nz_dim'], self._zgrid)
        output_dataset['gridpos'] = (['xyz', 'npts_dim'], self._gridpos[:, :self._npts])
        output_dataset['gridptr'] = (['8points', 'ncells_dim'], self._gridptr[:, :self._ncells])
        output_dataset['neighptr'] = (['6neighbours', 'ncells_dim'],
                                      self._neighptr[:, :self._ncells])
        output_dataset['treeptr'] = (['parent_child', 'ncells_dim'],
                                     self._treeptr[:, :self._ncells])
        output_dataset['cellflags'] = (['ncells_dim'], self._cellflags[:self._ncells])

        if save_radiances:
            output_dataset['fluxes'] = (['updown', 'npts_dim'], self._fluxes[:, :self._npts])
            output_dataset['shptr'] = (['nptsand1'], self._shptr[:self._npts+1])
            output_dataset['rshptr'] = (['nptsand2'], self._rshptr[:self._npts+2])
            output_dataset['source'] = (['nstokes_dim', 'sourcesize'],
                                        self._source[:, :self._shptr[self._npts]])
            output_dataset['radiance'] = (['nstokes_dim', 'radsize'],
                                          self._radiance[:, :self._rshptr[self._npts]])
        return output_dataset


    def _setup_medium(self, medium):
        """
        A utility function that does all the checks on the validity
        of the medium input.

        Checks that all scatterers in `medium` are on the same (valid) grid
        and that all required variable and dimension names are present with
        specified shapes (if required).

        Parameters
        ----------
        medium : dict, Tuple, List
            Contains the gridded optical properties of scatterers. This is the
             same `medium` variable used in __init__.

        Returns
        -------
        medium_dict : OrderedDict
            The gridded optical properties of all scatterers after filtering.
        """

        if isinstance(medium, OrderedDict):
            medium_dict = medium
        elif isinstance(medium, typing.Dict):
            medium_dict = OrderedDict(medium)
        elif isinstance(medium, (typing.Tuple, typing.List)):
            medium_dict = OrderedDict()
            for i, val in enumerate(medium):
                medium_dict['scatterer_{:03d}'.format(i)] = val
        elif isinstance(medium, xr.Dataset):
            medium_dict = OrderedDict({'scatterer_{:03d}'.format(0): medium})
        else:
            raise TypeError("'medium' argument to RTE must be an OrderedDict, dict or iterable"
                            "of xr.Dataset or an xr.Dataset")

        for name, dataset in medium_dict.items():
            pyshdom.checks.check_optical_properties(dataset, name=name)

        #check that all scatterers are on the same grid
        pyshdom.checks.check_grid_consistency(
            *medium_dict.values(),
            names=list(medium_dict.keys())
        )

        first_scatterer = list(medium_dict.values())[0]
        grid = xr.Dataset({'x': first_scatterer.coords['x'],
                           'y': first_scatterer.coords['y'],
                           'z': first_scatterer.coords['z'],
                           'delx': first_scatterer.delx,
                           'dely': first_scatterer.dely})
        for optional in ('nx', 'ny', 'nz'):
            if optional in first_scatterer:
                grid[optional] = first_scatterer[optional]

        scatterer_wavelengths = [scatterer.attrs['wavelength_center']
                                 for scatterer in medium_dict.values() if
                                 'wavelength_center' in scatterer.attrs]

        if scatterer_wavelengths:
            if not np.allclose(scatterer_wavelengths, self.wavelength):
                ValueError("Scatterers in medium have different wavelengths."
                           " An error has likely occurred in the calculation of mie properties."
                           "To override this exception make sure that "
                           "optical_properties.attrs['wavelength_center'] "
                           "either doesn't exist or are all identical.")
        else:
            warnings.warn("No wavelength specified in optical properties. "
                          "Make sure that the wavelength used is"
                          "consistent with {}".format(self.wavelength),
                          category=RuntimeWarning)

        return medium_dict, grid

    def _setup_atmosphere(self, atmosphere):
        """
        A utility function that does all the checks on the validity
        of the atmosphere input and preprocesses any temperature and
        gas absorption data.

        Chooses whether to use the 1D gas absorption array or to create a new
        3D scatterer object with the gas absorption.

        Parameters
        ----------
        atmosphere : xr.Dataset
            Contains temperature and possibly gas absorption data. Should be consistent
            with any rayleigh scatterer that is added.
        """
        self._pa.tempp = np.zeros(shape=(self._maxpg,), dtype=np.float32) + 273.0
        if atmosphere is None:
            if self._srctype in ('T', 'B'):
                raise KeyError("'temperature' variable was not specified in "
                               "`atmosphere` despite using thermal source.")
        elif isinstance(atmosphere, xr.Dataset):
            #atmosphere should be on the rte_grid.
            pyshdom.checks.check_grid(atmosphere)
            if not np.all(atmosphere.x == self._grid.x) & np.all(atmosphere.y == self._grid.y) & \
                np.all(atmosphere.z == self._grid.z):
                raise pyshdom.exceptions.GridError("`atmosphere` does not have a "
                                                   "consistent grid with medium.")

            if 'temperature' in atmosphere.data_vars:
                pyshdom.checks.check_positivity(atmosphere, 'temperature')
                pyshdom.checks.check_hasdim(atmosphere, temperature=['x', 'y', 'z'])
                if np.any(atmosphere.temperature >= 350.0) or \
                   np.any(atmosphere.temperature <= 150.0):
                    warnings.warn("Temperatures in `atmosphere` are out of "
                                  "Earth's range [150.0, 350.0].")
                self._pa.tempp[:] = atmosphere.temperature.data.ravel()
            elif self._srctype in ('T', 'B'):
                raise KeyError("'temperature' variable was not specified in "
                               "`atmosphere` despite using thermal source.")

            if 'gas_absorption' in atmosphere.data_vars:
                pyshdom.checks.check_positivity(atmosphere, 'gas_absorption')
                pyshdom.checks.check_hasdim(atmosphere, gas_absorption=['x', 'y', 'z'])
                if np.all(atmosphere.gas_absorption[0, 0] == atmosphere.gas_absorption):
                    self._pa.nzckd = atmosphere.sizes['z']
                    self._pa.zckd = atmosphere.z.data
                    self._pa.gasabs = atmosphere.gas_absorption[0, 0].data
                else:
                    warnings.warn("'gas_absorption' does not collapse to 1D so is being added "
                                  "to `medium`.")
                    if 'gas_absorption' in self.medium:
                        KeyError("'gas_absorption' key was already in `medium`.")
                    else:
                        gas_absorption_scatterer = xr.Dataset(
                            data_vars={
                                'extinction': (['x', 'y', 'z'], atmosphere.gas_absorption.data),
                                'ssalb': (['x', 'y', 'z'],
                                           np.zeros(atmosphere.gas_absorption.shape)),
                                'table_index': (['num_micro', 'x', 'y', 'z'],
                                                np.zeros((1,)+atmosphere.gas_absorption.shape,
                                                dtype=np.int)),
                                'phase_weights': (['num_micro', 'x', 'y', 'z'],
                                                np.ones((1,)+atmosphere.gas_absorption.shape,
                                                dtype=np.float32)),
                                'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
                                            np.zeros((6, 0, 0)))
                            },
                            coords={
                                'x': atmosphere.x,
                                'y': atmosphere.y,
                                'z': atmosphere.z,
                            }
                        )
                        self.medium['gas_absorption'] = gas_absorption_scatterer
            else:
                warnings.warn("No gas absorption found in `atmosphere` dataset."
                              " Name should be 'gas_absorption'.")
        else:
            raise TypeError("`atmosphere` should be an xr.Dataset or None")
        return atmosphere

    def _setup_source(self, source):
        """
        A utility function for preprocessing the information about the source
        for the RTE solution. (see source.py).

        Parameters
        ----------
        source : xr.Dataset
            Describes the source type. See source.py for details.
        """
        pyshdom.checks.check_positivity(source, 'wavelength', 'solarflux',
                                        'skyrad')
        solarmu = source.solarmu.data
        if not (-1.0 <= solarmu) & (solarmu < 0.0):
            raise pyshdom.exceptions.OutOfRangeError(
                "solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                "The SHDOM convention for solar direction is that it points"
                "in the direction of the propagation of radiance.".format(solarmu))

        pyshdom.checks.check_range(source, solaraz=(-np.pi, np.pi))
        self.wavelength = source.wavelength.data
        self._srctype = source.srctype.data
        self._solarflux = source.solarflux.data
        self._solarmu = source.solarmu.data
        self._solaraz = source.solaraz.data
        self._skyrad = source.skyrad.data
        self._units = source.units.data
        self._waveno = source.wavenumber.data

        if self._srctype not in ('T', 'S', 'B'):
            raise ValueError("Invalid source type '{}'. Make sure to choose "
                             "a source from source.py".format(self._srctype))
        if self._units not in ('R', 'T'):
            raise ValueError("Invalid units '{}'. Make sure to choose "
                             "a source from source.py".format(self._units))
        return source

    def _setup_surface(self, surface):
        """
        A utility function for preprocessing and checks on the surface properties
        that are input.

        Parameters
        ----------
        surface : xr.Dataset
            See surface.py for details.
        """
        # Surface parameters
        pyshdom.checks.check_range(surface, gndalbedo=(0.0, 1.0))
        pyshdom.checks.check_positivity(
            surface, 'gndtemp', 'delxsfc', 'delysfc', 'maxsfcpars',
            'nsfcpar', 'nxsfc', 'nxsfc'
        )

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
        self._sfcgridparms = np.zeros((self._nsfcpar,
                                       self._maxnbc),
                                      dtype=np.float32,
                                      order='F')

        if self._sfctype not in ('FL', 'VL', 'VW', 'VD', 'VO', 'VR'):
            raise ValueError("surface type '{}' not recognized. "
                             "Make sure to choose a supported surface from"
                             " surface.py".format(self._sfctype))

        #I don't know if this causes uniformly unwanted behaviour. At present
        #it is only warned against. - JRLoveridge 2021/02/21
        if (self._nxsfc > 1) & \
            (self._delxsfc * (self._nxsfc + 1) < self._grid.x[-1] - self._grid.x[0]):
            warnings.warn('Specified surface does not span RTE grid in the x direction.')
        if (self._nysfc > 1) & \
            (self._delysfc * (self._nysfc + 1) < self._grid.y[-1] - self._grid.y[0]):
            warnings.warn('Specified surface does not span RTE grid in the y direction.')

        if (self._sfctype[-1] in ('R', 'O')) and (self._nstokes > 1):
            raise ValueError("surface brdf '{}' is only"
                             " supported for unpolarized radiative"
                             " transfer (num_stokes=1)".format(surface.name.data))

        if self._sfctype.endswith('L'):
            self._maxbcrad = 2 * self._maxnbc
        else:
            self._maxbcrad = int((2 + self._nmu * self._nphi0max / 2) * self._maxnbc)

        if (self._gndtemp <= 150.0) or (self._gndtemp >= 350.0):
            warnings.warn("Ground temperature '{}' out of Earth range. "
                          "[150.0, 350.0]".format(self._gndtemp))

        return surface

    def _setup_numerical_params(self, numerical_params):
        """
        Set the numerical parameters of the SHDOM forward solver.

        Parameters
        ----------
        numerical_params : xr.Dataset
            An object which encapsulate numerical parameters such as number of azimuthal and
            zenith bins. See config.py & default_config.json
        """
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
        self._iterfixsh = int(numerical_params.iterfixsh.data)
        self._tautol = numerical_params.tautol.data
        self._angle_set = numerical_params.angle_set.data
        self._transcut = numerical_params.transcut.data

        if self._deltam.dtype != np.bool:
            raise TypeError("numerical_params.deltam should be of boolean type.")
        if self._highorderrad.dtype != np.bool:
            raise TypeError("numerical_params.high_order_radiance should be of boolean type.")

        if self._nphi < self._nmu:
            warnings.warn("Usually want NPHI > NMU. NMU={}, NPHI={}".format(self._nmu, self._nphi))

        for boundary in ('x_boundary_condition', 'y_boundary_condition'):
            if numerical_params[boundary] not in ('open', 'periodic'):
                raise ValueError("{} should be either ('open', 'periodic')".format(boundary))

        self._ipflag = numerical_params.ip_flag.data
        if self._ipflag not in np.arange(8):
            raise ValueError("numerical_params.ip_flag should be an integer "
                             "in (0, 1, 2, 3, 4, 5, 6, 7, 8) not '{}'".format(self._ipflag))
        # bcflag is set in _setup_grid and ipflag may be modified there to handle the
        # nx/ny = 1 special case.

        if self._angle_set not in (1, 2, 3):
            raise ValueError(
                "Numerical Parameter 'angle_set' must be in the set (1, 2, 3). See "
                "default_config.json or shdom.txt for more details."
            )
        if not (self._transcut >= 0.0 or self._transcut <= 5e-5):
            warnings.warn("TRANSCUT should be a small number in (0, 5e-5) not {}".format(self._transcut))

        return numerical_params

    def _setup_grid(self, grid):
        """
        A utility function to set the base grid and related grid structures for SHDOM.

        Parameters
        ----------
        grid: xr.Dataset
            The grid containing x, y, z coordinates. Must be a valid
            SHDOM grid. see grid.py for details.

        Notes
        -----
        This is a private method, it doesn't perform a check that `grid` is a valid
        SHDOM grid.
        """
        def ibits(val, bit, ret_val):
            if val & 2 ** bit:
                return ret_val
            return 0

        # Set shdom property array
        self._pa.npx = grid.dims['x']
        self._pa.npy = grid.dims['y']
        self._pa.npz = grid.dims['z']
        #All grids start from 0.0, xstart should be used only
        #for MPI as each worker will have different starting positions.
        assert np.allclose(grid.x[0], 0.0), 'X-dimension of property grid should start from 0.0'
        assert np.allclose(grid.y[0], 0.0), 'Y-dimension of property grid should start from 0.0'
        self._pa.xstart = grid.x[0]
        self._pa.ystart = grid.y[0]

        self._pa.delx = grid.delx.data
        self._pa.dely = grid.dely.data
        self._pa.zlevels = grid.z.data

        # Initialize shdom internal grid sizes to property array grid
        # If SHDOM grid sizes aren't specified. If they are then the
        # override.
        self._nx = self._pa.npx if 'nx' not in grid.data_vars else grid.nx.data
        self._ny = self._pa.npy if 'ny' not in grid.data_vars else grid.ny.data
        self._nz = self._pa.npz if 'nz' not in grid.data_vars else grid.nz.data

        # gridtype = 'P': Z levels taken from property file
        self._gridtype = 'E' if 'nz' in grid.data_vars else 'P'
        # If a vertical resolution is specified we use an equispaced vertical
        # grid for SHDOM overriding the property grid vertical levels.

        if self._nz < self._pa.npz:
            warnings.warn(
                "SHDOM vertical resolution nz={}, is less than property grid"
                " npz={}".format(self._nz, self._pa.npz)
                )
        if self._nx < self._pa.npx:
            warnings.warn(
                "SHDOM X resolution nz={}, is less than property grid"
                " npz={}".format(self._nx, self._pa.npx)
                )
        if self._ny < self._pa.npy:
            warnings.warn(
                "SHDOM Y resolution nz={}, is less than property grid"
                " npz={}".format(self._ny, self._pa.npy)
                )

        self._maxpg = grid.dims['x'] * grid.dims['y'] * grid.dims['z']
        self._maxnz = grid.dims['z']

        # Set the full domain grid sizes (variables end in t)
        self._nxt, self._nyt, self._npxt, self._npyt = \
            self._nx, self._ny, self._pa.npx, self._pa.npy
        self._nx, self._ny, self._nz = \
            max(1, self._nx), max(1, self._ny), max(2, self._nz)

        # if single plane or column then force independent pixel mode.
        if self._nx == 1:
            self._ipflag = self._ipflag | (1<<0) #first bit for X-dim.
        if self._ny == 1:
            self._ipflag = self._ipflag | (1<<1) # second bit for Y-dim.

        # set boundary condition flag based on updated information about ipflag.
        self._bcflag = 0
        if (self.numerical_params.x_boundary_condition == 'open') & (self._ipflag in (0, 2, 4, 6)):
            self._bcflag += 1
        if (self.numerical_params.y_boundary_condition == 'open')& (self._ipflag in (0, 1, 4, 5)):
            self._bcflag += 2

        # Set up base grid point actual size (NX1xNY1xNZ)
        self._nx1, self._ny1 = self._nx + 1, self._ny + 1
        if self._bcflag & 5 or ibits(self._ipflag, 0, 1):
            self._nx1 -= 1
        if self._bcflag & 7 or ibits(self._ipflag, 1, 1):
            self._ny1 -= 1
        self._nbpts = self._nx1 * self._ny1 * self._nz

        # Calculate the number of base grid cells depending on the BCFLAG
        self._nbcells = (self._nz - 1) * (self._nx + ibits(self._bcflag, 0, 1) - \
                                          ibits(self._bcflag, 2, 1)) * \
                                (self._ny + ibits(self._bcflag, 1, 1) - \
                                ibits(self._bcflag, 3, 1))

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

        self._setup_memory()

        self._npts, self._ncells, self._gridpos, self._gridptr, self._neighptr, \
        self._treeptr, self._cellflags = pyshdom.core.init_cell_structure(
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
        """A utility function to initialize internal memory parameters.

        Raises
        ------
        pyshdom.exceptions.SHDOMError
            If there are some errors in the memory
            allocation of the arrays that would cause the process to be exited in
            original SHDOM.
        """

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

        #max_legendre and numphase are recalculated
        max_legendre = max([scatterer.sizes['legendre_index'] for
                            scatterer in self.medium.values()])
        numphase = sum([scatterer.sizes['table_index'] for
                        scatterer in self.medium.values()])
        self._memword = self._nmu * (2 + 2 * self._nphi + 2 * self._nlm + 2 * 33 * 32) \
                        + 4.5*self._maxpg + 2*numphase*self._nstleg*max_legendre
        # Guess maximum number of grid points, cells, SH vector size needed
        # but don't let MAX_TOTAL_MB be exceeded

        self._max_int32_size = 2147483647
        #This is the typical output of HUGE(INTEGER) in Fortran which is what SHDOM uses
        #in the equivalent section of code. This may be inappropriate for
        #other computer architectures where the default INTEGER precision
        #is 64, I don't know anything about that sort of stuff.
        #Previously in pyshdom this was sys.maxsize - but that is max value for a int64
        #which is way too big. For very large grids there was an error attempting
        #to execute the fortran subroutines as the integer describing
        #array sizes overflowed. Hopefully using this
        #will cause the error to be caught here when setting sizes.
        # -JRLoveridge 2021/02/17.

        # This was done in SHDOM and would reduce the maximum memory to < 4GB.
        # This is probably legacy from 32bit systems.(?)
        # This is a different issue to the array allocation issue when
        # int32 integers overflow.
        # if self._max_total_mb * 1024 ** 2 > 1.8 * self._max_int32_size:
        #     self._max_total_mb = 1.8 * self._max_int32_size / 1024 ** 2
        #     warnings.warn("MAX_TOTAL_MB reduced to fit memory model: {}".format(
        #         self._max_total_mb))
        #It seems a better idea would be to know the machine's memory.
        #and not take more than 90% of it. This hasn't been tested so the 90%
        #value may need to be adjusted. - JRLoveridge 2021/02/20
        if self._max_total_mb * 1024**2 > 0.9*psutil.virtual_memory().total:
            self._max_total_mb = 0.9*psutil.virtual_memory().total
            warnings.warn("MAX_TOTAL_MB reduced to fit memory model: {}".format(
                self._max_total_mb))

        if self._splitacc < 0.0:
            self._adapt_grid_factor = 1.0

        self._big_arrays = 3
        if not self._accelflag:
            self._big_arrays = 2

        wantmem = self._adapt_grid_factor * self._nbpts * (
            28 + 16.5 * self._cell_to_point_ratio + self._nphi0max * self._nstokes + \
            self._num_sh_term_factor * self._nstokes * self._nlm * self._big_arrays)
        reduce = min(1.0, ((self._max_total_mb * 1024 ** 2) / 4 - self._memword) / wantmem)

        self._adapt_grid_factor *= reduce
        if self._adapt_grid_factor < 1.0:
            raise pyshdom.exceptions.SHDOMError(
                "max_total_mb memory limit ({}) exceeded with just base grid "
                "(final adapt_grid_factor={}). Either increase"
                "memory or reduce accuracy parameters (num_mu_bins={}, num_phi_bins={},"
                " num_sh_term_factor={})".format(
                    self._max_total_mb, self._adapt_grid_factor,
                    self.numerical_params.num_mu_bins.data,
                    self.numerical_params.num_phi_bins.data,
                    self._num_sh_term_factor
                )
            )
        if reduce < 1.0:
            print('adapt_grid_factor reduced to ', self._adapt_grid_factor)

        wantmem *= reduce
        # below is only necessary for 32 bit systems, I guess. - JRLoveridge 2021/03/02
        # if wantmem > self._max_int32_size:
        #     raise pyshdom.exceptions.SHDOMError(
        #         "Number of words of memory ({}) exceeds max integer size: {}".format(
        #             wantmem, self._max_int32_size
        #         ))

        #These are the sizes of the main arrays.
        self._maxig = int(self._adapt_grid_factor * self._nbpts)
        self._maxic = max(int(self._cell_to_point_ratio * self._maxig), self._nbcells)
        maxiv = int(self._num_sh_term_factor * self._nlm * self._maxig)
        if maxiv < self._nbpts*4:
            warnings.warn("User specified MAXIV={} is smaller than the minimum needed "
                "to initialize the radiance fields. It has been increased to {}.".format(
                    maxiv, self._nbpts*4
                    )
                )
        self._maxiv = max(maxiv, self._nbpts*4) #at minimum we need enough to allocate the
                                                # L=1 radiances on the base grid
                                                # or we will segfault in COMPUTE_SOURCE
                                                # as pointers are made in INIT_RADIANCE
                                                # assuming there is enough room.
                                                # For balancing of adaptive grid points and
                                                # spherical harmonics you are on your own but
                                                # should trigger informative exceptions.

        self._maxido = self._maxig * self._nphi0max

        #important that these numbers are smaller than sys._max_int32_size or
        #the f2py Fortran routines will fail to run. This sets an upper bound
        #on the number of gridpoints that can be used.
        #To upgrade this, it will be necessary to upgrade the precision in
        #the Fortran subroutines.

        if self._maxiv > self._max_int32_size:
            raise pyshdom.exceptions.SHDOMError(
                "Size of largest array (RADIANCE) with shape (NSTOKES, MAXIV+MAXIG) "
                "likely exceeds the max integer number of bytes. This will result in array probably"
                " failing to allocate. Current Size: {} > Max Size: {}".format(
                    self._maxiv, self._max_int32_size)
                )

        if 4.0*8.0*self._maxic > self._max_int32_size:
            raise pyshdom.exceptions.SHDOMError(
                "Size of gridptr array with shape (8*MAXIC) "
                "likely exceeds the max integer number of bytes. This will result in array probably"
                " failing to allocate. {} > {}".format(
                    4.0*8.0*self._maxic, self._max_int32_size))

        self._maxnbc = int(self._maxig * 3 / self._nz)

    def _prepare_optical_properties(self):
        """
        A utility function that prepares the optical properties for the initialization
        and solution of SHDOM.

        This forms the SHDOM property arrays (which hold the raw optical properties)
        as well as the arrays that will be used in the solution of SHDOM.
        The SHDOM optical properties includes adaptive points and may be delta-M scaled.
        """
        # Iterate over all particle types and aggregate the legendre scattering table

        # sefl._maxpg is different to self._nbpts as that has the extra boundary points.
        self._pa.extinctp = np.zeros(shape=[self._maxpg, len(self.medium)], dtype=np.float32)
        self._pa.albedop = np.zeros(shape=[self._maxpg, len(self.medium)], dtype=np.float32)

        # find maximum number of phase pointers across all species.
        self._pa.max_num_micro = max(
            [scatterer.num_micro.size for scatterer in self.medium.values()]
        )

        self._pa.iphasep = np.zeros(
            shape=[self._pa.max_num_micro, self._maxpg, len(self.medium)],
            dtype=np.int32
        )
        self._pa.phasewtp = np.zeros(
            shape=[self._pa.max_num_micro, self._maxpg, len(self.medium)],
            dtype=np.float32
        )

        for i, scatterer in enumerate(self.medium.values()):
            self._pa.extinctp[:, i] = scatterer.extinction.data.ravel()
            self._pa.albedop[:, i] = scatterer.ssalb.data.ravel()
            self._pa.iphasep[..., i] = scatterer.table_index.pad(
                {'num_micro': (0, self._pa.max_num_micro-scatterer.num_micro.size)},
                 constant_values=1
                ).data.reshape(self._pa.max_num_micro, -1) + self._pa.iphasep.max()
            self._pa.phasewtp[..., i] = scatterer.phase_weights.pad(
                {'num_micro': (0, self._pa.max_num_micro-scatterer.num_micro.size)},
                constant_values=0.0
                ).data.reshape(self._pa.max_num_micro, -1)

        #In regions which are not covered by any optical scatterer they have an iphasep of 0.
        #In original SHDOM these would be pointed to the rayleigh phase function
        #(which is always included in the legendre table even if there is no rayleigh
        # extinction.)
        #Here, instead we set them to whatever the first phase function is.
        #An arbitrary valid choice can be made as the contribution from these grid points is zero.
        self._pa.iphasep[np.where(self._pa.iphasep == 0)] = 1

        # Concatenate all scatterer tables into one large table
        max_legendre = max([scatterer.sizes['legendre_index'] for
                            scatterer in self.medium.values()])
        padded_legcoefs = [scatterer.legcoef.pad(
            {'legendre_index':
             (0, max_legendre - scatterer.legcoef.sizes['legendre_index'])
             }, constant_values=0.0)
                           for scatterer in self.medium.values()]

        legendre_table = xr.concat(padded_legcoefs, dim='table_index')

        self._pa.numphase = legendre_table.sizes['table_index']

        if np.any(self._pa.iphasep < 1) or np.any(self._pa.iphasep > self._pa.numphase):
            raise pyshdom.exceptions.OutOfRangeError("Phase function indices are out of bounds.")

        # Determine the number of legendre coefficient for a given angular resolution
        self._nleg = self._ml + 1 if self._deltam else self._ml

        self._pa.nlegp = max(legendre_table.sizes['legendre_index'] - 1, self._nleg)
        self._nscatangle = max(36, min(721, 2 * self._pa.nlegp))

        # Check if legendre table needs padding. It will only need
        # padding if angular resolution is larger than the number of
        # non-zero phase function legendre coefficients.
        if self._pa.nlegp > legendre_table.sizes['legendre_index']:
            legendre_table = legendre_table.pad(
                {'legendre_index':
                 (0, 1 + self._pa.nlegp - legendre_table.sizes['legendre_index'])
                }, constant_values=0.0
            )

        # Check if scalar or vector RTE as they are treated differently
        # in core.transfer_pa_to_grid.
        if self._nstokes == 1:
            self._pa.legenp = legendre_table.isel(
                stokes_index=0).data.ravel(order='F').astype(np.float32)
        else:
            self._pa.legenp = legendre_table.data.ravel(order='F').astype(np.float32)

        self._maxasym = np.max(legendre_table.isel(stokes_index=0, legendre_index=1).data / 3.0)
        self._maxpgl = self._nstleg * self._pa.numphase * (legendre_table.sizes['legendre_index']+1)

        if legendre_table.sizes['table_index'] > 0:
            self._maxigl = self._nstleg*legendre_table.sizes['table_index'] * \
            (self._pa.nlegp + 1)
        else:
            self._maxigl = self._maxig * (legendre_table.sizes['legendre_index'] + 1)

        self._npart = len(self.medium)
        self._nstphase = min(self._nstleg, 2)

        # Transfer property arrays into internal grid structures
        # also does the delta-M scaling if necessary.
        # ._extinct etc are the arrays used to solve the RTE. They still
        # have shape=(npts, nparticles) as the relative contribution of each
        # particle to the total extinction is needed for the mixing of the phase
        # functions during the source function calculation.
        # Note that this is not perfectly optimized in terms of space. Especially
        # the phaseinterpwt and iphase.
        # Technically, the albedo/extinct of each
        # species is only needed for the source computation and could be absorbed
        # into the phaseinterpwt which would remove the need for a 'particle' dimension
        # and wasted space when some 'particles' have much more micro dimensions
        # than others. Important thing is making sure the derivatives work similarly.
        self._temp, self._planck, self._extinct, self._albedo, self._legen, \
        self._iphase, self._total_ext, self._extmin, self._scatmin,         \
        self._albmax, ierr, errmsg, self._phaseinterpwt                  \
         = pyshdom.core.transfer_pa_to_grid(
            phasewtp=self._pa.phasewtp,
            maxnmicro=self._pa.max_num_micro,
            nlegp=self._pa.nlegp,
            phasemax=self._phasemax,
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
            maxpg=self._maxpg,
            ml=self._ml,
            mm=self._mm,
            numphase=self._pa.numphase,
            deltam=self._deltam,
            units=self._units,
            waveno=self._waveno,
            wavelen=self.wavelength,
            gridpos=self._gridpos,
            nleg=self._nleg,
            maxig=self._maxig,
            npx=self._pa.npx,
            npy=self._pa.npy,
            npz=self._pa.npz,
            srctype=self._srctype,
            npts=self._npts,
            interpmethod=self._interpmethod,
            )
        pyshdom.checks.check_errcode(ierr, errmsg)

        #calculate cell averaged extinctions so that warnings can be raised
        #about optical thickness of cells. High optical thickness across a cell
        #leads to lower accuracy for SHDOM.
        reshaped_ext = self._total_ext[:self._nbpts].reshape(self._nx1, self._ny1, self._nz)
        if self._nx1 == 1:
            cell_averaged_extinct = (reshaped_ext[:, 1:, 1:] + reshaped_ext[:, 1:, :-1] +   \
                                     reshaped_ext[:, :-1, 1:] + reshaped_ext[:, :-1, :-1])/4.0
        elif self._ny1 == 1:
            cell_averaged_extinct = (reshaped_ext[:1, :, 1:] + reshaped_ext[:1, :, :-1] + \
                                 reshaped_ext[:-1, :, 1:] + reshaped_ext[:-1, :, :-1])/4.0
        else:
            cell_averaged_extinct = (reshaped_ext[1:, 1:, 1:] + reshaped_ext[1:, 1:, :-1] +   \
                                     reshaped_ext[1:, :-1, 1:] + reshaped_ext[1:, :-1, :-1] + \
                                     reshaped_ext[:-1, 1:, 1:] + reshaped_ext[:-1, 1:, :-1] + \
                                     reshaped_ext[:-1, :-1, 1:] + reshaped_ext[:-1, :-1, :-1])/8.0
        cell_volume = (np.diff(self._zgrid)*np.diff(self._xgrid.data)[0]*np.diff(self._ygrid.data)[0])**(1/3)
        cell_tau_approx = cell_volume[np.newaxis, np.newaxis, :]*cell_averaged_extinct
        number_thick_cells = np.sum(cell_tau_approx >= 2.0)

        # compute something for the gradient later.
        self._maxsubgridints = int(np.nanmean(cell_tau_approx) / self._tautol)

        if number_thick_cells > 0:
            warnings.warn("Number of SHDOM grid cells with optical depth greater than 2: '{}'. "
                          "Max cell optical depth: '{}'".format(
                              number_thick_cells, np.max(cell_tau_approx)
                              )
                         )

        if (not self._deltam) & (self._srctype != 'T') & (self._maxasym > 0.5):
            warnings.warn("Delta-M should be use for solar radiative transfer problems with highly "
                          "peaked phase functions.")

    def _release_big_arrays(self):
        """
        Destroy the big arrays if they exist to free memory.
        This could be called after solving to save memory
        in a large number of sequential but independent RTE
        solutions.
        """
        self._source = None
        self._radiance = None
        self._delsource = None
        self._work = None
        self._work1 = None
        self._work2 = None
        self._sfc_brdf_do = None
        self._fluxes = None
        self._dirflux = None
        self._bcrad = None

    def _init_solution(self, setup_grid=True):
        """
        Initialize the solution (I, J fields) from the direct transmission and
        a 1D Delta-Eddington two-stream model.
        Many arrays are initialized including the dirflux (the direct solar transmission)
        and the 'big' source, radiance and delsource (if acceleration_flag) arrays.

        Parameters
        ----------
        setup_grid : bool
            If True, then the SHDOM grid is reinitialized to the base grid.
            If you want to reuse an SHDOM grid but reperform a solution then you
            can set this to False. Note that this will be overwritten by any
            adaptive grid provided through RTE.load_solution.
        """
        if not self._setup_grid_flag:
            warnings.warn(
                "`setup_grid` flag has been overwritten as an adaptive grid"
                " has been loaded through self.load_solution."
                )
            setup_grid = False
        if setup_grid:
            self._setup_grid(self._grid)

        # Restart solution criteria
        self._oldnpts = 0
        self._solcrit = 1.0
        self._iters = 0

        #initialize the zero-ed large arrays.
        #Doing this here rather than in Fortran releases the old memory of
        #self._radiance etc if it exists which prevents doubling of memory if
        # running RTE.solve more than one time.
        self._source = np.zeros((self._nstokes, self._maxiv), dtype=np.float32, order='F')
        if self._accelflag:
            self._delsource = np.zeros((self._nstokes, self._maxiv), dtype=np.float32, order='F')
            self._ndelsource = self._maxiv
        else:
            self._delsource = np.zeros((self._nstokes, 1), dtype=np.float32, order='F')
            self._ndelsource = 1
        self._radiance = np.zeros((self._nstokes, self._maxiv+self._maxig), dtype=np.float32, order='F')
        self._dirflux = np.zeros((self._maxig), dtype=np.float32, order='F')
        self._work1 = np.zeros((8*self._maxig), dtype=np.int32, order='F')
        self._work = np.zeros((self._maxido*self._nstokes), dtype=np.float32, order='F')
        self._work2_size = max((self._maxig*self._nstokes, self._maxic))
        self._work2 = np.zeros((self._work2_size), dtype=np.float32, order='F')
        self._bcrad = np.zeros((self._nstokes, self._maxbcrad), dtype=np.float32, order='F')
        self._fluxes = np.zeros((2, self._maxig), dtype=np.float32, order='F')

        self._shptr = np.zeros((self._maxig+1), dtype=np.int32, order='F')
        self._rshptr = np.zeros((self._maxig+2), dtype=np.int32, order='F')

        if self._restore_data is not None:
            #if data is made available to use instead of
            #initializing the radiance/source/fluxes with 1D delta-Eddington then
            #this is done here.
            self._inradflag = False #this flag turns off the 1D radiance initialization
                                    #in core.init_solution.
            self._shptr[:self._npts+1] = self._restore_data.shptr.data
            self._rshptr[:self._npts+2] = self._restore_data.rshptr.data
            self._source[:, :self._shptr[self._npts]] = self._restore_data.source.data
            self._radiance[:, :self._rshptr[self._npts]] = self._restore_data.radiance.data
            self._fluxes[:, :self._npts] = self._restore_data.fluxes.data
        else:
            self._inradflag = True

        #set the cached spherical_harmonics and net flux divergence to None
        #as a new solution has been formed.
        self._netfluxdiv = None
        self._shterms = None

        self._nang, self._nphi0, self._mu, self._phi, self._wtdo, \
        self._sfcgridparms, self._ntoppts, self._nbotpts, self._bcptr, \
        self._rshptr, self._shptr, self._oshptr, self._source, self._delsource, \
        self._radiance, self._fluxes, self._dirflux, self._ylmsun, self._uniformzlev, \
        self._pa.extdirp, self._bcrad, self._cx, self._cy, self._cz, self._cxinv, \
        self._cyinv, self._czinv, self._ipdirect, self._di, self._dj, self._dk, \
        self._epss, self._epsz, self._xdomain, self._ydomain, self._delxd, \
        self._delyd, self._deljdot, self._deljold, self._deljnew, self._jnorm, \
        self._fftflag, self._cmu1, self._cmu2, self._wtmu, self._cphi1, \
        self._cphi2, self._wphisave, self._work, self._work1, self._work2, \
        self._uniform_sfc_brdf, self._sfc_brdf_do, ierr, errmsg \
         = pyshdom.core.init_solution(
            ordinateset=self._angle_set,
            phasewtp=self._pa.phasewtp,
            maxnmicro=self._pa.max_num_micro,
            adjflag=self._adjflag,
            nlegp=self._pa.nlegp,
            phasemax=self._phasemax,
            interpmethod=self._interpmethod,
            phaseinterpwt=self._phaseinterpwt,
            work=self._work,
            work2_size=self._work2_size,
            maxpg=self._maxpg,
            ndelsource=self._ndelsource,
            work1=self._work1,
            work2=self._work2,
            fluxes=self._fluxes,
            bcrad=self._bcrad,
            dirflux=self._dirflux,
            delsource=self._delsource,
            source=self._source,
            shptr=self._shptr,
            rshptr=self._rshptr,
            inradflag=self._inradflag,
            radiance=self._radiance,
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
            cellflags=self._cellflags,
            treeptr=self._treeptr,
            neighptr=self._neighptr,
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
        pyshdom.checks.check_errcode(ierr, errmsg)

    def _precompute_phase(self):
        """
        Precompute angular scattering for the entire legendre table.
        Perform a negativity check. (negcheck=True).
        """
        # Use the property phase table as LEGEN is truncated to
        # only include the L indices needed for the solver
        # whereas legenp holds the full phase functions.
        legen_reshape = self._pa.legenp.reshape(
            (self._nstleg, self._pa.nlegp+1,self._pa.numphase),
            order='F'
        )
        self._phasetab, errmsg, ierr = pyshdom.core.precompute_phase_check(
            negcheck=True,
            nscatangle=self._nscatangle,
            numphase=self._pa.numphase,
            nstphase=self._nstphase,
            nstokes=self._nstokes,
            nstleg=self._nstleg,
            nleg=self._pa.nlegp,
            ml=self._ml,
            nlm=self._nlm,
            legen=legen_reshape,
            deltam=self._deltam,
        )
        pyshdom.checks.check_errcode(ierr, errmsg)

    def _make_direct(self):
        """
        Compute the direct transmission from the solar beam.
        """
        self._dirflux, self._pa.extdirp, self._cx, self._cy, self._cz, \
        self._cxinv, self._cyinv, self._czinv, self._di, self._dj, self._dk, \
        self._ipdirect, self._delxd, self._delyd, self._xdomain, self._ydomain, \
        self._epss, self._epsz, self._uniformzlev = \
            pyshdom.core.make_direct(
                nstleg=self._nstleg,
                dirflux=self._dirflux,
                extdirp=self._pa.extdirp,
                npts=self._npts,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                ml=self._ml,
                nlegp=self._pa.nlegp,
                maxnmicro=self._pa.max_num_micro,
                maxpg=self._maxpg,
                npart=self._npart,
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
                phasewtp=self._pa.phasewtp,
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs
            )

    @property
    def num_iterations(self):
        return self._iters
