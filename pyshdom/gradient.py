"""
Contains the implementation of the Levis approximation to the Frechet derivatives which
does all of the hardwork of a) preparing uncertainties and b) the mechanistic
part of the evaluation of the cost function and the gradient.
To implement a slightly different cost function / gradient that still makes
use of the Levis approximation, you will need to a) create a new cost function
in src/shdomsub4.f UPDATE_COSTFUNCTION, b) make sure an uncertainty model is
compatible (uncertainties.py) and possibly, c) add a class that inherits from
LevisApproxGradient (e.g. LevisApproxGradientUncorrelated).
"""
import copy
import warnings

import numpy as np
import xarray as xr
import pandas as pd

import pyshdom.core
import pyshdom.parallel

class LevisApproxGradient:
    """
    A base class for different cost function evaluations that does
    the mechanical part of preparing the gradient and the
    inverse error covariance matrix.
    """
    def __init__(self, measurements, solvers, forward_sensors,
                 unknown_scatterers, parallel_solve_kwargs, gradient_kwargs,
                 uncertainty_kwargs):
        #TODO do type checks here.
        self.measurements = measurements
        self.solvers = solvers
        self.forward_sensors = forward_sensors
        self.unknown_scatterers = unknown_scatterers
        self.parallel_solve_kwargs = parallel_solve_kwargs
        self.gradient_kwargs = gradient_kwargs
        self.uncertainty_kwargs = uncertainty_kwargs
        self._rte_sensors = None
        self._sensor_mapping = None

        #These are variables that are not used in standard gradient calculations.
        #They were used to debug the radiance calculation in the gradient
        #function. They are still in place for tests but are just left
        #with these defaults. DO NOT CHANGE self._uselongrad unless you know
        #what you are doing.
        self._uselongrad = 'Q'
        self._longradiance = None
        # tautol is the subgrid interval size for the radiance integration.
        # this was added as an option for testing. But it is hardcoded
        # to the original value used in SHDOM radiance integration of 0.2.
        self._tautol = 0.2

        for name, instrument in self.measurements.items():
            if instrument['uncertainty_model'] is None:
                warnings.warn(
                    "No uncertainty model supplied for instrument '{}'. "
                    "Using pyshdom.uncertainties.NullUncertainty which is"
                    "equivalent to unweighted least squares.".format(name))
                self.measurements.add_uncertainty_model(
                    name,
                    pyshdom.uncertainties.NullUncertainty(self.gradient_kwargs['cost_function'])
                    )
            if instrument['uncertainty_model'].cost_function != self.gradient_kwargs['cost_function']:
                raise ValueError(
                    "Uncertainty model's assumed cost_function '{}' "
                    "is inconsistent with the one being used '{}'".format(
                        instrument['uncertainty_model'].cost_function,
                        self.gradient_kwargs['cost_function']
                        )
                    )
            for sensor in instrument['sensor_list']:
                instrument['uncertainty_model'].calculate_uncertainties(sensor)
                if self.uncertainty_kwargs['add_noise']:
                    self.measurements.add_noise(sensor)

    def _prep_gradient(self):

        # for solver in self.solvers.values():
        #     if solver._srctype != 'S':
        #         raise NotImplementedError(
        #             "Only Solar Source is supported for gradient calculations.")

        self.solvers.parallel_solve(**self.parallel_solve_kwargs)
        #does some preprocessing for calculating the sensitivity of a gridpoint's
        #solar source to the optical properties along the path to the sun.
        self.solvers.add_direct_beam_derivatives()
        #adds the _dext/_dleg/_dalb/_diphase etc to the solvers.
        self.solvers.add_microphysical_partial_derivatives(self.unknown_scatterers)
        #prepare the sensors for the fortran subroutine for calculating gradient.
        rte_sensors, sensor_mapping = self.forward_sensors.sort_sensors(
            self.solvers, self.measurements
            )
        self._rte_sensors = rte_sensors
        self._sensor_mapping = sensor_mapping
        mpi_comm = self.parallel_solve_kwargs['mpi_comm'] if 'mpi_comm' in self.parallel_solve_kwargs else None
        n_jobs = self.parallel_solve_kwargs['n_jobs'] if 'n_jobs' in self.parallel_solve_kwargs else None

        #The treatment of the gradient_kwargs is quite clumsy here as they are known in self
        #but are sent, instead of redefining gradient_fun to be self.levis_approximation_grad
        #WITH the kwargs set.
        outputs = pyshdom.parallel.parallel_gradient(
            self.solvers, rte_sensors, sensor_mapping, self.forward_sensors,
            gradient_fun=self.levis_approximation_grad,
            mpi_comm=mpi_comm,
            n_jobs=n_jobs, **self.gradient_kwargs
            )
        return outputs

    def levis_approximation_grad(self, rte_solver, sensor, cost_function='L2',
                                 indices_for_jacobian=None, exact_single_scatter=True):
        """
        Calculates the gradient of a cost function according to the Levis approximation to the Frechet
        derivatives of the RTE equation.

        Calls the main fortran subroutine in src/polarized/shdomsub4.f.

        Parameters
        ----------
        rte_solver : pyshdom.solver.RTE
            Should be solved, so that this it holds a solution to an RTE.
        sensor : xr.Dataset
            A sensor which should contain pixel-level uncertainties and measurement data
            as well as the ray & pixel geometry for calculating the forward model
            pixel values for evaluation of the cost function and its gradient.
        Returns
        -------
        loss: float64
            The value of the cost function accumulated over all pixels.
        gradient: np.array(shape=(rte_solver._nbpts, rte_solver.num_derivatives), dtype=np.float64)
            The gradient with respect to all parameters at every grid base point
        integrated_rays : xr.Dataset
            The forward model output used to evaluate the cost function against
            the measurements.
        """
        #This function could also be a method of solver.RTE just like
        #calculate_microphysical_partial_derivatives and calculate_direct_beam_derivative.
        #both of which are preparatory for the inverse problem.
        #However, this function doesn't modify the solver in place.
        #It is expected that if the API for pyshdom.core.levisapprox_gradient
        #changes it will be todo with cost function / etc changes rather than changes
        #in how solver.RTE is defined so we are happy accessing all of the 'private'
        #variables of solver.RTE here.
        if not isinstance(rte_solver, pyshdom.solver.RTE):
            raise TypeError(
                "`rte_solver` must be of type pyshdom.solver.RTE for this gradient"
                "calculation. "
            )

        if not isinstance(sensor, xr.Dataset):
            raise TypeError("`sensor` should be an xr.Dataset not "
                            "of type '{}''".format(type(sensor)))
        pyshdom.checks.check_hasdim(sensor, ray_mu='nrays', ray_phi='nrays',
                                    ray_x='nrays', ray_y='nrays', ray_z='nrays',
                                    stokes='stokes_index')

        #Update this code here if new cost functions are implemented in
        #src/shdomsub4.f UPDATE_COSTFUNCTION that require a larger number of scalar
        #or vector quantities (e.g. normalized cross correlation.)
        if cost_function in ('L2', 'LL'):
            cost_size = 1
            gradient_size = 1
        else:
            raise NotImplementedError("`cost_function` '{}' is not valid.".format(cost_function))

        camx = sensor['ray_x'].data
        camy = sensor['ray_y'].data
        camz = sensor['ray_z'].data
        cammu = sensor['ray_mu'].data
        camphi = sensor['ray_phi'].data

        total_pix = sensor.sizes['npixels']
        measurement_data = sensor['measurement_data'].data
        stokes_weights = sensor['stokes_weights'].data
        ray_weights = sensor['ray_weight'].data
        rays_per_pixel = sensor['rays_per_pixel'].data
        uncertainties = sensor['uncertainties'].data
        num_uncertainty = sensor['num_uncertainty'].size

        if indices_for_jacobian is None:
            jacobian = np.empty(
                (rte_solver._nstokes, rte_solver._num_derivatives, 1, 1),
                order='F',
                dtype=np.float32
            )
            num_jacobian_pts = 1
            jacobian_ptr = np.zeros(num_jacobian_pts)
            jacobian_flag = False
        else:
            jacobian_ptr = np.ravel_multi_index(
                (indices_for_jacobian),
                (rte_solver._pa.npx, rte_solver._pa.npy, rte_solver._pa.npz)
                ) + 1
            num_jacobian_pts = np.size(jacobian_ptr)
            jacobian = np.empty(
                (rte_solver._nstokes, rte_solver._num_derivatives, num_jacobian_pts, total_pix),
                order='F',
                dtype=np.float32
            )
            jacobian_flag = True
        if self._longradiance is None:
            self._longradiance = np.zeros((rte_solver._nstokes, rte_solver._npts), dtype=np.float32,
                                          order='F')
        gradient, loss, images, jacobian = pyshdom.core.levisapprox_gradient(
            camx=camx,
            camy=camy,
            camz=camz,
            cammu=cammu,
            camphi=camphi,
            npix=total_pix,
            costfunc=cost_function,
            ncost=cost_size,
            ngrad=gradient_size,
            nuncertainty=num_uncertainty,
            uncertainties=uncertainties,
            rays_per_pixel=rays_per_pixel,
            ray_weights=ray_weights,
            stokes_weights=stokes_weights,
            exact_single_scatter=exact_single_scatter,
            measurements=measurement_data,
            jacobian=jacobian,
            jacobianptr=jacobian_ptr,
            num_jacobian_pts=num_jacobian_pts,
            makejacobian=jacobian_flag,
            longradiance=self._longradiance,
            uselongrad=self._uselongrad,
            tautol=self._tautol,
            diphaseind=rte_solver._diphaseind,
            nstphase=rte_solver._nstphase,
            dpath=rte_solver._direct_derivative_path,
            dptr=rte_solver._direct_derivative_ptr,
            npx=rte_solver._pa.npx,
            npy=rte_solver._pa.npy,
            npz=rte_solver._pa.npz,
            delx=rte_solver._pa.delx,
            dely=rte_solver._pa.dely,
            xstart=rte_solver._pa.xstart,
            ystart=rte_solver._pa.ystart,
            zlevels=rte_solver._pa.zlevels,
            extdirp=rte_solver._pa.extdirp,
            uniformzlev=rte_solver._uniformzlev,
            partder=rte_solver._unknown_scatterer_indices,
            numder=rte_solver._num_derivatives,
            dext=rte_solver._dext,
            dalb=rte_solver._dalb,
            diphase=rte_solver._diphase,
            dleg=rte_solver._dleg,
            dphasetab=rte_solver._dphasetab,
            dnumphase=rte_solver._dnumphase,
            nscatangle=rte_solver._nscatangle,
            phasetab=rte_solver._phasetab,
            ylmsun=rte_solver._ylmsun,
            nstokes=rte_solver._nstokes,
            nstleg=rte_solver._nstleg,
            nx=rte_solver._nx,
            ny=rte_solver._ny,
            nz=rte_solver._nz,
            bcflag=rte_solver._bcflag,
            ipflag=rte_solver._ipflag,
            npts=rte_solver._npts,
            nbpts=rte_solver._nbpts,
            ncells=rte_solver._ncells,
            nbcells=rte_solver._nbcells,
            ml=rte_solver._ml,
            mm=rte_solver._mm,
            ncs=rte_solver._ncs,
            nlm=rte_solver._nlm,
            numphase=rte_solver._pa.numphase,
            nmu=rte_solver._nmu,
            nphi0max=rte_solver._nphi0max,
            nphi0=rte_solver._nphi0,
            maxnbc=rte_solver._maxnbc,
            ntoppts=rte_solver._ntoppts,
            nbotpts=rte_solver._nbotpts,
            nsfcpar=rte_solver._nsfcpar,
            gridptr=rte_solver._gridptr,
            neighptr=rte_solver._neighptr,
            treeptr=rte_solver._treeptr,
            shptr=rte_solver._shptr,
            bcptr=rte_solver._bcptr,
            cellflags=rte_solver._cellflags,
            iphase=rte_solver._iphase[:rte_solver._npts],
            deltam=rte_solver._deltam,
            solarflux=rte_solver._solarflux,
            solarmu=rte_solver._solarmu,
            solaraz=rte_solver._solaraz,
            gndtemp=rte_solver._gndtemp,
            gndalbedo=rte_solver._gndalbedo,
            skyrad=rte_solver._skyrad,
            waveno=rte_solver._waveno,
            wavelen=rte_solver.wavelength,
            mu=rte_solver._mu,
            phi=rte_solver._phi,
            wtdo=rte_solver._wtdo,
            xgrid=rte_solver._xgrid,
            ygrid=rte_solver._ygrid,
            zgrid=rte_solver._zgrid,
            gridpos=rte_solver._gridpos,
            sfcgridparms=rte_solver._sfcgridparms,
            bcrad=copy.deepcopy(rte_solver._bcrad),
            extinct=rte_solver._extinct[:rte_solver._npts],
            albedo=rte_solver._albedo[:rte_solver._npts],
            legen=rte_solver._legen,
            dirflux=rte_solver._dirflux[:rte_solver._npts],
            fluxes=rte_solver._fluxes,
            source=rte_solver._source,
            srctype=rte_solver._srctype,
            sfctype=rte_solver._sfctype,
            units=rte_solver._units,
            rshptr=rte_solver._rshptr,
            radiance=rte_solver._radiance,
            total_ext=rte_solver._total_ext[:rte_solver._npts],
            planck=rte_solver._planck[:rte_solver._npts]
        )

        integrated_rays = sensor.copy(deep=True)
        data = {}
        if rte_solver._nstokes == 1:
            data['I'] = (['npixels'], images[0])
        elif rte_solver._nstokes > 1:
            data['I'] = (['npixels'], images[0])
            data['Q'] = (['npixels'], images[1])
            data['U'] = (['npixels'], images[2])
        if rte_solver._nstokes == 4:
            data['V'] = (['npixels'], images[3])
        for key, val in data.items():
            integrated_rays[key] = val

        if not jacobian_flag:
            jacobian = None
        return gradient, loss, integrated_rays, jacobian



    def __call__(self):
        """
        A method to be overwritten in inheritance.
        """
        outputs = self._prep_gradient()
        return outputs, None, None

class LevisApproxGradientUncorrelated(LevisApproxGradient):
    """
    LevisApproxGradient that assumes different wavelengths
    are uncorrelated.

    This is the default method to use.
    """
    def __call__(self):

        loss, gradient, other_outputs = self._prep_gradient()
        #uncorrelated among the output of all (possibly parallel) workers.
        #if wanting to impose error correlations between wavelengths etc
        #then make a new class with a modified __call__ method.
        loss = np.sum(loss) / self.forward_sensors.nmeasurements
        gradient = np.sum(gradient, axis=-1) / self.forward_sensors.nmeasurements
        #turn gradient into a gridded dataset for use in project_gradient_to_state
        gradient_dataset = make_gradient_dataset(gradient, self.unknown_scatterers, self.solvers)
        if other_outputs:
            jacobian_dataset = make_jacobian_dataset(
                other_outputs[0], self.unknown_scatterers,
                self.gradient_kwargs['indices_for_jacobian'], self.solvers, self._rte_sensors
                )
        else:
            jacobian_dataset = None

        return loss, gradient_dataset, jacobian_dataset

def make_gradient_dataset(gradient, unknown_scatterers, solvers):
    """
    A utility function that forms an xr.Dataset for the gradient
    for downstream ease of use when postprocessing at the
    script level.

    Parameters
    ----------
    gradient : np.ndarray, shape=(nbpts, numder)
        The gradient of the cost function.
    unknown_scatterers : pyshdom.containers.UnknownScatterers
        Contains the information for defining the names and variables
        that derivatives have been calculated for.
    solvers : pyshdom.containers.SolversDict
        Contains the solver.RTE objects that were used to calculate the gradient.
        This is used to verify consistency between `unknown_scatterers`
        and `solvers` and to provide the RTE grid used.

    Returns
    -------
    gradient_dataset : xr.Dataset
        The dataset containing the gradient.

    Raises
    ------
    AssertionError
        If `solvers` and `unknown_scatterers` have inconsistent representations
        of which scattering species gradients have been calculated for.
    """
    derivative_names = []
    unknown_scatterer_names2 = []
    for name, values in unknown_scatterers.items():
        for variable_name in values['variable_name_list']:
            derivative_names.append(variable_name)
            unknown_scatterer_names2.append(name)
    unknown_scatterer_indices = list(solvers.values())[0]._unknown_scatterer_indices - 1
    unknown_scatterer_names = np.array(list(list(solvers.values())[0].medium.keys()))[unknown_scatterer_indices]
    assert np.all(unknown_scatterer_names == np.atleast_1d(unknown_scatterer_names2)), 'Two different ways of listing unknown scatterer names do not match.'
    derivative_index = pd.MultiIndex.from_arrays([unknown_scatterer_names, derivative_names], names=("scatterer_name", "variable_name"))

    grid = list(solvers.values())[0]._grid
    gradient_dataset = xr.Dataset(
                        data_vars={
                            'gradient': (['x', 'y', 'z', 'derivative_index'],
                                 gradient.reshape((grid.sizes['x'], grid.sizes['y'], grid.sizes['z'], -1)))
                        },
        coords={
            'x': grid.x,
            'y': grid.y,
            'z': grid.z,
            'derivative_index': derivative_index
        }
    )
    return gradient_dataset

def make_jacobian_dataset(jacobian_list, unknown_scatterers, indices_for_jacobian, solvers, rte_sensors):
    """
    A utility function that forms an xr.Dataset for the Frechet derivatives
    under the Levis approximation for downstream ease of use when postprocessing
    at the script level.

    Parameters
    ----------
    jacobian_list : List
        List of Jacobian arrays (possibly from parallel workers).
    unknown_scatterers : pyshdom.containers.UnknownScatterers
        Contains the information for defining the names and variables
        that derivatives have been calculated for.
    indices_for_jacobian : Tuple
        The indices for which grid points to save the Jacobian for.
    solvers : pyshdom.containers.SolversDict
        Contains the solver.RTE objects that were used to calculate the gradient.
        This is only used to verify consistency between `unknown_scatterers`
        and `solvers`.
    rte_sensors : OrderedDict
        Sensors grouped by key in `solvers`. Contains the information on how
        to split and merge jacobian_list.

    Returns
    -------
    jacobian_dataset : xr.Dataset
        The dataset containing the Frechet derivatives.

    Raises
    ------
    AssertionError
        If `solvers` and `unknown_scatterers` have inconsistent representations
        of which scattering species gradients have been calculated for.
    """
    merged_jacobian = np.concatenate(jacobian_list, axis=-1)
    split_indices = []
    for rte_sensor in rte_sensors.values():
        size = rte_sensor.sizes['npixels']
        split_indices.append(size)
    split_indices = np.cumsum(split_indices) #TODO verify that no +1 needs to be added etc.
    split_jacobian = np.split(merged_jacobian, split_indices, axis=-1)[:-1] #TODO verify the [:-1]

    grid = list(solvers.values())[0]._grid
    grid_index = pd.MultiIndex.from_arrays([grid.x.data[indices_for_jacobian[0]],
                                            grid.y.data[indices_for_jacobian[1]],
                                            grid.z.data[indices_for_jacobian[2]]],
                                           names=("x", "y", "z"))

    derivative_names = []
    unknown_scatterer_names2 = []
    for name, values in unknown_scatterers.items():
        for variable_name in values['variable_name_list']:
            derivative_names.append(variable_name)
            unknown_scatterer_names2.append(name)
    unknown_scatterer_indices = list(solvers.values())[0]._unknown_scatterer_indices - 1
    unknown_scatterer_names = np.array(list(list(solvers.values())[0].medium.keys()))[unknown_scatterer_indices]
    assert np.all(unknown_scatterer_names == np.atleast_1d(unknown_scatterer_names2)), 'Two different ways of listing unknown scatterer names do not match.'
    derivative_index = pd.MultiIndex.from_arrays([unknown_scatterer_names, derivative_names], names=("scatterer_name", "variable_name"))

    jacobian_dataset = xr.Dataset(
                                data_vars ={
                                'jacobian_{:1.3f}'.format(wavelength): (['nstokes', 'derivative_index', 'grid_index', 'npixels_{:1.3f}'.format(wavelength)], jacobian)
                                    for wavelength, jacobian in zip(solvers.keys(), split_jacobian)
                                },
                    coords={
                        'grid_index': grid_index,
                        'derivative_index': derivative_index
                    }
    )
    return jacobian_dataset

# def grad_l2_old(rte_solver, sensor, exact_single_scatter=True,
#     jacobian_flag=False, stokes_weights=[1.0,1.0,1.0,0.0]):
#     """
#     The core l2 gradient method.
#
#     Parameters
#     ----------
#     rte_solver: pyshdom.RteSolver
#         A solver with all the associated parameters and the solution to the RTE
#     projection: pyshdom.Projection
#         A projection model which specified the position and direction of each and every pixel
#     pixels: np.array(shape=(projection.npix), dtype=np.float32)
#         The acquired pixels driving the error and optimization.
#     uncertainties: np.array(shape=(projection.npix), dtype=np.float32)
#         The pixel uncertainties.
#
#     Returns
#     -------
#     gradient: np.array(shape=(rte_solver._nbpts, self.num_derivatives), dtype=np.float64)
#         The gradient with respect to all parameters at every grid base point
#     loss: float64
#         The total loss accumulated over all pixels
#     images: np.array(shape=(rte_solver._nstokes, projection.npix), dtype=np.float32)
#         The rendered (synthetic) images.
#     """
#     camx = sensor['ray_x'].data
#     camy = sensor['ray_y'].data
#     camz = sensor['ray_z'].data
#     cammu = sensor['ray_mu'].data
#     camphi = sensor['ray_phi'].data
#     #TODO
#     #Some checks on the dimensions: this kind of thing.
#     #assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
#     #nrays = sensor.sizes['nrays']
#     total_pix = sensor.sizes['npixels']
#     measurement_data = sensor['measurement_data'].data
#     stokes_weights = sensor['stokes_weights'].data
#     ray_weights = sensor['ray_weight'].data
#     rays_per_pixel = sensor['rays_per_pixel'].data
#     uncertainties = sensor['uncertainties'].data
#     #TODO fix jacobian.
#     if jacobian_flag:
#         #maximum size is hard coded - possible source of seg faults.
#         #(8*MAX_DOMAIN_LENGTH)**2 based on the beam from sensor to voxel and from voxel to sun.
#         largest_dim = int(64*(rte_solver._pa.npx**2+rte_solver._pa.npy**2+rte_solver._pa.npz**2))
#         jacobian = np.empty((rte_solver._nstokes,self.num_derivatives,total_pix*largest_dim),order='F',dtype=np.float32)
#         jacobian_ptr = np.empty((2, total_pix*largest_dim), order='F', dtype=np.int32)
#     else:
#         jacobian = np.empty(
#             (rte_solver._nstokes, rte_solver._num_derivatives, 1),
#             order='F',
#             dtype=np.float32)
#         jacobian_ptr = np.empty((2, 1), order='F', dtype=np.int32)
#
#     gradient, loss, images, jacobian, jacobian_ptr, counter = pyshdom.core.gradient_l2_old(
#         uncertainties=uncertainties,
#         #rays_per_pixel=rays_per_pixel,
#         #ray_weights=ray_weights,
#         weights=np.array(stokes_weights)[:rte_solver._nstokes],
#         exact_single_scatter=exact_single_scatter,
#         nstphase=rte_solver._nstphase,
#         dpath=rte_solver._direct_derivative_path,
#         dptr=rte_solver._direct_derivative_ptr,
#         npx=rte_solver._pa.npx,
#         npy=rte_solver._pa.npy,
#         npz=rte_solver._pa.npz,
#         delx=rte_solver._pa.delx,
#         dely=rte_solver._pa.dely,
#         xstart=rte_solver._pa.xstart,
#         ystart=rte_solver._pa.ystart,
#         zlevels=rte_solver._pa.zlevels,
#         extdirp=rte_solver._pa.extdirp,
#         uniformzlev=rte_solver._uniformzlev,
#         partder=rte_solver._unknown_scatterer_indices,
#         numder=rte_solver._num_derivatives,
#         dext=rte_solver._dext,
#         dalb=rte_solver._dalb,
#         diphase=rte_solver._diphase,
#         dleg=rte_solver._dleg,
#         dphasetab=rte_solver._dphasetab,
#         dnumphase=rte_solver._dnumphase,
#         nscatangle=rte_solver._nscatangle,
#         phasetab=rte_solver._phasetab,
#         ylmsun=rte_solver._ylmsun,
#         nstokes=rte_solver._nstokes,
#         nstleg=rte_solver._nstleg,
#         nx=rte_solver._nx,
#         ny=rte_solver._ny,
#         nz=rte_solver._nz,
#         bcflag=rte_solver._bcflag,
#         ipflag=rte_solver._ipflag,
#         npts=rte_solver._npts,
#         nbpts=rte_solver._nbpts,
#         ncells=rte_solver._ncells,
#         nbcells=rte_solver._nbcells,
#         ml=rte_solver._ml,
#         mm=rte_solver._mm,
#         ncs=rte_solver._ncs,
#         nlm=rte_solver._nlm,
#         numphase=rte_solver._pa.numphase,
#         nmu=rte_solver._nmu,
#         nphi0max=rte_solver._nphi0max,
#         nphi0=rte_solver._nphi0,
#         maxnbc=rte_solver._maxnbc,
#         ntoppts=rte_solver._ntoppts,
#         nbotpts=rte_solver._nbotpts,
#         nsfcpar=rte_solver._nsfcpar,
#         gridptr=rte_solver._gridptr,
#         neighptr=rte_solver._neighptr,
#         treeptr=rte_solver._treeptr,
#         shptr=rte_solver._shptr,
#         bcptr=rte_solver._bcptr,
#         cellflags=rte_solver._cellflags,
#         iphase=rte_solver._iphase[:rte_solver._npts],
#         deltam=rte_solver._deltam,
#         solarflux=rte_solver._solarflux,
#         solarmu=rte_solver._solarmu,
#         solaraz=rte_solver._solaraz,
#         gndtemp=rte_solver._gndtemp,
#         gndalbedo=rte_solver._gndalbedo,
#         skyrad=rte_solver._skyrad,
#         waveno=rte_solver._waveno,
#         wavelen=rte_solver.wavelength,
#         mu=rte_solver._mu,
#         phi=rte_solver._phi,
#         wtdo=rte_solver._wtdo,
#         xgrid=rte_solver._xgrid,
#         ygrid=rte_solver._ygrid,
#         zgrid=rte_solver._zgrid,
#         gridpos=rte_solver._gridpos,
#         sfcgridparms=rte_solver._sfcgridparms,
#         bcrad=copy.deepcopy(rte_solver._bcrad),
#         extinct=rte_solver._extinct[:rte_solver._npts],
#         albedo=rte_solver._albedo[:rte_solver._npts],
#         legen=rte_solver._legen,
#         dirflux=rte_solver._dirflux[:rte_solver._npts],
#         fluxes=rte_solver._fluxes,
#         source=rte_solver._source,
#         camx=camx,
#         camy=camy,
#         camz=camz,
#         cammu=cammu,
#         camphi=camphi,
#         npix=len(camx),
#         srctype=rte_solver._srctype,
#         sfctype=rte_solver._sfctype,
#         units=rte_solver._units,
#         measurements=measurement_data,
#         rshptr=rte_solver._rshptr,
#         radiance=rte_solver._radiance,
#         total_ext=rte_solver._total_ext[:rte_solver._npts],
#         jacobian=jacobian,
#         jacobianptr=jacobian_ptr,
#         makejacobian=jacobian_flag
#     )
#     jacobian = jacobian[:, :, :counter-1]
#     jacobian_ptr = jacobian_ptr[:, :counter-1]
#     if jacobian_flag:
#         return gradient, loss, images, jacobian, jacobian_ptr, np.array([counter-1])
#     else:
#         integrated_rays = sensor.copy(deep=True)
#         data = {}
#         if rte_solver._nstokes == 1:
#             data['I'] = (['npixels'], images[0])
#         elif rte_solver._nstokes > 1:
#             data['I'] = (['npixels'], images[0])
#             data['Q'] = (['npixels'], images[1])
#             data['U'] = (['npixels'], images[2])
#         if rte_solver._nstokes == 4:
#             data['V'] = (['npixels'], images[3])
#         for key, val in data.items():
#             integrated_rays[key] = val
#
#         return gradient, loss, integrated_rays
