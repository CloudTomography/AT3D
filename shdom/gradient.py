import numpy as np
import xarray as xr
from collections import OrderedDict
import shdom
from joblib import Parallel, delayed
import inspect
import copy
import pandas as pd

def create_derivative_tables(unknown_scatterers):
    """
    TODO
    """
    partial_derivative_tables = OrderedDict()
    for key in unknown_scatterers[0][1].keys(): #the wavelengths in the poly tables.
        partial_derivative_dict = OrderedDict()
        for unknown_scatterer, tables, variable_names in unknown_scatterers:
            table = tables[key]
            variable_names = np.atleast_1d(variable_names)
            #For each unknown_scatterer compute the tables of derivatives
            valid_microphysics_coords = {name:table[name] for name in table.coords if name not in ('table_index', 'stokes_index')}
            derivatives_tables = OrderedDict()
            for variable_name in variable_names:

                if variable_name in valid_microphysics_coords.keys():
                    differentiated = table.differentiate(coord=variable_name)

                elif variable_name == 'extinction':
                    differentiated = table.copy(data={
                        'extinction': np.ones(table.extinction.shape),
                        'legcoef': np.zeros(table.legcoef.shape),
                        'ssalb': np.zeros(table.ssalb.shape)
                    })
                elif variable_name == 'ssalb':
                    differentiated = table.copy(data={
                        'extinction': np.zeros(table.extinction.shape),
                        'legcoef': np.zeros(table.legcoef.shape),
                        'ssalb': np.ones(table.ssalb.shape)
                    })
                elif 'legendre_' in variable_name:
                    leg_index = int(variable_name[len('legendre_X_'):])
                    stokes_index = int(variable_name[len('legendre_')])
                    legcoef = np.zeros(table.legcoef.shape)
                    legcoef[stokes_index,leg_index,...] = 1.0
                    differentiated = table.copy(data={
                        'extinction': np.zeros(table.extinction.shape),
                        'legcoef': legcoef,
                        'ssalb': np.zeros(table.ssalb.shape)
                    })
                elif variable_name == 'density':
                    differentiated = table.copy(data={
                        'extinction': table.extinction,
                        'legcoef': np.zeros(table.legcoef.shape),
                        'ssalb': np.zeros(table.ssalb.shape)
                    })
                else:
                    raise ValueError('Invalid unknown scatterer name', variable_name)

                derivatives_tables[variable_name] = differentiated
            partial_derivative_dict[unknown_scatterer] = derivatives_tables
        partial_derivative_tables[key] = partial_derivative_dict
    return partial_derivative_tables

# def get_derivatives(solvers, all_derivative_tables):
#     """
#     TODO
#     Calculates partial derivatives on the rte_grid.
#     The solvers are modified in place to contain the partial_derivatives
#
#     This requires density to be defined in scatterer.
#     """
#     #TODO more detailed checks that derivative tables match solver scatterers.
#     #At the moment it is assumed that they are defined consistently.
#     for key, rte_solver in solvers.items():
#
#         rte_solver._precompute_phase()
#         solver_derivative_table = all_derivative_tables[key]
#         num_derivatives = sum([len(scatterer_derivative_table.values()) for name,scatterer_derivative_table in \
#                               solver_derivative_table.items()
#                               ])
#         rte_solver._num_derivatives = np.array(num_derivatives, dtype=np.int32)
#         unknown_scatterer_indices = []
#
#         dext = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.float32)
#         dalb = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.float32)
#         diphase = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.int32)
#
#         #one loop through to find max_legendre and unkonwn_scatterer_indices
#         max_legendre = []
#         i = 0
#         for name,scatterer_derivative_table in solver_derivative_table.items():
#             scatterer = rte_solver.medium[name]
#             for variable_derivative_table in scatterer_derivative_table.values():
#                 derivative_on_grid = shdom.medium.table_to_grid(scatterer,variable_derivative_table, inverse_mode=True)
#                 max_legendre.append(derivative_on_grid.sizes['legendre_index'])
#                 unknown_scatterer_indices.append(i+1)
#             i += 1
#         max_legendre = max(max_legendre)
#         rte_solver._unknown_scatterer_indices = np.array(unknown_scatterer_indices).astype(np.int32)
#
#         #second loop to assign everything else.
#         padded_legcoefs = []
#         count = 0
#         for name,scatterer_derivative_table in solver_derivative_table.items():
#             scatterer = rte_solver.medium[name]
#             for variable_derivative_table in scatterer_derivative_table.values():
#                 derivative_on_grid = shdom.medium.table_to_grid(scatterer,variable_derivative_table, inverse_mode=True)
#
#                 dext[:, count] = derivative_on_grid.extinction.data.ravel()
#                 dalb[:, count] = derivative_on_grid.ssalb.data.ravel()
#                 diphase[:, count] = derivative_on_grid.table_index.data.ravel() + diphase.max()
#
#                 padded_legcoefs.append(derivative_on_grid.legcoef.pad(
#                     {'legendre_index': (0, max_legendre - derivative_on_grid.legcoef.sizes['legendre_index'])},
#                     constant_values=0.0
#                 ))
#
#                 count += 1
#
#         #COPIED FROM LINE 427 OF solver.py
#         #In regions which are not covered by any optical scatterer they have an iphasep of 0.
#         #In original SHDOM these would be pointed to the rayleigh phase function (which is always included
#         #in the legendre table even if there is no rayleigh extinction.)
#         #Here, instead we set them to whatever the first phase function is.
#         #An arbitrary valid choice can be made as the contribution from these grid points is zero.
#         diphase[np.where(diphase == 0)] = 1
#
#         # Concatenate all legendre tables into one table
#         legendre_table = xr.concat(padded_legcoefs, dim='table_index')
#         dnumphase = legendre_table.sizes['table_index']
#         dleg = legendre_table.data
#         #TODO make sure the shaping of dleg and consistent with legenp/._legen. (._nleg and numphase)
#         #For a single scatterer/variable dnumphase == numphase, but this is not always the case.
#
#         # zero the first term of the first component of the phase function
#         # gradient. Pre-scale the legendre moments by 1/(2*l+1) which
#         # is done in the forward problem in TRILIN_INTERP_PROP
#         scaling_factor = np.atleast_3d(np.array([2.0*i+1.0 for i in range(0,rte_solver._nleg+1)]))
#         dleg[0,0,:] = 0.0
#         dleg = dleg[:rte_solver._nstleg] / scaling_factor
#
#         dphasetab = shdom.core.precompute_phase_check_grad(
#                                                      negcheck=False,
#                                                      nstphase=rte_solver._nstphase,
#                                                      nstleg=rte_solver._nstleg,
#                                                      nscatangle=rte_solver._nscatangle,
#                                                      nstokes=rte_solver._nstokes,
#                                                      dnumphase=dnumphase,
#                                                      ml=rte_solver._ml,
#                                                      nlm=rte_solver._nlm,
#                                                      nleg=rte_solver._nleg,
#                                                      dleg=dleg,
#                                                      deltam=rte_solver._deltam
#                                                      )
#
#         rte_solver._dphasetab, rte_solver._dext, rte_solver._dalb, rte_solver._diphase, rte_solver._dleg, rte_solver._dnumphase = \
#         dphasetab, dext, dalb, diphase, dleg, dnumphase
#
#     return solvers


def grad_l2(rte_solver, sensor, exact_single_scatter=True):
    """
    The core l2 gradient method.

    Parameters
    ----------
    rte_solver: shdom.RteSolver
        A solver with all the associated parameters and the solution to the RTE
    projection: shdom.Projection
        A projection model which specified the position and direction of each and every pixel
    pixels: np.array(shape=(projection.npix), dtype=np.float32)
        The acquired pixels driving the error and optimization.
    uncertainties: np.array(shape=(projection.npix), dtype=np.float32)
        The pixel uncertainties.

    Returns
    -------
    gradient: np.array(shape=(rte_solver._nbpts, self.num_derivatives), dtype=np.float64)
        The gradient with respect to all parameters at every grid base point
    loss: float64
        The total loss accumulated over all pixels
    images: np.array(shape=(rte_solver._nstokes, projection.npix), dtype=np.float32)
        The rendered (synthetic) images.
    """
    camx = sensor['ray_x'].data
    camy = sensor['ray_y'].data
    camz = sensor['ray_z'].data
    cammu = sensor['ray_mu'].data
    camphi = sensor['ray_phi'].data
    #TODO
    #Some checks on the dimensions: this kind of thing.
    #assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
    #nrays = sensor.sizes['nrays']
    total_pix = sensor.sizes['npixels']
    measurement_data = sensor['measurement_data'].data
    stokes_weights = sensor['stokes_weights'].data
    ray_weights = sensor['ray_weight'].data
    rays_per_pixel = sensor['rays_per_pixel'].data
    uncertainties = sensor['uncertainties'].data

    #unused jacobian stuff.
    jacobian = np.empty((rte_solver._nstokes,rte_solver._num_derivatives,1,1),order='F',dtype=np.float32)
    num_jacobian_pts = 1
    jacobian_ptr = np.zeros(num_jacobian_pts)
    jacobian_flag=False

    gradient, loss, images, jacobian = shdom.core.gradient_l2(
        uncertainties=uncertainties,
        rays_per_pixel=rays_per_pixel,
        ray_weights=ray_weights,
        stokes_weights=stokes_weights,
        exact_single_scatter=exact_single_scatter,
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
        wavelen=rte_solver._wavelen,
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
        camx=camx,
        camy=camy,
        camz=camz,
        cammu=cammu,
        camphi=camphi,
        npix=total_pix,
        srctype=rte_solver._srctype,
        sfctype=rte_solver._sfctype,
        units=rte_solver._units,
        measurements=measurement_data,
        rshptr=rte_solver._rshptr,
        radiance=rte_solver._radiance,
        total_ext=rte_solver._total_ext[:rte_solver._npts],
        jacobian=jacobian,
        jacobianptr=jacobian_ptr,
        num_jacobian_pts=num_jacobian_pts,
        makejacobian=jacobian_flag
    )

    integrated_rays = sensor.copy(deep=True)
    data = {}
    if rte_solver._nstokes == 1:
        data['I'] = (['npixels'],images[0])
    elif rte_solver._nstokes > 1:
        data['I'] = (['npixels'],images[0])
        data['Q'] = (['npixels'],images[1])
        data['U'] = (['npixels'],images[2])
    if rte_solver._nstokes == 4:
        data['V'] = (['npixels'],images[3])
    for key, val in data.items():
        integrated_rays[key] = val

    return gradient, loss, integrated_rays

def jacobian(rte_solver, sensor, indices_for_jacobian, exact_single_scatter=True):
    """
    The core l2 gradient method.

    Parameters
    ----------
    rte_solver: shdom.RteSolver
        A solver with all the associated parameters and the solution to the RTE
    projection: shdom.Projection
        A projection model which specified the position and direction of each and every pixel
    pixels: np.array(shape=(projection.npix), dtype=np.float32)
        The acquired pixels driving the error and optimization.
    uncertainties: np.array(shape=(projection.npix), dtype=np.float32)
        The pixel uncertainties.

    Returns
    -------
    gradient: np.array(shape=(rte_solver._nbpts, self.num_derivatives), dtype=np.float64)
        The gradient with respect to all parameters at every grid base point
    loss: float64
        The total loss accumulated over all pixels
    images: np.array(shape=(rte_solver._nstokes, projection.npix), dtype=np.float32)
        The rendered (synthetic) images.
    """
    camx = sensor['ray_x'].data
    camy = sensor['ray_y'].data
    camz = sensor['ray_z'].data
    cammu = sensor['ray_mu'].data
    camphi = sensor['ray_phi'].data
    #TODO
    #Some checks on the dimensions: this kind of thing.
    #assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
    #nrays = sensor.sizes['nrays']
    total_pix = sensor.sizes['npixels']
    measurement_data = sensor['measurement_data'].data
    stokes_weights = sensor['stokes_weights'].data
    ray_weights = sensor['ray_weight'].data
    rays_per_pixel = sensor['rays_per_pixel'].data
    uncertainties = sensor['uncertainties'].data

    jacobian_ptrs = np.ravel_multi_index((indices_for_jacobian),(rte_solver._pa.npx,
                                        rte_solver._pa.npy,
                                        rte_solver._pa.npz)) + 1
    num_jacobian_pts = np.size(jacobian_ptrs)
    jacobian = np.empty((rte_solver._nstokes, rte_solver._num_derivatives, num_jacobian_pts, total_pix), order='F', dtype=np.float32)
    jacobian_flag = True
    gradient, loss, images, jacobian = shdom.core.gradient_l2(
        uncertainties=uncertainties,
        rays_per_pixel=rays_per_pixel,
        ray_weights=ray_weights,
        stokes_weights=stokes_weights,
        exact_single_scatter=exact_single_scatter,
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
        wavelen=rte_solver._wavelen,
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
        camx=camx,
        camy=camy,
        camz=camz,
        cammu=cammu,
        camphi=camphi,
        npix=total_pix,
        srctype=rte_solver._srctype,
        sfctype=rte_solver._sfctype,
        units=rte_solver._units,
        measurements=measurement_data,
        rshptr=rte_solver._rshptr,
        radiance=rte_solver._radiance,
        total_ext=rte_solver._total_ext[:rte_solver._npts],
        jacobian=jacobian,
        jacobianptr=jacobian_ptrs,
        num_jacobian_pts=num_jacobian_pts,
        makejacobian=jacobian_flag
    )

    integrated_rays = sensor.copy(deep=True)
    data = {}
    if rte_solver._nstokes == 1:
        data['I'] = (['npixels'],images[0])
    elif rte_solver._nstokes > 1:
        data['I'] = (['npixels'],images[0])
        data['Q'] = (['npixels'],images[1])
        data['U'] = (['npixels'],images[2])
    if rte_solver._nstokes == 4:
        data['V'] = (['npixels'],images[3])
    for key, val in data.items():
        integrated_rays[key] = val

    return gradient, loss, integrated_rays, jacobian

def grad_l2_old(rte_solver, sensor, exact_single_scatter=True,
    jacobian_flag=False, stokes_weights=[1.0,1.0,1.0,0.0]):
    """
    The core l2 gradient method.

    Parameters
    ----------
    rte_solver: shdom.RteSolver
        A solver with all the associated parameters and the solution to the RTE
    projection: shdom.Projection
        A projection model which specified the position and direction of each and every pixel
    pixels: np.array(shape=(projection.npix), dtype=np.float32)
        The acquired pixels driving the error and optimization.
    uncertainties: np.array(shape=(projection.npix), dtype=np.float32)
        The pixel uncertainties.

    Returns
    -------
    gradient: np.array(shape=(rte_solver._nbpts, self.num_derivatives), dtype=np.float64)
        The gradient with respect to all parameters at every grid base point
    loss: float64
        The total loss accumulated over all pixels
    images: np.array(shape=(rte_solver._nstokes, projection.npix), dtype=np.float32)
        The rendered (synthetic) images.
    """
    camx = sensor['ray_x'].data
    camy = sensor['ray_y'].data
    camz = sensor['ray_z'].data
    cammu = sensor['ray_mu'].data
    camphi = sensor['ray_phi'].data
    #TODO
    #Some checks on the dimensions: this kind of thing.
    #assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
    #nrays = sensor.sizes['nrays']
    total_pix = sensor.sizes['npixels']
    measurement_data = sensor['measurement_data'].data
    stokes_weights = sensor['stokes_weights'].data
    ray_weights = sensor['ray_weight'].data
    rays_per_pixel = sensor['rays_per_pixel'].data
    uncertainties = sensor['uncertainties'].data
    #TODO fix jacobian.
    if jacobian_flag:
        #maximum size is hard coded - possible source of seg faults.
        #(8*MAX_DOMAIN_LENGTH)**2 based on the beam from sensor to voxel and from voxel to sun.
        largest_dim = int(64*(rte_solver._pa.npx**2+rte_solver._pa.npy**2+rte_solver._pa.npz**2))
        jacobian = np.empty((rte_solver._nstokes,self.num_derivatives,total_pix*largest_dim),order='F',dtype=np.float32)
        jacobian_ptr = np.empty((2,total_pix*largest_dim),order='F',dtype=np.int32)
    else:
        jacobian = np.empty((rte_solver._nstokes,rte_solver._num_derivatives,1),order='F',dtype=np.float32)
        jacobian_ptr = np.empty((2,1),order='F',dtype=np.int32)

    gradient, loss, images, jacobian, jacobian_ptr, counter = shdom.core.gradient_l2_old(
        uncertainties=uncertainties,
        #rays_per_pixel=rays_per_pixel,
        #ray_weights=ray_weights,
        weights=np.array(stokes_weights)[:rte_solver._nstokes],
        exact_single_scatter=exact_single_scatter,
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
        wavelen=rte_solver._wavelen,
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
        camx=camx,
        camy=camy,
        camz=camz,
        cammu=cammu,
        camphi=camphi,
        npix=len(camx),
        srctype=rte_solver._srctype,
        sfctype=rte_solver._sfctype,
        units=rte_solver._units,
        measurements=measurement_data,
        rshptr=rte_solver._rshptr,
        radiance=rte_solver._radiance,
        total_ext=rte_solver._total_ext[:rte_solver._npts],
        jacobian=jacobian,
        jacobianptr=jacobian_ptr,
        makejacobian=jacobian_flag
    )
    jacobian = jacobian[:,:,:counter-1]
    jacobian_ptr = jacobian_ptr[:,:counter-1]
    if jacobian_flag:
        return gradient, loss, images, jacobian, jacobian_ptr, np.array([counter-1])
    else:
        integrated_rays = sensor.copy(deep=True)
        data = {}
        if rte_solver._nstokes == 1:
            data['I'] = (['npixels'],images[0])
        elif rte_solver._nstokes > 1:
            data['I'] = (['npixels'],images[0])
            data['Q'] = (['npixels'],images[1])
            data['U'] = (['npixels'],images[2])
        if rte_solver._nstokes == 4:
            data['V'] = (['npixels'],images[3])
        for key, val in data.items():
            integrated_rays[key] = val

        return gradient, loss, integrated_rays


def parallel_gradient(solvers, rte_sensors, sensor_mappings, forward_sensors, gradient_fun=grad_l2,
                     mpi_comm=None, n_jobs=1, **kwargs):
    """
    TODO
    """
    #organize **kwargs safely.
    grad_kwargs = {}
    grad_params = inspect.signature(gradient_fun).parameters
    for name,value in kwargs.items():
        if name in grad_params.keys():
            if grad_params[name].kind in (grad_params[name].POSITIONAL_OR_KEYWORD, grad_params[name].KEYWORD_ONLY,
                    grad_params[name].VAR_KEYWORD):
                grad_kwargs[name] = value
            else:
                warnings.warn("kwarg '{}' passed to shdom.gradient.parallel_gradient is unused by \
                                gradient_fun '{}''".format(name, gradient_fun.__name__))
        else:
            warnings.warn("kwarg '{}' passed to shdom.gradient.parallel_gradient is unused by \
                            gradient_fun '{}''".format(name, gradient_fun.__name__))

    if mpi_comm is not None:
        raise NotImplementedError

    else:
        if n_jobs == 1 or n_jobs >= forward_sensors.npixels:
            out = [gradient_fun(solvers[key],rte_sensors[key], **grad_kwargs) for key in solvers]
            keys = list(solvers.keys())
        else:
            #decide on the division of n_jobs among solvers based on total number of rays.
            keys, ray_start_end, pixel_start_end = shdom.util.subdivide_raytrace_jobs(rte_sensors, n_jobs)

            out = Parallel(n_jobs=n_jobs, backend='threading')(
                delayed(gradient_fun, check_pickle=False)(solvers[key],rte_sensors[key].sel(nrays=slice(ray_start,ray_end),
                                        npixels=slice(pix_start, pix_end)),**grad_kwargs)
                                for key, (ray_start,ray_end),(pix_start,pix_end) in zip(keys, ray_start_end, pixel_start_end))

        gradient = np.stack([i[0] for i in out],axis=-1)
        loss = np.array([i[1] for i in out])
        forward_model_output = [i[2] for i in out]
        #modify forward sensors in place to contain updated forward model estimates.
        forward_sensors.add_measurements_inverse(sensor_mappings, forward_model_output, keys)

    #special treatment of jacobian out.
    if gradient_fun == shdom.gradient.jacobian:
        jacobian_list = [i[-1] for i in out]
        return loss, gradient, jacobian_list
    else:
        return loss, gradient

def levis_approx_jacobian(measurements, solvers, forward_sensors, unknown_scatterers, indices_for_jacobian, n_jobs=1,mpi_comm=None,verbose=False,
                            maxiter=100, init_solution=True, setup_grid=True,
                            exact_single_scatter=True):

    #note this division of solving and 'raytracing' is not optimal for distributed memory parallelization
    #where tasks take varying amounts of time.
    for solver in solvers.values():
        if solver._srctype != 'S':
            raise NotImplementedError('Only Solar Source is supported for gradient calculations.')

    solvers.parallel_solve(n_jobs=n_jobs, mpi_comm=mpi_comm,verbose=verbose,maxiter=maxiter,
                            init_solution=init_solution, setup_grid=setup_grid)

    #These are called after the solution because they require at least _init_solution to be run.
    if not hasattr(list(solvers.values())[0], '_direct_derivative_path'):
        #only needs to be called once.
        solvers.add_direct_beam_derivatives()

    solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method,
                                                unknown_scatterers.table_data) #adds the _dext/_dleg/_dalb/_diphase etc to the solvers.
    #prepare the sensors for the fortran subroutine for calculating gradient.
    rte_sensors, sensor_mapping = forward_sensors.sort_sensors(solvers, measurements)
    loss, gradient, jacobian = parallel_gradient(solvers, rte_sensors, sensor_mapping, forward_sensors,
                        gradient_fun=shdom.gradient.jacobian,
                     mpi_comm=mpi_comm, n_jobs=n_jobs, exact_single_scatter=exact_single_scatter,indices_for_jacobian=indices_for_jacobian)

    #uncorrelated l2.
    loss = np.sum(loss) / forward_sensors.nmeasurements
    gradient =  np.sum(gradient, axis=-1)/ forward_sensors.nmeasurements

    #turn gradient into a gridded dataset for use in project_gradient_to_state
    gradient_dataset = make_gradient_dataset(gradient, unknown_scatterers, solvers)
    #turn jacobian into a dataset (minimal postprocessing)
    jacobian_dataset = make_jacobian_dataset(jacobian, unknown_scatterers, indices_for_jacobian, solvers,rte_sensors)

    return np.array(loss), gradient_dataset, jacobian_dataset


def levis_approx_uncorrelated_l2(measurements, solvers, forward_sensors, unknown_scatterers, n_jobs=1,
                                 mpi_comm=None,verbose=False, maxiter=100, init_solution=True,
                                 exact_single_scatter=True, setup_grid=True):
    """TODO"""
    #note this division of solving and 'raytracing' is not optimal for distributed memory parallelization
    #where tasks take varying amounts of time.
    for solver in solvers.values():
        if solver._srctype != 'S':
            raise NotImplementedError('Only Solar Source is supported for gradient calculations.')

    solvers.parallel_solve(n_jobs=n_jobs, mpi_comm=mpi_comm,verbose=verbose,maxiter=maxiter,
                            init_solution=init_solution, setup_grid=setup_grid)

    #These are called after the solution because they require at least _init_solution to be run.
    if not hasattr(list(solvers.values())[0], '_direct_derivative_path'):
        #only needs to be called once.
        solvers.add_direct_beam_derivatives()

    solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method,
                                                unknown_scatterers.table_data) #adds the _dext/_dleg/_dalb/_diphase etc to the solvers.

    #prepare the sensors for the fortran subroutine for calculating gradient.
    rte_sensors, sensor_mapping = forward_sensors.sort_sensors(solvers, measurements)

    loss, gradient = parallel_gradient(solvers, rte_sensors, sensor_mapping, forward_sensors, gradient_fun=grad_l2,
                     mpi_comm=mpi_comm, n_jobs=n_jobs, exact_single_scatter=exact_single_scatter)

    #uncorrelated l2.
    loss = np.sum(loss)/ forward_sensors.nmeasurements
    gradient =  np.sum(gradient, axis=-1)/ forward_sensors.nmeasurements

    #turn gradient into a gridded dataset for use in project_gradient_to_state
    gradient_dataset = make_gradient_dataset(gradient, unknown_scatterers, solvers)

    return np.array(loss), gradient_dataset


def make_gradient_dataset(gradient, unknown_scatterers, solvers):
    """
    TODO
    """
    derivative_names = []
    unknown_scatterer_names2  = []
    for name,values in unknown_scatterers.items():
        for variable_name in values['variable_name_list']:
            derivative_names.append(variable_name)
            unknown_scatterer_names2.append(name)
    unknown_scatterer_indices = list(solvers.values())[0]._unknown_scatterer_indices - 1
    unknown_scatterer_names = np.array(list(list(solvers.values())[0].medium.keys()))[unknown_scatterer_indices]
    assert np.all(unknown_scatterer_names == np.atleast_1d(unknown_scatterer_names2)), 'Two different ways of listing unknown scatterer names do not match.'
    derivative_index = pd.MultiIndex.from_arrays([unknown_scatterer_names, derivative_names], names=("scatterer_name","variable_name"))

    grid = list(solvers.values())[0]._grid
    gradient_dataset = xr.Dataset(
                        data_vars = {
                            'gradient': (['x','y','z','derivative_index'],
                                         gradient.reshape((grid.sizes['x'], grid.sizes['y'], grid.sizes['z'],-1)))
                        },
        coords = {
            'x': grid.x,
            'y': grid.y,
            'z': grid.z,
            'derivative_index': derivative_index
        }
    )
    return gradient_dataset

def make_jacobian_dataset(jacobian_list, unknown_scatterers, indices_for_jacobian, solvers,rte_sensors):
    """
    TODO
    """
    merged_jacobian = np.concatenate(jacobian_list, axis=-1)
    split_indices = []
    for rte_sensor in rte_sensors.values():
        size = rte_sensor.sizes['npixels']
        split_indices.append(size)
    split_indices = np.cumsum(split_indices) #TODO verify that no +1 needs to be added etc.
    split_jacobian = np.split(merged_jacobian, split_indices, axis=-1)[:-1] #TODO verify the [:-1]

    grid = list(solvers.values())[0]._grid
    grid_index = pd.MultiIndex.from_arrays([grid.x[indices_for_jacobian[0]], grid.y[indices_for_jacobian[1]],grid.z[indices_for_jacobian[2]]], names=("x", "y","z"))

    derivative_names = []
    unknown_scatterer_names2  = []
    for name,values in unknown_scatterers.items():
        for variable_name in values['variable_name_list']:
            derivative_names.append(variable_name)
            unknown_scatterer_names2.append(name)
    unknown_scatterer_indices = list(solvers.values())[0]._unknown_scatterer_indices - 1
    unknown_scatterer_names = np.array(list(list(solvers.values())[0].medium.keys()))[unknown_scatterer_indices]
    assert np.all(unknown_scatterer_names == np.atleast_1d(unknown_scatterer_names2)), 'Two different ways of listing unknown scatterer names do not match.'
    derivative_index = pd.MultiIndex.from_arrays([unknown_scatterer_names, derivative_names], names=("scatterer_name","variable_name"))

    jacobian_dataset = xr.Dataset(
                                data_vars ={
                                'jacobian_{:1.3f}'.format(wavelength): (['nstokes','derivative_index','grid_index', 'npixels_{:1.3f}'.format(wavelength)], jacobian)
                                        for wavelength, jacobian in zip(solvers.keys(), split_jacobian)
                                },
                    coords={
                    'grid_index': grid_index,
                    'derivative_index': derivative_index
                    }
    )
    return jacobian_dataset

# def gradient_one_solver(rte_solver, sensor, n_jobs, sensor_mapping, forward_sensors, gradient_fun, exact_single_scatter=True,
#                         **kwargs):
#     """
#     TODO
#     """
#     if n_jobs == 1:
#         out = [gradient_fun(rte_solver, sensor)]
#     else:
#         #split by rays
#         split = np.array_split(np.arange(sensor.sizes['nrays']),n_jobs)
#         start_end = [(i.min(),i.max()) for i in split]
#         #adjust start and end indices so that rays are grouped by pixel.
#         index_diffs = sensor.pixel_index.diff(dim='nrays')
#         transitions = np.where(index_diffs.data==1)[0] + 1
#         updated_start_end = []
#         updated_start_end.append((0, transitions[np.abs(transitions - start_end[0][1]).argmin()]))
#         for i in range(1,len(start_end)):
#             new_start = updated_start_end[i-1][1]
#             if i < len(start_end) - 1:
#                 new_end = transitions[np.abs(transitions - start_end[i][1]).argmin()]
#             else:
#                 new_end = start_end[i][1] + 1
#             updated_start_end.append((new_start, new_end))
#
#         pixel_start_end = []
#         pixel_inds = np.cumsum(np.concatenate([np.array([0]), sensor.rays_per_pixel.data])).astype(np.int)
#         for start,end in updated_start_end:
#             a = np.where(pixel_inds == start)[0][0]
#             b=np.where(pixel_inds == end)[0][0]
#             pixel_start_end.append((a,b))
#         out = Parallel(n_jobs=n_jobs, backend="threading")(delayed(gradient_fun, check_pickle=False)(rte_solver,sensor.sel(nrays=slice(start,end),
#         npixels=slice(pix_start, pix_end)),exact_single_scatter=exact_single_scatter, **kwargs)
#          for (start,end),(pix_start,pix_end) in zip(updated_start_end,pixel_start_end))
#
#     #DOESN"T suppport CORRELATED ERRORS.
#     loss = np.sum(np.array([i[1] for i in out]))
#     gradient = np.sum(np.stack([i[0] for i in out],axis=-1),axis=-1)
#
#     #Modify the forward_sensors in-place to include the latest observations.
#     var_list_nray = [str(name) for name in out[0][2].data_vars if str(name) not in ('rays_per_image', 'stokes',
#                                                                                    'rays_per_pixel', 'uncertainties',
#                                                                                    'measurement_data','stokes_weights',
#                                                                                    'I','Q','U','V')]
#     var_list_npixel = [str(name) for name in out[0][2].data_vars if str(name) in ('rays_per_pixel', 'uncertainties',
#                                                                                    'measurement_data','stokes_weights',
#                                                                                    'I','Q','U','V')]
#     merged = {}
#     for var in var_list_nray:
#         concatenated= xr.concat([data[2][var] for data in out], dim='nrays')
#         merged[var] = concatenated
#     merged2 = {}
#     for var in var_list_npixel:
#         concatenated= xr.concat([data[2][var] for data in out], dim='npixels')
#         merged2[var] = concatenated
#     merged.update(merged2)
#     merged['stokes'] = out[0][2].stokes
#     merged['rays_per_image'] = out[0][2].rays_per_image
#     integrated_rays = xr.Dataset(merged)
#
#     rendered_rays = shdom.sensor.split_sensor_pixels(integrated_rays)
#     for i,rendered_ray in enumerate(rendered_rays):
#         mapping = sensor_mapping[i]
#         forward_sensor = forward_sensors[mapping[0]]['sensor_list'][mapping[1]]
#         for stokes in forward_sensor.stokes_index:
#             if rendered_ray['stokes'].sel({'stokes_index':stokes}):
#                 #modification in-place.
#                 forward_sensor[str(stokes.data)] = rendered_ray[str(stokes.data)]
#
#     return gradient, loss
