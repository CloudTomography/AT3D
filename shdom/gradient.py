def calculate_direct_beam_derivative(solvers):
    """
    Calculate the geometry of the direct beam at each point and solver.
    Solver is modified in-place.
    TODO
    """
    #Intheory due to setting ._npts to ._nbpts this function can now be called
    #after solver solution. TODO
    #If it can't then ._init_solution needs to be called separately.

    #Only use the first solver.
    rte_solver = solvers.values()[0]

    #calculate the solar direct beam on the base grid
    #which ensures the solver has the required information to
    #calculate the derivative.
    #TODO check that only one solver needs to have the direct_beam calcualted.
    rte_solver._make_direct()

    direct_derivative_path, direct_derivative_ptr = \
        core.make_direct_derivative(
            npts=rte_solver._nbpts, #changed from ._npts to ._nbpts as this should only be calculated for base grid.
            bcflag=rte_solver._bcflag,
            gridpos=rte_solver._gridpos,
            npx=rte_solver._pa.npx,
            npy=rte_solver._pa.npy,
            npz=rte_solver._pa.npz,
            delx=rte_solver._pa.delx,
            dely=rte_solver._pa.dely,
            xstart=rte_solver._pa.xstart,
            ystart=rte_solver._pa.ystart,
            zlevels=rte_solver._pa.zlevels,
            ipdirect=rte_solver._ipdirect,
            di=rte_solver._di,
            dj=rte_solver._dj,
            dk=rte_solver._dk,
            epss=rte_solver._epss,
            epsz=rte_solver._epsz,
            xdomain=rte_solver._xdomain,
            ydomain=rte_solver._ydomain,
            cx=rte_solver._cx,
            cy=rte_solver._cy,
            cz=rte_solver._cz,
            cxinv=rte_solver._cxinv,
            cyinv=rte_solver._cyinv,
            czinv=rte_solver._czinv,
            uniformzlev=rte_solver._uniformzlev,
            delxd=rte_solver._delxd,
            delyd=rte_solver._delyd
        )


    for solver in solvers.values():
        solver._direct_derivative_path = direct_derivative_path
        solver._direct_derivative_ptr = direct_derivative_ptr


def create_derivative_tables(solvers,unknown_scatterers):
    """
    TODO
    """
    partial_derivative_tables = OrderedDict()

    for key in solvers.keys():
        medium = solvers[key].medium
        partial_derivative_list = []
        for unknown_scatterer, table, variable_names in unknown_scatterers[key]:
            variable_names = np.atleast_1d(variable_names)
            #For each unknown_scatterer compute the tables of derivatives
            valid_microphysics_coords = {name:table[name] for name in table.coords if name not in ('table_index', 'stokes_index')}
            derivatives_tables = OrderedDict()
            for variable_name in variable_names:

                if variable_name is in valid_microphysics_coords.keys():
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

                derivatives_tables[variable_name] = differentiated

            #find the index which each unknown scatterer corresponds to in the scatterer list.
            for i,scatterer in enumerate(medium):
                if unknown_scatterer is scatterer:
                    partial_derivative_list.append((i, derivatives_tables))
        unknown_scatterer_indices[key] = partial_derivative_list

    return partial_derivative_tables


def get_derivatives(solvers, derivative_tables):
    """
    TODO
    Calculates partial derivatives on the rte_grid.
    The solvers are modified in place to now contain the partial_derivatives
    """
    for key, rte_solver in solvers.items():

        solver._precompute_phase()
        derivative_table = derivative_tables[key]
        num_derivatives = len(derivative_table)

        dext = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.float32)
        dalb = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.float32)
        diphase = np.zeros(shape=[rte_solver._nbpts, num_derivatives], dtype=np.int32)

        for count, (i,table) in enumerate(derivative_table):
            scatterer = rte_solver.medium[i]
            derivative_on_grid = shdom.medium.table_to_grid(scatterer,table)

            dext[:, count] = derivative_on_grid.extinction.data.ravel()
            dalb[:, count] = derivative_on_grid.albedo.data.ravel()
            diphase[:, count] = derivative_on_grid.table_index.data.ravel() + diphase.max()

        # Concatenate all scatterer tables into one table
        max_legendre = max([table.sizes['legendre_index'] for i,table in derivative_table])
        padded_legcoefs = [table.legcoef.pad({'legendre_index': (0, max_legendre - table.legcoef.sizes['legendre_index'])}) for i,table in derivative_table]
        legendre_table = xr.concat(padded_legcoefs, dim='table_index')

        dnumphase = legendre_table.sizes['table_index']
        dleg = legendre_table.legcoef.data

        # zero the first term of the first component of the phase function
        # gradient. Pre-scale the legendre moments by 1/(2*l+1) which
        # is done in the forward problem in TRILIN_INTERP_PROP
        scaling_factor =np.array([2.0*i+1.0 for i in range(0,rte_solver._nleg+1)])
        dleg[0,0,:] = 0.0
        dleg = dleg[:rte_solver._nstleg] / scaling_factor

        dphasetab = core.precompute_phase_check_grad(
                                                     negcheck=False,
                                                     nstphase=rte_solver._nstphase,
                                                     nstleg=rte_solver._nstleg,
                                                     nscatangle=rte_solver._nscatangle,
                                                     nstokes=rte_solver._nstokes,
                                                     dnumphase=dnumphase,
                                                     ml=rte_solver._ml,
                                                     nlm=rte_solver._nlm,
                                                     nleg=rte_solver._nleg,
                                                     dleg=dleg,
                                                     deltam=rte_solver._deltam
                                                     )
        rte_solver._dphasetab, rte_solver._dext, rte_solver._dalb, rte_solver._diphase, rte_solver._dleg = \
        dphasetab, dext, dalb, diphase, dleg

def levis_approx_uncorrelated_l2(measurements, solvers, forward_sensors, unknown_scatterers,
                                                       table_derivatives, n_jobs=1,
                                 mpi_comm=None,verbose=False):
    """TODO"""
    #note this division of solving and 'raytracing' is not optimal for distributed memory parallelization
    #where tasks take varying amounts of time.
    shdom.script_util.parallel_solve(solvers, n_jobs=n_jobs, mpi_comm=mpi_comm,verbose=verbose)

    #These are called after the solution because they require at least _init_solution to be run.
    if not hasattr(solver.values()[0], '_direct_derivative_path'):
        #only needs to be called once.
        #If this doesn't work after 'solve' and requires only '_init_solution' then it will need to be
        #called inside parallel_solve. . .
        calculate_direct_beam_derivative(solvers)

    get_derivatives(solvers, table_derivatives)  #adds the _dext/_dleg/_dalb etc to the solvers.

    #prepare the sensors for the fortran subroutine for calculating gradient.
    merged_sensors = prepare_sensor_inverse(forward_sensors, solvers)

    #parallelized calculation of the gradient
    #TODO

    #merge gradients
    if mpi_comm is not None:
        #mpi_sum the gradients.

    return loss, gradient


def grad_l2(self, rte_solver, sensor, uncertainties,
            jacobian_flag=False):
    """
    TODO INPROGRESS/NONFUNCTIONAL
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
    assert camx.ndim == camy.ndim==camz.ndim==cammu.ndim==camphi.ndim==1
    total_pix = sensor.sizes['nrays']

    #TODO construct

    if jacobian_flag:
        #maximum size is hard coded - should be roughly an order of magnitude
        #larger than the maximum necssary size but is possible source of seg faults.
        largest_dim = int(100*np.sqrt(rte_solver._nx**2+rte_solver._ny**2+rte_solver._nz**2))
        jacobian = np.empty((rte_solver._nstokes,self.num_derivatives,total_pix*largest_dim),order='F',dtype=np.float32)
        jacobian_ptr = np.empty((2,total_pix*largest_dim),order='F',dtype=np.int32)
    else:
        jacobian = np.empty((rte_solver._nstokes,self.num_derivatives,1),order='F',dtype=np.float32)
        jacobian_ptr = np.empty((2,1),order='F',dtype=np.int32)

    gradient, loss, images, jacobian, jacobian_ptr, counter = core.gradient_l2(
        uncertainties=uncertainties,
        weights=self._stokes_weights[:rte_solver._nstokes],
        exact_single_scatter=self._exact_single_scatter,
        nstphase=rte_solver._nstphase,
        dpath=self._direct_derivative_path,
        dptr=self._direct_derivative_ptr,
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
        partder=self.unknown_scatterers_indices,
        numder=self.num_derivatives,
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
        bcrad=rte_solver._bcrad,
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
        measurements=pixels,
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
        return gradient, loss, images
