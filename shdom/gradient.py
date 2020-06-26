def calculate_direct_beam_derivative(solvers):
    """
    Calculate the geometry of the direct beam at each point and solver.
    TODO
    """
    direct_beam_derivative = OrderedDict()

    #All solvers are unsolved so they should have the same base grid.
    #npts = number of base grid points (nbpts) as no solutions have been performed.

    #Only use the first solver.
    rte_solver = solvers.values()[0]

    #calculate the solar direct beam on the base grid
    #which ensures the solver has the required information to
    #calculate the derivative.

    #TODO check that only one solver needs to have the direct_beam calcualted.
    rte_solver._make_direct()

    direct_derivative_path, direct_derivative_ptr = \
        core.make_direct_derivative(
            npts=rte_solver._npts,
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

    direct_beam_derivative['direct_derivative_path'] = direct_derivative_path
    direct_beam_derivative['direct_derivative_ptr'] = direct_derivative_ptr

    return direct_beam_derivative

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
