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
