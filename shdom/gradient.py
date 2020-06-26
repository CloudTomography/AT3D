def calculate_direct_beam_derivative(solvers):
    """
    Calculate the geometry of the direct beam at each point and solver.
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

    return solar_paths
