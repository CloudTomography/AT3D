import numpy as np

def check_positivity(dataset,*names):
    """
    TODO
    """
    check_exists(dataset, *names)
    fail_list = []
    for name in names:
        variable = dataset[name]
        if not np.all(variable.data >= 0.0):
            fail_list.append(name)
    if len(fail_list) > 0:
        raise ValueError("variables must be >= 0.0", *fail_list)

def check_hasdim(dataset, *argchecks):
    name_list = [a for a,b in argchecks]
    check_exists(dataset, *name_list)
    fail_list = []
    for (name,dim_names) in argchecks:
        for dim_name in dim_names:
            if dim_name not in dataset[name].dims:
                fail_list.append((name, dim_name))
    if len(fail_list) > 0:
        raise KeyError('Variables do not have the expected dimensions', *fail_list)

def check_exists(dataset, *names):
    fail_list = []
    for name in names:
        try:
            variable = dataset[name]
        except KeyError:
            fail_list.append(name)
    if len(fail_list) > 0:
        raise ValueError("variables were expected in dataset", *fail_list)

def check_grid(dataset):
    """
    TODO
    """
    check_exists(dataset, 'x', 'y', 'z')

    #check grid is 3D
    # fail_list = []
    # if len(dataset.x) == 1:
    #     fail_list.append('x')
    # if len(dataset.y) == 1:
    #     fail_list.append('y')
    # if len(fail_list) > 0:
    #     raise ValueError("Grid coordinates with size of 1 are not currently supported.", *fail_list)

    #check equispacing of horizontal grid 'x','y'
    x = dataset.x.diff('x')
    y = dataset.y.diff('y')

    fail_list = []
    if len(x) > 1:
        if not np.allclose(x,x[0],atol=1e-6):
            fail_list.append('x')
    if len(y) > 1:
        if not np.allclose(y,y[0], atol=1e-6):
            fail_list.append('y')
    if len(fail_list) > 0:
        raise ValueError("Grid coordinates must be equispaced for the RTE solver.", *fail_list)
    check_positivity(dataset,'z')

def check_legendre_poly_table(dataset):
    """TODO"""
    check_exists(dataset, ['legcoef', 'table_index', 'extinction', 'ssalb'])
    check_hasdim(dataset, ('legcoef', ['stokes_index', 'legendre_index']))

    if not np.all(dataset.table_index.data.ravel()== np.arange(1,dataset.extinction.size+1).astype(np.int)):
        raise ValueError("'table_index' must act as a multi_index for microphysical dims")

    if dataset['legcoef'].sizes['stokes_index'] != 6:
        raise ValueError("'stokes_index' dimension of 'legcoef' must have 6 components.")

    legendre_table = dataset['legcoef']
    if not np.allclose(legendre_table[0,0], 1.0):
        raise ValueError("0th Legendre Coefficients must be normalized to 1.0")

def check_legendre_rte(dataset):
    """TODO"""
    check_exists(dataset, ['legcoef', 'table_index', 'extinction', 'ssalb'])
    check_hasdim(dataset, ('legcoef', ['stokes_index', 'legendre_index']))

    if dataset['legcoef'].sizes['stokes_index'] != 6:
        raise ValueError("'stokes_index' dimension of 'legcoef' must have 6 components.")

    legendre_table = dataset['legcoef']
    if not np.allclose(legendre_table[0,0], 1.0):
        raise ValueError("0th Legendre Coefficients must be normalized to 1.0")

def check_range(dataset, *argchecks):
    """TODO"""
    name_list = [a for a,b,c in argchecks]
    check_exists(dataset, *name_list)
    fail_list = []
    for (name,low,high) in argchecks:
        data = dataset[name]
        if not np.all((low <=data) & (data<= high)):
            fail_list.append((name, low, high))

    if len(fail_list) > 0:
        raise ValueError('Variables not in required ranges', *fail_list)

def dataset_checks(**argchecks):
    """
    TODO
    A complicated decorator that allows
    tests to be applied to datasets or lists of datasets.
    """

    def decorator(func):
        code = func.__code__
        allargs  = code.co_varnames[:code.co_argcount]

        def tests(*pargs,**kargs):
            positionals = list(allargs)

            for argname in kargs:
                # for all passed by name, remove from expected
                positionals.remove(argname)
            test_fail = False
            #loop through all dataset variable names assigned to test.
            for argname,test_list in argchecks.items():
                #tidy up test_list
                test_list = test_list if isinstance(test_list,list) else [test_list]

                if argname in kargs:
                    #arg was passed as a keyword argument

                    #If the arg is not a list of (presumed) datasets
                    #then make it a list
                    if isinstance(kargs[argname], dict):
                        karg_list = list(kargs[argname].values())
                        list_flag = True
                    elif isinstance(kargs[argname], list):
                        karg_list = kargs[argname]
                        list_flag = True
                    else:
                        karg_list = [kargs[argname]]
                        list_flag = False
                    #iterate over all datasets in the list of datasets.
                    for i,single_arg in enumerate(karg_list):
                        #iterate over all tests.
                        for test in test_list:
                            #tidy up test_list args
                            if isinstance(test,tuple):
                                test,test_args = test
                                test_args = test_args if isinstance(test_args,list) else [test_args]
                                test_args = [single_arg] + test_args
                            else:
                                test = test
                                test_args = [single_arg]

                            #do the test
                            try:
                                test(*test_args)
                            except ValueError as err:
                                test_fail = True
                                #print an error message for all failed tests
                                if list_flag:
                                    print("The '{}' argument of '{}' failed with {} from '{}'".format(i,argname,err.args,test.__name__))
                                else:
                                    print("'{}' failed with {} from '{}'".format(argname,err.args,test.__name__))
                else:
                    #arg was passed as a positional argument
                    position = positionals.index(argname)
                    if isinstance(pargs[position], dict):
                        arg_list = list(pargs[position].values())
                        list_flag = True
                    elif isinstance(pargs[position], list):
                        arg_list = pargs[position]
                        list_flag = True
                    else:
                        arg_list = [pargs[position]]
                        list_flag = False
                    #iterate over all datasets in the list of datasets.
                    for i,single_arg in enumerate(arg_list):
                        for test in test_list:
                            #tidy up test_args
                            if isinstance(test,tuple):
                                test,test_args = test
                                test_args = test_args if isinstance(test_args,list) else [test_args]
                                test_args = [single_arg] + test_args
                            else:
                                test = test
                                test_args = [single_arg]

                            #do the test
                            try:
                                test(*test_args)
                            except ValueError as err:
                                test_fail=True
                                #print an error message for all failed tests
                                if list_flag:
                                    print("The '{}' argument of '{}' failed with {} from '{}'".format(i,argname,err.args,test.__name__))
                                else:
                                    print("'{}' failed with {} from '{}'".format(argname,err.args,test.__name__))
            if test_fail:
                raise ValueError("tests specified by dataset_checks failed as input to '{}'".format(func.__name__))
            return func(*pargs, **kargs)
        return tests
    return decorator

#------ example usage ---------------------
#We specifiy which variables are datasets/lists_of_datasets that we want to apply checks to ('microphysics, poly_table').
#We then specify which tests we want to apply and what variable names we want to apply them to (if applicable)
#A LIST of checks is supplied for each dataset/list_of_dataset name,
#For each check a TUPLE of the function, and a LIST of variable names is supplied.

# @dataset_checks(microphysics=(check_positivity,'density'),
#                poly_table=[(check_positivity, 'big_banananas')],
#                legendre=check_legendre)
# def fake_function_that_does_nothing(microphysics, legendre, poly_table):
#     return

# import xarray as xr
# import numpy as np
# microphysics = xr.Dataset(data_vars={'density': ('a',np.ones(50)*-1)})
# poly_table = xr.Dataset(data_vars={'big_banananas': ('a',np.ones(50)*-1)})
#legendre = xr.Dataset(data_vars={'legcoefs': (['a','b'], np.ones((10,10))+0.01)})
#fake_function_that_does_nothing([microphysics,microphysics],legendre=legendre,poly_table=[poly_table,poly_table])
