import numpy as np
import pyshdom.exceptions
import typing

def check_exists(dataset, *names):
    """TODO"""
    for name in names:
        try:
            dataset[name]
        except KeyError as err:
            raise type(err)(str("Expected variable with name '{}' in dataset".format(name)))

def check_positivity(dataset, *names, precision=7):
    """
    TODO
    """
    check_exists(dataset, *names)
    for name in names:
        variable = dataset[name]
        if not np.all(variable.data >= 0.0):
            dataset[name][:] = np.round(variable, decimals=precision)
            if not np.all(variable.data >= 0.0):
                raise pyshdom.exceptions.NegativeValueError("Negative values found in '{}'".format(name))

def check_range(dataset, **checkkwargs):
    """
    TODO
    """
    check_exists(dataset, *checkkwargs)
    for name, (low, high) in checkkwargs.items():
        data = dataset[name]
        if not np.all((low <= data) & (data <= high)):
            raise pyshdom.exceptions.OutOfRangeError(
                "Values outside of range '[{}, {}]' found in variable '{}'".format(
                low, high, name)
            )

def check_hasdim(dataset, **checkkwargs):
    """
    TODO
    """
    check_exists(dataset, *checkkwargs)
    for name, dim_names in checkkwargs.items():
        if not isinstance(dim_names, typing.Iterable):
            dim_names = [dim_names]
        for dim_name in dim_names:
            if dim_name not in dataset[name].dims:
                raise pyshdom.exceptions.MissingDimensionError(
                    "Expected '{}' to have dimension '{}'".format(
                    name, dim_name
                ))

def check_grid(dataset):
    """
    TODO
    """
    pyshdom.checks.check_exists(dataset, 'x', 'y', 'z', 'delx', 'dely')
    for dimension, deldim in zip(('x', 'y'), ('delx', 'dely')):

        if dataset[dimension][0] != 0.0:
            raise pyshdom.exceptions.GridError("Grid dimension '{}' should start from 0.0".format(dimension))

        if dataset[dimension].size > 1:
            diffx = dataset[dimension].diff(dimension).data
            if not (np.allclose(diffx, diffx[0], atol=1e-6)):
                raise pyshdom.exceptions.GridError("Grid dimension '{}' is not equispaced.".format(dimension))
            if not (np.allclose(diffx[0], dataset[deldim])):
                raise pyshdom.exceptions.GridError("'{a}' is not consistent with '{b}'. "
                                                   "'{a}' should be set based on '{b}', see grid.make_grid for details.".format(
                                                       a=deldim, b=dimension))
            if not np.all(diffx > 0.0):
                raise pyshdom.exceptions.GridError("Grid dimension '{}' is not strictly increasing.".format(dimension))

        if dataset[deldim] <= 0.0:
            raise pyshdom.exceptions.GridError("Grid dimension '{}' is not strictly increasing.".format(dimension))

    if not np.all(dataset.z >= 0.0) & np.all(dataset.z.diff('z') > 0.0) & (dataset.z.size >= 2):
        raise pyshdom.exceptions.GridError("Grid dimension 'z' should be positive, strictly increasing and "
                                           "have 2 or more elements.")

def check_legendre(dataset):
    """TODO"""
    check_exists(dataset, 'legcoef', 'table_index', 'extinction', 'ssalb')
    check_hasdim(dataset, legcoef=['stokes_index', 'legendre_index'])

    if dataset['legcoef'].sizes['stokes_index'] != 6:
        raise pyshdom.exceptions.LegendreTableError("'stokes_index' dimension of 'legcoef' must have 6 components.")
    legendre_table = dataset['legcoef']
    if not np.allclose(legendre_table[0, 0], 1.0, atol=1e-7):
        raise pyshdom.exceptions.LegendreTableError("0th Legendre/Wigner Coefficients must be normalized to 1.0")
    if not np.all(( legendre_table[0, 1]/3.0 >= -1.0) & (legendre_table[0, 1]/3.0 <= 1.0)):
        raise pyshdom.exceptions.LegendreTableError("Asymmetry Parameter (1st Legendre coefficient divided by 3)"
                                                   "is not in the range [-1.0, 1.0]")

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
