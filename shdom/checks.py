import shdom

def check_positivity(dataset,*names):
    """
    TODO
    """
    fail_list = []
    for name in names:
        variable = dataset[name]
        if not np.all(variable.data >= 0.0):
            fail_list.append(name)
    if len(fail_list) > 0:
        raise ValueError("variables are not positive semi-definite ('>= 0.0')", *fail_list)

def check_grid(dataset):
    """
    TODO
    """
    #check equispacing of horizontal grid 'x','y'
    x = dataset.x.diff('x')
    y = dataset.y.diff('y')

    fail_list = []
    if not np.all(x==x[0]):
        fail_list.append('x')
    if not np.all(y==y[0]):
        fail_list.append('y')
    if len(fail_list) > 0:
        raise ValueError("Grid coordinates must be equispaced for the RTE solver.", *fail_list)
    check_positivity(dataset,'z')

def check_legendre(dataset):
    """TODO"""
    legendre_table = dataset['legcoefs']
    if not np.all(legendre_table[0] == 1.0):
        raise ValueError("0th Legendre Coefficients must be normalized to 1.0")

def check_range(dataset, *argchecks):
    """TODO"""
    fail_list = []
    for (name,low,high) in argchecks:
        data = dataset[name]
        if not np.all(low <=data<= high):
            fail_list.append((name, low, high))

    if len(fail_list) > 0:
        raise ValueError('Variables not in required ranges', *fail_list)

def dataset_checks(**argchecks):
    """
    TODO
    A complicated decorator that allows
    tests to be applied to datasets or lists of datasets.
    See below for example usage.
    """

    def decorator(func):
        code = func.__code__
        allargs  = code.co_varnames[:code.co_argcount]

        def tests(*pargs,**kargs):
            positionals = list(allargs)
            for argname in kargs:
                # for all passed by name, remove from expected
                positionals.remove(argname)

            #loop through all dataset variable names assigned to test.
            for argname,test_list in argchecks.items():
                #tidy up test_list
                test_list = test_list if isinstance(test_list,list) else [test_list]

                if argname in kargs:
                    #arg was passed as a keyword argument

                    #If the arg is not a list of (presumed) datasets
                    #then make it a list
                    if isinstance(kargs[argname], list):
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
                                #print an error message for all failed tests
                                if list_flag:
                                    print("The '{}' argument of '{}' failed with {} from '{}'".format(i,argname,err.args,test.__name__))
                                else:
                                    print("'{}' failed with {} from '{}'".format(argname,err.args,test.__name__))
                else:
                    #arg was passed as a positional argument
                    position = positionals.index(argname)
                    if isinstance(pargs[position], list):
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
                                #print an error message for all failed tests
                                if list_flag:
                                    print("The '{}' argument of '{}' failed with {} from '{}'".format(i,argname,err.args,test.__name__))
                                else:
                                    print("'{}' failed with {} from '{}'".format(argname,err.args,test.__name__))

            raise ValueError("tests specified by dataset_checks failed as input to '{}'".format(func.__name__))
            return func(*pargs, **kargs)
        return tests
    return decorator

#------ example usage ---------------------
#We specifiy which variables are datasets/lists_of_datasets that we want to apply checks to ('microphysics, poly_table').
#We then specify which tests we want to apply and what variable names we want to apply them to (if applicable)
#A LIST of checks is supplied for each dataset/list_of_dataset name,
#For each check a TUPLE of the function, and a LIST of variable names is supplied.

#checks to apply ('check_positivity') to which datasets/lists_of_datasets ('microphysics')
@dataset_checks(microphysics=(check_positivity,'density'),
               poly_table=[(check_positivity, 'big_banananas')],
               legendre=check_legendre)
def fake_function_that_does_nothing(microphysics, legendre, poly_table):
    return

# import xarray as xr
# import numpy as np
# microphysics = xr.Dataset(data_vars={'density': ('a',np.ones(50)*-1)})
# poly_table = xr.Dataset(data_vars={'big_banananas': ('a',np.ones(50)*-1)})
#legendre = xr.Dataset(data_vars={'legcoefs': (['a','b'], np.ones((10,10))+0.01)})
#fake_function_that_does_nothing([microphysics,microphysics],legendre=legendre,poly_table=[poly_table,poly_table])
