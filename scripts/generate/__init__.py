"""
Generation scripts generate the atmospheric optical properties and grid. 
They can be subsequently input into a rendering script for forward modeling using SHDOM. 

A generation script should have all of the following components:


1. update_parser(parser) method
-------------------------------
This is a method to update the argument parser. It recieves the main parser as input 
and adds the necessary arguments for the generate script.

Parameters:

parser: argparse.ArgumentParser()
    The main parser to update.

Returns:

parser: argparse.ArgumentParser()
    The updated parser.



2. grid(args) method 
--------------------
This method generates the atmospheric grid.

Parameters:

args: args
    The arguments used by the method. They are defined in the update_parser method.

Returns:

grid: shdom.Grid
    The medium grid.  
  
  
  
3. extinction(grid, args) method 
---------------------------
This method generates the optical extinction field.

Parameters:

grid: shdom.Grid
    The arguments used by the method which are defined in the update_parser method.
args: args
    The arguments used by the method. They are defined in the update_parser method.


Returns:

extinction: shdom.GridData object
    This object contains the extinction on the pre-defined grid.
  
  
  
4. albedo(grid, args) method 
---------------------------
This method generates the single scattering albedo field.

Parameters:

grid: shdom.Grid
    The arguments used by the method which are defined in the update_parser method.
args: args
    The arguments used by the method. They are defined in the update_parser method.


Returns:

albedo: shdom.GridData object
    This object contains the albdeo on the pre-defined grid. 
  
  
  
5. phase(grid, args) method 
---------------------------
This method generates the phase function.

Parameters:

grid: shdom.Grid
    The arguments used by the method which are defined in the update_parser method.
args: args
    The arguments used by the method. They are defined in the update_parser method.


Returns:

phase: shdom.Phase object
    This could either be a shdom.TabulatedPhase or shdom.GridPhase object. 
  
"""


