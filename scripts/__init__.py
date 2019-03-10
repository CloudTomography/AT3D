"""
Scripts are divided into three categories:

1. Generate
Generate the medium grid and optical properties: extinction, single-scattering albedo and phase function.

2. Render
Render scripts solve the forward problem. The input is a generated medium (a generate script).
The outputs are measurements, i.e. synthetic images of the medium using SHDOM.

3. Optimize
Optimize scripts solve the inverse problem. The inputs are measurements and the outputs are the medium parameters.
"""

import generate, render, optimize