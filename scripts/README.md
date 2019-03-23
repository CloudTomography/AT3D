# Scripts

Scripts are divided into two categories:

2. Render
Rendering scripts solve the forward problem of simulating images from atmospheric optical properties.
A rendering script takes a generator (see: shdom/generate.py for more info) as input and outputs radiance images of the medium using SHDOM.
The resulting measurements and ground-truth are saved for subsequent use by an optimize script.


3. Optimize
Optimization scripts solve the inverse problem of recovering medium properties from radiance measurements. 
The inputs are measurements and the outputs are estimated medium parameters.