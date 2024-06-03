# Quick example of Reaction-diffusion-advection

The 1-dimensional reaction-diffusion-advection system.


## Current version

* The diffusion equation and two advection equations (up and down) that are uncoupled, no reactions.

* Matlab (consider making a Julia version)
* Forward-Euler for diffusion (consider changing to Backward Euler)

$$\frac{\partial u_D}{\partial t} = D \frac{\partial^2 u_D}{\partial z^2}$$

$$\frac{\partial u_v}{\partial t} = -v \frac{\partial u_v}{\partial z}$$

$$\frac{\partial u_v}{\partial t} = +v \frac{\partial u_v}{\partial z}$$