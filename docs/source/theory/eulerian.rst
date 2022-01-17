Eulerian simulation of drift trajectories
===========================================

Convection
-----------

Convection consists of advection and diffusion.

Diffusion is given by:

.. math::

  \frac{\partial U}{\partial t} = D \left( \frac{\partial^2 U}{\partial x^2} + \frac{\partial^2 U}{\partial y^2} \right)


The convection equation is (`wiki
<https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation>`_):

.. math::

  \frac{\partial c}{\partial t} = ...

with the assumptions that:

  * the diffusion constant `D` is constant for the field,
  * and that the flow `u` is incompressible (i.e. has *no divergence*).

the equation simplifies to:

.. math::

  \frac{\partial c}{\partial t} = D \nabla^2 c - \mathbf{v} \cdot \nabla T

where :math:`\nabla^2 = \triangle` is the Laplacian.


Diffusion
---------
Diffusivity (:math:`m^2/s`). E.g. between 0.01 and 0.1 for oil on the
surface of the ocean (`Matsuzakia et. al., 2017
<https://www.sciencedirect.com/science/article/pii/S0025326X16308426>`_).

Decreasing diffusivity places stricter stability criteria on time step.

Porosity
________

Porosity, rate of liquid volume to total volume (fraction of flux)

Numerical schemes
-----------------

Explicit simulation
___________________

A simple explicit scheme for integrating the convection-equation.

* Forward difference in time
* `ndimage.laplace` and `np.gradient` for spatial differences.

https://en.wikipedia.org/wiki/Numerical_solution_of_the_convection%E2%80%93diffusion_equation#Solving_the_convection%E2%80%93diffusion_equation_using_the_finite_difference_method

