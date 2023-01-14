# fSRK

Fortran implementation of the 2nd Order Weak Stochastic Runge-Kutta algorithm for simulating stochastic differential equations in Ito form, 
$$dx = f(x,t)dt + g(x,t)dW, $$ 
where $dW$ is a scalar Weiner noise process, $f$ is the drift term and $g$ is the diffusion term.
Equations are edited inside the Fortran code then imported into Python using `f2py`.

Instructions:
-------------

1. In the Fortran code, create functions (subroutines) for the drift and diffusion terms. 

2. Edit the existing subroutines `drift` and `diffus` so that they call upon these subroutines. 

3. Create python .so module by running command: `f2py -c -m modulename filename.f90`

4. The package may then used in python with `from modulename import *`.

Examples
--------
The examples show the application of the code to noisy dynamical systems and open quantum systems. These examples should clarify how parameters such as the initial condition, length of simulation and model parameters are passed from the Python code to Fortran.

<img src="https://github.com/alexeistepa/fSRK/blob/main/zeno_figa.png?raw=true" width="400" height="400">

References
----------
1. Breuer and Petruccione - The Theory of Open Quantum Systems pg 362 eq (7.47) - eq. (7.49).  

2. Kloeden and Platen - Numerical solutions of stochastic differential equations pg 486 eq. (1.1). 

3. Milstein and  Tretyakov - Stochastic Numerics for Mathematical Physics pg 104 eq. 2.20 .
