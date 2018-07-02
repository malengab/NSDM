# NSDM
Numerical steepest descent method

NSDM is used for a fast evaluation of oscillatory integrals of the type $\int f(x)e^{ig(x)/\varepsilon} dx$ where $\varepsilon\ll 1$.
We use this method in particular for $g$ quadratic complex polynomial.

Files: 

Gaussian_like_compare_eps.m

Setting up A, B, C in $g(x) = Ax^2 + Bx + C$, we can compute the integral on [-10,10] (can be changed manually) by employing the steepest descent solvers (Gauss-Laguerre quadrature in steepest_descent.m or Hermite quadrature in steepest_descent_H.m). Function f is chosen as $f(x) = e^{-dx^2}$.
The numerical solution is compared to the exact one and plotted as a function of $\epsilon$.
Needs GaussLaguerre.m to produce Laguerre points and weights.
