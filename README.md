# numerical-analysis-toolbox
basic numerical analysis algorithms implemented in MATLAB for coursework.

---

## algorithms list

### Linear System

- [x] matrix functions
    - [x] norm (1,2,infinity,Fibonacci)
    - [x] condition number
- [x] Gaussian elimination
    - [x] Gaussian elimination with partial pivoting
    - [x] LU (Doolittle) decomposition with partial pivoting
    - [x] Thomas method
    - [x] inverse (Gauss-Jordan method with pivoting)
- [x] decomposition for symmetric positive definite matrix
    - [x] Cholesky decomposition
    - [x] LDL decomposition
- [x] QR decomposition
    - [x] Gram-Schimidt orthogonalization
    - [x] Givens transformation
    - [x] Householder transformation
- [x] Jacobi iterative method
- [x] Gauss-Seidel iterative method
- [x] succesise over - relaxation method (SOR)
- [ ] conjugate gradient method
    - [ ] steepest descent
    - [x] conjugate gradient
    - [ ] preconditioned conjugate gradient 
- [x] least square problem

### Nonlinear Equation and System

- [x] bisection method
- [x] Steffensen acceleration for fixed-point iteration
- [x] Newton method
    - [x] Newton iteration
    - [x] Newton iteration for multiplicity > 1
- [x] secant method
- [x] Brent method
- [ ] fixed-point iteration for nonlinear system
- [ ] Newton method for nonlinear system
- [ ] quasi-Newton method

### Eigenvalue and Singular Value

- [x] QR method
- [x] Rayleigh quotient iteration
- [x] Jacobi method
- [x] power method
- [x] inverse power method
- [ ] singular value decomposition

### Interpolation

- [x] polynomial interpolation
    - [x] Lagrange polynomial
    - [x] Newton polynomial
    - [x] Neville algorithm
    - [x] Chebyshev interpolation
- [x] Hermite interpolation
    - [x] Lagrange polynomial
    - [x] Lagrange polynomial
- [x] nearest interpolation
- [x] linear interpolation
- [x] cubic spline interpolation
    -[x] natural end conditions
    -[x] not-a-knot end conditions

### Function Approximation

- [x] orthogonal polynomials
    - [x] Legendre
    - [x] Laguerre
    - [x] Hermite
    - [x] Chebyshev
- [x] polynomial fit
- [x] fast Fourier transform
    - [x] base 2
    - [x] base 3

### Numerical Integration and Differentiation

- [x] trapezoid
- [x] Simpson 
- [x] Gaussian quadrature
    - [x] Gauss-Legendre
    - [x] Gauss-Chebyshev
    - [x] Gauss-Hermite
    - [x] Gauss-Laguerre
- [x] Romberg

### Ordinary Differential Equation

- [x] Euler method
- [x] Runge - Kutta method
    - [x] 2 order RK (midepoint method)
    - [x] 4 order RK
    - [x] Runge - Kutta - Fehlberg
- [x] Adams-Bashforth four-step method with predictor-corrector

### Partial Differential Equation

- [ ] Poisson equation
- [ ] heat equation
- [ ] wave equation

### Stochastic Method

- [x] linear congruential generator

### Optimization

- [x] Nelder - Mead

---

## demo

- Runge effect
- Minimize error of polynomial interpolation using Chebyshev interpolation
- optimal relaxation factor for SOR

---

## reference

- 数值分析基础(第二版),关治,陆金甫. Fundamentals of Numerical Analysis (Second Edition), Zhi Guan, Jinfu Lu.
- Numerical Analysis, Rainer Kress
- An Introduction to Numerical Analysis, Endre suli and David F. Mayers
- Applied Numerical Linear Algebra, James W. Demmel
- Numerical Analysis (2nd edition), Timothy Sauer
- Numerical Analysis (9th edition), Richard L. Burden, J. Douglas Faires
- Matrix Computation (4th edition), Gene H. Golub, Charles F. Van Loan