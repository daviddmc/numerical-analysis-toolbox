# numerical-analysis-toolbox
basic numerical analysis algorithms implemented in MATLAB for coursework.

---

## algorithms list

### 1.Linear System

#### matrix functions

- [x] norm (1,2,infinity,Fibonacci)
- [x] condition number
- [ ] rank

#### direct methods for solving linear system
 
- [x] Gaussian elimination
    - [x] without pivoting
    - [x] with partial pivoting
    - [x] with complete pivoting
- [x] Thomas method

#### matrix decomposition

- [x] LU (Doolittle) decomposition 
    - [x] with partial pivoting
    - [x] with complete pivoting
- [x] decomposition for symmetric positive definite matrix
    - [x] Cholesky decomposition
    - [x] LDL decomposition
- [x] QR decomposition
    - [x] Gram-Schimidt orthogonalization
    - [x] Givens transformation
    - [x] Householder transformation
    - [x] column pivoting

#### inverse matrix
    - [x] Gauss-Jordan method with complete pivoting
    - [x] triangle
    - [x] QR
    - [ ] SVD
    - [x] LU

#### iterative method for solving linear system

- [x] Jacobi iterative method
- [x] Gauss-Seidel iterative method
- [x] succesise over - relaxation method
    - [x] SOR
    - [x] symmetric SOR
- [ ] conjugate gradient method
    - [ ] steepest descent
    - [x] conjugate gradient
    - [ ] preconditioned conjugate gradient 

### 2.Nonlinear Equation and System

#### single variable

- [x] bisection method
- [x] Steffensen acceleration for fixed-point iteration
- [x] Newton method
- [x] secant method
- [x] Brent method

#### multivariable

- [ ] Newton method
- [x] Broyden method

#### polynomail



### 3.Eigenvalue and Singular Value

#### reduction

- [x] Hessenberg
- [x] bidiagonal
- [x] tridiagonal

#### iterative method

- [x] power iteration
- [x] inverse iteration (with shift)
- [x] Rayleigh quotient iteration
- [x] deflation method

#### QR iteration

- [x] single shift
- [x] double shift
- [x] symmetric

#### others

- [x] Jacobi method
- [x] balances
- [x] singular value decomposition

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
- [x] least square problem

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

## Demo

- Runge effect
- Minimize error of polynomial interpolation using Chebyshev interpolation
- optimal relaxation factor for SOR

---

## Bibliography

- 数值分析基础(第二版),关治,陆金甫. Fundamentals of Numerical Analysis (Second Edition), Zhi Guan, Jinfu Lu.
- Numerical Analysis, Rainer Kress
- An Introduction to Numerical Analysis, Endre suli and David F. Mayers
- Applied Numerical Linear Algebra, James W. Demmel
- Numerical Analysis (2nd edition), Timothy Sauer
- Numerical Analysis (9th edition), Richard L. Burden, J. Douglas Faires
- Matrix Computation (4th edition), Gene H. Golub, Charles F. Van Loan
- Numerical Optimization (2nd edition), Jorge Nocedal, Stephen J. Wright