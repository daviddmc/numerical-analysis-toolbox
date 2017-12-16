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

#### inverse and pseudo-inverse of matrix

- [x] Newton iteration
- [x] Gauss-Jordan method with complete pivoting
- [x] triangle inverse
- [x] QR decomposition
- [x] SVD (Pseudo Inverse)
- [x] LU decomposition

#### iterative method for solving linear system

- [x] Jacobi iterative method
- [x] Gauss-Seidel iterative method
- [x] succesise over - relaxation method
    - [x] SOR
    - [x] symmetric SOR
- [ ] Krylov methods
    - [ ] CG
    - [ ] CGNE
    - [ ] CGNR
    - [ ] GMRES
    - [ ] MINRES
    - [ ] SYMMLQ
    - [ ] QMR
    - [ ] BiCG
    - [ ] CGS
    - [ ] BiCGstab
    - [ ] LSQR

### 2.Nonlinear Equation and System

#### single variable

- [x] bisection method
- [x] Steffensen acceleration for fixed-point iteration
- [x] Newton method
- [x] secant method
- [x] Brent method
- [x] Muller method

#### multivariable

- [ ] Newton method
- [ ] Broyden method
- [ ] solving nonlinear LS problem

#### polynomial

- [x] Horner method
- [x] finding roots using deflation
- [x] finding roots using QR iteration

### 3.Eigenvalue and Singular Value

#### reduction

- [x] Hessenberg
- [x] bidiagonal
- [x] tridiagonal

#### iterative method

- [x] power iteration
- [x] inverse iteration (with shift)
- [x] Rayleigh quotient iteration
- [x] Jacobi method

#### QR iteration

- [x] single shift
- [x] double shift
- [x] symmetric

#### pratical methods

- [x] eigenvalue
- [x] eigenpair of symmetric matrix
- [x] deflation method
- [x] matrix balance
- [x] singular value decomposition

### Interpolation

#### 1D

- [x] polynomial interpolation
    - [x] Lagrange polynomial
    - [x] Newton polynomial
    - [x] Neville algorithm
- [x] Hermite interpolation
    - [x] Lagrange polynomial
    - [x] Newton polynomial
- [x] Chebyshev interpolation
- [x] nearest interpolation
- [x] linear interpolation
- [x] cubic spline interpolation
    -[x] natural end conditions
    -[x] not-a-knot end conditions

#### 2D

- [ ] nearest interpolation
- [ ] linear interpolation

### Function Approximation

#### orthogonal polynomials

- [x] Legendre
- [x] Laguerre
- [x] Hermite
- [x] Chebyshev

#### least square problem

- [x] full rank least square
    - [x] normal equation
    - [x] QR
    - [x] SVD
- [x] rank deficient least square
    - [x] QR
    - [x] SVD
- [x] weighted least square
- [x] Tikhonov regularization
- [x] polynomial fit

#### fast Fourier transform

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

#### initail value problems

- [ ] Runge - Kutta method
    - [ ] classic
    - [ ] embedded
    - [ ] implicit
- [ ] Linear multistep method
    - [ ] Adams methods
    - [ ] implicit Adams methods
    - [ ]
    - [ ] Adams-Bashforth four-step method with predictor-corrector

#### boundary value problems

### Partial Differential Equation

- [ ] Poisson equation
- [ ] heat equation
- [ ] wave equation

### Stochastic Method

- [x] linear congruential generator

### Optimization

#### nonlinear least square

- [ ] Gauss-Newton
- [ ] Levenberg-Marquardt

#### derivative-free method
- [ ] Nelder - Mead

---

## Demo

- Runge phenomenon
- minimize error of polynomial interpolation using Chebyshev interpolation
- optimal relaxation factor for SOR
- round-off error in numerical differentiation
- common orthogonal polynomials
- Chebyshev acceleration in linear iteration
- Gershgorin disc theorem

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