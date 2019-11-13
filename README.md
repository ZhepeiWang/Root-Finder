# Root-Finder
This is an $\color{red}{univariate polynomial solver}$, which finds all $\color{red}{REAL} roots of any given polynomial.

# Feature

1. As for low order polynomials (linear, quadratic, cubic and quartic polynomials), the solver use $\color{blue}{closed form} solution.
In this case, the solver only takes about 4.0E-7 sec.
2. As for high order polynomials (order >= 5), the solver implements 2 different methods to find all roots. The recommended 
one, named $\color{blue}{Real Roots Isolation Method}, is based on Sturm's Theorem as well as other geometric property of polynomials. 
The other one is based on $\color{blue}{Companion Matrix Method}. The example is the comparision between these two methods.

3. The Real Roots Isolation Methods uses Cauchy’s bound as well as $\color{blue}{Kojima’s bound} to bracket all roots. Normally, the latter is tigher 
than the former for about 10^12 magnitude. Technically, Fujiwara’s bound is always better than Kojima's bound, while Kojima's bound 
is more numerically friendly and tight enough.

4. The Real Roots Isolation Methods uses $\color{blue}{Sturm's theorem} to determine the number of roots inside any given interval.

5. The Real Roots Isolation Methods is faster and more stable than the Companion Matrix Method. However, due to the truncation error 
of float point number, the former is only recommended for at most order-32 polynomials, while the latter is only recommended for 
at-most order-20 polynomials.

6. The solver can be furtherly optimized.

# Interface

Only the function below is needed.

Function:
std::set<double> RootFinder::solvePolyInterval(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol = 1e-3, bool isolation = true);

Inputs:
    coeffs: 
        Eigen VectorXd for coefficients of a polynomial. For example, the polynomial a(n)*x^n+a(n-1)*x^(n-1)+...+a(1)*x+a(0) can be expressed by 
        coeffs=[a(n), a(n-1), ..., a(1), a(0)].

    lbound and ubound:
        Open interval (lbound, ubound). If lbound = -INFINITY and ubound = INFINITY, then all real roots can be found by the solver.
    
    tol:
        tolerance for precision, only works when order is higher than 4.
    
    isolation:
        switch for Method used. Default one is Real Roots Isolation Method.

Outputs:
    std::set<double> which stores all searched real roots.


# Compile

    mkdir build
    cd build
    cmake ..
    make

