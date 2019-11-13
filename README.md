# Root-Finder
This is an univariate polynomial solver, which finds all REAL roots or just all real roots in a given interval.

# Feature

1. As for low order polynomials (linear, quadratic, cubic and quartic polynomials), the solver use closed form solution.
In this case, the solver only takes about 0.4E-6 sec.
2. As for high order polynomials (order >= 5), the solver implements 2 different methods to find all roots. The recommended 
one, names Real Roots Isolation Method, is based on Sturm's Theorem as well as other geometric property of polynomials. 
The other one is based on Companion Matrix Method. The example is the comparision between these two methods. Some results is shown below:

    --------------------- Speed Test between Roots Isolation Method and Companion Matrix Method ---------------------
    Ground Truth: 
    -48.0729 -47.7288 -35.9961 -24.952 -24.8682 -22.5804 -15.533 -9.50979 2.06432 2.19157 4.38561 6.10321 21.0671 28.6821 35.7077 44.6668 

    Average Execution Time for Real Roots Isolation Method: 
    1.7217e-05 secs

    Average Execution Time for Companion Matrix Method: 
    4.0405e-05 secs

    Roots Found by Real Roots Isolation Method: 
    -48.0729 -47.7288 -35.9961 -24.952 -24.8682 -22.5804 -15.533 -9.50979 2.06432 2.19157 4.38561 6.10321 21.0671 28.6821 35.7077 44.6668 

    Roots Found by  Companion Matrix Method: 
    -48.0729 -47.7288 -35.9961 -24.952 -24.8682 -22.5804 -15.533 -9.50979 2.06432 2.19157 4.38561 6.10321 21.0671 28.6821 35.7077 44.6668 

3. The Real Roots Isolation Methods uses Cauchy’s bound as well as Kojima’s bound to bracket all roots. Normally, the latter is tigher 
than the former for about 10^12 magnitude. Technically, Fujiwara’s bound is always better than Kojima's bound, while Kojima's bound 
is more numerically friendly and tight enough.

4. The Real Roots Isolation Methods uses Sturm's theorem to determine the number of roots inside any given interval.

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
