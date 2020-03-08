# Root-Finder
_Root-Finder_ is a header-only univariate polynomial solver, which finds/counts all REAL roots of any polynomial inside a given interval.

# Feature

0. The solver is a header-only C++11 lib, which is highly optimized on the premise of instruction set independence.

1. Only one header file "root_finder.hpp".

2. The interface of this lib only contains two functions. One is for _roots finding_. The other is for _roots counting_.

3. As for low order polynomials (linear, quadratic, cubic and quartic polynomials), the solver use closed form solution.
In this case, the solver only takes about 0.4E-6 sec.

4. As for high order polynomials (order >= 5), the solver implements 2 different methods to find all roots. The recommended 
one, named _Real Roots Isolation Method_. The other one is based on _Companion Matrix Method_. The example is the comparision 
between these two methods.

5. The _Real Roots Isolation Method_ uses Cauchy’s bound as well as Kojima’s bound to bracket all roots. Normally, the latter 
can be tighter than the former for several magnitude in most cases. Technically, Fujiwara’s bound is always better than 
Kojima's bound, while Kojima's bound is more numerically friendly and is tight enough.

6. In _Real Roots Isolation Method_, Sturm theory is employed to bracket each single root. Then _TOMS748_ is employed to shrink 
the interval efficiently

7. The _Real Roots Isolation Method_ is much faster and much more stable than the _Companion Matrix Method_. However, due to 
truncation error of float point number, the former is recommended for at most 32-order polynomials, while the latter is only 
recommended for at-most 20-order polynomials.

8. We provide benchmark example between our _Real Roots Isolation Method_ and _TOMS493: Jenkins–Traub Algorithm_ on two different 
platforms. The latter one is commonly known as the "RPOLY" algorithm. For 8-order polynomials, our method is about 27% faster than 
"RPOLY" under Intel i5-5200U CPU, while 13% slower under Intel i7-8700 CPU. In general, out lib has comparably low time comsumption 
as the widely employed "RPOLY" algorithm, in terms of real roots finding. 

9. Out lib is also capable of efficiently counting the number of roots inside an interval, of which "RPOLY" is incapable. For 8-order 
polynomials, our roots counter only takes about 0.2E-6 sec under Intel i7-8700 CPU and 0.8E-6 under Intel i5-5200U CPU. The time 
consumption for high order polynomial is extremely low as for those low order polynomials which have closed-form solutions.

10. Furthur instruction set independent optimization can be done to get better performance.

# Interface

Only two functions below is needed.

_Roots Finding Function_:

    std::set<double> RootFinder::solvePolynomial(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol, bool isolation = true)

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

Example:
    
    Eigen::VectorXd coeffs(6);
    coeffs << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0;
    std::set<double> allRoots = RootFinder::solvePolynomial(coeffs, -100.0, 100.0);

_Roots Counting Function_:

    int RootFinder::countRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound)

Inputs:

    coeffs: 
        Eigen VectorXd for coefficients of a polynomial. For example, the polynomial a(n)*x^n+a(n-1)*x^(n-1)+...+a(1)*x+a(0) can be expressed by 
        coeffs=[a(n), a(n-1), ..., a(1), a(0)].

    lbound and ubound:
        Open interval (lbound, ubound). Note that polynomial cannot be zeros at these two boundaries.

Outputs:

    The number of distinct roots inside (lbound, ubound).

Example:
    
    Eigen::VectorXd coeffs(6);
    coeffs << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0;
    int numRoots = RootFinder::countRoots(coeffs, -10, 0.5);

# Compile

    sudo apt install libeigen3-dev libboost-dev
    mkdir build
    cd build
    cmake ..
    make

# Note

0. Thanks to The "RPOLY" algorithm benchmarked comes a modified version by [Helen Oleynikova](https://github.com/helenol), 
which originates in http://www.netlib.org/toms/.

1. This lib is employed in a trajectory optimization/motion planning project of [ZJU-FAST-Lab](https://github.com/ZJU-FAST-Lab). 
See [AM-Traj](https://github.com/ZJU-FAST-Lab/am_traj) for details.
