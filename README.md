# Root-Finder
__Root-Finder__ is a header-only univariate polynomial solver, which finds/counts all REAL roots of any polynomial inside a given interval.

# Feature

0. The solver is a __C++11__ [__header-only__](https://en.wikipedia.org/wiki/Header-only) library, which is highly optimized on the 
premise of __instruction set independence__.

1. It only contains one header file "[__root_finder.hpp__](https://github.com/ZhepeiWang/Root-Finder/blob/master/root_finder/include/root_finder/root_finder.hpp)".

2. The interface only contains two functions. One is for __roots finding__ while the other one is for __roots counting__.

3. As for low order polynomials (linear, quadratic, cubic, and quartic polynomials), the solver calculates their closed-form solutions. 
In this case, the solver only takes about 0.4 microsecond.

4. As for high order polynomials (order >= 5), the solver implements 2 different methods to find all roots. The recommended 
one is __Real Roots Isolation Method__. The other one is based on [__Companion Matrix Method__](https://en.wikipedia.org/wiki/Companion_matrix). 

5. The __Real Roots Isolation Method__ uses [geometrical properties](https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots) 
of polynomial roots to tighten a given interval. Both Cauchy’s bound as well as Kojima’s bound are adopted. Normally, the latter 
can be tighter than the former for several magnitude in most cases. Technically, [Fujiwara’s bound is always better than 
Kojima's bound](https://doi.org/10.1016/j.cam.2003.10.019), while Kojima's bound is more numerically friendly and is tight enough.

6. In __Real Roots Isolation Method__, [__Sturm theory__](https://link.springer.com/book/10.1007%2F978-3-662-05355-3) is employed to bracket each 
single root. Then [__TOMS748__](https://doi.org/10.1145/210089.210111) is employed to shrink the interval efficiently.

7. The __Real Roots Isolation Method__ is much faster and much more stable than the __Companion Matrix Method__. However, due to 
truncation error of float point number, the former is recommended for at most 32-order polynomials, while the latter is only 
recommended for at most 20-order polynomials.

8. We perform benchmark example between our __Real Roots Isolation Method__ and 
[__TOMS493: Jenkins–Traub Algorithm__](https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm) on two different platforms. The latter one is commonly 
known as the "RPOLY" algorithm. In the benchmark, all roots of a polynomial are required. For 8-order polynomials, 
our method is about 27% faster than RPOLY under Intel i5-5200U CPU, while 13% slower under Intel i7-8700 CPU. 
In general, out library has comparably low time comsumption as the widely employed RPOLY algorithm, in terms of all real 
roots finding. Moreover, RPOLY is designed to find all roots of a polynomial while ours can handle any given interval. 
Therefore, when roots in an interval are required, our method performs far better.

9. The library is also capable of efficiently counting the number of roots inside an interval, of which RPOLY is incapable. For 8-order 
polynomials, our roots counter only takes about 0.2 microsecond under Intel i7-8700 CPU and 0.8 microsecond under Intel i5-5200U CPU. The time 
consumption for high order polynomial is extremely low as for those low order polynomials which have closed-form solutions.

10. Furthur instruction set independent optimization can be done to get better performance.

# Interface

Only two functions below are needed.

__Roots Finding Function__:

    std::set<double> RootFinder::solvePolynomial(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol, bool isolation = true)

Inputs:

    coeffs: 
        Eigen VectorXd for coefficients of a polynomial. For example, the polynomial a(n)*x^n+a(n-1)*x^(n-1)+...+a(1)*x+a(0) can be expressed by 
        coeffs=[a(n), a(n-1), ..., a(1), a(0)].

    lbound and ubound:
        Open interval (lbound, ubound). If lbound = -INFINITY and ubound = INFINITY, then all real roots can be found by the solver. Note that polynomial cannot be zero at these two boundaries.
    
    tol:
        tolerance for precision, only works when order is higher than 4.
    
    isolation:
        switch for Method used. Default one is Real Roots Isolation Method.

Outputs:

    std::set<double> which stores all searched real roots.

Example:
    
    Eigen::VectorXd coeffs(6);
    coeffs << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0;
    std::set<double> allRoots = RootFinder::solvePolynomial(coeffs, -100.0, 100.0, 0.00001);

__Roots Counting Function__:

    int RootFinder::countRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound)

Inputs:

    coeffs: 
        Eigen VectorXd for coefficients of a polynomial. For example, the polynomial a(n)*x^n+a(n-1)*x^(n-1)+...+a(1)*x+a(0) can be expressed by 
        coeffs=[a(n), a(n-1), ..., a(1), a(0)].

    lbound and ubound:
        Open interval (lbound, ubound). Note that polynomial cannot be zero at these two boundaries.

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

0. The "RPOLY" algorithm benchmarked comes from a modified version by [Helen Oleynikova](https://github.com/helenol), 
which originates in http://www.netlib.org/toms/.

1. This lib is employed in a trajectory optimization/motion planning project of [ZJU-FAST-Lab](https://github.com/ZJU-FAST-Lab). 
See [AM-Traj](https://github.com/ZJU-FAST-Lab/am_traj) for details.
