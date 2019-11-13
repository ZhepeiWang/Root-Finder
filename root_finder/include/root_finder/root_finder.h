#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <set>
#include <cmath>
#include <boost/math/tools/toms748_solve.hpp>
#include <Eigen/Eigen>

class RootFinder
{
public:
    // Note that tol only applies for high-order poly
    static std::set<double> solvePolyInterval(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol = 1e-3, bool isolation = true);
    static double solveGeneral(const std::function<double(double)> &func, double lbound, double ubound, double tol = 1e-3, int maxiters = -1);

    static Eigen::VectorXd polyConv(Eigen::VectorXd &lCoef, Eigen::VectorXd &rCoef);
    static double polyVal(Eigen::VectorXd &coeffs, double x);
    static Eigen::VectorXd polyModulus(const Eigen::VectorXd &u, const Eigen::VectorXd &v);

private:
    static std::set<double> solveCub(double a, double b, double c, double d);
    static std::set<double> solveQuart(double a, double b, double c, double d, double e);
    static std::set<double> solveQuartMonic(double a, double b, double c, double d);
    static int solveResolvent(double *x, double a, double b, double c);
    static std::set<double> isolateRealRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol);
    static std::set<double> eigenSolveRealRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol);
};

#endif
