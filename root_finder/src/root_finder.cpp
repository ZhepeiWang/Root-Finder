#include "root_finder/root_finder.h"

using namespace std;
using namespace Eigen;

set<double> RootFinder::solveCub(double a, double b, double c, double d)
{
    set<double> results;

    constexpr double cos120 = -0.50;
    constexpr double sin120 = 0.8660254037844386;

    if (fabs(d) < DBL_EPSILON)
    {
        // first solution is x = 0
        results.insert(0.0);

        // converting to quadratic equation
        d = c;
        c = b;
        b = a;
        a = 0.0;
    }

    if (fabs(a) < DBL_EPSILON)
    {
        if (fabs(b) < DBL_EPSILON)
        {
            // linear equation
            if (fabs(c) > DBL_EPSILON)
                results.insert(-d / c);
        }
        else
        {
            // quadratic equation
            double yy = c * c - 4 * b * d;
            if (yy >= 0)
            {
                double inv2b = 1 / (2 * b);
                double y = sqrt(yy);
                results.insert((-c + y) * inv2b);
                results.insert((-c - y) * inv2b);
            }
        }
    }
    else
    {
        // cubic equation
        double inva = 1 / a;
        double invaa = inva * inva;
        double bb = b * b;
        double bover3a = b * (1 / 3.0) * inva;
        double p = (3 * a * c - bb) * (1 / 3.0) * invaa;
        double halfq = (2 * bb * b - 9 * a * b * c + 27 * a * a * d) * (0.5 / 27) * invaa * inva;
        double yy = p * p * p / 27 + halfq * halfq;

        if (yy > DBL_EPSILON)
        {
            // sqrt is positive: one real solution
            double y = sqrt(yy);
            double uuu = -halfq + y;
            double vvv = -halfq - y;
            double www = fabs(uuu) > fabs(vvv) ? uuu : vvv;
            double w = (www < 0) ? -pow(fabs(www), 1 / 3.0) : pow(www, 1 / 3.0);
            results.insert(w - p / (3 * w) - bover3a);
        }
        else if (yy < -DBL_EPSILON)
        {
            // sqrt is negative: three real solutions
            double x = -halfq;
            double y = sqrt(-yy);
            double theta;
            double r;
            double ux;
            double uyi;

            if (fabs(x) > DBL_EPSILON)
            {
                theta = (x > 0) ? atan(y / x) : (atan(y / x) + M_PI);
                r = sqrt(x * x - yy);
            }
            else
            {
                theta = M_PI / 2;
                r = y;
            }

            theta /= 3.0;
            r = pow(r, 1 / 3.0);
            ux = cos(theta) * r;
            uyi = sin(theta) * r;
            // first solution
            results.insert(ux + ux - bover3a);
            // second solution, rotate +120 degrees
            results.insert(2 * (ux * cos120 - uyi * sin120) - bover3a);
            // third solution, rotate -120 degrees
            results.insert(2 * (ux * cos120 + uyi * sin120) - bover3a);
        }
        else
        {
            // sqrt is zero: two real solutions
            double www = -halfq;
            double w = (www < 0) ? -pow(fabs(www), 1 / 3.0) : pow(www, 1 / 3.0);
            // first solution
            results.insert(w + w - bover3a);
            // second solution, rotate +120 degrees
            results.insert(2 * w * cos120 - bover3a);
        }
    }
    return results;
}

set<double> RootFinder::solveQuart(double a, double b, double c, double d, double e)
{
    if (fabs(a) < DBL_EPSILON)
    {
        return solveCub(b, c, d, e);
    }
    else
    {
        return solveQuartMonic(b / a, c / a, d / a, e / a);
    }
}

set<double> RootFinder::solvePolyInterval(const VectorXd &coeffs, double lbound, double ubound, double tol, bool isolation)
{
    set<double> rts;

    int valid = coeffs.size();
    for (int i = 0; i < coeffs.size(); i++)
    {
        if (fabs(coeffs(i)) < DBL_EPSILON)
        {
            valid--;
        }
        else
        {
            break;
        }
    }

    int offset = 0;
    int nonzeros = valid;
    if (valid > 0)
    {
        for (int i = 0; i < valid; i++)
        {
            if (fabs(coeffs(coeffs.size() - i - 1)) < DBL_EPSILON)
            {
                nonzeros--;
                offset++;
            }
            else
            {
                break;
            }
        }
    }

    if (nonzeros == 0)
    {
        rts.insert(INFINITY);
        rts.insert(-INFINITY);
    }
    else if (nonzeros == 1 && offset == 0)
    {
        rts.clear();
    }
    else
    {
        VectorXd ncoeffs(std::max(5, nonzeros));
        ncoeffs.setZero();
        ncoeffs.tail(nonzeros) << coeffs.segment(coeffs.size() - valid, nonzeros);

        if (nonzeros <= 5)
        {
            rts = solveQuart(ncoeffs(0), ncoeffs(1), ncoeffs(2), ncoeffs(3), ncoeffs(4));

            for (auto it = rts.begin(); it != rts.end();)
            {
                if (*it >= lbound && *it <= ubound)
                {
                    it++;
                }
                else
                {
                    it = rts.erase(it);
                }
            }
        }
        else
        {
            if (isolation)
            {
                rts = isolateRealRoots(ncoeffs, lbound, ubound, tol);
            }
            else
            {
                rts = eigenSolveRealRoots(ncoeffs, lbound, ubound, tol);
            }
        }

        if (offset > 0)
        {
            rts.insert(0);
        }
    }

    return rts;
}

double RootFinder::solveGeneral(const function<double(double)> &func, double lbound, double ubound, double tol, int maxiters)
{
    double rt = NAN;

    using namespace boost::math::tools;

    auto tolf = [&](double a, double b) { return fabs(a - b) < tol; };

    boost::uintmax_t maxit = maxiters >= 0 ? maxiters : 2 * std::max(1, 1 + (int)log2(fabs((ubound - lbound) / tol)));

    double fl = func(lbound);
    double fu = func(ubound);

    if (fl == 0)
    {
        rt = lbound;
    }

    if (fu == 0)
    {
        rt = ubound;
    }

    if (fl * fu < 0)
    {
        pair<double, double> result = toms748_solve(func, lbound, ubound, tolf, maxit);
        rt = (result.first + result.second) / 2.0;
    }

    return rt;
}

VectorXd RootFinder::polyConv(VectorXd &lCoef, VectorXd &rCoef)
{
    VectorXd result(lCoef.size() + rCoef.size() - 1);
    result.setZero();
    for (int i = 0; i < result.size(); i++)
    {
        for (int j = 0; j <= i; j++)
        {
            result(i) += (j < lCoef.size() && (i - j) < rCoef.size()) ? (lCoef(j) * rCoef(i - j)) : 0;
        }
    }

    return result;
}

double RootFinder::polyVal(VectorXd &coeffs, double x)
{
    double xn = 1.0;
    double retVal = 0.0;

    for (int i = (int)coeffs.size() - 1; i >= 0; i--)
    {
        retVal += coeffs(i) * xn;
        xn *= x;
    }

    return retVal;
}

set<double> RootFinder::solveQuartMonic(double a, double b, double c, double d)
{
    set<double> results;

    double a3 = -b;
    double b3 = a * c - 4. * d;
    double c3 = -a * a * d - c * c + 4. * b * d;

    // cubic resolvent
    // y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0

    double x3[3];
    int iZeroes = solveResolvent(x3, a3, b3, c3);

    double q1, q2, p1, p2, D, sqD, y;

    y = x3[0];
    // choosing Y with maximal absolute value.
    if (iZeroes != 1)
    {
        if (fabs(x3[1]) > fabs(y))
            y = x3[1];
        if (fabs(x3[2]) > fabs(y))
            y = x3[2];
    }

    // h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

    D = y * y - 4 * d;
    if (fabs(D) < DBL_EPSILON) //in other words - D==0
    {
        q1 = q2 = y * 0.5;
        // g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
        D = a * a - 4 * (b - y);
        if (fabs(D) < DBL_EPSILON) //in other words - D==0
            p1 = p2 = a * 0.5;

        else
        {
            sqD = sqrt(D);
            p1 = (a + sqD) * 0.5;
            p2 = (a - sqD) * 0.5;
        }
    }
    else
    {
        sqD = sqrt(D);
        q1 = (y + sqD) * 0.5;
        q2 = (y - sqD) * 0.5;
        // g1+g2 = a && g1*h2 + g2*h1 = c   ( && g === p )  Krammer
        p1 = (a * q1 - c) / (q1 - q2);
        p2 = (c - a * q2) / (q1 - q2);
    }

    // solving quadratic eq. - x^2 + p1*x + q1 = 0
    D = p1 * p1 - 4 * q1;
    if (fabs(D) < DBL_EPSILON)
    {
        results.insert(-p1 * 0.5);
    }
    else if (D > 0.0)
    {
        sqD = sqrt(D);
        results.insert((-p1 + sqD) * 0.5);
        results.insert((-p1 - sqD) * 0.5);
    }

    // solving quadratic eq. - x^2 + p2*x + q2 = 0
    D = p2 * p2 - 4 * q2;
    if (fabs(D) < DBL_EPSILON)
    {
        results.insert(-p2 * 0.5);
    }
    else if (D > 0.0)
    {
        sqD = sqrt(D);
        results.insert((-p2 + sqD) * 0.5);
        results.insert((-p2 - sqD) * 0.5);
    }

    return results;
}

int RootFinder::solveResolvent(double *x, double a, double b, double c)
{
    double a2 = a * a;
    double q = (a2 - 3 * b) / 9;
    double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
    double r2 = r * r;
    double q3 = q * q * q;
    double A, B;
    if (r2 < q3)
    {
        double t = r / sqrt(q3);
        if (t < -1)
            t = -1;
        if (t > 1)
            t = 1;
        t = acos(t);
        a /= 3;
        q = -2 * sqrt(q);
        x[0] = q * cos(t / 3) - a;
        x[1] = q * cos((t + M_PI * 2) / 3) - a;
        x[2] = q * cos((t - M_PI * 2) / 3) - a;
        return 3;
    }
    else
    {
        A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3);
        if (r < 0)
            A = -A;
        B = (0 == A ? 0 : q / A);

        a /= 3;
        x[0] = (A + B) - a;
        x[1] = -0.5 * (A + B) - a;
        x[2] = 0.5 * sqrt(3.) * (A - B);
        if (fabs(x[2]) < DBL_EPSILON)
        {
            x[2] = x[1];
            return 2;
        }

        return 1;
    }
}

set<double> RootFinder::isolateRealRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol)
// leading coefficient must be nonzero
{
    set<double> rts;

    int order = (int)coeffs.size() - 1;
    Eigen::VectorXd monicCoeffs(order + 1);
    monicCoeffs << 1.0, coeffs.tail(order) / coeffs(0);

    // calculate Cauchy’s bound for the roots of a polynomial
    double rho_c = 1 + monicCoeffs.tail(order).cwiseAbs().maxCoeff();

    // calculate Kojima’s bound for the roots of a polynomial
    VectorXd nonzeroCoeffs(order + 1);
    nonzeroCoeffs.setZero();
    int nonzeros = 0;
    double tempEle;
    for (int i = 0; i < order + 1; i++)
    {
        tempEle = monicCoeffs(i);
        if (fabs(tempEle) >= DBL_EPSILON)
        {
            nonzeroCoeffs(nonzeros++) = tempEle;
        }
    }
    nonzeroCoeffs = nonzeroCoeffs.head(nonzeros).eval();
    VectorXd kojimaVec = nonzeroCoeffs.tail(nonzeros - 1).cwiseQuotient(nonzeroCoeffs.head(nonzeros - 1)).cwiseAbs();
    kojimaVec.tail(1) /= 2.0;
    double rho_k = 2.0 * kojimaVec.maxCoeff();

    // compare the bound and then loosen it by 1.0
    double rho = std::min(rho_c, rho_k) + 1.0;

    // reasonable tries for iterations
    int maxits = 2 * std::max(1, 1 + (int)log2(rho / tol));

    // tighten the bound to search in
    lbound = std::max(lbound, -rho);
    ubound = std::min(ubound, rho);

    // build sturm sequence
    vector<VectorXd> sturmSeqs;
    sturmSeqs.emplace_back(order + 1);
    sturmSeqs.back() << monicCoeffs;

    sturmSeqs.emplace_back(order);
    for (int i = 0; i < order; i++)
    {
        sturmSeqs.back()(i) = sturmSeqs.front()(i) * (order - i);
    }
    sturmSeqs.back() = (sturmSeqs.back() / fabs(sturmSeqs.back()(0))).eval();

    bool remainderConstant = false;
    VectorXd remainder;
    while (!remainderConstant)
    {
        remainder = polyModulus(*(sturmSeqs.end() - 2), sturmSeqs.back());
        remainderConstant = remainder.size() == 1;
        sturmSeqs.emplace_back(remainder / (-fabs(remainder(0))));
    }

    // useful functors
    auto func = [&](double t) {
        double tn = 1.0, retVal = 0.0;
        for (int i = (int)monicCoeffs.size() - 1; i >= 0; i--)
        {
            retVal += tn * monicCoeffs(i);
            tn *= t;
        }
        return retVal;
    };

    auto dfunc = [&](double t) {
        double tn = 1.0, retVal = 0.0;
        for (int i = (int)sturmSeqs[1].size() - 1; i >= 0; i--)
        {
            retVal += tn * sturmSeqs[1](i);
            tn *= t;
        }
        return retVal;
    };

    auto nchanges = [&](double x) {
        double y, lasty;
        int num = 0;
        lasty = polyVal(sturmSeqs[0], x);
        for (auto it = ++sturmSeqs.begin(); it != sturmSeqs.end(); it++)
        {
            y = polyVal(*it, x);
            if (lasty == 0.0 || lasty * y < 0.0)
            {
                ++num;
            }
            lasty = y;
        }
        return num;
    };

    // isolate only handles open interval, two end points are not considered
    function<void(double, double, int, int)> isolate = [&](double l, double r, int lchanges, int rchanges) {
        int nrts = lchanges - rchanges;

        if (nrts == 0)
        {
            return;
        }
        else if (nrts == 1)
        {
            if (func(l) * func(r) < 0)
            {
                rts.insert(solveGeneral(func, l, r, tol));
                return;
            }
            else
            {
                // binsect when non of above works
                double m;
                for (int i = 0; i < maxits; i++)
                {
                    // root with even multiplicity
                    if (dfunc(l) * dfunc(r) < 0)
                    {
                        rts.insert(solveGeneral(dfunc, l, r, tol));
                        return;
                    }

                    m = (l + r) / 2.0;
                    if (func(m) == 0 || fabs(r - l) < tol)
                    {
                        rts.insert(m);
                        return;
                    }
                    else
                    {
                        if (lchanges == nchanges(m))
                        {
                            l = m;
                        }
                        else
                        {
                            r = m;
                        }
                    }
                }
                rts.insert(m);
                return;
            }
        }
        else if (nrts > 1)
        {
            // more than one root in the interval
            double m;
            int mchanges;
            for (int i = 0; i < maxits; i++)
            {
                m = (l + r) / 2.0;
                mchanges = nchanges(m);
                if (func(m) == 0)
                {
                    rts.insert(m);
                    return;
                }
                else if (lchanges != mchanges && rchanges != mchanges)
                {
                    isolate(l, m, lchanges, mchanges);
                    isolate(m, r, mchanges, rchanges);
                    return;
                }
                else if (lchanges == mchanges)
                {
                    l = m;
                }
                else
                {
                    r = m;
                }
            }
            rts.insert(m);
            return;
        }
    };

    isolate(lbound, ubound, nchanges(lbound), nchanges(ubound));

    return rts;
}

set<double> RootFinder::eigenSolveRealRoots(const Eigen::VectorXd &coeffs, double lbound, double ubound, double tol)
{
    set<double> rts;

    int order = (int)coeffs.size() - 1;
    Eigen::VectorXd monicCoeffs(order + 1);
    monicCoeffs << 1.0, coeffs.tail(order) / coeffs(0);

    MatrixXd companionMat(order, order);
    companionMat.setZero();
    companionMat(0, order - 1) = -monicCoeffs(order);
    for (int i = 1; i < order; i++)
    {
        companionMat(i, i - 1) = 1.0;
        companionMat(i, order - 1) = -monicCoeffs(order - i);
    }
    VectorXcd eivals = companionMat.eigenvalues();
    double real;
    for (int i = 0; i < (int)eivals.size(); i++)
    {
        real = eivals(i).real();
        if (eivals(i).imag() < tol && real > lbound && real < ubound)
            rts.insert(real);
    }

    return rts;
}

VectorXd RootFinder::polyModulus(const VectorXd &u, const VectorXd &v)
// modulus of u(x) / v(x), the leading coefficient of v must be 1 or -1
{
    int orderu = (int)u.size() - 1;
    int orderv = (int)v.size() - 1;

    VectorXd r = u;

    if (v(0) < 0)
    {
        for (int i = orderv + 1; i <= orderu; i += 2)
        {
            r(i) = -r(i);
        }
        for (int i = 0; i <= orderu - orderv; i++)
        {
            for (int j = i + 1; j <= orderv + i; j++)
            {
                r(j) = -r(j) - r(i) * v(j - i);
            }
        }
    }
    else
    {
        for (int i = 0; i <= orderu - orderv; i++)
        {
            for (int j = i + 1; j <= orderv + i; j++)
            {
                r(j) = r(j) - r(i) * v(j - i);
            }
        }
    }

    int k = orderv - 1;
    while (k >= 0 && fabs(r(orderu - k)) < DBL_EPSILON)
    {
        r(orderu - k) = 0.0;
        k--;
    }

    return r.tail((k <= 0) ? 1 : (k + 1));
}
