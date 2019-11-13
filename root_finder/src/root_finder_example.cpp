#include "root_finder/root_finder.h"

using namespace std;
using namespace Eigen;

// use Mersenne Twister as PRNG
std::mt19937_64 gen;
std::uniform_real_distribution<double> uniformReal(0.0, 1.0);

void testRootFinder(int order)
{
    // randomly geenrate polynomial whose all real roots lie in interval [-50.0, 50.0]
    set<double> realRoots;
    VectorXd coeffs(1);
    coeffs << 1;
    double tempRoot;
    for (int i = 1; i <= order; i++)
    {
        VectorXd b(2);
        tempRoot = (uniformReal(gen) - 0.5) * 100;
        b << 1, -tempRoot;
        coeffs = RootFinder::polyConv(coeffs, b);
        realRoots.insert(tempRoot);
    }

    // show ground truth
    cout << "--------------------- Speed Test between Roots Isolation Method and Companion Matrix Method ---------------------" << endl;

    cout << "Ground Truth: " << endl;
    for (auto it = realRoots.begin(); it != realRoots.end(); it++)
    {
        cout << *it << " ";
    }
    cout << endl
         << endl;

    // running root finder for 1000 rounds
    int numRounds = 1000;
    set<double> rts1, rts2;
    clock_t ini = clock();
    for (int i = 0; i < numRounds; i++)
    {
        //find all roots by real roots isolation, i.e., in (-inf, inf)
        rts1 = RootFinder::solvePolyInterval(coeffs, -INFINITY, INFINITY, 1e-7);
    }
    clock_t mid = clock();
    for (int i = 0; i < numRounds; i++)
    {
        //find all roots by eigen values of companion matrix
        rts2 = RootFinder::solvePolyInterval(coeffs, -INFINITY, INFINITY, 1e-7, false);
    }
    clock_t fin = clock();

    // results
    cout << "Average Execution Time for Real Roots Isolation Method: " << endl
         << (mid - ini) * 1.0 / numRounds / CLOCKS_PER_SEC << " secs" << endl
         << endl;

    cout << "Average Execution Time for Companion Matrix Method: " << endl
         << (fin - mid) * 1.0 / numRounds / CLOCKS_PER_SEC << " secs" << endl
         << endl;

    cout << "Roots Found by Real Roots Isolation Method: " << endl;
    for (auto it = rts1.begin(); it != rts1.end(); it++)
    {
        cout << setprecision(6) << *it << " ";
    }
    cout << endl
         << endl;
    cout << "Roots Found by  Companion Matrix Method: " << endl;
    for (auto it = rts2.begin(); it != rts2.end(); it++)
    {
        cout << setprecision(6) << *it << " ";
    }
    cout << endl
         << endl;
}

int main(int argc, char **argv)
{

    testRootFinder(16);

    return 0;
}
