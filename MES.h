#ifndef MES_H
#define MES_H

#include <functional>
#include <vector>
#include <unordered_map>

namespace MES{
    using F = std::function<double (double)>;

    struct BasisF{
        F f;
        F fp;

        double a;
        double b;

        BasisF(const F& f, const F& fp, double a, double b): f(f), fp(fp), a(a), b(b) {}
    };
    using Basis = std::vector<BasisF>;

    /**
     * Generate basis:
     *       |(x - x_(k-1))/(x_k - x_(k-1)), x_(k-1) < x <=   x_k
     * e_k = |(x_(k+1) - x)/(x_(k+1) - x_k),   x_k   < x <= x_(k+1)
     *       | 0                           ,  otherwise
     */
    Basis generate_basis(int n, double a, double b);

    /**
     * Boundry value conditions
     * 
     * au' + bu + c = 0
     */
    struct BV{
        double a;
        double b;
        double c;

        enum Type{
            DIRICHLET,
            NEUMANN,
            ROBIN
        };

        Type get_type() const;
    };

    /**
     * MES:
     * (a(x)*u'(x))' + b(x)*u'(x) + c(x)*u(x) = f(x)
     */
    class Problem{
    public:
        Problem(F& a, const F& b, const F& c, const F& f): 
                a(a), b(b), c(c), f(f) {}

        void solve(int n, double x0, double xn, Basis& basis, const BV& left, const BV& right);
        void solve(int n, double x0, double xn, const BV& left, const BV& right);

    private:
        F a;
        F b;
        F c;
        F f;
    };
}

#endif