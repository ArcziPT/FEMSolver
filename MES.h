#ifndef MES_H
#define MES_H

#include <functional>
#include <vector>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
using float4 = boost::multiprecision::cpp_bin_float_quad;

#include <iostream>

namespace MES{
    template<typename T>
    using F = std::function<T (T)>;

    template <typename T>
    struct BasisF{
        F<T> f;
        F<T> fp;

        T a;
        T b;

        BasisF(const F<T>& f, const F<T>& fp, T a, T b): f(f), fp(fp), a(a), b(b) {}
    };

    template <typename T>
    using Basis = std::vector<BasisF<T>>;

    /**
     * Generate basis:
     *       |(x - x_(k-1))/(x_k - x_(k-1)), x_(k-1) < x <=   x_k
     * e_k = |(x_(k+1) - x)/(x_(k+1) - x_k),   x_k   < x <= x_(k+1)
     *       | 0                           ,  otherwise
     */
    template <typename T>
    Basis<T> generate_basis(int n, T x0, T xn){
        T h = (xn-x0)/n;

        auto basis_function = [h](T a, T b, T c) -> auto{
            return [a, b, c, h](T x) -> T{
                if(x <= a)
                    return 0;
                if(x >= c)
                    return 0;
                if(x <= b)
                    return (x-a)/h;
                else
                    return (c-x)/h;
            };
        };

        auto basis_function_p = [n, x0, xn](T a, T b, T c) -> auto{
            return [a, b, c, n, x0, xn](T x) -> T{
                if(x <= a)
                    return 0;
                if(x >= c)
                    return 0;
                if(x <= b)
                    return n/(xn-x0);
                else
                    return -n/(xn-x0);
            };
        };

        MES::Basis<T> basis;

        basis.push_back(BasisF<T>([x0, h](T x) -> T{
                if(x < x0)
                    return 0;
                if(x >= x0+h)
                    return 0;
                else
                    return (x0+h-x)/h;
            },[x0, h](T x) -> T{
                if(x <= x0)
                    return 0;
                if(x >= x0+h)
                    return 0;
                else
                    return -1/h;
            }, x0, x0+h));
            
        for(int i=1; i<n; i++){
            basis.push_back(BasisF<T>(basis_function(x0+(i-1)*h, x0+i*h, x0+(i+1)*h), basis_function_p(x0+(i-1)*h, x0+i*h, x0+(i+1)*h), x0+(i-1)*h, x0+(i+1)*h));
        }

        basis.push_back(BasisF<T>([xn, h](T x) -> T{
                if(x <= xn-h)
                    return 0;
                if(x > xn)
                    return 0;
                else
                    return (x-xn+h)/h;
            }, [xn, h](T x) -> T{
                if(x <= xn-h)
                    return 0;
                if(x >= xn)
                    return 0;
                else
                    return 1/h;
            }, xn-h, xn));

        return basis;
    }

    /**
     * Boundry value conditions
     * 
     * au' + bu + c = 0
     */
    template <typename T>
    struct BV{
        T a;
        T b;
        T c;

        enum Type{
            DIRICHLET,
            NEUMANN,
            ROBIN
        };

        Type get_type() const{
            if(a != 0 && b != 0)
                return ROBIN;
            if(a != 0 && b == 0)
                return NEUMANN;
            else
                return DIRICHLET;
        }
    };

    namespace ublas = boost::numeric::ublas;
    template <typename T>
    bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        // create a working copy of the input
        matrix<T> A(input);
        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());
        // perform LU-factorization
        int res = lu_factorize(A,pm);
            if( res != 0 ) return false;
        // create identity matrix of "inverse"
        inverse.assign(ublas::identity_matrix<T>(A.size1()));
        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);
        return true;
    }

    /**
     * MES:
     * (a(x)*u'(x))' + b(x)*u'(x) + c(x)*u(x) = f(x)
     */
    template <typename T>
    class Problem{
    public:
        Problem(F<T>& a, const F<T>& b, const F<T>& c, const F<T>& f): 
                a(a), b(b), c(c), f(f) {}

        void solve(int n, T a, T b, const BV<T>& left, const BV<T>& right){
            auto basis = generate_basis<T>(n, a, b);

            auto mul = [](const F<T>& u, const F<T>& v) -> auto{
                return [u, v](T x) -> T{
                    return u(x)*v(x);
                };
            };

            /**
             * MES:
             * (a(x)*u'(x))' + b(x)*u'(x) + c(x)*u(x) = f(x)
             * a(xn)*u'(xn)*v(xn) - a(x0)*u'(x0)*v(x0) + Int(...) = Int(...)
             * C1*a(xn)*u(xn)*v(xn) + C2*a(x0)*u(x0)*v(x0) + Int(...) = Int(...) + C3*a(xn)*v(xn) + C4*a(x0)*v(x0)
             */
            T C1 = 0;
            T C2 = 0;
            T C3 = 0;
            T C4 = 0;

            bool l_dirichlet = false;
            bool r_dirichlet = false;

            switch (left.get_type())
            {
            case BV<T>::DIRICHLET:
                C4 = 0;
                C2 = 0;
                l_dirichlet = true;
                break;
            case BV<T>::NEUMANN:
                basis[0].fp = [left, basis](T x) -> T{
                    return (-left.c/left.a)*basis[0].fp(x);
                };
                C4 = (left.c/left.a);
                C2 = 0;
                break;
            case BV<T>::ROBIN:
                C2 = (left.b/left.a);
                C4 = (-left.c/left.a);
                break;
            }

            switch (right.get_type())
            {
            case BV<T>::DIRICHLET:
                C3 = 0;
                C1 = 0;
                r_dirichlet = true;
                break;
            case BV<T>::NEUMANN:
                basis[n].fp = [right, basis, n](T x) -> T{
                    return (-right.c/right.a)*basis[n].fp(x);
                };
                C3 = (right.c/right.a);
                C1 = 0;
                break;
            case BV<T>::ROBIN:
                C1 = (right.b/right.a);
                C3 = (-right.c/right.a);
                break;
            }

            //std::cout<<left.get_type()<<" "<<right.get_type()<<std::endl;
            //std::cout<<C1<<" "<<C2<<" "<<C3<<" "<<C4<<std::endl;
            
            auto B = [a, b, this, mul, C1, C2, n](const BasisF<T>& u, const BasisF<T>& v) -> T{
                auto p1 = std::max(v.a, u.a);
                auto p2 = std::min(v.b, u.b);

                if(p2 <= p1)
                    return 0;

                auto i1 = boost::math::quadrature::gauss<T, 30>::integrate(mul(this->a, mul(u.fp, v.fp)), p1, p2);
                auto i2 = boost::math::quadrature::gauss<T, 30>::integrate(mul(this->b, mul(u.fp, v.f)), p1, p2);
                auto i3 = boost::math::quadrature::gauss<T, 30>::integrate(mul(this->c, mul(u.f, v.f)), p1, p2);

                auto r1 = C2 * this->a(a) * mul(u.f, v.f)(a);
                auto r2 = C1 * this->a(b) * mul(u.f, v.f)(b);

                //std::cout<<i1<<" "<<i2<<" "<<i3<<" "<<r1<<" "<<r2<<std::endl;

                return i1 + i2 + i3 + r1 + r2;
            };

            auto ud = BasisF<T>([left, right, basis, n, l_dirichlet, r_dirichlet](T x) -> T{
                if(l_dirichlet && r_dirichlet)
                    return (-left.c/left.b)*basis[0].f(x) + (-right.c/right.b)*basis[n].f(x);
                if(l_dirichlet)
                    return (-left.c/left.b)*basis[0].f(x);
                if(r_dirichlet)
                    return (-right.c/right.b)*basis[n].f(x);
            }, [left, right, basis, n, l_dirichlet, r_dirichlet](T x) -> T{
                if(l_dirichlet && r_dirichlet)
                    return (-left.c/left.b)*basis[0].fp(x) + (-right.c/right.b)*basis[n].fp(x);
                if(l_dirichlet)
                    return (-left.c/left.b)*basis[0].fp(x);
                if(r_dirichlet)
                    return (-right.c/right.b)*basis[n].fp(x);
            }, a, b);

            auto L = [a, b, this, mul, C3, C4, l_dirichlet, r_dirichlet, B, ud](const BasisF<T>& v) -> T{
                auto i = boost::math::quadrature::gauss<T, 30>::integrate(mul(this->f, v.f), v.a, v.b);

                auto r1 = C4 * this->a(a) * v.f(a);
                auto r2 = C3 * this->a(b) * v.f(b);

                if(l_dirichlet || r_dirichlet)
                    return i + r1 + r2 - B(ud, v);

                return i + r1 + r2;
            };

            boost::numeric::ublas::matrix<T> M(n+1, n+1);
            for(int i=0; i<=n; i++){
                for(int j=0; j<=n; j++){
                    if((i == 0 || j == 0) && l_dirichlet){
                        if(i == j)
                            M(j, i) = 1;
                        else
                            M(j, i) = 0;
                    }
                    else if((i == n || j == n) && r_dirichlet){
                        if(i == j)
                            M(j, i) = 1;
                        else
                            M(j, i) = 0;
                    }else
                        M(j, i) = B(basis[i], basis[j]);
                }
            }
            std::cout<<M<<std::endl;

            boost::numeric::ublas::vector<T> K(n+1);
            for(int i=0; i<=n; i++){
                if(i == 0 && l_dirichlet)
                    K(i) = 0;
                else if(i == n && r_dirichlet)
                    K(i) = 0;
                else
                    K(i) = L(basis[i]);
            }
            std::cout<<K<<std::endl;

            boost::numeric::ublas::matrix<T> MI(n+1, n+1);
            InvertMatrix(M, MI);

            std::cout<<std::setprecision(std::numeric_limits<T>::max_digits10)<<boost::numeric::ublas::prod(MI, K)<<std::endl;
        }

    private:
        F<T> a;
        F<T> b;
        F<T> c;
        F<T> f;
    };
}

MES::F<float4> F4;
MES::Basis<float4> B4;

template class MES::BasisF<float4>;
template class MES::Problem<float4>;

#endif