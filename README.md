# FEM Solver

Finite Element Method implementation for variational characterization of equation:
![alt text](eq.png?raw=true "equation")
Where:
```
u - unknown function
a, b, c, f - functions provided by user
A1, A2 - additional terms resulting from boundry conditions (automaticly calculated)
```
With boundry conditions:
```
au' + bu + c = 0
```
Where:
```
a, b, c - user provided coefficients
```

### Dependencies

* [Boost C++](https://www.boost.org/) - integral calculations, (math::quadrature::gauss, multiprecision::cpp_bin_float (might be removed manually))
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - matrix/vector calculations
* [matplotlib-cpp](https://github.com/lava/matplotlib-cpp) - plotting

### Usage
This example solves equation:
![alt text](ex_eq.png?raw=true "equation")

Include header `"FEM.h"`.

```c++
//functions
auto a = FEM::F<double>([](double x) -> double{return x;});
auto b = FEM::F<double>([](double x) -> double{return 1;});
auto c = FEM::F<double>([](double x) -> double{return 0;});
auto f = FEM::F<double>([](double x) -> double{return -std::exp(x);});

//boundry conditions
FEM::BV<float4> left = {
        .a = 0,
        .b = 1,
        .c = -2
    };

FEM::BV<float4> right = {
    .a = 1,
    .b = -1,
    .c = -1
};

MES::Problem<double> prob(a, b, c, f);
auto res = prob.solve(100, 0, 1); //solve

res.plot("plot.png", 10000); //open and save plot with 10000 points
```