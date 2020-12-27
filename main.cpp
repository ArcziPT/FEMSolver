#include "MES.h"

int main(){
    auto a = MES::F([](double x) -> double{return -1;});
    auto b = MES::F([](double x) -> double{return 0;});
    auto c = MES::F([](double x) -> double{return 1;});
    auto f = MES::F([](double x) -> double{return 1;});

    MES::BV left = {
        .a = 0,
        .b = 1,
        .c = -2
    };

    MES::BV right = {
        .a = -1,
        .b = 1,
        .c = -1
    };

    MES::Problem prob(a, b, c, f);
    prob.solve(3, 0, 1, left, right);

    return 0;
}