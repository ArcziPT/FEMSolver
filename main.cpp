#include "MES.h"

int main(){
    auto a = MES::F<float4>([](float4 x) -> float4{return -1;});
    auto b = MES::F<float4>([](float4 x) -> float4{return 0;});
    auto c = MES::F<float4>([](float4 x) -> float4{return 1;});
    auto f = MES::F<float4>([](float4 x) -> float4{return 1;});

    MES::BV<float4> left = {
        .a = 0,
        .b = 1,
        .c = -2
    };

    MES::BV<float4> right = {
        .a = -1,
        .b = 1,
        .c = -1
    };

    MES::Problem<float4> prob(a, b, c, f);
    prob.solve(3, 0, 1, left, right);

    return 0;
}