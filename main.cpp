#include "FEM.h"

int main(){
    auto a = FEM::F<float4>([](float4 x) -> float4{return -1;});
    auto b = FEM::F<float4>([](float4 x) -> float4{return 0;});
    auto c = FEM::F<float4>([](float4 x) -> float4{return -1;});
    auto f = FEM::F<float4>([](float4 x) -> float4{return 1;});

    FEM::BV<float4> left = {
        .a = 0,
        .b = 1,
        .c = 0
    };

    FEM::BV<float4> right = {
        .a = 0,
        .b = 1,
        .c = 0
    };

    FEM::Problem<float4> prob(a, b, c, f);
    auto res = prob.solve(100, 0, 1, left, right);

    //std::cout<<res;

    res.plot("file", 1000);

    return 0;
}