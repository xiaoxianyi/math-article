#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cstdlib>

struct Individual {
    std::vector<double> vars;  // T_in, t
    std::vector<double> objs;  // 舒适度, 能耗
    int rank;
    double crowding;
};

// 自定义clamp函数（兼容C++11）
double clamp(double value, double min, double max) {
    return (value < min) ? min : (value > max) ? max : value;
}

void evalObj(Individual& ind) {
    double T_in = ind.vars[0];
    double t = ind.vars[1];
    double T_out = 29.5;
    double k = 0.35;

    // 修改舒适度公式使最大值接近 1
    double comfort = 1.0 / (std::abs(T_in - 25) + 1e-6); // +1e-6避免除以零
    comfort = clamp(comfort, 0.0, 1.0); // 使用自定义clamp

    double energy = k * std::abs(T_in - T_out) * t;

    ind.objs.clear();
    ind.objs.push_back(-comfort); // 仍需要最小化 -comfort（即最大化 comfort）
    ind.objs.push_back(energy);
}

void crossover(Individual& p1, Individual& p2, Individual& c1, Individual& c2) {
    double rate = 0.9;
    if (rand() / (double)RAND_MAX < rate) {
        double alpha = rand() / (double)RAND_MAX;
        for (size_t i = 0; i < p1.vars.size(); ++i) {
            c1.vars[i] = alpha * p1.vars[i] + (1 - alpha) * p2.vars[i];
            c2.vars[i] = alpha * p2.vars[i] + (1 - alpha) * p1.vars[i];
        }
    } else {
        c1 = p1;
        c2 = p2;
    }

    // 添加变量范围约束（使用自定义clamp）
    c1.vars[0] = clamp(c1.vars[0], 20.0, 30.0);
    c1.vars[1] = clamp(c1.vars[1], 14.0, 15.0);
    c2.vars[0] = clamp(c2.vars[0], 20.0, 30.0);
    c2.vars[1] = clamp(c2.vars[1], 14.0, 15.0);
}

void mutate(Individual& ind) {
    double rate = 0.1;
    for (size_t i = 0; i < ind.vars.size(); ++i) {
        if (rand() / (double)RAND_MAX < rate) {
            ind.vars[i] += (rand() / (double)RAND_MAX - 0.5) * 0.1;
        }
    }

    // 严格约束范围
    ind.vars[0] = clamp(ind.vars[0], 20.0, 30.0);
    ind.vars[1] = clamp(ind.vars[1], 14.0, 15.0);
}

std::vector<Individual> nsga2(size_t popSize, size_t maxGen) {
    std::vector<Individual> pop(popSize); // 修复pop声明

    // 初始化种群
    for (size_t i = 0; i < popSize; ++i) {
        pop[i].vars.push_back(20 + (rand() / (double)RAND_MAX) * 10); // T_in ∈ [20,30]
        pop[i].vars.push_back(14 + (rand() / (double)RAND_MAX) * 1);  // t ∈ [14,15]
        evalObj(pop[i]);
    }

    for (size_t gen = 0; gen < maxGen; ++gen) {
        // ... 原有nsga2逻辑保持不变 ...
    }

    return pop;
}

int main() {
    srand(time(NULL));  
    size_t popSize = 1000;
    size_t maxGen = 100;

    std::vector<Individual> result = nsga2(popSize, maxGen);

std::cout << "\nPareto 前沿（共 " << result.size() << " 个解）：" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        double comfort = -result[i].objs[0];
        std::cout << "解 " << i+1 << ": T_in = " << result[i].vars[0] 
                  << "°C, t = " << result[i].vars[1] << "h, "
                  << "舒适度 = " << comfort << ", 能耗 = " << result[i].objs[1] << std::endl;
    }

    // 输出舒适度接近1的解
    std::cout << "\n舒适度接近1的解：" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        double comfort = -result[i].objs[0];
        if (comfort > 0.9) {
            std::cout << "解 " << i+1 << ": T_in = " << result[i].vars[0] 
                      << "°C, t = " << result[i].vars[1] << "h, "
                      << "舒适度 = " << comfort << ", 能耗 = " << result[i].objs[1] << std::endl;
        }
    }
    return 0;
}
