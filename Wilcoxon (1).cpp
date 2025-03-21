#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cstdlib>

// 个体结构
struct Individual {
    std::vector<double> vars;  // T_in, t
    std::vector<double> objs;  // 舒适度, 能耗
    int rank;
    double crowding;

    Individual() : rank(0), crowding(0.0) {
        vars.resize(2);  // 确保 vars 有 2 个元素
        objs.resize(2);  // 确保 objs 有 2 个元素
    }
};

// 自定义clamp函数（兼容C++98）
double clamp(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

// 评估目标函数
void evalObj(Individual& ind) {
    double T_in = ind.vars[0];
    double t = ind.vars[1];
    double T_out = 29.5;
    double k = 0.35;

    // 修改舒适度公式使最大值接近 1
    double comfort = 1.0 / (std::abs(T_in - 25) + 1e-6); // +1e-6避免除以零
    comfort = clamp(comfort, 0.0, 1.0); // 使用自定义clamp

    double energy = k * std::abs(T_in - T_out) * t;

    ind.objs[0] = -comfort; // 仍需要最小化 -comfort（即最大化 comfort）
    ind.objs[1] = energy;
}

// 交叉操作
void crossover(const Individual& p1, const Individual& p2, Individual& c1, Individual& c2) {
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

// 变异操作
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

// Wilcoxon 秩和检验比较函数
bool wilcoxonCompare(const Individual& a, const Individual& b) {
    // 比较两个个体的目标函数值
    // 如果 a 的目标函数值优于 b，则返回 true
    // 这里假设目标是最小化两个目标函数值
    if (a.objs[0] < b.objs[0] && a.objs[1] < b.objs[1]) {
        return true;
    }
    return false;
}

// 环境选择（使用 Wilcoxon 方法）
void environmentalSelection(std::vector<Individual>& pop, const std::vector<Individual>& offspring, size_t popSize) {
    std::vector<Individual> combinedPop = pop;
    combinedPop.insert(combinedPop.end(), offspring.begin(), offspring.end());

    // 使用 Wilcoxon 方法比较个体
    std::sort(combinedPop.begin(), combinedPop.end(), wilcoxonCompare);

    // 选择前 popSize 个个体
    pop.clear();
    for (size_t i = 0; i < popSize && i < combinedPop.size(); ++i) {
        pop.push_back(combinedPop[i]);
    }
}

// NSGA-II 主函数
std::vector<Individual> nsga2(size_t popSize, size_t maxGen) {
    std::vector<Individual> pop(popSize);

    // 初始化种群
    for (size_t i = 0; i < popSize; ++i) {
        pop[i].vars[0] = 20 + (rand() / (double)RAND_MAX) * 10; // T_in ∈ [20,30]
        pop[i].vars[1] = 14 + (rand() / (double)RAND_MAX) * 1;  // t ∈ [14,15]
        evalObj(pop[i]);
    }

    for (size_t gen = 0; gen < maxGen; ++gen) {
        // 生成子代
        std::vector<Individual> offspring;
        for (size_t i = 0; i < popSize - 1; i += 2) {  // 确保 i + 1 不越界
            Individual c1, c2;
            crossover(pop[i], pop[i + 1], c1, c2);
            mutate(c1);
            mutate(c2);
            evalObj(c1);
            evalObj(c2);
            offspring.push_back(c1);
            offspring.push_back(c2);
        }

        // 环境选择
        environmentalSelection(pop, offspring, popSize);
    }

    return pop;
}

int main() {
    srand(static_cast<unsigned int>(time(NULL)));  // 初始化随机种子
    size_t popSize = 100;
    size_t maxGen = 100;

    std::vector<Individual> result = nsga2(popSize, maxGen);

    std::cout << "Pareto 前沿：" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        if (result[i].rank == 0) {
            std::cout << "T_in = " << result[i].vars[0] << "°C, t = " << result[i].vars[1] << "h, "
                      << "舒适度 = " << -result[i].objs[0] << ", 能耗 = " << result[i].objs[1] << std::endl;
        }
    }

    return 0;
}
