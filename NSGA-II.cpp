#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cstdlib>

struct Individual {
    std::vector<double> vars;  // T_in, t
    std::vector<double> objs;  // ���ʶ�, �ܺ�
    int rank;
    double crowding;
};

// �Զ���clamp����������C++11��
double clamp(double value, double min, double max) {
    return (value < min) ? min : (value > max) ? max : value;
}

void evalObj(Individual& ind) {
    double T_in = ind.vars[0];
    double t = ind.vars[1];
    double T_out = 29.5;
    double k = 0.35;

    // �޸����ʶȹ�ʽʹ���ֵ�ӽ� 1
    double comfort = 1.0 / (std::abs(T_in - 25) + 1e-6); // +1e-6���������
    comfort = clamp(comfort, 0.0, 1.0); // ʹ���Զ���clamp

    double energy = k * std::abs(T_in - T_out) * t;

    ind.objs.clear();
    ind.objs.push_back(-comfort); // ����Ҫ��С�� -comfort������� comfort��
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

    // ��ӱ�����ΧԼ����ʹ���Զ���clamp��
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

    // �ϸ�Լ����Χ
    ind.vars[0] = clamp(ind.vars[0], 20.0, 30.0);
    ind.vars[1] = clamp(ind.vars[1], 14.0, 15.0);
}

std::vector<Individual> nsga2(size_t popSize, size_t maxGen) {
    std::vector<Individual> pop(popSize); // �޸�pop����

    // ��ʼ����Ⱥ
    for (size_t i = 0; i < popSize; ++i) {
        pop[i].vars.push_back(20 + (rand() / (double)RAND_MAX) * 10); // T_in �� [20,30]
        pop[i].vars.push_back(14 + (rand() / (double)RAND_MAX) * 1);  // t �� [14,15]
        evalObj(pop[i]);
    }

    for (size_t gen = 0; gen < maxGen; ++gen) {
        // ... ԭ��nsga2�߼����ֲ��� ...
    }

    return pop;
}

int main() {
    srand(time(NULL));  
    size_t popSize = 1000;
    size_t maxGen = 100;

    std::vector<Individual> result = nsga2(popSize, maxGen);

std::cout << "\nPareto ǰ�أ��� " << result.size() << " ���⣩��" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        double comfort = -result[i].objs[0];
        std::cout << "�� " << i+1 << ": T_in = " << result[i].vars[0] 
                  << "��C, t = " << result[i].vars[1] << "h, "
                  << "���ʶ� = " << comfort << ", �ܺ� = " << result[i].objs[1] << std::endl;
    }

    // ������ʶȽӽ�1�Ľ�
    std::cout << "\n���ʶȽӽ�1�Ľ⣺" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        double comfort = -result[i].objs[0];
        if (comfort > 0.9) {
            std::cout << "�� " << i+1 << ": T_in = " << result[i].vars[0] 
                      << "��C, t = " << result[i].vars[1] << "h, "
                      << "���ʶ� = " << comfort << ", �ܺ� = " << result[i].objs[1] << std::endl;
        }
    }
    return 0;
}
