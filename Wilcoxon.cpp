#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cstdlib>

// ������
struct Individual {
    std::vector<double> vars;  // T_in, t
    std::vector<double> objs;  // ���ʶ�, �ܺ�
    double fitness;            // ��Ӧ��ֵ
    int strength;              // ǿ��ֵ
    double crowding;           // ӵ����

    Individual() : fitness(0.0), strength(0), crowding(0.0) {
        vars.resize(2);
        objs.resize(2);
    }
};

// �Զ���clamp����
double clamp(double value, double min, double max) {
    return (value < min) ? min : (value > max) ? max : value;
}

// ����Ŀ�꺯��
void evalObj(Individual& ind) {
    double T_in = ind.vars[0];
    double t = ind.vars[1];

    // Լ������
    if (T_in < 20 || T_in > 30 || t < 14 || t > 15) {
        ind.objs[0] = std::numeric_limits<double>::max();
        ind.objs[1] = std::numeric_limits<double>::max();
        return;
    }

    // ���ʶȹ�ʽ
    double comfort = 1.0 / (std::abs(T_in - 25) + 1e-6); // +1e-6���������
    comfort = clamp(comfort, 0.0, 1.0); // ʹ���Զ���clamp

    // �ܺĹ�ʽ
    double T_out = 35;
    double k = 0.1;
    double energy = k * std::abs(T_in - T_out) * t;

    ind.objs[0] = -comfort; // ��С�� -comfort������� comfort��
    ind.objs[1] = energy;

    // �������
    std::cout << "T_in = " << T_in << ", t = " << t << ", comfort = " << comfort << ", energy = " << energy << std::endl;
}

// �жϸ��� a �Ƿ�֧����� b
bool dominates(const Individual& a, const Individual& b) {
    bool not_worse = true;
    bool better = false;

    for (size_t i = 0; i < a.objs.size(); ++i) {
        if (a.objs[i] > b.objs[i]) {
            not_worse = false;
            break;
        }
        if (a.objs[i] < b.objs[i]) {
            better = true;
        }
    }

    return not_worse && better;
}

// ����ǿ��ֵ
void calculate_strength(std::vector<Individual>& pop) {
    for (size_t i = 0; i < pop.size(); ++i) {
        pop[i].strength = 0;
        for (size_t j = 0; j < pop.size(); ++j) {
            if (i == j) continue;
            if (dominates(pop[i], pop[j])) {
                pop[i].strength++;
            }
        }
    }
}

// ������Ӧ��ֵ
void calculate_fitness(std::vector<Individual>& pop) {
    for (size_t i = 0; i < pop.size(); ++i) {
        double raw_fitness = 0.0;
        for (size_t j = 0; j < pop.size(); ++j) {
            if (i == j) continue;
            if (dominates(pop[j], pop[i])) {
                raw_fitness += pop[j].strength;
            }
        }
        pop[i].fitness = raw_fitness;
    }
}

// ����ӵ����
void calculate_crowding(std::vector<Individual>& front) {
    for (size_t i = 0; i < front.size(); ++i) {
        front[i].crowding = 0.0;
    }

    for (size_t m = 0; m < front[0].objs.size(); ++m) {
        // ����
        for (size_t i = 0; i < front.size(); ++i) {
            for (size_t j = i + 1; j < front.size(); ++j) {
                if (front[i].objs[m] > front[j].objs[m]) {
                    std::swap(front[i], front[j]);
                }
            }
        }

        front[0].crowding = front.back().crowding = std::numeric_limits<double>::max();

        for (size_t i = 1; i < front.size() - 1; ++i) {
            front[i].crowding += (front[i + 1].objs[m] - front[i - 1].objs[m]);
        }
    }
}

// �ȽϺ�����������
bool less_fitness(const Individual& a, const Individual& b) {
    return a.fitness < b.fitness;
}

// ����ѡ��
void environmental_selection(std::vector<Individual>& pop, std::vector<Individual>& archive, size_t archive_size) {
    std::vector<Individual> combined_pop = pop;
    combined_pop.insert(combined_pop.end(), archive.begin(), archive.end());

    calculate_strength(combined_pop);
    calculate_fitness(combined_pop);

    // ���������Ӧ��ֵ
    std::cout << "Fitness values before sorting:" << std::endl;
    for (size_t i = 0; i < combined_pop.size(); ++i) {
        std::cout << "Individual " << i << ": fitness = " << combined_pop[i].fitness << std::endl;
    }

    // ѡ����Ӧ��ֵ�ϵ͵ĸ���
    std::sort(combined_pop.begin(), combined_pop.end(), less_fitness);

    archive.clear();
    for (size_t i = 0; i < archive_size && i < combined_pop.size(); ++i) {
        archive.push_back(combined_pop[i]);
    }

    // ���������Ӧ��ֵ
    std::cout << "Fitness values after sorting:" << std::endl;
    for (size_t i = 0; i < archive.size(); ++i) {
        std::cout << "Archive " << i << ": fitness = " << archive[i].fitness << std::endl;
    }

    // ����ӵ������ά��������
    calculate_crowding(archive);
}

// �������
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

    // ��ӱ�����ΧԼ��
    c1.vars[0] = clamp(c1.vars[0], 20.0, 30.0);
    c1.vars[1] = clamp(c1.vars[1], 14.0, 15.0);
    c2.vars[0] = clamp(c2.vars[0], 20.0, 30.0);
    c2.vars[1] = clamp(c2.vars[1], 14.0, 15.0);

    // �������
    std::cout << "Crossover: c1 = (" << c1.vars[0] << ", " << c1.vars[1] << "), c2 = (" << c2.vars[0] << ", " << c2.vars[1] << ")" << std::endl;
}

// �������
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

    // �������
    std::cout << "Mutate: vars = (" << ind.vars[0] << ", " << ind.vars[1] << ")" << std::endl;
}

// SPEA2 ������
std::vector<Individual> spea2(size_t pop_size, size_t archive_size, size_t max_gen) {
    std::vector<Individual> pop(pop_size);
    std::vector<Individual> archive(archive_size);

    // ��ʼ����Ⱥ
    for (size_t i = 0; i < pop_size; ++i) {
        pop[i].vars[0] = 20 + (rand() / (double)RAND_MAX) * 10; // T_in �� [20,30]
        pop[i].vars[1] = 14 + (rand() / (double)RAND_MAX) * 1;  // t �� [14,15]
        evalObj(pop[i]);

        // �������
        std::cout << "Initial T_in = " << pop[i].vars[0] << ", t = " << pop[i].vars[1] << std::endl;
    }

    for (size_t gen = 0; gen < max_gen; ++gen) {
        // ������Ⱥ
        for (size_t i = 0; i < pop.size(); ++i) {
            evalObj(pop[i]);
        }

        // ����ѡ��
        environmental_selection(pop, archive, archive_size);

        // ��������Ⱥ
        std::vector<Individual> offspring;
        for (size_t i = 0; i < pop_size; i += 2) {
            Individual c1, c2;
            crossover(archive[i], archive[i + 1], c1, c2);
            mutate(c1);
            mutate(c2);
            evalObj(c1);
            evalObj(c2);
            offspring.push_back(c1);
            offspring.push_back(c2);
        }

        pop = offspring;
    }

    return archive;
}

int main() {
    srand(time(NULL));  // ��ʼ���������
    size_t pop_size = 100;    // ��Ⱥ��С
    size_t archive_size = 100; // �浵��С
    size_t max_gen = 100;     // ����������

    std::vector<Individual> result = spea2(pop_size, archive_size, max_gen);

    std::cout << "Pareto ǰ�أ�" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        std::cout << "T_in = " << result[i].vars[0] << "��C, t = " << result[i].vars[1] << "h, "
                  << "���ʶ� = " << -result[i].objs[0] << ", �ܺ� = " << result[i].objs[1] << std::endl;
    }

    return 0;
}
