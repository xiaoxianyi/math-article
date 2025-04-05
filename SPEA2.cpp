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
    double T_out = 29.5;
    double k = 0.35;

    // ����ԭʼ���ʶȹ�ʽ
    double comfort = 1.0 / (std::abs(T_in - 25) + 1e-6);
    comfort = clamp(comfort, 0.0, 1.0);

    double energy = k * std::abs(T_in - T_out) * t;

    ind.objs[0] = -comfort; // ��С�� -comfort������� comfort��
    ind.objs[1] = energy;
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

// �ȽϺ����������򣨰�Ŀ�꺯��ֵ����
struct ObjComparator {
    size_t m;
    ObjComparator(size_t m) : m(m) {}
    bool operator()(const Individual& a, const Individual& b) const {
        return a.objs[m] < b.objs[m];
    }
};

// ����ӵ����
void calculate_crowding(std::vector<Individual>& front) {
    if (front.size() <= 2) {
        for (size_t i = 0; i < front.size(); ++i) {
            front[i].crowding = std::numeric_limits<double>::max();
        }
        return;
    }

    for (size_t i = 0; i < front.size(); ++i) {
        front[i].crowding = 0.0;
    }

    for (size_t m = 0; m < front[0].objs.size(); ++m) {
        std::sort(front.begin(), front.end(), ObjComparator(m));
        
        double min_obj = front.front().objs[m];
        double max_obj = front.back().objs[m];
        double range = max_obj - min_obj;
        
        if (range < 1e-6) continue;
        
        front[0].crowding = front.back().crowding = std::numeric_limits<double>::max();
        
        for (size_t i = 1; i < front.size() - 1; ++i) {
            front[i].crowding += (front[i+1].objs[m] - front[i-1].objs[m]) / range;
        }
    }
}

// �ȽϺ����������򣨰���Ӧ��ֵ����
bool less_fitness(const Individual& a, const Individual& b) {
    if (a.fitness != b.fitness) {
        return a.fitness < b.fitness;
    }
    return a.crowding > b.crowding; // ӵ���ȴ������
}

// ����ѡ��
void environmental_selection(std::vector<Individual>& pop, std::vector<Individual>& archive, size_t archive_size) {
    std::vector<Individual> combined_pop = pop;
    combined_pop.insert(combined_pop.end(), archive.begin(), archive.end());

    calculate_strength(combined_pop);
    calculate_fitness(combined_pop);
    calculate_crowding(combined_pop);

    // ����Ӧ�Ⱥ�ӵ��������
    std::sort(combined_pop.begin(), combined_pop.end(), less_fitness);

    // ѡ����Ӧ����õĸ���
    archive.clear();
    size_t count = 0;
    for (size_t i = 0; i < combined_pop.size() && count < archive_size; ++i) {
        if (combined_pop[i].fitness < 1.0) { // ֻѡ���֧���
            archive.push_back(combined_pop[i]);
            count++;
        }
    }

    // �����֧��ⲻ�㣬����������
    if (archive.size() < archive_size) {
        for (size_t i = 0; i < combined_pop.size() && archive.size() < archive_size; ++i) {
            bool exists = false;
            for (size_t j = 0; j < archive.size(); ++j) {
                if (combined_pop[i].vars[0] == archive[j].vars[0] && 
                    combined_pop[i].vars[1] == archive[j].vars[1]) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                archive.push_back(combined_pop[i]);
            }
        }
    }
}

// �������
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

    // ��ӱ�����ΧԼ��
    c1.vars[0] = clamp(c1.vars[0], 20.0, 30.0);
    c1.vars[1] = clamp(c1.vars[1], 14.0, 15.0);
    c2.vars[0] = clamp(c2.vars[0], 20.0, 30.0);
    c2.vars[1] = clamp(c2.vars[1], 14.0, 15.0);
}

// �������
void mutate(Individual& ind) {
    double rate = 0.2; // ��߱�����
    for (size_t i = 0; i < ind.vars.size(); ++i) {
        if (rand() / (double)RAND_MAX < rate) {
            ind.vars[i] += (rand() / (double)RAND_MAX - 0.5) * 1.0; // ����������
        }
    }

    // �ϸ�Լ����Χ
    ind.vars[0] = clamp(ind.vars[0], 20.0, 30.0);
    ind.vars[1] = clamp(ind.vars[1], 14.0, 15.0);
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
    }

    for (size_t gen = 0; gen < max_gen; ++gen) {
        std::cout << "Generation " << gen + 1 << "/" << max_gen << std::endl;

        // ������Ⱥ
        for (size_t i = 0; i < pop.size(); ++i) {
            evalObj(pop[i]);
        }

        // ����ѡ��
        environmental_selection(pop, archive, archive_size);

        // ��������Ⱥ
        std::vector<Individual> offspring;
        while (offspring.size() < pop_size) {
            size_t i1 = rand() % archive.size();
            size_t i2 = rand() % archive.size();
            if (i1 == i2) continue;
            
            Individual c1, c2;
            crossover(archive[i1], archive[i2], c1, c2);
            mutate(c1);
            mutate(c2);
            evalObj(c1);
            evalObj(c2);
            offspring.push_back(c1);
            if (offspring.size() < pop_size) {
                offspring.push_back(c2);
            }
        }

        pop = offspring;
    }

    // ���ջ���ѡ��
    environmental_selection(pop, archive, archive_size);
    return archive;
}

int main() {
    srand(time(NULL));  // ��ʼ���������
    
    // ��������
    size_t pop_size = 100;     // ������Ⱥ��С
    size_t archive_size = 100;  // ����浵��С
    size_t max_gen = 50;       // ���ӵ�������

    std::vector<Individual> result = spea2(pop_size, archive_size, max_gen);

    // ���Paretoǰ��
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
