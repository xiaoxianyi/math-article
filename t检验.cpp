#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// �Զ���ȽϺ������󣨷º�����
struct CompareIndices {
    const std::vector<double>& diff;
    CompareIndices(const std::vector<double>& d) : diff(d) {}
    bool operator()(size_t i, size_t j) const {
        return std::abs(diff[i]) < std::abs(diff[j]);
    }
};

// ���� Wilcoxon �����ȼ���
double wilcoxon_signed_rank_test(const std::vector<double>& x, const std::vector<double>& y, double& W_plus, double& W_minus) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("������С������ͬ");
    }

    // �����ֵ
    std::vector<double> diff(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        diff[i] = x[i] - y[i];
    }

    // ȥ�����ֵ
    diff.erase(std::remove(diff.begin(), diff.end(), 0.0), diff.end());

    // ������в�ֵ��Ϊ�㣬���� p ֵΪ 1.0
    if (diff.empty()) {
        return 1.0;
    }

    // �Բ�ֵ�ľ���ֵ�������򲢷����ȴ�
    std::vector<size_t> indices(diff.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // ʹ���Զ���ȽϺ�����������
    std::sort(indices.begin(), indices.end(), CompareIndices(diff));

    std::vector<double> ranks(diff.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        ranks[indices[i]] = i + 1;
    }

    // ��������ֵ�͸���ֵ���Ⱥ�
    W_plus = 0.0;
    W_minus = 0.0;
    for (size_t i = 0; i < diff.size(); ++i) {
        if (diff[i] > 0) {
            W_plus += ranks[i];
        } else {
            W_minus += ranks[i];
        }
    }

    // ͳ���� W
    double W = std::min(W_plus, W_minus);

    // ���� p ֵ��������̬�ֲ���
    size_t n = diff.size();
    if (n < 10) {
        std::cerr << "���棺��������С����̬���ƿ��ܲ�׼ȷ��" << std::endl;
    }

    double mean = n * (n + 1) / 4.0;
    double std_dev = std::sqrt(n * (n + 1) * (2 * n + 1) / 24.0);
    double z = (W - mean) / std_dev;
    double p_value = 2 * (1 - 0.5 * (1 + erf(z / std::sqrt(2)))); // ʹ�� erf ������ std::erf

    return p_value;
}

int main() {
    // ʾ������
    std::vector<double> x;
    x.push_back(10.5);
    x.push_back(12.3);
    x.push_back(9.8);
    x.push_back(11.2);
    x.push_back(13.1);

    std::vector<double> y;
    y.push_back(10.0);
    y.push_back(11.8);
    y.push_back(9.5);
    y.push_back(10.9);
    y.push_back(12.5);

    try {
        double W_plus, W_minus;
        // ���� Wilcoxon �����ȼ���
        double p_value = wilcoxon_signed_rank_test(x, y, W_plus, W_minus);

        // ������
        std::cout << "Wilcoxon �����ȼ�������" << std::endl;
        std::cout << "ͳ���� W: " << std::min(W_plus, W_minus) << std::endl;
        std::cout << "������ n: " << x.size() << std::endl;
        std::cout << "p ֵ: " << std::setprecision(4) << p_value << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "����: " << e.what() << std::endl;
    }

    return 0;
}
