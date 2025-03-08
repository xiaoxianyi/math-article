#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// 自定义比较函数对象（仿函数）
struct CompareIndices {
    const std::vector<double>& diff;
    CompareIndices(const std::vector<double>& d) : diff(d) {}
    bool operator()(size_t i, size_t j) const {
        return std::abs(diff[i]) < std::abs(diff[j]);
    }
};

// 计算 Wilcoxon 符号秩检验
double wilcoxon_signed_rank_test(const std::vector<double>& x, const std::vector<double>& y, double& W_plus, double& W_minus) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("样本大小必须相同");
    }

    // 计算差值
    std::vector<double> diff(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        diff[i] = x[i] - y[i];
    }

    // 去除零差值
    diff.erase(std::remove(diff.begin(), diff.end(), 0.0), diff.end());

    // 如果所有差值都为零，返回 p 值为 1.0
    if (diff.empty()) {
        return 1.0;
    }

    // 对差值的绝对值进行排序并分配秩次
    std::vector<size_t> indices(diff.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // 使用自定义比较函数对象排序
    std::sort(indices.begin(), indices.end(), CompareIndices(diff));

    std::vector<double> ranks(diff.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        ranks[indices[i]] = i + 1;
    }

    // 计算正差值和负差值的秩和
    W_plus = 0.0;
    W_minus = 0.0;
    for (size_t i = 0; i < diff.size(); ++i) {
        if (diff[i] > 0) {
            W_plus += ranks[i];
        } else {
            W_minus += ranks[i];
        }
    }

    // 统计量 W
    double W = std::min(W_plus, W_minus);

    // 计算 p 值（近似正态分布）
    size_t n = diff.size();
    if (n < 10) {
        std::cerr << "警告：样本量较小，正态近似可能不准确。" << std::endl;
    }

    double mean = n * (n + 1) / 4.0;
    double std_dev = std::sqrt(n * (n + 1) * (2 * n + 1) / 24.0);
    double z = (W - mean) / std_dev;
    double p_value = 2 * (1 - 0.5 * (1 + erf(z / std::sqrt(2)))); // 使用 erf 而不是 std::erf

    return p_value;
}

int main() {
    // 示例数据
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
        // 计算 Wilcoxon 符号秩检验
        double p_value = wilcoxon_signed_rank_test(x, y, W_plus, W_minus);

        // 输出结果
        std::cout << "Wilcoxon 符号秩检验结果：" << std::endl;
        std::cout << "统计量 W: " << std::min(W_plus, W_minus) << std::endl;
        std::cout << "样本量 n: " << x.size() << std::endl;
        std::cout << "p 值: " << std::setprecision(4) << p_value << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }

    return 0;
}
