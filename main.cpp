#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

using namespace std;

/***
 * j бля всегда столбцы !!!
 * @param A
 * @param b
 * @return
 */
vector<double> operator*(vector<vector<double>> A, vector<double> b) {
    vector<double> c(A[0].size(), 0);
    for (int i = 0; i < A[0].size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            c[i] += A[j][i] * b[j];
        }
    }
    return c;
}

vector<double> operator*(double c, vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] = c * b[i];
    }
    return b;
}

vector<double> operator+(vector<double> c, vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] = c[i] + b[i];
    }
    return b;
}

double max_difference(vector<vector<double>> A, vector<vector<double>> B) {
    int max = 0;
    for (int j = 0; j < A.size(); j++) {
        for (int i = 0; i < A[j].size(); i++) {
            if (max < abs(A[j][i] - B[j][i])) {
                max = abs(A[j][i] - B[j][i]);
            }
        }
    }
    return max;
}

vector<vector<double>> K(double time, double h, vector<double> c, vector<double> y, vector<vector<double>> A,
                         function<vector<double>(double, vector<double>)> f) {
    vector<vector<double>> K0(c.size());
    for (int i = 0; i < K0.size(); i++) {
        K0[i].resize(y.size());
    }
    vector<vector<double>> K(c.size());
    for (int i = 0; i < K.size(); i++) {
        K[i].resize(y.size());
    }
    double error = 1e10;
    while (error > 1e-7) {
        for (int j = 0; j < K.size(); j++) {
            K[j] = f(time + h * c[j], y + h * (K * A[j]));
        }
        error = max_difference(K, K0);
        K0 = K;

    }
    return K;
}

vector<double> step(double time, double h, vector<double> c, vector<double> y, vector<vector<double>> A,
                    function<vector<double>(double, vector<double>)> f, vector<double> b) {
    return y + h * (K(time, h, c, y, A, f) * b);
}

vector<pair<double, vector<double>>>
solve(double start, double end, double h, vector<double> y0, vector<double> c, vector<vector<double>> A,
      vector<double> b,
      function<vector<double>(double, vector<double>)> f) {
    double time = start;
    vector<pair<double, vector<double>>> result;
    vector<double> y = y0;
    result.emplace_back(time, y);
    while (time < end) {
        y = step(time, h, c, y, A, f, b);
        time += h;
        result.emplace_back(time, y);
    }
    return result;
}

vector<double> f(double time, vector<double> y) {
    return {y[1], time * sqrt(y[0])};
}


double starting_condition(vector<double> c, vector<vector<double>> A, vector<double> b, double y_last) {
    double first = 1;
    double last = 2;
    auto result_a = solve(0, 1, 0.01, {0, first}, c, A, b, f);
    auto result_b = solve(0, 1, 0.01, {0, last}, c, A, b, f);
    while (abs(last-first) > 1e-6) {
        double median = (first + last) / 2;
        auto result_c = solve(0, 1, 0.01, {0, median}, c, A, b, f);
        if (result_c.back().second[0] <= y_last) {
            first = median;
            auto result_a = solve(0, 1, 0.01, {0, first}, c, A, b, f);
        }
        if (result_c.back().second[0] > y_last) {
            last = median;
            auto result_b = solve(0, 1, 0.01, {0, last}, c, A, b, f);
        }

    }
    return first;
}

int main() {
    vector<vector<double>> A = {{0,   0,   0, 0},
                                {0.5, 0,   0, 0},
                                {0,   0.5, 0, 0},
                                {0,   0,   1, 0}};
    vector<double> b = {1. / 6, 1. / 3, 1. / 3, 1. / 6};
    std::vector<double> c = {0, 0.5, 0.5, 1};
    auto result = solve(0, 1, 0.01, {0, starting_condition(c, A, b,2)}, c, A, b, f);
    cout << result.back().first << ' ' << result.back().second[0] << endl;
    return 0;
}