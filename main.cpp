#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

struct Circle {
    double x;
    double y;
    double r;
};

bool insideCircle(double px, double py, const Circle& circle) {
    return (px - circle.x) * (px - circle.x) + (py - circle.y) * (py - circle.y) <= circle.r * circle.r;
}

bool insideIntersection(double px, double py, const Circle& c1, const Circle& c2, const Circle& c3) {
    return insideCircle(px, py, c1) && insideCircle(px, py, c2) && insideCircle(px, py, c3);
}

double monteCarlo(const Circle& c1, const Circle& c2, const Circle& c3, int n, double& error, double scale = 1.0) {
    double minX = std::min({c1.x - c1.r, c2.x - c2.r, c3.x - c3.r}) * scale;
    double maxX = std::max({c1.x + c1.r, c2.x + c2.r, c3.x + c3.r}) * scale;
    double minY = std::min({c1.y - c1.r, c2.y - c2.r, c3.y - c3.r}) * scale;
    double maxY = std::max({c1.y + c1.r, c2.y + c2.r, c3.y + c3.r}) * scale;

    double rectangle = (maxX - minX) * (maxY - minY);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distX(minX, maxX);
    std::uniform_real_distribution<> distY(minY, maxY);

    int count = 0;
    for (int i = 0; i < n; ++i) {
        double px = distX(gen);
        double py = distY(gen);

        if (insideIntersection(px / scale, py / scale, c1, c2, c3)) {
            ++count;
        }
    }

    double estimated_res = static_cast<double>(count) / n * rectangle;
    error = std::abs(estimated_res - rectangle) / rectangle * 100.0;
    return estimated_res;
}

int main() {
    Circle c1 = {1.0, 1.0, 1.0};
    Circle c2 = {1.5, 2.0, 1.0};
    Circle c3 = {2.0, 1.5, 1.0};

    std::vector<int> N_values;
    for (int N = 100; N <= 100000; N += 500) {
        N_values.push_back(N);
    }

    std::vector<double> scales = {1.0, 1.5, 2.0};

    std::ofstream outFile("monte_carlo_results.txt");
    outFile << std::fixed << std::setprecision(6);
    outFile << "Scale\tN\tEstimated Area\tRelative Error\n";

    for (double scale : scales) {
        for (int N : N_values) {
            double error = 0.0;
            double estimated_res = monteCarlo(c1, c2, c3, N, error, scale);
            outFile << scale << "\t" << N << "\t" << estimated_res << "\t" << error << "\n";
        }
    }

    outFile.close();
    std::cout << "Done" << "\n";

    return 0;
}