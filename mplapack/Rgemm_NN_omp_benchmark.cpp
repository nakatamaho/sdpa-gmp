#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <gmpxx_mkII.h>
#include <mpblas_gmp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void Rgemm_NN_omp(mplapackint m, mplapackint n, mplapackint k, mpf_class alpha, mpf_class *A, mplapackint lda, mpf_class *B, mplapackint ldb, mpf_class beta, mpf_class *C, mplapackint ldc);

void generate_random_matrix(mplapackint rows, mplapackint cols, mpf_class *matrix, gmp_randclass &rng, unsigned int prec) {
    for (mplapackint i = 0; i < rows; ++i) {
        for (mplapackint j = 0; j < cols; ++j) {
            matrix[i + j * rows] = rng.get_f(prec);
        }
    }
}

template <typename Func> double benchmark(Func func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

std::pair<double, double> calculate_mean_and_variance(const std::vector<double> &values) {
    double mean = 0.0;
    double variance = 0.0;

    for (double value : values) {
        mean += value;
    }
    mean /= values.size();

    for (double value : values) {
        variance += (value - mean) * (value - mean);
    }
    variance /= values.size();

    return {mean, variance};
}

int main() {
#ifdef _OPENMP
    std::cout << "OpenMP is enabled.\n";
    std::cout << "Number of threads: " << omp_get_max_threads() << "\n";
#else
    std::cout << "OpenMP is not enabled.\n";
#endif

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    unsigned int prec = 512;

    // ベンチマーク用の行列サイズ (1から1024まで)
    //    std::vector<mplapackint> sizes = {8, 10, 15, 16, 32,48, 64, 128, 256, 512, 1024};
    std::vector<mplapackint> sizes = {128, 256, 512, 768};
    // スレッド数を指定
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 32};

    const int num_trials = 10; // ベンチマーク回数

    for (auto m : sizes) {
        for (auto n : sizes) {
            for (auto k : sizes) {
                // FLOP 数の計算: m * n * (2 * k + 1)
                double flop_count = static_cast<double>(m) * n * (2.0 * k + 1);

                // 行列とスカラーを初期化
                std::vector<mpf_class> A(m * k);
                std::vector<mpf_class> B(k * n);
                std::vector<mpf_class> C(m * n);

                // スカラー alpha, beta を初期化
                mpf_class alpha = rng.get_f(prec);
                mpf_class beta = rng.get_f(prec);

                // ランダム行列を生成
                generate_random_matrix(m, k, A.data(), rng, prec);
                generate_random_matrix(k, n, B.data(), rng, prec);
                generate_random_matrix(m, n, C.data(), rng, prec);

                std::cout << "Benchmarking m=" << m << ", n=" << n << ", k=" << k << ":\n";

                for (auto threads : thread_counts) {
#ifdef _OPENMP
                    omp_set_num_threads(threads);
#endif
                    std::vector<double> flops_results;

                    for (int trial = 0; trial < num_trials; ++trial) {
                        double elapsed = benchmark([&]() { Rgemm_NN_omp(m, n, k, alpha, A.data(), m, B.data(), k, beta, C.data(), m); });

                        double flops = flop_count / elapsed / 1e6; // MFLOPS
                        flops_results.push_back(flops);
                    }

                    // 平均と分散を計算
                    auto [mean_flops, variance_flops] = calculate_mean_and_variance(flops_results);

                    std::cout << "Threads: " << threads << ", Mean FLOPS: " << mean_flops << " MFLOPS"
                              << ", Variance: " << variance_flops << "\n";
                }
                std::cout << "---------------------------------\n";
            }
        }
    }

    return 0;
}
