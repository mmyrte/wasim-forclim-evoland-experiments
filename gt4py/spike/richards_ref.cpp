// Richards proxy, C++ reference: adaptive per-cell sub-steps and iterations, as WaSiM loops.
// The algorithm is specified in README.md; proxy_case.py writes the case.
//
// build: c++ -O3 -fopenmp -o richards_ref richards_ref.cpp
// usage: richards_ref CASE_DIR OUT_DIR [repeat [steps]]   (threads: OMP_NUM_THREADS)
// writes OUT_DIR/{theta,drain,excess,eta,substeps,iters}.f64 and prints ms per step
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

enum { KS, THS, THR, ALPHA, INV_N, M, INV_M, DZ, ZC, NPAR };

static std::vector<double> read(const std::string &path, size_t n) {
    std::vector<double> v(n);
    FILE *f = std::fopen(path.c_str(), "rb");
    if (!f || std::fread(v.data(), sizeof(double), n, f) != n) { std::perror(path.c_str()); std::exit(1); }
    std::fclose(f);
    return v;
}

static void write(const std::string &path, const std::vector<double> &v) {
    FILE *f = std::fopen(path.c_str(), "wb");
    std::fwrite(v.data(), sizeof(double), v.size(), f);
    std::fclose(f);
}

struct Column { const double *p; int L; };

static inline double sat(const Column &c, int k, double th) {
    double s = (th - c.p[THR * c.L + k]) / (c.p[THS * c.L + k] - c.p[THR * c.L + k]);
    if (s < 1e-6) s = 1e-6;
    if (s > 1.0) s = 1.0;
    return s;
}
static inline double suction(const Column &c, int k, double th) {
    double s = sat(c, k, th);
    return std::pow(std::pow(s, -c.p[INV_M * c.L + k]) - 1.0, c.p[INV_N * c.L + k]) / c.p[ALPHA * c.L + k];
}
static inline double krel(const Column &c, int k, double th) {
    double s = sat(c, k, th);
    double t = 1.0 - std::pow(1.0 - std::pow(s, c.p[INV_M * c.L + k]), c.p[M * c.L + k]);
    return std::sqrt(s) * t * t;
}

int main(int argc, char **argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: richards_ref CASE_DIR OUT_DIR [repeat [steps]]\n"); return 1; }
    std::string dir = argv[1], out = argv[2];
    int repeat = argc > 3 ? std::atoi(argv[3]) : 1;
    int N, L, T, max_iter, root_levels;
    double DT, max_err, min_step;
    FILE *f = std::fopen((dir + "/meta.txt").c_str(), "r");
    if (!f || std::fscanf(f, "%d %d %d %lf %lf %d %lf %d", &N, &L, &T, &DT, &max_err, &max_iter, &min_step,
                          &root_levels) != 8) { std::fprintf(stderr, "bad meta.txt\n"); return 1; }
    std::fclose(f);
    int T_all = T;
    if (argc > 4 && std::atoi(argv[4]) > 0 && std::atoi(argv[4]) < T) T = std::atoi(argv[4]);
    auto params = read(dir + "/params.f64", (size_t)N * NPAR * L);
    auto theta0 = read(dir + "/theta.f64", (size_t)N * L);
    auto rain = read(dir + "/rain.f64", (size_t)T_all * N);
    auto et = read(dir + "/et.f64", (size_t)T_all * N);
    std::vector<double> theta, drain((size_t)T * N), excess((size_t)T * N), eta((size_t)T * N),
        substeps((size_t)T * N), iters((size_t)T * N);
    double best = 1e30;
    for (int r = 0; r < repeat; r++) {
        theta = theta0;
        auto t0 = std::chrono::steady_clock::now();
        for (int t = 0; t < T; t++) {
#pragma omp parallel for schedule(dynamic, 64)
            for (int c = 0; c < N; c++) {
                Column col{&params[(size_t)c * NPAR * L], L};
                const double *dz = col.p + DZ * L, *zc = col.p + ZC * L, *thr = col.p + THR * L,
                             *ths = col.p + THS * L, *ks = col.p + KS * L;
                double *th = &theta[(size_t)c * L];
                double P = rain[(size_t)t * N + c], E = et[(size_t)t * N + c];
                double old[64], avg[64], h1[64], hh[64], K[64], q[65];
                // Courant sub-steps (WaSiM: unsatzon.cc, run, methode > 2)
                double dtmin = DT;
                for (int k = 0; k < L - 1; k++) {
                    double kr = P > 0 ? 1.0 : krel(col, k, th[k]);
                    if (kr > 0) {
                        double d = dz[k] / (kr * ks[k]);
                        if (d < dtmin) dtmin = d;
                    }
                }
                if (dtmin < min_step) dtmin = min_step;
                int n = (int)std::ceil(DT / dtmin);
                double dt = DT / n, Psub = P / n, Esub = E / n;
                double sum_drain = 0, sum_excess = 0, sum_eta = 0;
                int sum_iter = 0;
                for (int s = 0; s < n; s++) {
                    for (int k = 0; k < root_levels; k++) {
                        double take = Esub / root_levels / 1000.0 / dz[k];
                        double avail = th[k] - thr[k];
                        if (avail < 0) avail = 0;
                        if (take > avail) take = avail;
                        th[k] -= take;
                        sum_eta += take * dz[k] * 1000.0;
                    }
                    for (int k = 0; k < L; k++) {
                        old[k] = avg[k] = th[k];
                        h1[k] = suction(col, k, old[k]);
                    }
                    int z = 0;
                    long weight = 0;
                    double maxe = max_err, err;
                    do {
                        for (int k = 0; k < L; k++) {
                            double h2 = suction(col, k, avg[k]);
                            hh[k] = zc[k] - 0.5 * (h1[k] + h2);
                            K[k] = ks[k] * krel(col, k, avg[k]);
                        }
                        q[0] = Psub / (1000.0 * dt);
                        for (int k = 0; k < L - 1; k++)
                            q[k + 1] = std::sqrt(K[k] * K[k + 1]) * (hh[k] - hh[k + 1]) / (0.5 * (dz[k] + dz[k + 1]));
                        q[L] = K[L - 1];
                        for (int k = 0; k < L; k++) {  // no level below its residual water content
                            double allowed = q[k] + (old[k] - thr[k]) * dz[k] / dt;
                            if (q[k + 1] > allowed) q[k + 1] = allowed;
                        }
                        for (int k = L - 1; k >= 0; k--) {  // no level above saturation
                            double w = old[k] + (q[k] - q[k + 1]) * dt / dz[k];
                            if (w > ths[k]) q[k] = q[k + 1] + (ths[k] - old[k]) * dz[k] / dt;
                        }
                        z++;
                        if (z == 21) maxe = max_err * 5.0;
                        if (z > 3) weight = z - 3;
                        err = 0;
                        for (int k = 0; k < L; k++) {  // WaSiM new_thetas
                            double last = th[k];
                            if (z <= 3) avg[k] = th[k];
                            th[k] = old[k] + (q[k] - q[k + 1]) * dt / dz[k];
                            if (th[k] < 0.0) th[k] = 0.0;
                            else if (th[k] > 1.0) th[k] = 1.0;
                            double e = std::fabs(last - th[k]);
                            if (e > err) err = e;
                            if (z > 3) {
                                avg[k] *= weight;
                                avg[k] += th[k];
                                avg[k] /= (weight + 1);
                            } else {
                                avg[k] += th[k];
                                avg[k] /= 2.0;
                            }
                        }
                    } while (err > maxe && z < max_iter);
                    sum_excess += Psub - q[0] * dt * 1000.0;
                    sum_drain += q[L] * dt * 1000.0;
                    sum_iter += z;
                }
                size_t i = (size_t)t * N + c;
                drain[i] = sum_drain; excess[i] = sum_excess; eta[i] = sum_eta;
                substeps[i] = n; iters[i] = sum_iter;
            }
        }
        double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / T;
        if (ms < best) best = ms;
    }
    write(out + "/theta.f64", theta);
    write(out + "/drain.f64", drain);
    write(out + "/excess.f64", excess);
    write(out + "/eta.f64", eta);
    write(out + "/substeps.f64", substeps);
    write(out + "/iters.f64", iters);
    double it = 0, sub = 0;
    for (size_t i = 0; i < iters.size(); i++) { it += iters[i]; sub += substeps[i]; }
    std::printf("%d cells, %d levels, %d steps: %.3f ms per step (best of %d); "
                "mean sub-steps %.3f, mean iterations per cell-step %.2f\n",
                N, L, T, best, repeat, sub / iters.size(), it / iters.size());
    return 0;
}
