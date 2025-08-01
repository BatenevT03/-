#include "CGM.h"
#include <vector>
#include <cmath>
#include <iostream>

void CGM::init(const std::vector<int>& gi_s, const std::vector<int>& gj_s, const std::vector<double>& di_s,
    const std::vector<double>& gg_s, const std::vector<double>& rp_s, int n_s) {
    
 

    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    rp = rp_s;
    n = n_s;

    int m = gi[n];
    r.resize(n, 0.0);
    x0.resize(n, 0.0);
    z.resize(n, 0.0);
    p.resize(n, 0.0);
    s.resize(n, 0.0);

    L_di.resize(n, 0.0);
    L_gg.resize(m, 0.0);

    for (int i = 0; i < n; i++) {
        L_di[i] = di[i];
    }

    for (int i = 0; i < m; i++) {
        L_gg[i] = gg[i];
    }
}

void CGM::make_LLT_decomposition() {
    for (int k = 0; k < n; k++) {
        double sum_d = 0.0;
        int i_s = gi[k], i_e = gi[k + 1];

        for (int i = i_s; i < i_e; i++) {
            double sum_l = 0.0;
            int j_s = gi[gj[i]], j_e = gi[gj[i] + 1];

            for (int m = i_s; m < i; m++) {
                for (int j = j_s; j < j_e; j++) {
                    if (gj[m] == gj[j]) {
                        sum_l += L_gg[m] * L_gg[j];
                        j_s++;
                    }
                }
            }
            L_gg[i] = (gg[i] - sum_l) / L_di[gj[i]];
            sum_d += L_gg[i] * L_gg[i];
        }
        L_di[k] = std::sqrt(di[k] - sum_d);
    }
}

double CGM::dot_prod(const std::vector<double>& a, const std::vector<double>& b) {
    double d_p = 0.0;
    for (int i = 0; i < n; i++) {
        d_p += a[i] * b[i];
    }
    return d_p;
}

void CGM::mul_matrix(const std::vector<double>& f, std::vector<double>& x) {
    for (int i = 0; i < n; i++) {
        x[i] = di[i] * f[i];
        for (int k = gi[i]; k < gi[i + 1]; k++) {
            int j = gj[k];
            x[i] += gg[k] * f[j];
            x[j] += gg[k] * f[i];
        }
    }
}

void CGM::solve_L(const std::vector<double>& f, std::vector<double>& x) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = gi[i]; k < gi[i + 1]; k++) {
            sum += L_gg[k] * x[gj[k]];
        }
        x[i] = (f[i] - sum) / L_di[i];
    }
}

void CGM::solve_LT(std::vector<double>& f, std::vector<double>& x) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = f[i] / L_di[i];
        for (int k = gi[i]; k < gi[i + 1]; k++) {
            f[gj[k]] -= L_gg[k] * x[i];
        }
    }
}

void CGM::solve_LLT(const std::vector<double>& f, std::vector<double>& x) {
    solve_L(f, x);
    solve_LT(x, x);
}

void CGM::solve(std::vector<double>& solution) {
    int max_iter = 10000;
    double eps = 1e-70;

    mul_matrix(x0, r);
    for (int i = 0; i < n; i++) {
        r[i] = rp[i] - r[i];
    }

    make_LLT_decomposition();
    solve_LLT(r, z);

    for (int i = 0; i < n; i++) {
        p[i] = z[i];
    }

    double rp_norm = std::sqrt(dot_prod(rp, rp));
    double alpha, beta;

    for (int iter = 0; iter < max_iter; iter++) {
        double discr = std::sqrt(dot_prod(r, r));
        if (discr / rp_norm < eps) {
            break;
        }

        mul_matrix(z, s);
        alpha = dot_prod(p, r) / dot_prod(s, z);

        for (int i = 0; i < n; i++) {
            x0[i] += alpha * z[i];
            r[i] -= alpha * s[i];
        }

        solve_LLT(r, p);
        beta = dot_prod(p, r) / dot_prod(p, r);

        for (int i = 0; i < n; i++) {
            z[i] = p[i] + beta * z[i];
        }
    }

    solution = x0;
}