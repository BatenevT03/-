#pragma once

#include <vector>
#include <cmath>

class CGM {
public:
    void init(const std::vector<int>& gi_s, const std::vector<int>& gj_s, const std::vector<double>& di_s,
        const std::vector<double>& gg_s, const std::vector<double>& rp_s, int n_s);
    void solve(std::vector<double>& solution);

private:
    void make_LLT_decomposition();
    void mul_matrix(const std::vector<double>& f, std::vector<double>& x);
    void solve_L(const std::vector<double>& f, std::vector<double>& x);
    void solve_LT( std::vector<double>& f, std::vector<double>& x);
    void solve_LLT(const std::vector<double>& f, std::vector<double>& x);
    double dot_prod(const std::vector<double>& a, const std::vector<double>& b);

    int n;
    std::vector<int> gi, gj;
    std::vector<double> di, gg, rp, r, x0, z, p, s;
    std::vector<double> L_di, L_gg;
};