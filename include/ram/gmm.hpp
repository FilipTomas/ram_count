// gmm.hpp ────────────────────────────────────────────────────────────────
#pragma once
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cstdint>

namespace gmm1d {

struct Result {
    std::vector<double> means;
    std::vector<double> variances;
    std::vector<double> weights;
};

// Forward declaration (vector version)
Result fit(const std::vector<std::size_t>& hist,
           std::size_t n_components = 3,
           std::size_t max_iter     = 100,
           double      eps          = 1e-3);

// ────────────────────────────────────────────────────────────────────────
// Wrapper: build histogram only for depths 2 … 200 (inclusive)
// ────────────────────────────────────────────────────────────────────────
inline Result
fit(const std::unordered_map<std::uint64_t, std::size_t>& kmer_counts,
    std::size_t n_components = 3,
    std::size_t max_iter     = 100,
    double      eps          = 1e-3)
{
    constexpr std::size_t MIN_D = 2;
    constexpr std::size_t MAX_D = 200;

    std::vector<std::size_t> hist(MAX_D + 1, 0);       // 0…200

    for (const auto& kv : kmer_counts) {
        std::size_t d = kv.second;
        if (d >= MIN_D && d <= MAX_D)
            ++hist[d];
    }
    return fit(hist, n_components, max_iter, eps);
}

// ────────────────────────────────────────────────────────────────────────
// Vector-based EM solver (ignores depths 0-1 and >200)
// ────────────────────────────────────────────────────────────────────────
inline Result
fit(const std::vector<std::size_t>& src_hist,
    std::size_t n_components,
    std::size_t max_iter,
    double      eps)
{
    constexpr std::size_t MIN_D = 2;
    constexpr std::size_t MAX_D = 200;

    // build a trimmed copy so bins 0-1 and >200 are zeroed
    std::vector<std::size_t> hist(MAX_D + 1, 0);
    for (std::size_t d = MIN_D; d <= MAX_D && d < src_hist.size(); ++d)
        hist[d] = src_hist[d];

    const std::size_t D = hist.size();                 // = 201
    if (n_components == 0)
        return {{}, {}, {}};

    const double TWO_PI = 6.2831853071795864769;
    auto gaussian = [&](double x, double mu, double var) {
        double diff = x - mu;
        return std::exp(-0.5 * diff * diff / var) /
               std::sqrt(TWO_PI * var);
    };

    double total = 0.0;
    for (double c : hist) total += c;
    if (total == 0.0) return {{}, {}, {}};             // nothing in range

    // ── initial params (unchanged) ──────────────────────────────────────
    std::vector<double> means(n_components);
    std::vector<double> vars (n_components, 100.0);
    std::vector<double> mix  (n_components, 1.0 / n_components);

    {   std::vector<double> cdf(D);
        double acc = 0.0;
        for (std::size_t d = 0; d < D; ++d) {
            acc += hist[d];
            cdf[d] = acc / total;
        }
        for (std::size_t k = 0; k < n_components; ++k) {
            double q = (k + 0.5) / n_components;
            auto it  = std::lower_bound(cdf.begin(), cdf.end(), q);
            means[k] = static_cast<double>(std::distance(cdf.begin(), it));
        }
    }

    std::vector<std::vector<double>> resp(n_components,
                                          std::vector<double>(D));

    // ── EM iterations (unchanged) ───────────────────────────────────────
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        for (std::size_t d = 0; d < D; ++d) {
            if (hist[d] == 0) { for (auto& r : resp) r[d] = 0.0; continue; }
            double denom = 0.0;
            for (std::size_t k = 0; k < n_components; ++k) {
                resp[k][d] = mix[k] * gaussian(d, means[k], vars[k]);
                denom += resp[k][d];
            }
            if (denom < std::numeric_limits<double>::min()) {
                for (auto& r : resp) r[d] = 1.0 / n_components;
            } else {
                for (auto& r : resp) r[d] /= denom;
            }
        }

        bool converged = true;
        for (std::size_t k = 0; k < n_components; ++k) {
            double Nk = 0.0, mu = 0.0, var = 0.0;
            for (std::size_t d = 0; d < D; ++d) {
                double r = resp[k][d] * hist[d];
                Nk += r; mu += r * d;
            }
            if (Nk == 0.0) Nk = 1e-6;
            mu /= Nk;
            for (std::size_t d = 0; d < D; ++d) {
                double diff = d - mu;
                var += resp[k][d] * hist[d] * diff * diff;
            }
            var /= Nk; if (var < 1.0) var = 1.0;

            converged &= (std::fabs(mu - means[k]) < eps &&
                          std::fabs(var - vars[k]) < eps);

            means[k] = mu; vars[k] = var; mix[k] = Nk / total;
        }
        if (converged) break;
    }

    std::vector<std::size_t> idx(n_components);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b){ return means[a] < means[b]; });

    Result res;
    res.means.resize(n_components);
    res.variances.resize(n_components);
    res.weights.resize(n_components);
    for (std::size_t i = 0; i < n_components; ++i) {
        res.means[i]     = means[idx[i]];
        res.variances[i] = vars [idx[i]];
        res.weights[i]   = mix  [idx[i]];
    }
    return res;
}

} // namespace gmm1d
