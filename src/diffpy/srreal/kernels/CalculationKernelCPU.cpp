/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 *****************************************************************************/

#include <algorithm>
#include <future>
#include <thread>
#include <vector>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/PeakProfile.hpp>
#include <diffpy/srreal/kernels/PairGeometry.hpp>
#include <diffpy/srreal/kernels/PairTable.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

// Each worker groups a bounded number of pairs before running the SIMD
// windows. Symmetric structures often need only a few distinct windows.
constexpr long long GROUP_BATCH = 65536;
constexpr long long PAIRS_PER_WORKER = 8192;

int workerCount(long long pairs) {
  if (pairs < 2 * PAIRS_PER_WORKER) return 1;
  const long long hardware = std::max(1u, std::thread::hardware_concurrency());
  return int(std::min(pairs / PAIRS_PER_WORKER, hardware));
}

struct Worker {
  PairGroups groups;
  std::vector<double> sum;
  PairSumStats stats;
};

struct CPUState : PairSumWorkspace {
  std::vector<Worker> workers;
};

// Only immutable geometry, scattering arrays and mask data are read by
// the workers. Publish the result in range order after all workers succeed.
template <class GroupPair, class AddPeak>
PairSumStats sumPairs(const PairSumArgs& a, int npts,
                      const GroupPair& groupPair, const AddPeak& addPeak) {
  const long long count = pairCount(a.g);
  if (!count || !npts) return {};
  const int nworkers = workerCount(count);
  auto& workers = static_cast<CPUState&>(*a.workspace).workers;
  workers.resize(nworkers);
  auto run = [&](int i) {
    auto& w = workers[i];
    w.sum.assign(npts, 0.0);
    // Count locally while neighboring workers append to their group vectors.
    PairSumStats stats;
    const auto geometry = a.g;
    // Divide before multiplying so even a large implicit triangle fits.
    const long long stride = count / nworkers;
    const long long rest = count % nworkers;
    const long long first = i * stride + std::min<long long>(i, rest);
    const long long last = first + stride + (i < rest);
    for (long long begin = first; begin < last; begin += GROUP_BATCH) {
      w.groups.clear();
      const long long end = std::min(begin + GROUP_BATCH, last);
      ClusterPair cluster =
        geometry.translations ? ClusterPair{} : ClusterPair::at(begin);
      for (long long k = begin; k < end; ++k) {
        const Pair p = pairAt(geometry, k, cluster);
        if (!geometry.translations) cluster.next();
        if (!p.ok) continue;
        ++stats.pairs;
        if (a.mask && !a.mask->getPairMask(p.site0, p.site1)) continue;
        groupPair(w.groups, p, widthOf(a.w, p.dist, p.msd));
      }
      for (const auto& g : w.groups.groups()) addPeak(w.sum.data(), g);
      stats.peaks += w.groups.groups().size();
    }
    w.stats = stats;
  };
  std::vector<std::future<void>> pending;
  for (int i = 1; i < nworkers; ++i) {
    pending.push_back(std::async(std::launch::async, run, i));
  }
  run(0);
  for (auto& task : pending) task.get();

  PairSumStats stats;
  for (const Worker& w : workers) {
    stats.pairs += w.stats.pairs;
    stats.peaks += w.stats.peaks;
    for (int i = 0; i < npts; ++i) a.out[i] += w.sum[i];
  }
  return stats;
}

}  // namespace

void addProfilePeak(const PeakProfile& profile, double* out, int npts,
                    int rlosteps, double rstep, double dist, double fwhm,
                    double scale) {
  const int ilo =
    gridIndex(dist + profile.xboundlo(fwhm), rstep, rlosteps, 0, npts);
  const int ihi =
    gridIndex(dist + profile.xboundhi(fwhm), rstep, rlosteps, 1, npts);
  profile.addRDFPeak(out, ilo, ihi, rlosteps, rstep, dist, fwhm, scale);
}

PairSumStats cpuPDFPairSum(const PDFPairSumArgs& a, GaussianWindow gaussian) {
  auto group = [&](PairGroups& groups, const Pair& p, double width) {
    groups.add(p.dist, width, 0,
               a.sfsite[p.site0] * a.sfsite[p.site1] * p.weight);
  };
  auto peak = [&](double* out, const PairGroups::Group& g) {
    if (g.width <= 0) return;
    const double half = a.halfboundrel * g.width;
    const int ilo = gridIndex(g.dist - half, a.rstep, a.rlosteps, 0, a.npts);
    const int ihi = gridIndex(g.dist + half, a.rstep, a.rlosteps, 1, a.npts);
    gaussian({out, ilo, ihi, a.rlosteps, a.rstep, g.dist, g.width, g.scale,
              a.height / g.width, bool(a.cropped), a.halfboundrel});
  };
  return sumPairs(a, a.npts, group, peak);
}

PairSumStats cpuDebyePairSum(const DebyePairSumArgs& a, DebyeSum debye) {
  if (a.kqlo >= a.kqhi) return {};
  auto group = [&](PairGroups& groups, const Pair& p, double width) {
    const int t0 = a.typeofsite[p.site0], t1 = a.typeofsite[p.site1];
    const auto typepair = std::int64_t(t0) * a.ntypes + t1;
    groups.add(p.dist, fwhmToSigma(width), typepair, 1.0, p.weight);
  };
  auto peak = [&](double* out, const PairGroups::Group& g) {
    // Test the cutoff at the scale of one pair before multiplying by count.
    const int t0 = int(g.tag / a.ntypes), t1 = int(g.tag % a.ntypes);
    debye({out, a.kqlo, a.kqhi, a.qstep, g.dist, g.width, double(g.weight),
           a.prec, a.sf + size_t(a.curve[t0]) * a.nq,
           a.sf + size_t(a.curve[t1]) * a.nq, g.scale, a.occupancy[t0],
           a.occupancy[t1]});
  };
  return sumPairs(a, a.nq, group, peak);
}

PairSumWorkspace* makeCPUWorkspace() { return new CPUState; }

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
