/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <chrono>
#include <random>
#include <string>
#include <utility>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>

#include "korder.hh"
#include "SylvException.hh"

struct Rand
{
  static std::mt19937 mtgen;
  static std::uniform_real_distribution<> dis;
  static void init(int n1, int n2, int n3, int n4, int n5);
  static double get(double m);
  static int get(int m);
  static bool discrete(double prob); // answers true with given probability
};

std::mt19937 Rand::mtgen;
std::uniform_real_distribution<> Rand::dis;

ConstTwoDMatrix
make_matrix(int rows, int cols, const double *p)
{
  return ConstTwoDMatrix{rows, cols, ConstVector{p, rows*cols}};
}

void
Rand::init(int n1, int n2, int n3, int n4, int n5)
{
  decltype(mtgen)::result_type seed = n1;
  seed = 256*seed+n2;
  seed = 256*seed+n3;
  seed = 256*seed+n4;
  seed = 256*seed+n5;
  mtgen.seed(seed);
}

double
Rand::get(double m)
{
  return 2*m*(dis(mtgen)-0.5);
}

int
Rand::get(int m)
{
  return static_cast<int>(get(0.9999*m));
}

bool
Rand::discrete(double prob)
{
  return dis(mtgen) < prob;
}

struct SparseGenerator
{
  static std::unique_ptr<FSSparseTensor> makeTensor(int dim, int nv, int r,
                                                    double fill, double m);
  static void fillContainer(TensorContainer<FSSparseTensor> &c,
                            int maxdim, int nv, int r, double m);
};

std::unique_ptr<FSSparseTensor>
SparseGenerator::makeTensor(int dim, int nv, int r,
                            double fill, double m)
{
  auto res = std::make_unique<FSSparseTensor>(dim, nv, r);
  FFSTensor dummy(0, nv, dim);
  for (Tensor::index fi = dummy.begin(); fi != dummy.end(); ++fi)
    for (int i = 0; i < r; i++)
      if (Rand::discrete(fill))
        {
          double x = Rand::get(m);
          res->insert(fi.getCoor(), i, x);
        }
  return res;
}

void
SparseGenerator::fillContainer(TensorContainer<FSSparseTensor> &c,
                               int maxdim, int nv, int r,
                               double m)
{
  Rand::init(maxdim, nv, r, static_cast<int>(5*m), 0);
  double fill = 0.5;
  for (int d = 1; d <= maxdim; d++)
    {
      c.insert(makeTensor(d, nv, r, fill, m));
      fill *= 0.3;
    }
}

const double vdata[] =
  { // 3x3
   0.1307870268, 0.1241940078, 0.1356703123,
   0.1241940078, 0.1986920419, 0.2010160581,
   0.1356703123, 0.2010160581, 0.2160336975
  };

const double gy_data[] =
  { // 8x4
   0.3985178619, -0.5688233582, 0.9572900437, -0.6606847776, 0.1453004017,
   0.3025310675, -0.8627437750, -0.6903410191, 0.4751910580, -0.7270018589,
   -0.0939612498, -0.1463831989, 0.6742110220, 0.6046671043, 0.5215893126,
   -1.0412969986, -0.3524898417, -1.0986703430, 0.8006531522, 0.8879776376,
   -0.1037608317, -0.5587378073, -0.1010366945, 0.9462411248, -0.2439199881,
   1.3420621236, -0.7820285935, 0.3205293447, 0.3606124791, 0.2975422208,
   -0.5452861965, 1.6320340279
  };

const double gu_data[] =
  { // just some numbers, no structure
   1.8415286914, -0.2638743845, 1.7690713274, 0.9668585956, 0.2303143646,
   -0.2229624279, -0.4381991822, 1.0082401405, -0.3186555860, -0.0624691529,
   -0.5189085756, 1.4269672156, 0.1163282969, 1.4020183445, -0.0952660426,
   0.2099097124, 0.6912400502, -0.5180935114, 0.5288316624, 0.2188053448,
   0.5715516767, 0.7813893410, -0.6385073106, 0.8335131513, 0.3605202168,
   -1.1167944865, -1.2263750934, 0.6113636081, 0.6964915482, -0.6451217688,
   0.4062810500, -2.0552251116, -1.6383406284, 0.0198915095, 0.0111014458,
   -1.2421792262, -1.0724161722, -0.4276904972, 0.1801494950, -2.0716473264
  };

const double vdata2[] =
  { // 10×10 positive definite
   0.79666, -0.15536, 0.05667, -0.21026, 0.20262, 0.28505, 0.60341, -0.09703, 0.32363, 0.13299,
   -0.15536, 0.64380, -0.01131, 0.00980, 0.03755, 0.43791, 0.21784, -0.31755, -0.55911, -0.29655,
   0.05667, -0.01131, 0.56165, -0.34357, -0.40584, 0.20990, 0.28348, 0.20398, -0.19856, 0.35820,
   -0.21026, 0.00980, -0.34357, 0.56147, 0.10972, -0.34146, -0.49906, -0.19685, 0.21088, -0.31560,
   0.20262, 0.03755, -0.40584, 0.10972, 0.72278, 0.02155, 0.04089, -0.19696, 0.03446, -0.12919,
   0.28505, 0.43791, 0.20990, -0.34146, 0.02155, 0.75867, 0.77699, -0.31125, -0.55141, -0.02155,
   0.60341, 0.21784, 0.28348, -0.49906, 0.04089, 0.77699, 1.34553, -0.18613, -0.25811, -0.19016,
   -0.09703, -0.31755, 0.20398, -0.19685, -0.19696, -0.31125, -0.18613, 0.59470, 0.08386, 0.41750,
   0.32363, -0.55911, -0.19856, 0.21088, 0.03446, -0.55141, -0.25811, 0.08386, 0.98917, -0.12992,
   0.13299, -0.29655, 0.35820, -0.31560, -0.12919, -0.02155, -0.19016, 0.41750, -0.12992, 0.89608
  };

const double gy_data2[] =
  { // 600 items make gy 30×20, whose gy(6:25,:) has spectrum within unit
   0.39414, -0.29766, 0.08948, -0.19204, -0.00750, 0.21159, 0.05494, 0.06225, 0.01771, 0.21913,
   -0.01373, 0.20086, -0.06086, -0.10955, 0.14424, -0.08390, 0.03948, -0.14713, 0.11674, 0.05091,
   0.24039, 0.28307, -0.11835, 0.13030, 0.11682, -0.27444, -0.19311, -0.16654, 0.12867, 0.25116,
   -0.19781, 0.45242, -0.15862, 0.24428, -0.11966, 0.11483, -0.32279, 0.29727, 0.20934, -0.18190,
   -0.15080, -0.09477, -0.30551, -0.02672, -0.26919, 0.11165, -0.06390, 0.03449, -0.26622, 0.22197,
   0.45141, -0.41683, 0.09760, 0.31094, -0.01652, 0.05809, -0.04514, -0.05645, 0.00554, 0.47980,
   0.11726, 0.42459, -0.13136, -0.30902, -0.14648, 0.11455, 0.02947, -0.03835, -0.04044, 0.03559,
   -0.26575, -0.01783, 0.31243, -0.14412, -0.13218, -0.05080, 0.18576, 0.13840, -0.05560, 0.35530,
   -0.25573, -0.11560, 0.15187, -0.18431, 0.08193, -0.32278, 0.17560, -0.05529, -0.10020, -0.23088,
   -0.20979, -0.49245, 0.09915, -0.16909, -0.03443, 0.19497, 0.18473, 0.25662, 0.29605, -0.20531,
   -0.39244, -0.43369, 0.05588, 0.24823, -0.14236, -0.08311, 0.16371, -0.19975, 0.30605, -0.17087,
   -0.01270, 0.00123, -0.22426, -0.13810, 0.05079, 0.06971, 0.01922, -0.09952, -0.23177, -0.41962,
   -0.41991, 0.41430, -0.04247, -0.13706, -0.12048, -0.28906, -0.22813, -0.25057, -0.18579, -0.20642,
   -0.47976, 0.25490, -0.05138, -0.30794, 0.31651, 0.02034, 0.12954, -0.20110, 0.13336, -0.40775,
   -0.30195, -0.13704, 0.12396, 0.28152, 0.02986, 0.27669, 0.24623, 0.08635, -0.11956, -0.02949,
   0.37401, 0.20838, 0.24801, -0.26872, 0.11195, 0.00315, -0.19069, 0.12839, -0.23036, -0.48228,
   0.08434, -0.39872, -0.28896, -0.28754, 0.24668, 0.23285, 0.25437, 0.10456, -0.14124, 0.20483,
   -0.19117, -0.33836, -0.24875, 0.08207, -0.03930, 0.20364, 0.15384, -0.15270, 0.24372, -0.11199,
   -0.46591, 0.30319, 0.05745, 0.09084, 0.06058, 0.31884, 0.05071, -0.28899, -0.30793, -0.03566,
   0.02286, 0.28178, 0.00736, -0.31378, -0.18144, -0.22346, -0.27239, 0.31043, -0.26228, 0.22181,
   -0.15096, -0.36953, -0.06032, 0.21496, 0.29545, -0.13112, 0.16420, -0.07573, -0.43111, -0.43057,
   0.26716, -0.31209, -0.05866, -0.29101, -0.27437, -0.18727, 0.28732, -0.19014, 0.08837, 0.30405,
   0.06103, -0.35612, 0.00173, 0.25134, -0.08987, -0.22766, -0.03254, -0.18662, -0.08491, 0.49401,
   -0.12145, -0.02961, -0.03668, -0.30043, -0.08555, 0.01701, -0.12544, 0.10969, -0.48202, 0.07245,
   0.20673, 0.11408, 0.04343, -0.01815, -0.31594, -0.23632, -0.06258, -0.27474, 0.12180, 0.16613,
   -0.37931, 0.30219, 0.15765, 0.25489, 0.17529, -0.17020, -0.30060, 0.22058, -0.02450, -0.42143,
   0.49642, 0.46899, -0.28552, -0.22549, -0.01333, 0.21567, 0.22251, 0.21639, -0.19194, -0.19140,
   -0.24106, 0.10952, -0.11019, 0.29763, -0.02039, -0.25748, 0.23169, 0.01357, 0.09802, -0.19022,
   0.37604, -0.40777, 0.18131, -0.10258, 0.29573, -0.31773, 0.09069, -0.02198, -0.26594, 0.48302,
   -0.10041, 0.20210, -0.05609, -0.01169, -0.17339, 0.17862, -0.22502, 0.29009, -0.45160, 0.19771,
   0.27634, 0.31695, -0.09993, 0.17167, 0.12394, 0.28088, -0.12502, -0.16967, -0.06296, -0.17036,
   0.27320, 0.01595, 0.16955, 0.30146, -0.15173, -0.29807, 0.08178, -0.06811, 0.21655, 0.26348,
   0.06316, 0.45661, -0.29756, -0.05742, -0.14715, -0.03037, -0.16656, -0.08768, 0.38078, 0.40679,
   -0.32779, -0.09106, 0.16107, -0.07301, 0.07700, -0.22694, -0.15692, -0.02548, 0.38749, -0.12203,
   -0.02980, -0.22067, 0.00680, -0.23058, -0.29112, 0.23032, -0.16026, 0.23392, -0.09990, 0.03628,
   -0.42592, -0.33474, -0.09499, -0.17442, -0.20110, 0.24618, -0.06418, -0.06715, 0.40754, 0.29377,
   0.29543, -0.16832, -0.08468, 0.06491, -0.01410, 0.19988, 0.24950, 0.14626, -0.27851, 0.06079,
   0.48134, -0.13475, 0.25398, 0.11738, 0.23369, -0.00661, -0.16811, -0.04557, -0.12030, -0.39527,
   -0.35760, 0.01840, -0.15941, 0.03290, 0.09988, -0.08307, 0.06644, -0.24637, 0.34112, -0.08026,
   0.00951, 0.27656, 0.16247, 0.28217, 0.17198, -0.16389, -0.03835, -0.02675, -0.08032, -0.21045,
   -0.38946, 0.23207, 0.10987, -0.31674, -0.28653, -0.27430, -0.29109, -0.00648, 0.38431, -0.38478,
   -0.41195, -0.19364, -0.20977, -0.05524, 0.05558, -0.20109, 0.11803, -0.19884, 0.43318, -0.39255,
   0.26612, -0.21771, 0.12471, 0.12856, -0.15104, -0.11676, 0.17582, -0.25330, 0.00298, -0.31712,
   0.21532, -0.20319, 0.14507, -0.04588, -0.22995, -0.06470, 0.18849, -0.13444, 0.37107, 0.07387,
   -0.14008, 0.09896, 0.13727, -0.28417, -0.09461, -0.18703, 0.04080, 0.02343, -0.49988, 0.17993,
   0.23189, -0.30581, -0.18334, -0.09667, -0.27699, -0.05998, 0.09118, -0.32453, 0.46251, 0.41500,
   -0.45314, -0.00544, 0.08529, 0.29099, -0.00937, -0.31650, 0.26163, 0.14506, 0.37498, -0.16454,
   0.35215, 0.31642, -0.09161, -0.31452, -0.04792, -0.04677, -0.19523, 0.27998, 0.05491, 0.44461,
   -0.01258, -0.27887, 0.18361, -0.04539, -0.02977, 0.30821, 0.29454, -0.17932, 0.16193, 0.23934,
   0.47923, 0.25373, 0.23258, 0.31484, -0.17958, -0.01136, 0.17681, 0.12869, 0.03235, 0.43762,
   0.13734, -0.09433, -0.03735, 0.17949, 0.14122, -0.17814, 0.06359, 0.16044, 0.12249, -0.22314,
   0.40775, 0.05147, 0.12389, 0.04290, -0.01642, 0.00082, -0.18056, 0.02875, 0.32690, 0.17712,
   0.34001, -0.21581, -0.01086, -0.18180, 0.17480, -0.17774, -0.07503, 0.28438, -0.19747, 0.29595,
   -0.28002, -0.02073, -0.16522, -0.18234, -0.20565, 0.29620, 0.07502, 0.01429, -0.31418, 0.43693,
   -0.12212, 0.11178, -0.28503, 0.04683, 0.00072, 0.05566, 0.18857, 0.26101, -0.38891, -0.21216,
   -0.21850, -0.15147, -0.30749, -0.23762, 0.14984, 0.03535, -0.02862, -0.00105, -0.39907, -0.06909,
   -0.36094, 0.21717, 0.15930, -0.18924, 0.13741, 0.01039, 0.13613, 0.00659, 0.07676, -0.13711,
   0.24285, -0.07564, -0.28349, -0.15658, 0.03135, -0.30909, -0.22534, 0.17363, -0.19376, 0.26038,
   0.05546, -0.22607, 0.32420, -0.02552, -0.05400, 0.13388, 0.04643, -0.31535, -0.06181, 0.30237,
   -0.04680, -0.29441, 0.12231, 0.03960, -0.01188, 0.01406, 0.25402, 0.03315, 0.25026, -0.10922
  };

const double gu_data2[] =
  { // raw data 300 items
   0.26599, 0.41329, 0.31846, 0.92590, 0.43050, 0.17466, 0.02322, 0.72621, 0.37921, 0.70597,
   0.97098, 0.14023, 0.57619, 0.09938, 0.02281, 0.92341, 0.72654, 0.71000, 0.76687, 0.70182,
   0.88752, 0.49524, 0.42549, 0.42806, 0.57615, 0.76051, 0.15341, 0.47457, 0.60066, 0.40880,
   0.20668, 0.41949, 0.97620, 0.94318, 0.71491, 0.56402, 0.23553, 0.94387, 0.78567, 0.06362,
   0.85252, 0.86262, 0.25190, 0.03274, 0.93216, 0.37971, 0.08797, 0.14596, 0.73871, 0.06574,
   0.67447, 0.28575, 0.43911, 0.92133, 0.12327, 0.87762, 0.71060, 0.07141, 0.55443, 0.53310,
   0.91529, 0.25121, 0.07593, 0.94490, 0.28656, 0.82174, 0.68887, 0.67337, 0.99291, 0.03316,
   0.02849, 0.33891, 0.25594, 0.90071, 0.01248, 0.67871, 0.65953, 0.65369, 0.97574, 0.31578,
   0.23678, 0.39220, 0.06706, 0.80943, 0.57694, 0.08220, 0.18151, 0.19969, 0.37096, 0.37858,
   0.70153, 0.46816, 0.76511, 0.02520, 0.39387, 0.25527, 0.39050, 0.60141, 0.30322, 0.46195,
   0.12025, 0.33616, 0.04174, 0.00196, 0.68886, 0.74445, 0.15869, 0.18994, 0.95195, 0.62874,
   0.82874, 0.53369, 0.34383, 0.50752, 0.97023, 0.22695, 0.62407, 0.25840, 0.71279, 0.28785,
   0.31611, 0.20391, 0.19702, 0.40760, 0.85158, 0.68369, 0.63760, 0.09879, 0.11924, 0.32920,
   0.53052, 0.15900, 0.21229, 0.84080, 0.33933, 0.93651, 0.42705, 0.06199, 0.50092, 0.47192,
   0.57152, 0.01818, 0.31404, 0.50173, 0.87725, 0.50530, 0.10717, 0.04035, 0.32901, 0.33538,
   0.04780, 0.40984, 0.78216, 0.91288, 0.11314, 0.25248, 0.23823, 0.74001, 0.48089, 0.55531,
   0.82486, 0.01058, 0.05409, 0.44357, 0.52641, 0.68188, 0.94629, 0.61627, 0.33037, 0.11961,
   0.57988, 0.19653, 0.91902, 0.59838, 0.52974, 0.28364, 0.45767, 0.65836, 0.63045, 0.76140,
   0.27918, 0.27256, 0.46035, 0.77418, 0.92918, 0.14095, 0.89645, 0.25146, 0.21172, 0.47910,
   0.95451, 0.34377, 0.29927, 0.79220, 0.97654, 0.67591, 0.44385, 0.38434, 0.44860, 0.28170,
   0.90712, 0.20337, 0.00292, 0.55046, 0.62255, 0.45127, 0.80896, 0.43965, 0.59145, 0.23801,
   0.33601, 0.30119, 0.89935, 0.40850, 0.98226, 0.75430, 0.68318, 0.65407, 0.68067, 0.32942,
   0.11756, 0.27626, 0.83879, 0.72174, 0.75430, 0.13702, 0.03402, 0.58781, 0.07393, 0.23067,
   0.92537, 0.29445, 0.43437, 0.47685, 0.54548, 0.66082, 0.23805, 0.60208, 0.94337, 0.21363,
   0.72637, 0.57181, 0.77679, 0.63931, 0.72860, 0.38901, 0.94920, 0.04535, 0.12863, 0.40550,
   0.90095, 0.21418, 0.13953, 0.99639, 0.02526, 0.70018, 0.21828, 0.20294, 0.20191, 0.30954,
   0.39490, 0.68955, 0.11506, 0.15748, 0.40252, 0.91680, 0.61547, 0.78443, 0.19693, 0.67630,
   0.56552, 0.58556, 0.53554, 0.53507, 0.09831, 0.21229, 0.83135, 0.26375, 0.89287, 0.97069,
   0.70615, 0.42041, 0.43117, 0.21291, 0.26086, 0.26978, 0.77340, 0.43833, 0.46179, 0.54418,
   0.67878, 0.42776, 0.61454, 0.55915, 0.36363, 0.31999, 0.42442, 0.86649, 0.62513, 0.02047
  };

class TestRunnable
{
public:
  const std::string name;
  int dim; // dimension of the solved problem
  int nvar; // number of variable of the solved problem
  TestRunnable(std::string n, int d, int nv)
    : name{std::move(n)}, dim(d), nvar(nv)
  {
  }
  virtual ~TestRunnable() = default;
  bool test() const;
  virtual bool run() const = 0;
protected:
  static double korder_unfold_fold(int maxdim, int unfold_dim,
                                   int nstat, int npred, int nboth, int forw,
                                   const TwoDMatrix &gy, const TwoDMatrix &gu,
                                   const TwoDMatrix &v);
};

bool
TestRunnable::test() const
{
  std::cout << "Running test <" << name << ">" << std::endl;
  clock_t start = clock();
  auto start_real = std::chrono::steady_clock::now();
  bool passed = run();
  clock_t end = clock();
  auto end_real = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = end_real - start_real;
  std::cout << "CPU time  " << std::setprecision(4) << std::setw(8)
            << static_cast<double>(end-start)/CLOCKS_PER_SEC << " (CPU seconds)\n"
            << "Real time " << std::setw(8) << duration.count() << " (seconds).....................";
  if (passed)
    std::cout << "passed\n\n";
  else
    std::cout << "FAILED\n\n";
  return passed;
}

double
TestRunnable::korder_unfold_fold(int maxdim, int unfold_dim,
                                 int nstat, int npred, int nboth, int nforw,
                                 const TwoDMatrix &gy, const TwoDMatrix &gu,
                                 const TwoDMatrix &v)
{
  TensorContainer<FSSparseTensor> c(1);
  int ny = nstat+npred+nboth+nforw;
  int nu = v.nrows();
  int nz = nboth+nforw+ny+nboth+npred+nu;
  SparseGenerator::fillContainer(c, maxdim, nz, ny, 5.0);
  for (int d = 1; d <= maxdim; d++)
    std::cout << "\ttensor fill for dim=" << d << " is:   "
              << std::setprecision(2) << std::setw(6) << std::fixed
              << c.get(Symmetry{d}).getFillFactor()*100.0 << " %\n"
              << std::defaultfloat;
  Journal jr("out.txt");
  KOrder kord(nstat, npred, nboth, nforw, c, gy, gu, v, jr);
  // Perform unfolded steps until unfold_dim
  double maxerror = 0.0;
  for (int d = 2; d <= unfold_dim; d++)
    {
      clock_t pertime = clock();
      kord.performStep<Storage::unfold>(d);
      pertime = clock()-pertime;
      std::cout << "\ttime for unfolded step dim=" << d << ": " << std::setprecision(4)
                << static_cast<double>(pertime)/CLOCKS_PER_SEC << std::endl;
      clock_t checktime = clock();
      double err = kord.check<Storage::unfold>(d);
      checktime = clock()-checktime;
      std::cout << "\ttime for step check dim=" << d << ":    " << std::setprecision(4)
                << static_cast<double>(checktime)/CLOCKS_PER_SEC << '\n'
                << "\tmax error in step dim=" << d << ":      " << std::setprecision(6) << err
                << std::endl;
      maxerror = std::max(err, maxerror);
    }
  // Perform folded steps until maxdim
  if (unfold_dim < maxdim)
    {
      clock_t swtime = clock();
      kord.switchToFolded();
      swtime = clock()-swtime;
      std::cout << "\ttime for switching dim=" << unfold_dim << ":     " << std::setprecision(4)
                << static_cast<double>(swtime)/CLOCKS_PER_SEC << std::endl;

      for (int d = unfold_dim+1; d <= maxdim; d++)
        {
          clock_t pertime = clock();
          kord.performStep<Storage::fold>(d);
          pertime = clock()-pertime;
          std::cout << "\ttime for folded step dim=" << d << ":   " << std::setprecision(4)
                    << static_cast<double>(pertime)/CLOCKS_PER_SEC << std::endl;
          clock_t checktime = clock();
          double err = kord.check<Storage::fold>(d);
          checktime = clock()-checktime;
          std::cout << "\ttime for step check dim=" << d << ":    " << std::setprecision(4)
                    << static_cast<double>(checktime)/CLOCKS_PER_SEC << '\n'
                    << "\tmax error in step dim=" << d << ":      " << std::setprecision(6) << err
                    << std::endl;
          maxerror = std::max(err, maxerror);
        }
    }
  return maxerror;
}

class UnfoldKOrderSmall : public TestRunnable
{
public:
  UnfoldKOrderSmall()
    : TestRunnable("unfold-3 fold-4 korder (stat=2,pred=3,both=1,forw=2,u=3,dim=4)",
                   4, 18)
  {
  }

  bool
  run() const override
  {
    TwoDMatrix gy{make_matrix(8, 4, gy_data)};
    TwoDMatrix gu{make_matrix(8, 3, gu_data)};
    TwoDMatrix v{make_matrix(3, 3, vdata)};
    double err = korder_unfold_fold(4, 3, 2, 3, 1, 2,
                                    gy, gu, v);

    return err < 5e-7;
  }
};

// Same dimension as Smets & Wouters
class UnfoldKOrderSW : public TestRunnable
{
public:
  UnfoldKOrderSW()
    : TestRunnable("unfold S&W korder (stat=5,pred=12,both=8,forw=5,u=10,dim=4)",
                   4, 73)
  {
  }

  bool
  run() const override
  {
    TwoDMatrix gy{make_matrix(30, 20, gy_data2)};
    TwoDMatrix gu{make_matrix(30, 10, gu_data2)};
    TwoDMatrix v{make_matrix(10, 10, vdata2)};
    v.mult(0.001);
    gu.mult(.01);
    double err = korder_unfold_fold(4, 4, 5, 12, 8, 5,
                                    gy, gu, v);

    return err < 0.5;
  }
};

class UnfoldFoldKOrderSW : public TestRunnable
{
public:
  UnfoldFoldKOrderSW()
    : TestRunnable("unfold-2 fold-3 S&W korder (stat=5,pred=12,both=8,forw=5,u=10,dim=3)",
                   4, 73)
  {
  }

  bool
  run() const override
  {
    TwoDMatrix gy{make_matrix(30, 20, gy_data2)};
    TwoDMatrix gu{make_matrix(30, 10, gu_data2)};
    TwoDMatrix v{make_matrix(10, 10, vdata2)};
    v.mult(0.001);
    gu.mult(.01);
    double err = korder_unfold_fold(4, 3, 5, 12, 8, 5,
                                    gy, gu, v);

    return err < 0.5;
  }
};

int
main()
{
  std::vector<std::unique_ptr<TestRunnable>> all_tests;
  // Fill in vector of all tests
  all_tests.push_back(std::make_unique<UnfoldKOrderSmall>());
  all_tests.push_back(std::make_unique<UnfoldKOrderSW>());
  all_tests.push_back(std::make_unique<UnfoldFoldKOrderSW>());

  // Find maximum dimension and maximum nvar
  int dmax = 0;
  int nvmax = 0;
  for (const auto &test : all_tests)
    {
      if (dmax < test->dim)
        dmax = test->dim;
      if (nvmax < test->nvar)
        nvmax = test->nvar;
    }
  TLStatic::init(dmax, nvmax); // initialize library

  // Launch the tests
  int success = 0;
  for (const auto &test : all_tests)
    {
      try
        {
          if (test->test())
            success++;
        }
      catch (const TLException &e)
        {
          std::cout << "Caught TL exception in <" << test->name << ">:" << std::endl;
          e.print();
        }
      catch (SylvException &e)
        {
          std::cout << "Caught Sylv exception in <" << test->name << ">:" << std::endl;
          e.printMessage();
        }
    }

  int nfailed = all_tests.size() - success;
  std::cout << "There were " << nfailed << " tests that failed out of "
            << all_tests.size() << " tests run." << std::endl;

  if (nfailed)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
