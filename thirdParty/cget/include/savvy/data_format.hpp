/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_DATA_FORMAT_HPP
#define LIBSAVVY_DATA_FORMAT_HPP

#include <cstdint>

namespace savvy
{
  enum class fmt : std::uint8_t
  {
    gt = 1,
    ac,
    gp,
    gl,
    pl,
    ds,
    hds,

//    gt = genotype,
//    gp = genotype_probability,
//    gl = genotype_likelihood,
//    pl = phred_scaled_genotype_likelihood,
//    ds = dosage
//    ec = dosage
  };

  inline std::uint64_t sample_stride(fmt format, std::uint64_t ploidy)
  {
    switch (format)
    {
      case fmt::gt: return ploidy;
      case fmt::ac: return 1;
      case fmt::gp: return ploidy + 1;
      case fmt::gl: return ploidy + 1;
      case fmt::pl: return ploidy + 1;
      case fmt::ds: return 1;
      case fmt::hds: return ploidy;
    }
  }
}
#endif //LIBSAVVY_DATA_FORMAT_HPP
