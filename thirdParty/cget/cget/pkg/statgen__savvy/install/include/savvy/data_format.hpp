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
    allele = 1,
    genotype,
    genotype_probability,
    genotype_likelihood,
    phred_scaled_genotype_likelihood,
    dosage,
    haplotype_dosage,
    //phase,
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
      case fmt::allele: return ploidy;
      case fmt::genotype: return 1;
      case fmt::genotype_probability: return ploidy + 1;
      case fmt::genotype_likelihood: return ploidy + 1;
      case fmt::phred_scaled_genotype_likelihood: return ploidy + 1;
      case fmt::dosage: return 1;
      case fmt::haplotype_dosage: return ploidy;
    }
  }
}
#endif //LIBSAVVY_DATA_FORMAT_HPP
