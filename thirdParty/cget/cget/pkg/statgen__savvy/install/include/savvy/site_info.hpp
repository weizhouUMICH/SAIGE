/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SITE_INFO_HPP
#define LIBSAVVY_SITE_INFO_HPP

#include "compressed_vector.hpp"
#include "data_format.hpp"

#include <string>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <ostream>

namespace savvy
{
  class site_info
  {
  public:

    site_info()
    {
    }

    site_info(
      std::string&& chromosome,
      std::uint64_t pos,
      std::string&& ref,
      std::string&& alt,
      std::unordered_map<std::string, std::string>&& properties)
      :
      properties_(std::move(properties)),
      chromosome_(std::move(chromosome)),
      ref_(std::move(ref)),
      alt_(std::move(alt)),
      position_(pos)
    {

    }

    virtual ~site_info() {}

    const std::string& chromosome() const { return chromosome_; }
    const std::string& ref() const { return ref_; }
    const std::string& alt() const { return alt_; }
    //[[deprecated]] std::uint64_t locus() const { return position_; }
    std::uint64_t position() const { return position_; }
    const std::string& prop(const std::string& key) const
    {
      auto it = properties_.find(key);
      if (it == properties_.end())
        return empty_string;
      return it->second;
    }

    void prop(const std::string& key, std::string value)
    {
      properties_[key] = std::move(value);
    }
  private:
    std::unordered_map<std::string, std::string> properties_;
    std::string chromosome_;
    std::string ref_;
    std::string alt_;
    std::uint64_t position_;
    static const std::string empty_string;
  };

  template <typename T>
  class variant : public site_info
  {
  public:
    T& data() { return data_; }
    const T& data() const { return data_; }
  private:
    T data_;
  };

//  template <typename T>
//  class allele_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_probabilities_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class dosage_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class genotype_likelihoods_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };
//
//  template <typename T>
//  class phred_genotype_likelihoods_vector : public variant_vector<T>
//  {
//  public:
//    using variant_vector<T>::variant_vector;
//  };


//  template <typename T>
//  using dense_allele_vector = allele_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_allele_vector = allele_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_genotype_vector = genotype_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_genotype_vector = genotype_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_genotype_probabilities_vector = genotype_probabilities_vector<std::vector<T>>;
//  template <typename T>
//  using dense_genotype_likelihoods_vector = genotype_likelihoods_vector<std::vector<T>>;
//  template <typename T>
//  using dense_phred_genotype_likelihoods_vector = phred_genotype_likelihoods_vector<std::vector<T>>;
////  template <typename T>
////  using sparse_genotype_probabilities_vector = genotype_probabilities_vector<compressed_vector<T>>;
//
//  template <typename T>
//  using dense_dosage_vector = dosage_vector<std::vector<T>>;
//  template <typename T>
//  using sparse_dosage_vector = dosage_vector<compressed_vector<T>>;

  namespace detail
  {
    void print_vcf_site_info(std::ostream& out, const site_info& in, const std::vector<std::string>& info_fields);
  }

  template <typename T>
  void print_vcf_record(std::ostream& out, const site_info& in, const std::vector<T>& in_data, const std::vector<std::string>& info_fields, std::uint32_t ploidy = 2, bool phased = false)
  {
    detail::print_vcf_site_info(out, in, info_fields);

    std::ostreambuf_iterator<char> out_it(out);

    std::uint32_t a = 0;
    for (auto it = in_data.begin(); it != in_data.end(); ++it)
    {
      if (a == 0)
        out_it = '\t';

      if (*it <= 5.0)
        out_it = '0';
      else
        out_it = '1';

      ++a;
      if (a == ploidy)
        a = 0;
      else
        out_it = (phased ? '|' : '/');
    }

    out_it = '\n';
  }
}
#endif //LIBSAVVY_SITE_INFO_HPP
