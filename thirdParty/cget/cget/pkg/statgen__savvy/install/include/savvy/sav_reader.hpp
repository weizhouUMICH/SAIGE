/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SAV_READER_HPP
#define LIBSAVVY_SAV_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "site_info.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"
#include "utility.hpp"
#include "data_format.hpp"
#include "compressed_vector.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <fstream>
#include <tuple>
#include <cmath>
#include <unordered_map>
#include <type_traits>
#include <memory>
#include <set>
#include <unordered_set>
#include <random>
#include <chrono>

namespace savvy
{
  namespace sav
  {
    namespace detail
    {
      template<std::uint8_t BitWidth>
      struct allele_decoder
      {
        static const std::uint8_t denom = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static std::tuple<T, std::uint64_t> decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value);
      };

      template<std::uint8_t BitWidth>
      struct allele_encoder
      {
        static const std::uint8_t multiplier = std::uint8_t(~(std::uint8_t(0xFF) << BitWidth)) + std::uint8_t(1);
        template <typename T>
        static void encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it);
        template <typename T>
        static std::int8_t encode(const T& allele);
      };
    }

//    namespace detail
//    {
//      template <std::uint8_t Exp>
//      struct static_base2_pow; //              : public std::integral_constant<std::uint8_t, 0> {};
//
//      template <> struct static_base2_pow<0> : public std::integral_constant<std::uint8_t, 1>   {};
//      template <> struct static_base2_pow<1> : public std::integral_constant<std::uint8_t, 2>   {};
//      template <> struct static_base2_pow<2> : public std::integral_constant<std::uint8_t, 4>   {};
//      template <> struct static_base2_pow<3> : public std::integral_constant<std::uint8_t, 8>   {};
//      template <> struct static_base2_pow<4> : public std::integral_constant<std::uint8_t, 16>  {};
//      template <> struct static_base2_pow<5> : public std::integral_constant<std::uint8_t, 32>  {};
//      template <> struct static_base2_pow<6> : public std::integral_constant<std::uint8_t, 64>  {};
//      template <> struct static_base2_pow<7> : public std::integral_constant<std::uint8_t, 128> {};
//    }

    std::vector<std::string> query_chromosomes(const std::string& file_path);

    //################################################################//
    class reader_base
    {
    public:
      reader_base(const std::string& file_path);
      reader_base(const std::string& file_path, savvy::fmt data_format);

      reader_base(reader_base&& source);
      reader_base& operator=(reader_base&& source);

      //reader(const reader&) = delete;
      //reader& operator=(const reader&) = delete;
      virtual ~reader_base() {}

//      template <typename T>
//      bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
//      {
//        read_variant_details(destination);
//        read_genotypes(destination, missing_value);
//
//        return good();
//      }

      explicit operator bool() const { return input_stream_->good(); }
      bool good() const { return input_stream_->good(); }
      bool fail() const { return input_stream_->fail(); }
      bool bad() const { return input_stream_->bad(); }
      bool eof() const { return input_stream_->eof(); }
      const std::vector<std::string>& samples() const { return sample_ids_; }
//      std::vector<std::string>::const_iterator prop_fields_begin() const { return metadata_fields_.begin(); }
//      std::vector<std::string>::const_iterator prop_fields_end() const { return metadata_fields_.end(); }

      const std::vector<std::string>& info_fields() const { return metadata_fields_; }
      const std::vector<std::pair<std::string,std::string>>& headers() const { return headers_; }
      savvy::fmt data_format() const { return file_data_format_; }
      std::uint32_t ploidy() const { return ploidy_; }
      const std::array<std::uint8_t, 16>& uuid() const { return uuid_; }

      /**
       *
       * @param subset IDs to include if they exist in file.
       * @return intersect of subset and samples IDs in file.
       */
      std::vector<std::string> subset_samples(const std::set<std::string>& subset);

      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_->tellg(); }
    protected:
      void read_variant_details(site_info& annotations)
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::eofbit); // No more markers to read.
          }
          else
          {
            std::uint64_t sz;
            if (varint_decode(in_it, end_it, sz) == end_it)
            {
              this->input_stream_->setstate(std::ios::badbit);
            }
            else
            {
              ++in_it;
              std::string chrom;
              chrom.resize(sz);
              if (sz)
                input_stream_->read(&chrom[0], sz);

              std::uint64_t locus;
              if (varint_decode(in_it, end_it, locus) == end_it)
              {
                this->input_stream_->setstate(std::ios::badbit);
              }
              else
              {
                ++in_it;
                if (varint_decode(in_it, end_it, sz) == end_it)
                {
                  this->input_stream_->setstate(std::ios::badbit);
                }
                else
                {
                  ++in_it;
                  std::string ref;
                  ref.resize(sz);
                  if (sz)
                    input_stream_->read(&ref[0], sz);

                  if (varint_decode(in_it, end_it, sz) == end_it)
                  {
                    this->input_stream_->setstate(std::ios::badbit);
                  }
                  else
                  {
                    ++in_it;
                    std::string alt;
                    alt.resize(sz);
                    if (sz)
                      input_stream_->read(&alt[0], sz);

                    std::unordered_map<std::string, std::string> props;
                    props.reserve(this->metadata_fields_.size());
                    std::string prop_val;
                    for (const std::string& key : metadata_fields_)
                    {
                      if (varint_decode(in_it, end_it, sz) == end_it)
                      {
                        this->input_stream_->setstate(std::ios::badbit);
                        break;
                      }
                      else
                      {
                        ++in_it;
                        if (sz)
                        {
                          prop_val.resize(sz);
                          input_stream_->read(&prop_val[0], sz);
                          props[key] = prop_val;
                        }
                      }
                    }

                    annotations = site_info(std::move(chrom), locus, std::move(ref), std::move(alt), std::move(props));

                    if (!this->input_stream_->good())
                      this->input_stream_->setstate(std::ios::badbit);
                  }
                }
              }
            }
          }
        }
      }

      template <std::uint8_t BitWidth>
      void discard_genotypes_impl()
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            for (std::size_t i = 0; i < sz && in_it != end_it; ++i)
            {
              std::uint8_t allele;
              std::uint64_t offset;
              in_it = prefixed_varint<BitWidth>::decode(++in_it, end_it, allele, offset);
            }

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      void discard_genotypes()
      {
        if (this->file_data_format_ == fmt::gt)
          this->discard_genotypes_impl<1>();
        else
          this->discard_genotypes_impl<7>();
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_al(site_info& annotations, T& destination)
      {
        if (good())
        {
          const auto missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN();
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            std::size_t an = samples().size() * ploidy_level;
            std::size_t ac = 0;

            if (subset_size_ != samples().size())
            {
              an = subset_size_ * ploidy_level;
              destination.resize(an);

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (BitWidth != 1)
                  {
                    allele = std::round(allele);
                    if (allele != typename T::value_type())
                      destination[subset_map_[sample_index] * ploidy_level + (total_offset % ploidy_level)] = allele;
                  }
                  else
                  {
                    destination[subset_map_[sample_index] * ploidy_level + (total_offset % ploidy_level)] = allele;
                  }

                  if (std::isnan(allele))
                    --an;
                  else if (allele)
                    ++ac;
                }
              }
            }
            else
            {
              destination.resize(an);

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                if (BitWidth != 1)
                {
                  allele = std::round(allele);
                  if (allele != typename T::value_type())
                    destination[total_offset] = allele;
                }
                else
                {
                  destination[total_offset] = allele;
                }

                if (std::isnan(allele))
                  --an;
                else if (allele)
                  ++ac;
              }
            }

            annotations.prop("AC", std::to_string(ac));
            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(static_cast<float>(ac) / static_cast<float>(an)));

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_gt(site_info& annotations, T& destination)
      {
        if (good())
        {
          const auto missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN();
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            std::size_t an = samples().size() * ploidy_level;
            std::size_t ac = 0;

            if (subset_size_ != samples().size())
            {
              destination.resize(subset_size_);

              an = subset_size_ * ploidy_level;

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (BitWidth != 1)
                  {
                    allele = std::round(allele);
                    if (allele != typename T::value_type())
                      destination[subset_map_[sample_index]] += allele;
                  }
                  else
                  {
                    destination[subset_map_[sample_index]] += allele;
                  }

                  if (std::isnan(allele))
                    --an;
                  else if (allele)
                    ++ac;
                }
              }
            }
            else
            {
              destination.resize(samples().size());

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                if (BitWidth != 1)
                {
                  allele = std::round(allele);
                  if (allele != typename T::value_type())
                    destination[total_offset / ploidy_level] += allele;
                }
                else
                {
                  destination[total_offset / ploidy_level] += allele;
                }

                if (std::isnan(allele))
                  --an;
                else if (allele)
                  ++ac;
              }
            }

            annotations.prop("AC", std::to_string(ac));
            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(static_cast<float>(ac) / static_cast<float>(an)));

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_gp(site_info& annotations, T& destination)
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            const std::size_t stride = ploidy_level + 1;

            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            {
              std::size_t num_haps = samples().size() * ploidy_level;
              destination.resize(subset_size_ * stride);

              std::vector<typename T::value_type> hap_tmp;
              hap_tmp.reserve(ploidy_level);

              auto write_gp_to_dest = [this, stride, ploidy_level](std::size_t hap_index, const std::vector<typename T::value_type>& hap_probs, T& destination)
              {
                const std::uint64_t sample_index = hap_index / ploidy_level;
                if (this->subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  typename T::value_type gp = hds_to_gp<typename T::value_type>::get_first_prob(hap_probs);
                  if (gp != typename T::value_type(0))
                    destination[this->subset_map_[sample_index] * stride] = gp;

                  if (ploidy_level == 2)
                  {
                    gp = hap_probs[0] * (typename T::value_type(1) - hap_probs[1]) + hap_probs[1] * (typename T::value_type(1) - hap_probs[0]);
                    if (gp != typename T::value_type(0))
                      destination[this->subset_map_[sample_index] * stride + 1] = gp;
                  }
                  else
                  {
                    for (std::size_t g = 1; g < ploidy_level; ++g)
                    {
                      gp = hds_to_gp<typename T::value_type>::get_prob(hap_probs, g);
                      if (gp != typename T::value_type(0))
                      {
                        destination[this->subset_map_[sample_index] * stride + g] = gp;
                      }
                    }
                  }

                  gp = hds_to_gp<typename T::value_type>::get_last_prob(hap_probs);
                  if (gp != typename T::value_type(0))
                    destination[this->subset_map_[sample_index] * stride + ploidy_level] = gp;
                }
              };

              std::size_t h = 0;
              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());
                total_offset += offset;


                assert(total_offset < num_haps);
                for ( ; h < total_offset; ++h)
                {
                  hap_tmp.push_back(typename T::value_type(0));
                  if (hap_tmp.size() == ploidy_level)
                  {
                    write_gp_to_dest(h, hap_tmp, destination);
                    hap_tmp.resize(0);
                  }
                }

                hap_tmp.push_back(allele);
                if (hap_tmp.size() == ploidy_level)
                {
                  write_gp_to_dest(h, hap_tmp, destination);
                  hap_tmp.resize(0);
                }
                ++h;
              }

              // TODO: This section can be optimized. After hap_tmp is cleared out, the rest will be 1,0,0.
              for ( ; h < num_haps; ++h)
              {
                hap_tmp.push_back(typename T::value_type(0));
                if (hap_tmp.size() == ploidy_level)
                {
                  write_gp_to_dest(h, hap_tmp, destination);
                  hap_tmp.resize(0);
                }
              }
            }

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_hds(site_info& annotations, T& destination)
      {
        if (good())
        {
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            std::size_t an = samples().size() * ploidy_level;
            float dose_sum = 0.f;

            if (subset_size_ != samples().size())
            {
              an = subset_size_ * ploidy_level;
              destination.resize(an);


              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());

                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  destination[subset_map_[sample_index] * ploidy_level + (total_offset % ploidy_level)] = allele;

                  if (std::isnan(allele))
                    --an;
                  else
                    dose_sum += allele;
                }
              }
            }
            else
            {
              destination.resize(an);

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, std::numeric_limits<typename T::value_type>::quiet_NaN());

                total_offset += offset;

                assert(total_offset < (samples().size() * ploidy_level));
                destination[total_offset] = allele;

                if (std::isnan(allele))
                  --an;
                else
                  dose_sum += allele;
              }
            }

            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(dose_sum / static_cast<float>(an)));

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      template <std::size_t BitWidth, typename T>
      void read_genotypes_ds(site_info& annotations, T& destination)
      {
        if (good())
        {
          const typename T::value_type missing_value(std::numeric_limits<typename T::value_type>::quiet_NaN());
          std::istreambuf_iterator<char> in_it(*input_stream_);
          std::istreambuf_iterator<char> end_it;

          std::uint64_t ploidy_level;
          if (ploidy_ == 0)
          {
            if (varint_decode(in_it, end_it, ploidy_level) != end_it)
              ++in_it;
          }
          else
          {
            ploidy_level = ploidy_;
          }

          if (in_it == end_it)
          {
            this->input_stream_->setstate(std::ios::badbit);
          }
          else
          {
            std::uint64_t sz;
            varint_decode(in_it, end_it, sz);
            std::uint64_t total_offset = 0;

            std::size_t an = samples().size() * ploidy_level;
            float dose_sum = 0.f;

            if (subset_size_ != samples().size())
            {
              destination.resize(subset_size_);
              an = subset_size_ * ploidy_level;

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;

                const std::uint64_t sample_index = total_offset / ploidy_level;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  destination[subset_map_[sample_index]] += allele;

                  if (std::isnan(allele))
                    --an;
                  else
                    dose_sum += allele;
                }
              }
            }
            else
            {
              destination.resize(samples().size());

              for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
              {
                typename T::value_type allele;
                std::uint64_t offset;
                std::tie(allele, offset) = detail::allele_decoder<BitWidth>::decode(++in_it, end_it, missing_value);
                total_offset += offset;
                destination[total_offset / ploidy_level] += allele;

                if (std::isnan(allele))
                  --an;
                else
                  dose_sum += allele;
              }
            }

            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(dose_sum / static_cast<float>(an)));

            if (input_stream_->get() == std::char_traits<char>::eof())
            {
              assert(!"Truncated file");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
        }
      }

      template <typename T>
      void read_genotypes(site_info& annotations, T& destination)
      {
        destination.resize(0);
        if (true) //requested_data_formats_[idx] == file_data_format_)
        {
          if (requested_data_format_ == fmt::gt)
            file_data_format_ == fmt::gt ? read_genotypes_al<1>(annotations, destination) : read_genotypes_al<7>(annotations, destination);
          else if (requested_data_format_== fmt::ac)
            file_data_format_ == fmt::gt ? read_genotypes_gt<1>(annotations, destination) : read_genotypes_gt<7>(annotations, destination);
          else if (requested_data_format_ == fmt::gp)
            file_data_format_ == fmt::gt ? read_genotypes_gp<1>(annotations, destination) : read_genotypes_gp<7>(annotations, destination);
          else if (requested_data_format_ == fmt::ds)
            file_data_format_ == fmt::gt ? read_genotypes_ds<1>(annotations, destination) : read_genotypes_ds<7>(annotations, destination);
          else if (requested_data_format_ == fmt::hds)
            file_data_format_ == fmt::gt ? read_genotypes_hds<1>(annotations, destination) : read_genotypes_hds<7>(annotations, destination);
          else
            input_stream_->setstate(std::ios::failbit);
        }
        else
        {
          discard_genotypes();
        }
      }
    private:
      void parse_header();
      void init_subset_map();
    protected:
      std::vector<std::string> sample_ids_;
      std::vector<std::uint64_t> subset_map_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> metadata_fields_;
      std::string file_path_;
      std::uint64_t subset_size_;
      std::unique_ptr<std::istream> input_stream_;
      fmt file_data_format_;
      fmt requested_data_format_;
      std::uint32_t ploidy_ = 0;
      std::array<std::uint8_t, 16> uuid_;
    };
    //################################################################//

    //################################################################//
    class reader : public reader_base
    {
    public:
      using reader_base::reader_base;

      template <typename T>
      reader& operator>>(variant<T>& destination)
      {
        return this->read(destination, destination.data());
      }

      template <typename T>
      reader& read(site_info& annotations, T& destination)
      {
        this->read_variant_details(annotations);
        this->read_genotypes(annotations, destination);
        return *this;
      }

      template <typename Pred, typename T>
      reader& read_if(Pred fn, site_info& annotations, T& destination)
      {
        while (good())
        {
          this->read_variant_details(annotations);
          this->read_genotypes(annotations, destination);
          if (fn(annotations))
            break;
        }
        return *this;
      }
    };

    class indexed_reader : public reader_base
    {
    public:
      template <typename T>
      indexed_reader(const std::string& file_path, const std::string& index_file_path, const region& reg, bounding_point bound_type, T data_format)  :
        reader_base(file_path, data_format),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(reg)),
        i_(query_.begin()),
        reg_(reg),
        bounding_type_(bound_type),
        current_offset_in_block_(0),
        total_in_block_(0)
      {
        if (!index_.good())
          this->input_stream_->setstate(std::ios::badbit);
      }

      indexed_reader(const std::string& file_path, const region& reg, savvy::fmt data_format)  :
        indexed_reader(file_path, std::string(""), reg, bounding_point::beg, data_format)
      {
      }

      indexed_reader(const std::string& file_path, const std::string& index_file_path, const region& reg, savvy::fmt data_format)  :
        indexed_reader(file_path, index_file_path, reg, bounding_point::beg, data_format)
      {
      }

      indexed_reader(const std::string& file_path, const region& reg, bounding_point bounding_type, savvy::fmt data_format)  :
        indexed_reader(file_path, std::string(""), reg, bounding_type, data_format)
      {
      }

      std::vector<std::string> chromosomes() const
      {
        return index_.tree_names();
      }

      template <typename T>
      indexed_reader& operator>>(variant<T>& destination)
      {
        return this->read(destination, destination.data());
      }

      template <typename T>
      indexed_reader& read(site_info& annotations, T& destination)
      {
        while (this->good())
        {
          if (current_offset_in_block_ >= total_in_block_)
          {
            if (i_ == query_.end())
              this->input_stream_->setstate(std::ios::eofbit);
            else
            {
              total_in_block_ = std::uint32_t(0x000000000000FFFF & i_->value()) + 1;
              current_offset_in_block_ = 0;
              this->input_stream_->seekg(std::streampos((i_->value() >> 16) & 0x0000FFFFFFFFFFFF));
              ++i_;
            }
          }

          this->read_variant_details(annotations);
          if (!this->good())
          {
            if (current_offset_in_block_ < total_in_block_)
            {
              assert(!"Truncated block");
              this->input_stream_->setstate(std::ios::badbit);
            }
          }
          else
          {
            ++current_offset_in_block_;
            if (region_compare(bounding_type_, annotations, reg_))
            {
              this->read_genotypes(annotations, destination);
              break;
            }
            else
            {
              this->discard_genotypes();
            }
          }
        }
        return *this;
      }

      template <typename Pred, typename T>
      indexed_reader& read_if(Pred fn, site_info& annotations, T& destination)
      {
        while (this->good())
        {
          read(annotations, destination);
          if (fn(annotations))
            break;
        }

        return *this;
      }

      void reset_region(const region& reg)
      {
        current_offset_in_block_ = 0;
        total_in_block_ = 0;
        reg_ = reg;
        this->input_stream_->clear();
        query_ = index_.create_query(reg);
        i_ = query_.begin();
        if (!index_.good())
          this->input_stream_->setstate(std::ios::badbit);
      }
    private:
      s1r::reader index_;
      s1r::reader::query query_;
      s1r::reader::query::iterator i_;
      region reg_; //TODO: make this a default template argument when vector type is also a reader template.
      bounding_point bounding_type_;
      std::uint32_t current_offset_in_block_;
      std::uint32_t total_in_block_;
    };
    //################################################################//

    class writer
    {
    public:
      struct options
      {
        std::int8_t compression_level;
        std::uint16_t block_size;
        std::string index_path;
        options() :
          compression_level(3),
          block_size(2048)
        {
        }
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, fmt data_format) :
        writer(file_path, options(), samples_beg, samples_end, headers_beg, headers_end, data_format)
      {
      }

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator>
      writer(const std::string& file_path, const options& opts, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, fmt data_format) :
        rng_(std::chrono::high_resolution_clock::now().time_since_epoch().count() ^ std::clock() ^ (std::uint64_t)this),
        output_buf_(create_out_streambuf(file_path, opts.compression_level)), //opts.compression == compression_type::zstd ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path)) : std::unique_ptr<std::streambuf>(new std::filebuf(file_path, std::ios::binary))),
        output_stream_(output_buf_.get()),
        samples_(samples_beg, samples_end),
        file_path_(file_path),
        uuid_(::savvy::detail::gen_uuid(rng_)),
        index_file_(opts.index_path.size() ? ::savvy::detail::make_unique<s1r::writer>(opts.index_path, uuid_) : nullptr),
        current_block_min_(std::numeric_limits<std::uint32_t>::max()),
        current_block_max_(0),
        allele_count_(0),
        record_count_(0),
        record_count_in_block_(0),
        block_size_(opts.block_size),
        data_format_(data_format)
      {
        headers_.resize(std::distance(headers_beg, headers_end));
        auto copy_res = std::copy_if(headers_beg, headers_end, headers_.begin(), [](const std::pair<std::string,std::string>& kvp) { return kvp.first != "FORMAT" && kvp.first != "fileformat"; });
        headers_.resize(std::distance(headers_.begin(), copy_res));
      }


      template <typename RandAccessStringIterator>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, options opts = options()) :
        writer(file_path, std::forward<RandAccessStringIterator>(samples_beg), std::forward<RandAccessStringIterator>(samples_end), empty_string_pair_array.end(), empty_string_pair_array.end(), opts)
      {

      }

      ~writer()
      {
        // TODO: This is only a temp solution.
        if (index_file_)
        {
          if (record_count_in_block_)
          {
            auto file_pos = std::uint64_t(output_stream_.tellp());
            if (record_count_in_block_ > 0x10000) // Max records per block: 64*1024
            {
              assert(!"Too many records in zstd frame to be indexed!");
            }

            if (file_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
            {
              assert(!"File size to large to be indexed!");
            }

            s1r::entry e(current_block_min_, current_block_max_, (file_pos << 16) | std::uint16_t(record_count_in_block_ - 1));
            index_file_->write(current_chromosome_, e);
          }
        }
      }

      void write_header(std::int32_t ploidy)
      {
        if (ploidy_ == 0 && ploidy > 0 && good())
        {
          ploidy_ = ploidy;

          std::string version_string("sav\x00\x01\x00\x00", 7);
          output_stream_.write(version_string.data(), version_string.size());

          output_stream_.write((char*)uuid_.data(), uuid_.size());

          std::ostreambuf_iterator<char> out_it(output_stream_);


          // TODO: Handle unsupported formats.
          std::string fmt_str;
          if (data_format_ == fmt::hds)
            fmt_str = "<ID=HDS,Type=Float,Number=" + std::to_string(ploidy_) + ",Description=\"Haplotype dosages\">";
          else
            fmt_str = "<ID=GT,Type=Integer,Number=" + std::to_string(ploidy_) + ",Description=\"Genotype\">";
          headers_.push_back(std::make_pair(std::string("FORMAT"), fmt_str));

          std::unordered_set<std::string> unique_info_fields;

          varint_encode(headers_.size(), out_it);
          for (auto it = headers_.begin(); it != headers_.end(); ++it)
          {
            std::size_t str_sz = get_string_size(it->first);
            varint_encode(str_sz, out_it);
            if (str_sz)
            {
              output_stream_.write(it->first.data(), str_sz);

              str_sz = get_string_size(it->second);
              varint_encode(str_sz, out_it);
              if (str_sz)
                output_stream_.write(it->second.data(), str_sz);
            }

            if (it->first == "INFO")
            {
              std::string info_id = parse_header_sub_field(it->second, "ID");
              if (unique_info_fields.emplace(info_id).second)
                this->property_fields_.emplace_back(info_id);
            }
          }

          varint_encode(samples_.size(), out_it);
          for (auto it = samples_.begin(); it != samples_.end(); ++it)
          {
            std::size_t str_sz = it->size();
            varint_encode(str_sz, out_it);
            if (str_sz)
              output_stream_.write(&(*it)[0], str_sz);
          }
        }
      }

      template <typename T>
      writer& operator<<(const savvy::variant<T>& v)
      {
        write(v, v.data());
        return *this;
      }

//#define NO_LEB128 1
#ifdef NO_LEB128
      template <typename T>
      void write(const allele_vector<T>& m)
      {
        const typename T::value_type ref_value = typename T::value_type();
        //std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
        std::uint64_t sz = m.locus();

        output_stream_.write((char*)&sz, 8); //varint_encode(m.locus(), os_it);

        sz = (m.ref().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          output_stream_.write(m.ref().data(), m.ref().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        sz = (m.alt().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.alt().size())
          output_stream_.write(m.alt().data(), m.alt().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        for (const std::string& key : property_fields_)
        {
          std::string value(m.prop(key));
          sz = (m.ref().size() << 48);
          output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
          if (value.size())
            output_stream_.write(value.data(), value.size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
          //os.write(&source.ref_[0], source.ref_.size());
        }

        struct sparse_geno
        {
          std::uint32_t v: 1, offset: 31;
        };

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        output_stream_.write((char*)&non_zero_count, 8);//varint_encode(non_zero_count, os_it);

        std::vector<sparse_geno> tmp(non_zero_count);

        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        std::size_t non_ref_counter = 0;
        for (auto it = beg; it != m.end(); ++it)
        {
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            tmp[non_ref_counter].v = (std::isnan(*it)  ? std::uint8_t(0x80) : std::uint8_t(0x00));
            tmp[non_ref_counter].offset = offset;
            ++non_ref_counter;
          }
        }
        output_stream_.write((char*)tmp.data(), tmp.size() * 4);
      }
#else
      template <typename VecT>
      void write(const site_info& annotations, const VecT& data)
      {
        if (this->good())
        {
          if (data.size() % samples_.size() != 0)
          {
            output_stream_.setstate(std::ios::failbit);
          }
          else
          {
            if (ploidy_ == 0)
            {
              write_header(std::uint32_t((data.size() / samples_.size()) & 0xFFFFFFFF));
            }

            if (data.size() / samples_.size() != ploidy_)
            {
              output_stream_.setstate(std::ios::failbit);
            }
            else
            {
              // 1024*1024 non-ref GTs or 64*1024 records
              //if (allele_count_ >= 0x100000 || (record_count_ % 0x10000) == 0 || annotations.chromosome() != current_chromosome_)
              if (block_size_ != 0 && ((record_count_ % block_size_) == 0 || annotations.chromosome() != current_chromosome_))
              {
                if (index_file_ && record_count_in_block_)
                {
                  auto file_pos = std::uint64_t(output_stream_.tellp());
                  if (record_count_in_block_ > 0x10000) // Max records per block: 64*1024
                  {
                    assert(!"Too many records in zstd frame to be indexed!");
                    output_stream_.setstate(std::ios::badbit);
                  }

                  if (file_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
                  {
                    assert(!"File size to large to be indexed!");
                    output_stream_.setstate(std::ios::badbit);
                  }

                  s1r::entry e(current_block_min_, current_block_max_, (file_pos << 16) | std::uint16_t(record_count_in_block_ - 1));
                  index_file_->write(current_chromosome_, e);
                }
                output_stream_.flush();
                allele_count_ = 0;
                current_chromosome_ = annotations.chromosome();
                record_count_in_block_ = 0;
                current_block_min_ = std::numeric_limits<std::uint32_t>::max();
                current_block_max_ = 0;
              }

              std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

              varint_encode(annotations.chromosome().size(), os_it);
              std::copy(annotations.chromosome().begin(), annotations.chromosome().end(), os_it);

              varint_encode(annotations.position(), os_it);

              varint_encode(annotations.ref().size(), os_it);
              if (annotations.ref().size())
                std::copy(annotations.ref().begin(), annotations.ref().end(), os_it);
              //os.write(&source.ref_[0], source.ref_.size());

              varint_encode(annotations.alt().size(), os_it);
              if (annotations.alt().size())
                std::copy(annotations.alt().begin(), annotations.alt().end(), os_it);
              //os.write(&source.alt_[0], source.alt_.size());

              for (const std::string& key : property_fields_)
              {
                std::string value(annotations.prop(key));
                varint_encode(value.size(), os_it);
                if (value.size())
                  std::copy(value.begin(), value.end(), os_it);
              }

              if (data_format_ == fmt::hds)
              {
                write_hap_dosages(data);
              }
//            else if (data_format_ == fmt::genotype_probability)
//            {
//              write_probs(data);
//            }
              else
              {
                write_alleles(data);
              }

              current_block_min_ = std::min(current_block_min_, std::uint32_t(annotations.position()));
              current_block_max_ = std::max(current_block_max_, std::uint32_t(annotations.position() + std::max(annotations.ref().size(), annotations.alt().size())) - 1);
              ++record_count_in_block_;
              ++record_count_;
            }
          }
        }
      }
#endif
      explicit operator bool() const { return good(); }
      bool good() const { return output_stream_.good() && (!index_file_ || index_file_->good()); }
      bool fail() const { return output_stream_.fail(); }
      bool bad() const { return output_stream_.bad(); }
      bool eof() const { return output_stream_.eof(); }

      static bool create_index(const std::string& input_file_path, std::string output_file_path = "");
    protected:
      template <typename T>
      static std::size_t get_string_size(T str);

      static std::filebuf* create_std_filebuf(const std::string& file_path, std::ios::openmode mode)
      {
        std::filebuf* ret = new std::filebuf();
        ret->open(file_path.c_str(), mode);
        return ret;
      }

      static std::unique_ptr<std::streambuf> create_out_streambuf(const std::string& file_path, std::int8_t compression_level);

      template <std::size_t BitWidth, typename T, typename OutIt>
      static void serialize_alleles(const std::vector<T>& m, OutIt os_it)
      {
        std::uint64_t last_pos = 0;
        const auto beg = m.begin();
        for (auto it = beg; it != m.end(); ++it)
        {
          std::int8_t signed_allele = detail::allele_encoder<BitWidth>::encode(*it);
          if (signed_allele >= 0)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            prefixed_varint<BitWidth>::encode((std::uint8_t)(signed_allele), offset, os_it);
          }
        }

      }

      template <std::size_t BitWidth, typename T, typename OutIt>
      static void serialize_alleles(const savvy::compressed_vector<T>& m, OutIt os_it)
      {
        std::uint64_t last_pos = 0;
        auto end = m.end();
        for (auto it = m.begin(); it != end; ++it)
        {
          std::int8_t signed_allele = detail::allele_encoder<BitWidth>::encode(*it);
          if (signed_allele >= 0)
          {
            std::uint64_t dist = it.offset();
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            prefixed_varint<BitWidth>::encode((std::uint8_t)(signed_allele), offset, os_it);
          }
        }

      }

      template <typename T>
      void write_alleles(const std::vector<T>& m)
      {
        const T ref_value = T();

        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<1>(m, os_it);
      }

      template <typename T>
      void write_alleles(const savvy::compressed_vector<T>& m)
      {
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        allele_count_ += m.non_zero_size();
        varint_encode(m.non_zero_size(), os_it);

        serialize_alleles<1>(m, os_it);
      }

      template <typename T>
      void write_hap_dosages(const std::vector<T>& m)
      {
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint64_t non_zero_count = 0;
        for (auto it = m.begin(); it != m.end(); ++it)
        {
          if (detail::allele_encoder<7>::encode(*it) >= 0)
            ++non_zero_count;
        }

        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<7>(m, os_it);
      }

      template <typename T>
      void write_hap_dosages(const savvy::compressed_vector<T>& m)
      {
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());

        std::uint64_t non_zero_count = 0;
        for (auto it = m.begin(); it != m.end(); ++it)
        {
          if (detail::allele_encoder<7>::encode(*it) >= 0)
            ++non_zero_count;
        }

        allele_count_ += non_zero_count;
        varint_encode(non_zero_count, os_it);

        serialize_alleles<7>(m, os_it);
      }

//      template <typename T>
//      void write_probs(const std::vector<T>& m)
//      {
//        const T ref_value = T();
//
//        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
//
//        std::uint32_t ploidy = std::uint32_t((m.size() / sample_size_) & 0xFFFFFFFF) - 1;
//        std::uint32_t stride = ploidy + 1;
//
//        // TODO: check modulus and set error if needed.
//        varint_encode(ploidy, os_it);
//
//        auto beg = m.begin();
//        std::uint64_t non_zero_count = 0;
//        std::size_t c = 0;
//        for (auto it = m.begin(); it != m.end(); ++it,++c)
//        {
//          if (c % stride != 0)
//          {
//            if (allele_encoder<7>::encode(*it) >= 0)
//              ++non_zero_count;
//          }
//        }
//
//
//        allele_count_ += non_zero_count;
//        varint_encode(non_zero_count, os_it);
//        std::uint64_t last_pos = 0;
//        c = 0;
//        for (auto it = beg; it != m.end(); ++it,++c)
//        {
//          if (c % stride != 0)
//          {
//            //std::int8_t signed_allele = std::round((std::isnan(*it) ? T::value_type(0.5) : *it) * type_multiplier) - T::value_type(1);
//            std::int8_t signed_allele = allele_encoder<7>::encode(*it);
//            if (signed_allele >= 0)
//            {
//              std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
//              std::uint64_t offset = dist - last_pos;
//              last_pos = dist + 1;
//              prefixed_varint<7>::encode((std::uint8_t)(signed_allele), offset, os_it);
//            }
//          }
//        }
//      }
    private:
      static const std::array<std::string, 0> empty_string_array;
      static const std::array<std::pair<std::string, std::string>, 0> empty_string_pair_array;
    protected:
      std::mt19937_64 rng_;
      std::unique_ptr<std::streambuf> output_buf_;
      std::ostream output_stream_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> property_fields_;
      std::vector<std::string> samples_;
      std::string file_path_;
      std::array<std::uint8_t, 16> uuid_;
      std::unique_ptr<s1r::writer> index_file_;
      std::string current_chromosome_;
      std::uint32_t current_block_min_;
      std::uint32_t current_block_max_;
      std::uint32_t metadata_fields_cnt_;
      std::size_t allele_count_;
      std::size_t record_count_;
      std::size_t record_count_in_block_;
      std::uint16_t block_size_;
      fmt data_format_;
      std::int32_t ploidy_ = 0;
    };


    template <>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<0>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret{T(1), 0};
      in_it = varint_decode(in_it, end_it, std::get<1>(ret));
      return ret;
    }

    template<>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<1>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<1>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (allele ? T(1) : missing_value);
      return ret;
    }

    template<std::uint8_t BitWidth>
    template <typename T>
    inline std::tuple<T, std::uint64_t> detail::allele_decoder<BitWidth>::decode(std::istreambuf_iterator<char>& in_it, const std::istreambuf_iterator<char>& end_it, const T& missing_value)
    {
      std::tuple<T, std::uint64_t> ret;
      std::uint8_t allele;
      in_it = prefixed_varint<BitWidth>::decode(in_it, end_it, allele, std::get<1>(ret));
      std::get<0>(ret) = (static_cast<T>(allele) + T(1)) / denom;
      return ret;
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<0>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      varint_encode(offset, os_it);
    }

    template<>
    template <typename T>
    inline void detail::allele_encoder<1>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<1>::encode(std::uint8_t(std::isnan(allele) ? 0 : 1), offset, os_it);
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline void detail::allele_encoder<ByteWidth>::encode(const T& allele, std::uint64_t offset, std::ostreambuf_iterator<char>& os_it)
    {
      prefixed_varint<ByteWidth>::encode(std::uint8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1)), offset, os_it);
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<0>::encode(const T& allele)
    {
      return -1;
    }

    template<>
    template <typename T>
    inline std::int8_t detail::allele_encoder<1>::encode(const T& allele)
    {
      return std::int8_t(std::isnan(allele) ? 0 : (allele == T() ? -1 : 1));
    }

    template<std::uint8_t ByteWidth>
    template <typename T>
    inline std::int8_t detail::allele_encoder<ByteWidth>::encode(const T& allele)
    {
      return std::int8_t(std::round((std::isnan(allele) ? T(0.5) : allele) * multiplier) - T(1));
    }


    template <typename T>
    std::size_t writer::get_string_size(T str)
    {
      return str.size();
    }

    template <>
    inline std::size_t writer::get_string_size<const char*>(const char* str)
    {
      return std::strlen(str);
    }
  }
}

#endif //LIBSAVVY_SAV_READER_HPP
