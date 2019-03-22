/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_VCF_READER_HPP
#define LIBSAVVY_VCF_READER_HPP

#include "allele_status.hpp"
#include "site_info.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"
#include "utility.hpp"
#include "data_format.hpp"
#include "savvy.hpp"

#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <ctime>
#include <cassert>

namespace savvy
{
  namespace vcf
  {
    //################################################################//
    enum class compression_type : std::uint8_t
    {
      none = 0,
      bgzip
    };

    enum class empty_vector_policy : std::uint8_t
    {
      fail = 0,
      skip,
      skip_with_warning
    };
    //################################################################//

    std::vector<std::string> query_chromosomes(const std::string& file_path);

    namespace detail
    {
      enum class hts_info_type
      {
        bcf_bt_null  = 0,
        bcf_bt_int8  = 1,
        bcf_bt_int16 = 2,
        bcf_bt_int32 = 3,
        bcf_bt_float = 5,
        bcf_bt_char  = 7
      };

      class hts_file_base
      {
      public:
        static std::unique_ptr<hts_file_base> create_file(const std::string& file_path);
        static std::unique_ptr<hts_file_base> create_indexed_file(const std::string& file_path, const region& reg);

        virtual ~hts_file_base() {}

        virtual void init_headers(std::vector<std::pair<std::string,std::string>>& destination) = 0;
        virtual void init_sample_ids(std::vector<std::string>& destination) = 0;
        virtual void init_info_fields(std::vector<std::string>& destination) = 0;
        virtual std::size_t cur_fmt_field_size() const = 0;
        virtual const char*const cur_fmt_field(std::size_t idx) const = 0;
        virtual std::size_t cur_num_alleles() const = 0;
        virtual bool read_next_record() = 0;
        virtual site_info cur_site_info(std::size_t allele_index) const = 0;
        virtual bool get_cur_format_values_int32(const char* tag, int**buf, int*sz) const = 0;
        virtual bool get_cur_format_values_float(const char* tag, int**buf, int*sz) const = 0;
      };

      std::unique_ptr<std::ostream> create_out_stream(const std::string& file_path, compression_type type);
    }

    //################################################################//
    template <std::size_t VecCnt>
    class reader_base
    {
    public:
      reader_base() :
        subset_size_(0),
        state_(std::ios::goodbit),
        gt_(nullptr),
        gt_sz_(0),
        allele_index_(0)
      {}

      reader_base(reader_base&& source) :
        headers_(std::move(source.headers_)),
        subset_map_(std::move(source.subset_map_)),
        hts_file_(std::move(source.hts_file_)),
        subset_size_(source.subset_size_),
        state_(source.state_),
        empty_vector_policy_(source.empty_vector_policy_),
        gt_(source.gt_),
        gt_sz_(source.gt_sz_),
        allele_index_(source.allele_index_),
        property_fields_(std::move(source.property_fields_)),
        requested_data_formats_(source.requested_data_formats_)
      {
        source.gt_ = nullptr;
      }

      virtual ~reader_base()
      {
        if (gt_)
          free(gt_);
      }

      //virtual reader_base& operator>>(block& destination) = 0;

      explicit operator bool() const { return good(); }
      bool good() const { return state_ == std::ios::goodbit; }
      bool fail() const { return (state_ & std::ios::failbit) != 0; }
      bool bad() const { return (state_ & std::ios::badbit) != 0; }
      bool eof() const { return (state_ & std::ios::eofbit) != 0; }

      std::vector<std::string> subset_samples(const std::set<std::string>& subset);

      const std::vector<std::string>& samples() const;
      const std::vector<std::string>& info_fields() const;
      const std::vector<std::pair<std::string, std::string>>& headers() const { return headers_; }
      void set_policy(enum empty_vector_policy policy) { empty_vector_policy_ = policy; }
    protected:
      static const int bcf_gt_missing = 0;
      void init_sample_ids();
      void init_property_fields();
      void init_headers();

      template <typename T>
      void clear_destinations(T& destination);
      template <typename T, typename... T2>
      void clear_destinations(T& destination, T2&... other_destinations);

      void read_variant_details(site_info& destination);

      template <typename... T>
      std::size_t read_requested_genos(site_info& annotations, T&... vec);
      template <std::size_t Idx, typename T1>
      bool read_genos_to(fmt data_format, site_info& annotations, T1& vec);
      template <std::size_t Idx, typename T1, typename... T2>
      bool read_genos_to(fmt data_format, site_info& annotations, T1& vec, T2&... others);

      template <typename T>
      void read_genotypes_al(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_gt(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_gp(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_ds(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_hds(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_gl(site_info& annotations, T& destination);
      template <typename T>
      void read_genotypes_pl(site_info& annotations, T& destination);

      void init_requested_formats(fmt f);
      template <typename... T2>
      void init_requested_formats(fmt f, T2... args);
    protected:
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> sample_ids_;
      std::vector<std::string> property_fields_;
      std::array<fmt, VecCnt> requested_data_formats_;
      std::vector<std::uint64_t> subset_map_;
      std::unique_ptr<detail::hts_file_base> hts_file_;
      std::uint64_t subset_size_;
      empty_vector_policy empty_vector_policy_ = empty_vector_policy::fail;
      std::ios::iostate state_;
      int* gt_;
      int gt_sz_;
      int allele_index_;
    };
    //################################################################//

    //################################################################//
    template <std::size_t VecCnt>
    class reader : public reader_base<VecCnt>
    {
    public:
      //typedef reader_base::input_iterator input_iterator;
      template <typename... T>
      reader(const std::string& file_path, T... data_formats);
      reader(reader&& other);
      reader(const reader&) = delete;
      reader& operator=(const reader&) = delete;
      ~reader();

      template <typename T>
      reader& operator>>(variant<T>& destination)
      {
        //static_assert(VecCnt == 1, "Extraction operator only supported with single format field reader");
        return this->read(destination, destination.data());
      }

      template <typename... T>
      reader& read(site_info& annotations, T&... destinations);

//      static std::string get_chromosome(const reader& rdr, const marker& mkr);
    private:
//      template <typename T>
//      bool read_block(allele_vector<T>& destination)
//      {
//        bool ret = true;
//
//        ++allele_index_;
//        if (allele_index_ >= hts_rec_->n_allele)
//        {
//          ret = read_hts_record();
//        }
//
//        if (ret)
//        {
//          destination = allele_vector<T>(
//            std::string(bcf_hdr_id2name(hts_hdr_, hts_rec_->rid)),
//            static_cast<std::uint64_t>(hts_rec_->pos + 1),
//            std::string(hts_rec_->d.allele[0]),
//            std::string(hts_rec_->n_allele > 1 ? hts_rec_->d.allele[allele_index_] : ""),
//            sample_count(),
//            (gt_sz_ / hts_rec_->n_sample),
//            std::move(destination));
//
//          for (std::size_t i = 0; i < gt_sz_; ++i)
//          {
//            if (gt_[i] == bcf_gt_missing)
//              destination[i] = std::numeric_limits<typename T::value_type>::quiet_NaN();
//            else if (bcf_gt_allele(gt_[i]) == allele_index_)
//              destination[i] = 1.0;
//          }
//        }
//
//        return ret;
//      }

//      bcf_hdr_t* hts_hdr() const { return hts_hdr_; }
//      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
//      htsFile* hts_file_;
//      bcf_hdr_t* hts_hdr_;
//      bcf1_t* hts_rec_;
    };
    //################################################################//

    //################################################################//
    template <std::size_t VecCnt>
    class indexed_reader : public reader_base<VecCnt>
    {
    public:
      template <typename... T>
      indexed_reader(const std::string& file_path, const region& reg, T... data_formats);
      template <typename... T>
      indexed_reader(const std::string& file_path, const region& reg, bounding_point bounding_type, T... data_formats);
      //template <typename PathItr, typename RegionItr>
      //region_reader(PathItr file_paths_beg, PathItr file_paths_end, RegionItr regions_beg, RegionItr regions_end);
      ~indexed_reader();
      void reset_region(const region& reg);

      std::vector<std::string> chromosomes() const;

      template <typename T>
      indexed_reader& operator>>(variant<T>& destination);
      template <typename... T>
      indexed_reader<VecCnt>& read(site_info& annotations, T&... destinations);

      template <typename Pred, typename... T>
      indexed_reader<VecCnt>& read_if(Pred fn, site_info& annotations, T&... destinations);
    private:
//      index_reader& seek(const std::string& chromosome, std::uint64_t position);
//
//      static std::string get_chromosome(const index_reader& rdr, const marker& mkr);
//      bcf_hdr_t* hts_hdr() const { return bcf_sr_get_header(synced_readers_, 0); }
//      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      region region_;
      std::string file_path_;
//      bcf_srs_t* synced_readers_;
//      bcf1_t* hts_rec_;
      bounding_point bounding_type_;
    };
    //################################################################//

    //################################################################//
    template <std::size_t VecCnt>
    class writer
    {
    public:
      struct options
      {
        compression_type compression;
        options() :
          compression(compression_type::none)
        {}
      };

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator, typename... Fmt>
      writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, Fmt... data_formats);

      template <typename RandAccessStringIterator, typename RandAccessKVPIterator, typename... Fmt>
      writer(const std::string& file_path, const options& opts, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, Fmt... data_formats);

      void init_format_fields(savvy::fmt f);

      template <typename... Fmt>
      void init_format_fields(savvy::fmt f, Fmt... other);

      template <typename... T>
      writer& write(const site_info& anno, const T&... data);

      template <typename T>
      writer& operator<<(const variant<std::vector<T>>& v);

      bool good() { return output_stream_->good(); }
    private:
      template <typename T>
      static const std::vector<T>& get_vec(std::size_t offset, const std::vector<T>& m);

      template <typename T, typename... T2>
      static const std::vector<T>& get_vec(std::size_t offset, const std::vector<T>& m, const T2&... other);

      template <typename... T>
      void write_multi_sample_level_data(const std::size_t ploidy, const T&... data);
    private:
      std::vector<std::string> info_fields_;
      std::vector<int> info_field_types_;
      std::vector<fmt> format_fields_;
      std::unique_ptr<std::ostream> output_stream_;
      std::size_t sample_size_;
      char phase_character_ = '|';
    };
    //################################################################//






    //################################################################//
    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_headers()
    {
      if (hts_file_)
      {
        hts_file_->init_headers(headers_);
      }
    }

    template <std::size_t VecCnt>
    std::vector<std::string> reader_base<VecCnt>::subset_samples(const std::set<std::string>& subset)
    {
      std::vector<std::string> ret;
//      const char** beg = hts_hdr() ? (const char**) hts_hdr()->samples : nullptr;
//      const char** end = hts_hdr() ? (const char**) (hts_hdr()->samples + bcf_hdr_nsamples(hts_hdr())) : nullptr;
      ret.reserve(std::min(subset.size(), (std::size_t)std::distance(sample_ids_.begin(), sample_ids_.end())));

      subset_map_.clear();
      subset_map_.resize((std::size_t)std::distance(sample_ids_.begin(), sample_ids_.end()), std::numeric_limits<std::uint64_t>::max());
      std::uint64_t subset_index = 0;
      for (auto it = sample_ids_.begin(); it != sample_ids_.end(); ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map_[std::distance(sample_ids_.begin(), it)] = subset_index;
          ret.push_back(*it);
          ++subset_index;
        }
      }

      subset_size_ = subset_index;

      return ret;
    }

    template <std::size_t VecCnt>
    const std::vector<std::string>& reader_base<VecCnt>::samples() const
    {
      return sample_ids_;
    }

    template <std::size_t VecCnt>
    const std::vector<std::string>& reader_base<VecCnt>::info_fields() const
    {
      return property_fields_;
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_sample_ids()
    {
      if (hts_file_)
        hts_file_->init_sample_ids(sample_ids_);
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_property_fields()
    {
      if (hts_file_)
        hts_file_->init_info_fields(property_fields_);
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_requested_formats(fmt f)
    {
      requested_data_formats_[requested_data_formats_.size() - 1] = f;
    }

    template <std::size_t VecCnt>
    template <typename... T>
    void reader_base<VecCnt>::init_requested_formats(fmt f, T... args)
    {
      requested_data_formats_[requested_data_formats_.size() - (sizeof...(T) + 1)] = f;
      init_requested_formats(args...);
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::clear_destinations(T& destination)
    {
      destination.resize(0);
    }

    template <std::size_t VecCnt>
    template <typename T, typename... T2>
    void reader_base<VecCnt>::clear_destinations(T& destination, T2&... other_destinations)
    {
      clear_destinations(destination);
      clear_destinations(other_destinations...);
    }

    template <std::size_t VecCnt>
    template <typename... T>
    std::size_t reader_base<VecCnt>::read_requested_genos(savvy::site_info& annotations, T&... destinations)
    {
      std::size_t cnt = 0;
      clear_destinations(destinations...);

      for (int i = 0; i < hts_file_->cur_fmt_field_size(); ++i)
      {
        std::string fmt_key = hts_file_->cur_fmt_field(i);

        std::int64_t gt_idx     = std::distance(requested_data_formats_.begin(), std::find(requested_data_formats_.begin(), requested_data_formats_.end(), fmt::ac));
        std::int64_t allele_idx = std::distance(requested_data_formats_.begin(), std::find(requested_data_formats_.begin(), requested_data_formats_.end(), fmt::gt));

        if (fmt_key == "GT" && (gt_idx < VecCnt || allele_idx < VecCnt))
        {
          if (gt_idx < VecCnt)
          {
            cnt += read_genos_to<0>(fmt::ac, annotations, destinations...);
          }
          else // allele_idx < VecCnt
          {
            cnt += read_genos_to<0>(fmt::gt, annotations, destinations...);
          }
        }
        else if (fmt_key == "DS")
        {
          cnt += read_genos_to<0>(fmt::ds, annotations, destinations...);
        }
        else if (fmt_key == "HDS")
        {
          cnt += read_genos_to<0>(fmt::hds, annotations, destinations...);
        }
        else if (fmt_key == "GP")
        {
          cnt += read_genos_to<0>(fmt::gp, annotations, destinations...);
        }
        else if (fmt_key == "GL")
        {
          cnt += read_genos_to<0>(fmt::gl, annotations, destinations...);
        }
        else if (fmt_key == "PL")
        {
          cnt += read_genos_to<0>(fmt::pl, annotations, destinations...);
        }
        else
        {
          // Discard
        }
      }

      if (cnt != VecCnt && empty_vector_policy_ != empty_vector_policy::skip)
      {
        std::cerr << "Variant is missing requested data format type." << std::endl;
        if (empty_vector_policy_ == empty_vector_policy::fail)
          state_ = std::ios::badbit | state_;
      }

      return cnt;
    }

    template <std::size_t VecCnt>
    template <std::size_t Idx, typename T1>
    bool reader_base<VecCnt>::read_genos_to(fmt data_format, site_info& annotations, T1& destination)
    {
      bool ret = true;
      if (requested_data_formats_[Idx] == data_format)
      {
        switch (data_format)
        {
          case fmt::gt:
            read_genotypes_al(annotations, destination);
            break;
          case fmt::ac:
            read_genotypes_gt(annotations, destination);
            break;
          case fmt::gp:
            read_genotypes_gp(annotations, destination);
            break;
          case fmt::ds:
            read_genotypes_ds(annotations, destination);
            break;
          case fmt::hds:
            read_genotypes_hds(annotations, destination);
            break;
          case fmt::gl:
            read_genotypes_gl(annotations, destination);
            break;
          case fmt::pl:
            read_genotypes_pl(annotations, destination);
            break;
        }
      }
      else
      {
        // Discard Genotypes
        ret = false;
      }
      return ret;
    }

    template <std::size_t VecCnt>
    template <std::size_t Idx, typename T1, typename... T2>
    bool reader_base<VecCnt>::read_genos_to(fmt data_format, site_info& annotations, T1& vec, T2&... others)
    {
      if (requested_data_formats_[Idx] == data_format)
      {
        return read_genos_to<Idx>(data_format, annotations, vec);
      }
      else
      {
        return read_genos_to<Idx + 1>(data_format, annotations, others...);
      }
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::read_variant_details(site_info& destination)
    {
      if (good())
      {
        bool res = true;
        ++allele_index_;
        if (allele_index_ >= hts_file_->cur_num_alleles())
        {
          res = hts_file_->read_next_record();


          this->allele_index_ = 1;
        }
        
        if (res)
        {
          destination = hts_file_->cur_site_info(allele_index_);
        }

        if (!res)
          this->state_ = std::ios::failbit | std::ios::eofbit; // TODO: Find better way to deal with failed reads vs. eof.

      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_al(site_info& annotations, T& destination)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        if (allele_index_ > 1 || hts_file_->get_cur_format_values_int32("GT", &(gt_), &(gt_sz_)))
        {
          if (gt_sz_ % samples().size() != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const int allele_index_plus_one = allele_index_ + 1;
            const std::uint64_t ploidy(gt_sz_ / samples().size());

            std::size_t an = samples().size() * ploidy;
            std::size_t ac = 0;

            if (subset_map_.size())
            {
              an = subset_size_ * ploidy;
              destination.resize(an);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (gt_[i] == bcf_gt_missing)
                  {
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = std::numeric_limits<typename T::value_type>::quiet_NaN();
                    --an;
                  }
                  else if ((gt_[i] >> 1) == allele_index_plus_one)
                  {
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = alt_value;
                    ++ac;
                  }
                }
              }
            }
            else
            {
              destination.resize(an);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (gt_[i] == bcf_gt_missing)
                {
                  destination[i] = std::numeric_limits<typename T::value_type>::quiet_NaN();
                  --an;
                }
                else if ((gt_[i] >> 1) == allele_index_plus_one)
                {
                  destination[i] = alt_value;
                  ++ac;
                }
              }
            }

            annotations.prop("AC", std::to_string(ac));
            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(static_cast<float>(ac) / static_cast<float>(an)));

            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_gt(site_info& annotations, T& destination)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        if (allele_index_ > 1 || hts_file_->get_cur_format_values_int32("GT", &(gt_), &(gt_sz_)))
        {
          if (gt_sz_ % samples().size() != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / samples().size());
            const int allele_index_plus_one = allele_index_ + 1;

            std::size_t an = samples().size() * ploidy;
            std::size_t ac = 0;

            if (subset_map_.size())
            {
              destination.resize(subset_size_);

              an = subset_size_ * ploidy;

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (gt_[i] == bcf_gt_missing)
                  {
                    destination[subset_map_[sample_index]] += std::numeric_limits<typename T::value_type>::quiet_NaN();
                    --an;
                  }
                  else if ((gt_[i] >> 1) == allele_index_plus_one)
                  {
                    destination[subset_map_[sample_index]] += alt_value;
                    ++ac;
                  }
                }
              }
            }
            else
            {
              destination.resize(samples().size());

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (gt_[i] == bcf_gt_missing)
                {
                  destination[i / ploidy] += std::numeric_limits<typename T::value_type>::quiet_NaN();
                  --an;
                }
                else if ((gt_[i] >> 1) == allele_index_plus_one)
                {
                  destination[i / ploidy] += alt_value;
                  ++ac;
                }
              }
            }

            annotations.prop("AC", std::to_string(ac));
            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(static_cast<float>(ac) / static_cast<float>(an)));

            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_ds(site_info& annotations, T& destination)
    {
      if (good())
      {
        if (hts_file_->cur_num_alleles() > 2)
        {
          std::cerr << "multi-allelic GP not supported" << std::endl;
          state_ = std::ios::badbit;
          return;
        }


        if (hts_file_->get_cur_format_values_float("DS", &(gt_), &(gt_sz_)))
        {
          float* ds = (float*) (void*) (gt_);
          const std::size_t num_samples = sample_ids_.size();
          if (gt_sz_ % num_samples != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / num_samples);
            const typename T::value_type zero_value{0};

            std::size_t an = samples().size() * ploidy;
            float dose_sum = 0.f;

            if (subset_map_.size())
            {
              destination.resize(subset_size_);
              an = subset_size_ * ploidy;

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  float cur_ds = ds[i];
                  if (cur_ds != zero_value)
                  {
                    destination[subset_map_[sample_index]] = cur_ds;
                    if (std::isnan(cur_ds))
                      an -= ploidy;
                    else
                      dose_sum += cur_ds;
                  }
                }
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                float cur_ds = ds[i];
                if (cur_ds != zero_value)
                {
                  destination[i] = cur_ds;
                  if (std::isnan(cur_ds))
                    an -= ploidy;
                  else
                    dose_sum += cur_ds;
                }
              }
            }

            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(dose_sum / static_cast<float>(an)));

            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_hds(site_info& annotations, T& destination)
    {
      if (good())
      {
        if (hts_file_->cur_num_alleles() > 2)
        {
          std::cerr << "multi allelic HDS not supported" << std::endl;
          state_ = std::ios::badbit;
          return;
        }


        if (hts_file_->get_cur_format_values_float("HDS", &(gt_), &(gt_sz_)))
        {
          float *hds = (float *) (void *) (gt_);
          const std::size_t num_samples = sample_ids_.size();
          if (gt_sz_ % num_samples != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / num_samples);
            const typename T::value_type zero_value{0};

            std::size_t an = samples().size() * ploidy;
            float dose_sum = 0.f;

            if (subset_map_.size())
            {
              an = subset_size_ * ploidy;
              destination.resize(an);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  float cur_hds = hds[i];
                  if (cur_hds != zero_value)
                  {
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = cur_hds;
                    if (std::isnan(cur_hds))
                      --an;
                    else
                      dose_sum += cur_hds;
                  }
                }
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                float cur_hds = hds[i];
                if (cur_hds != zero_value)
                {
                  destination[i] = cur_hds;
                  if (std::isnan(cur_hds))
                    --an;
                  else
                    dose_sum += cur_hds;
                }
              }
            }

            annotations.prop("AN", std::to_string(an));
            if (an)
              annotations.prop("AF", std::to_string(dose_sum / static_cast<float>(an)));

            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_gp(site_info& annotations, T& destination)
    {
      if (good())
      {
        if (hts_file_->cur_num_alleles() > 2)
        {
          std::cerr << "multi allelic GP not supported" << std::endl;
          state_ = std::ios::badbit;
          return;
        }


        if (hts_file_->get_cur_format_values_float("GP", &(gt_), &(gt_sz_)))
        {
          float* gp = (float*) (void*) (gt_);
          const std::size_t num_samples = sample_ids_.size();
          if (gt_sz_ % num_samples != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / num_samples - 1);
            const std::uint64_t ploidy_plus_one = ploidy + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_plus_one);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy_plus_one;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index] * ploidy_plus_one + (i % ploidy_plus_one)] = gp[i];
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                destination[i] = gp[i];
              }
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_gl(site_info& annotations, T& destination)
    {
      if (good())
      {
        if (hts_file_->cur_num_alleles() > 2)
        {
          std::cerr << "multi allelic GL not supported" << std::endl;
          state_ = std::ios::badbit;
          return;
        }


        if (hts_file_->get_cur_format_values_float("GL", &(gt_), &(gt_sz_)))
        {
          float* gl = (float*) (void*) (gt_);
          const std::size_t num_samples = sample_ids_.size();
          if (gt_sz_ % num_samples != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / num_samples - 1);
            const std::uint64_t ploidy_plus_one = ploidy + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_plus_one);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy_plus_one;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index] * ploidy + (i % ploidy_plus_one)] = gl[i];
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                destination[i] = gl[i];
              }
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_pl(site_info& annotations, T& destination)
    {
      if (good())
      {
        if (hts_file_->cur_num_alleles() > 2)
        {
          std::cerr << "multi allelic PL not supported" << std::endl;
          state_ = std::ios::badbit;
          return;
        }

        if (hts_file_->get_cur_format_values_int32("PL", &(gt_), &(gt_sz_)))
        {
          const std::size_t num_samples = sample_ids_.size();
          if (gt_sz_ % num_samples != 0)
          {
            std::cerr << "ERROR: mixed ploidy at site" << std::endl;
            state_ = std::ios::badbit;
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / num_samples - 1);
            const std::uint64_t ploidy_plus_one = ploidy + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_plus_one);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy_plus_one;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index] * ploidy + (i % ploidy_plus_one)] = gt_[i];
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                destination[i] = gt_[i];
              }
            }
            return;
          }
        }

        this->state_ = std::ios::failbit;
      }
    }
    //################################################################//

    //################################################################//
    template <std::size_t VecCnt>
    template <typename... T>
    reader<VecCnt>::reader(const std::string& file_path, T... data_formats)
    {
      static_assert(VecCnt == sizeof...(T), "Number of requested format fields do not match VecCnt template parameter");
      this->init_requested_formats(data_formats...);


      this->hts_file_ = detail::hts_file_base::create_file(file_path);
      if (!this->hts_file_)
      {
        this->state_ = std::ios::badbit;
      }
      else
      {
        this->init_property_fields();
        this->init_headers();
        this->init_sample_ids();
      }
    }

    template <std::size_t VecCnt>
    reader<VecCnt>::reader(reader<VecCnt>&& source) :
      reader_base<VecCnt>(std::move(source))
    {
      source.gt_ = nullptr;
      source.state_ = std::ios::badbit;
    }

    template <std::size_t VecCnt>
    reader<VecCnt>::~reader()
    {

    }

//    std::vector<std::string> reader::chromosomes() const
//    {
//      std::vector<std::string> ret(hts_hdr_->n[BCF_DT_CTG] > 0 ? (unsigned)hts_hdr_->n[BCF_DT_CTG] : 0);
//      for (int i = 0; i < ret.size(); ++i)
//      {
//        ret[i] = hts_hdr()->id[BCF_DT_CTG][i].key;
//      }
//      return ret;
//    }



    template <std::size_t VecCnt>
    template <typename... T>
    reader<VecCnt>& reader<VecCnt>::read(site_info& annotations, T&... destinations)
    {
      static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");
      std::size_t vecs_read = 0;
      while (vecs_read == 0 && this->good())
      {
        this->read_variant_details(annotations);
        vecs_read = this->read_requested_genos(annotations, destinations...);
      }

      return *this;
    }
    //################################################################//



    //################################################################//
    template <std::size_t VecCnt>
    template <typename T>
    indexed_reader<VecCnt>& indexed_reader<VecCnt>::operator>>(variant<T>& destination)
    {
      //static_assert(VecCnt == 1, "Extraction operator only supported with one format field");
      return this->read(destination, destination.data());
    }

    template <std::size_t VecCnt>
    template <typename... T>
    indexed_reader<VecCnt>& indexed_reader<VecCnt>::read(site_info& annotations, T&... destinations)
    {
      static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");
      while (this->good())
      {
        this->read_variant_details(annotations);
        if (this->good() && region_compare(bounding_type_, annotations, region_))
        {
          if (this->read_requested_genos(annotations, destinations...) > 0)
            break;
        }
      }
      return *this;
    }

    template <std::size_t VecCnt>
    template <typename Pred, typename... T>
    indexed_reader<VecCnt>& indexed_reader<VecCnt>::read_if(Pred fn, site_info& annotations, T&... destinations)
    {
      static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");
      while (this->good())
      {
        this->read_variant_details(annotations);
        if (this->good())
        {
          if (fn(annotations) && region_compare(bounding_type_, annotations, region_))
          {
            if (this->read_requested_genos(annotations, destinations...) > 0)
              break;
          }
        }
      }

      return *this;
    }

    template <std::size_t VecCnt>
    template <typename... T>
    indexed_reader<VecCnt>::indexed_reader(const std::string& file_path, const region& reg, T... data_formats) :
      indexed_reader<VecCnt>(file_path, reg, bounding_point::beg, data_formats...)
    {
    }

    template <std::size_t VecCnt>
    template <typename... T>
    indexed_reader<VecCnt>::indexed_reader(const std::string& file_path, const region& reg, bounding_point bounding_type, T... data_formats) :
      region_(reg),
      file_path_(file_path),
      bounding_type_(bounding_type)
    {
      static_assert(VecCnt == sizeof...(T), "Number of requested format fields do not match VecCnt template parameter");

      this->init_requested_formats(data_formats...);

      this->hts_file_ = detail::hts_file_base::create_indexed_file(file_path, reg);
      if (this->hts_file_)
      {
        this->init_property_fields();
        this->init_headers();
        this->init_sample_ids();
      }
      else
        this->state_ = std::ios::badbit;
    }

    template <std::size_t VecCnt>
    indexed_reader<VecCnt>::~indexed_reader()
    {

    }

    template <std::size_t VecCnt>
    std::vector<std::string> indexed_reader<VecCnt>::chromosomes() const
    {
      std::vector<std::string> ret;

      if (this->good())
      {
        return query_chromosomes(file_path_);
      }

      return ret;
    }

    template <std::size_t VecCnt>
    void indexed_reader<VecCnt>::reset_region(const region& reg)
    {
      region_ = reg;
      this->state_ = std::ios::goodbit;
      this->hts_file_ = detail::hts_file_base::create_indexed_file(file_path_, reg);

      if (!this->hts_file_)
        this->state_ = std::ios::failbit;
    }
    //################################################################//



    //################################################################//
    template <std::size_t VecCnt>
    template <typename RandAccessStringIterator, typename RandAccessKVPIterator, typename... Fmt>
    writer<VecCnt>::writer(const std::string& file_path, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, Fmt... data_formats) :
      writer(file_path, options(), samples_beg, samples_end, headers_beg, headers_end, data_formats...)
    {
    }

    template <std::size_t VecCnt>
    template <typename RandAccessStringIterator, typename RandAccessKVPIterator, typename... Fmt>
    writer<VecCnt>::writer(const std::string& file_path, const options& opts, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessKVPIterator headers_beg, RandAccessKVPIterator headers_end, Fmt... data_formats) :
      output_stream_(detail::create_out_stream(file_path, opts.compression)),
      sample_size_(0)
    {
      static_assert(VecCnt == sizeof...(Fmt), "Number of requested format fields do not match VecCnt template parameter");

      this->format_fields_.reserve(sizeof...(Fmt));
      this->init_format_fields(data_formats...);
      std::unordered_set<std::string> unique_info_fields;

      (*output_stream_) << "##fileformat=VCFv4.2\n";

      for (auto it = headers_beg; it != headers_end; ++it)
      {
        if (it->first == "phasing" && (it->second == "partial" || it->second == "none"))
          phase_character_ = '/';

        if (it->first != "FORMAT" && it->first != "fileformat")
        {
          if (it->first == "fileDate")
          {
            std::time_t t = std::time(nullptr);
            char datestr[11];
            if (std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)))
            {
              (*output_stream_) << (std::string("##") + it->first + "=" + std::string(datestr)) << "\n";
            }
          }
          else
          {
            (*output_stream_) << (std::string("##") + it->first + "=" + it->second) << "\n";

            if (it->first == "INFO")
            {
              if (it->second.size() && it->second.front() == '<' && it->second.back() == '>')
              {
                auto header_info = parse_header_value(it->second);
                if (unique_info_fields.emplace(header_info.id).second)
                {
                  info_fields_.emplace_back(header_info.id);
                  int field_int_type = -1;
                  if (header_info.type == "Flag")
                    field_int_type = (int)detail::hts_info_type::bcf_bt_null;
                  info_field_types_.emplace_back(field_int_type);
                }
              }
            }
          }
        }
      }

      for (auto f : this->format_fields_)
      {
        if (f == savvy::fmt::gt)
          (*output_stream_) << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        else if (f == savvy::fmt::hds)
          (*output_stream_) << "##FORMAT=<ID=HDS,Number=.,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage\">\n";
        else if (f == savvy::fmt::ds)
          (*output_stream_) << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage\">\n";
        else if (f == savvy::fmt::gp)
          (*output_stream_) << "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes\">\n";
      }

      (*output_stream_) << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
      for (auto it = samples_beg; it != samples_end; ++it)
      {
        (*output_stream_) << (std::string("\t") + *it);
        ++sample_size_;
      }

      (*output_stream_) << "\n";
    }

    template <std::size_t VecCnt>
    void writer<VecCnt>::init_format_fields(savvy::fmt f)
    {
      this->format_fields_.push_back(f);
    }

    template <std::size_t VecCnt>
    template <typename... Fmt>
    void writer<VecCnt>::init_format_fields(savvy::fmt f, Fmt... other)
    {
      this->format_fields_.push_back(f);
      this->init_format_fields(other...);
    }

    template <std::size_t VecCnt>
    template <typename T>
    writer<VecCnt>& writer<VecCnt>::operator<<(const variant<std::vector<T>>& v)
    {
      return this->write(v, v.data());
    }

    template <std::size_t VecCnt>
    template <typename... T>
    writer<VecCnt>& writer<VecCnt>::write(const site_info& anno, const T&... data)
    {
      static_assert(VecCnt == sizeof...(T), "The number of source vectors must match class template size");

      if (good())
      {
        // VALIDATE VECTOR SIZES
        std::size_t ploidy = 0;
        for (std::size_t i = 0; i < format_fields_.size(); ++i)
        {
          fmt f = format_fields_[i];
          if (f == fmt::gt || f == fmt::hds)
          {
            if (ploidy)
            {
              if ((get_vec(i, data...).size() / sample_size_) != ploidy || (get_vec(i, data...).size() % sample_size_) != 0)
                this->output_stream_->setstate(std::ios::failbit);
            }
            else
            {
              ploidy = (get_vec(i, data...).size() / sample_size_);
              if ((get_vec(i, data...).size() % sample_size_) != 0 || ploidy == 0)
                this->output_stream_->setstate(std::ios::failbit);
            }
          }
          else if (f == fmt::gp)
          {
            if (ploidy)
            {
              if ((get_vec(i, data...).size() / sample_size_) - 1 != ploidy || (get_vec(i, data...).size() % sample_size_) != 0)
                this->output_stream_->setstate(std::ios::failbit);
            }
            else
            {
              ploidy = (get_vec(i, data...).size() / sample_size_) - 1;
              if ((get_vec(i, data...).size() % sample_size_) != 0 || ploidy == 0)
                this->output_stream_->setstate(std::ios::failbit);
            }
          }
          else if (f == fmt::ds)
          {
            if (sample_size_ != get_vec(i, data...).size())
              this->output_stream_->setstate(std::ios::failbit);
          }
        }

        if (good())
        {
          (*output_stream_) << anno.chromosome()
            << "\t" << anno.position()
            << "\t" << std::string(anno.prop("ID").size() ? anno.prop("ID") : ".")
            << "\t" << anno.ref()
            << "\t" << anno.alt()
            << "\t" << std::string(anno.prop("QUAL").size() ? anno.prop("QUAL") : ".")
            << "\t" << std::string(anno.prop("FILTER").size() ? anno.prop("FILTER") : ".");

          std::size_t i = 0;
          for (auto it = info_fields_.begin(); it != info_fields_.end(); ++it)
          {
            if (anno.prop(*it).size())
            {
              if (i == 0)
                (*output_stream_) << "\t";
              else
                (*output_stream_) << ";";

              if (info_field_types_[std::distance(info_fields_.begin(), it)] == (int)detail::hts_info_type::bcf_bt_null)
                (*output_stream_) << (*it); // Flag
              else
                (*output_stream_) << (*it + "=" + anno.prop(*it));

              ++i;
            }
          }

          if (i == 0)
            (*output_stream_) << "\t.";

          for (auto it = format_fields_.begin(); it != format_fields_.end(); ++it)
          {
            if (std::distance(format_fields_.begin(), it) > 0)
              output_stream_->put(':');
            if (*it == fmt::ds)
              (*output_stream_) << "\tDS";
            else if (*it == fmt::hds)
              (*output_stream_) << "\tHDS";
            else if (*it == fmt::gp)
              (*output_stream_) << "\tGP";
            else
              (*output_stream_) << "\tGT";
          }

//          if (VecCnt == 1)
//          {
//            this->write_single_sample_level_data(ploidy, data...);
//          }
//          else
          {
            this->write_multi_sample_level_data(ploidy, data...);
          }
        }
      }

      return *this;
    }

    template <std::size_t VecCnt>
    template <typename... T>
    void writer<VecCnt>::write_multi_sample_level_data(const std::size_t ploidy, const T&... data)
    {
      if (this->good())
      {
        const std::size_t ploidy_plus_one = ploidy + 1;
        std::ostreambuf_iterator<char> out_it(*output_stream_);
        for (std::size_t sample_index = 0; sample_index < sample_size_; ++sample_index)
        {
          for (std::size_t format_index = 0; format_index < format_fields_.size(); ++format_index)
          {
            const auto& v = get_vec(format_index, data...);
            fmt f = format_fields_[format_index];
            if (f == fmt::gt)
            {
              out_it = '\t';

              std::size_t i = sample_index * ploidy;
              if (std::isnan(v[i]))
                out_it = '.';
              else if (v[i] == 0)
                out_it = '0';
              else
                out_it = '1';

              std::size_t end = ploidy + i;
              ++i;
              for ( ; i < end; ++i)
              {
                out_it = phase_character_;

                if (std::isnan(v[i]))
                  out_it = '.';
                else if (v[i] == 0)
                  out_it = '0';
                else
                  out_it = '1';
              }
            }
            else if (f == fmt::gp)
            {
              out_it = '\t';

              std::size_t i = sample_index * ploidy_plus_one;
              if (v[i] == 0)
                out_it = '0';
              else if (v[i] == 1)
                out_it = '1';
              else if (std::isnan(v[i]))
                out_it = '.';
              else
              {
                for (const char c : std::to_string(v[i]))
                  out_it = c;
              }

              std::size_t end = ploidy_plus_one + i;
              ++i;
              for ( ; i < end; ++i)
              {
                out_it = ',';

                if (v[i] == 0)
                  out_it = '0';
                else if (v[i] == 1)
                  out_it = '1';
                else if (std::isnan(v[i]))
                  out_it = '.';
                else
                {
                  for (const char c : std::to_string(v[i]))
                    out_it = c;
                }
              }
            }
            else if (f == fmt::hds)
            {
              out_it = '\t';

              std::size_t i = sample_index * ploidy;
              if (v[i] == 0)
                out_it = '0';
              else if (std::isnan(v[i]))
                out_it = '.';
              else
              {
                for (const char c : std::to_string(v[i]))
                  out_it = c;
              }

              std::size_t end = ploidy + i;
              ++i;
              for ( ; i < end; ++i)
              {
                out_it = ',';

                if (v[i] == 0)
                  out_it = '0';
                else if (std::isnan(v[i]))
                  out_it = '.';
                else
                {
                  for (const char c : std::to_string(v[i]))
                    out_it = c;
                }
              }
            }
            else //if (f == fmt::dosage)
            {
              out_it = '\t';

              if (v[sample_index] == 0)
                out_it = '0';
              else if (std::isnan(v[sample_index]))
                out_it = '.';
              else
              {
                for (const char c : std::to_string(v[sample_index]))
                  out_it = c;
              }
            }
          }
        }

        (*output_stream_) << "\n";
      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    const std::vector<T>& writer<VecCnt>::get_vec(std::size_t offset, const std::vector<T>& m)
    {
      assert(offset == 0);
      return m;
    }

    template <std::size_t VecCnt>
    template <typename T, typename... T2>
    const std::vector<T>& writer<VecCnt>::get_vec(std::size_t offset, const std::vector<T>& m, const T2&... other)
    {
      if ((sizeof...(T2) + 1) - offset == 0)
        return m;
      else
        return get_vec(offset - 1, other...);
    }
    //################################################################//
  }
}

// unqualified lookup error:
// 'operator+' should be declared prior to the call site or in namespace
//inline savvy::vcf::marker::const_iterator operator+(const savvy::vcf::marker::const_iterator& a, savvy::vcf::marker::const_iterator::difference_type n);
//inline savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a);

#endif //LIBSAVVY_VCF_READER_HPP
