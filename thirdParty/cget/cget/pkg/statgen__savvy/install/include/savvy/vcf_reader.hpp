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

//namespace savvy
//{
//namespace vcf
//{
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
//}
//}
#include <shrinkwrap/gz.hpp>

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
#include <ctime>
#include <htslib/hts.h>

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
    //################################################################//

    std::vector<std::string> query_chromosomes(const std::string& file_path);

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
        subset_size_(source.subset_size_),
        state_(source.state_),
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
    protected:
      virtual bcf_hdr_t* hts_hdr() const = 0;
      virtual bcf1_t* hts_rec() const = 0;
      virtual bool read_hts_record() = 0;
      void init_sample_ids();
      void init_property_fields();
      void init_headers();

      template <typename T>
      void clear_destinations(T& destination);
      template <typename T, typename... T2>
      void clear_destinations(T& destination, T2&... other_destinations);

      void read_variant_details(site_info& destination);

      template <typename... T>
      std::size_t read_requested_genos(T&... vec);
      template <std::size_t Idx, typename T1>
      bool read_genos_to(fmt data_format, T1& vec);
      template <std::size_t Idx, typename T1, typename... T2>
      bool read_genos_to(fmt data_format, T1& vec, T2&... others);

      template <typename T>
      void read_genotypes_al(T& destination);
      template <typename T>
      void read_genotypes_gt(T& destination);
      template <typename T>
      void read_genotypes_gp(T& destination);
      template <typename T>
      void read_genotypes_ds(T& destination);
      template <typename T>
      void read_genotypes_hds(T& destination);
      template <typename T>
      void read_genotypes_gl(T& destination);
      template <typename T>
      void read_genotypes_pl(T& destination);

      void init_requested_formats(fmt f);
      template <typename... T2>
      void init_requested_formats(fmt f, T2... args);
    protected:
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> sample_ids_;
      std::vector<std::string> property_fields_;
      std::array<fmt, VecCnt> requested_data_formats_;
      std::vector<std::uint64_t> subset_map_;
      std::uint64_t subset_size_;
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

      bool read_hts_record();

      bcf_hdr_t* hts_hdr() const { return hts_hdr_; }
      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      htsFile* hts_file_;
      bcf_hdr_t* hts_hdr_;
      bcf1_t* hts_rec_;
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
      bool read_hts_record();
//      index_reader& seek(const std::string& chromosome, std::uint64_t position);
//
//      static std::string get_chromosome(const index_reader& rdr, const marker& mkr);
      bcf_hdr_t* hts_hdr() const { return bcf_sr_get_header(synced_readers_, 0); }
      bcf1_t* hts_rec() const { return hts_rec_; }
    private:
      region region_;
      std::string file_path_;
      bcf_srs_t* synced_readers_;
      bcf1_t* hts_rec_;
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
      std::vector<fmt> format_fields_;
      std::unique_ptr<std::ostream> output_stream_;
      std::size_t sample_size_;
    };
    //################################################################//






    //################################################################//
    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_headers()
    {
      bcf_hdr_t* hdr = hts_hdr();
      if (hdr)
      {
        this->headers_.reserve(std::size_t(hdr->nhrec - 1));
        for (int i = 0; i < hdr->nhrec; ++i)
        {
          std::string key, val;
          if (hdr->hrec[i]->key && hdr->hrec[i]->value)
          {
            key = hdr->hrec[i]->key;
            val = hdr->hrec[i]->value;
          }
          else if (hdr->hrec[i]->key && hdr->hrec[i]->nkeys) // (hdr->hrec[i]->type == BCF_HL_INFO || hdr->hrec[i]->type == BCF_HL_FLT || hdr->hrec[i]->type == BCF_HL_STR))
          {
            bcf_hrec_t* r = hdr->hrec[i];
            key = r->key;
            std::stringstream ss_val;

            ss_val << "<";
            for (int j = 0; j < r->nkeys - 1; ++j) // minus 1 to remove IDX;
            {
              if (j > 0)
                ss_val << ",";
              if (r->keys[j])
                ss_val << r->keys[j];
              ss_val << "=";
              if (r->vals[j])
                ss_val << r->vals[j];
            }
            ss_val << ">";
            val = ss_val.str();
          }

          if (key.size())
            this->headers_.emplace_back(std::move(key), std::move(val));
          //ret.insert(std::upper_bound(ret.begin(), ret.end(), std::make_pair(key, std::string()), [](const auto& a, const auto& b) { return a.first < b.first; }), {std::move(key), std::move(val)});
        }
      }
    }

    template <std::size_t VecCnt>
    std::vector<std::string> reader_base<VecCnt>::subset_samples(const std::set<std::string>& subset)
    {
      std::vector<std::string> ret;
      const char** beg = hts_hdr() ? (const char**) hts_hdr()->samples : nullptr;
      const char** end = hts_hdr() ? (const char**) (hts_hdr()->samples + bcf_hdr_nsamples(hts_hdr())) : nullptr;
      ret.reserve(std::min(subset.size(), (std::size_t)std::distance(beg, end)));

      subset_map_.clear();
      subset_map_.resize((std::size_t)std::distance(beg, end), std::numeric_limits<std::uint64_t>::max());
      std::uint64_t subset_index = 0;
      for (auto it = beg; it != end; ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map_[std::distance(beg, it)] = subset_index;
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
      bcf_hdr_t* hdr = hts_hdr();
      if (hdr)
      {
        const char **beg = (const char **) (hdr->samples);
        const char **end = (const char **) (hdr->samples + bcf_hdr_nsamples(hdr));
        sample_ids_.reserve(end - beg);

        for (; beg != end; ++beg)
        {
          sample_ids_.emplace_back(*beg);
        }
      }
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::init_property_fields()
    {
      bcf_hdr_t* hdr = hts_hdr();
      if (hdr)
      {

        this->property_fields_ = {"ID", "QUAL", "FILTER"};
        for (int i = 0; i < hdr->nhrec; ++i)
        {
          if (hdr->hrec[i]->type == BCF_HL_INFO)
          {
            bcf_hrec_t* r = hdr->hrec[i];
            for (int j = 0; j < r->nkeys; ++j)
            {
              if (strcmp(r->keys[j], "ID") == 0)
              {
                const char* inf = r->vals[j];
                if (inf)
                  this->property_fields_.emplace_back(inf);
              }
            }
          }
        }
      }
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
    std::size_t reader_base<VecCnt>::read_requested_genos(T&... destinations)
    {
      std::size_t cnt = 0;
      clear_destinations(destinations...);
      bcf_unpack(hts_rec(), BCF_UN_ALL);
      for (int i = 0; i < hts_rec()->n_fmt; ++i)
      {
        int fmt_id = hts_rec()->d.fmt[i].id;
        std::string fmt_key = hts_hdr()->id[BCF_DT_ID][fmt_id].key;

        std::int64_t gt_idx     = std::distance(requested_data_formats_.begin(), std::find(requested_data_formats_.begin(), requested_data_formats_.end(), fmt::ac));
        std::int64_t allele_idx = std::distance(requested_data_formats_.begin(), std::find(requested_data_formats_.begin(), requested_data_formats_.end(), fmt::gt));

        if (fmt_key == "GT" && (gt_idx < VecCnt || allele_idx < VecCnt))
        {
          if (gt_idx < VecCnt)
          {
            cnt += read_genos_to<0>(fmt::ac, destinations...);
          }
          else // allele_idx < VecCnt
          {
            cnt += read_genos_to<0>(fmt::gt, destinations...);
          }
        }
        else if (fmt_key == "DS")
        {
          cnt += read_genos_to<0>(fmt::ds, destinations...);
        }
        else if (fmt_key == "HDS")
        {
          cnt += read_genos_to<0>(fmt::hds, destinations...);
        }
        else if (fmt_key == "GP")
        {
          cnt += read_genos_to<0>(fmt::gp, destinations...);
        }
        else if (fmt_key == "GL")
        {
          cnt += read_genos_to<0>(fmt::gl, destinations...);
        }
        else if (fmt_key == "PL")
        {
          cnt += read_genos_to<0>(fmt::pl, destinations...);
        }
        else
        {
          // Discard
        }
      }
      return cnt;
    }

    template <std::size_t VecCnt>
    template <std::size_t Idx, typename T1>
    bool reader_base<VecCnt>::read_genos_to(fmt data_format, T1& destination)
    {
      bool ret = true;
      if (requested_data_formats_[Idx] == data_format)
      {
        switch (data_format)
        {
          case fmt::gt:
            read_genotypes_al(destination);
            break;
          case fmt::ac:
            read_genotypes_gt(destination);
            break;
          case fmt::gp:
            read_genotypes_gp(destination);
            break;
          case fmt::ds:
            read_genotypes_ds(destination);
            break;
          case fmt::hds:
            read_genotypes_hds(destination);
            break;
          case fmt::gl:
            read_genotypes_gl(destination);
            break;
          case fmt::pl:
            read_genotypes_pl(destination);
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
    bool reader_base<VecCnt>::read_genos_to(fmt data_format, T1& vec, T2&... others)
    {
      if (requested_data_formats_[Idx] == data_format)
      {
        return read_genos_to<Idx>(data_format, vec);
      }
      else
      {
        return read_genos_to<Idx + 1>(data_format, others...);
      }
    }

    template <std::size_t VecCnt>
    void reader_base<VecCnt>::read_variant_details(site_info& destination)
    {
      if (good())
      {
        bool res = true;
        ++allele_index_;
        if (!hts_rec() || allele_index_ >= hts_rec()->n_allele)
        {
          res = read_hts_record();
          if (res)
            bcf_unpack(hts_rec(), BCF_UN_SHR);

          this->allele_index_ = 1;
        }
        
        if (res)
        {
          std::size_t n_info = hts_rec()->n_info;
          std::size_t n_flt = hts_rec()->d.n_flt;
          bcf_info_t* info = hts_rec()->d.info;
          std::unordered_map<std::string, std::string> props;
          props.reserve(n_info + 2);

          if (std::isnan(hts_rec()->qual))
          {
            props["QUAL"] = ".";
          }
          else
          {
            std::string qual(std::to_string(hts_rec()->qual));
            qual.erase(qual.find_last_not_of(".0") + 1); // rtrim zeros.
            props["QUAL"] = std::move(qual);
          }

          std::stringstream ss;
          for (std::size_t i = 0; i < n_flt; ++i)
          {
            if (i > 0)
              ss << ";";
            ss << bcf_hdr_int2id(hts_hdr(), BCF_DT_ID, hts_rec()->d.flt[i]);
          }
          std::string fltr(ss.str());
          if (fltr == "." || fltr == "PASS")
            fltr.clear();
          props["FILTER"] = std::move(fltr);
          props["ID"] = hts_rec()->d.id;


          for (std::size_t i = 0; i < n_info; ++i)
          {
            // bcf_hdr_t::id[BCF_DT_ID][$key].key
            const char* key = hts_hdr()->id[BCF_DT_ID][info[i].key].key;
            if (key)
            {
              switch (info[i].type)
              {
                case BCF_BT_NULL:
                  props[key] = "1"; // Flag present so should be true.
                  break;
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                  props[key] = std::to_string(info[i].v1.i);
                  break;
                case BCF_BT_FLOAT:
                  props[key] = std::to_string(info[i].v1.f);
                  props[key].erase(props[key].find_last_not_of(".0") + 1); // rtrim zeros.
                  break;
                case BCF_BT_CHAR:
                  props[key] = std::string((char*)info[i].vptr, info[i].vptr_len);
                  break;
              }
            }
          }

          destination = site_info(
            std::string(bcf_hdr_id2name(hts_hdr(), hts_rec()->rid)),
            static_cast<std::uint64_t>(hts_rec()->pos + 1),
            std::string(hts_rec()->d.allele[0]),
            std::string(hts_rec()->n_allele > 1 ? hts_rec()->d.allele[allele_index_] : ""),
            std::move(props));
        }

        if (!res)
          this->state_ = std::ios::failbit | std::ios::eofbit; // TODO: Find better way to deal with failed reads vs. eof.

      }
    }

    template <std::size_t VecCnt>
    template <typename T>
    void reader_base<VecCnt>::read_genotypes_al(T& destination)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        if (allele_index_ > 1 || bcf_get_genotypes(hts_hdr(), hts_rec(), &(gt_), &(gt_sz_)) >= 0)
        {
          if (gt_sz_ % samples().size() != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const int allele_index_plus_one = allele_index_ + 1;
            const std::uint64_t ploidy(gt_sz_ / samples().size());
            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (gt_[i] == bcf_gt_missing)
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = std::numeric_limits<typename T::value_type>::quiet_NaN();
                  else if ((gt_[i] >> 1) == allele_index_plus_one)
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = alt_value;
                }
              }
            }
            else
            {
              destination.resize(gt_sz_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (gt_[i] == bcf_gt_missing)
                  destination[i] = std::numeric_limits<typename T::value_type>::quiet_NaN();
                else if ((gt_[i] >> 1) == allele_index_plus_one)
                  destination[i] = alt_value;
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
    void reader_base<VecCnt>::read_genotypes_gt(T& destination)
    {
      if (good())
      {
        const typename T::value_type alt_value = typename T::value_type(1);
        if (allele_index_ > 1 || bcf_get_genotypes(hts_hdr(), hts_rec(), &(gt_), &(gt_sz_)) >= 0)
        {
          if (gt_sz_ % samples().size() != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / samples().size());
            const int allele_index_plus_one = allele_index_ + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_);

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (gt_[i] == bcf_gt_missing)
                    destination[subset_map_[sample_index]] += std::numeric_limits<typename T::value_type>::quiet_NaN();
                  else if ((gt_[i] >> 1) == allele_index_plus_one)
                    destination[subset_map_[sample_index]] += alt_value;
                }
              }
            }
            else
            {
              destination.resize(samples().size());

              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (gt_[i] == bcf_gt_missing)
                  destination[i / ploidy] += std::numeric_limits<typename T::value_type>::quiet_NaN();
                else if ((gt_[i] >> 1) == allele_index_plus_one)
                  destination[i / ploidy] += alt_value;
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
    void reader_base<VecCnt>::read_genotypes_ds(T& destination)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        if (bcf_get_format_float(hts_hdr(),hts_rec(),"DS", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample);
            const typename T::value_type zero_value{0};

            if (subset_map_.size())
            {
              destination.resize(subset_size_);

              float* ds = (float*) (void*) (gt_);
              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (ds[i] != zero_value)
                    destination[subset_map_[sample_index]] = ds[i];
                }
              }
            }
            else
            {
              destination.resize(gt_sz_);

              float* ds = (float*) (void*) (gt_);
              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (ds[i] != zero_value)
                  destination[i] = ds[i];
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
    void reader_base<VecCnt>::read_genotypes_hds(T& destination)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic HDS not supported.
          return;
        }

        if (bcf_get_format_float(hts_hdr(),hts_rec(),"HDS", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample);
            const typename T::value_type zero_value{0};

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy);

              float* hds = (float*) (void*) (gt_);
              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                {
                  if (hds[i] != zero_value)
                    destination[subset_map_[sample_index] * ploidy + (i % ploidy)] = hds[i];
                }
              }
            }
            else
            {
              destination.resize(gt_sz_);

              float* hds = (float*) (void*) (gt_);
              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                if (hds[i] != zero_value)
                  destination[i] = hds[i];
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
    void reader_base<VecCnt>::read_genotypes_gp(T& destination)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        if (bcf_get_format_float(hts_hdr(),hts_rec(),"GP", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
            const std::uint64_t ploidy_plus_one = ploidy + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_plus_one);

              float* gp = (float*) (void*) (gt_);
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

              float* gp = (float*) (void*) (gt_);
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
    void reader_base<VecCnt>::read_genotypes_gl(T& destination)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        if (bcf_get_format_float(hts_hdr(),hts_rec(),"GL", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
            const std::uint64_t ploidy_plus_one = ploidy + 1;

            if (subset_map_.size())
            {
              destination.resize(subset_size_ * ploidy_plus_one);

              float* gp = (float*) (void*) (gt_);
              for (std::size_t i = 0; i < gt_sz_; ++i)
              {
                const std::uint64_t sample_index = i / ploidy_plus_one;
                if (subset_map_[sample_index] != std::numeric_limits<std::uint64_t>::max())
                  destination[subset_map_[sample_index] * ploidy + (i % ploidy_plus_one)] = gp[i];
              }
            }
            else
            {
              destination.resize(gt_sz_);

              float* gp = (float*) (void*) (gt_);
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
    void reader_base<VecCnt>::read_genotypes_pl(T& destination)
    {
      if (good())
      {
        if (hts_rec()->n_allele > 2)
        {
          state_ = std::ios::failbit; // multi allelic GP not supported.
          return;
        }

        if (bcf_get_format_int32(hts_hdr(),hts_rec(),"PL", &(gt_), &(gt_sz_)) >= 0)
        {
          int num_samples = hts_hdr()->n[BCF_DT_SAMPLE];
          if (gt_sz_ % num_samples != 0)
          {
            // TODO: mixed ploidy at site error.
          }
          else
          {
            const std::uint64_t ploidy(gt_sz_ / hts_rec()->n_sample - 1);
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
    reader<VecCnt>::reader(const std::string& file_path, T... data_formats) :
      hts_file_(bcf_open(file_path.c_str(), "r")),
      hts_hdr_(nullptr),
      hts_rec_(bcf_init1())
    {
      static_assert(VecCnt == sizeof...(T), "Number of requested format fields do not match VecCnt template parameter");
      this->init_requested_formats(data_formats...);
      if (!hts_file_ || !hts_rec_)
      {
        this->state_ = std::ios::badbit;
      }
      else
      {
        hts_hdr_ = bcf_hdr_read(hts_file_);
        if (!hts_hdr_)
          this->state_ = std::ios::badbit;
        else
        {
          this->init_property_fields();
          this->init_headers();
          this->init_sample_ids();
        }
      }
    }

    template <std::size_t VecCnt>
    reader<VecCnt>::reader(reader<VecCnt>&& source) :
      reader_base<VecCnt>(std::move(source)),
      hts_file_(source.hts_file_),
      hts_hdr_(source.hts_hdr_),
      hts_rec_(source.hts_rec_)
    {
      source.hts_file_ = nullptr;
      source.hts_hdr_ = nullptr;
      source.hts_rec_ = nullptr;
      source.gt_ = nullptr;
      source.state_ = std::ios::badbit;
    }

    template <std::size_t VecCnt>
    reader<VecCnt>::~reader()
    {
      try
      {
        if (hts_hdr_)
          bcf_hdr_destroy(hts_hdr_);
        if (hts_file_)
          bcf_close(hts_file_);
        if (hts_rec_)
          bcf_destroy1(hts_rec_);
      }
      catch (...)
      {
        assert(!"bcf_hdr_destroy or bcf_close threw exception!");
      }
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
    bool reader<VecCnt>::read_hts_record()
    {
      if (bcf_read(hts_file_, hts_hdr_, hts_rec_) >= 0)
      {
        return true;
      }
      return false;
    }

    template <std::size_t VecCnt>
    template <typename... T>
    reader<VecCnt>& reader<VecCnt>::read(site_info& annotations, T&... destinations)
    {
      static_assert(VecCnt == sizeof...(T), "The number of destination vectors must match class template size");
      std::size_t vecs_read = 0;
      while (vecs_read == 0 && this->good())
      {
        this->read_variant_details(annotations);
        vecs_read = this->read_requested_genos(destinations...);
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
          if (this->read_requested_genos(destinations...) > 0)
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
            if (this->read_requested_genos(destinations...) > 0)
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
      synced_readers_(bcf_sr_init()),
      hts_rec_(nullptr),
      bounding_type_(bounding_type)
    {
      static_assert(VecCnt == sizeof...(T), "Number of requested format fields do not match VecCnt template parameter");

      this->init_requested_formats(data_formats...);
      std::stringstream contigs;
      contigs << reg.chromosome();
      if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
        contigs << ":" << reg.from() << "-" << reg.to();


      if (bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0) == 0 && bcf_sr_add_reader(synced_readers_, file_path_.c_str()) == 1)
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
      if (synced_readers_)
        bcf_sr_destroy(synced_readers_);
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
//      if (this->good())
//      {
        region_ = reg;

        if (synced_readers_)
          bcf_sr_destroy(synced_readers_);
        synced_readers_ = bcf_sr_init();
        hts_rec_ = nullptr;
        this->state_ = std::ios::goodbit;

        std::stringstream contigs;
        contigs << reg.chromosome();
        if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
          contigs << ":" << reg.from() << "-" << reg.to();

        if (bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0) != 0 || bcf_sr_add_reader(synced_readers_, file_path_.c_str()) != 1)
          this->state_ = std::ios::failbit;
//      }
    }

    template <std::size_t VecCnt>
    bool indexed_reader<VecCnt>::read_hts_record()
    {
      if (bcf_sr_next_line(synced_readers_) && (hts_rec_ = bcf_sr_get_line(synced_readers_, 0)))
      {
        return true;
      }
      return false;
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
      output_stream_(opts.compression == compression_type::none ? std::unique_ptr<std::ostream>(new std::ofstream(file_path)) : std::unique_ptr<std::ostream>(new shrinkwrap::bgzf::ostream(file_path))),
      sample_size_(0)
    {
      static_assert(VecCnt == sizeof...(Fmt), "Number of requested format fields do not match VecCnt template parameter");

      this->format_fields_.reserve(sizeof...(Fmt));
      this->init_format_fields(data_formats...);

      (*output_stream_) << "##fileformat=VCFv4.2\n";


      for (auto it = headers_beg; it != headers_end; ++it)
      {
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
                std::string header_value = it->second;
                header_value.resize(header_value.size() - 1);

                auto curr_pos = header_value.begin() + 1;
                auto comma_pos = std::find(curr_pos, header_value.end(), ',');

                while (comma_pos != header_value.end())
                {
                  auto equals_pos = std::find(curr_pos, comma_pos, '=');
                  if (equals_pos != comma_pos)
                  {
                    std::string key(curr_pos, equals_pos);
                    std::string val(equals_pos + 1, comma_pos);

                    if (key == "ID")
                      info_fields_.emplace_back(std::move(val));
                  }

                  curr_pos = comma_pos + 1;
                  comma_pos = std::find(curr_pos, header_value.end(), ',');
                }

                auto equals_pos = std::find(curr_pos, comma_pos, '=');
                if (equals_pos != comma_pos)
                {
                  std::string key(curr_pos, equals_pos);
                  std::string val(equals_pos + 1, comma_pos);

                  if (key == "ID")
                    info_fields_.emplace_back(std::move(val));

                  curr_pos = comma_pos + 1;
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
                out_it = '|';

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
