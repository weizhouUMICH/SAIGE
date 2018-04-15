/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_M3VCF_READER_HPP
#define LIBSAVVY_M3VCF_READER_HPP

#include "allele_status.hpp"
#include "portable_endian.hpp"

#include <list>
#include <string>
#include <cstring>
#include <cstdint>
#include <vector>
#include <fstream>
#include <functional>

namespace savvy
{
  namespace m3vcf
  {
    namespace detail
    {
      void deserialize_string(std::string& output, std::istream& input);
      void serialize_string(std::ostream& output, const std::string& input);
    }

    class block;

    class marker
    {
    public:
      class const_iterator
      {
      public:
        typedef const_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef allele_status value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::bidirectional_iterator_tag iterator_category;

        const_iterator(const marker& parent, std::uint64_t index) : parent_(&parent), cur_(index) {}
        self_type& operator--(){ --cur_; return *this; }
        self_type operator--(int) { self_type r = *this; --cur_; return r; }
        self_type& operator++(){ ++cur_; return *this; }
        self_type operator++(int) { self_type r = *this; ++cur_; return r; }
        reference operator*() { return (*parent_)[cur_]; }
        pointer operator->() { return &(*parent_)[cur_]; }
        bool operator==(const self_type& rhs) { return cur_ == rhs.cur_; }
        bool operator!=(const self_type& rhs) { return cur_ != rhs.cur_; }
      private:
        const marker* parent_;
        std::size_t cur_;
      };

      marker(block& parent, std::uint32_t offset, const std::string& chromosome, std::uint64_t position, const std::string& ref, const std::string& alt);
      const allele_status& operator[](std::size_t i) const;
      const allele_status& at(std::size_t i) const;
      std::uint64_t haplotype_count() const;
      std::uint64_t pos() const;
      std::string ref() const;
      std::string alt() const;
      const_iterator begin() const { return const_iterator(*this, 0); }
      const_iterator end() const { return const_iterator(*this, haplotype_count()); }


//      bool has_alt_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool has_ref_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      bool is_missing_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      allele_status operator()(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
//      void for_each_allele(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
//      void for_each_missing(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
//      void for_each_non_ref(const std::function<void(allele_status status, std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);

      double calculate_allele_frequency() const;
    private:
      block& parent_;
      std::uint32_t offset_;

      std::string chromosome_;
      std::uint64_t position_;
      std::string ref_;
      std::string alt_;
    };

    class block
    {
    public:
      class const_iterator
      {
      public:
        typedef const_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef marker value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::random_access_iterator_tag iterator_category;

        const_iterator(pointer ptr) : ptr_(ptr) { }
        self_type& operator--() { --ptr_; return *this; }
        self_type operator--(int) { self_type r = *this; --ptr_; return r; }
        self_type& operator++() { ++ptr_; return *this; }
        self_type operator++(int) { self_type r = *this; ++ptr_; return r; }
        reference operator*() { return *ptr_; }
        pointer operator->() { return ptr_; }
        bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
        bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
      private:
        pointer ptr_;
      };

      block() = default;
      block(std::uint64_t sample_size, std::uint8_t ploidy) :
        sample_size_(sample_size),
        unique_haplotype_cnt_(0),
        ploidy_level_(ploidy)
      {}
      const allele_status& sample_haplotype_at(std::uint32_t marker_offset, std::uint64_t haplotype_offset);
      std::uint32_t sample_mapping_at(std::uint64_t haplotype_offset) { return sample_mappings_[haplotype_offset]; }
      const allele_status& unique_haplotype_at(std::uint32_t marker_offset, std::uint64_t unique_haplotype_offset);
      std::uint64_t unique_haplotype_weight_at(std::uint64_t unique_haplotype_offset) { return haplotype_weights_[unique_haplotype_offset]; }
      double calculate_allele_frequency(std::uint32_t marker_off) const;

      const_iterator begin();
      const_iterator end();

      std::uint64_t sample_count() const { return sample_size_; }
      std::uint64_t haplotype_count() const { return sample_size_ * ploidy_level_; }
      std::uint32_t unique_haplotype_count() const { return unique_haplotype_cnt_; }
      std::size_t marker_count() const { return markers_.size(); }
      std::uint8_t ploidy_level() const { return ploidy_level_; }
      const marker& operator[](std::size_t i) const;
      const marker& at(std::size_t i) const;
      static bool read(block& destination, std::istream& source, std::uint64_t sample_count, std::uint8_t ploidy);
      static bool write(std::ostream& destination, const block& source);


      template <typename RandAccessAlleleStatusIt>
      bool add_marker(std::uint64_t position, const std::string& ref, const std::string& alt, RandAccessAlleleStatusIt hap_array_beg, RandAccessAlleleStatusIt hap_array_end);
    private:
      static const allele_status const_has_ref;
      static const allele_status const_has_alt;
      static const allele_status const_is_missing;
      std::vector<marker> markers_;

      //---- GT Data ----//
      std::vector<std::uint64_t> haplotype_weights_;
      std::vector<std::uint32_t> sample_mappings_;
      //std::vector<char> unique_haplotype_matrix_;
      std::vector<std::vector<char>> unique_haplotype_matrix_;
      std::uint64_t sample_size_;
      std::uint32_t unique_haplotype_cnt_;
      std::uint8_t ploidy_level_;
      //---- GT Data ----//
    };

    class reader
    {
    public:
      class input_iterator
      {
      public:
        typedef input_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef marker value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::input_iterator_tag iterator_category;
        typedef block buffer;

        input_iterator() : file_reader_(nullptr), buffer_(nullptr), i_(0) {}
        input_iterator(reader& file_reader, block& buffer) : file_reader_(&file_reader), buffer_(&buffer), i_(0)
        {
          *file_reader_ >> *buffer_;
        }
        void increment()
        {
          ++i_;
          if (i_ >= buffer_->marker_count())
          {
            i_ = 0;
            if (!(*file_reader_ >> *buffer_))
              file_reader_ = nullptr;
          }
        }
        self_type& operator++(){ increment(); return *this; }
        void operator++(int) { increment(); }
        reference operator*() { return (*buffer_)[i_]; }
        pointer operator->() { return &(*buffer_)[i_]; }
        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
      private:
        reader* file_reader_;
        block* buffer_;
        std::uint32_t i_;
      };

      reader(std::istream& input_stream);
      explicit operator bool() const { return input_stream_.good(); }
      bool good() const { return input_stream_.good(); }
      bool fail() const { return input_stream_.fail(); }
      bool bad() const { return input_stream_.bad(); }
      std::uint64_t sample_count() const { return sample_ids_.size(); }
      reader& operator>>(block& destination);
    private:
      const std::string file_path_;
      std::istream& input_stream_;
      std::uint32_t ploidy_level_;
      std::vector<std::string> sample_ids_;
    };

    class writer
    {
    public:
      template <typename RandAccessStringIterator>
      writer(std::ostream& output_stream, const std::string& chromosome, std::uint32_t ploidy, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end) :
        output_stream_(output_stream),
        sample_size_(samples_end - samples_beg),
        ploidy_level_(ploidy)
      {
        std::string version_string("m3vcf\x00\x01\x00\x00", 9);
        output_stream_.write(version_string.data(), version_string.size());


        detail::serialize_string(output_stream_, chromosome);

        std::uint32_t ploidy_nbo = htobe32(ploidy);
        output_stream_.write((char*)(&ploidy_nbo), 4);
        std::uint32_t sample_size_nbo = htobe32(sample_size_);
        output_stream_.write((char*)(&sample_size_nbo), 4);
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          detail::serialize_string(output_stream_, *it);
        }
      }

      writer& operator<<(const block& b)
      {
        if (output_stream_.good())
        {
          if (b.haplotype_count() != sample_size_ * ploidy_level_)
            output_stream_.setstate(std::ios::failbit);
          else
            block::write(output_stream_, b);
        }
        return *this;
      }

    private:
      std::ostream& output_stream_;
      std::uint32_t sample_size_;
      std::uint32_t ploidy_level_;
    };

    template <typename RandAccessAlleleStatusIt>
    bool block::add_marker(std::uint64_t position, const std::string& ref, const std::string& alt, RandAccessAlleleStatusIt hap_array_beg, RandAccessAlleleStatusIt hap_array_end)
    {
      bool ret = false;

      std::size_t hap_array_sz = hap_array_end - hap_array_beg;
      if (hap_array_sz)
      {
        if (markers_.size() == 0)
        {
          sample_mappings_.resize(hap_array_sz, 0xFFFFFFFF);

          char hap;
          switch ((allele_status)((int)(*hap_array_beg)))
          {
            case allele_status::has_ref:
              hap = '0';
              break;
            case allele_status::has_alt:
              hap = '1';
              break;
            case allele_status::is_missing:
              hap = '.';
              break;
          }

          unique_haplotype_matrix_.emplace_back(1, hap);
          sample_mappings_[0] = 0;

          for (std::size_t i = 1; i < hap_array_sz; ++i)
          {
            switch ((allele_status)(int)(*(hap_array_beg + i)))
            {
              case allele_status::has_ref:
                hap = '0';
                break;
              case allele_status::has_alt:
                hap = '1';
                break;
              case allele_status::is_missing:
                hap = '.';
                break;
            }

            std::size_t j = 0;
            for (; j < unique_haplotype_matrix_.size(); ++j)
            {
              if (unique_haplotype_matrix_[j].back() == hap)
                break;
            }

            if (j == unique_haplotype_matrix_.size())
              unique_haplotype_matrix_.emplace_back(1, hap);

            sample_mappings_[i] = static_cast<std::uint32_t>(j);
          }

          markers_.push_back(marker(*this, 0, "", position, ref, alt));
          unique_haplotype_cnt_ = unique_haplotype_matrix_.size();
          ret = true;
        }
        else
        {
          const std::size_t old_marker_size = markers_.size();
          const std::size_t new_marker_size = old_marker_size + 1;
          const std::size_t old_haplotype_matrix_column_count = unique_haplotype_matrix_.size();

          std::int64_t current_savings = static_cast<std::int64_t>(hap_array_sz * old_marker_size) - static_cast<std::int64_t>(old_haplotype_matrix_column_count * old_marker_size + hap_array_sz);

          if (hap_array_sz == sample_mappings_.size())
          {
            std::vector<std::uint32_t> tmp_sample_mappings = sample_mappings_;
            for (std::size_t i = 0; i < hap_array_sz; ++i)
            {
              char hap;
              switch ((allele_status)((int)(*(hap_array_beg + i))))
              {
                case allele_status::has_ref:
                  hap = '0';
                  break;
                case allele_status::has_alt:
                  hap = '1';
                  break;
                case allele_status::is_missing:
                  hap = '.';
                  break;
              }

              if (unique_haplotype_matrix_[tmp_sample_mappings[i]].size() == old_marker_size)
              {
                unique_haplotype_matrix_[tmp_sample_mappings[i]].push_back(hap);
              }
              else
              {
                if (unique_haplotype_matrix_[tmp_sample_mappings[i]][old_marker_size] != hap)
                {
                  std::size_t new_columns_index = old_haplotype_matrix_column_count;
                  for (; new_columns_index < unique_haplotype_matrix_.size(); ++new_columns_index)
                  {
                    std::vector<char>& new_column_to_compare = unique_haplotype_matrix_[new_columns_index];

                    if (std::memcmp(new_column_to_compare.data(), unique_haplotype_matrix_[tmp_sample_mappings[i]].data(), old_marker_size) == 0 && new_column_to_compare.back() == hap)
                      break;
                  }

                  if (new_columns_index == unique_haplotype_matrix_.size())
                  {
                    unique_haplotype_matrix_.emplace_back(unique_haplotype_matrix_[tmp_sample_mappings[i]]);
                    unique_haplotype_matrix_.back().back() = hap;
                  }
                  tmp_sample_mappings[i] = static_cast<std::uint32_t>(new_columns_index); //TODO: Check if in 32-bit range.
                }
              }
            }



            std::int64_t new_savings = static_cast<std::int64_t>(hap_array_sz * new_marker_size) - static_cast<std::int64_t>(unique_haplotype_matrix_.size() * new_marker_size + hap_array_sz);

            if (new_savings >= current_savings)
            {
              markers_.push_back(marker(*this, markers_.size(), "", position, ref, alt));
              sample_mappings_ = std::move(tmp_sample_mappings);
              unique_haplotype_cnt_ = unique_haplotype_matrix_.size();
              ret = true;
            }
            else
            {
              // Undo unique_haplotype_matrix_ edits.
              unique_haplotype_matrix_.erase(unique_haplotype_matrix_.begin() + old_haplotype_matrix_column_count, unique_haplotype_matrix_.end());
              for (auto it = unique_haplotype_matrix_.begin(); it != unique_haplotype_matrix_.end(); ++it)
                it->pop_back();
            }
          }
        }
      }

      return ret;
    }

//    template <typename RandAccessAlleleStatusIt>
//    bool block::add_marker(std::uint64_t position, const std::string& ref, const std::string& alt, RandAccessAlleleStatusIt hap_array_beg, RandAccessAlleleStatusIt hap_array_end)
//    {
//      bool ret = false;
//
//      std::size_t hap_array_sz = hap_array_end - hap_array_beg;
//      if (markers_.size() == 0)
//        sample_mappings_.resize(hap_array_sz, 0xFFFFFFFF);
//      if (hap_array_sz == sample_mappings_.size())
//      {
//        std::int64_t current_savings = static_cast<std::int64_t>(hap_array_sz * markers_.size()) - static_cast<std::int64_t>(unique_haplotype_cnt_ * markers_.size() + hap_array_sz);
//
//        std::list<std::string> new_unique_haps;
//
//        std::string column_string(markers_.size() + 1, '\0');
//        std::vector<std::uint32_t> new_sample_mappings(hap_array_sz);
//
//        for (std::size_t i = 0; i < hap_array_sz; ++i)
//        {
//          std::size_t j = 0;
//          for (; j < markers_.size(); ++j)
//          {
//            column_string[j] = unique_haplotype_matrix_[(j * unique_haplotype_cnt_) + sample_mappings_[i]];
//          }
//          switch (*(hap_array_beg + i))
//          {
//            case allele_status::has_ref: column_string[j] = '0'; break;
//            case allele_status::has_alt: column_string[j] = '1'; break;
//            case allele_status::is_missing: column_string[j] = '.'; break;
//          }
//
//          std::uint32_t unique_hap_offset = 0;
//          auto find_it = new_unique_haps.begin();
//          while (find_it != new_unique_haps.end())
//          {
//            if (*find_it == column_string)
//              break;
//            ++find_it;
//            ++unique_hap_offset;
//          }
//
//          if (find_it == new_unique_haps.end())
//            new_unique_haps.push_back(column_string);
//
//          new_sample_mappings[i] = unique_hap_offset;
//        }
//
//        std::int64_t new_savings = static_cast<std::int64_t>(hap_array_sz * (markers_.size() + 1)) - static_cast<std::int64_t>(new_unique_haps.size() * (markers_.size() + 1) + hap_array_sz);
//        if (new_savings >= current_savings)
//        {
//          markers_.push_back(marker(*this, markers_.size(), "", position, ref, alt));
//          std::vector<char> new_unique_haplotype_matrix(new_unique_haps.size() * markers_.size());
//
//          std::size_t i = 0;
//          for (auto it = new_unique_haps.begin(); it != new_unique_haps.end(); ++it,++i)
//          {
//            for (std::size_t j = 0; j < it->size(); ++j)
//            {
//              new_unique_haplotype_matrix[(j * new_unique_haps.size()) + i] = (*it)[j];
//            }
//          }
//
//          unique_haplotype_cnt_ = new_unique_haps.size();
//          unique_haplotype_matrix_ = std::move(new_unique_haplotype_matrix);
//          sample_mappings_ = std::move(new_sample_mappings);
//
//          ret = true;
//        }
//      }
//
//      return ret;
//    }
  }
}

#endif //VC_M3VCF_READER_HPP
