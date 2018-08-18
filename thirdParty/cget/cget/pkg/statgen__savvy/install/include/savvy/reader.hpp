/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_READER_HPP
#define LIBSAVVY_READER_HPP

#include "site_info.hpp"
#include "sav_reader.hpp"
#include "vcf_reader.hpp"
#include "savvy.hpp"

#include <string>
#include <memory>
#include <stdexcept>

namespace savvy
{
  //################################################################//
  class reader_base
  {
  public:
    class sample_iterator
    {
    public:
      typedef sample_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef std::string value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::bidirectional_iterator_tag iterator_category;
    public:
      sample_iterator() : cstring_itr_(nullptr), stdstring_itr_{} {}
      sample_iterator(const char** cstring_itr) : cstring_itr_(cstring_itr), stdstring_itr_{} {}
      sample_iterator(std::vector<std::string>::const_iterator stdstring_itr) : cstring_itr_(nullptr), stdstring_itr_(stdstring_itr) {}

      self_type& operator+=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ += n;
        else
          stdstring_itr_ += n;
        return *this;
      }
      self_type operator+(difference_type n) const { self_type ret(*this); return (ret += n); }
      self_type& operator-=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ -= n;
        else
          stdstring_itr_ -= n;
        return *this;
      }
      self_type operator-(difference_type n) const { self_type ret(*this); return (ret -= n); }
      difference_type operator-(const self_type& b) const
      {
        if (cstring_itr_)
          return cstring_itr_ - b.cstring_itr_;
        return stdstring_itr_ - b.stdstring_itr_;
      }

      self_type& operator--()
      {
        if (cstring_itr_)
          --cstring_itr_;
        else
          --stdstring_itr_;
        return *this;
      }
      self_type operator--(int) { self_type r = *this; --(*this); return r; }
      self_type& operator++()
      {
        if (cstring_itr_)
          ++cstring_itr_;
        else
          ++stdstring_itr_;
        return *this;
      }
      self_type operator++(int) { self_type r = *this; ++(*this); return r; }
      reference operator*()
      {
        if (cstring_itr_)
        {
          tmp_ = *cstring_itr_;
          return tmp_;
        }
        return *stdstring_itr_;
      }
      pointer operator->() { return &(operator*()); }
      bool operator==(const self_type& rhs)
      {
        if (cstring_itr_)
          return cstring_itr_ == rhs.cstring_itr_;
        return stdstring_itr_ == rhs.stdstring_itr_;
      }
      bool operator!=(const self_type& rhs) { return !(*this == rhs); }
    private:
      const char** cstring_itr_;
      std::vector<std::string>::const_iterator stdstring_itr_;
      std::string tmp_;
    };

    virtual ~reader_base() {}

    operator bool() const
    {
      return this->good();
    }

    bool good() const
    {
      if (sav_impl())
        return sav_impl()->good();
      else if (vcf_impl())
        return vcf_impl()->good();
      return false;
    }

    bool fail() const
    {
      if (sav_impl())
        return sav_impl()->fail();
      else if (vcf_impl())
        return vcf_impl()->fail();
      return true;
    }

    bool bad() const
    {
      if (sav_impl())
        return sav_impl()->bad();
      else if (vcf_impl())
        return vcf_impl()->bad();
      return true;
    }

//    template <typename T>
//    bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
//    {
//      if (sav_impl())
//        sav_impl()->read_variant(destination, missing_value);
//      else if (vcf_impl())
//        vcf_impl()->read_variant(destination, missing_value);
//
//      return good();
//    }

    const std::vector<std::string>& info_fields() const;
    const std::vector<std::string>& samples() const;
    const std::vector<std::pair<std::string, std::string>>& headers() const;

    std::vector<std::string> subset_samples(const std::set<std::string>& subset);
  private:
    static const std::vector<std::string> empty_string_vector;
    static const std::vector<std::pair<std::string, std::string>> empty_string_pair_vector;
  protected:
    virtual savvy::sav::reader_base* sav_impl() const = 0;
    virtual savvy::vcf::reader_base<1>* vcf_impl() const = 0;
  };
  //################################################################//

  //################################################################//
  class reader : public reader_base
  {
  public:
    reader() {}
    reader(const std::string& file_path, savvy::fmt data_format);
    ~reader() {}

    template <typename T>
    reader& operator>>(variant<T>& destination);

    template <typename T>
    reader& read(site_info& annotations, T& destination);
  private:
    savvy::sav::reader_base* sav_impl() const { return sav_reader_.get(); }
    savvy::vcf::reader_base<1>* vcf_impl() const { return vcf_reader_.get(); }
  private:
    std::unique_ptr<sav::reader> sav_reader_;
    std::unique_ptr<vcf::reader<1>> vcf_reader_;
  };
  //################################################################//

  //################################################################//
  class indexed_reader : public reader_base
  {
  public:
    indexed_reader() {}
    indexed_reader(const std::string& file_path, const region& reg, savvy::fmt data_format);
    indexed_reader(const std::string& file_path, const region& reg, bounding_point bounding_type, savvy::fmt data_format);
    void reset_region(const region& reg);

    std::vector<std::string> chromosomes() const;

    template <typename T>
    indexed_reader& operator>>(variant<T>& destination);

    template <typename T>
    indexed_reader& read(site_info& annotations, T& destination);

    template <typename Pred, typename T>
    indexed_reader& read_if(Pred fn, site_info& annotations, T& destination);
  private:
    savvy::sav::reader_base* sav_impl() const { return sav_reader_.get(); }
    savvy::vcf::reader_base<1>* vcf_impl() const { return vcf_reader_.get(); }
  private:
    std::unique_ptr<sav::indexed_reader> sav_reader_;
    std::unique_ptr<vcf::indexed_reader<1>> vcf_reader_;
  };
  //################################################################//


  //################################################################//
  template <typename T>
  reader& reader::operator>>(variant<T>& destination)
  {
    return this->read(destination, destination.data());
  }

  template <typename T>
  reader& reader::read(site_info& annotations, T& destination)
  {
    if (sav_impl())
      sav_reader_->read(annotations, destination);
    else if (vcf_impl())
      vcf_reader_->read(annotations, destination);
    return *this;
  }
  //################################################################//

  //################################################################//
  template <typename T>
  indexed_reader& indexed_reader::operator>>(variant<T>& destination)
  {
    return this->read(destination, destination.data());
  }

  template <typename T>
  indexed_reader& indexed_reader::read(site_info& annotations, T& destination)
  {
    if (sav_impl())
      sav_reader_->read(annotations, destination);
    else if (vcf_impl())
      vcf_reader_->read(annotations, destination);
    return *this;
  }

  template <typename Pred, typename T>
  indexed_reader& indexed_reader::read_if(Pred fn, site_info& annoations, T& destination)
  {
    if (sav_reader_)
      sav_reader_->read_if(fn, annoations, destination);
    else if (vcf_reader_)
      vcf_reader_->read_if(fn, annoations, destination);

    return *this;
  }
  //################################################################//
}

#endif //VC_READER_HPP
