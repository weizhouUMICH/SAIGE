/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SPARSE_VECTOR_HPP
#define LIBSAVVY_SPARSE_VECTOR_HPP

#include <vector>
#include <algorithm>

namespace savvy
{
  template<typename T>
  class compressed_vector
  {
  public:
    typedef T value_type;
    typedef compressed_vector<T> self_type;
    static const T const_value_type;

    class iterator
    {
    public:
      typedef iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef value_type& reference;
      typedef value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      iterator() : vec_(nullptr), beg_(nullptr), cur_(beg_) {}
      iterator(compressed_vector& parent, std::size_t off) :
        vec_(&parent),
        beg_(parent.values_.data()),
        cur_(parent.values_.data() + off)
      {

      }

      std::size_t offset() const
      {
        return vec_->offsets_[cur_ - beg_];
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() { return *cur_; }
      pointer operator->() { return cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      compressed_vector* vec_;
      const value_type*const beg_;
      value_type* cur_;
    };

    class const_iterator
    {
    public:
      typedef const_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::input_iterator_tag iterator_category;

      const_iterator() : vec_(nullptr), beg_(nullptr), cur_(beg_) {}
      const_iterator(const compressed_vector& parent, std::size_t off) :
        vec_(&parent),
        beg_(parent.values_.data()),
        cur_(beg_ + off)
      {

      }

      std::size_t offset() const
      {
        return vec_->offsets_[cur_ - beg_];
      }

      self_type operator++()
      {
        self_type ret = *this;
        ++cur_;
        return ret;
      }

      void operator++(int) { ++cur_; }
      reference operator*() const { return *cur_; }
      const pointer operator->() const { return cur_; }
      bool operator==(const self_type& rhs) const { return (cur_ == rhs.cur_); }
      bool operator!=(const self_type& rhs) const { return (cur_ != rhs.cur_); }
    private:
      const compressed_vector* vec_;
      const value_type*const beg_;
      const value_type* cur_;
    };

    compressed_vector(std::size_t sz = 0)
    {
      resize(sz);
    }

    value_type& operator[](std::size_t pos)
    {
      if (offsets_.size() && offsets_.back() < pos)
      {
        offsets_.emplace_back(pos);
        values_.emplace_back();
        return values_.back();
      }
      else
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
        if (it == offsets_.end() || *it != pos)
        {
          it = offsets_.insert(it, pos);
          return *(values_.insert(values_.begin() + std::distance(offsets_.begin(), it), value_type()));
        }
        return values_[it - offsets_.begin()];
      }
    }

    const_iterator cbegin() const  { return const_iterator(*this, 0); }
    const_iterator cend() const { return const_iterator(*this, this->values_.size()); }

    const_iterator begin() const  { return this->cbegin(); }
    const_iterator end() const { return this->cend(); }

    iterator begin() { return iterator(*this, 0); }
    iterator end() { return iterator(*this, this->values_.size()); }

    const value_type& operator[](std::size_t pos) const
    {
      auto it = std::lower_bound(offsets_.begin(), offsets_.end(), pos);
      if (it == offsets_.end() || *it != pos)
        return const_value_type;
      return values_[it - offsets_.begin()];
    }

    void resize(std::size_t sz, value_type val = value_type())
    {
      if (sz < size_)
      {
        auto it = std::lower_bound(offsets_.begin(), offsets_.end(), sz);
        offsets_.erase(it, offsets_.end());
        values_.resize(offsets_.size());
      }
      else if (val != value_type())
      {
        values_.resize(sz, val);
        offsets_.resize(sz);
        for (std::size_t i = size_; i < sz; ++i)
          offsets_[i] = i;
      }

      size_ = sz;
    }

    void clear()
    {
      resize(0);
    }

    const std::size_t* const index_data() const { return offsets_.data(); }
    const value_type* const value_data() const { return values_.data(); }
    std::size_t size() const { return size_; }
    std::size_t non_zero_size() const { return values_.size(); }
  private:
    std::vector<value_type> values_;
    std::vector<std::size_t> offsets_;
    std::size_t size_;
  };

  template <typename T>
  const T compressed_vector<T>::const_value_type = T();
}

#endif //LIBSAVVY_SPARSE_VECTOR_HPP
