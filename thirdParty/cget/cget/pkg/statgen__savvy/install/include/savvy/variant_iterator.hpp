/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_VARIANT_ITERATOR_HPP
#define LIBSAVVY_VARIANT_ITERATOR_HPP

//#include "site_info.hpp"
//
//#include <iterator>
//
//namespace savvy
//{
//  template <typename ReaderType, typename VectorType>
//  class basic_variant_iterator
//  {
//  public:
//    typedef basic_variant_iterator self_type;
//    typedef std::ptrdiff_t difference_type;
//    typedef VectorType value_type;
//    typedef const value_type& reference;
//    typedef const value_type* pointer;
//    typedef std::input_iterator_tag iterator_category;
//
//    basic_variant_iterator() : file_reader_(nullptr) {}
//    basic_variant_iterator(ReaderType& file_reader) :
//      file_reader_(&file_reader)
//    {
//      increment();
//    }
//
//    void increment()
//    {
//      file_reader_->read(m_);
//      if (!file_reader_->good())
//        file_reader_ = nullptr;
//    }
//    self_type& operator++(){ increment(); return *this; }
//    void operator++(int) { increment(); }
//    reference operator*() { return m_; }
//    pointer operator->() { return &m_; }
//    bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
//    bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
//  private:
//    ReaderType* file_reader_;
//    value_type m_;
//  };
//}

#endif //LIBSAVVY_VARIANT_ITERATOR_HPP
