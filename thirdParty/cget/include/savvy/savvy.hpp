/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_SAVVY_HPP
#define LIBSAVVY_SAVVY_HPP

#include <string>
#include <functional>
#include <memory>
#include "utility.hpp"

namespace savvy
{
  namespace detail
  {
    bool has_extension(const std::string& fullString, const std::string& ext);

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
//    template<typename F, typename Tuple, std::size_t... S>
//    decltype(auto) apply_impl(F&& fn, Tuple&& t, std::index_sequence<S...>)
//    {
//      return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
//    }
//
//    template<typename F, typename Tuple>
//    decltype(auto) apply(F&& fn, Tuple&& t)
//    {
//      std::size_t constexpr tuple_size
//        = std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
//      return apply_impl(std::forward<F>(fn),
//        std::forward<Tuple>(t),
//        std::make_index_sequence<tuple_size>());
//    }



//    struct variadic_file_opener
//    {
//      template<typename... TupleArgs, typename Fn, typename...>
//      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler)
//      {
//        apply(handler, std::move(readers));
//      }
//
//      template<typename... TupleArgs, typename Fn, typename File, typename... AddlFiles>
//      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler, const File& file_path, const AddlFiles& ... addl_file_paths)
//      {
//        if (detail::has_extension(file_path, ".sav"))
//        {
//          variadic_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::sav::reader(file_path))), std::ref(handler), addl_file_paths...);
//        }
////        else if (detail::has_extension(file_path, ".m3vcf"))
////        {
////          std::ifstream ifs(file_path);
////          variadic_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::m3vcf::reader(ifs))), std::ref(handler), addl_file_paths...);
////        }
//        else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
//        {
//          variadic_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::vcf::reader(file_path))), std::ref(handler), addl_file_paths...);
//        }
//      }
//    };

//    struct variadic_indexed_file_opener
//    {
//      template<typename... TupleArgs, typename Fn, typename...>
//      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler)
//      {
//        apply(handler, std::move(readers));
//      }
//
//      template<typename... TupleArgs, typename Fn, typename File, typename... AddlFiles>
//      void operator()(std::tuple<TupleArgs...>&& readers, Fn&& handler, const File& file_path, const AddlFiles& ... addl_file_paths)
//      {
//        if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
//        {
//          variadic_indexed_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::vcf::index_reader(file_path))), std::ref(handler), addl_file_paths...);
//        }
////        else if (detail::has_extension(file_path, ".sav"))
////        {
////          std::ifstream ifs(file_path);
////          variadic_indexed_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::sav::reader(ifs))), std::ref(handler), addl_file_paths...);
////        }
////        else if (detail::has_extension(file_path, ".m3vcf"))
////        {
////          std::ifstream ifs(file_path);
////          variadic_indexed_file_opener::operator()(std::tuple_cat(std::move(readers), std::make_tuple(savvy::m3vcf::reader(ifs))), std::ref(handler), addl_file_paths...);
////        }
//      }
//    };
#endif
  }
#if 0
  template <typename Fn, typename File, typename File2, typename... AddlFiles>
  void open_files(Fn&& handler, const File& file_path, const File2& file_path2, const AddlFiles&... addl_file_paths)
  {
    detail::variadic_file_opener()(std::tuple<>(), std::move(handler), file_path, file_path2, addl_file_paths...);
  }

  template <typename... TplArgs, typename Fn>
  void open_files(const std::tuple<TplArgs...>& file_paths, Fn&& handler)
  {
    detail::apply(detail::variadic_file_opener(), std::tuple_cat(std::make_tuple(std::tuple<>(), std::ref(handler)), file_paths));
  }

  template <typename Fn>
  void open_file(const std::string& file_path, Fn&& handler)
  {
    if (detail::has_extension(file_path, ".sav"))
    {
      savvy::sav::reader input(file_path);
      handler(std::move(input));
    }
//    else if (detail::has_extension(file_path, ".m3vcf"))
//    {
//      std::ifstream ifs(file_path);
//      savvy::m3vcf::reader input(ifs);
//      handler(std::move(input));
//    }
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
    {
      savvy::vcf::reader input(file_path);
      handler(std::move(input));
    }
  }
#endif
//  template <typename Fn, typename File, typename File2, typename... AddlFiles>
//  void open_indexed_files(Fn&& handler, const File& file_path, const File2& file_path2, const AddlFiles&... addl_file_paths)
//  {
//    detail::variadic_indexed_file_opener()(std::tuple<>(), std::move(handler), file_path, file_path2, addl_file_paths...);
//  }
//
//  template <typename... TplArgs, typename Fn>
//  void open_indexed_files(const std::tuple<TplArgs...>& file_paths, Fn&& handler)
//  {
//    detail::apply(detail::variadic_indexed_file_opener(), std::tuple_cat(std::make_tuple(std::tuple<>(), std::ref(handler)), file_paths));
//  }
//
//  template <typename Fn>
//  void open_indexed_file(const std::string& file_path, Fn&& handler)
//  {
//    if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, "vcf.gz") || detail::has_extension(file_path, ".bcf"))
//    {
//      savvy::vcf::index_reader input(file_path);
//      handler(std::move(input));
//    }
////    else if (detail::has_extension(file_path, ".sav"))
////    {
////      std::ifstream ifs(file_path);
////      savvy::sav::reader input(ifs);
////      handler(std::move(input));
////    }
////    else if (detail::has_extension(file_path, ".m3vcf"))
////    {
////      std::ifstream ifs(file_path);
////      savvy::m3vcf::reader input(ifs);
////      handler(std::move(input));
////    }
//  }

//  template <typename Fn>
//  void open_marker_file(const std::string& file_path, Fn&& handler)
//  {
//    open_marker_file(file_path, std::move(handler));
//  }

  namespace detail
  {
    template<typename Fn>
    class reader_iterator
    {
    public:
      reader_iterator(Fn& fn)
        :handler_(fn)
      {}

      template<typename Reader>
      void operator()(Reader&& r)
      {
        typename Reader::input_iterator::buffer buf;
        typename Reader::input_iterator eof;
        typename Reader::input_iterator it(r, buf);

        while (it != eof)
        {
          handler_(*it);
          ++it;
        }
      }

    private:
      Fn& handler_;
    };
  }

  template <typename Fn>
  void iterate_file(const std::string& file_path, Fn&& handler)
  {
    detail::reader_iterator<Fn> mrkr_iter_functor(handler);
    open_file(file_path, mrkr_iter_functor);
  }

  template <typename ReaderType, typename VecType>
  std::uint8_t get_ploidy(const ReaderType& r, const VecType& v)
  {
    return static_cast<std::uint8_t>(v.size() / (r.samples_end() - r.samples_begin()));
  }

  template <typename ReaderType, typename VecType>
  std::tuple<std::uint8_t, bool> get_validated_ploidy(const ReaderType& r, const VecType& v)
  {
    const std::uint64_t sample_count = (r.samples_end() - r.samples_begin());
    std::tuple<std::uint8_t, bool> ret(static_cast<std::uint8_t>(v.size() / sample_count), (v.size() % sample_count) == 0);
    return ret;
  }

//  template <typename R>
//  class region_query
//  {
//  public:
//    class const_iterator
//    {
//    public:
//      typedef const_iterator self_type;
//      typedef std::ptrdiff_t difference_type;
//      typedef typename R::input_iterator::value_type value_type;
//      typedef const value_type& reference;
//      typedef const value_type* pointer;
//      typedef std::input_iterator_tag iterator_category;
//
//      const_iterator(region_query& parent, typename R::input_iterator&& beg_itr) :
//        parent_query_(&parent),
//        itr_(std::move(beg_itr))
//      {
//      }
//
//      void increment()
//      {
//        ++itr_;
//        if (*this == parent_query_->end() || itr_->pos() > parent_query_->end_pos() || R::get_chromosome(parent_query_->reader(), *itr_) != parent_query_->chromosome())
//        {
//          itr_ = typename R::input_iterator();
//        }
//      }
//
//      self_type& operator++(){ increment(); return *this; }
//      void operator++(int) { increment(); }
//      reference operator*() { return *itr_; }
//      pointer operator->() { return &(*itr_); }
//      bool operator==(const self_type& rhs) { return (itr_ == rhs.itr_); }
//      bool operator!=(const self_type& rhs) { return (itr_ != rhs.itr_); }
//    private:
//      region_query* parent_query_;
//      typename R::input_iterator itr_;
//    };
//
//    region_query(R& index_reader, const std::string& chromosome, std::uint64_t start_pos, std::uint64_t end_pos) :
//      reader_(index_reader),
//      chromosome_(chromosome),
//      start_pos_(start_pos),
//      end_pos_(end_pos)
//    {
//
//    }
//
//    const_iterator begin()
//    {
//      reader_.seek(chromosome_, start_pos_);
//      const_iterator ret(*this, typename R::input_iterator(reader_, buf_));
//      return ret;
//    }
//
//    const_iterator end() { return const_iterator(*this, typename R::input_iterator()); }
//    const R& reader() const { return reader_; }
//    const std::string& chromosome() { return chromosome_; }
//    std::uint64_t end_pos() const { return end_pos_; }
//  private:
//    R& reader_;
//    typename R::input_iterator::buffer buf_;
//
//    std::string chromosome_;
//    std::uint64_t start_pos_;
//    std::uint64_t end_pos_;
//  };
//
//  template <typename R>
//  region_query<R> make_region_query(R& index_reader, const std::string& chromosome, std::uint64_t start_pos, std::uint64_t end_pos)
//  {
//    region_query<R> ret(index_reader, chromosome, start_pos, end_pos);
//    return ret;
//  }
}
#endif //LIBSAVVY_SAVVY_HPP