#ifndef SHRINKWRAP_ISTREAM_HPP
#define SHRINKWRAP_ISTREAM_HPP

#include "xz.hpp"
#include "gz.hpp"
#include "zstd.hpp"

#include <streambuf>
#include <memory>

namespace shrinkwrap
{
  namespace detail
  {
    template<typename T, typename ...Args>
    std::unique_ptr<T> make_unique(Args&& ...args)
    {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
  }

  class istream : public std::istream
  {
  public:
    istream(const std::string& file_path)
      :
      std::istream(nullptr)
    {
      FILE* fp = fopen(file_path.c_str(), "rb");

      int first_byte = fgetc(fp);
      ungetc(first_byte, fp);

      switch (char(first_byte))
      {
        case '\x1F':
          sbuf_ = detail::make_unique<::shrinkwrap::bgzf::ibuf>(fp);
          break;
        case char('\xFD'):
          sbuf_ = detail::make_unique<::shrinkwrap::xz::ibuf>(fp);
          break;
        case '\x28':
          sbuf_ = detail::make_unique<::shrinkwrap::zstd::ibuf>(fp);
          break;
        default:
          throw std::runtime_error("raw files not yet supported.");

      }

      this->rdbuf(sbuf_.get());
    }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
    istream(istream&& src)
      :
      std::istream(src.sbuf_.get()),
      sbuf_(std::move(src.sbuf_))
    {
    }

    istream& operator=(istream&& src)
    {
      if (&src != this)
      {
        std::istream::operator=(std::move(src));
        sbuf_ = std::move(src.sbuf_);
      }
      return *this;
    }
#endif
  private:
    std::unique_ptr<std::streambuf> sbuf_;
  };
}

#endif //SHRINKWRAP_ISTREAM_HPP