#ifndef SHRINKWRAP_XZ_HPP
#define SHRINKWRAP_XZ_HPP

#include <streambuf>
#include <array>
#include <vector>
#include <stdio.h>
#include <lzma.h>
#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

namespace shrinkwrap
{
  namespace xz
  {
    class ibuf : public std::streambuf
    {
    public:
      ibuf(FILE* fp)
        :
        decoded_position_(0),
        discard_amount_(0),
        fp_(fp),
        put_back_size_(0),
        lzma_index_(nullptr),
        at_block_boundary_(true),
        lzma_block_decoder_(LZMA_STREAM_INIT)
      {
        if (fp_)
        {
          fread(stream_header_.data(), stream_header_.size(), 1, fp_); // TODO: handle error.
          lzma_res_ = lzma_stream_header_decode(&stream_header_flags_, stream_header_.data());
          if (lzma_res_ != LZMA_OK)
          {
            // TODO: handle error.
          }
        }
        char* end = ((char*) decompressed_buffer_.data()) + decompressed_buffer_.size();
        setg(end, end, end);
      }

      ibuf(const std::string& file_path) :ibuf(fopen(file_path.c_str(), "rb")) {}

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      ibuf(ibuf&& src)
        :
        std::streambuf(std::move(src))
      {
        this->move(std::move(src));
      }

      ibuf& operator=(ibuf&& src)
      {
        if (&src != this)
        {
          std::streambuf::operator=(std::move(src));
          this->destroy();
          this->move(std::move(src));
        }

        return *this;
      }
#endif

      virtual ~ibuf()
      {
        this->destroy();
      }

    protected:
      virtual std::streambuf::int_type underflow()
      {
        if (!fp_)
          return traits_type::eof();
        if (gptr() < egptr()) // buffer not exhausted
          return traits_type::to_int_type(*gptr());

        while (gptr() >= egptr() && lzma_res_ == LZMA_OK)
        {
          lzma_block_decoder_.next_out = decompressed_buffer_.data();
          lzma_block_decoder_.avail_out = decompressed_buffer_.size();

          if (at_block_boundary_)
          {
            std::vector<std::uint8_t> block_header(LZMA_BLOCK_HEADER_SIZE_MAX);
            if (lzma_block_decoder_.avail_in == 0 && !feof(fp_) && !ferror(fp_))
            {
              replenish_compressed_buffer();
            }
            // TODO: make sure avail_in is greater than 0;
            std::memcpy(block_header.data(), lzma_block_decoder_.next_in, 1);
            ++(lzma_block_decoder_.next_in);
            --(lzma_block_decoder_.avail_in);

            if (block_header[0] == 0x00)
            {
              // Index indicator found
              lzma_res_ = LZMA_STREAM_END;
            }
            else
            {
              lzma_block_.version = 0;
              lzma_block_.check = stream_header_flags_.check;
              lzma_block_.filters = lzma_block_filters_buf_.data();
              lzma_block_.header_size = lzma_block_header_size_decode (block_header[0]);

              std::size_t bytes_already_copied = 0;
              if (lzma_block_decoder_.avail_in < (lzma_block_.header_size - 1))
              {
                bytes_already_copied = lzma_block_decoder_.avail_in;
                std::memcpy(&block_header[1], lzma_block_decoder_.next_in, bytes_already_copied);
                lzma_block_decoder_.avail_in -= bytes_already_copied;
                lzma_block_decoder_.next_in += bytes_already_copied;
                assert(lzma_block_decoder_.avail_in == 0);
                replenish_compressed_buffer();
              }

              // TODO: make sure avail_in is greater than (lzma_block_.header_size - 1) - bytes_already_copied.
              std::size_t bytes_left_to_copy = (lzma_block_.header_size - 1) - bytes_already_copied;
              std::memcpy(&block_header[1 + bytes_already_copied], lzma_block_decoder_.next_in, bytes_left_to_copy);
              lzma_block_decoder_.avail_in -= bytes_left_to_copy;
              lzma_block_decoder_.next_in += bytes_left_to_copy;

              lzma_res_ = lzma_block_header_decode(&lzma_block_, nullptr, block_header.data());
              if (lzma_res_ != LZMA_OK)
              {
                // TODO: handle error.
              }
              else
              {
                lzma_res_ = lzma_block_decoder(&lzma_block_decoder_, &lzma_block_);
                // TODO: handle error.
              }
            }
            at_block_boundary_ = false;
          }

          if (lzma_res_ == LZMA_OK)
          {
            if (lzma_block_decoder_.avail_in == 0 && !feof(fp_) && !ferror(fp_))
            {
              replenish_compressed_buffer();
            }

            assert(lzma_block_decoder_.avail_in > 0);

            lzma_ret r = lzma_code(&lzma_block_decoder_, LZMA_RUN);
            if (r == LZMA_STREAM_END)
            {
              // End of block.
              at_block_boundary_ = true;
              r = LZMA_OK;
            }
            lzma_res_ = r;
          }

          char* start = ((char*) decompressed_buffer_.data());
          setg(start, start, start + (decompressed_buffer_.size() - lzma_block_decoder_.avail_out));
          decoded_position_ += (egptr() - gptr());

          if (discard_amount_ > 0)
          {
            std::uint64_t advance_amount = discard_amount_;
            if ((egptr() - gptr()) < advance_amount)
              advance_amount = (egptr() - gptr());
            setg(start, gptr() + advance_amount, egptr());
            discard_amount_ -= advance_amount;
          }
        }

        if (lzma_res_ == LZMA_STREAM_END && gptr() >= egptr())
          return traits_type::eof();
        else if (lzma_res_ != LZMA_OK && lzma_res_ != LZMA_STREAM_END)
          return traits_type::eof();

        return traits_type::to_int_type(*gptr());
      }

      virtual std::streambuf::pos_type seekoff(std::streambuf::off_type off, std::ios_base::seekdir way, std::ios_base::openmode which)
      {
        std::uint64_t current_position = decoded_position_ - (egptr() - gptr());
        current_position += discard_amount_; // TODO: overflow check.

        pos_type pos{off_type(current_position)};

        if (off == 0 && way == std::ios::cur)
          return pos; // Supports tellg for streams that can't seek.

        if (way == std::ios::cur)
        {
          pos = pos + off;
        }
        else if (way == std::ios::end)
        {
          if (!lzma_index_)
          {
            if (!init_index())
              return pos_type(off_type(-1));
          }

          pos = pos_type(lzma_index_uncompressed_size(lzma_index_)) + off;
        }
        else
        {
          pos = off;
        }

        return seekpos(pos, which);
      }

      virtual std::streambuf::pos_type seekpos(std::streambuf::pos_type pos, std::ios_base::openmode which)
      {
        if (fp_ == 0 || sync())
          return pos_type(off_type(-1));

        if (!lzma_index_) //stream_flags_.backward_size == LZMA_VLI_UNKNOWN)
        {
          if (!init_index())
            return pos_type(off_type(-1));
        }

        if (lzma_index_iter_locate(&lzma_index_itr_, (std::uint64_t) off_type(pos))) // Returns true on failure.
          return pos_type(off_type(-1));

        long seek_amount = (lzma_index_itr_.block.compressed_file_offset > std::numeric_limits<long>::max() ? std::numeric_limits<long>::max() : static_cast<long>(lzma_index_itr_.block.compressed_file_offset));
        if (fseek(fp_, seek_amount, SEEK_SET))
          return pos_type(off_type(-1));

        discard_amount_ = off_type(pos) - lzma_index_itr_.block.uncompressed_file_offset;
        decoded_position_ = lzma_index_itr_.block.uncompressed_file_offset;

        at_block_boundary_ = true;
        lzma_block_decoder_.next_in = nullptr;
        lzma_block_decoder_.avail_in = 0;
        char* end = ((char*) decompressed_buffer_.data()) + decompressed_buffer_.size();
        setg(end, end, end);

        return pos;
      }

    private:
      //ixzbuf(const ixzbuf& src) = delete;
      //ixzbuf& operator=(const ixzbuf& src) = delete;

      void destroy()
      {
        if (lzma_block_decoder_.internal)
          lzma_end(&lzma_block_decoder_);
        if (lzma_index_)
          lzma_index_end(lzma_index_, nullptr);
        if (fp_)
          fclose(fp_);
      }

      void move(ibuf&& src)
      {
        stream_header_flags_ = src.stream_header_flags_;
        stream_footer_flags_ = src.stream_footer_flags_;
        lzma_block_decoder_ = src.lzma_block_decoder_;
        if (src.lzma_block_decoder_.internal)
          src.lzma_block_decoder_.internal = nullptr;
        lzma_block_ = src.lzma_block_;
        lzma_block_filters_buf_ = src.lzma_block_filters_buf_; // TODO: handle filter.options
        lzma_index_itr_ = src.lzma_index_itr_; // lzma_index_iter_init() doesn't allocate any memory, thus there is no lzma_index_iter_end().
        stream_header_ = src.stream_header_;
        stream_footer_ = src.stream_footer_;
        compressed_buffer_ = src.compressed_buffer_;
        decompressed_buffer_ = src.decompressed_buffer_;
        decoded_position_ = src.decoded_position_;
        discard_amount_ = src.discard_amount_;
        fp_ = src.fp_;
        if (src.fp_)
          src.fp_ = nullptr;
        put_back_size_ = src.put_back_size_;
        lzma_index_ = src.lzma_index_;
        if (src.lzma_index_)
          src.lzma_index_ = nullptr;
        lzma_res_ = src.lzma_res_;
        at_block_boundary_ = src.at_block_boundary_;
      }

      void replenish_compressed_buffer()
      {
        lzma_block_decoder_.next_in = compressed_buffer_.data();
        lzma_block_decoder_.avail_in = fread(compressed_buffer_.data(), 1, compressed_buffer_.size(), fp_);
      }

      bool init_index()
      {
        if (!fp_)
          return false;

        if (fseek(fp_, -12, SEEK_END) || !fread(stream_footer_.data(), 12, 1, fp_))
          return false;

        if (lzma_stream_footer_decode(&stream_footer_flags_, stream_footer_.data()) != LZMA_OK)
          return false;

        std::vector<std::uint8_t> index_raw(stream_footer_flags_.backward_size);
        if (fseek(fp_, -(stream_footer_flags_.backward_size + 12), SEEK_END) || !fread(index_raw.data(), index_raw.size(), 1, fp_))
          return false;

        std::uint64_t memlimit = UINT64_MAX;
        size_t in_pos = 0;
        auto res = lzma_index_buffer_decode(&lzma_index_, &memlimit, nullptr, index_raw.data(), &in_pos, index_raw.size());
        if (res != LZMA_OK)
          return false;

        lzma_index_iter_init(&lzma_index_itr_, lzma_index_);

        return true;
      }

    private:
      lzma_stream_flags stream_header_flags_;
      lzma_stream_flags stream_footer_flags_;
      lzma_stream lzma_block_decoder_;
      lzma_block lzma_block_;
      std::array<lzma_filter, LZMA_FILTERS_MAX + 1> lzma_block_filters_buf_;
      lzma_index_iter lzma_index_itr_;
      std::array<std::uint8_t, LZMA_STREAM_HEADER_SIZE> stream_header_;
      std::array<std::uint8_t, LZMA_STREAM_HEADER_SIZE> stream_footer_;
      std::array<std::uint8_t, (BUFSIZ >= LZMA_BLOCK_HEADER_SIZE_MAX ? BUFSIZ : LZMA_BLOCK_HEADER_SIZE_MAX)> compressed_buffer_;
      std::array<std::uint8_t, (BUFSIZ >= LZMA_BLOCK_HEADER_SIZE_MAX ? BUFSIZ : LZMA_BLOCK_HEADER_SIZE_MAX)> decompressed_buffer_;
      std::uint64_t decoded_position_;
      std::uint64_t discard_amount_;
      FILE* fp_;
      std::size_t put_back_size_;
      lzma_index* lzma_index_;
      lzma_ret lzma_res_;
      bool at_block_boundary_;
    };

    class obuf : public std::streambuf
    {
    public:
      obuf(FILE* fp)
        :
        fp_(fp)
      {
        if (!fp_)
        {
          char* end = ((char*) decompressed_buffer_.data()) + decompressed_buffer_.size();
          setp(end, end);
        }
        else
        {
          lzma_stream_encoder_ = LZMA_STREAM_INIT;

          lzma_res_ = lzma_easy_encoder(&lzma_stream_encoder_, LZMA_PRESET_DEFAULT, LZMA_CHECK_CRC64);
          if (lzma_res_ != LZMA_OK)
          {
            // TODO: handle error.
          }

          lzma_stream_encoder_.next_out = compressed_buffer_.data();
          lzma_stream_encoder_.avail_out = compressed_buffer_.size();

          char* end = ((char*) decompressed_buffer_.data()) + decompressed_buffer_.size();
          setp((char*) decompressed_buffer_.data(), end);
        }
      }

      obuf(const std::string& file_path) : obuf(fopen(file_path.c_str(), "wb")) {}

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      obuf(obuf&& src)
        :
        std::streambuf(std::move(src))
      {
        this->move(std::move(src));
      }

      obuf& operator=(obuf&& src)
      {
        if (&src != this)
        {
          std::streambuf::operator=(std::move(src));
          this->close();
          this->move(std::move(src));
        }

        return *this;
      }
#endif
      virtual ~obuf()
      {
        this->close();
      }

    protected:

      virtual int overflow(int c)
      {
        if (!fp_)
          return traits_type::eof();

        if ((epptr() - pptr()) > 0)
        {
          assert(!"Put buffer not empty, this should never happen");
          this->sputc(static_cast<char>(0xFF & c));
        }
        else
        {
          lzma_stream_encoder_.next_in = decompressed_buffer_.data();
          lzma_stream_encoder_.avail_in = decompressed_buffer_.size();
          while (lzma_res_ == LZMA_OK && lzma_stream_encoder_.avail_in > 0)
          {
            lzma_res_ = lzma_code(&lzma_stream_encoder_, LZMA_RUN);
            if (lzma_stream_encoder_.avail_out == 0 || lzma_res_ == LZMA_STREAM_END)
            {
              if (!fwrite(compressed_buffer_.data(), compressed_buffer_.size() - lzma_stream_encoder_.avail_out, 1, fp_))
              {
                // TODO: handle error.
                return traits_type::eof();
              }
              lzma_stream_encoder_.next_out = compressed_buffer_.data();
              lzma_stream_encoder_.avail_out = compressed_buffer_.size();
            }
          }

          if (lzma_res_ == LZMA_STREAM_END)
            lzma_res_ = LZMA_OK;

          assert(lzma_stream_encoder_.avail_in == 0);
          decompressed_buffer_[0] = reinterpret_cast<unsigned char&>(c);
          setp((char*) decompressed_buffer_.data() + 1, (char*) decompressed_buffer_.data() + decompressed_buffer_.size());
        }

        return (lzma_res_ == LZMA_OK ? traits_type::to_int_type(c) : traits_type::eof());
      }

      virtual int sync()
      {
        if (!fp_)
          return -1;

        lzma_stream_encoder_.next_in = decompressed_buffer_.data();
        lzma_stream_encoder_.avail_in = decompressed_buffer_.size() - (epptr() - pptr());
        if (lzma_stream_encoder_.avail_in)
        {
          while (lzma_res_ == LZMA_OK)
          {
            lzma_res_ = lzma_code(&lzma_stream_encoder_, LZMA_FULL_FLUSH);
            if (lzma_stream_encoder_.avail_out == 0 || (lzma_res_ == LZMA_STREAM_END && compressed_buffer_.size() != lzma_stream_encoder_.avail_out))
            {
              if (!fwrite(compressed_buffer_.data(), compressed_buffer_.size() - lzma_stream_encoder_.avail_out, 1, fp_))
              {
                // TODO: handle error.
                return -1;
              }
              lzma_stream_encoder_.next_out = compressed_buffer_.data();
              lzma_stream_encoder_.avail_out = compressed_buffer_.size();
            }
          }

          if (lzma_res_ == LZMA_STREAM_END)
            lzma_res_ = LZMA_OK;

          if (lzma_res_ != LZMA_OK)
            return -1;

          assert(lzma_stream_encoder_.avail_in == 0);
          setp((char*) decompressed_buffer_.data(), (char*) decompressed_buffer_.data() + decompressed_buffer_.size());
        }

        return 0;
      }

    private:

      void move(obuf&& src)
      {
        compressed_buffer_ = src.compressed_buffer_;
        decompressed_buffer_ = src.decompressed_buffer_;
        lzma_stream_encoder_ = src.lzma_stream_encoder_;
        if (src.lzma_stream_encoder_.internal)
          src.lzma_stream_encoder_.internal = nullptr;
        fp_ = src.fp_;
        if (src.fp_)
          src.fp_ = nullptr;
        lzma_res_ = src.lzma_res_;
      }

      void close()
      {
        if (lzma_stream_encoder_.internal)
        {
          lzma_stream_encoder_.next_in = decompressed_buffer_.data();
          lzma_stream_encoder_.avail_in = decompressed_buffer_.size() - (epptr() - pptr());
          while (lzma_res_ == LZMA_OK)
          {
            lzma_res_ = lzma_code(&lzma_stream_encoder_, LZMA_FINISH);
            if (lzma_stream_encoder_.avail_out == 0 || (lzma_res_ == LZMA_STREAM_END && compressed_buffer_.size() != lzma_stream_encoder_.avail_out))
            {
              if (!fwrite(compressed_buffer_.data(), compressed_buffer_.size() - lzma_stream_encoder_.avail_out, 1, fp_))
              {
                break;
              }
              lzma_stream_encoder_.next_out = compressed_buffer_.data();
              lzma_stream_encoder_.avail_out = compressed_buffer_.size();
            }
          }
          lzma_end(&lzma_stream_encoder_);
        }

        if (fp_)
          fclose(fp_);
      }

    private:
      static const std::size_t default_block_size = 1024; //4 * 1024 * 1024;
      std::array<std::uint8_t, (default_block_size >= LZMA_BLOCK_HEADER_SIZE_MAX ? default_block_size : LZMA_BLOCK_HEADER_SIZE_MAX)> compressed_buffer_;
      std::array<std::uint8_t, (default_block_size >= LZMA_BLOCK_HEADER_SIZE_MAX ? default_block_size : LZMA_BLOCK_HEADER_SIZE_MAX)> decompressed_buffer_;
      lzma_stream lzma_stream_encoder_;
      FILE* fp_;
      lzma_ret lzma_res_;
    };

    class istream : public std::istream
    {
    public:
      istream(const std::string& file_path)
        :
        std::istream(&sbuf_),
        sbuf_(file_path)
      {
      }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      istream(istream&& src)
        :
        std::istream(&sbuf_),
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
      ::shrinkwrap::xz::ibuf sbuf_;
    };

    class ostream : public std::ostream
    {
    public:
      ostream(const std::string& file_path)
        :
        std::ostream(&sbuf_),
        sbuf_(file_path)
      {
      }

#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
      ostream(ostream&& src)
        :
        std::ostream(&sbuf_),
        sbuf_(std::move(src.sbuf_))
      {
      }

      ostream& operator=(ostream&& src)
      {
        if (&src != this)
        {
          std::ostream::operator=(std::move(src));
          sbuf_ = std::move(src.sbuf_);
        }
        return *this;
      }
#endif
    private:
      ::shrinkwrap::xz::obuf sbuf_;
    };
  }
}

#endif //SHRINKWRAP_XZ_HPP