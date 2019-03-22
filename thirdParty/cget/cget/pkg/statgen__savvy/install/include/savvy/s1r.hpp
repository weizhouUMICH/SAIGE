/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_S1R_INDEX_READER_HPP
#define LIBSAVVY_S1R_INDEX_READER_HPP

#include "portable_endian.hpp"
#include "region.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cerrno>
#include <map>
#include <cassert>
#include <array>
#include <cstring>
#include <tuple>

namespace savvy
{
  namespace s1r
  {
    class internal_entry
    {
    public:
      std::uint32_t region_start() const
      {
        std::uint32_t ret = be32toh(region_start_);
        return ret;
      }
      std::uint32_t region_length() const
      {
        std::uint32_t ret = be32toh(region_length_);
        return ret;
      }
      std::uint32_t region_end() const
      {
        std::uint32_t ret = this->region_start() + this->region_length();
        return ret;
      }

      internal_entry() :
        region_start_(0),
        region_length_(0)
      { }

      internal_entry(std::uint32_t beg, std::uint32_t end) :
        region_start_(htobe32(beg)),
        region_length_(end > beg ? htobe32(end - beg) : 0)
      { }
    private:
      std::uint32_t region_start_;
      std::uint32_t region_length_;
    };

    static_assert(sizeof(internal_entry) == 8, "Alignment issue (internal_entry)");

    class entry : public internal_entry
    {
    public:
      std::uint64_t value() const
      {
        return be64toh(value_);
      }

      entry() :
        value_(0)
      { }

      entry(std::uint32_t beg, std::uint32_t end, std::uint64_t value) :
        internal_entry(beg, end),
        value_(htobe64(value))
      { }
    private:
      std::uint64_t value_;
    };
    static_assert(sizeof(entry) == 16, "Alignment issue (entry)");

    namespace detail
    {
      inline std::uint64_t ceil_divide(std::uint64_t x, std::uint64_t y)
      {
        return (x + y - 1) / y;
      }

      inline std::uint16_t ceil_divide(std::uint16_t x, std::uint16_t y)
      {
        return (x + y - std::uint16_t(1)) / y;
      }

      inline std::size_t entries_per_leaf_node(std::size_t block_size)
      {
        return std::uint16_t(block_size / sizeof(entry));
      }

      inline std::size_t entries_per_internal_node(std::size_t block_size)
      {
        return std::uint16_t(block_size / sizeof(internal_entry));
      }
    }

    class tree_base
    {
    public:
      struct node_position
      {
        node_position(std::uint64_t lev, std::uint64_t node_off) : level(lev), node_offset(node_off) {}
        std::uint64_t level;
        std::uint64_t node_offset;
      };

      struct tree_position : public node_position
      {
        tree_position(std::uint64_t lev, std::uint64_t node_off, std::uint64_t entry_off) : node_position(lev, node_off), entry_offset(entry_off) {}
        tree_position(node_position node_pos, std::uint64_t entry_off) : node_position(node_pos), entry_offset(entry_off) {}
        bool operator==(const tree_position& other) const { return level == other.level && node_offset == other.node_offset && entry_offset == other.entry_offset; }
        bool operator!=(const tree_position& other) const { return !(*this == other); }
        std::uint64_t entry_offset;
      };

      tree_base(std::uint8_t block_size_in_kib, std::uint64_t block_offset, std::uint64_t entry_count)
      {
        //const std::uint16_t header_data_size = 7 + 2 + 8;
        block_size_ = 1024u * (std::uint32_t(block_size_in_kib) + 1);
        block_offset_ = block_offset; //ceil_divide(header_data_size, block_size_);
        entry_count_ = entry_count;

        std::size_t node_count = 1;
        entry_counts_per_level_.clear();
        this->entry_counts_per_level_.insert(this->entry_counts_per_level_.begin(), entry_count_);
        for (std::uint64_t nodes_at_current_level = detail::ceil_divide(entry_count_, (std::uint64_t) this->entries_per_leaf_node());
             nodes_at_current_level > 1;
             nodes_at_current_level = detail::ceil_divide(nodes_at_current_level, (std::uint64_t) this->entries_per_internal_node()))
        {
          node_count += nodes_at_current_level;
          this->entry_counts_per_level_.insert(this->entry_counts_per_level_.begin(), nodes_at_current_level);
        }

        root_block_end_pos_ = (block_offset + node_count) * block_size_;
      }

      std::uint64_t entry_count() const { return entry_count_; }
      std::uint64_t bucket_size() const { return block_size_; }
      std::uint64_t entry_count_at_level(std::size_t level) const { return entry_counts_per_level_[level]; }
      std::uint16_t entries_per_leaf_node() const { return std::uint16_t(block_size_ / sizeof(entry)); }
      std::uint16_t entries_per_internal_node() const { return std::uint16_t(block_size_ / sizeof(internal_entry)); }

      std::uint64_t calculate_node_size(node_position input) const
      {
        std::uint64_t ret = 0;

        if (input.level < entry_counts_per_level_.size())
        {
          std::uint64_t entries_per_node = (input.level + 1 ==  entry_counts_per_level_.size() ? this->entries_per_leaf_node() : this->entries_per_internal_node());

          if (input.node_offset < entry_counts_per_level_[input.level] / entries_per_node)
            ret = entries_per_node;
          else
            ret = entry_counts_per_level_[input.level] % entries_per_node;
        }

        return ret;
      }

      tree_position calculate_parent_position(const node_position& input) const
      {
        return tree_position(input.level - 1, input.node_offset / this->entries_per_internal_node(), input.node_offset % this->entries_per_internal_node());
      }

      node_position calculate_child_position(const tree_position& input) const
      {
        return node_position(input.level + 1, this->entries_per_internal_node() * input.node_offset + input.entry_offset);
      }

      std::streampos calculate_file_position(const node_position& input) const
      {
        std::streampos ret = root_block_end_pos_ - std::streampos(block_size_);

        if (input.level > 0)
        {
          for (auto it = entry_counts_per_level_.begin(); it != entry_counts_per_level_.end(); ++it)
          {
            if (1 + std::distance(entry_counts_per_level_.begin(), it) < input.level)
            {
              ret -= (*it * block_size_);
            }
            else
            {
              ret -= (*it * block_size_ - input.node_offset * block_size_);
              break;
            }
          }
        }

        assert(ret >= std::streampos(block_offset_));

        return ret;
      }


      std::uint64_t tree_height() const
      {
        return entry_counts_per_level_.size();
      }


    private:
      std::streampos root_block_end_pos_;
      std::uint64_t block_offset_;
      std::uint32_t block_size_;
      std::vector<std::uint64_t> entry_counts_per_level_;
      std::uint64_t entry_count_;
    };

    class tree_reader : public tree_base
    {
    public:
      class query
      {
      public:
        class iterator
        {
        public:
          typedef iterator self_type;
          typedef std::ptrdiff_t difference_type;
          typedef entry value_type;
          typedef const value_type& reference;
          typedef const value_type* pointer;
          typedef std::bidirectional_iterator_tag iterator_category;

          iterator(tree_reader& rdr, std::istream& ifs, std::uint64_t beg, std::uint64_t end, tree_position pos = {0, 0, 0}) :
            reader_(&rdr),
            ifs_(&ifs),
            beg_(beg),
            end_(end),
            leaf_node_(reader_->entries_per_leaf_node()),
            position_(pos)
          {
            if (position_ != reader_->end_tree_position())
            {
              const std::uint64_t leaf_level = reader_->tree_height() - 1;
              ifs_->seekg(reader_->calculate_file_position(position_));

              if (position_.level == leaf_level)
                ifs_->read((char*) leaf_node_.data(), reader_->bucket_size());
              else
              {
                traversal_chain_.emplace(reader_->entries_per_internal_node());
                ifs_->read((char*) (traversal_chain_.top().data()), reader_->bucket_size());
              }

              traverse_right();
            }
          }

          void traverse_right()
          {
            const std::uint64_t leaf_level = reader_->tree_height() - 1;
            bool leaf_entry_found = false;
            while (!leaf_entry_found && position_ != reader_->end_tree_position())
            {
              if (position_.level == leaf_level)
              {
                auto entry_it = leaf_node_.begin() + position_.entry_offset;
                const auto entry_end_it = leaf_node_.begin() + reader_->calculate_node_size(position_);
                for (; entry_it != entry_end_it; ++entry_it)
                {
                  if (entry_it->region_start() <= end_ && entry_it->region_end() >= beg_)
                    break;
                }

                if (entry_it == entry_end_it)
                {
                  if (position_.level > 0)
                  {
                    position_ = reader_->calculate_parent_position(position_);
                    position_ = tree_position(position_, position_.entry_offset + 1);
                  }
                  else
                  {
                    position_ = reader_->end_tree_position();
                  }
                }
                else
                {
                  position_ = tree_position(position_.level, position_.node_offset, (std::uint64_t) std::distance(leaf_node_.begin(), entry_it));
                  leaf_entry_found = true;
                }
              }
              else
              {
                auto entry_it = traversal_chain_.top().begin() + position_.entry_offset;
                const auto entry_end_it = traversal_chain_.top().begin() + reader_->calculate_node_size(position_);
                for (; entry_it != entry_end_it; ++entry_it)
                {
                  if (entry_it->region_start() <= end_ && entry_it->region_end() >= beg_)
                    break;
                }

                if (entry_it == entry_end_it)
                {
                  if (position_.level > 0)
                  {
                    position_ = reader_->calculate_parent_position(position_);
                    position_ = tree_position(position_, position_.entry_offset + 1);
                    traversal_chain_.pop();
                  }
                  else
                  {
                    position_ = reader_->end_tree_position();
                  }
                }
                else
                {
                  position_ = tree_position(reader_->calculate_child_position(tree_position(position_, (std::uint64_t) std::distance(traversal_chain_.top().begin(), entry_it))), 0);

                  ifs_->seekg(reader_->calculate_file_position(position_));
                  if (position_.level == leaf_level)
                    ifs_->read((char*) leaf_node_.data(), reader_->bucket_size());
                  else
                  {
                    traversal_chain_.emplace(reader_->entries_per_internal_node());
                    ifs_->read((char*) (traversal_chain_.top().data()), reader_->bucket_size());
                  }
                }
              }
            }
          }

          self_type& operator++()
          {
            position_ = tree_position(position_, position_.entry_offset + 1);
            traverse_right();
            return *this;
          }

          self_type operator++(int)
          {
            self_type r = *this;
            ++(*this);
            return r;
          }

          reference operator*() { return leaf_node_[position_.entry_offset]; }
          pointer operator->() { return &(leaf_node_[position_.entry_offset]); }
          bool operator==(const self_type& other) { return position_ == other.position_; }
          bool operator!=(const self_type& other) { return position_ != other.position_; }

        private:
          tree_reader* reader_;
          std::istream* ifs_;
          std::uint64_t beg_;
          std::uint64_t end_;
          std::stack<std::vector<internal_entry>> traversal_chain_;
          std::vector<entry> leaf_node_;
          tree_position position_;
        };

        query(tree_reader& rdr, std::istream& ifs, std::uint64_t beg, std::uint64_t end) :
          reader_(&rdr),
          ifs_(&ifs),
          beg_(beg),
          end_(end)
        {
        }

        iterator begin() { return iterator(*reader_, *ifs_, beg_, end_); }
        iterator end() { return iterator(*reader_, *ifs_, beg_, end_, reader_->end_tree_position()); }
      private:
        tree_reader* reader_;
        std::istream* ifs_;
        std::uint64_t beg_;
        std::uint64_t end_;
      };

      tree_reader(std::ifstream& file, std::uint8_t block_size_in_kib, std::uint64_t block_offset, const std::string& name, std::uint64_t entry_count) :
        tree_base(block_size_in_kib, block_offset, entry_count),
        ifs_(file),
        name_(name)
      {
//        std::string version(7, '\0');
//        ifs_.read(&version[0], version.size());
//        std::uint8_t block_size_exponent;
//        ifs_.read((char*)(&block_size_exponent), 1);
//        block_size_ = static_cast<std::uint16_t>(std::pow(8.0, (0x03 & block_size_exponent) + 2));
//        ifs_.read((char*)(&entry_count_), 8);
//        entry_count_ = be64toh(entry_count_);
      }

      bool good() const { return ifs_.good(); }

      std::tuple<std::uint32_t, std::uint32_t> range()
      {
        std::tuple<std::uint32_t, std::uint32_t> ret{std::numeric_limits<std::uint32_t>::max(), 0};

        auto fn = [&ret](const internal_entry& e)
        {
          std::get<0>(ret) = std::min(std::get<0>(ret), e.region_start());
          std::get<1>(ret) = std::max(std::get<1>(ret), e.region_end());
        };

        node_position position{0, 0};

        const std::uint64_t leaf_level = this->tree_height() - 1;
        ifs_.seekg(this->calculate_file_position(position));

        if (position.level == leaf_level)
        {
          std::vector<entry> leaf_node(this->entries_per_leaf_node());
          ifs_.read((char*) leaf_node.data(), this->bucket_size());

          const auto entry_end_it = leaf_node.begin() + this->calculate_node_size(position);
          std::for_each(leaf_node.begin(), entry_end_it, fn);
        }
        else
        {
          std::vector<internal_entry> internal_node(this->entries_per_internal_node());
          ifs_.read((char*) internal_node.data(), this->bucket_size());
          const auto entry_end_it = internal_node.begin() + this->calculate_node_size(position);
          std::for_each(internal_node.begin(), entry_end_it, fn);
        }

        return ret;
      }

      const tree_position end_tree_position()
      {
        return tree_position(0, 0, calculate_node_size(tree_position(0, 0, 0)));
      }

      query create_query(std::uint64_t beg, std::uint64_t end)
      {
        query ret(*this, ifs_, beg, end);
        return ret;
      }

      const std::string& name() const { return name_; }
    private:
      std::ifstream& ifs_;
      std::string name_;
    };

    enum class sort_point : std::uint8_t
    {
      mid = 0,
      beg = 0x10,
      end = 0x01
    };

    class reader
    {
    public:
      typedef reader self_type;


      reader(const std::string& file_path) :
        file_path_(file_path),
        input_file_(file_path, std::ios::binary)
      {
        std::array<char, 26> footer;

        std::uint8_t block_size_byte = 0;
        std::uint64_t block_count = 0;

        input_file_.seekg(-(footer.size()), std::ios::end);
        input_file_.read(footer.data(), footer.size());
        if (!input_file_.good())
        {
          input_file_.setstate(std::ios::badbit);
        }
        else
        {
          std::size_t bytes_parsed = 0;
          std::string version(footer.begin() + (footer.size() - 7), footer.begin() + footer.size());
          bytes_parsed += version.size();

          if (version.substr(0, 3) != "s1r")
          {
            input_file_.setstate(std::ios::badbit);
          }
          else
          {
            std::string uuid(footer.begin() + (footer.size() - bytes_parsed - 16), footer.begin() + (footer.size() - bytes_parsed));
            bytes_parsed += uuid.size();

            std::uint16_t tree_details_size_be;
            std::memcpy((char*)(&tree_details_size_be), footer.data() + (footer.size() - bytes_parsed - sizeof(tree_details_size_be)), sizeof(tree_details_size_be));
            bytes_parsed += sizeof(tree_details_size_be);

            std::uint16_t tree_details_size = be16toh(tree_details_size_be);


            std::memcpy((char*)(&block_size_byte), footer.data() + (footer.size() - bytes_parsed - sizeof(block_size_byte)), sizeof(block_size_byte));
            bytes_parsed += sizeof(block_size_byte);

            const std::uint32_t block_size = 1024u * (std::uint32_t(block_size_byte) + 1);

            struct tree_details
            {
              std::string name;
              std::uint64_t entry_count = 0;
            };

            input_file_.seekg(-(std::int64_t(footer.size()) + tree_details_size), std::ios::end);
            std::vector<tree_details> tree_details_array;

            int tree_details_bytes_left = tree_details_size;
            while (tree_details_bytes_left > 0 && input_file_.good())
            {
              tree_details details;
              std::getline(input_file_, details.name, '\0');

              std::uint64_t entry_count_be = 0;
              input_file_.read((char*)(&entry_count_be), 8);
              details.entry_count = be64toh(entry_count_be);

              tree_details_bytes_left -= (details.name.size() + 1 + 8);


              tree_details_array.emplace_back(std::move(details));
            }

            assert(tree_details_bytes_left == 0);

            trees_.reserve(tree_details_array.size());
            for (auto it = tree_details_array.begin(); it != tree_details_array.end(); ++it)
            {
              trees_.emplace_back(input_file_, block_size_byte, block_count, it->name, it->entry_count);

              for (std::uint64_t nodes_at_current_level = detail::ceil_divide(it->entry_count, (std::uint64_t) detail::entries_per_leaf_node(block_size));
                nodes_at_current_level > 1;
                nodes_at_current_level = detail::ceil_divide(nodes_at_current_level, (std::uint64_t) detail::entries_per_internal_node(block_size)))
              {
                block_count += nodes_at_current_level;
              }

              block_count += 1;
            }
          }
        }

        trees_.emplace_back(input_file_, block_size_byte, block_count, "", 0); // empty tree (end marker).


//        std::uint8_t block_size_exponent;
//        ifs_.read((char*)(&block_size_exponent), 1);
//        block_size_ = static_cast<std::uint16_t>(std::pow(8.0, (0x03 & block_size_exponent) + 2));
//        ifs_.read((char*)(&entry_count_), 8);
//        entry_count_ = be64toh(entry_count_);
      }

      bool good() const { return input_file_.good(); }

      std::vector<std::string> tree_names() const
      {
        std::vector<std::string> ret;
        std::size_t sz = trees_.size();
        if (sz > 0)
        {
          --sz;
          ret.reserve(sz);
          for (auto it = trees_.begin(); it != std::prev(trees_.end()); ++it)
          {
            ret.push_back(it->name());
          }
        }

        return ret;
      }

      std::vector<tree_reader>::iterator trees_begin()
      {
        return trees_.begin();
      }

      std::vector<tree_reader>::iterator trees_end()
      {
        if (trees_.size())
          return trees_.begin() + trees_.size() - 1;
        return trees_.begin();
      }

      class query;
      query create_query(region reg);
      query create_query(std::vector<region> regs);
    private:
      std::string file_path_;
      std::ifstream input_file_;
      std::vector<tree_reader> trees_;
    };

    class reader::query
    {
    public:
      query(std::ifstream& ifs, std::vector<tree_reader>& trees, std::vector<region> regs) :
        ifs_(&ifs),
        regions_{regs}
      {
        for (auto i = regions_.begin(); i != regions_.end(); ++i)
        {
          for (auto j = trees.begin(); j != std::prev(trees.end()); ++j)
          {
            if (i->chromosome() == j->name())
            {
              auto tmp_query = j->create_query(i->from(), i->to());
              auto b = tmp_query.begin();
              auto e = tmp_query.end();
              if (b != e)
                tree_queries_.emplace_back(std::move(tmp_query));
            }
          }
        }
        tree_queries_.emplace_back(trees.back().create_query(0, 0)); // empty tree.
      }

      query(std::ifstream& ifs, std::vector<tree_reader>& trees, region reg) : query(ifs, trees, std::vector<region>({reg})) {}

      class iterator;
      iterator begin();
      iterator end();
    private:
      std::ifstream* ifs_;
      std::vector<tree_reader::query> tree_queries_;
      std::vector<region> regions_;
    };

    class reader::query::iterator
    {
    public:
      typedef iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef entry value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::bidirectional_iterator_tag iterator_category;

      iterator(std::istream& ifs, std::vector<tree_reader::query>::iterator tree_query_it, tree_reader::query::iterator tree_query_beg, tree_reader::query::iterator tree_query_end) :
        tree_it_(tree_query_it),
        tree_query_it_(tree_query_beg),
        tree_query_end_(tree_query_end),
        ifs_(&ifs)
      {

      }

      self_type& operator++()
      {
        ++tree_query_it_;
        if (tree_query_it_ == tree_query_end_)
        {
          ++tree_it_;
          tree_query_it_ = tree_it_->begin();
          tree_query_end_ = tree_it_->end();
        }

        return *this;
      }

      self_type operator++(int)
      {
        self_type r = *this;
        ++(*this);
        return r;
      }

      reference operator*() { return *tree_query_it_; }
      pointer operator->() { return &(*tree_query_it_); }
      bool operator==(const self_type& other) { return tree_it_ == other.tree_it_ && tree_query_it_ == other.tree_query_it_; }
      bool operator!=(const self_type& other) { return tree_it_ != other.tree_it_ || tree_query_it_ != other.tree_query_it_; }

    private:
      std::vector<tree_reader::query>::iterator tree_it_;
      tree_reader::query::iterator tree_query_it_;
      tree_reader::query::iterator tree_query_end_;
      std::istream* ifs_;
    };

    class writer
    {
    public:
      writer(const std::string& file_path, const std::array<std::uint8_t, 16>& uuid, std::uint8_t block_size_in_kib = 4 - 1) :
        file_path_(file_path),
        ofs_(file_path, std::ios::binary),
        uuid_(uuid)
      {
        this->block_size_ = 1024u * (std::uint32_t(block_size_in_kib) + 1);

        current_leaf_node_.reserve(detail::entries_per_leaf_node(block_size_));
      }

      ~writer()
      {
        if (chromosomes_.size())
        {
          this->write_internal_nodes();
        }

        std::uint16_t index_size = 0;
        for (auto it = chromosomes_.begin(); it != chromosomes_.end(); ++it)
        {
          ofs_.write(it->first.c_str(), it->first.size() + 1);

          std::uint64_t entry_count_be = htobe64(it->second);
          ofs_.write((char *) (&entry_count_be), 8);

          index_size += (1 + it->first.size() + 8);
        }

        std::string footer_block(26, '\0');
        std::size_t cur = 0;

        std::uint8_t block_size_in_kib = std::uint8_t(block_size_ / 1024 - 1);
        std::memcpy(&footer_block[cur], (char *) (&block_size_in_kib), 1);
        cur += 1;

        std::uint16_t index_size_be = htobe16(index_size);
        std::memcpy(&footer_block[cur], (char *) (&index_size_be), 2);
        cur += 2;

        std::memset(&footer_block[cur], '\0', 16); // uuid
        cur += 16;

        std::memcpy(&footer_block[cur], "s1r\x00\x01\x00\x00", 7);
        ofs_.write(footer_block.data(), footer_block.size());
      }


      writer& write(const std::string& chrom, const entry& e)
      {
        if (chromosomes_.empty())
          chromosomes_.emplace_back(chrom, 0);

        if (chrom != chromosomes_.back().first)
        {
          this->write_internal_nodes();
          this->chromosomes_.emplace_back(chrom, 0);
        }


        current_leaf_node_.emplace_back(e);
        if (current_leaf_node_.size() == detail::entries_per_leaf_node(block_size_))
        {
          ofs_.write((char *)(current_leaf_node_.data()), block_size_);

          current_leaf_node_.resize(0);
        }

        ++(chromosomes_.back().second);

        return *this;
      }

      void write_internal_nodes()
      {
        if (this->current_leaf_node_.size())
        {
          this->current_leaf_node_.resize(detail::entries_per_leaf_node(block_size_));
          ofs_.write((char *)(current_leaf_node_.data()), block_size_);

          current_leaf_node_.resize(0);
        }

        ofs_.flush();

        std::uint64_t num_leaf_nodes = detail::ceil_divide(this->chromosomes_.back().second, std::uint64_t(detail::entries_per_leaf_node(block_size_)));
        tree_base tree(std::uint8_t(block_size_ / 1024 - 1), std::uint64_t(ofs_.tellp()) - block_size_ * num_leaf_nodes, this->chromosomes_.back().second);

        std::vector<std::pair<std::vector<internal_entry>, tree_base::tree_position>> current_nodes_at_each_internal_level;
        current_nodes_at_each_internal_level.reserve(tree.tree_height() - 1);
        for (std::size_t i = 0; i < tree.tree_height() - 1; ++i)
          current_nodes_at_each_internal_level.emplace_back(std::make_pair(std::vector<internal_entry>(tree.entries_per_internal_node()), tree_base::tree_position(i, 0, 0)));

        std::vector<entry> current_leaf_node(detail::entries_per_leaf_node(block_size_));
        tree_base::tree_position current_leaf_position(tree.tree_height() - 1, 0, 0);

        std::ifstream ifs(file_path_, std::ios::binary);


        ifs.seekg(ofs_.tellp() - std::streamoff(block_size_ * num_leaf_nodes));

        for (std::size_t i = 0; i < num_leaf_nodes && ifs.good(); ++i)
        {
          bool last_leaf_node = (i + 1) == num_leaf_nodes;

          ifs.read((char*)current_leaf_node.data(), block_size_);

          std::uint32_t node_range_min = (std::uint32_t)-1;
          std::uint32_t node_range_max = 0;
          std::size_t current_leaf_node_end = last_leaf_node
            ? (this->chromosomes_.back().second % detail::entries_per_leaf_node(block_size_)) //detail::entries_per_leaf_node(block_size_) - (this->chromosomes_.back().second % detail::entries_per_leaf_node(block_size_))
            : detail::entries_per_leaf_node(block_size_);

          if (!current_leaf_node_end)
          {
            current_leaf_node_end = detail::entries_per_leaf_node(block_size_);
          }

          for (std::size_t i = 0; i < current_leaf_node_end; ++i)
          {
            node_range_min = std::min(node_range_min, current_leaf_node[i].region_start());
            node_range_max = std::max(node_range_max, current_leaf_node[i].region_end());
          }

          for (auto rit = current_nodes_at_each_internal_level.rbegin(); rit != current_nodes_at_each_internal_level.rend(); ++rit)
          {
            rit->first[rit->second.entry_offset] = internal_entry(node_range_min, node_range_max);

            if (rit->second.entry_offset + 1 == rit->first.size() || last_leaf_node)
            {
              ofs_.seekp(tree.calculate_file_position(rit->second));
              ofs_.write((char*)(rit->first.data()), block_size_);

              node_range_min = (std::uint32_t)-1;
              node_range_max = 0;
              for (std::size_t i = 0; i <= rit->second.entry_offset; ++i)
              {
                node_range_min = std::min(node_range_min, rit->first[i].region_start());
                node_range_max = std::max(node_range_max, rit->first[i].region_end());
              }

              rit->first = std::vector<internal_entry>(tree.entries_per_internal_node());
              ++(rit->second.node_offset);
              rit->second.entry_offset = 0;
            }
            else
            {
              ++(rit->second.entry_offset);
              break;
            }
          }

        }

        if (!ifs.good())
          ofs_.setstate(std::ios::badbit);
      }






//      writer(const std::string& file_path, std::map<std::string, std::vector<tree_base::entry>>& chrom_data, std::uint8_t block_size_in_kib) :
//        ofs_(file_path, std::ios::binary)
//      {
//        const std::uint32_t block_size = 1024u * (std::uint32_t(block_size_in_kib) + 1);
//        //std::uint64_t header_size = 7 + 16 + 2;
//
//        for (auto it = chrom_data.begin(); it != chrom_data.end(); ++it)
//        {
//          auto beg = it->second.begin();
//          auto end = it->second.end();
//          std::sort(beg, end, [](const tree_reader::entry& a, const tree_reader::entry& b)
//          {
//            double mid_a = static_cast<double>(a.region_start()) + (static_cast<double>(a.region_length()) / 2.0);
//            double mid_b = static_cast<double>(b.region_start()) + (static_cast<double>(b.region_length()) / 2.0);
//
//            return mid_a < mid_b;
//          });
//        }
//
//
//
//        //std::uint64_t block_offset = header_block.size() / std::uint64_t(1024u * (std::uint32_t(block_size_in_kib) + 1));
//        std::size_t block_offset = 0;
//        for (auto it = chrom_data.begin(); it != chrom_data.end(); ++it)
//        {
//          auto beg = it->second.begin();
//          auto end = it->second.end();
//
//          entry_count_ = (std::uint64_t) std::distance(beg, end);
//
//          tree_base::init(block_size_in_kib, block_offset, it->first);
//
//          std::vector<std::pair<std::vector<internal_entry>, tree_position>> current_nodes_at_each_internal_level;
//          current_nodes_at_each_internal_level.reserve(this->tree_height() - 1);
//          for (std::size_t i = 0; i < this->tree_height() - 1; ++i)
//            current_nodes_at_each_internal_level.emplace_back(std::make_pair(std::vector<internal_entry>(this->entries_per_internal_node()), tree_position(i, 0, 0)));
//
//          std::vector<entry> current_leaf_node(this->entries_per_leaf_node());
//          tree_position current_leaf_position(this->tree_height() - 1, 0, 0);
//
//          std::size_t node_counter = 0;
//          for (auto entry_it = beg; entry_it != end; ++entry_it)
//          {
//            const bool last_entry = entry_it + 1 == end;
//            const tree_reader::entry& e = *entry_it;
//
//            current_leaf_node[current_leaf_position.entry_offset] = e;
//
//            if (current_leaf_position.entry_offset + 1 == current_leaf_node.size() || last_entry)
//            {
//              ofs_.seekp(this->calculate_file_position(current_leaf_position));
//              ofs_.write((char*)(current_leaf_node.data()), block_size_);
//              ++node_counter;
//
//              std::uint32_t node_range_min = (std::uint32_t)-1;
//              std::uint32_t node_range_max = 0;
//              for (std::size_t i = 0; i <= current_leaf_position.entry_offset; ++i)
//              {
//                node_range_min = std::min(node_range_min, current_leaf_node[i].region_start());
//                node_range_max = std::max(node_range_max, current_leaf_node[i].region_end());
//              }
//
//              for (auto rit = current_nodes_at_each_internal_level.rbegin(); rit != current_nodes_at_each_internal_level.rend(); ++rit)
//              {
//                rit->first[rit->second.entry_offset] = internal_entry(node_range_min, node_range_max);
//
//
//
//                if (rit->second.entry_offset + 1 == rit->first.size() || last_entry)
//                {
//                  ofs_.seekp(this->calculate_file_position(rit->second));
//                  ofs_.write((char*)(rit->first.data()), block_size_);
//                  ++node_counter;
//
//                  node_range_min = (std::uint32_t)-1;
//                  node_range_max = 0;
//                  for (std::size_t i = 0; i <= rit->second.entry_offset; ++i)
//                  {
//                    node_range_min = std::min(node_range_min, rit->first[i].region_start());
//                    node_range_max = std::max(node_range_max, rit->first[i].region_end());
//                  }
//
//                  rit->first = std::vector<internal_entry>(this->entries_per_internal_node());
//                  ++(rit->second.node_offset);
//                  rit->second.entry_offset = 0;
//                }
//                else
//                {
//                  ++(rit->second.entry_offset);
//                  break;
//                }
//              }
//
//              current_leaf_node = std::vector<entry>(this->entries_per_leaf_node());
//              ++(current_leaf_position.node_offset);
//              current_leaf_position.entry_offset = 0;
//            }
//            else
//            {
//              ++(current_leaf_position.entry_offset);
//            }
//          }
//
//          block_offset += node_counter;
//        }
//
//        std::uint16_t index_size = 0;
//        for (auto it = chrom_data.begin(); it != chrom_data.end(); ++it)
//        {
//          ofs_.write(it->first.c_str(), it->first.size() + 1);
//
//          std::uint64_t entry_count_be = htobe64((std::uint64_t)std::distance(it->second.begin(), it->second.end()));
//          ofs_.write((char*)(&entry_count_be), 8);
//
//          index_size += (1 + it->first.size() + 8);
//        }
//
//        std::string footer_block(26, '\0');
//        std::size_t cur = 0;
//
//        std::memcpy(&footer_block[cur], (char*)(&block_size_in_kib), 1);
//        cur += 1;
//
//        std::uint16_t index_size_be = htobe16(index_size);
//        std::memcpy(&footer_block[cur], (char*)(&index_size_be), 2);
//        cur += 2;
//
//        std::memset(&footer_block[cur], '\0', 16); // uuid
//        cur += 16;
//
//        std::memcpy(&footer_block[cur], "s1r\x00\x01\x00\x00", 7);
//        ofs_.write(footer_block.data(), footer_block.size());
//      }

      explicit operator bool() const { return ofs_.good(); }
      bool good() { return ofs_.good(); }

      bool flush()
      {
        return ofs_.flush().good();
      }
    private:
      std::string file_path_;
      std::ofstream ofs_;
      const std::array<std::uint8_t, 16> uuid_;
      std::uint32_t block_size_;
      std::vector<entry> current_leaf_node_;
      std::vector<std::pair<std::string, std::uint64_t>> chromosomes_;
    };

//    inline bool create_file(const std::string& file_path, std::map<std::string, std::vector<s1r::tree_base::entry>>& entries, std::uint8_t block_size_in_kib)
//    {
//      writer w(file_path, block_size_in_kib);
//      return w.flush();
//    }


    inline reader::query reader::create_query(region reg)
    {
      query ret(input_file_, trees_, reg);
      return ret;
    }

    inline reader::query reader::create_query(std::vector<region> regs)
    {
      query ret(input_file_, trees_, regs);
      return ret;
    }

    inline reader::query::iterator reader::query::begin()
    {
      return reader::query::iterator(*ifs_, tree_queries_.begin(), tree_queries_.begin()->begin(), tree_queries_.begin()->end());
    }

    inline reader::query::iterator reader::query::end()
    {
      return reader::query::iterator(*ifs_, std::prev(tree_queries_.end()), std::prev(tree_queries_.end())->begin(), std::prev(tree_queries_.end())->end());
    }
  }
}

#endif //LIBSAVVY_S1R_INDEX_READER_HPP