#ifndef VARIANT_GROUP_ITERATOR_HPP
#define VARIANT_GROUP_ITERATOR_HPP

#include "savvy/reader.hpp"

#include <regex>
#include <tuple>
#include <iterator>
#include <cstddef>
#include <list>

class variant_group_iterator_util
{
public:
  static std::vector<savvy::genomic_region> merge_regions(const std::vector<savvy::genomic_region>& un_merged_regions)
  {
    std::vector<savvy::genomic_region> ret;

    for (auto it = un_merged_regions.begin(); it != un_merged_regions.end(); ++it)
    {
      if (ret.empty() || ret.back().chromosome() != it->chromosome())
      {
        ret.emplace_back(*it);
      }
      else
      {
        std::uint64_t from = std::min(ret.back().from(), it->from());
        std::uint64_t to = std::max(ret.back().to(), it->to());
        ret.back() = savvy::genomic_region(ret.back().chromosome(), from, to);
      }
    }

    return ret;
  }

  // [CHROM]:[POS]_[REF]/[ALT]
  static savvy::site_info marker_id_to_site_info(std::string::const_iterator beg, std::string::const_iterator end)
  {
    auto colon_it = std::find(beg, end, ':');
    std::string chrom(beg, colon_it);
    if (colon_it != end)
    {
      auto underscore_it = std::find(++colon_it, end, '_');
      std::uint64_t pos = static_cast<std::uint64_t>(std::atoll(std::string(colon_it, underscore_it).c_str()));
      if (underscore_it != end)
      {
        auto slash_it = std::find(++underscore_it, end, '/');
        std::string ref(underscore_it, slash_it);
        if (slash_it != end)
        {
          std::string alt(++slash_it, end);
          return savvy::site_info{std::move(chrom), std::uint32_t(pos), std::move(ref), {std::move(alt)}};
        }
      }
    }

    return savvy::site_info{};
  }

  // [CHROM]:[POS]_[REF]/[ALT]
  static savvy::genomic_region site_info_to_region(const savvy::site_info& site)
  {
    std::size_t length = site.ref().size();
    for (auto it = site.alts().begin(); it != site.alts().end(); ++it)
      length = std::max(length, it->size());

    if (length > 0)
    {
      return savvy::genomic_region{site.chromosome(), site.position(), site.position() + length - 1};
    }

    return savvy::genomic_region{""};
  }

  static std::tuple<std::string, std::list<savvy::site_info>> parse_marker_group_line(const std::string& input)
  {
    std::tuple<std::string, std::list<savvy::site_info>> ret;
    auto delim_it = std::find(input.begin(), input.end(), '\t');
    if (delim_it != input.end())
    {
      std::get<0>(ret) = std::string(input.begin(), delim_it);
      ++delim_it;

      std::string::const_iterator next_delim_it;
      while ((next_delim_it = std::find(delim_it, input.end(), '\t')) != input.end())
      {
        std::get<1>(ret).emplace_back(marker_id_to_site_info(delim_it, next_delim_it));
        delim_it = next_delim_it + 1;
      }

      std::get<1>(ret).emplace_back(marker_id_to_site_info(delim_it, input.end()));
    }

    return ret;
  }
};

class variant_group_iterator : public variant_group_iterator_util
{
public:
  typedef variant_group_iterator self_type;
  typedef std::ptrdiff_t difference_type;
  typedef savvy::variant value_type;
  typedef const value_type& reference;
  typedef const value_type *pointer;
  typedef std::input_iterator_tag iterator_category;

  variant_group_iterator& operator=(const variant_group_iterator& src)
  {
    if (&src != this)
    {
      rdr_ = src.rdr_;
      group_id_ = src.group_id_;
      sites_ = src.sites_;
      merged_regions_ = src.merged_regions_;
      region_offset_ = src.region_offset_;
      variant_ = src.variant_;
    }

    return *this;
  }

  variant_group_iterator& operator=(variant_group_iterator&& src)
  {
    if (&src != this)
    {
      rdr_ = src.rdr_;
      src.rdr_ = nullptr;
      group_id_ = std::move(src.group_id_);
      sites_ = std::move(src.sites_);
      merged_regions_ = std::move(src.merged_regions_);
      region_offset_ = src.region_offset_;
      variant_ = std::move(src.variant_);
    }

    return *this;
  }

  variant_group_iterator(const variant_group_iterator& src)
  {
    operator=(src);
  }

  variant_group_iterator(variant_group_iterator&& src)
  {
    operator=(std::move(src));
  }

  variant_group_iterator(savvy::reader& rdr, std::string marker_group_file_line) :
    rdr_(&rdr)
  {
    std::tie(group_id_, sites_) = parse_marker_group_line(marker_group_file_line);
    init();
  }

  variant_group_iterator(savvy::reader& rdr, std::string group_id, std::list<savvy::site_info> sites) :
    rdr_(&rdr),
    group_id_(std::move(group_id)),
    sites_(std::move(sites))
  {
    init();
  }

  variant_group_iterator() :
    rdr_(nullptr)
  {
  }

  self_type& operator++()
  {
    increment();
    return *this;
  }

  void operator++(int)
  { increment(); }

  reference operator*()
  { return variant_; }

  pointer operator->()
  { return &variant_; }

  bool operator==(const self_type& rhs)
  { return (rdr_ == rhs.rdr_); }

  bool operator!=(const self_type& rhs)
  { return (rdr_ != rhs.rdr_); }

  const std::string& group_id() const
  { return group_id_; }

  const std::list<savvy::site_info>& sites() const
  { return sites_; }

private:
  void init()
  {
    std::vector<savvy::genomic_region> un_merged_regions(sites_.size(), {""});
    auto out_it = un_merged_regions.begin();
    for (auto in_it = sites_.begin(); in_it != sites_.end(); ++in_it, ++out_it)
      *out_it = site_info_to_region(*in_it);

    merged_regions_ = merge_regions(un_merged_regions);
    region_offset_ = 0;

    rdr_->reset_bounds(merged_regions_[region_offset_]);

    increment();
  }

  void increment()
  {
    while (region_offset_ < merged_regions_.size())
    {
      while (!sites_.empty() && rdr_->read(variant_))
      {
        while (!sites_.empty())
        {
          if (sites_.front().position() >= variant_.position() || sites_.front().chromosome() != variant_.chromosome())
            break;
          sites_.pop_front();
        }

        if (sites_.empty() || sites_.front().chromosome() != variant_.chromosome())
          break;

        for (auto sites_it = sites_.begin(); sites_it != sites_.end() && sites_.front().position() == sites_it->position(); ++sites_it)
        {
          if (sites_it->alts().empty()) continue;

          for (const std::string& a : variant_.alts())
          {
//              std::string target_id = sites_it->chromosome() + ":" + std::to_string(sites_it->position()) + "_" + sites_it->ref() + "/" + sites_it->alts()[0];
//              std::string current_id = variant_.chromosome() + ":" + std::to_string(variant_.position()) + "_" + variant_.ref() + "/" + a;
            if (
              sites_it->chromosome() == variant_.chromosome() &&
              sites_it->position() == variant_.position() &&
              sites_it->ref() == variant_.ref() &&
              sites_it->alts()[0] == a)
            {
              return;
            }
          }
        }
      }

      ++region_offset_;
      if (region_offset_ < merged_regions_.size())
      {
        rdr_->reset_bounds(merged_regions_[region_offset_]);
        while (sites_.size() && sites_.front().chromosome() != merged_regions_[region_offset_].chromosome())
          sites_.pop_front();
      }
    }

    rdr_ = nullptr;
  }

  std::string join_vector_to_string(const std::vector<std::string>& vec, std::string delim)
  {
    std::string ret;
    for (auto it = vec.begin(); it != vec.end(); ++it)
    {
      if (it != vec.begin())
        ret += delim;
      ret += *it;
    }
    return ret;
  }

private:
  savvy::reader* rdr_;
  std::string group_id_;
  std::list<savvy::site_info> sites_;
  std::vector<savvy::genomic_region> merged_regions_;
  std::size_t region_offset_;
  value_type variant_;
};

#endif // VARIANT_GROUP_ITERATOR_HPP
