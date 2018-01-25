/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_REGION_HPP
#define LIBSAVVY_REGION_HPP

#include "site_info.hpp"

#include <cstdint>
#include <string>
#include <limits>

namespace savvy
{
  enum class coord_bound : std::uint8_t
  {
    any = 0,
    all,
    left,
    right
  };

  class region
  {
  public:
    region(const std::string& chromosome, std::uint64_t from = 1, std::uint64_t to = std::numeric_limits<std::uint64_t>::max()) :
      chromosome_(chromosome),
      from_(from),
      to_(to)
    {
    }
    const std::string& chromosome() const { return chromosome_; }
    std::uint64_t from() const { return from_; }
    std::uint64_t to() const { return to_; }
  private:
    std::string chromosome_;
    std::uint64_t from_;
    std::uint64_t to_;
  };

  namespace detail
  {
    struct any_coordinate_within_region
    {
      static bool compare(const site_info& var, const region& reg)
      {
        return (var.position() <= reg.to() && (var.position() + std::max(var.ref().size(), var.alt().size()) - 1) >= reg.from() && var.chromosome() == reg.chromosome());
      }
    };

    struct all_coordinates_within_region
    {
      static bool compare(const site_info& var, const region& reg)
      {
        return (var.position() >= reg.from() && (var.position() + std::max(var.ref().size(), var.alt().size()) - 1) <= reg.to() && var.chromosome() == reg.chromosome());
      }
    };

    struct leftmost_coordinate_within_region
    {
      static bool compare(const site_info& var, const region& reg)
      {
        return (var.position() >= reg.from() && var.position() <= reg.to() && var.chromosome() == reg.chromosome());
      }
    };

    struct rightmost_coordinate_within_region
    {
      static bool compare(const site_info& var, const region& reg)
      {
        std::uint64_t right = (var.position() + std::max(var.ref().size(), var.alt().size()) - 1);
        return (right >= reg.from() && right <= reg.to() && var.chromosome() == reg.chromosome());
      }
    };
  }

  bool region_compare(coord_bound bounding_type, const site_info& var, const region& reg);
}

#endif //LIBSAVVY_REGION_HPP
