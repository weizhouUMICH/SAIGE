/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_UTILITY_HPP
#define LIBSAVVY_UTILITY_HPP

#include <string>
#include <functional>
#include <memory>
#include <vector>
#include <iomanip>
#include <cstdint>
#include <array>
#include <sstream>
#include <cstring>

namespace savvy
{
  namespace detail
  {
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args)
    {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    template<typename Rng>
    static std::array<std::uint8_t, 16> gen_uuid(Rng& rng)
    {
      static_assert(sizeof(typename Rng::result_type) == 8, "gen_uuid requires a 64 bit PRNG");
      // xxxxxxxx-xxxx-4xxx-{8,9,A,B}xxx-xxxxxxxxxxxx
      // https://www.cryptosys.net/pki/uuid-rfc4122.html

      std::array<std::uint8_t, 16> ret;

      std::uint64_t r1 = rng();
      std::uint64_t r2 = rng();

      std::memcpy(ret.data(), &r1, 8);
      std::memcpy(ret.data() + 8, &r2, 8);

      ret[6] = static_cast<std::uint8_t>(ret[6] & 0x0F) | static_cast<std::uint8_t>(0x40);
      ret[8] = static_cast<std::uint8_t>(ret[8] & 0x3F) | static_cast<std::uint8_t>(0x80);

      return ret;
    }

    template<typename Rng>
    std::string gen_uuid_str(Rng& rng)
    {
      std::array<std::uint8_t, 16> tmp = gen_uuid(rng);
      std::stringstream ret;
      ret << std::hex << std::setfill('0');

      std::size_t i = 0;
      for ( ; i < 16; ++i)
      {
        if (i == 4 || i == 6 || i == 8 || i == 10)
          ret << "-";
        ret << std::setw(2) << (unsigned)tmp[i];
      }

      return ret.str();
    }
  }

  std::string savvy_version();

  struct header_value_details
  {
    std::string id;
    std::string type;
    std::string number;
    std::string description;
  };
  header_value_details parse_header_value(std::string header_value);
  std::string parse_header_sub_field(std::string header_value, std::string field_to_parse);

  template <typename T>
  class hds_to_gp
  {
  public:
    static T get_first_prob(const std::vector<T>& hap_probs)
    {
      T ret = (T(1) - hap_probs[0]);
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= (T(1) - hap_probs[i]);
      return ret;
    }

    static T get_prob(const std::vector<T>& hap_probs, std::size_t num_alleles)
    {
      return choose(hap_probs, num_alleles);
    }

    static T get_last_prob(const std::vector<T>& hap_probs)
    {
      T ret = hap_probs[0];
      for (std::size_t i = 1; i < hap_probs.size(); ++i)
        ret *= hap_probs[i];
      return ret;
    }
  private:
    static void choose(const std::vector<T>& input, std::vector<T>& buf, T& output, std::size_t k, std::size_t offset)
    {
      if (k == 0)
      {
        T product = T(1);
        std::size_t j = 0;
        for (std::size_t i = 0; i < input.size(); ++i)
        {
          if (j < buf.size() && i == buf[j])
          {
            product *= input[i];
            ++j;
          }
          else
          {
            product *= (T(1) - input[i]);
          }
        }
        output += product;
        return;
      }

      for (std::size_t i = offset; i <= input.size() - k; ++i)
      {
        buf.push_back(i);
        choose(input, buf, output, k-1, i+1);
        buf.pop_back();
      }
    }

    static T choose(const std::vector<T>& input, std::size_t k)
    {
      T ret = T(0);
      std::vector<T> buf;
      buf.reserve(k);
      choose(input, buf, ret, k, 0);
      return ret;
    }
  };

  namespace detail
  {
#if __cpp_decltype_auto >= 201304
    template<typename F, typename Tuple, std::size_t... S>
    decltype(auto) apply_impl(F&& fn, Tuple&& t, std::index_sequence<S...>)
    {
      return std::forward<F>(fn)(std::get<S>(std::forward<Tuple>(t))...);
    }

    template<typename F, typename Tuple>
    decltype(auto) apply(F&& fn, Tuple&& t)
    {
      std::size_t constexpr tuple_size
        = std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
      return apply_impl(std::forward<F>(fn),
        std::forward<Tuple>(t),
        std::make_index_sequence<tuple_size>());
    }
#endif
  }
}

#endif // LIBSAVVY_UTILITY_HPP
