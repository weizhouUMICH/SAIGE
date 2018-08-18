/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_ALLELE_STATUS_HPP
#define LIBSAVVY_ALLELE_STATUS_HPP

#include <cstdint>

namespace savvy
{
  enum class allele_status : std::int8_t { is_missing = -1, has_ref, has_alt };
}

#endif //LIBSAVVY_ALLELE_STATUS_HPP
