/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_UBLAS_VECTOR_HPP
#define LIBSAVVY_UBLAS_VECTOR_HPP

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace savvy
{
  namespace ublas
  {
    template <typename T>
    using sparse_vector = boost::numeric::ublas::compressed_vector<T>;

    template <typename T>
    using dense_vector = boost::numeric::ublas::vector<T>;
  }
}

#endif //LIBSAVVY_UBLAS_VECTOR_HPP