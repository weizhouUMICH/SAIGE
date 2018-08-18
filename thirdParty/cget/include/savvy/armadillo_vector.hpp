/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_ARMADILLO_VECTOR_HPP
#define LIBSAVVY_ARMADILLO_VECTOR_HPP

#include <armadillo>
#include <cstddef>

namespace savvy
{
  namespace armadillo
  {
    template <typename T>
    class sparse_vector : public arma::SpCol<T>
    {
    public:
      typedef T value_type;
      using arma::SpCol<T>::SpCol;

      void resize(std::size_t sz)
      {
        arma::SpCol<T>::resize(sz, 1);
      }
    };

    template <typename T>
    class dense_vector : public arma::Col<T>
    {
    public:
      typedef T value_type;
      using arma::Col<T>::Col;

      void resize(std::size_t sz)
      {
        std::size_t before_size = arma::Col<T>::size();
        arma::Col<T>::resize(sz);
        if (arma::Col<T>::size() > before_size)
        {
          for (std::size_t i = before_size; i < arma::Col<T>::size(); ++i)
            (*this)[i] = value_type();
        }
      }
    };
  }
}

#endif //LIBSAVVY_ARMADILLO_VECTOR_HPP