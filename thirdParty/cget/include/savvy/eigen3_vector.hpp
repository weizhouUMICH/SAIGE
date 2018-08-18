/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_EIGEN3_VECTOR_HPP
#define LIBSAVVY_EIGEN3_VECTOR_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cstddef>

namespace savvy
{
  namespace eigen3
  {
    template <typename T>
    class sparse_vector : public ::Eigen::SparseVector<T>
    {
    public:
      using ::Eigen::SparseVector<T>::SparseVector;

      T& operator[](std::size_t idx)
      {
        return ::Eigen::SparseVector<T>::coeffRef(idx);
      }

      T operator[](std::size_t idx) const
      {
        return ::Eigen::SparseVector<T>::coeff(idx);
      }
    };


    template <typename T>
    class dense_vector : public ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>
    {
    public:
      using ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::Matrix;

      T& operator[](std::size_t idx)
      {
        return ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::coeffRef(idx);
      }

      const T& operator[](std::size_t idx) const
      {
        return ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::coeffRef(idx);
      }

      void resize(std::size_t sz)
      {
        using value_type = typename ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::value_type;
        std::size_t before_size = ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::size();
        ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::resize(sz);
        if (::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::size() > before_size)
        {
          for (std::size_t i = before_size; i < ::Eigen::Matrix<T, 1, ::Eigen::Dynamic>::size(); ++i)
            (*this)[i] = value_type();
        }
      }
    };
  }
}
#endif //LIBSAVVY_EIGEN3_VECTOR_HPP
