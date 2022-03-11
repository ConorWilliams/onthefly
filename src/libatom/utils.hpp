#pragma once

#include "Eigen/Core"

namespace otf {

  // Compile time constants

#ifndef LIBATOM_SPATIAL_DIMS
  inline constexpr int spatial_dims = 3;
#else
  inline constexpr int spatial_dims = LIBATOM_SPATIAL_DIMS;
#  undef LIBATOM_SPATIAL_DIMS
#endif

  // Types aliases

  template <typename T> using Vec = Eigen::Vector<T, spatial_dims>;
  template <typename T> using Mat = Eigen::Matrix<T, spatial_dims, spatial_dims>;

  template <typename T> using VecN = Eigen::Matrix<T, spatial_dims, Eigen::Dynamic>;

}  // namespace otf