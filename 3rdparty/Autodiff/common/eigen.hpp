#pragma once

#include "../../Eigen/Core"

#include "vectortraits.hpp"

//=====================================================================================================================
//
// EIGEN MACROS FOR CREATING NEW TYPE ALIASES
//
//=====================================================================================================================

#define AUTODIFF_DEFINE_EIGEN_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)                 \
    using Array##SizeSuffix##SizeSuffix##TypeSuffix = Eigen::Array<Type, Size, Size>;      \
    using Array##SizeSuffix##TypeSuffix = Eigen::Array<Type, Size, 1>;                     \
    using Matrix##SizeSuffix##TypeSuffix = Eigen::Matrix<Type, Size, Size, 0, Size, Size>; \
    using Vector##SizeSuffix##TypeSuffix = Eigen::Matrix<Type, Size, 1, 0, Size, 1>;       \
    using RowVector##SizeSuffix##TypeSuffix = Eigen::Matrix<Type, 1, Size, 1, 1, Size>;

#define AUTODIFF_DEFINE_EIGEN_FIXED_TYPEDEFS(Type, TypeSuffix, Size)                \
    using Array##Size##X##TypeSuffix = Eigen::Array<Type, Size, -1>;                \
    using Array##X##Size##TypeSuffix = Eigen::Array<Type, -1, Size>;                \
    using Matrix##Size##X##TypeSuffix = Eigen::Matrix<Type, Size, -1, 0, Size, -1>; \
    using Matrix##X##Size##TypeSuffix = Eigen::Matrix<Type, -1, Size, 0, -1, Size>;

#define AUTODIFF_DEFINE_EIGEN_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
    AUTODIFF_DEFINE_EIGEN_TYPEDEFS(Type, TypeSuffix, 2, 2)         \
    AUTODIFF_DEFINE_EIGEN_TYPEDEFS(Type, TypeSuffix, 3, 3)         \
    AUTODIFF_DEFINE_EIGEN_TYPEDEFS(Type, TypeSuffix, 4, 4)         \
    AUTODIFF_DEFINE_EIGEN_TYPEDEFS(Type, TypeSuffix, -1, X)        \
    AUTODIFF_DEFINE_EIGEN_FIXED_TYPEDEFS(Type, TypeSuffix, 2)      \
    AUTODIFF_DEFINE_EIGEN_FIXED_TYPEDEFS(Type, TypeSuffix, 3)      \
    AUTODIFF_DEFINE_EIGEN_FIXED_TYPEDEFS(Type, TypeSuffix, 4)

namespace autodiff
{
    namespace detail
    {

        //=====================================================================================================================
        //
        // DEFINE VECTOR TRAITS FOR EIGEN TYPES
        //
        //=====================================================================================================================

        template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
        struct VectorTraits<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>
        {
            using ValueType = Scalar;

            template <typename NewValueType>
            using ReplaceValueType = Eigen::Matrix<NewValueType, Rows, Cols, Options, MaxRows, MaxCols>;
        };

        template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
        struct VectorTraits<Eigen::Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>
        {
            using ValueType = Scalar;

            template <typename NewValueType>
            using ReplaceValueType = Eigen::Array<NewValueType, Rows, Cols, Options, MaxRows, MaxCols>;
        };

        template <typename VectorType, int Size>
        struct VectorTraits<Eigen::VectorBlock<VectorType, Size>>
        {
            using ValueType = typename PlainType<VectorType>::Scalar;

            template <typename NewValueType>
            using ReplaceValueType = VectorReplaceValueType<VectorType, NewValueType>;
        };

        // #if EIGEN_VERSION_AT_LEAST(3, 3, 90)

        //         template <typename VectorType, typename IndicesType>
        //         struct VectorTraits<Eigen::IndexedView<VectorType, IndicesType, Eigen::internal::SingleRange>>
        //         {
        //             using ValueType = typename PlainType<VectorType>::Scalar;

        //             template <typename NewValueType>
        //             using ReplaceValueType = VectorReplaceValueType<VectorType, NewValueType>;
        //         };

        //         template <typename VectorType, typename IndicesType>
        //         struct VectorTraits<Eigen::IndexedView<VectorType, Eigen::internal::SingleRange, IndicesType>>
        //         {
        //             using ValueType = typename PlainType<VectorType>::Scalar;

        //             template <typename NewValueType>
        //             using ReplaceValueType = VectorReplaceValueType<VectorType, NewValueType>;
        //         };

        // #endif

        template <typename MatrixType>
        struct VectorTraits<Eigen::Ref<MatrixType>>
        {
            using ValueType = VectorValueType<MatrixType>;

            template <typename NewValueType>
            using ReplaceValueType = VectorReplaceValueType<MatrixType, NewValueType>;
        };

        template <typename VectorType, int MapOptions, typename StrideType>
        struct VectorTraits<Eigen::Map<VectorType, MapOptions, StrideType>>
        {
            using ValueType = VectorValueType<VectorType>;

            template <typename NewValueType>
            using ReplaceValueType = Eigen::Map<VectorReplaceValueType<VectorType, NewValueType>, MapOptions, StrideType>;
        };

        //=====================================================================================================================
        //
        // AUXILIARY TEMPLATE TYPE ALIASES
        //
        //=====================================================================================================================

        template <typename Scalar>
        using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, 0, Eigen::Dynamic, 1>;

        template <typename Scalar>
        using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, Eigen::Dynamic>;

    } // namespace detail
} // namespace autodiff
