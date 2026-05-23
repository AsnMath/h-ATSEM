#pragma once

// C++ includes
#include <vector>

// autodiff includes
#include "meta.hpp"

namespace autodiff
{
    namespace detail
    {

        /// An auxiliary template type to indicate VectorTraits has not been defined for a type.
        template <typename V>
        struct VectorTraitsNotDefinedFor
        {
        };

        /// An auxiliary template type to indicate VectorTraits::ReplaceValueType is not supported for a type.
        template <typename V>
        struct VectorReplaceValueTypeNotSupportedFor
        {
        };

        /// A vector traits to be defined for each autodiff number.
        template <typename V, class Enable = void>
        struct VectorTraits
        {
            /// The value type of each entry in the vector.
            using ValueType = VectorTraitsNotDefinedFor<V>;

            /// The template alias to replace the value type of a vector type with another value type.
            using ReplaceValueType = VectorReplaceValueTypeNotSupportedFor<V>;
        };

        /// A template alias used to get the type of the values in a vector type.
        template <typename V>
        using VectorValueType = typename VectorTraits<PlainType<V>>::ValueType;

        /// A template alias used to get the type of a vector that is equivalent to another but with a different value type.
        template <typename V, typename NewValueType>
        using VectorReplaceValueType = typename VectorTraits<PlainType<V>>::template ReplaceValueType<NewValueType>;

        /// A compile-time constant that indicates with a type is a vector type.
        template <typename V>
        constexpr bool isVector = !std::is_same_v<VectorValueType<PlainType<V>>, VectorTraitsNotDefinedFor<PlainType<V>>>;

        /// Implementation of VectorTraits for std::vector.
        template <typename T, template <class> typename Allocator>
        struct VectorTraits<std::vector<T, Allocator<T>>>
        {
            using ValueType = T;

            template <typename NewValueType>
            using ReplaceValueType = std::vector<NewValueType, Allocator<NewValueType>>;
        };

    } // namespace detail
} // namespace autodiff
