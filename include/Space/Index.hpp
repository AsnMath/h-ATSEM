#ifndef ASN_INDEX_HPP
#define ASN_INDEX_HPP

#include "Index.h"

namespace Asn::Math
{
    template <Int DIM>
    inline std::strong_ordering operator<=>(const Index<DIM> &index_1, const Index<DIM> &index_2)
    {
        for (Int i = 0; i < DIM; i++)
        {
            if (index_1(i) < index_2(i))
            {
                return std::strong_ordering::less;
            }
            if (index_1(i) > index_2(i))
            {
                return std::strong_ordering::greater;
            }
        }
        return std::strong_ordering::equal;
    }

    template <Int DIM>
    inline Index<DIM>::Index()
    {
        data.fill(NPOS);
    };

    template <Int DIM>
    inline Index<DIM>::Index(const std::initializer_list<Int> &data)
    {
        assert(data.size() == DIM);
        Int k = 0;
        for (auto it = data.begin(); it != data.end() && k < DIM; it++, k++)
        {
            this->data[k] = *it;
        }
    }

    template <Int DIM>
    inline Index<DIM>::Index(const std::array<Int, DIM> &data)
    {
        this->data = data;
    }

    template <Int DIM>
    template <typename ...IndexTypes>
    inline Index<DIM>::Index(const Int &first_index, const IndexTypes & ...other_indices) : Index({first_index, other_indices ...})
    {
        static_assert((sizeof...(other_indices) + 1) == DIM);
    }

    template <Int DIM>
    inline Int &Index<DIM>::operator()(const Int &n)
    {
        return this->data[n];
    }

    template <Int DIM>
    inline const Int &Index<DIM>::operator()(const Int &n) const
    {
        return this->data[n];
    }

    template <Int DIM>
    inline Int Index<DIM>::sum() const
    {
        Int s = 0;
        for (Int i = 0; i < data.size(); i++)
        {
            s += data[i];
        }
        return s;
    }

    template <Int DIM, typename DataType>
    inline IndexedData<DIM, DataType>::IndexedData(const Int &order, const DataType &default_value)
    {
        if (order >= 0)
        {
            this->init(order, default_value);
        }
    }

    template <Int DIM, typename DataType>
    inline void IndexedData<DIM, DataType>::init(const Int &order, const DataType &default_value)
    {
        this->data.resize(order + 1);
        if constexpr (DIM == 1)
        {
            std::fill(this->data.begin(), this->data.end(), default_value);
        }
        if constexpr (DIM > 1)
        {
            for (Int i = 0; i <= order; i++)
            {
                this->data[i].init(order - i, default_value);
            }
        }
    }

    template <Int DIM, typename DataType>
    inline DataType &IndexedData<DIM, DataType>::operator()(const Index<DIM> &index)
    {
        if constexpr (DIM == 1)
        {
            return this->data[index(0)];
        }
        if constexpr (DIM > 1)
        {
            Index<DIM - 1> new_index;
            for (Int i = 1; i < DIM; i++)
            {
                new_index(i - 1) = index(i);
            }
            return data[index(0)](new_index);
        }
        assert(false);
    }

    template <Int DIM, typename DataType>
    inline const DataType &IndexedData<DIM, DataType>::operator()(const Index<DIM> &index) const
    {
        if constexpr (DIM == 1)
        {
            return this->data[index(0)];
        }
        if constexpr (DIM > 1)
        {
            Index<DIM - 1> new_index;
            for (Int i = 1; i < DIM; i++)
            {
                new_index(i - 1) = index(i);
            }
            return data[index(0)](new_index);
        }
        assert(false);
    }
};

#endif
