#ifndef ASN_INDEX_H
#define ASN_INDEX_H

#include "../Config.h"

namespace Asn::Math
{
    template <Int DIM>
    class Index;
    template <Int DIM, typename DataType>
    class IndexedData;

    template <Int DIM>
    class Index
    {
        static_assert(DIM > 0);

    private:
        std::array<Int, DIM> data;

    public:
        Index();
        Index(const std::initializer_list<Int> &data);
        explicit Index(const std::array<Int, DIM> &data);
        template <typename ...IndexTypes>
        explicit Index(const Int &first_index, const IndexTypes & ...other_indices);
        Index(const Index &indexed_data) = default;
        ~Index() = default;

        Int &operator()(const Int &n);
        const Int &operator()(const Int &n) const;

        Index &operator=(const Index &index) = default;

        Int sum() const;
    };

    template <Int DIM, typename DataType>
    class IndexedData
    {
        static_assert(DIM > 0);

    private:
        std::conditional_t<DIM == 1, std::vector<DataType>, std::vector<IndexedData<DIM - 1, DataType>>> data;

    public:
        IndexedData() = default;
        IndexedData(const Int &order, const DataType &default_value);
        IndexedData(const IndexedData &indexed_data) = default;
        ~IndexedData() = default;

        void init(const Int &order, const DataType &default_value);

        DataType &operator()(const Index<DIM> &index);
        const DataType &operator()(const Index<DIM> &index) const;

        IndexedData &operator=(const IndexedData &indexed_data) = default;
    };
};

#include "Index.hpp"

#endif
