#ifndef ASN_MESH_UTILS_H
#define ASN_MESH_UTILS_H

#include "../Config.h"

namespace Asn::Mesh
{
    enum MeshElemType
    {
        Unknown,
        Interior,
        Boundary,
    };

    String to_lower(const String &str);

    String get_file_extension(const String &s);

    template <typename Type>
    void switch_data(Type &data, const Int &data1, const Int &data2);
    template <Int N>
    void switch_data_array(Array<Int, N> &data, const Int &data1, const Int &data2);

    template <Int N>
    Bool contain(Array<Int, N> &list, const Int &data);

    Bool seek_keyword(std::ifstream &fp, const String &keyword);

    template <Int DIM, Real TOL = DEFAULT_REL_TOL>
    struct CmpPoint
    {
        Bool operator()(const Point<DIM> &x, const Point<DIM> &y) const;
    };
};

#endif
