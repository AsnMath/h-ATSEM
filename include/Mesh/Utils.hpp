#ifndef ASN_MESH_UTILS_HPP
#define ASN_MESH_UTILS_HPP

#include <iomanip>
#include <chrono>
#include <fstream>

#include "Utils.h"

namespace Asn::Mesh
{
    inline String to_lower(const String &str)
    {
        String res = str;
        for (Int i = 0; i < static_cast<Int>(res.length()); i++)
        {
            res[i] = std::tolower(res[i]);
        }
        return res;
    }

    inline String get_file_extension(const String &s)
    {
        const size_t pos = s.find_last_of(".");
        if (pos >= s.size())
        {
            return "";
        }
        return s.substr(pos + 1);
    }

    template <typename Type>
    inline void switch_data(Type &data, const Int &data1, const Int &data2)
    {
        if (data1 == data2)
        {
            return;
        }
        auto it1 = data.find(data1);
        auto it2 = data.find(data2);
        const Bool contain_1 = (it1 != data.end());
        const Bool contain_2 = (it2 != data.end());
        if (contain_1 && contain_2)
        {
            return;
        }
        if (contain_1)
        {
            data.erase(it1);
            data.insert(data2);
        }
        if (contain_2)
        {
            data.erase(it2);
            data.insert(data1);
        }
        return;
    }

    template <Int N>
    inline void switch_data_array(Array<Int, N> &data, const Int &data1, const Int &data2)
    {
        if (data1 == data2)
        {
            return;
        }
        for (Int i = 0; i < N; i++)
        {
            if (data[i] == data1)
            {
                data[i] = data2;
            }
            else if (data[i] == data2)
            {
                data[i] = data1;
            }
        }
        return;
    }

    template <Int N>
    Bool contain(Array<Int, N> &list, const Int &data)
    {
        for (Int i = 0; i < N; i++)
        {
            if (list[i] == data)
            {
                return true;
            }
        }
        return false;
    }

    inline Bool seek_keyword(std::ifstream &fp, const String &keyword)
    {
        String s;
        do
        {
            fp >> s;
        } while (s != keyword && fp.good());
        return !fp.eof();
    }

    template <Int DIM, Real TOL>
    inline Bool CmpPoint<DIM, TOL>::operator()(const Point<DIM> &x, const Point<DIM> &y) const
    {
        for (Int i = 0; i < DIM; i++)
        {
            const Real tol = TOL * std::max(TOL, std::max(std::abs(x(i)), std::abs(y(i))));
            if (x(i) < (y(i) - tol))
            {
                return true;
            }
            if (x(i) > (y(i) + tol))
            {
                return false;
            }
        }
        return false;
    }
};

#endif
