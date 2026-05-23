#ifndef ASN_MATH_QUADRATURE_HPP
#define ASN_MATH_QUADRATURE_HPP

#include <fstream>
#include <iostream>

#include "Quadrature.h"

namespace Asn::Math
{
    template <>
    inline constexpr Real Quadrature::simplex_volume<1>()
    {
        return 1.0;
    }

    template <>
    inline constexpr Real Quadrature::simplex_volume<2>()
    {
        return 1.0 / 2.0;
    }

    template <>
    inline constexpr Real Quadrature::simplex_volume<3>()
    {
        return 1.0 / 6.0;
    }

    template <Int DIM>
    inline constexpr Real Quadrature::simplex_volume()
    {
        static_assert(DIM > 3);
        Real res = 1.0;
        for (Int i = 2; i <= DIM; i++)
        {
            res = res / static_cast<Real>(i);
        }
        return res;
    }

    inline Quadrature::Quadrature() : _order_line(0), _order_tri(0), _order_tet(0),
                                      order_line(_order_line), order_tri(_order_tri), order_tet(_order_tet),
                                      line(_line), tri(_tri), tet(_tet) {}

    inline Quadrature::Quadrature(const Int &order) : Quadrature()
    {
        init(order);
    }

    inline Quadrature::Quadrature(const Quadrature &quadrature) : Quadrature()
    {
        this->_order_line = quadrature.order_line;
        this->_order_tri = quadrature.order_tri;
        this->_order_tet = quadrature.order_tet;
        this->_line = quadrature.line;
        this->_tri = quadrature.tri;
        this->_tet = quadrature.tet;
    }

    inline Quadrature &Quadrature::init(const Int &order)
    {
        init_line(order);
        init_tri(order);
        init_tet(order);
        return *this;
    }

    inline Quadrature &Quadrature::init_line(const Int &order)
    {
        assert(1 <= order && order <= MAX_ORDER);
        this->_order_line = order;
        _line.clear();
        const String filename = QUAD_PATH + "line/o" + std::to_string(order / 2 + 1) + ".txt";
        std::ifstream fp(filename, std::ios::in);
        if (!fp.is_open())
        {
            ASN_ERROR("Fail to load quadrature data");
        }
        while (!fp.eof())
        {
            Point<1> p;
            Real w = 0.0;
            fp >> p(0);
            fp >> w;
            if (fp.good())
            {
                _line.push_back(std::make_pair(p, w));
            }
        }
        fp.close();
        return *this;
    }

    inline Quadrature &Quadrature::init_tri(const Int &order)
    {
        assert(1 <= order && order <= MAX_ORDER);
        this->_order_tri = order;
        _tri.clear();
        const String filename = QUAD_PATH + "tri/o" + std::to_string(order) + ".txt";
        std::ifstream fp(filename, std::ios::in);
        if (!fp.is_open())
        {
            ASN_ERROR("Fail to load quadrature data");
        }
        while (!fp.eof())
        {
            Point<2> p;
            Real w = 0.0;
            fp >> p(0);
            fp >> p(1);
            fp >> w;
            p(0) = (p(0) + 1.0) / 2.0;
            p(1) = (p(1) + 1.0) / 2.0;
            w = w / 2.0;
            if (fp.good())
            {
                _tri.push_back(std::make_pair(p, w));
            }
        }
        fp.close();
        return *this;
    }

    inline Quadrature &Quadrature::init_tet(const Int &order)
    {
        assert(1 <= order && order <= MAX_ORDER);
        this->_order_tet = order;
        _tet.clear();
        const String filename = QUAD_PATH + "tet/o" + std::to_string(order) + ".txt";
        std::ifstream fp(filename, std::ios::in);
        if (!fp.is_open())
        {
            ASN_ERROR("Fail to load quadrature data");
        }
        while (!fp.eof())
        {
            Point<3> p;
            Real w = 0.0;
            fp >> p(0);
            fp >> p(1);
            fp >> p(2);
            fp >> w;
            p(0) = (p(0) + 1.0) / 2.0;
            p(1) = (p(1) + 1.0) / 2.0;
            p(2) = (p(2) + 1.0) / 2.0;
            w = 3.0 * w / 4.0;
            if (fp.good())
            {
                _tet.push_back(std::make_pair(p, w));
            }
        }
        fp.close();
        return *this;
    }

    template <typename FuncType>
    inline Real Quadrature::int_line(const FuncType &f)
    {
        Real res = 0.0;
        for (Int i = 0; i < static_cast<Int>(this->line.size()); i++)
        {
            const auto &[p, w] = line[i];
            res += w * f(p(0));
        }
        return res * simplex_volume<1>();
    }

    template <typename FuncType>
    inline Real Quadrature::int_tri(const FuncType &f)
    {
        Real res = 0.0;
        for (Int i = 0; i < static_cast<Int>(this->tri.size()); i++)
        {
            const auto &[p, w] = tri[i];
            res += w * f(p(0), p(1));
        }
        return res * simplex_volume<2>();
    }

    template <typename FuncType>
    inline Real Quadrature::int_tet(const FuncType &f)
    {
        Real res = 0.0;
        for (Int i = 0; i < static_cast<Int>(this->tet.size()); i++)
        {
            const auto &[p, w] = tet[i];
            res += w * f(p(0), p(1), p(2));
        }
        return res * simplex_volume<3>();
    }

    inline Quadrature &Quadrature::operator=(const Quadrature &quadrature)
    {
        if (this != &quadrature)
        {
            this->_order_line = quadrature.order_line;
            this->_order_tri = quadrature.order_tri;
            this->_order_tet = quadrature.order_tet;
            this->_line = quadrature.line;
            this->_tri = quadrature.tri;
            this->_tet = quadrature.tet;
        }
        return *this;
    }
};

#endif
