#ifndef ASN_MATH_QUADRATURE_H
#define ASN_MATH_QUADRATURE_H

#include <filesystem>
#include <source_location>

#include "../Config.h"

namespace Asn::Math
{
    inline const String QUAD_PATH = std::filesystem::path(std::source_location::current().file_name()).parent_path().append("./QuadInfo/").string();

    class Quadrature;

    class Quadrature
    {
    public:
        template <Int DIM>
        static constexpr Real simplex_volume();

    public:
        static constexpr Int MAX_ORDER = 40;

    private:
        Int _order_line;
        Int _order_tri;
        Int _order_tet;
        List<Pair<Point<1>, Real>> _line;
        List<Pair<Point<2>, Real>> _tri;
        List<Pair<Point<3>, Real>> _tet;

    public:
        const Int &order_line;
        const Int &order_tri;
        const Int &order_tet;
        const List<Pair<Point<1>, Real>> &line;
        const List<Pair<Point<2>, Real>> &tri;
        const List<Pair<Point<3>, Real>> &tet;

    public:
        Quadrature();
        explicit Quadrature(const Int &order);
        Quadrature(const Quadrature &quadrature);
        ~Quadrature() = default;

        Quadrature &init(const Int &order);
        Quadrature &init_line(const Int &order);
        Quadrature &init_tri(const Int &order);
        Quadrature &init_tet(const Int &order);

        template <typename FuncType>
        Real int_line(const FuncType &f);
        template <typename FuncType>
        Real int_tri(const FuncType &f);
        template <typename FuncType>
        Real int_tet(const FuncType &f);

        Quadrature &operator=(const Quadrature &quadrature);
    };
};

#endif
