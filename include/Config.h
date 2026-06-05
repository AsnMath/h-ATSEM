#ifndef ASN_CONFIG_H
#define ASN_CONFIG_H

#include <complex>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <array>
#include <string>
#include <vector>
#include <cassert>
#include <iomanip>

#include "../3rdparty/Eigen/Eigen"
#include "../3rdparty/unsupported/Eigen/CXX11/Tensor"
#include "../3rdparty/Autodiff/reverse/var.hpp"
#include "../3rdparty/Autodiff/reverse/eigen.hpp"

static_assert(__cplusplus >= 202000, "C++ 20 or higher is required.");

#define ASN_ALERT printf("\n\033[38;2;0;255;0mAsn: Alert at %s: %d.\033[0m\n", __FILE__, __LINE__);
#define ASN_ERROR(text) printf("\n\033[38;2;255;0;0mAsn: Error at %s: %d. %s.\033[0m\n", __FILE__, __LINE__, text); exit(0);

namespace Asn
{
    using Bool = bool;
    using Int = long long int;
    using Real = double;
    using Var = autodiff::Variable<Real>;
    using String = std::string;
    using Order = std::strong_ordering;

    template <typename FirstType, typename SecondType>
    using Pair = std::pair<FirstType, SecondType>;
    template <typename ...Args>
    using Tuple = std::tuple<Args ...>;

    template <typename Type, Int SIZE>
    using Array = std::array<Type, SIZE>;
    template <typename Type>
    using List = std::vector<Type>;

    template <typename Type, typename Cmp = std::less<Type>>
    using Set = std::set<Type, Cmp>;
    template <typename Type, typename Hash = std::hash<Type>, typename Equal = std::equal_to<Type>>
    using HashSet = std::unordered_set<Type, Hash, Equal>;
    template <typename Key, typename Value, typename Cmp = std::less<Key>>
    using Map = std::map<Key, Value, Cmp>;
    template <typename Key, typename Value, typename Cmp = std::less<Key>>
    using MultiMap = std::multimap<Key, Value, Cmp>;
    template <typename Key, typename Value, typename Hash = std::hash<Key>, typename Equal = std::equal_to<Key>>
    using HashMap = std::unordered_map<Key, Value, Hash, Equal>;

    template <typename Type>
    using Func = std::function<Type>;

    template <typename Type, Int ROWS = Eigen::Dynamic>
    using Vector = Eigen::Vector<Type, ROWS>;
    template <Int ROWS, typename Type = Real>
    using Point = Vector<Type, ROWS>;
    template <typename Type, Int ROWS = Eigen::Dynamic, Int COLS = Eigen::Dynamic>
    using Matrix = Eigen::Matrix<Type, ROWS, COLS>;
    template <typename Type, Int Mode = Eigen::RowMajor>
    using SparseMatrix = Eigen::SparseMatrix<Type, Mode>;
    template <typename Type, Int DIM>
    using Tensor = Eigen::Tensor<Type, DIM>;

    inline constexpr Int NPOS = -1;
    inline constexpr Real DEFAULT_REL_TOL = std::numeric_limits<float>::epsilon();
    inline constexpr Real DEFAULT_ABS_TOL = std::numeric_limits<float>::min();
    inline constexpr Real DEFAULT_MAX_TOL = std::numeric_limits<float>::max();
    inline constexpr Int DEFAULT_MAX_ITER = std::numeric_limits<short>::max();
    inline constexpr Int MAXIMUM_MAX_ITER = std::numeric_limits<int>::max() - 1;

    static_assert(DEFAULT_REL_TOL >= 0.0);
    static_assert(DEFAULT_ABS_TOL >= 0.0);
    static_assert(DEFAULT_MAX_TOL >= 0.0);
    static_assert(DEFAULT_MAX_ITER >= 0);
    static_assert(MAXIMUM_MAX_ITER >= 0);

    template <typename InputType, typename OutputType>
    inline OutputType convert(const InputType &input)
    {
        OutputType res(input.begin(), input.end());
        return res;
    }

    inline String to_fixed_precision(const Real &value, const Int &precision)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << value;
        return oss.str();
    }

    inline Real timer(const String &tag = "")
    {
        static HashMap<String, std::chrono::steady_clock::time_point> time_table;
        const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        Real duration = -1.0;
        if (time_table.contains(tag))
        {
            duration = std::chrono::duration_cast<std::chrono::duration<Real>>(now - time_table[tag]).count();
        }
        time_table[tag] = now;
        return duration;
    }
};

#endif
