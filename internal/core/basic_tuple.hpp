#ifndef VDS_SLICE_BASIC_TUPLE_HPP
#define VDS_SLICE_BASIC_TUPLE_HPP

#include <array>
#include <iostream>

namespace Coords {
namespace detail {

template < typename Base, typename T, std::size_t Dims >
class basic_tuple : private std::array< T, Dims > {
    using base_type = std::array< T, Dims >;
    static_assert(Dims > 0, "0-tuples are non-sensical, is this a bug?");

    /*
     * Dimensionalities and coordinates are all structurally identical, but
     * semantically different. It makes perfect sense for them all to be
     * different types, but it's quite tedious and difficult to maintain
     * multiple identical implementations.
     *
     * Which means it's time to bring out the mixins with CRTP.
     */

public:

    static constexpr const auto dimensions = Dims;

    basic_tuple() noexcept (true) { this->fill(0); }
    basic_tuple(const base_type& t) noexcept (true) : base_type(t) {}

    using reference  = typename base_type::reference;
    using value_type = typename base_type::value_type;

    template < typename... Vs >
    basic_tuple(T v, Vs... vs) noexcept (true) :
    basic_tuple(base_type{ v, static_cast< T >(vs) ... }) {
        /*
         * A super-duper hack to support brace-initialization, emplace back and
         * similar.
         *
         * Really, this boils down to something along the lines of:
         *
         *  basic_tuple(std::size_t x, std::size_t y, std::size_t z) :
         *      std::array<>({ x, y, z })
         *  {}
         *
         * but for N dimensions. The first argument is fixed to size_t (and not
         * variadic), otherwise arbitrary overloads are easily picked up on,
         * and a compile error quickly follows.
         *
         * It delegates to the array constructor, in case more behaviour should
         * be added to it.
         */
        static_assert(
            sizeof...(vs) + 1 == Dims,
            "constructor must have exactly Dims arguments"
        );
    }

    /* inherit methods from std::array */
    using base_type::operator [];
    using base_type::begin;
    using base_type::end;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::size;
    using base_type::front;
    using base_type::back;
    using base_type::at;

    std::string string() const;

    /*
     * Comparisons, but only within same type - no conversion!
     */
    friend
    bool operator != (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs != rhs;
    }

    friend
    bool operator == (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs == rhs;
    }

    friend
    bool operator < (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs < rhs;
    }

    friend
    bool operator <= (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs <= rhs;
    }

    friend
    bool operator > (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs > rhs;
    }

    friend
    bool operator >= (const Base& left, const Base& right) noexcept (true) {
        const auto& lhs = static_cast< const base_type& >(left);
        const auto& rhs = static_cast< const base_type& >(right);
        return lhs >= rhs;
    }

    friend std::ostream&
    operator << (std::ostream& o, const Base& p) {
        static_assert(
            Dims > 1,
            "ostream << is only implemented for Dims > 1, "
            "fix it by writing a better join");

        /*
         * C++ :------------)
         *
         * '(' + ', '.join(*this) + ')'
         */
        o << '(';
        for (auto x = p.begin(); x != p.end() - 1; ++x)
            o << *x << ", ";
        return o << *(p.end() - 1) << ')';
    }
};

template < typename Base, typename T, std::size_t Dims >
std::string basic_tuple< Base, T, Dims >::string() const {
    const auto& self = static_cast< const base_type& >(*this);
    std::string fmt = "(";
    for (auto x = self.begin(); x != self.end() - 1; ++x) {
        fmt += (std::to_string(*x) + ", ");
    }
    fmt += (std::to_string(*(self.end() - 1)) + ")");
    return fmt;
}

} // namespace detail

} // namespace Coords

#endif // VDS_SLICE_BASIC_TUPLE_HPP
