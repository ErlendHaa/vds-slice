#ifndef VDS_SLICE_TYPEDEFS_HPP
#define VDS_SLICE_TYPEDEFS_HPP

#include <cstddef>

#include "basic_tuple.hpp"

namespace Coords {

using namespace detail;

/* World (cdp) coordinate system - defined at the upper-left corner of the voxel */
using World = double;

struct WorldPoint : public basic_tuple< WorldPoint, World, 2 > {
    using base_type = basic_tuple< WorldPoint, World, 2 >;
    using base_type::base_type;

    World x() const noexcept (false) { return this->at(0); }
    World y() const noexcept (false) { return this->at(1); }
};

/* Annotation coordinate system - defined at the upper-left corner of the voxel */
using Annotation = double;

template < std::size_t Dims >
struct AnnotationPoint : public basic_tuple< AnnotationPoint< Dims >, Annotation, Dims > {
    using base_type = basic_tuple< AnnotationPoint, Annotation, Dims >;
    using base_type::base_type;
};

// Inline /crossline
struct TracePoint : public basic_tuple< TracePoint, Annotation, 2 > {
    using base_type = basic_tuple< TracePoint, Annotation, 2 >;
    using base_type::base_type;

    Annotation iline() const noexcept (false) { return this->at(0); }
    Annotation xline() const noexcept (false) { return this->at(1); }

    AnnotationPoint< 3 > extend(Annotation sample) const noexcept (true) {
        return { this->iline(), this->xline(), sample };
    }
};

namespace detail {

template < typename T, std::size_t Dims >
struct IndexBase : public basic_tuple< IndexBase< T, Dims >, T, Dims > {
    using base_type = basic_tuple< IndexBase< T, Dims >, T, Dims >;
    using base_type::base_type;

    T i() const noexcept (false) {
        static_assert(Dims >= 1, "i() only defined for Dims >= 1");
        return this->at(0);
    }
    T j() const noexcept (false) {
        static_assert(Dims >= 2, "j() only defined for Dims >= 2");
        return this->at(1);
    }
    T k() const noexcept (false) {
       static_assert(Dims >= 3, "k() only defined for Dims >= 3");
       return this->at(2);
    }
};

} // namespace

/* IJK (index) coordinate system - defined in the upper-left corner of the voxel*/
using Index = std::size_t;

template < std::size_t Dims >
struct IndexPoint : public IndexBase< Index, Dims > {
    using base_type = IndexBase< Index, Dims>;
    using base_type::base_type;
};

struct IJPoint : public IndexBase< Index, 2 > {
    using base_type = IndexBase< Index, 2 >;
    using base_type::base_type;
};

using Center = double;

template < std::size_t Dims >
struct IndexCenter : public IndexBase< Center, Dims > {
    using base_type = IndexBase< Center, Dims >;
    using base_type::base_type;
};

} // namespace Coords

#endif // VDS_SLICE_TYPEDEFS_HPP
