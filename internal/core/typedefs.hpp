#ifndef VDS_SLICE_TYPEDEFS_HPP
#define VDS_SLICE_TYPEDEFS_HPP

#include <cstddef>

#include "basic_tuple.hpp"

namespace Coords {

using namespace detail;

/* World (cdp) coordinate system - defined at the center of the voxel */
using World = double;

struct WorldPoint : public basic_tuple< WorldPoint, World, 2 > {
    using base_type = basic_tuple< WorldPoint, World, 2 >;
    using base_type::base_type;

    World x() const noexcept (false) { return this->at(0); }
    World y() const noexcept (false) { return this->at(1); }
};

/* IJK (index) coordinate system - defined at the center of the voxel */
using Index = double;

template < std::size_t Dims >
struct IndexPoint : public basic_tuple< IndexPoint< Dims >, Index, Dims > {
    using base_type = basic_tuple< IndexPoint, Index, Dims >;
    using base_type::base_type;

    Index i() const noexcept (false) {
        static_assert(Dims >= 1, "i() only defined for Dims >= 1");
        return this->at(0);
    }
    Index j() const noexcept (false) {
        static_assert(Dims >= 2, "j() only defined for Dims >= 2");
        return this->at(1);
    }
    Index k() const noexcept (false) {
        static_assert(Dims >= 3, "k() only defined for Dims >= 3");
       return this->at(2);
    }
};

struct IJPoint : public IndexPoint< 2 > {
    using base_type = IndexPoint< 2 >;
    using base_type::base_type;
};

/* Annotation coordinate system - defined at the upper left corner of the voxel */
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

} // namespace Coords

#endif // VDS_SLICE_TYPEDEFS_HPP
