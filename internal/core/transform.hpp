#ifndef VDS_SLICE_TRANSFORM_HPP
#define VDS_SLICE_TRANSFORM_HPP

#include <OpenVDS/IJKCoordinateTransformer.h>

#include "typedefs.hpp"

namespace Coords {

struct Transformer : private OpenVDS::IJKCoordinateTransformer {

    Transformer(OpenVDS::IJKCoordinateTransformer transformer)
        : OpenVDS::IJKCoordinateTransformer(transformer)
    {}

    /* Conversion to annotations */
    TracePoint to_annotation(WorldPoint const& world) const noexcept (false) {
        auto annotation = this->WorldToAnnotation({world.x(), world.y(), 0});
        return {annotation[0], annotation[1]};
    }

    TracePoint to_annotation(IJPoint const& ij) const noexcept (false) {
        auto annotation = this->IJKPositionToAnnotation({ij.i(), ij.j(), 0});
        return {annotation[0], annotation[1]};
    }


    TracePoint to_annotation(TracePoint const& point) const noexcept (false) {
        return point;
    }

    /* Conversion to index */

    IndexCenter< 2 > to_center(TracePoint const& trace) const noexcept (false) {
        auto voxel = this->AnnotationToIJKPosition({trace[0], trace[1], 0});
        return { voxel[0] + 0.5, voxel[1] + 0.5 };
    }

    IndexCenter< 3 > to_center(AnnotationPoint< 3 > const& annotation) const noexcept (false) {
        auto voxel = this->AnnotationToIJKPosition(
            {annotation[0], annotation[1], annotation[2]}
        );
        return { voxel[0] + 0.5, voxel[1] + 0.5, voxel[2] + 0.5 };
    }

    template< std::size_t Dims >
    IndexCenter< Dims > to_center(
        IndexPoint< Dims > const& index
    ) const noexcept (false) {
        return index += 0.5;
    }

    IndexCenter< 2 > to_center(
        AnnotationPoint< 2 > const& annotation
    ) const noexcept (false) {
        auto voxel = this->AnnotationToIJKPosition(
            {annotation[0], annotation[1], 0 }
        );
        return { voxel[0] + 0.5, voxel[1] + 0.5 };
    }

    bool is_out_of_range(TracePoint const& point) const noexcept (false) {
        auto const& start = this->IJKAnnotationStart();
        auto const& end   = this->IJKAnnotationEnd();

        TracePoint min(start[0], start[1]);
        TracePoint max(  end[0],   end[1]);

        return ((point < min) or (point > max));
    }

};

} // namespace Coords

#endif // VDS_SLICE_TRANSFORM_HPP
