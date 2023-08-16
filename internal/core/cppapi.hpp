#ifndef VDS_SLICE_CPPAPI_HPP
#define VDS_SLICE_CPPAPI_HPP

#include <vector>

#include "ctypes.h"

#include "attribute.hpp"
#include "datahandle.hpp"
#include "direction.hpp"
#include "regularsurface.hpp"

namespace cppapi {

void slice(
    DataHandle& handle,
    Direction const direction,
    int lineno,
    std::vector< Bound > const& bounds,
    response* out
) noexcept (false);

void fence(
    DataHandle& handle,
    enum coordinate_system coordinate_system,
    const float* coordinates,
    size_t npoints,
    enum interpolation_method interpolation_method,
    response* out
) noexcept (false);

void horizon_size(
    DataHandle& handle,
    RegularSurface const& surface,
    float above,
    float below,
    std::size_t* out
) noexcept (false);

void horizon(
    DataHandle& handle,
    RegularSurface const& surface,
    float above,
    float below,
    enum interpolation_method interpolation,
    std::size_t from,
    std::size_t to,
    void* out
) noexcept (false);

void attributes(
    DataHandle& handle,
    Horizon const& horizon,
    RegularSurface const& surface,
    VerticalWindow const& src_window,
    VerticalWindow const& dst_window,
    enum attribute* attributes,
    std::size_t nattributes,
    std::size_t from,
    std::size_t to,
    void** out
) noexcept (false);

/**
 * Given two input surfaces, primary and secondary, updates third surface,
 * aligned, which is expected to be shaped as primary surface, with data
 * belonging to secondary surface.
 *
 * For each point on primary surface a nearest point on the secondary surface
 * will be found and its value will be written to the resulting aligned surface.
 *
 * If according to the algorithm described above surfaces appear to intersect,
 * exception would be thrown.
 *
 * If the resulting point appears to be out of secondary surface bounds, aligned
 * surface fillvalue will be stored at the position. If for the primary or
 * secondary surface at the point the value of the data is surface's
 * corresponding fillvalue, aligned surface fillvalue will be stored at the
 * postion.
 *
 * Additionally a parameter primary_is_top would be set determining whether
 * primary or resulting aligned surface appeared on top of another.
 */
void align_surfaces(
    RegularSurface const& primary,
    RegularSurface const& secondary,
    RegularSurface &aligned,
    bool* primary_is_top
) noexcept (false);

void slice_metadata(
    DataHandle& handle,
    Direction const direction,
    int lineno,
    std::vector< Bound > const& bounds,
    response* out
) noexcept (false);


void fence_metadata(
    DataHandle& handle,
    size_t npoints,
    response* out
) noexcept (false);

void metadata(
    DataHandle& handle,
    response* out
) noexcept (false);

void attributes_metadata(
    DataHandle& handle,
    std::size_t nrows,
    std::size_t ncols,
    response* out
) noexcept (false);

/* implementations */

namespace {

void to_response(
    std::unique_ptr< char[] > data,
    std::int64_t const size,
    response* response
) {
    /* The data should *not* be free'd on success, as it's returned to CGO */
    response->data = data.release();
    response->size = static_cast<unsigned long>(size);
}

bool equal(const char* lhs, const char* rhs) {
    return std::strcmp(lhs, rhs) == 0;
}

} // namespace

namespace detail {

template< typename Coordinate >
void fetch_fence(
    DataHandle& handle,
    const float* coordinates,
    size_t npoints,
    enum interpolation_method interpolation_method,
    response* out
) {
    MetadataHandle const& metadata = handle.get_metadata();

    std::unique_ptr< voxel[] > coords(new voxel[npoints]{{0}});

    auto transformer = metadata.transformer();

    Axis inline_axis = metadata.iline();
    Axis crossline_axis = metadata.xline();

    for (size_t i = 0; i < npoints; i++) {
        const float x = *(coordinates++);
        const float y = *(coordinates++);
        Coordinate coordinate(x, y);

        auto annotation = transformer.to_annotation(coordinate);

        if (transformer.is_out_of_range(annotation)) {
            throw std::runtime_error(
                "Coordinate " + coordinate.string() + " is out of range"
            );
        }

        auto const center = transformer.to_center(annotation);

        coords[i][   inline_axis.dimension()] = center.i();
        coords[i][crossline_axis.dimension()] = center.j();
    }

    std::int64_t const size = handle.traces_buffer_size(npoints);

    std::unique_ptr< char[] > data(new char[size]);

    handle.read_traces(
        data.get(),
        size,
        coords.get(),
        npoints,
        interpolation_method
    );

    return to_response(std::move(data), size, out);
}

} // namespace detail

} // namespace cppapi

#endif // VDS_SLICE_CPPAPI_HPP
