#include "vds.h"

#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <utility>
#include <cmath>

#include "nlohmann/json.hpp"

#include "vdshandle.h"

using namespace std;

namespace internal {

requestdata requestdata_from_dump( const nlohmann::json::string_t& dump ) {
    requestdata tmp{ new char[dump.size()], nullptr, dump.size() };
    std::copy(dump.begin(), dump.end(), tmp.data);
    return tmp;
}

struct requestdata fetch_slice(
    const std::string url,
    const std::string credentials,
    const Axis ax,
    const int lineno
) {
    VDSHandle vds_handle(url, credentials);
    return vds_handle.get_slice_of( ax, lineno );
}

struct requestdata fetch_slice_metadata(
    const std::string url,
    const std::string credentials,
    const Axis ax
) {
    VDSHandle vds_handle(url, credentials);

    nlohmann::json meta;
    meta["format"] = vds_handle.get_format_string_of_seismic_channel();

    meta["x"] = vds_handle.get_axis_metadata( ax, 0 );
    meta["y"] = vds_handle.get_axis_metadata( ax, 1 );

    return requestdata_from_dump( meta.dump() );
}

struct requestdata fetch_fence(
    const std::string& url,
    const std::string& credentials,
    const enum CoordinateSystem coordinate_system,
    const float* coordinates,
    const size_t npoints,
    const enum InterpolationMethod interpolation_method
) {
    VDSHandle vds_handle(url, credentials);
    return vds_handle.get_fence_of( coordinate_system, coordinates, npoints, interpolation_method);
}

struct requestdata handle_error(
    const std::exception& e
) {
    requestdata buf {};
    buf.err = new char[std::strlen(e.what()) + 1];
    std::strcpy(buf.err, e.what());
    return buf;
}

struct requestdata metadata(
    const std::string& url,
    const std::string& credentials
) {
    VDSHandle vds_handle(url, credentials);

    nlohmann::json meta;
    meta["format"] = vds_handle.get_format_string_of_seismic_channel();
    meta["crs"] = vds_handle.get_crs_string();

    const auto bbox = vds_handle.get_bounding_box();
    meta["boundingBox"]["ij"]   = bbox.index;
    meta["boundingBox"]["cdp"]  = bbox.world;
    meta["boundingBox"]["ilxl"] = bbox.annotation;

    for (int i = 2; i >= 0 ; i--) {
        meta["axis"].push_back( vds_handle.get_axis_metadata( i ) );
    }
    return internal::requestdata_from_dump( meta.dump() );
}

struct requestdata fetch_fence_metadata(
    const std::string url,
    const std::string credentials,
    const size_t npoints
) {
    VDSHandle vds_handle(url, credentials);

    nlohmann::json meta;
    meta["shape"] = vds_handle.get_shape_metadata(npoints);
    meta["format"] = vds_handle.get_format_string_of_seismic_channel();

    return internal::requestdata_from_dump( meta.dump() );
}

} // namespace internal

/*
 * Below here the interface methods are defined.
 */


struct requestdata slice(
    char const * const vds,
    char const * const credentials,
    const int lineno,
    const Axis ax
) {
    const std::string cube(vds);
    const std::string cred(credentials);

    try {
        return internal::fetch_slice(cube, cred, ax, lineno);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata slice_metadata(
    char const * const vds,
    char const * const credentials,
    const Axis ax
) {
    const std::string cube(vds);
    const std::string cred(credentials);

    try {
        return internal::fetch_slice_metadata(cube, cred, ax);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata fence(
    char const * const vds,
    char const * const credentials,
    const enum CoordinateSystem coordinate_system,
    float const * const coordinates,
    const size_t npoints,
    const enum InterpolationMethod interpolation_method
) {
    const std::string cube(vds);
    const std::string cred(credentials);

    try {
        return internal::fetch_fence(
            cube, cred, coordinate_system, coordinates, npoints,
            interpolation_method);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata fence_metadata(
    char const * const vds,
    char const * const credentials,
    const size_t npoints
) {
    const std::string cube(vds);
    const std::string cred(credentials);

    try {
        return internal::fetch_fence_metadata(cube, cred, npoints);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata metadata(
    char const * const vds,
    char const * const credentials
) {
    try {
        const std::string cube(vds);
        const std::string cred(credentials);
        return internal::metadata(cube, cred);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

void requestdata_delete(struct requestdata* buf) {
    if (!buf)
        return;

    delete[] buf->data;
    delete[] buf->err;
    *buf = requestdata {};
}
