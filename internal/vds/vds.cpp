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

#include <OpenVDS/OpenVDS.h>
#include <OpenVDS/KnownMetadata.h>
#include <OpenVDS/IJKCoordinateTransformer.h>

using namespace std;

namespace internal {

struct VoxelBounds {
    int lower[OpenVDS::VolumeDataLayout::Dimensionality_Max]{0, 0, 0, 0, 0, 0};
    int upper[OpenVDS::VolumeDataLayout::Dimensionality_Max]{1, 1, 1, 1, 1, 1};
};

enum VDSChannelID {
    Amplitude = 0,
    Trace = 1,
    SegyTraceHeader = 2
};

enum VDSLevelOfDetailID {
    Level0 = 0
};

enum VDSAxisID {
    DepthSampleTime=0,
    Crossline=1,
    Inline=2,
};

struct AxisUnitCombination {
    const char* label;
    const std::vector<const char*> units;

    AxisUnitCombination( const char* l, const std::vector<const char*> vu) : label(l), units(std::move(vu))
    {}
};

class ValidZAxisCombinations {

    private:
        enum Index{
            depth=0,
            time=1,
            sample=2
        };

        static const std::array<AxisUnitCombination, 3 > label_unit_combinations_;
        static const std::array<const char*, 3> axis_labels_;

    public:
        static const std::array<const char*, 3>& axis_labels(){
            return axis_labels_;
        }

        static const std::vector<const char*>& depth_units() {
            return label_unit_combinations_[Index::depth].units;
        }

        static const std::vector<const char*>& time_units() {
            return label_unit_combinations_[Index::time].units;
        }

        static const std::vector<const char*>& sample_units() {
            return label_unit_combinations_[Index::sample].units;
        }
};

/* Define some convenient lookup tables for labels and units */
const std::array<AxisUnitCombination, 3 > ValidZAxisCombinations::label_unit_combinations_ {
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Depth(),
        std::vector<const char*>{
            OpenVDS::KnownUnitNames::Meter(),
            OpenVDS::KnownUnitNames::Foot(),
            OpenVDS::KnownUnitNames::USSurveyFoot()
        }
    ),
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Time(),
        std::vector<const char*>{
            OpenVDS::KnownUnitNames::Millisecond(),
            OpenVDS::KnownUnitNames::Second()
        }
    ),
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Sample(),
        std::vector<const char*>{
            OpenVDS::KnownUnitNames::Unitless()
        }
    )
};

const std::array<const char*, 3> ValidZAxisCombinations::axis_labels_ {
    label_unit_combinations_[Index::depth].label,
    label_unit_combinations_[Index::time].label,
    label_unit_combinations_[Index::sample].label,
};

requestdata requestdata_from_dump( const nlohmann::json::string_t& dump ) {
    requestdata tmp{ new char[dump.size()], nullptr, dump.size() };
    std::copy(dump.begin(), dump.end(), tmp.data);
    return tmp;
}

requestdata requestdata_from_requested_data( std::unique_ptr< char[] >  &data, std::size_t size ) {
    requestdata tmp{ data.get(), nullptr, size };
    data.release();
    return tmp;
}

template<typename REQUEST_TYPE>
requestdata finalize_request( std::shared_ptr<REQUEST_TYPE>& request,
                            const std::string message,
                            std::unique_ptr< char[] >& data,
                            const std::size_t size ) {

    const bool success = request.get()->WaitForCompletion();
    if( not success ) {
        throw std::runtime_error(message);
    }

    return internal::requestdata_from_requested_data( data, size );
}

class VDSHandle {

    private:
        OpenVDS::ScopedVDSHandle handle_;
        OpenVDS::Error error_;
        OpenVDS::VolumeDataAccessManager access_manager_;
        const OpenVDS::VolumeDataLayout *layout_;
        OpenVDS::IJKCoordinateTransformer ijk_coordinate_transformer_;

        const std::string seismic_channel_name_{"Amplitude"};

        void validate_dimension() {
            if (layout_->GetDimensionality() != 3) {
                throw std::runtime_error(
                    "Unsupported VDS, expected 3 dimensions, got " +
                    std::to_string(layout_->GetDimensionality())
                );
            }
        }

    public:

        VDSHandle( std::string url, std::string credentials) {
            this->handle_ = OpenVDS::Open(url, credentials, this->error_);
            if(this->error_.code != 0)
                throw std::runtime_error("Could not open VDS: " + this->error_.string);

            this->access_manager_ = OpenVDS::GetAccessManager(handle_);
            this->layout_ = access_manager_.GetVolumeDataLayout();

            if (this->layout_ == nullptr)
                throw std::runtime_error("VDS does not contain valid data layout");

            this->ijk_coordinate_transformer_ = OpenVDS::IJKCoordinateTransformer(this->layout_);

            validate_dimension();
        }

        OpenVDS::VolumeDataAccessManager& access_manager() {
            return this->access_manager_;
        }

        const OpenVDS::VolumeDataLayout &layout() const {
            return *this->layout_;
        }

        const OpenVDS::IJKCoordinateTransformer &ijk_coordinate_transformer() const {
            return this->ijk_coordinate_transformer_;
        }

        std::string get_channel_format_string( VDSChannelID id ) {
            using namespace OpenVDS;
            VolumeDataFormat format = layout_->GetChannelFormat(id);
            switch (format) {
                case OpenVDS::VolumeDataFormat::Format_U8:  return "<u1";
                case OpenVDS::VolumeDataFormat::Format_U16: return "<u2";
                case OpenVDS::VolumeDataFormat::Format_R32: return "<f4";
                default: {
                    throw std::runtime_error("unsupported VDS format type");
                }
            }
        }

        std::string get_crs_string() const {
            auto crs = OpenVDS::KnownMetadata::SurveyCoordinateSystemCRSWkt();
            return layout_->GetMetadataString(crs.GetCategory(), crs.GetName());
        }

        nlohmann::json get_axis_metadata( const int voxel_dim ) const {
            return nlohmann::json {
                { "annotation", layout_->GetDimensionName(voxel_dim)       },
                { "min",        layout_->GetDimensionMin(voxel_dim)        },
                { "max",        layout_->GetDimensionMax(voxel_dim)        },
                { "samples",    layout_->GetDimensionNumSamples(voxel_dim) },
                { "unit",       layout_->GetDimensionUnit(voxel_dim)       },
            };
        }

        nlohmann::json get_shape_metadata( const int npoints ) const {
            return nlohmann::json::array({ npoints,
                                           layout_->GetDimensionNumSamples(VDSAxisID::DepthSampleTime)});
        }


        int convert_ijk_to_voxel_axis_id( const int ijk_axis_id ) const {
            const OpenVDS::IntVector3& mapping = this->ijk_coordinate_transformer_.IJKToVoxelDimensionMap();
            if ( ijk_axis_id > -1 && ijk_axis_id < 3 ) {
                return mapping[ijk_axis_id];
            }
            else {
                throw std::runtime_error("Unhandled axis");
            }
        }

        int convert_voxel_to_ijk_axis_id( const int voxel_axis_id ) const {
            return convert_ijk_to_voxel_axis_id( voxel_axis_id );
        }

        int get_max_voxel_index( const int voxel_dimension ) const {
            const int world_axis_id = convert_voxel_to_ijk_axis_id(voxel_dimension);
            return layout_->GetDimensionNumSamples(world_axis_id);
        }

};


int axis_todim(Axis ax) {
    switch (ax) {
        case I:
        case INLINE:
            return 0;
        case J:
        case CROSSLINE:
            return 1;
        case K:
        case DEPTH:
        case TIME:
        case SAMPLE:
            return 2;
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
}

CoordinateSystem axis_tosystem(Axis ax) {
    switch (ax) {
        case I:
        case J:
        case K:
            return INDEX;
        case INLINE:
        case CROSSLINE:
        case DEPTH:
        case TIME:
        case SAMPLE:
            return ANNOTATION;
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
}

const std::string axis_tostring(Axis ax) {
    switch (ax) {
        case I:         return std::string( OpenVDS::KnownAxisNames::I()         );
        case J:         return std::string( OpenVDS::KnownAxisNames::J()         );
        case K:         return std::string( OpenVDS::KnownAxisNames::K()         );
        case INLINE:    return std::string( OpenVDS::KnownAxisNames::Inline()    );
        case CROSSLINE: return std::string( OpenVDS::KnownAxisNames::Crossline() );
        case DEPTH:     return std::string( OpenVDS::KnownAxisNames::Depth()     );
        case TIME:      return std::string( OpenVDS::KnownAxisNames::Time()      );
        case SAMPLE:    return std::string( OpenVDS::KnownAxisNames::Sample()    );
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
}

OpenVDS::InterpolationMethod to_interpolation(InterpolationMethod interpolation) {
    switch (interpolation)
    {
        case NEAREST: return OpenVDS::InterpolationMethod::Nearest;
        case LINEAR: return OpenVDS::InterpolationMethod::Linear;
        case CUBIC: return OpenVDS::InterpolationMethod::Cubic;
        case ANGULAR: return OpenVDS::InterpolationMethod::Angular;
        case TRIANGULAR: return OpenVDS::InterpolationMethod::Triangular;
        default: {
            throw std::runtime_error("Unhandled interpolation method");
        }
    }
}

/*
 * Unit validation of Z-slices
 *
 * Verify that the units of the VDS' Z axis matches the requested slice axis.
 * E.g. a Time slice is only valid if the units of the Z-axis in the VDS is
 * "Seconds" or "Milliseconds"
 */
bool unit_validation(Axis ax, const char* zunit) {

    auto isoneof = [zunit](const char* x) {
        return !std::strcmp(x, zunit);
    };

    switch (ax) {
        case I:
        case J:
        case K:
        case INLINE:
        case CROSSLINE:
            return true;
        case DEPTH:
            return std::any_of(ValidZAxisCombinations::depth_units().begin(), ValidZAxisCombinations::depth_units().end(), isoneof);
        case TIME:
            return std::any_of(ValidZAxisCombinations::time_units().begin(), ValidZAxisCombinations::time_units().end(), isoneof);
        case SAMPLE:
            return std::any_of(ValidZAxisCombinations::sample_units().begin(), ValidZAxisCombinations::sample_units().end(), isoneof);
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
};

/*
 * Until we know more about how VDS' are constructed w.r.t. to axis ordering
 * we're gonna assume that all VDS' are ordered like so:
 *
 *     voxel[0] -> depth/time/sample
 *     voxel[1] -> crossline
 *     voxel[2] -> inline
 *
 * This function will return 0 if that's not the case
 */
bool axis_order_validation(const OpenVDS::VolumeDataLayout &layout) {
    if (std::strcmp(layout.GetDimensionName(VDSAxisID::Inline), OpenVDS::KnownAxisNames::Inline())) {
        return false;
    }

    if (std::strcmp(layout.GetDimensionName(VDSAxisID::Crossline), OpenVDS::KnownAxisNames::Crossline())) {
        return false;
    }

    auto z = layout.GetDimensionName(VDSAxisID::DepthSampleTime);

    auto isoneof = [z](const char* x) {
        return !std::strcmp(x, z);
    };

    return std::any_of(ValidZAxisCombinations::axis_labels().begin(), ValidZAxisCombinations::axis_labels().end(), isoneof);
}


void axis_validation(Axis ax, const OpenVDS::VolumeDataLayout &layout) {
    if (not axis_order_validation(layout)) {
        std::string msg = "Unsupported axis ordering in VDS, expected ";
        msg += "Depth/Time/Sample, Crossline, Inline";
        throw std::runtime_error(msg);
    }

    auto zaxis = layout.GetAxisDescriptor(VDSAxisID::DepthSampleTime);
    const char* zunit = zaxis.GetUnit();
    if (not unit_validation(ax, zunit)) {
        std::string msg = "Unable to use " + axis_tostring(ax);
        msg += " on cube with depth units: " + std::string(zunit);
        throw std::runtime_error(msg);
    }
}

int lineno_annotation_to_voxel(
    int lineno,
    int vdim,
    const OpenVDS::VolumeDataLayout &layout
) {
    /* Assume that annotation coordinates are integers */
    int min      = layout.GetDimensionMin(vdim);
    int max      = layout.GetDimensionMax(vdim);
    int nsamples = layout.GetDimensionNumSamples(vdim);

    auto stride = (max - min) / (nsamples - 1);

    if (lineno < min || lineno > max || (lineno - min) % stride) {
        throw std::runtime_error(
            "Invalid lineno: " + std::to_string(lineno) +
            ", valid range: [" + std::to_string(min) +
            ":" + std::to_string(max) +
            ":" + std::to_string(stride) + "]"
        );
    }

    int voxelline = (lineno - min) / stride;
    return voxelline;
}

int lineno_index_to_voxel(
    int lineno,
    int vdim,
    const OpenVDS::VolumeDataLayout &layout
) {
    /* Line-numbers in IJK match Voxel - do bound checking and return*/
    int min = 0;
    int max = layout.GetDimensionNumSamples(vdim) - 1;

    if (lineno < min || lineno > max) {
        throw std::runtime_error(
            "Invalid lineno: " + std::to_string(lineno) +
            ", valid range: [" + std::to_string(min) +
            ":" + std::to_string(max) +
            ":1]"
        );
    }

    return lineno;
}

/*
 * Convert target dimension/axis + lineno to VDS voxel coordinates.
 */
VoxelBounds get_voxel_bounds(
    Axis ax,
    int lineno,
    const VDSHandle& vds_handle
) {
    VoxelBounds voxel_bounds;
    for (std::size_t i = 0; i < 3; ++i)
        voxel_bounds.upper[i] = vds_handle.layout().GetDimensionNumSamples(i);

    int voxelline;
    const int dimension = axis_todim(ax);
    auto vdim   = vds_handle.convert_ijk_to_voxel_axis_id(dimension);

    const int system = axis_tosystem(ax);
    switch (system) {
        case ANNOTATION: {
            if (not vds_handle.ijk_coordinate_transformer().AnnotationsDefined()) {
                throw std::runtime_error("VDS doesn't define annotations");
            }
            voxelline = lineno_annotation_to_voxel(lineno, vdim, vds_handle.layout());
            break;
        }
        case INDEX: {
            voxelline = lineno_index_to_voxel(lineno, vdim, vds_handle.layout());
            break;
        }
        case CDP:
        default: {
            throw std::runtime_error("Unhandled coordinate system");
        }
    }

    voxel_bounds.lower[vdim] = voxelline;
    voxel_bounds.upper[vdim] = voxelline + 1;

    return voxel_bounds;
}


class BoundingBox {
public:
    explicit BoundingBox(
        const OpenVDS::VolumeDataLayout &layout
    ) : layout(layout)
    {
        transformer = OpenVDS::IJKCoordinateTransformer(&layout);
    }

    std::vector< std::pair<int, int> >       index()      noexcept (true);
    std::vector< std::pair<int, int> >       annotation() noexcept (true);
    std::vector< std::pair<double, double> > world()      noexcept (true);
private:
    OpenVDS::IJKCoordinateTransformer transformer;
    const OpenVDS::VolumeDataLayout &layout;
};


std::vector< std::pair<int, int> > BoundingBox::index() noexcept (true) {
    auto ils = layout.GetDimensionNumSamples(VDSAxisID::Inline) - 1;
    auto xls = layout.GetDimensionNumSamples(VDSAxisID::Crossline) - 1;

    return { {0, 0}, {ils, 0}, {ils, xls}, {0, xls} };
}

std::vector< std::pair<double, double> > BoundingBox::world() noexcept (true) {
    std::vector< std::pair<double, double> > world_points;

    auto points = this->index();
    std::for_each(points.begin(), points.end(),
        [&](const std::pair<int, int>& point) {
            auto p = this->transformer.IJKIndexToWorld(
                { point.first, point.second, 0 }
            );
            world_points.emplace_back(p[0], p[1]);
        }
    );

    return world_points;
};

std::vector< std::pair<int, int> > BoundingBox::annotation() noexcept (true) {
    auto points = this->index();
    std::transform(points.begin(), points.end(), points.begin(),
        [this](std::pair<int, int>& point) {
            auto anno = this->transformer.IJKIndexToAnnotation({
                point.first,
                point.second,
                0
            });
            return std::pair<int, int>{anno[0], anno[1]};
        }
    );

    return points;
};

struct requestdata fetch_slice(
    std::string url,
    std::string credentials,
    Axis ax,
    int lineno
) {
    VDSHandle vds_handle(url, credentials);

    axis_validation(ax, vds_handle.layout());

    VoxelBounds voxel_bounds = get_voxel_bounds(ax, lineno, vds_handle);

    auto format = vds_handle.layout().GetChannelFormat(VDSChannelID::Amplitude);
    const int size = vds_handle.access_manager().GetVolumeSubsetBufferSize(
        voxel_bounds.lower,
        voxel_bounds.upper,
        format,
        VDSLevelOfDetailID::Level0,
        VDSChannelID::Amplitude);

    std::unique_ptr< char[] > data(new char[size]());
    auto request = vds_handle.access_manager().RequestVolumeSubset(
        data.get(),
        size,
        OpenVDS::Dimensions_012,
        VDSLevelOfDetailID::Level0,
        VDSChannelID::Amplitude,
        voxel_bounds.lower,
        voxel_bounds.upper,
        format
    );

    return finalize_request( request, "Failed to fetch slice from VDS", data, size );
}

struct requestdata fetch_slice_metadata(
    std::string url,
    std::string credentials,
    Axis ax
) {
    VDSHandle vds_handle(url, credentials);

    axis_validation(ax, vds_handle.layout());

    auto dimension = axis_todim(ax);
    auto vdim = vds_handle.convert_ijk_to_voxel_axis_id(dimension);

    nlohmann::json meta;
    meta["format"] = vds_handle.get_channel_format_string(VDSChannelID::Amplitude);

    /*
     * SEGYImport always writes annotation 'Sample' for axis K. We, on the
     * other hand, decided that we base the valid input direction on the units
     * of said axis. E.g. ms/s -> Time, etc. This leads to an inconsistency
     * between what we require as input for axis K and what we return as
     * metadata. In the ms/s case we require the input to be asked for in axis
     * 'Time', but the return metadata can potentially say 'Sample'.
     *
     * TODO: Either revert the 'clever' unit validation, or patch the
     * K-annotation here. IMO the later is too clever for it's own good and
     * would be quite suprising for people that use this API in conjunction
     * with the OpenVDS library.
     */
    std::vector< int > dims;
    for (int i = 2; i >= 0; --i) {
        if (i == vdim) continue;
        dims.push_back(i);
    }

    meta["x"] = vds_handle.get_axis_metadata( dims[0] );
    meta["y"] = vds_handle.get_axis_metadata( dims[1] );

    return requestdata_from_dump( meta.dump() );
}

struct requestdata fetch_fence(
    const std::string& url,
    const std::string& credentials,
    enum CoordinateSystem coordinate_system,
    const float* coordinates,
    size_t npoints,
    enum InterpolationMethod interpolation_method
) {
    VDSHandle vds_handle(url, credentials);

    const auto dimension_map =
            vds_handle.layout().GetVDSIJKGridDefinitionFromMetadata().dimensionMap;

    unique_ptr< float[][OpenVDS::Dimensionality_Max] > coords(
        new float[npoints][OpenVDS::Dimensionality_Max]{{0}}
    );

    auto transform_coordinate = [&] (const float x, const float y) {
        switch (coordinate_system) {
            case INDEX:
                return OpenVDS::Vector<double, 3> {x, y, 0};
            case ANNOTATION:
                return vds_handle.ijk_coordinate_transformer().AnnotationToIJKPosition({x, y, 0});
            case CDP:
                return vds_handle.ijk_coordinate_transformer().WorldToIJKPosition({x, y, 0});
            default: {
                throw std::runtime_error("Unhandled coordinate system");
            }
        }
    };

    for (size_t i = 0; i < npoints; i++) {
        const float x = *(coordinates++);
        const float y = *(coordinates++);

        auto coordinate = transform_coordinate(x, y);

        auto validate_boundary = [&] (const int voxel) {
            const auto min = -0.5;
            const auto max = vds_handle.get_max_voxel_index( voxel ) - 0.5;
            if(coordinate[voxel] < min || coordinate[voxel] >= max) {
                const std::string coordinate_str =
                    "(" +std::to_string(x) + "," + std::to_string(y) + ")";
                throw std::runtime_error(
                    "Coordinate " + coordinate_str + " is out of boundaries "+
                    "in dimension "+ std::to_string(voxel)+ "."
                );
            }
        };

        for (size_t dim = 0; dim < 2; ++dim) {
            validate_boundary(dim);

            /* openvds uses rounding down for Nearest interpolation.
             * As it is counterintuitive, we fix it by snapping to nearest index
             * and rounding half-up.
             */
            if (interpolation_method == NEAREST) {
                coordinate[dim] = std::round(coordinate[dim] + 1) - 1;
            }

            coords[i][dimension_map[dim]] = coordinate[dim];
        }

    }

    // TODO: Verify that trace dimension is always 0
    auto size = vds_handle.access_manager().GetVolumeTracesBufferSize(npoints, VDSLevelOfDetailID::Level0);

    std::unique_ptr< char[] > data(new char[size]());

    auto request = vds_handle.access_manager().RequestVolumeTraces(
            (float*)data.get(),
            size,
            OpenVDS::Dimensions_012,
            VDSLevelOfDetailID::Level0,
            VDSChannelID::Amplitude,
            coords.get(),
            npoints,
            to_interpolation(interpolation_method),
            0
    );

    return finalize_request( request, "Failed to fetch fence from VDS", data, size );
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
    internal::VDSHandle vds_handle(url, credentials);

    nlohmann::json meta;
    meta["format"] = vds_handle.get_channel_format_string(internal::VDSChannelID::Amplitude);
    meta["crs"] = vds_handle.get_crs_string();

    auto bbox = internal::BoundingBox(vds_handle.layout());
    meta["boundingBox"]["ij"]   = bbox.index();
    meta["boundingBox"]["cdp"]  = bbox.world();
    meta["boundingBox"]["ilxl"] = bbox.annotation();

    for (int i = 2; i >= 0 ; i--) {
        meta["axis"].push_back( vds_handle.get_axis_metadata( i ) );
    }
    return internal::requestdata_from_dump( meta.dump() );
}

struct requestdata fetch_fence_metadata(
    std::string url,
    std::string credentials,
    size_t npoints
) {
    internal::VDSHandle vds_handle(url, credentials);

    nlohmann::json meta;
    meta["shape"] = vds_handle.get_shape_metadata(npoints);
    meta["format"] = vds_handle.get_channel_format_string(internal::VDSChannelID::Amplitude);

    return internal::requestdata_from_dump( meta.dump() );
}

} // namespace internal

/*
 * Below here the interface methods are defined.
 */


struct requestdata slice(
    const char* vds,
    const char* credentials,
    int lineno,
    Axis ax
) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return internal::fetch_slice(cube, cred, ax, lineno);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata slice_metadata(
    const char* vds,
    const char* credentials,
    Axis ax
) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return internal::fetch_slice_metadata(cube, cred, ax);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata fence(
    const char* vds,
    const char* credentials,
    enum CoordinateSystem coordinate_system,
    const float* coordinates,
    size_t npoints,
    enum InterpolationMethod interpolation_method
) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return internal::fetch_fence(
            cube, cred, coordinate_system, coordinates, npoints,
            interpolation_method);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata fence_metadata(
    const char* vds,
    const char* credentials,
    size_t npoints
) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return internal::fetch_fence_metadata(cube, cred, npoints);
    } catch (const std::exception& e) {
        return internal::handle_error(e);
    }
}

struct requestdata metadata(
    const char* vds,
    const char* credentials
) {
    try {
        std::string cube(vds);
        std::string cred(credentials);
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
