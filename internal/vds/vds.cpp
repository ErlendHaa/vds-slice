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

#include "coordinatetransformers.h"

using namespace std;

namespace internal {

enum VDSLevelOfDetailID {
    Level0 = 0
};


static const std::unique_ptr<CoordinateToVDSVoxelTransformer> get_coordinate_system( const Axis axis,
                                                                                OpenVDS::VolumeDataLayout const * vds_layout_)  {
    switch (axis) {
        case I:
        case J:
        case K:
            return std::make_unique<IndexToVDSVoxel>(axis, vds_layout_);
        case INLINE:
        case CROSSLINE:
        case DEPTH:
        case TIME:
        case SAMPLE:
            return std::make_unique<AnnotationToVDSVoxel>(axis, vds_layout_);
    }
    return std::make_unique<CDPToVDSVoxel>(axis, vds_layout_);
}


requestdata requestdata_from_dump( const nlohmann::json::string_t& dump ) {
    requestdata tmp{ new char[dump.size()], nullptr, dump.size() };
    std::copy(dump.begin(), dump.end(), tmp.data);
    return tmp;
}


class BoundingBox {
public:
    explicit BoundingBox(
        const OpenVDS::VolumeDataLayout &layout
    ) : layout(layout)
    {
        transformer = OpenVDS::IJKCoordinateTransformer(&layout);
    }

    std::vector< std::pair<int, int> >       index()      const noexcept (true);
    std::vector< std::pair<int, int> >       annotation() const noexcept (true);
    std::vector< std::pair<double, double> > world()      const noexcept (true);
private:
    OpenVDS::IJKCoordinateTransformer transformer;
    const OpenVDS::VolumeDataLayout &layout;
};


std::vector< std::pair<int, int> > BoundingBox::index() const noexcept (true) {
    const auto ils = layout.GetDimensionNumSamples(VDSAxisID::Inline) - 1;
    const auto xls = layout.GetDimensionNumSamples(VDSAxisID::Crossline) - 1;

    return { {0, 0}, {ils, 0}, {ils, xls}, {0, xls} };
}

std::vector< std::pair<double, double> > BoundingBox::world() const noexcept (true) {
    std::vector< std::pair<double, double> > world_points;

    const auto points = this->index();
    std::for_each(points.begin(), points.end(),
        [&](const std::pair<int, int>& point) {
            const auto p = this->transformer.IJKIndexToWorld(
                { point.first, point.second, 0 }
            );
            world_points.emplace_back(p[0], p[1]);
        }
    );

    return world_points;
};

std::vector< std::pair<int, int> > BoundingBox::annotation() const noexcept (true) {
    auto points = this->index();
    std::transform(points.begin(), points.end(), points.begin(),
        [this](const std::pair<int, int>& point) {
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

class VDSHandle {

    private:
        OpenVDS::ScopedVDSHandle handle_;
        OpenVDS::Error error_;
        OpenVDS::VolumeDataAccessManager access_manager_;
        OpenVDS::VolumeDataLayout const * layout_;
        OpenVDS::IJKCoordinateTransformer ijk_coordinate_transformer_;
        int seismic_channel_id_;

        const std::string seismic_channel_name_{"Amplitude"};

        static OpenVDS::InterpolationMethod to_interpolation( const InterpolationMethod interpolation) {
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

        static requestdata requestdata_from_requested_data( std::unique_ptr< char[] >  &data,
                                                            const std::size_t size ) {
            requestdata tmp{ data.get(), nullptr, size };
            data.release();
            return tmp;
        }

        template<typename REQUEST_TYPE>
        static requestdata finalize_request( const std::shared_ptr<REQUEST_TYPE>& request,
                                             const std::string message,
                                             std::unique_ptr< char[] >& data,
                                             const std::size_t size ) {

            const bool success = request.get()->WaitForCompletion();
            if( not success ) {
                throw std::runtime_error(message);
            }

            return requestdata_from_requested_data( data, size );
        }

        void validate_dimension() const {
            if (layout_->GetDimensionality() != 3) {
                throw std::runtime_error(
                    "Unsupported VDS, expected 3 dimensions, got " +
                    std::to_string(layout_->GetDimensionality())
                );
            }
        }

        void check_axis( const int dimension, char const * const  expected_name ) const
        {
            char const * const actual_name = this->layout_->GetDimensionName(dimension);
            if (std::strcmp(actual_name , expected_name)) {
                const std::string msg = std::string("Unsupported axis ordering in VDS for axis nr. ")
                                        + std::to_string(dimension)
                                        + " named "
                                        + actual_name
                                        + ", expected "
                                        + expected_name;
                throw std::runtime_error(msg);
            }
        }

        void validate_axes_order() const {
            check_axis( VDSAxisID::Inline, OpenVDS::KnownAxisNames::Inline() );
            check_axis( VDSAxisID::Crossline, OpenVDS::KnownAxisNames::Crossline() );

            const auto z = this->layout_->GetDimensionName(VDSAxisID::DepthSampleTime);

            auto isoneof = [z](const char* x) {
                return !std::strcmp(x, z);
            };

            using namespace constants::annotation;
            if ( not std::any_of(axis_labels.begin(),
                                 axis_labels.end(),
                                 isoneof) )
            {
                const char* actual_name = this->layout_->GetDimensionName(VDSAxisID::DepthSampleTime);
                const std::string msg = std::string("Unsupported axis ordering in VDS for axis nr. 0 named ")
                                        + actual_name
                                        + ", expected "
                                        + "Depth, Time, or Sample";
                throw std::runtime_error(msg);
            }
        }

        void validate_annotations_are_defined() const {
            if (not this->ijk_coordinate_transformer_.AnnotationsDefined()) {
                throw std::runtime_error("VDS doesn't define annotations");
            }
        }

        void validate_data_store() const {
            validate_dimension();
            validate_axes_order();
            validate_annotations_are_defined();
        }

        int get_voxel_axis_id_of( const Axis axis ) const {
            //this->validate_request_axis(axis);
            const auto coordinate_system = get_coordinate_system(axis, this->layout_);
            return this->convert_ijk_to_voxel_axis_id(coordinate_system->axis_to_vds_dimension());
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
            seismic_channel_id_ = layout_->GetChannelIndex( seismic_channel_name_.c_str() );

            validate_data_store();
        }

        std::string get_format_string_of_seismic_channel() const {
            using namespace OpenVDS;
            const VolumeDataFormat format = layout_->GetChannelFormat(seismic_channel_id_);
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
            const auto crs = OpenVDS::KnownMetadata::SurveyCoordinateSystemCRSWkt();
            return layout_->GetMetadataString(crs.GetCategory(), crs.GetName());
        }

        nlohmann::json get_axis_metadata( const Axis axis, const int requested_dimension ) const {
            const int vdim = this->get_voxel_axis_id_of(axis);

            //TODO: Make this nicer and more expressive
            //We assemple an array of all dimensions that
            std::vector< int > non_voxel_dimensions;
            for (int i = 2; i >= 0; --i) {
                if (i == vdim) continue;
                non_voxel_dimensions.push_back(i);
            }

            const int mapped_id = non_voxel_dimensions[requested_dimension];

            return get_axis_metadata( mapped_id );
        }

        nlohmann::json get_axis_metadata( const int requested_dimension ) const {
            return nlohmann::json {
                { "annotation", layout_->GetDimensionName(requested_dimension)       },
                { "min",        layout_->GetDimensionMin(requested_dimension)        },
                { "max",        layout_->GetDimensionMax(requested_dimension)        },
                { "samples",    layout_->GetDimensionNumSamples(requested_dimension) },
                { "unit",       layout_->GetDimensionUnit(requested_dimension)       },
            };
        }

        nlohmann::json get_shape_metadata( const int npoints ) const {
            return nlohmann::json::array({ npoints,
                                           layout_->GetDimensionNumSamples(VDSAxisID::DepthSampleTime)});
        }

        BoundingBox get_bounding_box() {
            return BoundingBox(*this->layout_);
        }

        requestdata get_slice_of( const Axis axis, const int line_number ) {
            const auto coordinate_system = get_coordinate_system(axis, this->layout_);
            VoxelBounds vb = coordinate_system->get_slice_bounds_from( line_number, this->get_voxel_axis_id_of(axis) );
            return this->request_volume_subset( std::move(vb) );
        }

        unique_ptr< float[][OpenVDS::Dimensionality_Max] >
        get_fence_request_coordinates( const enum CoordinateSystem coordinate_system,
                                       const float* coordinates,
                                       const size_t npoints,
                                       const enum InterpolationMethod interpolation_method  ) const {
            const OpenVDS::IntVector3 dimension_map =
                this->layout_->GetVDSIJKGridDefinitionFromMetadata().dimensionMap;

            unique_ptr< float[][OpenVDS::Dimensionality_Max] > coords(
                new float[npoints][OpenVDS::Dimensionality_Max]{{0}}
            );

            auto transform_coordinate = [&] (const float x, const float y) {
                switch (coordinate_system) {
                    case INDEX:
                        return OpenVDS::Vector<double, 3> {x, y, 0};
                    case ANNOTATION:
                        return this->ijk_coordinate_transformer_.AnnotationToIJKPosition({x, y, 0});
                    case CDP:
                        return this->ijk_coordinate_transformer_.WorldToIJKPosition({x, y, 0});
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
                    const auto max = this->get_max_voxel_index( voxel ) - 0.5;
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
            return coords;
        }

        requestdata request_volume_trace( const unique_ptr<float[][OpenVDS::Dimensionality_Max]> &coordinates,
                                          const std::size_t npoints,
                                          const InterpolationMethod interpolation_method ) {
            // TODO: Verify that trace dimension is always 0
            const std::size_t size = this->access_manager_.GetVolumeTracesBufferSize(npoints, VDSLevelOfDetailID::Level0);

            std::unique_ptr< char[] > data(new char[size]());

            const auto request = access_manager_.RequestVolumeTraces(
                    (float*)data.get(),
                    size,
                    OpenVDS::Dimensions_012,
                    VDSLevelOfDetailID::Level0,
                    seismic_channel_id_,
                    coordinates.get(),
                    npoints,
                    to_interpolation(interpolation_method),
                    0
            );

            return finalize_request( request, "Failed to fetch fence from VDS", data, size );
        }

        requestdata request_volume_subset( const VoxelBounds& voxel_bounds ) {

            const auto format = layout_->GetChannelFormat(this->seismic_channel_id_);
            const int size = access_manager_.GetVolumeSubsetBufferSize(
                voxel_bounds.lower,
                voxel_bounds.upper,
                format,
                VDSLevelOfDetailID::Level0,
                this->seismic_channel_id_);

            std::unique_ptr< char[] > data(new char[size]());
            const auto request = access_manager_.RequestVolumeSubset(
                data.get(),
                size,
                OpenVDS::Dimensions_012,
                VDSLevelOfDetailID::Level0,
                this->seismic_channel_id_,
                voxel_bounds.lower,
                voxel_bounds.upper,
                format
            );

            return finalize_request( request, "Failed to fetch slice from VDS", data, size );
        }

};

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

    const auto request_coordinates =
        vds_handle.get_fence_request_coordinates( coordinate_system,
                                                  coordinates,
                                                  npoints,
                                                  interpolation_method);

    return vds_handle.request_volume_trace( request_coordinates, npoints, interpolation_method );
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
    meta["format"] = vds_handle.get_format_string_of_seismic_channel();
    meta["crs"] = vds_handle.get_crs_string();

    const auto bbox = vds_handle.get_bounding_box();
    meta["boundingBox"]["ij"]   = bbox.index();
    meta["boundingBox"]["cdp"]  = bbox.world();
    meta["boundingBox"]["ilxl"] = bbox.annotation();

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
    internal::VDSHandle vds_handle(url, credentials);

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
