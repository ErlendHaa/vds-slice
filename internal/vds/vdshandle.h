#ifndef VDSHANDLE_H
#define VDSHANDLE_H

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "nlohmann/json.hpp"

#include <OpenVDS/OpenVDS.h>
#include <OpenVDS/KnownMetadata.h>
#include <OpenVDS/IJKCoordinateTransformer.h>

#include "coordinatetransformers.h"
#include "vds.h"

// Aliased types
using IntPointList = std::vector<std::pair<int, int > >;
using DoublePointList = std::vector<std::pair<double, double> >;

struct BoundingBox {
    IntPointList       index;
    IntPointList       annotation;
    DoublePointList world;
};

class AxisDescr {
    int dimension;

    int min;
    int max;
    int stride;
    int samples;

    const char * name;
    const char * unit;
};

// Implementation of Axis that knows its own string representation
// and which coordinate system it belongs too.
// TODO naming conflict with current Axis
struct Axis {
    Axis(CAxis axis) : axis(axis) {}

    CordinateSystem system() const {
        switch (this->value) {
            /* ... */
        }
    }

    std::string string() const { /* ... */ }
};

private:
    CAxis axis;
}

struct SeismicAxisMap {
    virtual int iline()  const = 0;
    virtual int xline()  const = 0;
    virtual int sample() const = 0;
    virtual int offset() const = 0; // Planing for a future of Prestack support
};

struct PostStackAxisMap : public SeismicAxisMap {
public:
    PostStackAxisMap(int i, int x, int s) : i(i), x(x), s(s) {};

    int iline()  const override { return this->i };
    int xline()  const override { return this->x };
    int sample() const override { return this->s };
    int offset() const override { throw std::not_implemented(""); };
private:
    int i;
    int x;
    int s;
};

/** Verifies e.g. annotations is present */
class SeismicVDSHandle {
public:
    SeismicVDSHandle(
        std::string           url,
        std::string           connection,
        struct Channel        default_channel,
        struct LOD            default_lod,
        struct SeismicAxisMap axis_map
    );

    struct SeismicVDSValidator { /* ... */ };

    // Maps from our Axis to a VDS axisDescriptor.
    AxisDescr   get_axis(Axis axis) const;
    BoundingBox get_bounding_box() const;
    std::string get_crs() const;
    std::string get_format(Channel ch = Channel::None) const;
private:
    struct SeismicAxisMap axis_map;
};

AxisDescr SeismicVDSHandle::get_axis(Axis axis) const {
    int axis;
    switch (axis.Value()) {
        case CAxis::I:
        case CAxis::Inline: {
            axis = this->axis_map.iline;
            break;
        }
        // ...
    }

    return {
        axis,
        this->layout.GetDimensionMin(axis),
        this->layout.GetDimensionMax(axis),
        // TODO stride
        this->layout.GetDimensionNumSamples(axis),
        // TODO add units and name
    }
}


// Does the actual validation for the PostStack class
struct PostStackValidator {
public:
    // Verify dimensions
    // Verify axis order
    bool validate(const PostStack& handle) {
    }

    bool units(const PostStack& handle) {
    }

private:
    // ...
};

/** Verifies dimension is 3D and axis order */
class PostStack : SeismicVDSHandle {
public:
    PostStack(std::string url, std::string conn)
        : SeismicVDSHandle(/* ... */, PostStackAxisMap( 2, 1, 0 ))  // hardcode axis map, this assumption will be validated
    {
        // TODO hardcode defaults for lod and channel and send to SeismicVDSHandle
        if not PostStackValidator().validate(this) throw std::runtime_error("");
    }

    requestdata slice(
        Axis axis,
        int lineno,
        LOD lod = LOD::None,
        Channel channel = Channel::None
    );

private:
    SubVolume slice_as_subvolume(int axis, int lineno);
    buffer    get_subvolume(SubVolume volume);
};

requestbuffer PostStack::slice(
    Axis axis,
    int lineno,
    LOD lod = LOD::Default,
    Channel channel = Channel::Default
) {
    AxisDescr description = this->get_axis(axis);
    CoordinateSystem = Axis(axis).system();

    if (system = CoordinateSystem::Annotation) {
        lineno = (lineno - axis_description.min) / axis_description.stride;
    }

    if (not PostStackValidator(this).units())
        throw std::invalid_argument("");

    if (not this->slice_bounds_check(description.axis, lineno))
        throw std::out_of_bounds("");

    SubVolume subvolume = this->slice_as_subvolume(axis, lineno);
    return this->get_subvolume(subvolume, lod, channel);
}

requestbuffer slice(Axis axis, int lineno) {
    PostStackHandle poststack(url, conn);
    return poststack.slice(Axis, lineno);
}

class VDSHandle {

    private:

        enum VDSLevelOfDetailID {
            Level0 = 0
        };

        OpenVDS::ScopedVDSHandle handle_;
        OpenVDS::Error error_;
        OpenVDS::VolumeDataAccessManager access_manager_;
        OpenVDS::VolumeDataLayout const * layout_;
        OpenVDS::IJKCoordinateTransformer ijk_coordinate_transformer_;
        int seismic_channel_id_;

        const std::string seismic_channel_name_{"Amplitude"};

        static OpenVDS::InterpolationMethod to_interpolation( const InterpolationMethod interpolation);

        int convert_ijk_to_voxel_axis_id( const int ijk_axis_id ) const;

        int convert_voxel_to_ijk_axis_id( const int voxel_axis_id ) const;

        int get_max_voxel_index( const int voxel_dimension ) const;

        static requestdata requestdata_from_requested_data( std::unique_ptr< char[] >  &data,
                                                            const std::size_t size );

        template<typename REQUEST_TYPE>
        static requestdata finalize_request( const std::shared_ptr<REQUEST_TYPE>& request,
                                             const std::string message,
                                             std::unique_ptr< char[] >& data,
                                             const std::size_t size );

        void validate_dimension() const;

        void check_axis( const int dimension, char const * const  expected_name ) const;

        void validate_axes_order() const;

        void validate_annotations_are_defined() const;

        void validate_data_store() const;

    IntPointList get_bounding_box_in_index_coordinates() const;

    DoublePointList convert_ijk_to_world(const IntPointList& points_as_ijk) const;

    IntPointList convert_ijk_to_annotation(const IntPointList& points_as_ijk) const;

        int get_voxel_axis_id_of( const Axis axis ) const;

        requestdata request_volume_trace( const std::unique_ptr<float[][OpenVDS::Dimensionality_Max]> &coordinates,
                                          const std::size_t npoints,
                                          const InterpolationMethod interpolation_method );

        std::unique_ptr< float[][OpenVDS::Dimensionality_Max] >
        get_fence_request_coordinates( const enum CoordinateSystem coordinate_system,
                                       const float* coordinates,
                                       const size_t npoints,
                                       const enum InterpolationMethod interpolation_method  ) const;

    public:

        VDSHandle( std::string url, std::string credentials);

        std::string get_format_string_of_seismic_channel() const;

        std::string get_crs_string() const;

        nlohmann::json get_axis_metadata( const Axis axis, const int requested_dimension ) const;

        nlohmann::json get_axis_metadata( const int requested_dimension ) const;

        nlohmann::json get_shape_metadata( const int npoints ) const;

        BoundingBox get_bounding_box() const;

        requestdata get_slice_of( const Axis axis, const int line_number );

        requestdata get_fence_of( const enum CoordinateSystem coordinate_system,
                                  const float* coordinates,
                                  const size_t npoints,
                                  const enum InterpolationMethod interpolation_method );

        requestdata request_volume_subset( const internal::VoxelBounds& voxel_bounds );

};

#endif /* VDSHANDLE_H */
