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

        IntPointList get_bounding_box_in_index_coordinates();

        DoublePointList convert_ijk_to_world( const IntPointList& points_as_ijk );

        IntPointList convert_ijk_to_annotation( const IntPointList& points_as_ijk );

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

        BoundingBox get_bounding_box();

        requestdata get_slice_of( const Axis axis, const int line_number );

        requestdata get_fence_of( const enum CoordinateSystem coordinate_system,
                                  const float* coordinates,
                                  const size_t npoints,
                                  const enum InterpolationMethod interpolation_method );

        requestdata request_volume_subset( const internal::VoxelBounds& voxel_bounds );

};

#endif /* VDSHANDLE_H */