#ifndef COORDINATETRANSFORMERS_H
#define COORDINATETRANSFORMERS_H

#include <array>
#include <vector>

#include <OpenVDS/OpenVDS.h>

#include "vds.h"

namespace internal {

enum VDSAxisID {
    DepthSampleTime=0,
    Crossline=1,
    Inline=2,
};

struct VoxelBounds {
    int lower[OpenVDS::VolumeDataLayout::Dimensionality_Max]{0, 0, 0, 0, 0, 0};
    int upper[OpenVDS::VolumeDataLayout::Dimensionality_Max]{1, 1, 1, 1, 1, 1};
};


namespace constants {
namespace annotation {

enum Index{
    depth=0,
    time=1,
    sample=2
};

struct AxisUnitCombination {
    const char* label;
    const std::vector<char const * const> units;

    AxisUnitCombination( char const * const l, const std::vector<char const * const> u) : label(l), units(std::move(u))
    {}
};

extern const std::array<const AxisUnitCombination, 3> label_unit_combinations;
extern const std::array<char const* const, 3> axis_labels;

} //namespace annotation
} //namespace constants

class CoordinateToVDSVoxelTransformer {

    const CoordinateSystem coordinate_system_label_;

    protected:
    const Axis axis_;
    VoxelBounds voxel_limits_;
    const std::vector<char const * const> legal_units_;

    public:
    CoordinateToVDSVoxelTransformer( const Axis axis,
                                     OpenVDS::VolumeDataLayout const * vds_layout_,
                                     const CoordinateSystem coordinate_system_label );

    CoordinateSystem get_coordinate_system_label() const;

    virtual const std::string get_name() const = 0;

    virtual VoxelBounds get_slice_bounds_from( const int lineno,
                                               const int vdim ) const = 0;

    virtual int axis_to_vds_dimension() const = 0;

    virtual ~CoordinateToVDSVoxelTransformer();
};



class IndexToVDSVoxel : public CoordinateToVDSVoxelTransformer {

    public:
    IndexToVDSVoxel(const Axis axis,
                    OpenVDS::VolumeDataLayout const * vds_layout_);

    const std::string get_name() const override;

    VoxelBounds get_slice_bounds_from( const int lineno,
                                       const int vdim ) const override;

    int axis_to_vds_dimension() const override;

    ~IndexToVDSVoxel() override;
};


class AnnotationToVDSVoxel : public CoordinateToVDSVoxelTransformer {

    std::array<int, 3> dimensions_min_;
    std::array<int, 3> dimensions_max_;
    const char* unit_of_depth_sample_time_axis_;

    static const std::string axis_tostring(const Axis ax);

    protected:

    //TODO: Why is this only checking units if axis is the third coordinate axis?
    void validate_units_of_z_axis() const;

    public:
    AnnotationToVDSVoxel(const Axis axis,
                         OpenVDS::VolumeDataLayout const * vds_layout_);

    const std::string get_name() const override;

    virtual VoxelBounds get_slice_bounds_from( const int lineno,
                                               const int vdim ) const override;

    int axis_to_vds_dimension() const override;

    ~AnnotationToVDSVoxel() override;
};


class CDPToVDSVoxel : public CoordinateToVDSVoxelTransformer {
    public:

    CDPToVDSVoxel(const Axis axis, OpenVDS::VolumeDataLayout const * vds_layout_);

    const std::string get_name() const override;

    VoxelBounds get_slice_bounds_from(  [[maybe_unused]] const int lineno,
                                        [[maybe_unused]] const int vdim ) const override;

    int axis_to_vds_dimension() const override;

    ~CDPToVDSVoxel() override;
};


} // namespace internal

#endif
