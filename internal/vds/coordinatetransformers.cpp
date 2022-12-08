#include "coordinatetransformers.h"

#include <OpenVDS/KnownMetadata.h>

#include "vds.h"

namespace internal {

namespace constants {
namespace annotation {

const std::array<const AxisUnitCombination, 3> label_unit_combinations{
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Depth(),
        std::vector<char const * const>{
            OpenVDS::KnownUnitNames::Meter(),
            OpenVDS::KnownUnitNames::Foot(),
            OpenVDS::KnownUnitNames::USSurveyFoot()
        }
    ),
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Time(),
        std::vector<char const * const>{
            OpenVDS::KnownUnitNames::Millisecond(),
            OpenVDS::KnownUnitNames::Second()
        }
    ),
    AxisUnitCombination(
        OpenVDS::KnownAxisNames::Sample(),
        std::vector<char const * const>{
            OpenVDS::KnownUnitNames::Unitless()
        }
    )
};

const std::array<char const* const, 3> axis_labels{
    label_unit_combinations[Index::depth].label,
    label_unit_combinations[Index::time].label,
    label_unit_combinations[Index::sample].label,
};

} //namespace annotation
} //namespace constants


/*
 * class CoordinateToVDSVoxelTransformer
 */
CoordinateToVDSVoxelTransformer::CoordinateToVDSVoxelTransformer(
    const Axis axis,
    OpenVDS::VolumeDataLayout const* vds_layout_,
    const CoordinateSystem coordinate_system_label)
    : coordinate_system_label_(coordinate_system_label), axis_(axis) {
    for (std::size_t i = 0; i < 3; ++i) {
        this->voxel_limits_.upper[i] = vds_layout_->GetDimensionNumSamples(i);
    }
}

CoordinateSystem CoordinateToVDSVoxelTransformer::get_coordinate_system_label() const {
    return coordinate_system_label_;
}

CoordinateToVDSVoxelTransformer::~CoordinateToVDSVoxelTransformer() {}


/*
 * class IndexToVDSVoxel
 */
IndexToVDSVoxel::IndexToVDSVoxel(const Axis axis,
                                 OpenVDS::VolumeDataLayout const * vds_layout_) :
    CoordinateToVDSVoxelTransformer(axis, vds_layout_, CoordinateSystem::INDEX)  {
}

const std::string IndexToVDSVoxel::get_name() const {
    return std::string("INDEX");
}

VoxelBounds IndexToVDSVoxel::get_slice_bounds_from( const int lineno,
                                    const int vdim ) const {

    const int min = 0;
    const int max = this->voxel_limits_.upper[vdim] - 1;

    if (lineno < min || lineno > max) {
        throw std::runtime_error(
            "Invalid lineno: " + std::to_string(lineno) +
            ", valid range: [" + std::to_string(min) +
            ":" + std::to_string(max) +
            ":1]"
        );
    }

    VoxelBounds bounds = this->voxel_limits_;
    bounds.lower[vdim] = lineno;
    bounds.upper[vdim] = lineno+1;
    return bounds;
}

int IndexToVDSVoxel::axis_to_vds_dimension() const {
    switch (this->axis_) {
        case Axis::I:
            return 0;
        case Axis::J:
            return 1;
        case Axis::K:
            return 2;
        case Axis::INLINE:
        case Axis::CROSSLINE:
        case Axis::DEPTH:
        case Axis::TIME:
        case Axis::SAMPLE:
            throw std::runtime_error("IndexToVDSVoxel: Unhandled axis  " +
                                     std::to_string(this->axis_));
    }
}

IndexToVDSVoxel::~IndexToVDSVoxel() {}


/*
 * class AnnotationToVDSVoxel
 */
const std::string AnnotationToVDSVoxel::axis_tostring(const Axis ax) {
    switch (ax) {
        case Axis::INLINE:
            return std::string(OpenVDS::KnownAxisNames::Inline());
        case Axis::CROSSLINE:
            return std::string(OpenVDS::KnownAxisNames::Crossline());
        case Axis::DEPTH:
            return std::string(OpenVDS::KnownAxisNames::Depth());
        case Axis::TIME:
            return std::string(OpenVDS::KnownAxisNames::Time());
        case Axis::SAMPLE:
            return std::string(OpenVDS::KnownAxisNames::Sample());
        case Axis::I:
        case Axis::J:
        case Axis::K:
            throw std::runtime_error("Unhandled axis");
    }
}

void AnnotationToVDSVoxel::validate_units_of_z_axis() const {

    auto isoneof = [this](char const * const x) {
        return !std::strcmp(x, this->unit_of_depth_sample_time_axis_);
    };


    bool is_valid_axis_unit = false;
    using namespace constants::annotation;
    switch (this->axis_) {
        case Axis::DEPTH:
            is_valid_axis_unit = std::any_of(
                label_unit_combinations[Index::depth].units.begin(),
                label_unit_combinations[Index::depth].units.end(), isoneof);
            break;
        case Axis::TIME:
            is_valid_axis_unit = std::any_of(
                label_unit_combinations[Index::time].units.begin(),
                label_unit_combinations[Index::time].units.end(), isoneof);
            break;
        case Axis::SAMPLE:
            is_valid_axis_unit = std::any_of(
                label_unit_combinations[Index::sample].units.begin(),
                label_unit_combinations[Index::sample].units.end(), isoneof);
            break;
        case Axis::CROSSLINE:
        case Axis::INLINE:
            is_valid_axis_unit = true;
            break;
        case Axis::I:
        case Axis::J:
        case Axis::K:
            // TODO: This case should never happen. Verify thins and potentially
            // update/improve code.
            throw std::runtime_error("AnnotationToVDSVoxel: Unhandled axis");
    }

    if (not is_valid_axis_unit) {
        std::string msg = "Unable to use " + axis_tostring(this->axis_);
        msg += " on cube with depth units: " + std::string(this->unit_of_depth_sample_time_axis_);
        throw std::runtime_error(msg);
    }
}

AnnotationToVDSVoxel::AnnotationToVDSVoxel(
    const Axis axis,
    OpenVDS::VolumeDataLayout const* vds_layout_)
    : CoordinateToVDSVoxelTransformer(axis,
                                      vds_layout_,
                                      CoordinateSystem::ANNOTATION) {
    for (std::size_t i = 0; i < 3; ++i) {
        dimensions_min_[i] = vds_layout_->GetDimensionMin(i);
        dimensions_max_[i] = vds_layout_->GetDimensionMax(i);
    }
    unit_of_depth_sample_time_axis_ = vds_layout_->GetDimensionUnit(VDSAxisID::DepthSampleTime);
    validate_units_of_z_axis();
}


const std::string AnnotationToVDSVoxel::get_name() const {
    return std::string("ANNOTATION");
}

VoxelBounds AnnotationToVDSVoxel::get_slice_bounds_from( const int lineno,
                                            const int vdim ) const {
    const int min      = this->dimensions_min_[vdim];
    const int max      = this->dimensions_max_[vdim];
    const int nsamples = this->voxel_limits_.upper[vdim];

    const int stride = (max - min) / (nsamples - 1);

    if (lineno < min || lineno > max || (lineno - min) % stride) {
        throw std::runtime_error(
            "Invalid lineno: " + std::to_string(lineno) +
            ", valid range: [" + std::to_string(min) +
            ":" + std::to_string(max) +
            ":" + std::to_string(stride) + "]"
        );
    }

    const int voxelline = (lineno - min) / stride;
    VoxelBounds bounds = this->voxel_limits_;
    bounds.lower[vdim] = voxelline;
    bounds.upper[vdim] = voxelline+1;
    return bounds;
}

int AnnotationToVDSVoxel::axis_to_vds_dimension() const {
    switch (this->axis_) {
        case Axis::INLINE:
            return 0;
        case Axis::CROSSLINE:
            return 1;
        case Axis::DEPTH:
        case Axis::TIME:
        case Axis::SAMPLE:
            return 2;
        case Axis::I:
        case Axis::J:
        case Axis::K:
            throw std::runtime_error("AnnotationToVDSVoxel: Unhandled axis " +
                                     std::to_string(this->axis_));
    }
}

AnnotationToVDSVoxel::~AnnotationToVDSVoxel() {

}


/*
 * class CDPToVDSVoxel
 */
CDPToVDSVoxel::CDPToVDSVoxel(const Axis axis, OpenVDS::VolumeDataLayout const * vds_layout_) :
    CoordinateToVDSVoxelTransformer(axis, vds_layout_, CoordinateSystem::CDP ) {
}

const std::string CDPToVDSVoxel::get_name() const {
    return std::string("CDP");
}

VoxelBounds CDPToVDSVoxel::get_slice_bounds_from(  [[maybe_unused]] const int lineno,
                                    [[maybe_unused]] const int vdim ) const {
    throw std::runtime_error("Unhandled coordinate system");
}

int CDPToVDSVoxel::axis_to_vds_dimension() const {
    throw std::runtime_error("CDPToVDSVoxel: Unhandled axis " + std::to_string(this->axis_));
}

CDPToVDSVoxel::~CDPToVDSVoxel() {}


} // namespace internal
