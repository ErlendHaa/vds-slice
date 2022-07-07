#include "vds.h"

#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <algorithm>
#include <string>
#include <stdexcept>

#include <OpenVDS/OpenVDS.h>
#include <OpenVDS/KnownMetadata.h>
#include <OpenVDS/IJKCoordinateTransformer.h>

using namespace std;

void vdsbuffer_delete(struct vdsbuffer* buf) {
    if (!buf)
        return;

    delete[] buf->data;
    delete[] buf->err;
    *buf = vdsbuffer {};
}

enum coord_system {
    INDEX      = 0,
    ANNOTATION = 1,
};

int todimension(axis ax) {
    switch (ax) {
        case X:
        case ILINE:
            return 0;
        case Y:
        case XLINE:
            return 1;
        case Z:
        case DEPTH:
        case TIME:
        case SAMPLE:
            return 2;
    }
}

coord_system tosystem(axis ax) {
    switch (ax) {
        case X:
        case Y:
        case Z:
            return INDEX;
        case ILINE:
        case XLINE:
        case DEPTH:
        case TIME:
        case SAMPLE:
            return ANNOTATION;
    }
}

const std::string axis_tostring(axis ax) {
    switch (ax) {
        case X:      return std::string( OpenVDS::KnownAxisNames::X()    );
        case Y:      return std::string( OpenVDS::KnownAxisNames::Y()    );
        case Z:      return std::string( OpenVDS::KnownAxisNames::Z()    );
        case ILINE:  return std::string( OpenVDS::KnownAxisNames::Inline()    );
        case XLINE:  return std::string( OpenVDS::KnownAxisNames::Crossline() );
        case DEPTH:  return std::string( OpenVDS::KnownAxisNames::Depth()     );
        case TIME:   return std::string( OpenVDS::KnownAxisNames::Time()      );
        case SAMPLE: return std::string( OpenVDS::KnownAxisNames::Sample()    );
        default: {
            throw std::runtime_error("Unhandled axis");
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
bool zaxisvalidation(axis ax, const char* zunit) {
    /* Define some convenient lookup tables for units */
    static const std::array< const char*, 3 > depthunits = {
        OpenVDS::KnownUnitNames::Meter(),
        OpenVDS::KnownUnitNames::Foot(),
        OpenVDS::KnownUnitNames::USSurveyFoot()
    };

    static const std::array< const char*, 2 > timeunits = {
        OpenVDS::KnownUnitNames::Millisecond(),
        OpenVDS::KnownUnitNames::Second()
    };

    static const std::array< const char*, 1 > sampleunits = {
        OpenVDS::KnownUnitNames::Unitless(),
    };

    auto isoneof = [zunit](const char* x) {
        return !std::strcmp(x, zunit);
    };

    switch (ax) {
        case X:
        case Y:
        case Z:
        case ILINE:
        case XLINE:
            return true;
        case DEPTH:
            return std::any_of(depthunits.begin(), depthunits.end(), isoneof);
        case TIME:
            return std::any_of(timeunits.begin(), timeunits.end(), isoneof);
        case SAMPLE:
            return std::any_of(sampleunits.begin(), sampleunits.end(), isoneof);
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
};

void axisvalidation(axis ax, const OpenVDS::VolumeDataLayout* layout) {
    /* This assumes that the Z-axis is always the first one in the VDS */
    auto zaxis = layout->GetAxisDescriptor(0);
    const char* zunit = zaxis.GetUnit();
    if (not zaxisvalidation(ax, zunit)) {
        std::string msg = "Unable to use " + axis_tostring(ax);
        msg += " on cube with depth units: " + std::string(zunit);
        throw std::runtime_error(msg);
    }
}

/*
 * Convert target dimension/axis + lineno to VDS voxel coordinates.
 */
void set_voxels(
    coord_system sys,
    int dimension,
    int lineno,
    const OpenVDS::VolumeDataLayout *layout,
    int (&voxelmin)[OpenVDS::VolumeDataLayout::Dimensionality_Max],
    int (&voxelmax)[OpenVDS::VolumeDataLayout::Dimensionality_Max]
) {

    /*
     * This seams to be a very roundabout way of converting from
     * dimension/slice to voxel coordinates which might be improved upon in the
     * future.
     *
     * 1) Init min/max to span from 0 - samplemax (Voxel)
     * 2) Convert to index or annotation
     * 3) Update the target dimension
     * 4) Convert back to voxel
     */
    auto vmin = OpenVDS::IntVector3 { 0, 0, 0 };
    auto vmax = OpenVDS::IntVector3 {
        layout->GetDimensionNumSamples(0) - 1,
        layout->GetDimensionNumSamples(1) - 1,
        layout->GetDimensionNumSamples(2) - 1
    };

    auto transformer = OpenVDS::IJKCoordinateTransformer(layout);

    switch (sys) {
        case ANNOTATION: {
            auto amin = transformer.VoxelIndexToAnnotation(vmin);
            auto amax = transformer.VoxelIndexToAnnotation(vmax);
            amin[dimension] = lineno;
            amax[dimension] = lineno + 1;
            vmin = transformer.AnnotationToVoxelIndex(amin);
            vmax = transformer.AnnotationToVoxelIndex(amax);
            break;
        }
        case INDEX: {
            auto imin = transformer.VoxelIndexToIJKIndex(vmin);
            auto imax = transformer.VoxelIndexToIJKIndex(vmax);
            imin[dimension] = lineno;
            imax[dimension] = lineno + 1;
            vmin = transformer.IJKIndexToVoxelIndex(imin);
            vmax = transformer.IJKIndexToVoxelIndex(imax);
            break;
        }
    }

    for (int i = 0; i < 3; i++) {
        voxelmin[i] = vmin[i];
        voxelmax[i] = vmax[i];
    }
}

struct vdsbuffer fetch_slice(
    std::string url,
    std::string credentials,
    axis ax,
    int lineno
) {
    OpenVDS::Error error;
    OpenVDS::VDSHandle handle = OpenVDS::Open(url, credentials, error);

    if(error.code != 0) {
        throw std::runtime_error("Could not open VDS: " + error.string);
    }

    auto access = OpenVDS::GetAccessManager(handle);
    auto const *layout = access.GetVolumeDataLayout();

    axisvalidation(ax, layout);

    int vmin[OpenVDS::Dimensionality_Max] = { 0, 0, 0, 0, 0, 0};
    int vmax[OpenVDS::Dimensionality_Max] = { 1, 1, 1, 1, 1, 1};
    auto dimension = todimension(ax);
    auto sys = tosystem(ax);
    set_voxels(sys, dimension, lineno, layout, vmin, vmax);

    auto format = layout->GetChannelFormat(0);
    auto size = access.GetVolumeSubsetBufferSize(vmin, vmax, format, 0, 0);


    auto *data = new char[size];
    auto request = access.RequestVolumeSubset(
        data,
        size,
        OpenVDS::Dimensions_012,
        0,
        0,
        vmin,
        vmax,
        format
    );

    request.get()->WaitForCompletion();

    vdsbuffer buffer{};
    buffer.size = size;
    buffer.data = data;

    return buffer;
}

struct vdsbuffer slice(
    const char* vds,
    const char* credentials,
    int lineno,
    axis ax
) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return fetch_slice(cube, cred, ax, lineno);
    } catch (const std::exception& e) {
        vdsbuffer buf {};
        buf.err = new char[std::strlen(e.what()) + 1];
        std::strcpy(buf.err, e.what());
        return buf;
    }
}
