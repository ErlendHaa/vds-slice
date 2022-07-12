#include "vds.h"

#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <algorithm>
#include <string>
#include <stdexcept>

#include <OpenVDS/OpenVDS.h>
#include <OpenVDS/KnownMetadata.h>

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

int axis_todim(axis ax) {
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

coord_system axis_tosystem(axis ax) {
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
bool unitvalidation(axis ax, const char* zunit) {
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
bool axis_order_validation(axis ax, const OpenVDS::VolumeDataLayout *layout) {
    if (std::strcmp(layout->GetDimensionName(2), OpenVDS::KnownAxisNames::Inline())) {
        return false;
    }

    if (std::strcmp(layout->GetDimensionName(1), OpenVDS::KnownAxisNames::Crossline())) {
        return false;
    }

    auto z = layout->GetDimensionName(0);

    /* Define some convenient lookup tables for units */
    static const std::array< const char*, 3 > depth = {
        OpenVDS::KnownAxisNames::Depth(),
        OpenVDS::KnownAxisNames::Time(),
        OpenVDS::KnownAxisNames::Sample()
    };

    auto isoneof = [z](const char* x) {
        return !std::strcmp(x, z);
    };

    return std::any_of(depth.begin(), depth.end(), isoneof);
}


void axisvalidation(axis ax, const OpenVDS::VolumeDataLayout* layout) {
    if (not axis_order_validation(ax, layout)) {
        std::string msg = "Unsupported axis ordering in VDS, expected ";
        msg += "Depth/Time/Sample, Crossline, Inline";
        throw std::runtime_error(msg);
    }

    auto zaxis = layout->GetAxisDescriptor(0);
    const char* zunit = zaxis.GetUnit();
    if (not unitvalidation(ax, zunit)) {
        std::string msg = "Unable to use " + axis_tostring(ax);
        msg += " on cube with depth units: " + std::string(zunit);
        throw std::runtime_error(msg);
    }
}

int dim_tovoxel(int dimension, const OpenVDS::VolumeDataLayout *layout) {
    /*
     * For now assume that the axis order is always depth/time/sample,
     * crossline, inline. This should be checked elsewhere.
     */
    switch (dimension) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 0;
        default: {
            throw std::runtime_error("Unhandled axis");
        }
    }
}

int lineno_tovoxel(int lineno, int vdim, const OpenVDS::VolumeDataLayout *layout) {
    auto min      = layout->GetDimensionMin(vdim);
    auto max      = layout->GetDimensionMax(vdim);
    auto nsamples = layout->GetDimensionNumSamples(vdim);

    auto stride    = (max - min) / (nsamples - 1);
    auto voxelline = (lineno - min) / stride;

    return voxelline;
}

/*
 * Convert target dimension/axis + lineno to VDS voxel coordinates.
 */
void set_voxels(
    axis ax,
    int dimension,
    int lineno,
    const OpenVDS::VolumeDataLayout *layout,
    int (&voxelmin)[OpenVDS::VolumeDataLayout::Dimensionality_Max],
    int (&voxelmax)[OpenVDS::VolumeDataLayout::Dimensionality_Max]
) {
    auto vmin = OpenVDS::IntVector3 { 0, 0, 0 };
    auto vmax = OpenVDS::IntVector3 {
        layout->GetDimensionNumSamples(0) - 1,
        layout->GetDimensionNumSamples(1) - 1,
        layout->GetDimensionNumSamples(2) - 1
    };

    int voxelline;
    auto vdim   = dim_tovoxel(dimension, layout);
    auto system = axis_tosystem(ax);
    switch (system) {
        case ANNOTATION: {
            voxelline = lineno_tovoxel(lineno, vdim, layout);
            break;
        }
        case INDEX: {
            /* Line-numbers in IJK match Voxel */
            voxelline = lineno;
            break;
        }
    }

    vmin[vdim] = voxelline;
    vmax[vdim] = voxelline + 1;

    /* Commit */
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
    auto dimension = axis_todim(ax);
    set_voxels(ax, dimension, lineno, layout, vmin, vmax);

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
