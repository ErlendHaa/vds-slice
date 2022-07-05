#include "vds.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdexcept>

#include <OpenVDS/IJKCoordinateTransformer.h>
#include <OpenVDS/OpenVDS.h>

using namespace std;

void vdsbuffer_delete(struct vdsbuffer* buf) {
    if (!buf)
        return;

    delete[] buf->data;
    delete[] buf->err;
    *buf = vdsbuffer {};
}


int todimension(std::string direction, const OpenVDS::VolumeDataLayout *layout) {
    for (int i = 0; i < 3; i++) {
        std::string name(layout->GetDimensionName(i));

        if (direction.compare(name) == 0) {
            return i;
        }
    }

    throw std::runtime_error(
        "Unable to convert direction '" +
        direction +
        "' to VDS dimension"
    );
};

struct vdsbuffer fetch_slice(
    std::string url,
    std::string credentials,
    std::string direction,
    int lineno
) {
    OpenVDS::Error error;
    OpenVDS::VDSHandle handle = OpenVDS::Open(url, credentials, error);

    if(error.code != 0) {
        throw std::runtime_error("Could not open VDS: " + error.string);
    }

    auto access = OpenVDS::GetAccessManager(handle);
    auto const *layout = access.GetVolumeDataLayout();

    auto dimension = todimension(direction, layout);

    int voxelMin[OpenVDS::Dimensionality_Max] = { 0, 0, 0, 0, 0, 0};
    int voxelMax[OpenVDS::Dimensionality_Max] = { 1, 1, 1, 1, 1, 1};

    voxelMax[0] = layout->GetDimensionNumSamples(0);
    voxelMax[1] = layout->GetDimensionNumSamples(1);
    voxelMax[2] = layout->GetDimensionNumSamples(2);

    voxelMin[dimension] = lineno;
    voxelMax[dimension] = lineno + 1;

    auto format = layout->GetChannelFormat(0);
    auto size = access.GetVolumeSubsetBufferSize(
        voxelMin,
        voxelMax,
        format,
        0,
        0
    );

    auto *data = new char[size];
    auto request = access.RequestVolumeSubset(
        data,
        size,
        OpenVDS::Dimensions_012,
        0,
        0,
        voxelMin,
        voxelMax,
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
    const char* direction,
    int lineno
) {
    std::string cube(vds);
    std::string cred(credentials);
    std::string dir(direction);

    try {
        return fetch_slice(cube, cred, dir, lineno);
    } catch (const std::exception& e) {
        vdsbuffer buf {};
        buf.err = new char[std::strlen(e.what()) + 1];
        std::strcpy(buf.err, e.what());
        return buf;
    }
}
