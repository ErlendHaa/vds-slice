#include "vds.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

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

struct vdsbuffer getData(std::string url, std::string credentials, int dimension, int lineno)
{
    OpenVDS::Error error;
    OpenVDS::VDSHandle handle = OpenVDS::Open(url, credentials, error);

    if(error.code != 0) {
        throw std::runtime_error("Could not open VDS: " + error.string);
    }

    auto accessManager = OpenVDS::GetAccessManager(handle);
    auto const *layout = accessManager.GetVolumeDataLayout();

    int voxelMin[OpenVDS::Dimensionality_Max] = { 0, 0, 0, 0, 0, 0};
    int voxelMax[OpenVDS::Dimensionality_Max] = { 1, 1, 1, 1, 1, 1};

    voxelMax[0] = layout->GetDimensionNumSamples(0);
    voxelMax[1] = layout->GetDimensionNumSamples(1);
    voxelMax[2] = layout->GetDimensionNumSamples(2);

    voxelMin[dimension] = lineno;
    voxelMax[dimension] = lineno + 1;

    auto request = accessManager.RequestVolumeSubset<float>(OpenVDS::Dimensions_012, 0, 0, voxelMin, voxelMax);
    auto slice = request->Data();

    /* There doesn't seam to be any straight forward way of taking ownership
     * over a std::vector's internal buffer, which means that any pointer to
     * that buffer is invalid after the vector goes out of scope. Hence we
     * copy the data to a new buffer to keep the returned pointer valid until
     * go has had a chance to parse it.
     */
    char *data = new char[slice.size()];
    std::copy(slice.begin(), slice.end(), data);

    vdsbuffer buf{};
    buf.size = slice.size();
    buf.data = data;

    return buf;
}

struct vdsbuffer fetch_slice(const char* vds,
                             const char* credentials,
                             int dim,
                             int lineno) {
    std::string cube(vds);
    std::string cred(credentials);

    try {
        return getData(cube, cred, dim, lineno);
    } catch (const std::exception& e) {
        vdsbuffer buf {};
        buf.err = new char[std::strlen(e.what()) + 1];
        std::strcpy(buf.err, e.what());
        return buf;
    }
}
