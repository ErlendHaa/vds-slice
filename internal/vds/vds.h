#ifndef VDS_SLICE_CGO_VDS_H
#define VDS_SLICE_CGO_VDS_H

#ifdef __cplusplus
extern "C" {
#endif

struct vdsbuffer {
    char*         data;
    char*         err;
    unsigned long size;
};

enum axis {
    X      = 0,
    Y      = 1,
    Z      = 2,
    ILINE  = 3,
    XLINE  = 4,
    DEPTH  = 5,
    TIME   = 6,
    SAMPLE = 7,
};

struct vdsbuffer slice(
    const char* vds,
    const char* credentials,
    int lineno,
    enum axis direction
);

struct vdsbuffer slice_metadata(
    const char* vds,
    const char* credentials,
    int lineno,
    enum axis direction
);

void vdsbuffer_delete(struct vdsbuffer*);


#ifdef __cplusplus
}
#endif
#endif // VDS_SLICE_CGO_VDS_H
