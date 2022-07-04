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

struct vdsbuffer fetch_slice(const char* vds,
                             const char* credentials,
                             int dim,
                             int lineno);

void vdsbuffer_delete(struct vdsbuffer*);


#ifdef __cplusplus
}
#endif
#endif // VDS_SLICE_CGO_VDS_H
