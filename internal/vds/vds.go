package vds

/*
#cgo LDFLAGS: -lopenvds
#cgo CXXFLAGS: -std=c++11
#include <vds.h>
#include <stdlib.h>
*/
import "C"
import "unsafe"

import (
	"errors"
)

func Slice(vds, credentials, direction string, lineno int) ([]byte, error) {
	cvds := C.CString(vds)
	defer C.free(unsafe.Pointer(cvds))

	ccred := C.CString(credentials)
	defer C.free(unsafe.Pointer(ccred))

	cdirection := C.CString(direction)
	defer C.free(unsafe.Pointer(cdirection))

	result := C.slice(
		cvds,
		ccred,
		cdirection,
		C.int(lineno),
	)

	defer C.vdsbuffer_delete(&result)

	if result.err != nil {
		err := C.GoString(result.err)
		return nil, errors.New(err)
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}
