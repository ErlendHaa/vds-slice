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

const (
	AxisX      = C.X
	AxisY      = C.Y
	AxisZ      = C.Z
	AxisIline  = C.ILINE
	AxisXline  = C.XLINE
	AxisDepth  = C.DEPTH
	AxisTime   = C.TIME
	AxisSample = C.SAMPLE
)

func Slice(vds, credentials string, lineno, direction int) ([]byte, error) {
	cvds := C.CString(vds)
	defer C.free(unsafe.Pointer(cvds))

	ccred := C.CString(credentials)
	defer C.free(unsafe.Pointer(ccred))

	result := C.slice(
		cvds,
		ccred,
		C.int(lineno),
		C.enum_axis(direction),
	)

	defer C.vdsbuffer_delete(&result)

	if result.err != nil {
		err := C.GoString(result.err)
		return nil, errors.New(err)
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}
