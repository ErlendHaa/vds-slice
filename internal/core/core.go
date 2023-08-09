package core

/*
#cgo LDFLAGS: -lopenvds
#cgo CXXFLAGS: -std=c++17
#include <capi.h>
#include <ctypes.h>
#include <stdlib.h>
*/
import "C"
import (
	"errors"
	"fmt"
	"math"
	"strings"
	"unsafe"
)

const (
	AxisI         = C.I
	AxisJ         = C.J
	AxisK         = C.K
	AxisInline    = C.INLINE
	AxisCrossline = C.CROSSLINE
	AxisDepth     = C.DEPTH
	AxisTime      = C.TIME
	AxisSample    = C.SAMPLE
)

const (
	CoordinateSystemIndex      = C.INDEX
	CoordinateSystemAnnotation = C.ANNOTATION
	CoordinateSystemCdp        = C.CDP
)

// @Description Axis description
type Axis struct {
	// Name/Annotation of axis
	Annotation string `json:"annotation" example:"Sample"`

	// Minimum axis value
	Min float64 `json:"min" example:"4.0"`

	// Maximum axis value
	Max float64 `json:"max" example:"4500.0"`

	// Number of samples along the axis
	Samples int `json:"samples" example:"1600"`

	// Axis units
	Unit string `json:"unit" example:"ms"`
} // @name Axis

// @Description Geometrical plane with depth/time datapoints
type RegularSurface struct {
	// Values / height-map
	Values [][]float32 `json:"values" binding:"required"`

	// Rotation of the X-axis (East), counter-clockwise, in degrees
	Rotation *float32 `json:"rotation" binding:"required"`

	// X-coordinate of the origin
	Xori *float32 `json:"xori" binding:"required"`

	// Y-coordinate of the origin
	Yori *float32 `json:"yori" binding:"required"`

	// X-increment - The physical distance between height-map columns
	Xinc float32 `json:"xinc" binding:"required"`

	// Y-increment - The physical distance between height-map rows
	Yinc float32 `json:"yinc" binding:"required"`

	// Any sample in the input values with value == fillValue will be ignored
	// and the fillValue will be used in the amplitude map.
	// I.e. for any [i, j] where values[i][j] == fillValue then
	// output[i][j] == fillValue.
	// Additionally, the fillValue is used for any point of the surface that
	// falls outside the bounds of the seismic volume.
	FillValue *float32 `json:"fillValue" binding:"required"`
} // @name RegularSurface

// @Description The bounding box of the survey, defined by its 4 corner
// @Description coordinates. The bounding box is given in 3 different
// @Description coordinate systems. The points are sorted in the same order for
// @Description each system. E.g. [[1,2], [2,2], [2,3], [1,3]]
type BoundingBox struct {
	Cdp  [][]float64 `json:"cdp"`
	Ilxl [][]float64 `json:"ilxl"`
	Ij   [][]float64 `json:"ij"`
} //@name BoundingBox

type Array struct {
	// Data format is represented by numpy-style formatcodes. Currently the
	// format is always 4-byte floats, little endian (<f4).
	Format string `json:"format" example:"<f4"`

	// Shape of the returned data
	Shape []int `json:"shape" swaggertype:"array,integer" example:"10,50"`
}

type Bound struct {
	Direction string
	Lower     float32
	Upper     float32
}

// @Description Slice metadata
type SliceMetadata struct {
	Array

	// X-axis information
	X Axis `json:"x"`

	// Y-axis information
	Y Axis `json:"y"`

	/* Override shape for docs */

	// Shape of the returned slice. Equals to [Y.Samples, X.Samples]
	Shape []int `json:"shape" swaggertype:"array,integer" example:"10,50"`

	// Horizontal bounding box of the slice. For inline/crossline slices this
	// is a linestring, while for time/depth slices this is essentially the
	// bounding box of the volume.
	Geospatial [][]float64 `json:"geospatial"`
} // @name SliceMetadata

// @Description Metadata
type Metadata struct {
	// Coordinate reference system
	Crs string `json:"crs" example:"PROJCS[\"ED50 / UTM zone 31N\",..."`

	// The original input file name
	InputFileName string `json:"inputFileName" example:"file.segy"`

	// Import time stamp in ISO8601 format
	ImportTimeStamp string `json:"importTimeStamp" example:"2021-02-18T21:54:42.123Z"`

	// Bounding box
	BoundingBox BoundingBox `json:"boundingBox"`

	// Axis descriptions
	//
	// Describes the axes of the requested 2-dimensional slice.
	Axis []*Axis `json:"axis"`
} // @name Metadata

// @Description Fence metadata
type FenceMetadata struct {
	Array
} // @name FenceMetadata

// @Description Attribute metadata
type AttributeMetadata struct {
	Array
} // @name AttributeMetadata

func GetAxis(direction string) (int, error) {
	switch direction {
	case "i":
		return AxisI, nil
	case "j":
		return AxisJ, nil
	case "k":
		return AxisK, nil
	case "inline":
		return AxisInline, nil
	case "crossline":
		return AxisCrossline, nil
	case "depth":
		return AxisDepth, nil
	case "time":
		return AxisTime, nil
	case "sample":
		return AxisSample, nil
	default:
		options := "i, j, k, inline, crossline or depth/time/sample"
		msg := "invalid direction '%s', valid options are: %s"
		return -1, NewInvalidArgument(fmt.Sprintf(msg, direction, options))
	}
}

func GetCoordinateSystem(coordinateSystem string) (int, error) {
	switch strings.ToLower(coordinateSystem) {
	case "ij":
		return CoordinateSystemIndex, nil
	case "ilxl":
		return CoordinateSystemAnnotation, nil
	case "cdp":
		return CoordinateSystemCdp, nil
	default:
		options := "ij, ilxl, cdp"
		msg := "coordinate system not recognized: '%s', valid options are: %s"
		return -1, NewInvalidArgument(fmt.Sprintf(msg, coordinateSystem, options))
	}
}

func GetInterpolationMethod(interpolation string) (int, error) {
	switch strings.ToLower(interpolation) {
	case "":
		fallthrough
	case "nearest":
		return C.NEAREST, nil
	case "linear":
		return C.LINEAR, nil
	case "cubic":
		return C.CUBIC, nil
	case "angular":
		return C.ANGULAR, nil
	case "triangular":
		return C.TRIANGULAR, nil
	default:
		options := "nearest, linear, cubic, angular or triangular"
		msg := "invalid interpolation method '%s', valid options are: %s"
		return -1, NewInvalidArgument(fmt.Sprintf(msg, interpolation, options))
	}
}

func GetAttributeType(attribute string) (int, error) {
	switch strings.ToLower(attribute) {
	case "samplevalue":
		return C.VALUE, nil
	case "min":
		return C.MIN, nil
	case "max":
		return C.MAX, nil
	case "maxabs":
		return C.MAXABS, nil
	case "mean":
		return C.MEAN, nil
	case "meanabs":
		return C.MEANABS, nil
	case "meanpos":
		return C.MEANPOS, nil
	case "meanneg":
		return C.MEANNEG, nil
	case "median":
		return C.MEDIAN, nil
	case "rms":
		return C.RMS, nil
	case "var":
		return C.VAR, nil
	case "sd":
		return C.SD, nil
	case "sumpos":
		return C.SUMPOS, nil
	case "sumneg":
		return C.SUMNEG, nil
	case "":
		fallthrough
	default:
		options := []string{
			"samplevalue", "min", "max", "maxabs", "mean", "meanabs", "meanpos",
			"meanneg", "median", "rms", "var", "sd", "sumpos", "sumneg",
		}
		msg := "invalid attribute '%s', valid options are: %s"
		return -1, NewInvalidArgument(fmt.Sprintf(
			msg,
			attribute,
			strings.Join(options, ", "),
		))
	}
}

/** Translate C status codes into Go error types */
func toError(status C.int, ctx *C.Context) error {
	if status == C.STATUS_OK {
		return nil
	}

	msg := C.GoString(C.errmsg(ctx))

	switch status {
	case C.STATUS_NULLPTR_ERROR:
		fallthrough
	case C.STATUS_RUNTIME_ERROR:
		return NewInternalError(msg)
	case C.STATUS_BAD_REQUEST:
		return NewInvalidArgument(msg)
	default:
		return errors.New(msg)
	}
}

type CRegularSurface struct {
	cSurface *C.struct_RegularSurface
	cData    []C.float
}

func (r *CRegularSurface) get() *C.struct_RegularSurface {
	return r.cSurface
}

func (r *CRegularSurface) Close() error {
	var cCtx = C.context_new()
	defer C.context_free(cCtx)

	cErr := C.regular_surface_free(cCtx, r.cSurface)
	if err := toError(cErr, cCtx); err != nil {
		return err
	}

	return nil
}

func (surface *RegularSurface) ToCRegularSurface() (CRegularSurface, error) {
	var cCtx = C.context_new()
	defer C.context_free(cCtx)

	nrows := len(surface.Values)
	ncols := len(surface.Values[0])

	cdata := make([]C.float, nrows*ncols)
	for i := range surface.Values {
		if len(surface.Values[i]) != ncols {
			msg := fmt.Sprintf(
				"Surface rows are not of the same length. "+
					"Row 0 has %d elements. Row %d has %d elements",
				ncols, i, len(surface.Values[i]),
			)
			return CRegularSurface{}, NewInvalidArgument(msg)
		}

		for j := range surface.Values[i] {
			cdata[i*ncols+j] = C.float(surface.Values[i][j])
		}
	}

	var cSurface *C.struct_RegularSurface
	cErr := C.regular_surface_new(
		cCtx,
		&cdata[0],
		C.size_t(nrows),
		C.size_t(ncols),
		C.float(*surface.Xori),
		C.float(*surface.Yori),
		C.float(surface.Xinc),
		C.float(surface.Yinc),
		C.float(*surface.Rotation),
		C.float(*surface.FillValue),
		&cSurface,
	)

	if err := toError(cErr, cCtx); err != nil {
		C.regular_surface_free(cCtx, cSurface)
		return CRegularSurface{}, err
	}

	return CRegularSurface{cSurface: cSurface, cData: cdata}, nil
}

type VDSHandle struct {
	handle *C.struct_DataHandle
	ctx    *C.struct_Context
}

func (v VDSHandle) Handle() *C.struct_DataHandle {
	return v.handle
}

func (v VDSHandle) context() *C.struct_Context {
	return v.ctx
}

func (v VDSHandle) Error(status C.int) error {
	return toError(status, v.context())
}

func (v VDSHandle) Close() error {
	defer C.context_free(v.ctx)

	cerr := C.datahandle_free(v.ctx, v.handle)
	return toError(cerr, v.ctx)
}

func NewVDSHandle(conn Connection) (VDSHandle, error) {
	curl := C.CString(conn.Url())
	defer C.free(unsafe.Pointer(curl))

	ccred := C.CString(conn.ConnectionString())
	defer C.free(unsafe.Pointer(ccred))

	var cctx = C.context_new()
	var handle *C.struct_DataHandle

	cerr := C.datahandle_new(cctx, curl, ccred, &handle)

	if err := toError(cerr, cctx); err != nil {
		defer C.context_free(cctx)
		return VDSHandle{}, err
	}

	return VDSHandle{handle: handle, ctx: cctx}, nil
}

func (v VDSHandle) GetMetadata() ([]byte, error) {
	var result C.struct_response
	cerr := C.metadata(v.context(), v.Handle(), &result)

	defer C.response_delete(&result)

	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func newCSliceBounds(bounds []Bound) ([]C.struct_Bound, error) {
	var cBounds []C.struct_Bound
	for _, bound := range bounds {
		direction, err := GetAxis(bound.Direction)
		if err != nil {
			return nil, err
		}

		cBound := C.struct_Bound{
			C.float(bound.Lower),
			C.float(bound.Upper),
			C.enum_axis_name(direction),
		}
		cBounds = append(cBounds, cBound)
	}

	return cBounds, nil
}

func (v VDSHandle) GetSlice(lineno, direction int, bounds []Bound) ([]byte, error) {
	var result C.struct_response

	cBounds, err := newCSliceBounds(bounds)
	if err != nil {
		return nil, err
	}

	var bound *C.struct_Bound
	if len(cBounds) > 0 {
		bound = &cBounds[0]
	}

	cerr := C.slice(
		v.context(),
		v.Handle(),
		C.int(lineno),
		C.enum_axis_name(direction),
		bound,
		C.size_t(len(cBounds)),
		&result,
	)

	defer C.response_delete(&result)
	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func (v VDSHandle) GetSliceMetadata(
	lineno    int,
	direction int,
	bounds    []Bound,
) ([]byte, error) {
	var result C.struct_response

	cBounds, err := newCSliceBounds(bounds)
	if err != nil {
		return nil, err
	}

	var bound *C.struct_Bound
	if len(cBounds) > 0 {
		bound = &cBounds[0]
	}

	cerr := C.slice_metadata(
		v.context(),
		v.Handle(),
		C.int(lineno),
		C.enum_axis_name(direction),
		bound,
		C.size_t(len(cBounds)),
		&result,
	)

	defer C.response_delete(&result)

	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func (v VDSHandle) GetFence(
	coordinateSystem int,
	coordinates [][]float32,
	interpolation int,
) ([]byte, error) {
	coordinate_len := 2
	ccoordinates := make([]C.float, len(coordinates)*coordinate_len)
	for i := range coordinates {

		if len(coordinates[i]) != coordinate_len {
			msg := fmt.Sprintf(
				"invalid coordinate %v at position %d, expected [x y] pair",
				coordinates[i],
				i,
			)
			return nil, NewInvalidArgument(msg)
		}

		for j := range coordinates[i] {
			ccoordinates[i*coordinate_len+j] = C.float(coordinates[i][j])
		}
	}

	var result C.struct_response
	cerr := C.fence(
		v.context(),
		v.Handle(),
		C.enum_coordinate_system(coordinateSystem),
		&ccoordinates[0],
		C.size_t(len(coordinates)),
		C.enum_interpolation_method(interpolation),
		&result,
	)

	defer C.response_delete(&result)

	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func (v VDSHandle) GetFenceMetadata(coordinates [][]float32) ([]byte, error) {
	var result C.struct_response
	cerr := C.fence_metadata(
		v.context(),
		v.Handle(),
		C.size_t(len(coordinates)),
		&result,
	)

	defer C.response_delete(&result)

	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func (v VDSHandle) GetAttributeMetadata(data [][]float32) ([]byte, error) {
	var result C.struct_response
	cerr := C.attribute_metadata(
		v.context(),
		v.Handle(),
		C.size_t(len(data)),
		C.size_t(len(data[0])),
		&result,
	)

	defer C.response_delete(&result)

	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	buf := C.GoBytes(unsafe.Pointer(result.data), C.int(result.size))
	return buf, nil
}

func (v VDSHandle) GetAttributes(
	surface RegularSurface,
	above float32,
	below float32,
	stepsize float32,
	attributes []string,
	interpolation int,
) ([][]byte, error) {
	var targetAttributes []int
	for _, attr := range attributes {
		id, err := GetAttributeType(attr)
		if err != nil {
			return nil, err
		}
		targetAttributes = append(targetAttributes, id)
	}

	var nrows = len(surface.Values)
	var ncols = len(surface.Values[0])
	var hsize = nrows * ncols

	cAttributes := make([]C.enum_attribute, len(targetAttributes))
	for i := range targetAttributes {
		cAttributes[i] = C.enum_attribute(targetAttributes[i])
	}

	cSurface, err := surface.ToCRegularSurface()
	if err != nil {
		return nil, err
	}
	defer cSurface.Close()

	var horizonSize C.size_t

	cerr := C.horizon_size(
		v.context(),
		v.Handle(),
		cSurface.get(),
		C.float(above),
		C.float(below),
		&horizonSize,
	)
	if err := v.Error(cerr); err != nil {
		return nil, err
	}

	horizon, err := v.fetchHorizon(
		cSurface,
		nrows,
		ncols,
		above,
		below,
		interpolation,
		horizonSize,
	)
	if err != nil {
		return nil, err
	}

	return v.calculateAttributes(
		cSurface,
		hsize,
		horizon,
		horizonSize,
		cAttributes,
		above,
		below,
		stepsize,
	)
}

func (v VDSHandle) fetchHorizon(
	cSurface CRegularSurface,
	nrows int,
	ncols int,
	above float32,
	below float32,
	interpolation int,
	horizonSize C.size_t,
) ([]byte, error) {
	buffer := make([]byte, horizonSize)
	hsize := nrows * ncols
	// note that it is possible to hit go's own goroutines limit
	// but we do not deal with it here

	// max number of goroutines running at the same time
	// too low number doesn't utilize all CPU, too high overuses it
	// value should be experimented with
	maxConcurrentGoroutines := max(nrows/2, 1)
	guard := make(chan struct{}, maxConcurrentGoroutines)

	// the size of the data processed in one goroutine
	// decides how many parts data is split into
	// value should be experimented with
	chunkSize := max(nrows, 1)

	from := 0
	to := from + chunkSize

	errs := make(chan error, hsize/chunkSize+1)
	nRoutines := 0

	for from < hsize {
		guard <- struct{}{} // block if guard channel is filled
		go func(from, to int) {
			var cCtx = C.context_new()
			defer C.context_free(cCtx)

			cerr := C.horizon(
				cCtx,
				v.Handle(),
				cSurface.get(),
				C.float(above),
				C.float(below),
				C.enum_interpolation_method(interpolation),
				C.size_t(from),
				C.size_t(to),
				unsafe.Pointer(&buffer[0]),
			)
			errs <- toError(cerr, cCtx)
			<-guard
		}(from, to)

		nRoutines += 1

		from += chunkSize
		to = min(to+chunkSize, hsize)
	}

	// Wait for all gorutines to finish and collect any errors
	var computeErrors []error
	for i := 0; i < nRoutines; i++ {
		err := <-errs
		if err != nil {
			computeErrors = append(computeErrors, err)
		}
	}

	if len(computeErrors) > 0 {
		return nil, computeErrors[0]
	}
	return buffer, nil
}

func (v VDSHandle) calculateAttributes(
	cSurface CRegularSurface,
	hsize int,
	horizon []byte,
	horizonSize C.size_t,
	cAttributes []uint32,
	above float32,
	below float32,
	stepsize float32,
) ([][]byte, error) {
	nAttributes := len(cAttributes)
	var mapsize = hsize * 4
	buffer := make([]byte, mapsize*nAttributes)

	maxConcurrency := 32
	windowsPerRoutine := int(math.Ceil(float64(hsize) / float64(maxConcurrency)))

	errs := make(chan error, maxConcurrency)

	from := 0
	remaining := hsize
	nRoutines := 0
	for remaining > 0 {
		nRoutines++

		size := min(windowsPerRoutine, remaining)
		to := from + size

		go func(from, to int) {
			var cCtx = C.context_new()
			defer C.context_free(cCtx)

			cErr := C.attribute(
				cCtx,
				v.Handle(),
				cSurface.get(),
				unsafe.Pointer(&horizon[0]),
				horizonSize,
				&cAttributes[0],
				C.size_t(nAttributes),
				C.float(above),
				C.float(below),
				C.float(stepsize),
				C.size_t(from),
				C.size_t(to),
				unsafe.Pointer(&buffer[0]),
			)

			errs <- toError(cErr, cCtx)
		}(from, to)

		from = to
		remaining -= size
	}

	// Wait for all gorutines to finish and collect any errors
	var computeErrors []error
	for i := 0; i < nRoutines; i++ {
		err := <-errs
		if err != nil {
			computeErrors = append(computeErrors, err)
		}
	}

	if len(computeErrors) > 0 {
		return nil, computeErrors[0]
	}

	out := make([][]byte, nAttributes)
	for i := 0; i < nAttributes; i++ {
		out[i] = buffer[i*mapsize : (i+1)*mapsize]
	}

	return out, nil
}
