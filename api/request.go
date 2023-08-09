package api

import (
	"encoding/json"
	"fmt"
	"net/url"
	"strings"

	"github.com/equinor/vds-slice/internal/cache"
	"github.com/equinor/vds-slice/internal/core"
)

type RequestedResource struct {
	// The blob URL can be sent in either format:
	// - https://account.blob.core.windows.net/container/blob
	//    In the above case the user must provide a sas-token as a separate key.
	//
	// - https://account.blob.core.windows.net/container/blob?sp=r&st=2022-09-12T09:44:17Z&se=2022-09-12T17:44:17Z&spr=https&sv=2021-06-08&sr=c&sig=...
	//	  Instead of passing the sas-token explicitly in the sas field, you can
	//	  pass an sign url. If the sas-token is provided in both fields, the
	//	  sas-token in the sas field is prioritized.
	//
	// Note that your whole query string will be passed further down to
	// openvds. We expect query parameters to contain sas-token and sas-token
	// only and give no guarantee that Openvds/Azure returns you if you provide
	// any additional arguments.
	//
	// Warning: We do not expect storage accounts to have snapshots. If your
	// storage account has any, please contact System Admin, as due to caching
	// you might end up with incorrect data.
	Vds string `json:"vds" binding:"required" example:"https://account.blob.core.windows.net/container/blob"`

	// A valid sas-token with read access to the container specified in Vds
	Sas string `json:"sas,omitempty" example:"sp=r&st=2022-09-12T09:44:17Z&se=2022-09-12T17:44:17Z&spr=https&sv=2021-06-08&sr=c&sig=..."`
}

type Stringable interface {
	toString() (string, error)
}
type Normalizable interface {
	NormalizeConnection() error
}

func (r *RequestedResource) NormalizeConnection() error {
	url, err := url.Parse(r.Vds)
	if err != nil {
		return core.NewInvalidArgument(err.Error())
	}
	if strings.TrimSpace(r.Sas) == "" {
		if url.RawQuery == "" {
			return core.NewInvalidArgument("No valid Sas token is found in the request")
		}
		r.Sas = url.RawQuery
	}

	url.RawQuery = ""
	url.Host = url.Hostname()
	r.Vds = url.String()
	return nil
}

type MetadataRequest struct {
	RequestedResource
} //@name MetadataRequest

func (m MetadataRequest) toString() (string, error) {
	m.Sas = ""
	out, err := json.Marshal(m)
	if err != nil {
		return "", err
	}
	str := string(out)
	return str, nil
}

type FenceRequest struct {
	RequestedResource
	// Coordinate system for the requested fence
	// Supported options are:
	// ilxl : inline, crossline pairs
	// ij   : Coordinates are given as in 0-indexed system, where the first
	//        line in each direction is 0 and the last is number-of-lines - 1.
	// cdp  : Coordinates are given as cdpx/cdpy pairs. In the original SEGY
	//        this would correspond to the cdpx and cdpy fields in the
	//        trace-headers after applying the scaling factor.
	CoordinateSystem string `json:"coordinateSystem" binding:"required" example:"cdp"`

	// A list of (x, y) points in the coordinate system specified in
	// coordinateSystem, for example [[2000.5, 100.5], [2050, 200], [10, 20]].
	// All coordinates are expected to be inside the bounding box (with outer margin
	// of half a distance between consecutive lines).
	Coordinates [][]float32 `json:"coordinates" binding:"required"`

	// Interpolation method
	// Supported options are: nearest, linear, cubic, angular and triangular.
	// Defaults to nearest.
	// This field is passed on to OpenVDS, which does the actual interpolation.
	// Note: For nearest interpolation result will snap to the nearest point
	// as per "half up" rounding. This is different from openvds logic.
	Interpolation string `json:"interpolation" example:"linear"`
} //@name FenceRequest

func (f FenceRequest) toString() (string, error) {
	coordinates := func() string {
		var length = len(f.Coordinates)
		const halfPrintLength = 5
		const printLength = halfPrintLength * 2
		if length > printLength {
			return fmt.Sprintf("%v, ...[%d element(s) skipped]..., %v",
				f.Coordinates[0:halfPrintLength],
				length-printLength,
				f.Coordinates[length-halfPrintLength:length])
		} else {
			return fmt.Sprintf("%v", f.Coordinates)
		}
	}()

	return fmt.Sprintf("{vds: %s, coordinate system: %s, coordinates: %s, interpolation (optional): %s}",
		f.Vds,
		f.CoordinateSystem,
		coordinates,
		f.Interpolation,
	), nil
}

/** Compute a hash of the request that uniquely identifies the requested fence
 *
 * The hash is computed based on all fields that contribute toward a unique response.
 * I.e. every field except the sas token.
 */
func (f FenceRequest) Hash() (string, error) {
	// Strip the sas token before computing hash
	f.Sas = ""
	return cache.Hash(f)
}

// Query for slice endpoints
// @Description Query payload for slice endpoint /slice.
type SliceRequest struct {
	RequestedResource

	// Direction can be specified in two domains
	// - Annotation. Valid options: Inline, Crossline and Depth/Time/Sample
	// - Index. Valid options: i, j and k.
	//
	// When requesting z-slices using annotation it's recommended to use Depth
	// or Time as they include validation against the VDS itself. I.e. the
	// server returns an error if you request a time-slice from a depth cube
	// and visa-versa. Use Sample if you don't care for such guarantees or as a
	// fallback if the VDS file is wonky.
	//
	// i, j, k are zero-indexed and correspond to Inline, Crossline,
	// Depth/Time/Sample, respectively.
	//
	// All options are case-insensitive.
	Direction string `json:"direction" binding:"required" example:"inline"`

	// Line number of the slice
	Lineno *int `json:"lineno" binding:"required" example:"10000"`

	// Restrict the slice in the other dimensions (sub-slicing)
	//
	// Bounds can be used to retrieve sub-slices. For example: when requesting
	// a time slice you can set inline and/or crossline bounds to restrict the
	// area of the slice.
	//
	// Bounds are optional. Not specifying any bounds will produce the requested
	// slice through the entire volume.
	//
	// Any bounds in the same direction as the slice itself are ignored.
	//
	// Bounds are applied one after the other. If there are multiple bounds in the
	// same direction, the last one takes precedence.
	//
	// Bounds will throw out-of-bounds error if their range is outside the
	// cubes dimensions.
	//
	// Bounds can be set using both annotation and index. You are free to mix
	// and match as you see fit.
	Bounds []core.Bound `json:"bounds"`
} //@name SliceRequest

/** Compute a hash of the request that uniquely identifies the requested slice
 *
 * The hash is computed based on all fields that contribute toward a unique response.
 * I.e. every field except the sas token.
 */
func (s SliceRequest) Hash() (string, error) {
	// Strip the sas token before computing hash
	s.Sas = ""
	return cache.Hash(s)
}

func (s SliceRequest) toString() (string, error) {
	s.Sas = ""
	out, err := json.Marshal(s)
	if err != nil {
		return "", err
	}
	str := string(out)
	return str, nil
}

// Query for Attribute endpoints
// @Description Query payload for attribute endpoint.
type AttributeRequest struct {
	RequestedResource

	// Surface along which data must be retrieved
	Surface core.RegularSurface `json:"surface" binding:"required"`

	// Horizontal interpolation method
	// Supported options are: nearest, linear, cubic, angular and triangular.
	// Defaults to nearest.
	// This field is passed on to OpenVDS, which does the actual interpolation.
	//
	// This only applies to the horizontal plane. Traces are always
	// interpolated with cubic interpolation (algorithm: modified makima)
	// Note: For nearest interpolation result will snap to the nearest point
	// as per "half up" rounding. This is different from openvds logic.
	Interpolation string `json:"interpolation" example:"linear"`

	// Samples interval above the horizon to include in attribute calculation.
	// This value should be given in the VDS's vertical domain. E.g. if the
	// vertical domain is 'ms' then a value of 22 means to include samples up
	// to 22 ms above the horizon definition. The value is rounded down to the
	// nearest whole sample. I.e. if the cube is sampled at 4ms, the attribute
	// calculation will include samples 4, 8, 12, 16 and 20ms above the
	// horizon, while the sample at 24ms is excluded.
	//
	// Defaults to zero
	Above float32 `json:"above"`

	// Samples interval below the horizon to include in attribute calculation.
	// Implements the same behaviour as 'above'.
	//
	// Defaults to zero
	Below float32 `json:"below"`

	// Stepsize for samples within the window defined by above below
	//
	// Samples within the vertical window will be re-sampled to 'stepsize'
	// using cubic interpolation (modified makima) before the attributes are
	// calculated.
	//
	// This value should be given in the vertical domain of the traces. E.g.
	// 0.1 implies re-sample samples at an interval of 0.1 meter (if it's a
	// depth cube) or 0.1 ms (if it's a time cube with vertical units in
	// milliseconds).
	//
	// Setting this to zero, or omitting it will default it to the vertical
	// stepsize in the VDS volume.
	Stepsize float32 `json:"stepsize"`

	// Requested attributes. Multiple attributes can be calculated by the same
	// request. This is considerably faster than doing one request per
	// attribute.
	Attributes []string `json:"attributes" binding:"required"`
} //@name AttributeRequest

/** Compute a hash of the request that uniquely identifies the requested attributes
 *
 * The hash is computed based on all fields that contribute toward a unique response.
 * I.e. every field except the sas token.
 */
func (h AttributeRequest) Hash() (string, error) {
	// Strip the sas token before computing hash
	h.Sas = ""
	return cache.Hash(h)
}

func (h AttributeRequest) toString() (string, error) {
	msg := "{vds: %s, Horizon: (ncols: %d, nrows: %d), Rotation: %.2f, " +
		"Origin: [%.2f, %.2f], Increment: [%.2f, %.2f], FillValue: %.2f, " +
		"interpolation: %s, Above: %.2f, Below: %.2f, Stepsize: %.2f, " +
		"Attributes: %v}"
	return fmt.Sprintf(
		msg,
		h.Vds,
		len(h.Surface.Values[0]),
		len(h.Surface.Values),
		*h.Surface.Rotation,
		*h.Surface.Xori,
		*h.Surface.Yori,
		h.Surface.Xinc,
		h.Surface.Yinc,
		*h.Surface.FillValue,
		h.Interpolation,
		h.Above,
		h.Below,
		h.Stepsize,
		h.Attributes,
	), nil
}
