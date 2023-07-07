package api

import (
	"encoding/json"
	"fmt"
	"math"
	"net/http"
	"strings"

	"github.com/gin-gonic/gin"
	"github.com/gin-gonic/gin/binding"

	"github.com/equinor/vds-slice/internal/cache"
	"github.com/equinor/vds-slice/internal/vds"
)

func httpStatusCode(err error) int {
	switch err.(type){
	case *vds.InvalidArgument: return http.StatusBadRequest
	case *vds.InternalError:   return http.StatusInternalServerError
	default:
		return http.StatusInternalServerError
	}
}

/* Call abortOnError on the context in case of an error
 *
 * This function is designed specifically for our endpoint handler functions
 * and aims at making the errorhandling as short and concise as possible.
 *
 * If err != nil the error will be mapped to an appropriate http status code
 * through the httpStatusCode mapper, and ctx.AbortWithError will be called
 * with this status and the error itself. It then returns true to indicate that
 * the context have been aborted.
 *
 * If err == nil the ctx is left untouched and this function returns false,
 * indicating that the context was not aborted.
 *
 * The result is a oneline error handling:
 *
 *     err, _ := func()
 *     if abortOnError(ctx, err) { return }
 */
func abortOnError(ctx *gin.Context, err error) bool {
	if err == nil {
		return false
	}

	ctx.AbortWithError(httpStatusCode(err), err)

	return true
}

type Endpoint struct {
	MakeVdsConnection vds.ConnectionMaker
	Cache             cache.Cache
	Goroutines        int
}

func prepareRequestLogging(ctx *gin.Context, request Stringable) {
	// ignore possible errors as they should not change outcome for the user
	requestString, _ := request.toString()
	ctx.Set("request", requestString)
}

func (e *Endpoint) metadata(ctx *gin.Context, request MetadataRequest) {
	prepareRequestLogging(ctx, request)
	conn, err := e.MakeVdsConnection(request.Vds, request.Sas)
	if abortOnError(ctx, err) { return }

	handle, err := vds.NewVDSHandle(conn)
	if abortOnError(ctx, err) { return }
	defer handle.Close()

	buffer, err := handle.GetMetadata()
	if abortOnError(ctx, err) { return }

	ctx.Data(http.StatusOK, "application/json", buffer)
}

func (e *Endpoint) slice(ctx *gin.Context, request SliceRequest) {
	prepareRequestLogging(ctx, request)
	conn, err := e.MakeVdsConnection(request.Vds, request.Sas)
	if abortOnError(ctx, err) { return }

	cacheKey, err := request.Hash()
	if abortOnError(ctx, err) { return }

	cacheEntry, hit := e.Cache.Get(cacheKey)
	if hit && conn.IsAuthorizedToRead() {
		ctx.Set("cache-hit", true)
		writeResponse(ctx, cacheEntry.Metadata(), cacheEntry.Data())
		return
	}

	handle, err := vds.NewVDSHandle(conn)
	if abortOnError(ctx, err) { return }
	defer handle.Close()

	axis, err := vds.GetAxis(strings.ToLower(request.Direction))
	if abortOnError(ctx, err) { return }

	metadata, err := handle.GetSliceMetadata(*request.Lineno, axis)
	if abortOnError(ctx, err) { return }

	data, err := handle.GetSlice(*request.Lineno, axis)
	if abortOnError(ctx, err) { return }

	e.Cache.Set(cacheKey, cache.NewCacheEntry([][]byte{data}, metadata));

	writeResponse(ctx, metadata, [][]byte{data})
}

type FencePart struct {
	data []byte
	part int
}

func min(a, b int) int {
	if a < b {
		return a
	}

	return b
}

func (e *Endpoint) fence(ctx *gin.Context, request FenceRequest) {
	prepareRequestLogging(ctx, request)
	conn, err := e.MakeVdsConnection(request.Vds, request.Sas)
	if abortOnError(ctx, err) { return }

	cacheKey, err := request.Hash()
	if abortOnError(ctx, err) { return }

	cacheEntry, hit := e.Cache.Get(cacheKey)
	if hit && conn.IsAuthorizedToRead() {
		ctx.Set("cache-hit", true)
		writeResponse(ctx, cacheEntry.Metadata(), cacheEntry.Data())
		return
	}

	coordinateSystem, err := vds.GetCoordinateSystem(
		strings.ToLower(request.CoordinateSystem),
	)
	if abortOnError(ctx, err) { return }

	interpolation, err := vds.GetInterpolationMethod(request.Interpolation)
	if abortOnError(ctx, err) { return }
	
	handle, err := vds.NewVDSHandle(conn)
	if abortOnError(ctx, err) { return }
	defer handle.Close()

	shape, err := handle.Shape()
	if abortOnError(ctx, err) { return }

	nParts := e.Goroutines
	responses := make(chan FencePart, nParts)
	errors    := make(chan error,     nParts) 
	
	nTraces := len(request.Coordinates)
	nTracesPerPart := int(math.Ceil(float64(nTraces) / float64(nParts)))
	partSize := nTracesPerPart * shape.Samples * 4

	start := 0
	partN := 0
	remaining := nTraces
	for remaining > 0 { 
		size := min(nTracesPerPart, remaining)
		end := start + size

		go func(start, end, part int) {
			// TODO: The handle itself should be threadsafe (assuming
			// OpenVDS::Handle is threadsafe, which I'm not completly sure
			// about btw). But the context object it uses for error messageging
			// is not threadsafe. That should be fairly easy to fix, but I've
			// ignored it for now.
			data, err := handle.GetFence(
				coordinateSystem,
				request.Coordinates[start : end],
				interpolation,
			)
			if err != nil {
				errors <- err
				return
			}

			responses <- FencePart{ data: data, part: part }
		}(start, end, partN)

		start = end
		remaining -= size
		partN++
	}

	data := make([]byte, nTraces * shape.Samples * 4)
	for i := 0; i < partN; i++ {
		select {
		case err := <-errors:
			if abortOnError(ctx, err) { return }
		case response := <-responses: 
			start := response.part * partSize
			end   := start + len(response.data)
			copy(data[start : end], response.data)
		}
	}

	metadata, err := handle.GetFenceMetadata(request.Coordinates)
	if abortOnError(ctx, err) { return }

	e.Cache.Set(cacheKey, cache.NewCacheEntry([][]byte{data}, metadata));

	writeResponse(ctx, metadata, [][]byte{data})
}

func validateVerticalWindow(above float32, below float32, stepSize float32) error {
	const lowerBound = 0
	const upperBound = 250

	if above < lowerBound || above >= upperBound {
		return vds.NewInvalidArgument(fmt.Sprintf(
			"'above' out of range! Must be within [%d, %d], was %f",
			lowerBound,
			upperBound,
			above,
		))
	}
	
	if below < lowerBound || below >= upperBound {
		return vds.NewInvalidArgument(fmt.Sprintf(
			"'below' out of range! Must be within [%d, %d], was %f",
			lowerBound,
			upperBound,
			below,
		))
	}
	
	if stepSize < lowerBound {
		return vds.NewInvalidArgument(fmt.Sprintf(
			"'stepsize' out of range! Must be bigger than %d, was %f",
			lowerBound,
			stepSize,
		))
	}

	return nil
}

func (e *Endpoint) attributes(ctx *gin.Context, request AttributeRequest) {
	prepareRequestLogging(ctx, request)

	err := validateVerticalWindow(request.Above, request.Below, request.Stepsize)
	if abortOnError(ctx, err) { return }

	conn, err := e.MakeVdsConnection(request.Vds, request.Sas)
	if abortOnError(ctx, err) { return }

	cacheKey, err := request.Hash()
	if abortOnError(ctx, err) { return }

	handle, err := vds.NewVDSHandle(conn)
	if abortOnError(ctx, err) { return }
	defer handle.Close()

	cacheEntry, hit := e.Cache.Get(cacheKey)
	if hit && conn.IsAuthorizedToRead() {
		ctx.Set("cache-hit", true)
		writeResponse(ctx, cacheEntry.Metadata(), cacheEntry.Data())
		return
	}

	interpolation, err := vds.GetInterpolationMethod(request.Interpolation)
	if abortOnError(ctx, err) { return }

	metadata, err := handle.GetAttributeMetadata(request.Horizon)
	if abortOnError(ctx, err) { return }

	data, err := handle.GetAttributes(
		request.Horizon,
		*request.Xori,
		*request.Yori,
		request.Xinc,
		request.Yinc,
		*request.Rotation,
		*request.FillValue,
		request.Above,
		request.Below,
		request.Stepsize,
		request.Attributes,
		interpolation,
	)
	if abortOnError(ctx, err) { return }

	e.Cache.Set(cacheKey, cache.NewCacheEntry(data, metadata));

	writeResponse(ctx, metadata, data)
}

func parseGetRequest(ctx *gin.Context, v Normalizable) error {
	if err := json.Unmarshal([]byte(ctx.Query("query")), v); err != nil {
		return vds.NewInvalidArgument(err.Error())
	}

	if err := binding.Validator.ValidateStruct(v); err != nil {
		return vds.NewInvalidArgument(err.Error())
	}

	return v.NormalizeConnection()
}

func parsePostRequest(ctx *gin.Context, v Normalizable) error {
	if err := ctx.ShouldBind(v); err != nil {
		return vds.NewInvalidArgument(err.Error())
	}
	return v.NormalizeConnection()
}

func (e *Endpoint) Health(ctx *gin.Context) {
	ctx.String(http.StatusOK, "I am up and running")
}

// MetadataGet godoc
// @Summary  Return volumetric metadata about the VDS
// @Tags     metadata
// @Param    query  query  string  True  "Urlencoded/escaped MetadataRequest"
// @Produce  json
// @Success  200 {object} vds.Metadata
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /metadata  [get]
func (e *Endpoint) MetadataGet(ctx *gin.Context) {
	var request MetadataRequest
	err := parseGetRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.metadata(ctx, request)
}

// MetadataPost godoc
// @Summary  Return volumetric metadata about the VDS
// @Tags     metadata
// @Param    body  body  MetadataRequest  True  "Request parameters"
// @Produce  json
// @Success  200 {object} vds.Metadata
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /metadata  [post]
func (e *Endpoint) MetadataPost(ctx *gin.Context) {
	var request MetadataRequest
	err := parsePostRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.metadata(ctx, request)
}

// SliceGet godoc
// @Summary  Fetch a slice from a VDS
// @description.markdown slice
// @Tags     slice
// @Param    query  query  string  True  "Urlencoded/escaped SliceRequest"
// @Produce  multipart/mixed
// @Success  200 {object} vds.SliceMetadata "(Example below only for metadata part)"
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /slice  [get]
func (e *Endpoint) SliceGet(ctx *gin.Context) {
	var request SliceRequest
	err := parseGetRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.slice(ctx, request)
}

// SlicePost godoc
// @Summary  Fetch a slice from a VDS
// @description.markdown slice
// @Tags     slice
// @Param    body  body  SliceRequest  True  "Query Parameters"
// @Accept   application/json
// @Produce  multipart/mixed
// @Success  200 {object} vds.SliceMetadata "(Example below only for metadata part)"
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /slice  [post]
func (e *Endpoint) SlicePost(ctx *gin.Context) {
	var request SliceRequest
	err := parsePostRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.slice(ctx, request)
}

// FenceGet godoc
// @Summary  Returns traces along an arbitrary path, such as a well-path
// @description.markdown fence
// @Tags     fence
// @Param    query  query  string  True  "Urlencoded/escaped FenceResponse"
// @Accept   application/json
// @Produce  multipart/mixed
// @Success  200 {object} vds.FenceMetadata "(Example below only for metadata part)"
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /fence  [get]
func (e *Endpoint) FenceGet(ctx *gin.Context) {
	var request FenceRequest
	err := parseGetRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.fence(ctx, request)
}

// FencePost godoc
// @Summary  Returns traces along an arbitrary path, such as a well-path
// @description.markdown fence
// @Tags     fence
// @Param    body  body  FenceRequest  True  "Request Parameters"
// @Accept   application/json
// @Produce  multipart/mixed
// @Success  200 {object} vds.FenceMetadata "(Example below only for metadata part)"
// @Failure  400 {object} ErrorResponse "Request is invalid"
// @Failure  500 {object} ErrorResponse "openvds failed to process the request"
// @Router   /fence  [post]
func (e *Endpoint) FencePost(ctx *gin.Context) {
	var request FenceRequest
	err := parsePostRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.fence(ctx, request)
}

// AttributesPost godoc
// @Summary  Returns horizon attributes
// @description.markdown attribute
// @Tags     horizon
// @Param    body  body  AttributeRequest  True  "Request Parameters"
// @Accept   application/json
// @Produce  multipart/mixed
// @Router   /horizon  [post]
func (e *Endpoint) AttributesPost(ctx *gin.Context) {
	var request AttributeRequest
	err := parsePostRequest(ctx, &request)
	if abortOnError(ctx, err) { return }

	e.attributes(ctx, request)
}
