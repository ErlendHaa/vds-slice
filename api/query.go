package api

import (
	"errors"
	"fmt"
	"log"
	"net/http"
	"strings"

	"github.com/gin-gonic/gin"

	"github.com/equinor/vds-slice/internal/vds"
)

type SliceQuery struct {
	Vds       string `form:"vds"       json:"vds"       binding:"required"`
	Direction string `form:"direction" json:"direction" binding:"required"`
	Lineno    *int   `form:"lineno"    json:"lineno"    binding:"required"`
}

func getAxis(direction string) (int, error) {
	switch direction {
		case "x":      return vds.AxisX,      nil
		case "y":      return vds.AxisY,      nil
		case "z":      return vds.AxisZ,      nil
		case "iline":  return vds.AxisIline,  nil
		case "xline":  return vds.AxisXline,  nil
		case "depth":  return vds.AxisDepth,  nil
		case "time":   return vds.AxisTime,   nil
		case "sample": return vds.AxisSample, nil
		default:
			msg := fmt.Sprintf("Invalid direction %v", direction)
			return 0, errors.New(msg)
	}
}

type Endpoint struct {
	StorageURL string
}

func (e *Endpoint) Health(ctx *gin.Context) {
	ctx.String(http.StatusOK, "I am up and running")
}

func (e *Endpoint) Slice(ctx *gin.Context) {
	var query SliceQuery

	if err := ctx.ShouldBind(&query); err != nil {
		ctx.AbortWithError(http.StatusBadRequest, err)
		return
	}

	querystr := ctx.Request.URL.Query()

	delete(querystr, "vds")
	delete(querystr, "direction")
	delete(querystr, "lineno")

	url := fmt.Sprintf("azure://%v", query.Vds)

	/*
	 * This is super buggy as assumes that no other query-paramters are
	 * present.
	 */
	cred := fmt.Sprintf(
		"BlobEndpoint=%v;SharedAccessSignature=?%v",
		e.StorageURL,
		querystr.Encode(),
	)

	axis, err := getAxis(strings.ToLower(query.Direction))

	if err != nil {
		ctx.AbortWithError(http.StatusBadRequest, err)
		return
	}

	buffer, err := vds.Slice(url, cred, *query.Lineno, axis)
	if err != nil {
		log.Println(err)
		ctx.AbortWithStatus(http.StatusInternalServerError)
		return
	}

	/*
	 * TODO: How should the data be returned? Returning a raw bytestream like
	 * this does not provide the caller with enough information to re-create
	 * the 2D array.
	 */
	ctx.Data(http.StatusOK, "application/octet-stream", buffer)
}
