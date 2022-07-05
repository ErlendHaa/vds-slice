package api

import (
	"errors"
	"fmt"
	"log"
	"net/http"

	"github.com/gin-gonic/gin"

	"github.com/equinor/vds-slice/internal/vds"
)

type SliceQuery struct {
	Vds       string `form:"vds"       json:"vds"       binding:"required"`
	Direction string `form:"direction" json:"direction" binding:"required"`
	Lineno    int    `form:"lineno"    json:"lineno"   `
}

var directions = map[string]string{
	"iline": "Inline",
	"xline": "Crossline",
	"sampe": "Sample",
}

func vdsdirection(direction string) (string, error) {
	if val, ok := directions[direction]; ok {
		return val, nil
	}

	msg := "Invalid direction. Expected 'iline', 'xline' or 'sample', got: %v"
	return "", errors.New(fmt.Sprintf(msg, direction))
}

type Endpoint struct {
	StorageURL string
}

func (e *Endpoint) Health(ctx *gin.Context) {
	ctx.String(http.StatusOK, "I am up and running")
}

func (e *Endpoint) SliceGet(ctx *gin.Context) {
	var query SliceQuery

	if err := ctx.ShouldBind(&query); err != nil {
		ctx.AbortWithError(http.StatusBadRequest, err)
		return
	}

	dir, err := vdsdirection(query.Direction)
	if err != nil {
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

	buffer, err := vds.Slice(url, cred, dir, query.Lineno)
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
