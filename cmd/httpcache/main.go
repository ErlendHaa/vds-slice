package main

import (
	"bufio"
	"bytes"
	"crypto/sha1"
	"encoding/hex"
	"fmt"
	"io/ioutil"
	"net/http"
	"net/http/httputil"
	"net/url"
	"os"
	"strconv"
	"unsafe"

	"github.com/gin-gonic/gin"
	"github.com/pborman/getopt/v2"

	"github.com/equinor/vds-slice/internal/cache"
)

type opts struct {
    port      uint32
    cacheSize int64
}

func parseAsUint32(fallback uint32, value string) uint32 {
    if len(value) == 0 {
		return fallback
    }
    out, err := strconv.ParseUint(value, 10, 32)
    if err != nil {
		panic(err)
    }

    return uint32(out)
}

func parseAsInt64(fallback int64, value string) int64 {
    if len(value) == 0 {
		return fallback
    }
    out, err := strconv.ParseInt(value, 10, 64)
    if err != nil {
		panic(err)
    }

    return out
}

func parseopts() opts {
    help := getopt.BoolLong("help", 0, "print this help text")
    
    opts := opts{
	port:      parseAsUint32(8079, os.Getenv("VDSSLICE_PROXY_PORT")),
	cacheSize: parseAsInt64(1024,  os.Getenv("VDSSLICE_PROXY_CACHE_SIZE")),

    }

    getopt.FlagLong(
	&opts.port,
		"port",
		0,
		"Port to start server on. Defaults to 8079.\n" +
		"Can also be set by environment variable 'VDSSLICE_HTTP_CACHE_PORT'",
		"int",
    )

    getopt.FlagLong(
	&opts.cacheSize,
		"cache-size",
		0,
		"Max size of the response cache. In megabytes. Defaults to 1024.\n" +
		"Can also be set by environment variable 'VDSSLICE_HTTP_CACHE_SIZE'",
		"int",
    )

    getopt.Parse()
    if *help {
		getopt.Usage()
		os.Exit(0)
    }

    return opts
}

func clone(resp *http.Response) (*http.Response, error) {
    // Clone the response
    clone := new(http.Response)
    *clone = *resp

    bodyBytes, err := ioutil.ReadAll(resp.Body)
    if err != nil {
		return nil, err
    }

    resp.Body.Close()

    // The body of a response is read-once. In other words, we have already
    // consumed the body of our original response and hence it must be reset to
    // it's original state for furter use
    resp.Body  = ioutil.NopCloser(bytes.NewBuffer(bodyBytes))
    clone.Body = ioutil.NopCloser(bytes.NewBuffer(bodyBytes))

    return clone, nil
}

func hashURL(input *url.URL) string {
    // Don't include the query parameters in the hash (i.e. the sas-token)
    x := *input
    x.RawQuery = ""

    hasher := sha1.New()
    hasher.Write([]byte(x.String()))
    return hex.EncodeToString(hasher.Sum(nil))
}

type cacheItem struct {
	request  []byte
	response []byte
}

func (c *cacheItem) Size() int64 {
	return int64(len(c.request) + len(c.response) + int(unsafe.Sizeof(*c)))
}

// TODO proper cache
type httpRequestCache struct {
    cache *cache.RistrettoCache
}

func (h *httpRequestCache) get(key string) (*http.Response, bool) {
	value, hit := h.cache.Cache.Get(key)
	if !hit {
		return nil, false
	}

	item := value.(cacheItem)

	req, err := http.ReadRequest( 
		bufio.NewReader(bytes.NewReader(item.request)),
	)

	if err != nil {
		fmt.Println(err)
		return nil, false
	}

	resp, err := http.ReadResponse( 
		bufio.NewReader(bytes.NewReader(item.response)),
		req,
	)
	if err != nil {
		fmt.Println(err)
		return nil, false
	}

	return resp, true
}

func (h *httpRequestCache) set(key string, resp *http.Response) error {
	rawRequest, err := httputil.DumpRequest(resp.Request, true)
	if err != nil {
		return err
	}
	rawResponse, err := httputil.DumpResponse(resp, true)
	if err != nil {
		return err
	}

	item := cacheItem{
		request:  rawRequest,
		response: rawResponse,
	}

	h.cache.Cache.Set(key, item, item.Size())
	return nil
}

func (h *httpRequestCache) ginHandler(ctx *gin.Context) {
    host := ctx.Param("host")
    path := ctx.Param("path")
    sas  := ctx.Request.URL.RawQuery

    remote, _ := url.Parse(fmt.Sprintf("https://%s%s?%s", host, path, sas))
    proxy := httputil.NewSingleHostReverseProxy(remote)

	key := hashURL(remote)
	cachedResponse, hit := h.get(key)

	// Update request before forwarding
    proxy.Director = func(req *http.Request) {
		req.Host       = remote.Host
		req.URL.Host   = remote.Host
		req.URL.Scheme = remote.Scheme
		req.URL.Path   = remote.Path

		if hit {
			req.Header.Add("If-Modified-Since", cachedResponse.Header.Get("Last-Modified"))
		}
	}

	// Update or cache response before forwarding
    proxy.ModifyResponse = func(resp *http.Response) error {
		switch (resp.StatusCode) {
		case http.StatusNotModified:
			// fmt.Println("Blob not modified, returning cached response...")
			*resp = *cachedResponse
		case http.StatusOK:
			// fmt.Println("Request was not in cache")
			h.set(key, resp)
		default:
			// Request failed, so no need to cache nor update the response. Return as is
		}

		return nil
    }

    proxy.ServeHTTP(ctx.Writer, ctx.Request)
}

// TODO: might want to generalize the cache package to better fit this usecase
// too
func NewHTTPRequestCache(cacheSize int64) *httpRequestCache {
	if cacheSize == 0 {
		panic("Cachesize must be greater than 0")
	}

    return &httpRequestCache{
		cache: cache.NewRistrettoCache(cacheSize),
    }
}

func main() {
    opts := parseopts()

    fmt.Printf("%+v", opts)
    reqCache := NewHTTPRequestCache(opts.cacheSize * 1024 * 1024)

    app := gin.New()
    app.SetTrustedProxies(nil)

    app.GET("/:host/*path", reqCache.ginHandler)
    app.Run(fmt.Sprintf(":%d", opts.port))
}
