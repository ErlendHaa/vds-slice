package main

import (
	"fmt"
	"os"
	"strings"
	"strconv"

	"github.com/gin-gonic/gin"
	"github.com/pborman/getopt/v2"
	"github.com/swaggo/gin-swagger"
	"github.com/swaggo/files"
	"github.com/gin-contrib/gzip"

	_ "github.com/equinor/vds-slice/docs"
	"github.com/equinor/vds-slice/api"
	"github.com/equinor/vds-slice/internal/vds"
	"github.com/equinor/vds-slice/internal/cache"
	"github.com/equinor/vds-slice/internal/logging"
)

type opts struct {
	storageAccounts string
	port            string
	cacheSize       uint32
}

func parseCacheSize(cacheSize string) uint32 {
	if len(cacheSize) == 0 {
		return 0
	}
	out, err := strconv.ParseUint(cacheSize, 10, 32)
	if err != nil {
		panic(err)
	}

	return uint32(out)
}

func parseopts() opts {
	help := getopt.BoolLong("help", 0, "print this help text")
	
	opts := opts{
		storageAccounts: os.Getenv("STORAGE_ACCOUNTS"),
		port:            "8080",
		cacheSize:       parseCacheSize(os.Getenv("VDSSLICE_CACHE_SIZE")),
	}

	getopt.FlagLong(
		&opts.storageAccounts,
		"storage-accounts",
		0,
		"Comma-separated list of storage accounts that should be accepted by the API. " +
		"E.g. https://<account1>.blob.core.windows.net,https://<account2>.blob.core.windows.net",
		"string",
	)

	getopt.FlagLong(
		&opts.port,
		"port",
		0,
		"Port to start server on. Defaults to 8080",
	)

	getopt.FlagLong(
		&opts.cacheSize,
		"cache-size",
		0,
		"Max size of the response cache. In megabytes. A value of zero effectively\n" +
		"disables caching. Defaults to the value of the environment variable\n" +
		"VDSSLICE_CACHE_SIZE, or zero if the env var is not set.",
		"string",
	)

	getopt.Parse()
	if *help {
		getopt.Usage()
		os.Exit(0)
	}

	return opts
}

// @title        VDS-slice API
// @version      0.0
// @description  Serves seismic slices and fences from VDS files.
// @contact.name Equinor ASA
// @contact.url  https://github.com/equinor/vds-slice/issues
// @license.name GNU Affero General Public License
// @license.url  https://www.gnu.org/licenses/agpl-3.0.en.html
// @schemes      https
func main() {
	opts := parseopts()

	storageAccounts := strings.Split(opts.storageAccounts, ",")
	
	endpoint := api.Endpoint{
		MakeVdsConnection: vds.MakeAzureConnection(storageAccounts),
		Cache:             cache.NewCache(opts.cacheSize),
	}

	app := gin.New()
	app.Use(logging.FormattedLogger())
	app.Use(gin.Recovery())
	app.Use(gzip.Gzip(gzip.BestSpeed))

	app.GET("/", endpoint.Health)

	app.GET("metadata", api.ErrorHandler, endpoint.MetadataGet)
	app.POST("metadata", api.ErrorHandler, endpoint.MetadataPost)

	app.GET("slice", api.ErrorHandler, endpoint.SliceGet)
	app.POST("slice", api.ErrorHandler, endpoint.SlicePost)

	app.GET("fence", api.ErrorHandler, endpoint.FenceGet)
	app.POST("fence", api.ErrorHandler, endpoint.FencePost)

	app.GET("/swagger/*any", ginSwagger.WrapHandler(swaggerFiles.Handler))
	app.Run(fmt.Sprintf(":%s", opts.port))
}
