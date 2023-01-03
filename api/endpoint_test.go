package api

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"mime"
	"mime/multipart"
	"net/http"
	"net/http/httptest"
	"net/url"
	"reflect"
	"strings"
	"testing"

	"github.com/gin-gonic/gin"
	"github.com/stretchr/testify/require"

	"github.com/equinor/vds-slice/internal/vds"
)

const well_known = "../testdata/well_known/well_known_default.vds"

type sliceTest struct {
	name           string
	method         string
	slice          testSliceRequest
	sliceRequest   string
	expectedStatus int
	expectedError  string
}

type fenceTest struct {
	name           string
	method         string
	fence          testFenceRequest
	fenceRequest   string
	expectedStatus int
	expectedError  string
}

type metadataTest struct {
	name            string
	method          string
	metadata        testMetadataRequest
	metadataRequest string
	expectedStatus  int
	expectedError   string
}

// define own help types to assure separation between production and test code
type testErrorResponse struct {
	Error string `json:"error" binding:"required"`
}

type testSliceRequest struct {
	Vds       string `json:"vds"`
	Direction string `json:"direction"`
	Lineno    int    `json:"lineno"`
	Sas       string `json:"sas"`
}

type testFenceRequest struct {
	Vds              string      `json:"vds"`
	CoordinateSystem string      `json:"coordinateSystem"`
	Coordinates      [][]float32 `json:"coordinates"`
	Sas              string      `json:"sas"`
}

type testMetadataRequest struct {
	Vds string `json:"vds"`
	Sas string `json:"sas"`
}

type testSliceAxis struct {
	Annotation string  `json:"annotation" binding:"required"`
	Max        float32 `json:"max"        binding:"required"`
	Min        float32 `json:"min"        binding:"required"`
	Samples    int     `json:"samples"    binding:"required"`
	Unit       string  `json:"unit"       binding:"required"`
}

type testSliceMetadata struct {
	X      testSliceAxis `json:"x"      binding:"required"`
	Y      testSliceAxis `json:"y"      binding:"required"`
	Format string        `json:"format" binding:"required"`
}

func TestSliceHappyHTTPResponse(t *testing.T) {
	testcases := []sliceTest{
		{
			name:   "Valid GET Request",
			method: http.MethodGet,
			slice: testSliceRequest{
				Vds:       well_known,
				Direction: "i",
				Lineno:    0, //side-effect assurance that 0 is accepted
				Sas:       "n/a",
			},
			expectedStatus: http.StatusOK,
		},
		{
			name:   "Valid json POST Request",
			method: http.MethodPost,
			slice: testSliceRequest{
				Vds:       well_known,
				Direction: "crossline",
				Lineno:    10,
				Sas:       "n/a",
			},
			expectedStatus: http.StatusOK,
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareSliceRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		parts := readMultipartData(t, w)
		if len(parts) != 2 {
			msg := "Got %d parts in reply; want it to always contain 2 in case '%s'"
			t.Errorf(
				msg,
				len(parts),
				testcase.name,
			)
		}

		inlineAxis := testSliceAxis{
			Annotation: "Inline", Max: 5.0, Min: 1.0, Samples: 3, Unit: "unitless",
		}
		crosslineAxis := testSliceAxis{
			Annotation: "Crossline", Max: 11.0, Min: 10.0, Samples: 2, Unit: "unitless",
		}
		sampleAxis := testSliceAxis{
			Annotation: "Sample", Max: 16.0, Min: 4.0, Samples: 4, Unit: "ms",
		}
		expectedFormat := "<f4"

		var expectedMetadata *testSliceMetadata
		switch testcase.slice.Direction {
		case "i":
			expectedMetadata = &testSliceMetadata{
				X:      sampleAxis,
				Y:      crosslineAxis,
				Format: expectedFormat}
		case "crossline":
			expectedMetadata = &testSliceMetadata{
				X:      sampleAxis,
				Y:      inlineAxis,
				Format: expectedFormat}
		default:
			t.Fatalf("Unhandled direction %s in case %s", testcase.slice.Direction, testcase.name)
		}

		metadata := &testSliceMetadata{}
		err := json.Unmarshal(parts[0], metadata)

		if err != nil {
			msg := "Failed json metadata extraction in case '%s'"
			t.Fatalf(
				msg,
				testcase.name,
			)
		}

		if !reflect.DeepEqual(metadata, expectedMetadata) {
			msg := "Got %v as metadata; want it to be %v in case '%s'"
			t.Fatalf(
				msg,
				metadata,
				expectedMetadata,
				testcase.name,
			)
		}

		expectedDataLength := expectedMetadata.X.Samples *
			expectedMetadata.Y.Samples * 4 //4 bytes each
		if len(parts[1]) != expectedDataLength {
			msg := "Got %d bytes in data reply; want it to be %d in case '%s'"
			t.Errorf(
				msg,
				len(parts[2]),
				expectedDataLength,
				testcase.name,
			)
		}
	}
}

func TestSliceErrorHTTPResponse(t *testing.T) {
	testcases := []sliceTest{
		{
			name:           "Invalid json GET request",
			method:         http.MethodGet,
			sliceRequest:   "help I am a duck",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "invalid character",
		},
		{
			name:           "Invalid json POST request",
			method:         http.MethodPost,
			sliceRequest:   "help I am a duck",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "invalid character",
		},
		{
			name:   "Missing parameters GET request",
			method: http.MethodGet,
			sliceRequest: "{\"vds\":\"" + well_known +
				"\", \"direction\":\"i\", \"sas\": \"n/a\"}",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "Error:Field validation for 'Lineno'",
		},
		{
			name:   "Missing parameters POST Request",
			method: http.MethodPost,
			sliceRequest: "{\"vds\":\"" + well_known +
				"\", \"lineno\":1, \"sas\": \"n/a\"}",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "Error:Field validation for 'Direction'",
		},
		{
			name:   "Request with unknown axis",
			method: http.MethodPost,
			slice: testSliceRequest{
				Vds:       well_known,
				Direction: "unknown",
				Lineno:    1,
				Sas:       "n/a",
			},
			expectedStatus: http.StatusBadRequest,
			expectedError:  "invalid direction 'unknown', valid options are",
		},
		{
			name:   "Request which passed all input checks but still should fail",
			method: http.MethodPost,
			slice: testSliceRequest{
				Vds:       well_known,
				Direction: "i",
				Lineno:    10,
				Sas:       "n/a",
			},
			expectedStatus: http.StatusInternalServerError,
			expectedError:  "Invalid lineno: 10, valid range: [0:2:1]",
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareSliceRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		assertError(w, t, testcase.name, testcase.expectedError)
	}
}

func TestFenceHappyHTTPResponse(t *testing.T) {
	testcases := []fenceTest{
		{
			name:   "Valid GET Request",
			method: http.MethodGet,
			fence:  testFenceRequest{
				Vds:              well_known,
				CoordinateSystem: "ilxl",
				Coordinates:      [][]float32{{3, 11}, {2, 10}},
				Sas:              "n/a",
			},
			expectedStatus: http.StatusOK,
		},
		{
			name:   "Valid json POST Request",
			method: http.MethodPost,
			fence:  testFenceRequest{
				Vds:              well_known,
				CoordinateSystem: "ij",
				Coordinates:      [][]float32{{0, 1}, {1, 1}, {1, 0}},
				Sas:              "n/a",
			},
			expectedStatus: http.StatusOK,
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareFenceRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		parts := readMultipartData(t, w)
		if len(parts) != 2 {
			msg := "Got %d parts in reply; want it to always contain 3 in case '%s'"
			t.Errorf(
				msg,
				len(parts),
				testcase.name,
			)
		}

		metadata := string(parts[0])
		expectedMetadata := `{
			"shape": [` + fmt.Sprint(len(testcase.fence.Coordinates)) + `, 4],
			"format": "<f4"
		}`

		if metadata != expectedMetadata {
			msg := "Metadata not equal in case '%s'"
			require.JSONEq(t, expectedMetadata, metadata, fmt.Sprintf(msg, testcase.name))
		}

		expectedDataLength := len(testcase.fence.Coordinates) * 4 * 4 //4 bytes, 4 samples per each requested
		if len(parts[1]) != expectedDataLength {
			msg := "Got %d bytes in data reply; want it to be %d in case '%s'"
			t.Errorf(
				msg,
				len(parts[2]),
				expectedDataLength,
				testcase.name,
			)
		}
	}
}

func TestFenceErrorHTTPResponse(t *testing.T) {
	testcases := []fenceTest{
		{
			name:           "Invalid json GET request",
			method:         http.MethodGet,
			fenceRequest:   "help I am a duck",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "invalid character",
		},
		{
			name:           "Invalid json POST request",
			method:         http.MethodPost,
			fenceRequest:   "help I am a duck",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "invalid character",
		},
		{
			name:   "Missing parameters GET request",
			method: http.MethodGet,
			fenceRequest: "{\"vds\":\"" + well_known +
				"\", \"coordinateSystem\":\"ilxl\", \"coordinates\":[[0, 0]]}",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "Error:Field validation for 'Sas'",
		},
		{
			name:   "Missing parameters POST Request",
			method: http.MethodPost,
			fenceRequest: "{\"vds\":\"" + well_known +
				"\", \"coordinateSystem\":\"ilxl\", \"sas\": \"n/a\"}",
			expectedStatus: http.StatusBadRequest,
			expectedError:  "Error:Field validation for 'Coordinates'",
		},
		{
			name:   "Request with unknown coordinate system",
			method: http.MethodPost,
			fence:  testFenceRequest{
				Vds:              well_known,
				CoordinateSystem: "unknown",
				Coordinates:      [][]float32{{3, 12}, {2, 10}},
				Sas:              "n/a",
			},
			expectedStatus: http.StatusBadRequest,
			expectedError:  "coordinate system not recognized: 'unknown', valid options are",
		},
		{
			name:   "Request which passed all input checks but still should fail",
			method: http.MethodPost,
			fence:  testFenceRequest{
				Vds:              "unknown",
				CoordinateSystem: "ilxl",
				Coordinates:      [][]float32{{3, 12}, {2, 10}},
				Sas:              "n/a",
			},
			expectedStatus: http.StatusInternalServerError,
			expectedError:  "Could not open VDS",
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareFenceRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		assertError(w, t, testcase.name, testcase.expectedError)
	}
}

func TestMetadataHappyHTTPResponse(t *testing.T) {
	testcases := []metadataTest{
		{
			name:   "Valid GET Request",
			method: http.MethodGet,
			metadata: testMetadataRequest{
				Vds: well_known,
				Sas: "n/a",
			},
			expectedStatus: http.StatusOK,
		},
		{
			name:   "Valid json POST Request",
			method: http.MethodPost,
			metadata: testMetadataRequest{
				Vds: well_known,
				Sas: "n/a",
			},
			expectedStatus: http.StatusOK,
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareMetadataRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		metadata := w.Body.String()
		expectedMetadata := `{
			"axis": [
				{"annotation": "Inline", "max": 5.0, "min": 1.0, "samples" : 3, "unit": "unitless"},
				{"annotation": "Crossline", "max": 11.0, "min": 10.0, "samples" : 2, "unit": "unitless"},
				{"annotation": "Sample", "max": 16.0, "min": 4.0, "samples" : 4, "unit": "ms"}
			],
			"boundingBox": {
				"cdp": [[5,0],[9,8],[4,11],[0,3]],
				"ilxl": [[1, 10], [5, 10], [5, 11], [1, 11]],
				"ij": [[0, 0], [2, 0], [2, 1], [0, 1]]
			},
			"crs": "utmXX",
			"format": "<f4"
		}`

		if metadata != expectedMetadata {
			msg := "Metadata not equal in case '%s'"
			require.JSONEq(t, expectedMetadata, metadata, fmt.Sprintf(msg, testcase.name))
		}
	}
}

func TestMetadataErrorHTTPResponse(t *testing.T) {
	testcases := []metadataTest{
		{
			name:            "Invalid json GET request",
			method:          http.MethodGet,
			metadataRequest: "help I am a duck",
			expectedStatus:  http.StatusBadRequest,
			expectedError:   "invalid character",
		},
		{
			name:            "Invalid json POST request",
			method:          http.MethodPost,
			metadataRequest: "help I am a duck",
			expectedStatus:  http.StatusBadRequest,
			expectedError:   "invalid character",
		},
		{
			name:            "Missing parameters GET request",
			method:          http.MethodGet,
			metadataRequest: "{\"vds\":\"" + well_known + "\"}",
			expectedStatus:  http.StatusBadRequest,
			expectedError:   "Error:Field validation for 'Sas'",
		},
		{
			name:            "Missing parameters POST Request",
			method:          http.MethodPost,
			metadataRequest: "{\"sas\":\"somevalidsas\"}",
			expectedStatus:  http.StatusBadRequest,
			expectedError:   "Error:Field validation for 'Vds'",
		},
		{
			name:     "Request which passed all input checks but still should fail",
			method:   http.MethodPost,
			metadata: testMetadataRequest{
				Vds: "unknown",
				Sas: "n/a",
			},
			expectedStatus: http.StatusInternalServerError,
			expectedError:  "Could not open VDS",
		},
	}

	for _, testcase := range testcases {
		w := httptest.NewRecorder()
		ctx, r := gin.CreateTestContext(w)
		setupTestServer(r)

		prepareMetadataRequest(ctx, t, testcase)
		r.ServeHTTP(w, ctx.Request)

		if w.Result().StatusCode != testcase.expectedStatus {
			msg := "Got status %v; want %d %s in case '%s'"
			t.Errorf(
				msg,
				w.Result().Status,
				testcase.expectedStatus,
				http.StatusText(testcase.expectedStatus),
				testcase.name,
			)
		}
		assertError(w, t, testcase.name, testcase.expectedError)
	}
}

func MakeFileConnection() vds.ConnectionMaker {
	return func(path, sas string) (vds.ResourceConnection, error) {
		path = fmt.Sprintf("file://%s", path)
		return vds.NewFileConnection(path), nil
	}
}

func setupTestServer(r *gin.Engine) {
	r.Use(gin.Logger())
	endpoint := Endpoint{ MakeVdsConnection: MakeFileConnection() }

	r.POST("/slice", ErrorHandler, endpoint.SlicePost)
	r.GET("/slice", ErrorHandler, endpoint.SliceGet)

	r.POST("/fence", ErrorHandler, endpoint.FencePost)
	r.GET("/fence", ErrorHandler, endpoint.FenceGet)

	r.POST("/metadata", ErrorHandler, endpoint.MetadataPost)
	r.GET("/metadata", ErrorHandler, endpoint.MetadataGet)
}

func assertError(w *httptest.ResponseRecorder, t *testing.T, test_name string,
	expected_error string) {

	testErrorInfo := &testErrorResponse{}
	err := json.Unmarshal(w.Body.Bytes(), testErrorInfo)
	if err != nil {
		t.Fatalf("Couldn't unmarshal data")
	}

	if !strings.Contains(testErrorInfo.Error, expected_error) {
		msg := "Got error string %s; want it to contain '%s' in case '%s'"
		t.Errorf(
			msg,
			testErrorInfo.Error,
			expected_error,
			test_name,
		)
	}
}

func readMultipartData(t *testing.T, w *httptest.ResponseRecorder) [][]byte {
	_, params, err := mime.ParseMediaType(w.Result().Header.Get("Content-Type"))
	if err != nil {
		t.Fatalf("Cannot parse Content Type")
	}
	mr := multipart.NewReader(w.Body, params["boundary"])

	parts := [][]byte{}
	for {
		p, err := mr.NextPart()
		if err == io.EOF {
			return parts
		}
		if err != nil {
			t.Fatalf("Couldn't process part")
		}
		data, err := io.ReadAll(p)
		if err != nil {
			t.Fatalf("Couldn't read part")
		}
		parts = append(parts, data)
	}
}

func prepareRequest(ctx *gin.Context, t *testing.T, endpoint string,
	method string, request interface{}, jsonReq string) {
	jsonRequest := jsonReq
	if jsonRequest == "" {
		req, err := json.Marshal(request)
		if err != nil {
			t.Fatalf("Cannot marshal request")
		}
		jsonRequest = string(req)
	}

	switch method := method; method {
	case http.MethodGet:
		ctx.Request, _ = http.NewRequest(
			http.MethodGet,
			endpoint,
			nil,
		)

		q := url.Values{}
		q.Add("query", jsonRequest)
		ctx.Request.URL.RawQuery = q.Encode()

	case http.MethodPost:
		ctx.Request, _ = http.NewRequest(
			http.MethodPost,
			endpoint,
			bytes.NewBuffer([]byte(jsonRequest)),
		)
		ctx.Request.Header.Set("Content-Type", "application/json")
	default:
		t.Fatalf("Unknown method")
	}
}

func prepareSliceRequest(ctx *gin.Context, t *testing.T, testcase sliceTest) {
	prepareRequest(
		ctx,
		t,
		"/slice",
		testcase.method,
		testcase.slice,
		testcase.sliceRequest,
	)
}

func prepareFenceRequest(ctx *gin.Context, t *testing.T, testcase fenceTest) {
	prepareRequest(
		ctx,
		t,
		"/fence",
		testcase.method,
		testcase.fence,
		testcase.fenceRequest,
	)
}

func prepareMetadataRequest(ctx *gin.Context, t *testing.T, testcase metadataTest) {
	prepareRequest(
		ctx,
		t,
		"/metadata",
		testcase.method,
		testcase.metadata,
		testcase.metadataRequest,
	)
}
