import {
  getRandomIndexInDimension,
  getMaxIndexInDimension,
  convertDimension,
} from "./helpers/request-helpers.js";
import { sendSetupMetadataRequest } from "./helpers/metadata-helpers.js";
import { sendSliceRequest } from "./helpers/slice-helpers.js";
import { createSummary, thresholds, summaryTrendStats } from "./helpers/report-helpers.js";

import exec from "k6/execution";

export const options = {
  scenarios: {
    constantSlice: {
      executor: "constant-vus",
      vus: 1,
      duration: "2m",
    },
  },
  thresholds: thresholds(),
  summaryTrendStats: summaryTrendStats(),
};

export function setup() {
  const annotatedDimension = "Inline";
  const metadata = sendSetupMetadataRequest();
  const startIndex = getRandomIndexInDimension(metadata, annotatedDimension);
  const indexDimension = convertDimension(annotatedDimension);
  const maxIndex = getMaxIndexInDimension(metadata, annotatedDimension);
  return [startIndex, maxIndex, indexDimension];
}

export default function (params) {
  const [startIndex, maxIndex, indexDimension] = params;
  const index = (startIndex + exec.scenario.iterationInInstance) % maxIndex;
  sendSliceRequest(indexDimension, index);
}

export function handleSummary(data) {
  return createSummary(data);
}
