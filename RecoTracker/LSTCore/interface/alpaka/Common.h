#ifndef RecoTracker_LSTCore_interface_alpaka_Common_h
#define RecoTracker_LSTCore_interface_alpaka_Common_h

#include <numbers>

#include "FWCore/Utilities/interface/HostDeviceConstant.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoTracker/LSTCore/interface/Common.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::lst {

  using namespace ::lst;

  ALPAKA_FN_HOST ALPAKA_FN_INLINE void lstWarning(std::string_view warning) {
#ifdef LST_STANDALONE
    printf("%s\n", warning.data());
#else
    edm::LogWarning("LST") << warning;
#endif
  }

  // The constants below are usually used in functions like alpaka::math::min(),
  // expecting a reference (T const&) in the arguments. Hence,
  // HOST_DEVICE_CONSTANT needs to be used instead of constexpr.

  HOST_DEVICE_CONSTANT float kPi = std::numbers::pi_v<float>;
  // 15 MeV constant from the approximate Bethe-Bloch formula
  HOST_DEVICE_CONSTANT float kMulsInGeV = 0.015;
  HOST_DEVICE_CONSTANT float kMiniMulsPtScaleBarrel[6] = {0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};
  HOST_DEVICE_CONSTANT float kMiniMulsPtScaleEndcap[5] = {0.006, 0.006, 0.006, 0.006, 0.006};
  HOST_DEVICE_CONSTANT float kMiniRminMeanBarrel[6] = {
      25.007152356, 37.2186993757, 52.3104270826, 68.6658656666, 85.9770373007, 108.301772384};
  HOST_DEVICE_CONSTANT float kMiniRminMeanEndcap[5] = {
      130.992832231, 154.813883559, 185.352604327, 221.635123002, 265.022076742};
  HOST_DEVICE_CONSTANT float k2Rinv1GeVf = (2.99792458e-3 * 3.8) / 2;
  HOST_DEVICE_CONSTANT float kR1GeVf = 1. / (2.99792458e-3 * 3.8);
  HOST_DEVICE_CONSTANT float kSinAlphaMax = 0.95;
  HOST_DEVICE_CONSTANT float kDeltaZLum = 15.0;
  HOST_DEVICE_CONSTANT float kPixelPSZpitch = 0.15;
  HOST_DEVICE_CONSTANT float kStripPSZpitch = 2.4;
  HOST_DEVICE_CONSTANT float kStrip2SZpitch = 5.0;
  HOST_DEVICE_CONSTANT float kWidth2S = 0.009;
  HOST_DEVICE_CONSTANT float kWidthPS = 0.01;
  HOST_DEVICE_CONSTANT float kPt_betaMax = 7.0;
  // To be updated with std::numeric_limits<float>::infinity() in the code and data files
  HOST_DEVICE_CONSTANT float kVerticalModuleSlope = 123456789.0;

  HOST_DEVICE_CONSTANT float kMiniDeltaTilted[3] = {0.26f, 0.26f, 0.26f};
  HOST_DEVICE_CONSTANT float kMiniDeltaFlat[6] = {0.26f, 0.16f, 0.16f, 0.18f, 0.18f, 0.18f};
  HOST_DEVICE_CONSTANT float kMiniDeltaLooseTilted[3] = {0.4f, 0.4f, 0.4f};
  HOST_DEVICE_CONSTANT float kMiniDeltaEndcap[5][15] = {
      {0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, /*10*/ 0.18f, 0.18f, 0.18f, 0.18f, 0.18f},
      {0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, /*10*/ 0.18f, 0.18f, 0.18f, 0.18f, 0.18f},
      {0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.18f, 0.18f, /*10*/ 0.18f, 0.18f, 0.18f, 0.18f, 0.18f},
      {0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.18f, 0.18f, /*10*/ 0.18f, 0.18f, 0.18f, 0.18f, 0.18f},
      {0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.18f, /*10*/ 0.18f, 0.18f, 0.18f, 0.18f, 0.18f}};

  namespace dnn {

    // Common constants for both DNNs
    HOST_DEVICE_CONSTANT float kEta_norm = 2.5f;
    HOST_DEVICE_CONSTANT float kPhi_norm = kPi;
    HOST_DEVICE_CONSTANT float kEtaSize = 0.25f;  // Bin size in eta.
    constexpr unsigned int kPtBins = 2;
    constexpr unsigned int kEtaBins = 10;

    namespace t3dnn {
      HOST_DEVICE_CONSTANT float kZ_max = 224.149505f;
      HOST_DEVICE_CONSTANT float kR_max = 98.932365f;
      HOST_DEVICE_CONSTANT float kWp_prompt[kPtBins][kEtaBins] = {
          {0.4957, 0.5052, 0.5201, 0.5340, 0.4275, 0.4708, 0.4890, 0.4932, 0.5400, 0.5449},
          {0.0302, 0.0415, 0.0994, 0.1791, 0.1960, 0.2467, 0.3227, 0.3242, 0.2367, 0.2187}};
      HOST_DEVICE_CONSTANT float kWp_displaced[kPtBins][kEtaBins] = {
          {0.0334, 0.0504, 0.0748, 0.0994, 0.1128, 0.1123, 0.1118, 0.1525, 0.1867, 0.1847},
          {0.0091, 0.0075, 0.0350, 0.0213, 0.0435, 0.0676, 0.1957, 0.1649, 0.1080, 0.1046}};
    }  // namespace t3dnn

    namespace t5dnn {
      HOST_DEVICE_CONSTANT float kZ_max = 267.2349854f;
      HOST_DEVICE_CONSTANT float kR_max = 110.1099396f;
      HOST_DEVICE_CONSTANT float kWp[kPtBins][kEtaBins] = {
          {0.4493, 0.4939, 0.5715, 0.6488, 0.5709, 0.5938, 0.7164, 0.7565, 0.8103, 0.8593},
          {0.4488, 0.4448, 0.5067, 0.5929, 0.4836, 0.4112, 0.4968, 0.4403, 0.5597, 0.5067}};
    }  // namespace t5dnn

    namespace pt3dnn {
      HOST_DEVICE_CONSTANT float kWp[kEtaBins] = {
          0.189, 0.1805, 0.2267, 0.3104, 0.4719, 0.3159, 0.1372, 0.1571, 0.3198, 0.186};
      HOST_DEVICE_CONSTANT float kWpHigh = 0.0473f;
    }  // namespace pt3dnn

  }  // namespace dnn

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::lst
#endif
