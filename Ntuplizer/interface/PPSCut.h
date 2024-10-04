#include "aperture_param_v2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TF2.h"
#include <iostream>
#include <string>
bool IsWithinAperture(const std::string &period, const std::string &arm, double xangle, double xi, double th_x_star)
{
    TF2 aperture("a", GetApertureParametrisation_version2(period, arm).c_str());
    
    aperture.SetParameter("xangle", xangle);
    aperture.SetParameter("xi", xi);
    const double th_x_star_upper_limit = -aperture.EvalPar(nullptr);

    return (th_x_star < th_x_star_upper_limit);
}
/*
void fillFiducialCutsVectors(
    TString era, std::map<std::pair<int, int>, double> &fiducialXLow_,
    std::map<std::pair<int, int>, double> &fiducialXHigh_,
    std::map<std::pair<int, int>, double> &fiducialYLow_,
    std::map<std::pair<int, int>, double> &fiducialYHigh_) {
  if (era == "2017B") {
    fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.098;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.698;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.698;
  } else if (era == "2017C1") {
    fiducialXLow_[std::pair<int, int>(0, 2)] = 1.860;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.334;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.098;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.698;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.698;
  } else if (era == "2017E") {
    fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -10.098;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.998;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -9.698;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.398;
  } else if (era == "2017F1") {
    fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -10.098;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.998;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -9.698;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.398;
  } else if (era == "2018A") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.710;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.927;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -11.598;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 3.698;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -10.898;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.398;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 18.498;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -11.298;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.098;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.420;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 20.045;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.398;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.098;
  } else if (era == "2018B1") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.850;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.927;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -11.598;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 3.698;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.420;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -10.798;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 18.070;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -11.198;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.098;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.420;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 25.045;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.398;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.098;
  } else if (era == "2018B2") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.562;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.640;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.198;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.135;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.798;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.498;
  } else if (era == "2018C") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.564;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.930;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.198;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.698;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.398;
  } else if (era == "2018D1") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.847;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.930;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.098;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.698;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.398;
  } else if (era == "2018D2") {
    fiducialXLow_[std::pair<int, int>(0, 0)] = 2.847;
    fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.931;
    fiducialYLow_[std::pair<int, int>(0, 0)] = -10.598;
    fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.498;
    fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
    fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
    fiducialYLow_[std::pair<int, int>(0, 2)] = -11.598;
    fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.398;
    fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
    fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
    fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
    fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
    fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
    fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
    fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
    fiducialYHigh_[std::pair<int, int>(1, 2)] = 3.898;
  } else
    std::cout << "WARNING: Era not recognized!" << std::endl;
  return;
}*/
