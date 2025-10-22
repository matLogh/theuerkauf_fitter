#pragma once

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include <TCanvas.h>
#include <TCollection.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TRatioPlot.h>
#include <TVirtualFitter.h>
#pragma GCC diagnostic pop

#include "Util.hpp"
#include "tabulate.hpp"

class TheuerkaufPeak
{
    friend class TheuerkaufFitter;

  public:
    enum ParamState
    {
        FREE,
        FIXED,
        NONE,
        SAME
    };

  public:
    TheuerkaufPeak() : fFcn(nullptr){};
    TheuerkaufPeak(double min, double max, int id = 0, bool hasTL = false, bool hasTR = false, bool hasStep = false);
    TheuerkaufPeak(const TheuerkaufPeak &peak);
    TheuerkaufPeak(const TheuerkaufPeak *peak);
    TheuerkaufPeak(TheuerkaufPeak &&peak);
    ~TheuerkaufPeak();

    TheuerkaufPeak &operator=(const TheuerkaufPeak &other);

    void Print() const;

    int GetID() const noexcept
    {
        return fId;
    };

    int GetIndex_Volume() const noexcept;
    int GetIndex_Position() const noexcept;
    int GetIndex_Sigma() const noexcept;
    int GetIndex_TailLeft() const noexcept;
    int GetIndex_TailRight() const noexcept;
    int GetIndex_StepHeight() const noexcept;
    int GetIndex_StepWidth() const noexcept;

    void Draw(std::string options = "")
    {
        auto temp_param_index = fParamIndex;
        this->ResetIndexes(true);
        fFcn->Draw(options.c_str());
        fParamIndex = temp_param_index;
    };

    void SetProperty_Volume(const ParamState prop) noexcept
    {
        fParamState[0] = prop;
    };
    void SetProperty_Position(const ParamState prop) noexcept
    {
        fParamState[1] = prop;
    };
    void SetProperty_Sigma(const ParamState prop) noexcept
    {
        fParamState[2] = prop;
    };
    void SetProperty_TailLeft(const ParamState prop) noexcept
    {
        fParamState[3] = prop;
    };
    void SetProperty_TailRight(const ParamState prop) noexcept
    {
        fParamState[4] = prop;
    };
    void SetProperty_StepHeight(const ParamState prop) noexcept
    {
        fParamState[5] = prop;
    };
    void SetProperty_StepWidth(const ParamState prop) noexcept
    {
        fParamState[6] = prop;
    };

    /// @brief It is expected to get full range of parameters in this order:
    ///     p[fParamIndex[0]] = volume
    ///     p[fParamIndex[1]] = position
    ///     p[fParamIndex[2]] = sigma
    ///     p[fParamIndex[3]] = left tail
    ///     p[fParamIndex[4]] = right tail
    ///     p[fParamIndex[5]] = step height
    ///     p[fParamIndex[6]] = step width
    /// @param x
    /// @param p
    /// @return
    double Eval(const double *x, const double *p) const;
    TF1 *GetFunction() const noexcept
    {
        return fFcn;
    };

    /// @brief These functions are relevant only with Theuerkauf fitter, or when peak
    /// function is coupled with another function and you want to share some parameters
    /// between them, such as sigma etc
    /// @param index
    /// @return
    TheuerkaufPeak *SetIndex_Volume(const int index);
    TheuerkaufPeak *SetIndex_Position(const int index);
    TheuerkaufPeak *SetIndex_Sigma(const int index);
    TheuerkaufPeak *SetIndex_TailLeft(const int index);
    TheuerkaufPeak *SetIndex_TailRight(const int index);
    TheuerkaufPeak *SetIndex_StepHeight(const int index);
    TheuerkaufPeak *SetIndex_StepWidth(const int index);

    TheuerkaufPeak *SetParameter_Volume(double val, ParamState prop = FREE, double min = -1e9, double max = 1e9);
    TheuerkaufPeak *SetParameter_Position(double val, ParamState prop = FREE, double min = -1e9, double max = 1e9);
    TheuerkaufPeak *SetParameter_Sigma(double val, ParamState prop = FREE,
                                       double min = std::numeric_limits<double>::epsilon(), double max = 1e9);
    TheuerkaufPeak *SetParameter_TailLeft(double val = 10., ParamState prop = NONE, double min = 0, double max = 1e12);
    TheuerkaufPeak *SetParameter_TailRight(double val = 10., ParamState prop = NONE, double min = 0, double max = 1e12);
    TheuerkaufPeak *SetParameter_StepHeight(double val, ParamState prop = NONE, double min = 0, double max = 1e9);
    TheuerkaufPeak *SetParameter_StepWidth(double val, ParamState prop = NONE, double min = 0, double max = 1e9);

    bool HasStep() const noexcept
    {
        return fHasStep;
    };
    bool HasTL() const noexcept
    {
        return fHasLeftTail;
    };
    bool HasTR() const noexcept
    {
        return fHasRightTail;
    };

    TheuerkaufPeak *SetRange(double min, double max);

    ParamState GetState(int index) const;
    ParamState GetState_Volume() const noexcept
    {
        return fParamState[0];
    };
    ParamState GetState_Position() const noexcept
    {
        return fParamState[1];
    };
    ParamState GetState_Sigma() const noexcept
    {
        return fParamState[2];
    };
    ParamState GetState_TailLeft() const noexcept
    {
        return fParamState[3];
    };
    ParamState GetState_TailRight() const noexcept
    {
        return fParamState[4];
    };
    ParamState GetState_StepHeight() const noexcept
    {
        return fParamState[5];
    };
    ParamState GetState_StepWidth() const noexcept
    {
        return fParamState[6];
    };

    double GetVol() const noexcept
    {
        return fHistBinning_normalization <= 0 ? fFcn->GetParameter(0)
                                               : fFcn->GetParameter(0) * fHistBinning_normalization;
        // return fFcn->GetParameter(0);
    };
    double GetPos() const noexcept
    {
        return fFcn->GetParameter(1);
    };
    double GetSig() const noexcept
    {
        return fFcn->GetParameter(2);
    };
    double GetTL() const noexcept
    {
        return fFcn->GetParameter(3);
    };
    double GetTR() const noexcept
    {
        return fFcn->GetParameter(4);
    };
    double GetSH() const noexcept
    {
        return fFcn->GetParameter(5);
    };
    double GetSW() const noexcept
    {
        return fFcn->GetParameter(6);
    };

    double GetFWHM() const noexcept
    {
        return sig_to_fwhm * this->GetSig();
    };

    // get width at x-th maximum (x=2 for FWHM, x=10 for FW at tenth maximum etc.)
    double GetFWxM(const double width_multiple);

    double GetVolErr() const noexcept
    {
        return fHistBinning_normalization <= 0 ? fFcn->GetParError(0)
                                               : fFcn->GetParError(0) * fHistBinning_normalization;
    };
    double GetPosErr() const noexcept
    {
        return fFcn->GetParError(1);
    };
    double GetSigErr() const noexcept
    {
        return fFcn->GetParError(2);
    };
    double GetTLErr() const noexcept
    {
        return fFcn->GetParError(3);
    };
    double GetTRErr() const noexcept
    {
        return fFcn->GetParError(4);
    };
    double GetSHErr() const noexcept
    {
        return fFcn->GetParError(5);
    };
    double GetSWErr() const noexcept
    {
        return fFcn->GetParError(6);
    };

    double GetFWHMErr() const noexcept
    {
        return sig_to_fwhm * this->GetSigErr();
    };

    void ResetIndexes(bool disregard_id = false) const;
    void SetBinning(const double &energy_per_bin)
    {
        fHistBinning_normalization = 1. / energy_per_bin;
    };
    double GetBinning() const noexcept
    {
        return 1. / fHistBinning_normalization;
    };

  private:
    double EvalNoStep(const double *x, const double *p) const;
    double EvalStep(const double *x, const double *p) const;
    double GetNorm(const double sigma, const double tl, const double tr) const;
    // void SetMinMax(double min, double max);

    int fId;
    double fXMin, fXMax;
    bool fHasLeftTail, fHasRightTail, fHasStep;
    TF1 *fFcn;
    double fHistBinning_normalization{-1.}; // necessary to get the real volume, cleared of the

    static const double sig_to_fwhm;

    mutable double fCachedNorm{-1.}, fCachedSigma{-1.}, fCachedTL{-1.}, fCachedTR{-1.};
    mutable std::array<int, 7> fParamIndex;
    mutable std::array<ParamState, 7> fParamState;
};

//!  THEUERKAUF FITTER
//!  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TheuerkaufFitter
{
  public:
    TheuerkaufFitter(double min, double max);
    ~TheuerkaufFitter();

    /// @brief Add peak to fitter
    /// @param position rought position of the peak
    /// @param leftTail has left tail?
    /// @param rightTail has right tail?
    /// @param step has step?
    /// @return id of the peak
    int AddPeak(const double position, const bool leftTail, const bool rightTail, const bool step);

    /// @brief This function will create a copy of the peak and set id to one of the
    /// internal structure
    /// @param peak
    /// @return id of added peak
    int AddPeak(const TheuerkaufPeak &peak);
    TheuerkaufFitter *SetBackground(std::unique_ptr<TF1> &bcg_fcn) noexcept;

    /// @brief get peak based on its ID
    /// @param id
    /// @return
    std::shared_ptr<TheuerkaufPeak> GetPeak(const int id);

    /// @brief Get number of peaks
    /// @return
    int GetNPeaks();

    /// @brief Set internal polynomial background
    /// @param poly_order number of parameters (poly_order=0 no background, poly_order=3
    /// quadratic etc)
    /// @return
    TheuerkaufFitter *SetBackgroundPoly(unsigned int poly_order) noexcept;

    /// @brief remove function set by the TheuerkaufFitter::SetBackground, internal
    /// polynomial background remains
    /// @return
    TheuerkaufFitter *RemoveBackground() noexcept;

    /// @brief If true, allows peaks to have negative volume
    /// @param allow_negative_peaks
    /// @return
    TheuerkaufFitter *AllowNegativePeaks(bool allow_negative_peaks) noexcept;

    /// @brief Change range of the fitter
    /// @param min X min
    /// @param max X max
    /// @return
    TheuerkaufFitter *SetRange(double min, double max);

    /// @brief Get X range of the fitter
    void GetRange(double &min, double &max) const noexcept;
    double GetXmin() const noexcept
    {
        return fXMin;
    };

    double GetXmax() const noexcept
    {
        return fXMax;
    };

    double Eval(const double *x, const double *p) const;

    /// @brief Fit histogram. Only here the total function is being computed.
    /// @param histFit histogram to be fitted
    /// @param options
    ///        "OUTPUT_NONE" - do not print fit results,
    ///        "OUTPUT_STANDARD" - print standard output,
    ///        "OUTPUT_ROOT" - print ROOT fit output,
    ///        "OUTPUT_PLAIN" - plain output without table (id, pos, pos_err, fwhm,
    ///        fwhm_err, vol, vol_err) "LIKELIHOOD" - use likelihood in fitting, "POISSON"
    ///        - use chi2 in fitting
    void Fit(TH1 *histFit, std::string options = "");

    /// @brief Draw canvas with 2 pads: first with fitted spectrum, summed fit function
    /// and individual peaks; second pad with residuals, confidence intervals and 2 sigma
    /// bands
    /// @param histAna
    /// @param toPrint
    void Analyze(TH1 *histAna);

    /// @brief Draws the fit function and the background function.
    /// @param hist histogram to be drawn - if set tu nullptr it is assumed that histogram is already plotted. If
    /// provided, function creates copy of the histogram and draws it.
    /// @param toPrint pad to draw the fit function on - if set to nullptr, a new canvas will be drawn
    void DrawFit(TH1 *hist = nullptr, TVirtualPad *toPrint = nullptr);

    TFitResultPtr GetFitResults()
    {
        return fFitResults;
    };

    std::shared_ptr<TF1> GetTotalFunction()
    {
        return fSumFunc;
    };

    // std::vector<std::pair<double,double>> GetPeaksIntegral(); //Filippo
    void PrintFitResults() const;
    void PrintFitResults_plain(std::ostream &os = std::cout) const;

    ///@brief Set output verbosity
    ///@param verbose 0 - do not print anything, 1 - standard print, 2 - print all
    void SetVerbosity(const int verbose)
    {
        fVerbose = verbose;
    };

  private:
    double EvalTotalBackground(const double *x, const double *p);
    int GetNumParams() const noexcept
    {
        return static_cast<int>(fPeaks.size() * 7 + fPolyBcgDegree);
    };
    void HandleParameterStates();
    void DistributeParametersToPeaks();

  private:
    int fVerbose;
    double fXMin, fXMax;
    double fChiSquare;
    bool fOnlypositivepeaks;
    unsigned int fPolyBcgDegree;
    TFitResultPtr fFitResults;

    std::shared_ptr<TF1> fSumFunc;
    std::unique_ptr<TF1> fBcgFunc;
    std::shared_ptr<TH1> fTempHist;

    std::vector<std::shared_ptr<TheuerkaufPeak>> fPeaks;

    std::vector<std::shared_ptr<TheuerkaufPeak>> fTempPeaks;

    // handle graphical objects
  private:
    std::vector<std::shared_ptr<TObject>> fTempObjects;

  public:
    /// @brief Get graphical objects that were created during the fit - this is to ensure they will continue existing
    /// after fitter object is deleted
    std::vector<std::shared_ptr<TObject>> &GetGraphicalObjects()
    {
        return fTempObjects;
    }

  public:
    void GetConfidenceIntervals(unsigned int n, const double *x, double *ci, double cl = 0.95, bool norm = true);
    void GetConfidenceIntervals(TH1 *hfit, double cl = 0.95);

    static double GetMinimumInRange(const std::shared_ptr<TH1> h, double x_min, double x_max) noexcept;
    static double GetMaximumInRange(const std::shared_ptr<TH1> h, double x_min, double x_max) noexcept;
};