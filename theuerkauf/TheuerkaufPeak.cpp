#include "TheuerkaufPeak.h"

#include <TLegend.h>
#include <TStyle.h>
#include <iomanip>

const double TheuerkaufPeak::sig_to_fwhm = 2. * TMath::Sqrt(2. * TMath::Log(2.));

double TheuerkaufPeak::GetFWxM(const double width_multiple)
{
    auto temp_param_index = fParamIndex;
    this->ResetIndexes(true);
    const double fraction_max = fFcn->GetMaximum() / width_multiple;
    double left = fFcn->GetX(fraction_max, fXMin, this->GetPos(), 1E-5, 1000, false);
    double right = fFcn->GetX(fraction_max, this->GetPos(), fXMax, 1E-5, 1000, false);

    fParamIndex = temp_param_index;
    return right - left;
}

TheuerkaufPeak &TheuerkaufPeak::operator=(const TheuerkaufPeak &other)
{
    if (this != &other) // protect against self-assignment
    {
        // Directly assign each member from 'other' to 'this'
        this->fId = other.fId;
        this->fXMin = other.fXMin;
        this->fXMax = other.fXMax;
        this->fHasLeftTail = other.fHasLeftTail;
        this->fHasRightTail = other.fHasRightTail;
        this->fHasStep = other.fHasStep;
        this->fParamIndex = other.fParamIndex;
        this->fParamState = other.fParamState;
        this->fHistBinning_normalization = other.fHistBinning_normalization;
        this->fCachedNorm = other.fCachedNorm;
        this->fCachedSigma = other.fCachedSigma;
        this->fCachedTL = other.fCachedTL;
        this->fCachedTR = other.fCachedTR;

        if (fFcn != nullptr)
        {
            delete fFcn;
            fFcn = nullptr;
        }

        std::string fcn_name = "theurekauf_" + std::to_string(fId);
        this->fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, fXMin,
                             fXMax, 7, "TheuerkaufPeak", "Eval");
        fFcn->SetBit(kCanDelete);

        this->fFcn->SetNpx(10000);
        this->fFcn->SetLineColor(GetColor());
        for (int i = 0; i < 7; i++)
        {
            double min, max;
            other.fFcn->GetParLimits(i, min, max);
            this->fFcn->SetParName(i, other.fFcn->GetParName(i));
            this->fFcn->SetParLimits(i, min, max);
            this->fFcn->SetParameter(i, other.fFcn->GetParameter(i));
        }
    }
    return *this;
}

TheuerkaufPeak::TheuerkaufPeak(TheuerkaufPeak &&peak)
{
    this->fId = peak.fId;
    this->fXMin = peak.fXMin;
    this->fXMax = peak.fXMax;
    this->fHasLeftTail = peak.fHasLeftTail;
    this->fHasRightTail = peak.fHasRightTail;
    this->fHasStep = peak.fHasStep;
    this->fParamIndex = peak.fParamIndex;
    this->fParamState = peak.fParamState;

    this->fHistBinning_normalization = peak.fHistBinning_normalization;
    this->fCachedNorm = peak.fCachedNorm;
    this->fCachedSigma = peak.fCachedSigma;
    this->fCachedTL = peak.fCachedTL;
    this->fCachedTR = peak.fCachedTR;

    // we need to make a new function to reference the ::Eval of THIS member
    std::string fcn_name = "theurekauf_" + std::to_string(fId);
    this->fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, fXMin, fXMax,
                         7, "TheuerkaufPeak", "Eval");
    fFcn->SetBit(kCanDelete);

    this->fFcn->SetNpx(10000);
    this->fFcn->SetLineColor(GetColor());
    for (int i = 0; i < 7; i++)
    {
        double min, max;
        peak.fFcn->GetParLimits(i, min, max);
        this->fFcn->SetParName(i, peak.fFcn->GetParName(i));
        this->fFcn->SetParLimits(i, min, max);
        this->fFcn->SetParameter(i, peak.fFcn->GetParameter(i));
    }
}

TheuerkaufPeak::TheuerkaufPeak(const TheuerkaufPeak &peak)
{
    this->fId = peak.fId;
    this->fXMin = peak.fXMin;
    this->fXMax = peak.fXMax;
    this->fHasLeftTail = peak.fHasLeftTail;
    this->fHasRightTail = peak.fHasRightTail;
    this->fHasStep = peak.fHasStep;
    this->fParamIndex = peak.fParamIndex;
    this->fParamState = peak.fParamState;

    this->fHistBinning_normalization = peak.fHistBinning_normalization;
    this->fCachedNorm = peak.fCachedNorm;
    this->fCachedSigma = peak.fCachedSigma;
    this->fCachedTL = peak.fCachedTL;
    this->fCachedTR = peak.fCachedTR;

    // we need to make a new function to reference the ::Eval of THIS member
    std::string fcn_name = "theurekauf_" + std::to_string(fId);
    this->fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, fXMin, fXMax,
                         7, "TheuerkaufPeak", "Eval");
    fFcn->SetBit(kCanDelete);

    this->fFcn->SetNpx(10000);
    this->fFcn->SetLineColor(GetColor());
    for (int i = 0; i < 7; i++)
    {
        double min, max;
        peak.fFcn->GetParLimits(i, min, max);
        this->fFcn->SetParName(i, peak.fFcn->GetParName(i));
        this->fFcn->SetParLimits(i, min, max);
        this->fFcn->SetParameter(i, peak.fFcn->GetParameter(i));
    }
}

TheuerkaufPeak::TheuerkaufPeak(const TheuerkaufPeak *peak)
{
    this->fId = peak->fId;
    this->fXMin = peak->fXMin;
    this->fXMax = peak->fXMax;
    this->fHasLeftTail = peak->fHasLeftTail;
    this->fHasRightTail = peak->fHasRightTail;
    this->fHasStep = peak->fHasStep;
    this->fParamIndex = peak->fParamIndex;
    this->fParamState = peak->fParamState;

    this->fHistBinning_normalization = peak->fHistBinning_normalization;
    this->fCachedNorm = peak->fCachedNorm;
    this->fCachedSigma = peak->fCachedSigma;
    this->fCachedTL = peak->fCachedTL;
    this->fCachedTR = peak->fCachedTR;
    this->fParamIndex = peak->fParamIndex;
    this->fParamState = peak->fParamState;

    // we need to make a new function to reference the ::Eval of THIS member
    std::string fcn_name = "theurekauf_" + std::to_string(fId);
    this->fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, fXMin, fXMax,
                         7, "TheuerkaufPeak", "Eval");
    fFcn->SetBit(kCanDelete);

    this->fFcn->SetNpx(10000);
    this->fFcn->SetLineColor(GetColor());
    for (int i = 0; i < 7; i++)
    {
        double min, max;
        peak->fFcn->GetParLimits(i, min, max);
        this->fFcn->SetParName(i, peak->fFcn->GetParName(i));
        this->fFcn->SetParLimits(i, min, max);
        this->fFcn->SetParameter(i, peak->fFcn->GetParameter(i));
    }
}

TheuerkaufPeak::TheuerkaufPeak(double min, double max, int id, bool hasTL, bool hasTR, bool hasStep)
    : fId(id), fXMin(min), fXMax(max), fHasLeftTail(hasTL), fHasRightTail(hasTR), fHasStep(hasStep)
{
    assert(id >= 0);
    std::string fcn_name = "theurekauf_" + std::to_string(id);
    fFcn = new TF1(GetFuncUniqueName(fcn_name.c_str(), this).c_str(), this, &TheuerkaufPeak::Eval, min, max, 7,
                   "TheuerkaufFitter", "Eval");
    fFcn->SetBit(kCanDelete);

    fFcn->SetNpx(10000);
    fFcn->SetLineColor(GetColor());

    this->ResetIndexes();
    this->SetParameter_Volume(0., FREE)
        ->SetParameter_Position((min + max) / 2, FREE)
        ->SetParameter_Sigma(std::numeric_limits<double>::epsilon(), FREE)
        ->SetParameter_TailLeft(0., hasTL ? FREE : NONE)
        ->SetParameter_TailRight(0., hasTR ? FREE : NONE)
        ->SetParameter_StepHeight(0., hasStep ? FREE : NONE)
        ->SetParameter_StepWidth(0., hasStep ? FREE : NONE);
};

TheuerkaufPeak::~TheuerkaufPeak()
{
    if (fFcn)
    {
        delete fFcn;
        fFcn = nullptr;
    }
}

int TheuerkaufPeak::GetIndex_Volume() const noexcept
{
    return fParamIndex[0];
};
int TheuerkaufPeak::GetIndex_Position() const noexcept
{
    return fParamIndex[1];
};
int TheuerkaufPeak::GetIndex_Sigma() const noexcept
{
    return fParamIndex[2];
};
int TheuerkaufPeak::GetIndex_TailLeft() const noexcept
{
    return fParamIndex[3];
};
int TheuerkaufPeak::GetIndex_TailRight() const noexcept
{
    return fParamIndex[4];
};
int TheuerkaufPeak::GetIndex_StepHeight() const noexcept
{
    return fParamIndex[5];
};
int TheuerkaufPeak::GetIndex_StepWidth() const noexcept
{
    return fParamIndex[6];
};

TheuerkaufPeak *TheuerkaufPeak::SetIndex_Volume(const int index)
{
    fParamIndex[0] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_Position(const int index)
{
    fParamIndex[1] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_Sigma(const int index)
{
    fParamIndex[2] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_TailLeft(const int index)
{
    fParamIndex[3] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_TailRight(const int index)
{
    fParamIndex[4] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_StepHeight(const int index)
{
    fParamIndex[5] = index;
    return this;
};
TheuerkaufPeak *TheuerkaufPeak::SetIndex_StepWidth(const int index)
{
    fParamIndex[6] = index;
    return this;
};

TheuerkaufPeak *TheuerkaufPeak::SetRange(double min, double max)
{
    assert(min < max);

    fXMin = min;
    fXMax = max;
    fFcn->SetRange(fXMin, fXMax);

    double _min, _max;
    fFcn->GetParLimits(1, _min, _max);
    if (_min < min)
        _min = min;
    if (_max > max)
        _max = max;
    fFcn->SetParLimits(1, _min, _max);

    return this;
}

TheuerkaufPeak::ParamState TheuerkaufPeak::GetState(int index) const
{
    assert(index >= 0 && index < 7);
    return fParamState[index];
};

//  void TheuerkaufPeak::SetMinMax(double min, double max)
// {
//     assert(min < max);
//     fXMin = min;
//     fXMax = max;
//     fFcn->SetRange(min, max);
// };

void TheuerkaufPeak::Print() const
{
    std::cout << std::endl;
    std::cout << "TheuerkaufPeak model, peakID " << fId << std::endl;
    std::cout << "Function name " << fFcn->GetName() << std::endl << std::endl;
    std::cout << "PARAM   STATUS   VAL   MIN   MAX" << std::endl;

    for (int i = 0; i < 7; i++)
    {
        std::cout << fFcn->GetParName(i) << " ";
        std::cout << "index " << fParamIndex[i] << " ";
        if (GetState(i) == FREE)
            std::cout << "FREE  ";
        if (GetState(i) == FIXED)
            std::cout << "FIXED ";
        if (GetState(i) == NONE)
            std::cout << "NONE  ";
        if (GetState(i) == SAME)
            std::cout << "SAME  ";
        std::cout << fFcn->GetParameter(i) << " ";
        double min, max;
        fFcn->GetParLimits(i, min, max);
        std::cout << min << " " << max << " " << std::endl;
    }
    std::cout << std::endl;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_Volume(double val, ParamState prop, double min, double max)
{
    assert(prop != NONE);
    assert(val >= min && val <= max);

    int index = 0;
    std::string name = std::string("vol_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());

    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    default:
        throw std::runtime_error("It is not allowed to have volume set as NONE or SAME!");
    }
    fFcn->SetParLimits(index, min, max);

    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_Position(double val, ParamState prop, double min, double max)
{
    assert(prop != NONE && prop != SAME);
    if (min < fXMin)
        min = fXMin;
    if (max > fXMax)
        max = fXMax;
    if (val < min || val > max)
    {
        throw std::runtime_error("Position value " + std::to_string(val) + " is out of range [" + std::to_string(min) +
                                 ", " + std::to_string(max) + "]");
    }

    int index = 1;
    std::string name = std::string("pos_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    default:
        throw std::runtime_error("It is not allowed to have volume set as NONE or SAME!");
    }
    fFcn->SetParLimits(index, min, max);

    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_Sigma(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 2;
    std::string name = std::string("sig_") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        break;
    case NONE:
        std::cerr << "Setting sigma to NONE will turns peak to delta function, are you "
                     "sure this is what you wanted to do?"
                  << std::endl;
        fFcn->FixParameter(index, std::numeric_limits<double>::epsilon());
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_TailLeft(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 3;
    std::string name = std::string("TL__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fFcn->SetParLimits(index, min, max);
        fHasLeftTail = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasLeftTail = true;
        break;
    case NONE:
        fFcn->FixParameter(index, std::numeric_limits<double>::max());
        fHasLeftTail = false;
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fFcn->SetParLimits(index, min, max);
        fHasLeftTail = true;
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_TailRight(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 4;
    std::string name = std::string("TR__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasRightTail = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasRightTail = true;
        break;
    case NONE:
        fHasRightTail = false;
        fFcn->FixParameter(index, std::numeric_limits<double>::max());
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasRightTail = true;
        break;
    }

    fFcn->SetParLimits(index, min, max);
    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_StepHeight(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 5;
    std::string name = std::string("SH__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    bool hasStep = (prop != NONE);
    if (fHasStep != hasStep)
    {
        fHasStep = hasStep;
        this->SetParameter_StepWidth(1., prop);
    }

    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasStep = true;
        break;
    case NONE:
        fFcn->FixParameter(index, 0.);
        fHasStep = false;
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    }
    fFcn->SetParLimits(index, min, max);
    return this;
}

TheuerkaufPeak *TheuerkaufPeak::SetParameter_StepWidth(double val, ParamState prop, double min, double max)
{
    assert(val >= min && val <= max);

    int index = 6;
    std::string name = std::string("SW__") + "_" + std::to_string(fId);
    fFcn->SetParName(index, name.c_str());
    fParamState[index] = prop;
    bool hasStep = (prop != NONE);
    if (fHasStep != hasStep)
    {
        fHasStep = hasStep;
        this->SetParameter_StepHeight(1., prop);
    }

    switch (prop)
    {
    case FREE:
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    case FIXED:
        fFcn->FixParameter(index, val);
        fHasStep = true;
        break;
    case NONE:
        fHasStep = false;
        fFcn->FixParameter(index, 0.);
        break;
    case SAME: // not really applicable for single peak, but matters for TheuerkaufFitter
        fFcn->ReleaseParameter(index);
        fFcn->SetParameter(index, val);
        fHasStep = true;
        break;
    }
    fFcn->SetParLimits(index, min, max);
    return this;
}

double TheuerkaufPeak::Eval(const double *x, const double *p) const
{
    return EvalNoStep(x, p) + EvalStep(x, p);
}

double TheuerkaufPeak::EvalNoStep(const double *x, const double *p) const
{
    double vol = p[fParamIndex[0]];
    double dx = *x - p[fParamIndex[1]];
    double sigma = p[fParamIndex[2]];

    double tl = p[fParamIndex[3]];
    double tr = p[fParamIndex[4]];
    double norm = GetNorm(sigma, tl, tr);
    double _x;

    // Peak function
    if (dx < -tl && fHasLeftTail)
    {
        _x = tl / (sigma * sigma) * (dx + tl / 2.0);
    }
    else if (dx < tr || !fHasRightTail)
    {
        _x = -dx * dx / (2.0 * sigma * sigma);
    }
    else
    {
        _x = -tr / (sigma * sigma) * (dx - tr / 2.0);
    }

    return vol * norm * std::exp(_x);
}

double TheuerkaufPeak::EvalStep(const double *x, const double *p) const
{
    //! Step function

    if (fHasStep)
    {
        double vol = p[fParamIndex[0]];
        double dx = x[0] - p[fParamIndex[1]];
        double sigma = p[fParamIndex[2]];
        double tl = p[fParamIndex[3]];
        double tr = p[fParamIndex[4]];
        double sh = p[fParamIndex[5]];
        double sw = p[fParamIndex[6]];

        double norm = this->GetNorm(sigma, tl, tr);

        return vol * norm * sh * (M_PI / 2. + std::atan(sw * dx / (std::sqrt(2.) * sigma)));
    }
    else
    {
        return 0.0;
    }
}

double TheuerkaufPeak::GetNorm(const double sigma, const double tl, const double tr) const
{
    if (fCachedSigma == sigma && fCachedTL == tl && fCachedTR == tr)
    {
        return fCachedNorm;
    }

    double vol;

    // Contribution from left tail + left half of truncated gaussian
    if (fHasLeftTail)
    {
        vol = (sigma * sigma) / tl * std::exp(-(tl * tl) / (2.0 * sigma * sigma));
        vol += std::sqrt(M_PI / 2.0) * sigma * std::erf(tl / (std::sqrt(2.0) * sigma));
    }
    else
    {
        vol = std::sqrt(M_PI / 2.0) * sigma;
    }

    // Contribution from right tail + right half of truncated gaussian
    if (fHasRightTail)
    {
        vol += (sigma * sigma) / tr * std::exp(-(tr * tr) / (2.0 * sigma * sigma));
        vol += std::sqrt(M_PI / 2.0) * sigma * std::erf(tr / (std::sqrt(2.0) * sigma));
    }
    else
    {
        vol += std::sqrt(M_PI / 2.0) * sigma;
    }

    fCachedSigma = sigma;
    fCachedTL = tl;
    fCachedTR = tr;
    // fCachedNorm  = fHistBinning_normalization <= 0 ? 1. / vol :
    // fHistBinning_normalization / vol;
    fCachedNorm = 1. / vol;

    return fCachedNorm;
}

void TheuerkaufPeak::ResetIndexes(bool disregard_id) const
{
    if (disregard_id)
        std::iota(fParamIndex.begin(), fParamIndex.end(), 0);
    else
        std::iota(fParamIndex.begin(), fParamIndex.end(), fId * 7);
};

//!  THEUERKAUF FITTER
//!  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TheuerkaufFitter::~TheuerkaufFitter()
{
    fTempHist.reset(); // Explicitly reset before ROOT cleanup
    fSumFunc.reset();
    fBcgFunc.reset();
}

TheuerkaufFitter::TheuerkaufFitter(double min, double max) : fXMin(min), fXMax(max)
{
    fVerbose = 1;
    fSumFunc = nullptr;
    fBcgFunc = nullptr;
    fPolyBcgDegree = 2;
    fOnlypositivepeaks = true;
    fChiSquare = std::numeric_limits<double>::quiet_NaN();
}

TheuerkaufFitter *TheuerkaufFitter::SetBackground(std::unique_ptr<TF1> &bcg_fcn) noexcept
{
    fBcgFunc = std::move(bcg_fcn);
    return this;
};

void TheuerkaufFitter::PrintFitResults_plain(std::ostream &os) const
{
    for (auto &peak : fPeaks)
    {
        std::string _id = std::to_string(peak->GetID());
        std::string _pos = std::to_string(peak->GetPos());
        std::string _pos_e = std::to_string(peak->GetPosErr());
        std::string _fwhm = std::to_string(peak->GetFWHM());
        std::string _fwhm_err = std::to_string(peak->GetFWHMErr());
        std::string _vol = std::to_string(peak->GetVol());
        std::string _vol_err = std::to_string(peak->GetVolErr());
        std::string _tl = peak->HasTL() ? std::to_string(peak->GetTL()) : "--";
        std::string _tr = peak->HasTR() ? std::to_string(peak->GetTR()) : "--";
        std::string _sw = peak->HasStep() ? std::to_string(peak->GetSW()) : "--";
        std::string _sh = peak->HasStep() ? std::to_string(peak->GetSH()) : "--";

        os << _id << " " << _pos << " " << _pos_e << " " << _fwhm << " " << _fwhm_err << " " << _vol << " " << _vol_err
           << std::endl;
    }
}

void TheuerkaufFitter::PrintFitResults() const
{
    // using namespace tabulate;
    assert(fSumFunc && fFitResults && "Fit was likely not yet performed!");
    tabulate::Table fit_results;

    std::string _status;

    if (fFitResults->IsValid())
        _status = "Fit OK";
    else
        _status = "Fit not valid, minimizer status " + std::to_string(fFitResults->Status());
    std::string _chi2 = "chi2 " + std::to_string(fChiSquare);
    std::string _rchi2 = "rchi2 " + std::to_string(fChiSquare / fSumFunc.get()->GetNDF());
    std::string _errmat = "error matrix ";

    tabulate::Color _matrix_status;
    switch (fFitResults->CovMatrixStatus())
    {
    case 0:
        _errmat += "not calculated";
        _matrix_status = tabulate::Color::red;
        break;
    case 1:
        _errmat += "approximated";
        _matrix_status = tabulate::Color::yellow;
        break;
    case 2:
        _errmat += "pos def";
        _matrix_status = tabulate::Color::yellow;
        break;
    case 3:
        _errmat += "accurate";
        _matrix_status = tabulate::Color::green;
        break;
    default:
        _errmat += "unknown";
        _matrix_status = tabulate::Color::red;
    }

    std::string _hname = fTempHist != nullptr ? fTempHist->GetName() : "unknown hist";

    if (fFitResults->IsValid())
    {
        fit_results.add_row({_hname, _status, _chi2, _rchi2, _errmat});
    }
    else
    {
        std::string _edm = "edm " + std::to_string(fFitResults->Edm());
        std::string _ncalls = "ncalls " + std::to_string(fFitResults->NCalls());
        fit_results.add_row({_hname, _status, _chi2, _rchi2, _edm, _ncalls, _errmat});
    }

    // fit_results.format()
    //            .border_top("*")
    //            .border_bottom("-")
    //            .border_left(" ")
    //            .border_right(" ");
    fit_results.row(0)
        .format()
        .background_color(tabulate::Color::white)
        .font_style({tabulate::FontStyle::bold})
        .font_align(tabulate::FontAlign::center)
        .background_color(tabulate::Color::grey)
        .font_color(tabulate::Color::white);
    if (fFitResults->IsValid())
        fit_results.row(0).format().background_color(tabulate::Color::green);
    else
        fit_results.row(0).format().background_color(tabulate::Color::red);
    fit_results[0][4].format().background_color(_matrix_status);
    fit_results[0][0].format().background_color(tabulate::Color::magenta);

    tabulate::Table peak_parameters;
    // peak_parameters.format()
    //            .border_top("*")
    //            .border_bottom("-")
    //            .border_left(" ")
    //            .border_right(" ");
    peak_parameters.add_row(
        {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Volume", "VolErr", "TL", "TR", "SW", "SH"});
    // check if we indeed have Volumes or just normalizations by checking the binwidth
    // parameter of each peak
    bool binning_set = std::find_if(fPeaks.begin(), fPeaks.end(),
                                    [](const auto &peak) { return peak->GetBinning() <= 0; }) == fPeaks.end();
    if (binning_set)
    {
        peak_parameters.add_row(
            {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Volume", "VolErr", "TL", "TR", "SW", "SH"});
    }
    else
    {
        peak_parameters.add_row(
            {"Peak ID", "Position", "PosError", "FWHM", "FWHMErr", "Normalization", "NormErr", "TL", "TR", "SW", "SH"});
    }

    // bool has_tl = std::find(fPeaks.begin(), fPeaks.end(), [](auto& peak) {return
    // peak->HasTL();}) == fPeaks.end() ? false : true; bool has_tr =
    // std::find(fPeaks.begin(), fPeaks.end(), [](auto& peak) {return peak->HasTR();}) ==
    // fPeaks.end() ? false : true; bool has_step = std::find(fPeaks.begin(),
    // fPeaks.end(), [](auto& peak) {return peak->HasStep();}) == fPeaks.end() ? false :
    // true;

    peak_parameters.row(0)
        .format()
        .font_style({tabulate::FontStyle::bold})
        .font_align(tabulate::FontAlign::center)
        .font_color(tabulate::Color::yellow);

    for (auto &peak : fPeaks)
    {
        std::string _id = std::to_string(peak->GetID());
        std::string _pos = std::to_string(peak->GetPos());
        std::string _pos_e = std::to_string(peak->GetPosErr());
        std::string _fwhm = std::to_string(peak->GetFWHM());
        std::string _fwhm_err = std::to_string(peak->GetFWHMErr());
        std::string _vol = std::to_string(peak->GetVol());
        std::string _vol_err = std::to_string(peak->GetVolErr());
        std::string _tl = peak->HasTL() ? std::to_string(peak->GetTL()) : "--";
        std::string _tr = peak->HasTR() ? std::to_string(peak->GetTR()) : "--";
        std::string _sw = peak->HasStep() ? std::to_string(peak->GetSW()) : "--";
        std::string _sh = peak->HasStep() ? std::to_string(peak->GetSH()) : "--";

        peak_parameters.add_row({_id, _pos, _pos_e, _fwhm, _fwhm_err, _vol, _vol_err, _tl, _tr, _sw, _sh});
    }

    peak_parameters.column(0).format().font_style({tabulate::FontStyle::bold}).font_align(tabulate::FontAlign::center);

    std::cout << std::endl << fit_results << std::endl;
    std::cout << peak_parameters << std::endl << std::endl;
}

std::shared_ptr<TheuerkaufPeak> TheuerkaufFitter::GetPeak(const int id)
{
    auto it = std::find_if(fPeaks.begin(), fPeaks.end(), [&](const auto &pk) { return pk->GetID() == id; });
    assert(it != fPeaks.end());

    return (*it);
}

int TheuerkaufFitter::GetNPeaks()
{
    return fPeaks.size();
}

TheuerkaufFitter *TheuerkaufFitter::SetRange(double min, double max)
{
    fXMin = min;
    fXMax = max;

    for (auto &peak : fPeaks)
    {
        peak->SetRange(min, max);
    }
    return this;
}

TheuerkaufFitter *TheuerkaufFitter::SetBackgroundPoly(unsigned int poly_order) noexcept
{
    fPolyBcgDegree = poly_order;
    return this;
};

TheuerkaufFitter *TheuerkaufFitter::RemoveBackground() noexcept
{
    fBcgFunc.reset(nullptr);
    return this;
}

TheuerkaufFitter *TheuerkaufFitter::AllowNegativePeaks(bool allow_negative_peaks) noexcept
{
    fOnlypositivepeaks = !allow_negative_peaks;
    return this;
};

int TheuerkaufFitter::AddPeak(const double position, const bool leftTail, const bool rightTail, const bool step)
{
    fPeaks.emplace_back(std::make_shared<TheuerkaufPeak>(fXMin, fXMax, fPeaks.size(), leftTail, rightTail, step));
    fPeaks.back()->SetParameter_Sigma(1., TheuerkaufPeak::ParamState::SAME);
    fPeaks.back()->SetParameter_Position(position, TheuerkaufPeak::ParamState::FREE, fXMin, fXMax);
    if (leftTail)
        fPeaks.back()->SetParameter_TailLeft(10., TheuerkaufPeak::ParamState::FREE);
    if (rightTail)
        fPeaks.back()->SetParameter_TailRight(10., TheuerkaufPeak::ParamState::FREE);
    if (step)
    {
        fPeaks.back()->SetParameter_StepHeight(0., TheuerkaufPeak::ParamState::FREE);
        fPeaks.back()->SetParameter_StepWidth(1., TheuerkaufPeak::ParamState::FREE);
    }

    return fPeaks.back()->fId;
}

int TheuerkaufFitter::AddPeak(const TheuerkaufPeak &peak)
{
    fPeaks.emplace_back(std::make_shared<TheuerkaufPeak>(peak));
    fPeaks.back()->fId = fPeaks.size() - 1;
    return fPeaks.back()->fId;
}

double TheuerkaufFitter::Eval(const double *x, const double *p) const
{
    //! Private: evaluation function for fit
    auto num_params = GetNumParams();

    // Evaluate background function, if it has been given
    double sum = fBcgFunc ? fBcgFunc->Eval(*x) : 0.0;

    // Evaluate internal background
    sum += std::accumulate(std::reverse_iterator<const double *>(p + num_params),
                           std::reverse_iterator<const double *>(p + num_params - fPolyBcgDegree), 0.0,
                           [&x](double bg, double param) { return bg * *x + param; });

    // Evaluate peaks
    return std::accumulate(fPeaks.begin(), fPeaks.end(), sum,
                           [x, p](double _sum, const auto &peak) { return _sum + peak->Eval(x, p); });
}

double TheuerkaufFitter::EvalTotalBackground(const double *x, const double *p)
{
    int num_params = GetNumParams();

    // eval bcg fcn
    double sum = fBcgFunc ? fBcgFunc->Eval(*x) : 0.0;

    sum += std::accumulate(std::reverse_iterator<const double *>(p + num_params),
                           std::reverse_iterator<const double *>(p + num_params - fPolyBcgDegree), 0.0,
                           [&x](double bg, double param) { return bg * *x + param; });
    // eval peak steps
    sum += std::accumulate(fPeaks.begin(), fPeaks.end(), 0.,
                           [x, p](double _sum, const auto &peak) { return _sum + peak->EvalStep(x, p); });

    return sum;
}

void TheuerkaufFitter::DrawFit(TH1 *hist, TVirtualPad *toPrint)
{
    if (!fSumFunc)
        return;
    fTempObjects.clear();
    fTempPeaks.clear();

    // some function-local things
    const int num_params = this->GetNumParams();
    const auto color = GetColor();
    const int fcn_npx = fSumFunc->GetNpx();

    // prepare canvas
    std::shared_ptr<TPad> can(nullptr);
    if (toPrint == nullptr)
    {
        can = std::make_shared<TCanvas>(GetFuncUniqueName("fit_canvas", this).c_str(), "Fitted spectrum", 1400, 600);
        fTempObjects.emplace_back(can);
    }

    // prepare histogram
    std::shared_ptr<TH1> _h(nullptr);
    if (hist != nullptr)
    {
        // auto hname = GetFuncUniqueName("_hist", this);
        _h = std::shared_ptr<TH1>(dynamic_cast<TH1 *>(hist->Clone()));
        _h->SetName(GetFuncUniqueName("_hist", _h.get()).c_str());
        _h->SetTitle("Energy spectrum");
        _h->SetDirectory(0);
        _h->Draw();
        fTempObjects.emplace_back(_h);
    }

    // create histogram out of a background function
    auto tot_bcg_fcn = std::make_unique<TF1>(GetFuncUniqueName("total_bcg", this).c_str(), this,
                                             &TheuerkaufFitter::EvalTotalBackground, fXMin, fXMax, num_params,
                                             "TheuerkaufFitter", "EvalTotalBackground");
    tot_bcg_fcn->SetNpx(fcn_npx);
    tot_bcg_fcn->SetParameters(fSumFunc->GetParameters());
    std::shared_ptr<TH1> tot_bcg_hist(dynamic_cast<TH1 *>(tot_bcg_fcn->GetHistogram()->Clone()));
    tot_bcg_hist->SetName(GetFuncUniqueName("tot_bcg_hist", tot_bcg_hist.get()).c_str());
    tot_bcg_hist->SetDirectory(0);
    fTempObjects.emplace_back(tot_bcg_hist);
    tot_bcg_hist->SetLineColor(color);
    tot_bcg_hist->SetLineWidth(0);
    tot_bcg_hist->SetFillStyle(3002);
    tot_bcg_hist->SetFillColor(color);

    std::shared_ptr<TH1> sum_hist(dynamic_cast<TH1 *>(fSumFunc->GetHistogram()->Clone()));
    sum_hist->SetName(GetFuncUniqueName("sum_hist", sum_hist.get()).c_str());
    sum_hist->SetDirectory(0);
    fTempObjects.emplace_back(sum_hist);
    sum_hist->SetLineColor(color);
    sum_hist->SetLineWidth(2);

    // We consider the step to be a background, so we need to remove it from the peaks by zeroing the Parameters. In
    // order not to mess with the parameters of the peaks, we need to make a temporary deep copy of them and work our
    // magic there. These peak functions will be later coverted into histograms and drawn on the canvas.
    fTempPeaks.reserve(fPeaks.size());
    std::for_each(fPeaks.begin(), fPeaks.end(),
                  [this](const auto &peak) { fTempPeaks.emplace_back(std::make_shared<TheuerkaufPeak>(*peak)); });
    std::sort(fTempPeaks.begin(), fTempPeaks.end(),
              [](const auto &a, const auto &b) { return a->GetPos() < b->GetPos(); });

    // int j = 0;
    // for (auto &temp_peak : fTempPeaks)
    // {
    //     TF1 *fcn = temp_peak->GetFunction();
    //     fcn->SetLineColor(GetColor());
    //     fcn->SetLineStyle(4);
    //     fcn->SetNpx(10000.);
    //     j++;
    // }

    if (hist != nullptr)
    {
        std::vector<double> yMax;

        yMax.push_back(this->GetMaximumInRange(tot_bcg_hist, fXMin, fXMax));
        yMax.push_back(this->GetMaximumInRange(_h, fXMin, fXMax));

        _h->GetXaxis()->SetRangeUser(fXMin - 20, fXMax + 20);
        _h->GetYaxis()->SetRangeUser(0, (*std::max_element(yMax.begin(), yMax.end())) * 1.1);
    }

    sum_hist->Draw("SAME");
    tot_bcg_hist->Draw("SAME");
}

#include <TRandom3.h>

double TheuerkaufFitter::GetMinimumInRange(const std::shared_ptr<TH1> h, double x_min, double x_max) noexcept
{
    assert(x_min < x_max && "Minimum range must be valid!");
    assert(h && "Histogram must not be null!");

    double min = std::numeric_limits<double>::max();

    for (auto i = h->FindBin(x_min); i <= h->FindBin(x_max); i++)
    {
        if (h->GetBinContent(i) < min)
        {
            min = h->GetBinContent(i);
        }
    }
    return min;
}

double TheuerkaufFitter::GetMaximumInRange(const std::shared_ptr<TH1> h, double x_min, double x_max) noexcept
{
    assert(x_min < x_max && "Minimum range must be valid!");
    assert(h && "Histogram must not be null!");

    double max = std::numeric_limits<double>::min();

    for (auto i = h->FindBin(x_min); i <= h->FindBin(x_max); i++)
    {
        if (h->GetBinContent(i) > max)
        {
            max = h->GetBinContent(i);
        }
    }
    return max;
}

TCanvas *TheuerkaufFitter::Analyze(TH1 *histAna)
{
    assert(histAna && "Histogram for analysis must not be null!");
    if (!fSumFunc)
        return nullptr;
    fTempObjects.clear();
    fTempPeaks.clear();

    // some function-local things
    const int fcn_npx = fSumFunc->GetNpx();
    const int num_params = GetNumParams();

    // prepare canvas
    std::string canvas_name = GetFuncUniqueName(Form("analyze_canvas_%02.0f", this->GetPeak(0)->GetPos()), this);
    std::string canvas_title = "Peak fit analysis of ";
    canvas_title += this->GetNPeaks() > 1 ? "peaks " : "peak ";
    std::ostringstream oss;
    std::for_each(fPeaks.begin(), fPeaks.end(),
                  [&oss](const auto &peak) { oss << " " << std::fixed << std::setprecision(1) << peak->GetPos(); });
    canvas_title += oss.str();
    auto *can = new TCanvas(canvas_name.c_str(), canvas_title.c_str(), 1400, 800);
    gStyle->SetOptStat(0);
    // can->Divide(1, 2, 0, 0);
    can->ToggleEventStatus();
    can->Divide(1, 2);
    can->Update();
    can->cd(1);

    // prepare main histogram
    std::shared_ptr<TH1> _h(nullptr);
    {
        auto hname = GetFuncUniqueName(Form("_hist_d%f", this->GetPeak(0)->GetPos()), this);
        _h = std::shared_ptr<TH1>(dynamic_cast<TH1 *>(histAna->Clone(hname.c_str())));
        _h->SetDirectory(0);
        // subtracted histo is not resetted here, so we need to clean the list of functions attached to it
        ClearListOfFunctions(_h.get());
        // _h->Draw();
    }

    // create histogram out of a background function
    auto tot_bcg_fcn = std::make_unique<TF1>(GetFuncUniqueName("total_bcg", this).c_str(), this,
                                             &TheuerkaufFitter::EvalTotalBackground, fXMin, fXMax, num_params,
                                             "TheuerkaufFitter", "EvalTotalBackground");
    tot_bcg_fcn->SetNpx(fcn_npx);
    tot_bcg_fcn->SetParameters(fSumFunc->GetParameters());
    std::shared_ptr<TH1> tot_bcg_hist((TH1 *)tot_bcg_fcn->GetHistogram()->Clone());
    tot_bcg_fcn->SetName(GetFuncUniqueName("bcg_hist", tot_bcg_fcn.get()).c_str());
    tot_bcg_hist->SetDirectory(0);
    tot_bcg_hist->SetLineColor(kGray);
    tot_bcg_hist->SetFillStyle(3002);
    tot_bcg_hist->SetFillColor(kGray);

    // create historam out of the sum function (sum of all peaks)
    std::shared_ptr<TH1> sum_fcn_hist((TH1 *)fSumFunc->GetHistogram()->Clone());
    sum_fcn_hist->SetName(GetFuncUniqueName("_sum_fcn_hist", sum_fcn_hist.get()).c_str());
    sum_fcn_hist->SetDirectory(0);
    sum_fcn_hist->SetLineColor(kViolet);
    sum_fcn_hist->SetLineWidth(2);

    // We consider the step to be a background, so we need to remove it from the peaks by zeroing the Parameters. In
    // order not to mess with the parameters of the peaks, we need to make a temporary deep copy of them and work our
    // magic there. These peak functions will be later coverted into histograms and drawn on the canvas.
    fTempPeaks.reserve(fPeaks.size());
    std::for_each(fPeaks.begin(), fPeaks.end(),
                  [this](const auto &peak) { fTempPeaks.emplace_back(std::make_shared<TheuerkaufPeak>(*peak)); });
    std::sort(fTempPeaks.begin(), fTempPeaks.end(),
              [](const auto &a, const auto &b) { return a->GetPos() < b->GetPos(); });
    int j = 0;
    for (auto &temp_peak : fTempPeaks)
    {
        temp_peak->ResetIndexes(true);
        TF1 *fcn = temp_peak->GetFunction();
        fcn->SetParameter(5, 0);
        fcn->SetParameter(6, 0);

        fcn->SetLineColor(GetColor(j + 2));
        // fcn->SetLineStyle(2);
        fcn->SetLineWidth(1);
        fcn->SetNpx(fcn_npx);

        j++;
    }

    // calculate the confidence intervals and the subtracted histogram
    std::shared_ptr<TH1> fit_confidence_histo95((TH1 *)_h->Clone());
    fit_confidence_histo95->SetName(GetFuncUniqueName("fit_confidence", fit_confidence_histo95.get()).c_str());
    fit_confidence_histo95->SetTitle("95% confidence band");
    fit_confidence_histo95->SetDirectory(0);
    fit_confidence_histo95->SetFillColor(kRed);
    fit_confidence_histo95->SetFillStyle(3002);
    fit_confidence_histo95->SetLineWidth(0);
    fit_confidence_histo95->Reset();
    // legend->AddEntry(fit_confidence_histo95.get(), "95% confidence band", "f");
    // this->GetConfidenceIntervals(fit_confidence_histo95.get(), 0.95);
    this->GetConfidenceIntervals(fit_confidence_histo95.get(), 0.95);
    fit_confidence_histo95->ResetStats();

    // set XY ranges for the histogram and the canvas 1
    double yMin = std::numeric_limits<double>::max(), yMax = std::numeric_limits<double>::min();
    double xMinRange, xMaxRange;

    for (auto &temp_peak : fTempPeaks)
    {
        // auto fP = temp_peak->GetFunction();
        if (yMin > fSumFunc->Eval(temp_peak->GetPos()))
            yMin = fSumFunc->Eval(temp_peak->GetPos());
        if (yMax < fSumFunc->Eval(temp_peak->GetPos()))
            yMax = fSumFunc->Eval(temp_peak->GetPos());
    }

    fSumFunc->GetRange(xMinRange, xMaxRange);

    if (yMin > tot_bcg_fcn->Eval(xMinRange))
        yMin = tot_bcg_fcn->Eval(xMinRange);
    if (yMin > tot_bcg_fcn->Eval(xMaxRange))
        yMin = tot_bcg_fcn->Eval(xMaxRange);
    if (yMax < tot_bcg_fcn->Eval(xMinRange))
        yMax = tot_bcg_fcn->Eval(xMinRange);
    if (yMax < tot_bcg_fcn->Eval(xMaxRange))
        yMax = tot_bcg_fcn->Eval(xMaxRange);

    yMax = yMax > this->GetMaximumInRange(_h, fXMin, fXMax) ? yMax : this->GetMaximumInRange(_h, fXMin, fXMax);
    yMin = yMin < this->GetMinimumInRange(_h, fXMin, fXMax) ? yMin : this->GetMinimumInRange(_h, fXMin, fXMax);
    yMin *= 0.9;
    yMax *= 1.1;

    can->cd(1);
    // gPad->Clear();

    if (yMin >= 0 && yMax > 0)
    {
        gPad->SetLogy(true);
        yMin = yMin > 1 ? yMin : 0.1;
    }

    _h->GetXaxis()->SetRangeUser(fXMin, fXMax);
    _h->GetYaxis()->SetRangeUser(yMin, yMax);
    _h->SetTitle(";Energy [keV];Counts");

    // gPad->SetLogy(true);
    _h->Draw("HIST");
    sum_fcn_hist->Draw("SAME");
    tot_bcg_hist->Draw("SAME");
    // fit_confidence_histo95->Draw("SAME E3");
    can->cd(1);

    TLegend *legend1 = new TLegend(0.80, 0.7, 0.99, 0.99);
    legend1->AddEntry(sum_fcn_hist.get(), "Full fit function", "l");
    legend1->AddEntry(tot_bcg_hist.get(), "Background", "f");
    // legend1->AddEntry(fit_confidence_histo95.get(), "95% confidence band", "f");

    // draw function representations
    std::for_each(fTempPeaks.begin(), fTempPeaks.end(), [&](auto &peak) {
        // const auto &fcn = peak->GetFunction()->Draw("SAME");
        const auto &fcn = peak->GetFunction();
        std::shared_ptr<TH1> fcn_hist((TH1 *)fcn->GetHistogram()->Clone());
        fcn_hist->SetName(GetFuncUniqueName(Form("peak_%i", peak->GetID()), fcn_hist.get()).c_str());
        fcn_hist->SetDirectory(0);
        fTempObjects.emplace_back(fcn_hist);

        fcn_hist->Add(tot_bcg_hist.get());

        fcn_hist->SetLineColor(fcn->GetLineColor());
        fcn_hist->SetLineStyle(fcn->GetLineStyle());
        fcn_hist->SetLineWidth(fcn->GetLineWidth());
        fcn_hist->Draw("SAME HIST L");

        std::ostringstream oss;
        oss << "Peak " << std::fixed << std::setprecision(2) << peak->GetPos() << " keV";
        std::string peak_name = oss.str();
        // std::string peak_name = "Peak " + std::to_string(peak->GetPos()) + " keV";
        legend1->AddEntry(fcn_hist.get(), peak_name.c_str(), "l");
    });

    legend1->Draw("SAME");

    /// Switch focus to second pad for the analysis of the fit - prepare histograms for it
    TLegend *legend2 = new TLegend(0.80, 0.85, 0.99, 0.99);
    // create histogram of residuals
    std::shared_ptr<TH1> subtracted_histo((TH1 *)_h->Clone());
    subtracted_histo->SetName(GetFuncUniqueName("subtracted_histo", subtracted_histo.get()).c_str());
    subtracted_histo->SetTitle("Residuals");
    subtracted_histo->SetDirectory(0);
    legend2->AddEntry(subtracted_histo.get(), "Residuals", "l");

    // create 2 sigma interval around the residuals
    std::shared_ptr<TH1> sigma_histo_plus((TH1 *)_h->Clone());
    sigma_histo_plus->SetName(GetFuncUniqueName("_sigma_histo_plus", sigma_histo_plus.get()).c_str());
    sigma_histo_plus->SetTitle("95% statistical uncertainty band (2#sqrt{N})");
    sigma_histo_plus->SetDirectory(0);
    sigma_histo_plus->Reset();
    legend2->AddEntry(sigma_histo_plus.get(), "95% statistical uncertainty band (2#sqrt{N})", "l");

    std::shared_ptr<TH1> sigma_histo_minus((TH1 *)_h->Clone());
    sigma_histo_minus->SetName(GetFuncUniqueName("_sigma_histo_minus", sigma_histo_minus.get()).c_str());
    // sigma_histo_minus->SetTitle("95% statistical error band");
    sigma_histo_minus->SetDirectory(0);
    sigma_histo_minus->Reset();

    // calculate the residuals
    const double normalization = 1. / _h->GetXaxis()->GetBinWidth(1);
    for (int i = _h->FindBin(fXMin); i < _h->FindBin(fXMax); i++)
    {
        double low_edge = subtracted_histo->GetXaxis()->GetBinLowEdge(i);
        double high_edge = subtracted_histo->GetXaxis()->GetBinUpEdge(i);
        double fcn_integral = fSumFunc->Integral(low_edge, high_edge) * normalization;
        subtracted_histo->SetBinContent(i, subtracted_histo->GetBinContent(i) - fcn_integral);
    }

    // set range
    for (int i = _h->FindBin(fXMin); i < _h->FindBin(fXMax); i++)
    {
        sigma_histo_plus->SetBinContent(i, 2. * sqrt(_h->GetBinContent(i)));
        sigma_histo_minus->SetBinContent(i, -2. * sqrt(_h->GetBinContent(i)));
    }

    // calculate residuals over the 95% statistical error band
    int counter_tot = 0;
    int counter_out = 0;
    for (int i = _h->FindBin(fXMin); i < _h->FindBin(fXMax); i++)
    {
        counter_tot++;
        if (subtracted_histo->GetBinContent(i) > sigma_histo_plus->GetBinContent(i) ||
            subtracted_histo->GetBinContent(i) < sigma_histo_minus->GetBinContent(i))
        {
            counter_out++;
        }
    }
    float frac_in = 100 - 100. * (float)counter_out / (float)counter_tot;
    legend2->AddEntry((TObject *)0, Form("(%.1f%%) bins within 95%% band", frac_in), "");

    sigma_histo_plus->SetLineColor(kRed);
    sigma_histo_minus->SetLineColor(kRed);

    subtracted_histo->ResetStats();
    sigma_histo_plus->ResetStats();
    sigma_histo_minus->ResetStats();

    double maxY_sub = subtracted_histo->GetBinContent(subtracted_histo->GetMaximumBin());
    double minY_sub = subtracted_histo->GetBinContent(subtracted_histo->GetMinimumBin());
    double maxY_sig = sigma_histo_plus->GetBinContent(sigma_histo_plus->GetMaximumBin());
    double minY_sig = sigma_histo_minus->GetBinContent(sigma_histo_minus->GetMinimumBin());

    double maxY = maxY_sub > maxY_sig ? maxY_sub * 1.1 : maxY_sig * 1.1;
    double minY = minY_sub < minY_sig ? minY_sub * 1.1 : minY_sig * 1.1;
    // double range = std::fabs(maxY) > std::fabs(minY) ? std::fabs(maxY) : std::fabs(minY);
    subtracted_histo->SetMinimum(minY);
    subtracted_histo->SetMaximum(maxY);
    sigma_histo_minus->SetMinimum(minY);
    sigma_histo_minus->SetMaximum(maxY);
    subtracted_histo->SetMinimum(minY);
    subtracted_histo->SetMaximum(maxY);

    can->cd(2);
    // gPad->SetTopMargin(0.02);
    // gPad->SetBottomMargin(0.2);
    // gPad->SetLeftMargin(0.1);
    // gPad->SetRightMargin(0.01);

    sigma_histo_plus->SetTitle(";Energy (keV); Counts - Fit");
    sigma_histo_plus->GetXaxis()->SetRangeUser(fXMin, fXMax);
    sigma_histo_plus->GetYaxis()->SetLimits(minY, maxY);
    sigma_histo_plus->GetYaxis()->SetRangeUser(minY, maxY);

    sigma_histo_plus->Draw("HIST");
    subtracted_histo->Draw("SAME");
    sigma_histo_minus->Draw("SAME");

    can->SetCrosshair(true);
    can->Update();
    legend2->Draw("SAME");

    // gPad->RedrawAxis();
    // can->RangeChanged();
    // can->ForceUpdate();
    // can->Flush();

    gPad->Update();
    gPad->Modified();
    gPad->Update();

    // store all graphical objects to "private" container that is not managed by ROOT
    fTempObjects.emplace_back(can);
    fTempObjects.emplace_back(_h);
    fTempObjects.emplace_back(tot_bcg_hist);
    fTempObjects.emplace_back(sum_fcn_hist);
    fTempObjects.emplace_back(fit_confidence_histo95);
    fTempObjects.emplace_back(subtracted_histo);
    fTempObjects.emplace_back(sigma_histo_minus);
    fTempObjects.emplace_back(sigma_histo_plus);

    return can;
}

void TheuerkaufFitter::HandleParameterStates()
{
    if (!fSumFunc)
        return;
    // The fSumFunct has complete set of parameters for every TheuerkaufPeak + polynomial
    // background. Each TheuerkaufPeak knows which parameters it should take (based on the
    // peakID) However, some of them are not used by the TheuerkaufPeak function, e.g.
    // when NONE is set, but TheuerkaufFitter does not know about this. We need to fix
    // them otherwise minimizer will vary parameters that are not really used, which will
    // slow down minimization process, produce warnings and who knows what else. The
    // parameters with NONE and FIXED are fixed to its values with
    // fSumFunc->FixParameter() to solve this. If SAME state is observed, all parameters
    // of the same type (e.g. sigma) must point to the same index of the TheuerkaufFitter
    // parameter array, while first being set the FREE and rest fixed. For example we have
    // 3 peaks and linear background, the number of Fitter parameters is 21+2=23. The
    // TheuerkaufPeak internally knows which index to take: peak with ID=0 knows, that it
    // in Eval() it should take indexes no, 0-6 of the fitter, peak with ID=1 indexes 7-13
    // of Fitter etc. Let's say that SAME is set to all sigma parameters. Sigma is third
    // parameter of TheuerkaufPeak, so in Fitter the sigma's will have indexes
    // 2(peakID=0), 9(peakID=1) and 16(peakID=2). So, to apply the "SAME" option, we tell
    // TheuerkaufPeaks with ID 1 and 2, that they should use index 2 for their sigma. The
    // indexes 9 and 16 are set as Fix in the Fitter, thus only index 2 will remain free.
    // This index redirection is done using the SetIndex_Sigma() function of
    // TheuerkaufPeak

    using PAR_STATE = TheuerkaufPeak::ParamState;
    int first_same_index;

    // set volume
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_Volume();
        auto param_index = peak->GetIndex_Volume();
        // NONE is invalid state for volume
        if (state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_Volume(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set POSITION
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_Position();
        auto param_index = peak->GetIndex_Position();
        // NONE and SAME are invalid states for position
        if (state == PAR_STATE::FIXED) // NONE is invalid state for position
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SIGMA
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_Sigma();
        auto param_index = peak->GetIndex_Sigma();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_Sigma(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set TL
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_TailLeft();
        auto param_index = peak->GetIndex_TailLeft();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_TailLeft(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set TR
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_TailRight();
        auto param_index = peak->GetIndex_TailRight();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_TailRight(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SH
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_StepHeight();
        auto param_index = peak->GetIndex_StepHeight();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_StepHeight(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
    // set SW
    first_same_index = -1;
    for (auto &peak : fPeaks)
    {
        auto state = peak->GetState_StepWidth();
        auto param_index = peak->GetIndex_StepWidth();
        if (state == PAR_STATE::NONE || state == PAR_STATE::FIXED)
        {
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
        else if (state == PAR_STATE::SAME)
        {
            // find the first SAME index, keep it free
            if (first_same_index < 0)
            {
                first_same_index = param_index;
                continue;
            }
            // find other SAME parameters, redirect their internal indexes and set to
            // fixed
            peak->SetIndex_StepWidth(first_same_index);
            fSumFunc->FixParameter(param_index, fSumFunc->GetParameter(param_index));
        }
    }
}

void TheuerkaufFitter::DistributeParametersToPeaks()
{
    if (!fSumFunc)
        return;
    for (auto &peak : fPeaks)
    {
        int i_vol = peak->GetIndex_Volume();
        int i_pos = peak->GetIndex_Position();
        int i_sig = peak->GetIndex_Sigma();
        int i_tl = peak->GetIndex_TailLeft();
        int i_tr = peak->GetIndex_TailRight();
        int i_sh = peak->GetIndex_StepHeight();
        int i_sw = peak->GetIndex_StepWidth();

        peak->GetFunction()->SetParameter(0, fSumFunc->GetParameter(i_vol));
        peak->GetFunction()->SetParameter(1, fSumFunc->GetParameter(i_pos));
        peak->GetFunction()->SetParameter(2, fSumFunc->GetParameter(i_sig));
        peak->GetFunction()->SetParameter(3, fSumFunc->GetParameter(i_tl));
        peak->GetFunction()->SetParameter(4, fSumFunc->GetParameter(i_tr));
        peak->GetFunction()->SetParameter(5, fSumFunc->GetParameter(i_sh));
        peak->GetFunction()->SetParameter(6, fSumFunc->GetParameter(i_sw));

        peak->GetFunction()->SetParError(0, fSumFunc->GetParError(i_vol));
        peak->GetFunction()->SetParError(1, fSumFunc->GetParError(i_pos));
        peak->GetFunction()->SetParError(2, fSumFunc->GetParError(i_sig));
        peak->GetFunction()->SetParError(3, fSumFunc->GetParError(i_tl));
        peak->GetFunction()->SetParError(4, fSumFunc->GetParError(i_tr));
        peak->GetFunction()->SetParError(5, fSumFunc->GetParError(i_sh));
        peak->GetFunction()->SetParError(6, fSumFunc->GetParError(i_sw));
        // peak->ResetIndexes(true);
    }

    if (fBcgFunc != nullptr)
    {
        for (unsigned int i = 0; i < fPolyBcgDegree; ++i)
        {
            fBcgFunc->SetParameter(i, fSumFunc->GetParameter(fSumFunc->GetNpar() - fPolyBcgDegree - 1 + i));
            fBcgFunc->SetParError(i, fSumFunc->GetParError(fSumFunc->GetNpar() - fPolyBcgDegree - 1 + i));
        }
    }
}

void TheuerkaufFitter::Fit(TH1 *histFit, std::string options)
{
    std::shared_ptr<TH1> ptr((TH1 *)(histFit->Clone()));
    fTempHist = ptr;
    fTempHist->SetDirectory(0);
    // TH1 *hist = (TH1 *)histFit->Clone();

    assert(fPeaks.size() > 0);

    int num_params = GetNumParams();

    // set binning, becasue peak volume depends on the histogram binning: this is because the bin content represents
    // events that falls within some energy interval, but a function is evaluated at a single point. Technically, the
    // function should evaluate the integral of the function over the bin width...
    double hist_binning = fTempHist->GetXaxis()->GetBinWidth(1);
    std::for_each(fPeaks.begin(), fPeaks.end(), [&hist_binning](auto &peak) { peak->SetBinning(hist_binning); });

    // Create fit function
    fSumFunc = std::make_unique<TF1>(GetFuncUniqueName("fSumFunc", this).c_str(), this, &TheuerkaufFitter::Eval, fXMin,
                                     fXMax, num_params, "TheuerkaufFitter", "Eval");
    int nbins = fTempHist->FindBin(fXMax) - fTempHist->FindBin(fXMin);
    fSumFunc->SetNpx(nbins * 10);
    // set bcg params names
    {
        int polynomial_order = 0;
        for (int i = fPeaks.size() * 7; i < num_params; i++)
        {
            std::string par_name = "bcg_p" + std::to_string(polynomial_order);
            fSumFunc->SetParName(i, par_name.c_str());
            polynomial_order++;
        }
    }

    // *** Initial parameter estimation ***
    int b1 = fTempHist->FindBin(fXMin);
    int b2 = fTempHist->FindBin(fXMax);
    // Check if any of the peaks contain steps
    bool steps =
        std::find_if(fPeaks.begin(), fPeaks.end(), [](const auto &peak) { return peak->HasStep(); }) != fPeaks.end();
    // If there is internal background, we need to estimate it first. We will
    // estimate that the background is constant at the level of the bin with the
    // lowest content if there are no steps, or constant at the level of the
    // leftmost bin in the fit region otherwise (both after substraction of
    // possible external background).
    // NOTE: we generally assume that the step width is positive, so the step
    // function goes to zero at the far left side of the peak. This seems
    // reasonable, as the step width is usually fixed at 1.0.

    double intBg0 = 0.0;
    if (fPolyBcgDegree >= 1)
    {
        if (steps)
        {
            intBg0 = fTempHist->GetBinContent(b1);
            if (fBcgFunc != nullptr)
            {
                intBg0 -= fBcgFunc->Eval(fTempHist->GetBinCenter(b1));
            }
        }
        else
        {
            intBg0 = std::numeric_limits<double>::infinity();

            if (fBcgFunc != nullptr)
            {
                for (int b = b1; b <= b2; ++b)
                {
                    double bc = fTempHist->GetBinContent(b) - fBcgFunc->Eval(fTempHist->GetBinCenter(b));
                    if (bc < intBg0)
                    {
                        intBg0 = bc;
                    }
                }
            }
            else
            {
                for (int b = b1; b <= b2; ++b)
                {
                    double bc = fTempHist->GetBinContent(b);
                    if (bc < intBg0)
                    {
                        intBg0 = bc;
                    }
                }
            }
        }

        // Set background parameters of sum function
        fSumFunc->SetParameter(num_params - fPolyBcgDegree, intBg0);
        if (fPolyBcgDegree >= 2)
        {
            for (int i = num_params - fPolyBcgDegree + 1; i < num_params; ++i)
            {
                fSumFunc->SetParameter(i, 0.0);
            }
        }
    }
    // Next, we must estimate possible steps in the background. We estimate the
    // sum of all step heights as the difference between the first and the last
    // bin in the region (with background substracted). From this, we substract
    // all step heights that have been fixed, and evenly distribute the difference
    // among the others.
    // For the rest of the initial parameter estimation process, we will generally
    // assume the steps to be sharp (step width zero). This is because the step
    // width depends on the peak width, which will only be estimated in the end.
    double avgFreeStep = 0.0;
    if (steps)
    {
        struct Result
        {
            int nStepFree;
            double sumFixedStep;
        };

        auto result =
            std::accumulate(fPeaks.begin(), fPeaks.end(), Result{0, 0.0}, [](Result _result, const auto &peak) {
                if (peak->HasStep())
                {
                    if (peak->GetState_StepHeight() == TheuerkaufPeak::ParamState::FIXED)
                    {
                        _result.sumFixedStep += peak->GetSH();
                    }
                    else
                    {
                        ++_result.nStepFree;
                    }
                }
                return _result;
            });

        double sumStep = fTempHist->GetBinContent(b2) - fTempHist->GetBinContent(b1);
        if (result.nStepFree != 0)
        {
            avgFreeStep = (sumStep - result.sumFixedStep) / result.nStepFree;
        }
    }

    // Estimate peak amplitudes:
    // We assume that the peak positions provided are already a good estimate of
    // the peak centers. The peak amplitude is then estimated as the bin content
    // at the center, with possible external and internal background, and a
    // possible step, substracted. Note that our estimate gets bad if peaks
    // overlap a lot, but it seems hard to do something about this, because we do
    // not know the peak width yet.
    std::vector<double> amps;
    amps.reserve(fPeaks.size());
    double sumAmp = 0.0;

    // First: no steps
    std::transform(fPeaks.begin(), fPeaks.end(), std::back_inserter(amps), [&](const auto &peak) {
        double pos = peak->GetPos();
        double amp = fTempHist->GetBinContent(fTempHist->FindBin(pos)) - intBg0;
        if (fBcgFunc)
        {
            amp -= fBcgFunc->Eval(pos);
        }
        sumAmp += amp;
        return amp;
    });

    // Second: include steps
    using PeakVector_t = std::vector<TheuerkaufPeak>;
    using PeakID_t = PeakVector_t::size_type;

    if (steps)
    {
        // Generate a list of peak IDs sorted by position
        std::vector<PeakID_t> sortedPeakIDs(fPeaks.size());
        std::iota(sortedPeakIDs.begin(), sortedPeakIDs.end(), 0);
        std::sort(sortedPeakIDs.begin(), sortedPeakIDs.end(), [&](const PeakID_t &lhs, const PeakID_t &rhs) {
            return fPeaks[lhs]->GetPos() < fPeaks[rhs]->GetPos();
        });

        struct Sums
        {
            double step, amp;
        };
        auto sums =
            std::accumulate(sortedPeakIDs.begin(), sortedPeakIDs.end(), Sums{0.0, 0.0}, [&](Sums _sums, PeakID_t id) {
                const auto &peak = fPeaks[id];
                double curStep = 0.0;
                if (peak->HasStep())
                {
                    curStep = peak->HasStep() ? peak->GetSH() : avgFreeStep;
                }
                amps[id] -= _sums.step + curStep / 2.0;
                _sums.amp -= _sums.step + curStep / 2.0;
                _sums.step += curStep;
                return _sums;
            });
        sumAmp -= sums.amp;
    }

    // Estimate peak parameters
    //
    // Assuming that all peaks in the fit have the same width, their volume is
    // proportional to their amplitude. We thus calculate the total volume as a
    // sum over bin contents (with background substracted) and distribute it among
    // the peaks according to their amplitude. Last, we can estimate the common
    // width from the total volume and the sum of the amplitudes.
    //
    // NOTE: This assumes that the peaks are purely Gaussian, i.e. that there are
    // no tails.

    // First: calculate total volume
    double sumVol = 0.0;
    for (int b = b1; b <= b2; ++b)
    {
        sumVol += fTempHist->GetBinContent(b);
    }
    sumVol -= intBg0 * (b2 - b1 + 1.);
    if (fBcgFunc != nullptr)
    {
        for (int b = b1; b <= b2; ++b)
        {
            sumVol -= fBcgFunc->Eval(fTempHist->GetBinCenter(b));
        }
    }

    if (steps)
    {
        sumVol -= std::accumulate(fPeaks.begin(), fPeaks.end(), 0.0, [&](double sum, const auto &peak) {
            if (peak->HasStep())
            {
                double curStep = peak->HasStep() ? avgFreeStep : peak->GetSH();
                int b = fTempHist->FindBin(peak->GetPos());
                sum -= curStep * (b2 - std::min(b, b2) + 0.5);
            }
            return sum;
        });
    }

    // normalize on binning
    sumAmp = sumAmp / fTempHist->GetXaxis()->GetBinWidth(1);

    // Second: calculate average peak width (sigma)
    double avgSigma = std::fabs(sumVol / (sumAmp * std::sqrt(2. * M_PI)));
    if (fVerbose > 2)
        std::cout << "sumvVol " << sumVol << " sumAmp " << sumAmp << " " << std::endl;
    if (fVerbose > 2)
        std::cout << "avgSigma " << avgSigma << std::endl;

    // Third: calculate sum of free volumes and amplitudes
    double sumFreeAmp = sumAmp;
    double sumFreeVol = sumVol;
    auto ampIter = amps.begin();
    for (const auto &peak : fPeaks)
    {
        if (peak->GetState_Volume() != TheuerkaufPeak::ParamState::FREE)
        {
            sumFreeAmp -= *(ampIter++);
            sumFreeVol -= peak->GetVol();
        }
    }

    // Init fit parameters for peaks

    int par_n = 0;
    // first initialize all parameters to the their peak values
    for (const auto &peak : fPeaks)
    {
        for (int i = 0; i < 7; i++)
        {
            double _min, _max;
            peak->GetFunction()->GetParLimits(i, _min, _max);
            fSumFunc->SetParLimits(par_n, _min, _max);

            std::string par_name = peak->GetFunction()->GetParName(i);
            switch (peak->GetState(i))
            {
            case TheuerkaufPeak::ParamState::SAME:
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak->GetID() << " is set to \"SAME\""
                              << std::endl;
                par_name += "_SAME ";
                fSumFunc->FixParameter(par_n, peak->GetFunction()->GetParameter(i));
                break;
            case TheuerkaufPeak::ParamState::FIXED:
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak->GetID() << " is set to \"FIXED\""
                              << std::endl;
                par_name += "_FIXED";
                fSumFunc->FixParameter(par_n, peak->GetFunction()->GetParameter(i));
                break;
            case TheuerkaufPeak::ParamState::NONE:
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak->GetID() << " is set to \"NONE\""
                              << std::endl;
                par_name += "_NONE ";
                fSumFunc->FixParameter(par_n, peak->GetFunction()->GetParameter(i));
                break;
            default:
                if (fVerbose > 1)
                    std::cout << "Parameter " << par_name << " of peak " << peak->GetID() << " is set to \"FREE\""
                              << std::endl;
                fSumFunc->SetParameter(par_n, peak->GetFunction()->GetParameter(i));
                break;
            }
            fSumFunc->SetParName(par_n, par_name.c_str());
            par_n++;
        }
    }
    this->HandleParameterStates();

    ampIter = amps.begin();
    for (const auto &peak : fPeaks)
    {
        int i_vol = peak->GetIndex_Volume();
        int i_pos = peak->GetIndex_Position();
        int i_sig = peak->GetIndex_Sigma();
        int i_tl = peak->GetIndex_TailLeft();
        int i_tr = peak->GetIndex_TailRight();
        int i_sw = peak->GetIndex_StepWidth();
        int i_sh = peak->GetIndex_StepHeight();

        if (fVerbose > 1)
            std::cout << "peak " << peak->GetID() << " i_vol " << i_vol << " i_pos " << i_pos << " i_sig " << i_sig
                      << " i_tl " << i_tl << " i_tr " << i_tr << " i_sh " << i_sh << " i_sw " << i_sw << std::endl;

        double amp = *(ampIter++);

        if (fOnlypositivepeaks)
        {
            // The root fit algorithm produces strange results when fitting with extreme
            // limits/boundary conditions, hence set some sane upper limits
            fSumFunc->SetParameter(i_vol, std::max(sumFreeVol * amp / sumFreeAmp, 1.));
            fSumFunc->SetParLimits(i_vol, 0., std::max(100. * sumVol, 0.) + 1e9);
        }
        else
        {
            fSumFunc->SetParameter(i_vol, std::max(sumFreeVol * amp / sumFreeAmp, 1.));
        }

        fSumFunc->SetParameter(i_sig, avgSigma);
        fSumFunc->SetParLimits(i_sig, 0., 10 * (fXMax - fXMin) + 1e3);

        if (peak->HasTL())
        {
            fSumFunc->ReleaseParameter(i_tl);
            fSumFunc->SetParameter(i_tl, avgSigma * 3);
            if (peak->GetState_TailLeft() == TheuerkaufPeak::ParamState::FIXED)
            {
                fSumFunc->FixParameter(i_tl, peak->GetTL());
            }
            // fSumFunc->SetParLimits(i_tl, 0, 1E9);
            // fSumFunc->SetParLimits(i_tl, 0, 1e9);
        }
        if (peak->HasTR())
        {
            fSumFunc->ReleaseParameter(i_tr);
            fSumFunc->SetParameter(i_tr, avgSigma * 3);
            if (peak->GetState_TailRight() == TheuerkaufPeak::ParamState::FIXED)
            {
                fSumFunc->FixParameter(i_tr, peak->GetTR());
            }
            // fSumFunc->SetParLimits(i_tr, 0, 1E9);
            // fSumFunc->SetParameter(i_tr, peak->GetTR());
            // fSumFunc->SetParLimits(i_tr, 0, 100);
        }
        if (peak->HasStep())
        {
            fSumFunc->ReleaseParameter(i_sh);
            fSumFunc->ReleaseParameter(i_sw);

            fSumFunc->SetParameter(i_sw, 1);
            fSumFunc->SetParLimits(i_sw, 0, 5);
        }
    }

    if (fVerbose > 1)
        std::cout << "Sum function has " << fSumFunc->GetNpar() << " from which " << fSumFunc->GetNumberFreeParameters()
                  << " are set free" << std::endl;

    // Now, do the fit
    std::string fit_options = "0NRMIS";
    std::transform(options.begin(), options.end(), options.begin(), ::toupper);
    if (options.find("OUTPUT_ROOT") == std::string::npos)
        fit_options += "Q";
    if ((options.find("LIKELIHOOD") != std::string::npos) || options.find("POISSON") != std::string::npos)
        fit_options += "L";

    // fitter settings
    // TVirtualFitter::SetMaxIterations(1000000);
    // TVirtualFitter::SetMaxFunctionCalls(10000000);
    // TVirtualFitter::SetPrecision(1e-4);
    // TVirtualFitter::SetDefaultFitter("Minuit2");

    fFitResults = fTempHist->Fit(fSumFunc.get(), fit_options.c_str());
    // if fit is not valid or the covariance matrix is weird, lets give it another go
    if (!fFitResults->IsValid() || fFitResults->CovMatrixStatus() != 3)
    {
        fFitResults = fTempHist->Fit(fSumFunc.get(), fit_options.c_str());
    }

    this->DistributeParametersToPeaks();
    // Store Chi^2
    fChiSquare = fSumFunc->GetChisquare();
    if (fVerbose != 0 &&
        (options.find("OUTPUT_NONE") == std::string::npos || options.find("OUTPUT_STANDARD") != std::string::npos))
        this->PrintFitResults();
    if (fVerbose != 0 && options.find("OUTPUT_PLAIN") != std::string::npos)
        this->PrintFitResults_plain();
}

void TheuerkaufFitter::GetRange(double &min, double &max) const noexcept
{
    min = fXMin;
    max = fXMax;
}

#include <Math/OneDimFunctionAdapter.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/RichardsonDerivator.h>
#include <limits>

void TheuerkaufFitter::GetConfidenceIntervals(TH1 *hfit, double cl)

{
    TF1 *f = (TF1 *)fSumFunc.get();
    Int_t npar = f->GetNpar();
    Double_t *grad = new Double_t[npar];
    Double_t *sum_vector = new Double_t[npar];
    // Double_t x[3];
    Double_t x[1];

    Int_t hxfirst = hfit->GetXaxis()->GetFirst();
    Int_t hxlast = hfit->GetXaxis()->GetLast();
    // Int_t hyfirst = hfit->GetYaxis()->GetFirst();
    // Int_t hylast = hfit->GetYaxis()->GetLast();
    // Int_t hzfirst = hfit->GetZaxis()->GetFirst();
    // Int_t hzlast = hfit->GetZaxis()->GetLast();

    TAxis *xaxis = hfit->GetXaxis();
    // TAxis *yaxis = hfit->GetYaxis();
    // TAxis *zaxis = hfit->GetZaxis();
    Double_t t = TMath::StudentQuantile(0.5 + cl / 2, f->GetNDF());
    Double_t chidf = TMath::Sqrt(f->GetChisquare() / f->GetNDF());
    auto covmat = fFitResults->GetCovarianceMatrix();
    // Double_t *matr;
    Double_t c = 0;
    // for (Int_t binz = hzfirst; binz <= hzlast; binz++)
    // {
    //     x[2] = zaxis->GetBinCenter(binz);
    //     for (Int_t biny = hyfirst; biny <= hylast; biny++)
    //     {
    //         x[1] = yaxis->GetBinCenter(biny);
    for (Int_t binx = hxfirst; binx <= hxlast; binx++)
    {
        x[0] = xaxis->GetBinCenter(binx);
        f->GradientPar(x, grad);
        for (Int_t irow = 0; irow < npar; irow++)
        {
            sum_vector[irow] = 0;
            for (Int_t icol = 0; icol < npar; icol++)
                // sum_vector[irow] += matr[irow * npar + icol] * grad[icol];
                sum_vector[irow] += covmat(irow, icol) * grad[icol];
        }
        c = 0;
        for (Int_t i = 0; i < npar; i++)
            c += grad[i] * sum_vector[i];
        c = TMath::Sqrt(c);
        // hfit->SetBinContent(binx, biny, binz, f->EvalPar(x));
        hfit->SetBinContent(binx, f->EvalPar(x));
        // hfit->SetBinError(binx, biny, binz, c * t * chidf);
        hfit->SetBinError(binx, c * t * chidf);
    }
    //     }
    // }
    delete[] grad;
    delete[] sum_vector;
}

void TheuerkaufFitter::GetConfidenceIntervals(unsigned int n, const double *x, double *ci, double cl, bool norm)
{
    unsigned int stride1 = 1;
    unsigned int stride2 = 1;
    // stride1 stride in coordinate  stride2 stride in dimension space
    // i.e. i-th point in k-dimension is x[ stride1 * i + stride2 * k]
    // compute the confidence interval of the fit on the given data points
    // the dimension of the data points must match the dimension of the fit function
    // confidence intervals are returned in array ci

    if (!fSumFunc)
    {
        // check if model function exists
        throw std::runtime_error(
            "FitResult::GetConfidenceIntervals: Cannot compute Confidence Intervals without fit model function");
    }
    assert(fSumFunc);

    // use student quantile in case of normalized errors
    auto chi2 = fSumFunc->GetChisquare();
    auto ndf = fSumFunc->GetNDF();

    double corrFactor = 1;
    if (chi2 <= 0 || ndf == 0)
    {
        norm = false;
    }
    if (norm)
        corrFactor = TMath::StudentQuantile(0.5 + cl / 2, ndf) * std::sqrt(chi2 / ndf);
    else
        // correction to apply to the errors given a CL different than 1 sigma (cl=0.683)
        corrFactor = ROOT::Math::normal_quantile(0.5 + cl / 2, 1);

    unsigned int ndim = fSumFunc->GetNdim();
    unsigned int npar = fSumFunc->GetNpar();

    std::vector<double> xpoint(ndim);
    std::vector<double> grad(npar);
    std::vector<double> vsum(npar);

    auto IsParameterFixed = [&](unsigned int ipar) {
        double pmin, pmax;
        fSumFunc->GetParLimits(ipar, pmin, pmax);
        return pmin == pmax;
    };

    // loop on the points
    for (unsigned int ipoint = 0; ipoint < n; ++ipoint)
    {

        for (unsigned int kdim = 0; kdim < ndim; ++kdim)
        {
            unsigned int i = ipoint * stride1 + kdim * stride2;
            assert(i < ndim * n);
            xpoint[kdim] = x[i];
        }

        // calculate gradient of fitted function w.r.t the parameters
        ROOT::Math::RichardsonDerivator d;
        for (unsigned int ipar = 0; ipar < npar; ++ipar)
        {
            if (!IsParameterFixed(ipar))
            {
                ROOT::Math::OneDimParamFunctionAdapter<const ROOT::Math::IParamMultiFunction &> fadapter(
                    *fFitResults->FittedFunction(), &xpoint.front(), &fFitResults->Parameters().front(), ipar);
                d.SetFunction(fadapter);
                // compute step size as a small fraction of the error
                // (see numerical recipes in C 5.7.8)   1.E-5 is ~ (eps)^1/3
                if (fFitResults->Errors()[ipar] > 0)
                    d.SetStepSize(std::max(fFitResults->Errors()[ipar] * 1.E-5, 1.E-15));
                else
                    d.SetStepSize(std::min(std::max(fFitResults->Parameters()[ipar] * 1.E-5, 1.E-15), 0.0001));

                grad[ipar] = d(fFitResults->Parameters()[ipar]); // evaluate df/dp
            }
            else
                grad[ipar] = 0.; // for fixed parameters
        }

        // multiply covariance matrix with gradient
        vsum.assign(npar, 0.0);
        for (unsigned int ipar = 0; ipar < npar; ++ipar)
        {
            for (unsigned int jpar = 0; jpar < npar; ++jpar)
            {
                vsum[ipar] += fFitResults->GetCovarianceMatrix()(ipar, jpar) * grad[jpar];
            }
        }
        // multiply gradient by vsum
        double r2 = 0;
        for (unsigned int ipar = 0; ipar < npar; ++ipar)
        {
            r2 += grad[ipar] * vsum[ipar];
        }
        double r = std::sqrt(r2);
        ci[ipoint] = r * corrFactor;
    }
}