

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "TApplication.h"
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"

#include "TheuerkaufPeak.h"

Double_t gauss(Double_t *x, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0)
        arg = (x[0] - par[1]) / par[2];
    Double_t val = par[0] * TMath::Exp(-0.5 * arg * arg);
    return val;
}

int main(int argc, char **argv)
{
    TApplication app("app", &argc, argv);

    TH1F *h1 = new TH1F("h1", "histo from a gaussian", 1200, 0, 200);
    TF1 *fcn_gen = new TF1("fcn_gen", "gaus", 0, 200);
    fcn_gen->SetParameter(0, 1);
    fcn_gen->SetParameter(1, 30);
    fcn_gen->SetParameter(2, 5);

    fcn_gen->SetLineColor(kGreen);

    // generate random spectrum
    auto fcn_bcg = std::make_unique<TF1>("fcn_bcg", "pol2", 0, 200);
    fcn_bcg->SetParameter(0, 30);
    fcn_bcg->SetParameter(1, 3);
    fcn_bcg->SetParameter(2, 0.4);
    h1->FillRandom("fcn_bcg", 1000);

    h1->FillRandom("fcn_gen", 10000);
    fcn_gen->SetParameter(1, 70);
    h1->FillRandom("fcn_gen", 10000);
    fcn_gen->SetParameter(1, 50);
    h1->FillRandom("fcn_gen", 3000);
    fcn_gen->SetParameter(1, 150);
    h1->FillRandom("fcn_gen", 3000);

    h1->Sumw2();

    TheuerkaufPeak pk(0, 200, 0, false, false, false);
    pk.SetParameter_Volume(sqrt(2. * M_PI), TheuerkaufPeak::ParamState::FREE, 0, 1E9);
    pk.SetParameter_Position(150, TheuerkaufPeak::ParamState::FREE, 0, 200);
    pk.SetParameter_Sigma(5., TheuerkaufPeak::ParamState::FREE, 0, 10);
    pk.SetParameter_TailLeft(5., TheuerkaufPeak::ParamState::FREE);

    TF1 *tfcn = pk.GetFunction();
    h1->FillRandom(tfcn->GetName(), 10000);

    TheuerkaufFitter fitter(0, 200);
    fitter.SetBackgroundPoly(3);
    fitter.AddPeak(70, false, false, false);
    fitter.AddPeak(30, false, false, false);
    fitter.AddPeak(50, false, false, false);
    fitter.AddPeak(150, false, false, false);

    fitter.Fit(h1);

    for (int i = 0; i < fitter.GetNPeaks(); i++)
    {
        std::cout << "Peak " << i << ": " << std::endl;
        std::cout << "       FWHM: " << fitter.GetPeak(i)->GetFWxM(2.) << std::endl;
        std::cout << "       FWFM: " << fitter.GetPeak(i)->GetFWxM(5.) << std::endl;
        std::cout << "       FWTM: " << fitter.GetPeak(i)->GetFWxM(10.) << std::endl;
    }
    fitter.Analyze(h1);

    app.Run();

    return 0;
}

void test()
{
    int argc;
    char **argv;
    main(argc, argv);
}
