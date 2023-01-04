/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2022 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#include "Toolkit.h"

#include "TRandom3.h"

//--------------------------------------------------------------------------
// ++++++++++++++++++++++ Additional functions used in URANOS ++++++++++++++


//--------------------------------------------------------------------------
// ++++++++++++++++++++++ Fit functions ++++++++++++++++++++++++++++++++++++



//Lorentzian Peak Function with offset and slope
Double_t  lorentzianPeak(Double_t* x, Double_t* par)
{
    return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]) +par[3];
}

Double_t  gaussoffset(Double_t* x, Double_t* par)
{
    //two gaussians
    return (par[0]*( TMath::Exp(-0.5*TMath::Power(((x[0]-par[1])/par[2]),2)) ) +par[3]);
}

Double_t errf( Double_t *x, Double_t *par)
{
  return par[0]*TMath::Erf((x[0]-par[1])/TMath::Sqrt(2)/par[2])+ par[3];
}


// configures the ROOT output for plots

void  rootlogon()
{
    bool largeStyle = false;
    bool atlasStyle = true;

    bool tightMargins = false;
    bool largeMargins = true;

    TStyle *plain = new TStyle("Plain", "Plain Style");
    plain->SetCanvasBorderMode(0);
    plain->SetPadBorderMode(0);
    plain->SetPadColor(0);
    plain->SetCanvasColor(0);
    plain->SetTitleColor(1);
    //plain->SetTitleSize(0.001);
    plain->SetStatColor(0);
    plain->SetPalette(1) ;
    plain->SetAxisColor(1, "X");

    plain->SetGridWidth(0.1);
    plain->SetGridColor(13);
    plain->SetGridStyle(3);

    gROOT->SetStyle( "Plain" );

    gStyle->SetPaperSize(20, 26);//A4
    gStyle->SetPalette(1) ;

    gStyle->SetLineWidth(1.);
    //gStyle->SetHistLineWidth(0.5);

    //plain->SetTextSize(1.1);
    plain->SetLineStyleString(2,"[12 12]");

    plain->SetPadTickX(1);
    plain->SetPadTickY(1);

    gStyle->SetLabelSize(0.025,"xy");
    gStyle->SetTitleSize(0.025,"xy");
    //gStyle->SetTitleOffset(1.,"x");
    //gStyle->SetTitleOffset(1.,"y");
    //gStyle->SetPadTopMargin(0.1);
    //gStyle->SetPadRightMargin(0.1);
    //gStyle->SetPadBottomMargin(0.16);
    //gStyle->SetPadLeftMargin(0.12);
    //gStyle->SetOptFit(0111);
    gStyle->SetOptFit(0000);

    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat("");
    //gStyle->SetOptStat(1111);
    //gStyle->SetOptStat("e");

    if (largeStyle)
    {
        //gStyle->SetLineWidth(3.8);
        //gStyle->SetHistLineWidth(3.8);
        gStyle->SetTextSize(2.);

        gStyle->SetLabelSize(0.04,"xy");
        gStyle->SetTitleSize(0.042,"xy");
        //gStyle->SetTitleOffset(1.4,"x");
        //gStyle->SetTitleOffset(1.4,"y");
        //gStyle->SetTitleOffset(1.7,"y");
    }

    // ATLAS Style
    if (atlasStyle)
    {
        Int_t font=42;
        Double_t tsize=0.025;
        Double_t zsize=0.025;
        plain->SetTextFont(font);
        plain->SetTextSize(tsize);
        plain->SetLabelFont(font,"x");
        plain->SetTitleFont(font,"x");
        plain->SetLabelFont(font,"y");
        plain->SetTitleFont(font,"y");
        plain->SetLabelFont(font,"z");
        plain->SetTitleFont(font,"z");
        plain->SetLabelSize(tsize,"x");
        plain->SetTitleSize(tsize,"x");
        plain->SetLabelSize(tsize,"y");
        plain->SetTitleSize(tsize,"y");
        plain->SetLabelSize(tsize,"z");
        plain->SetTitleSize(tsize,"z");

        plain->SetLegendFont(font);

        //this is really a fucker
        plain->SetFuncWidth(1.5);

        if(tightMargins)
        {
            plain->SetPaperSize(20,26);
            plain->SetPadTopMargin(0.05);
            plain->SetPadRightMargin(0.05);
            plain->SetPadBottomMargin(0.13);
            plain->SetPadLeftMargin(0.12);
        }
        else
        {
            plain->SetPaperSize(20,26);
            plain->SetPadTopMargin(0.05);
            plain->SetPadRightMargin(0.07);
            plain->SetPadBottomMargin(0.10);
            plain->SetPadLeftMargin(0.07);
        }
        if (largeMargins)
        {
            plain->SetPaperSize(20,26);
            plain->SetPadTopMargin(0.06);
            plain->SetPadRightMargin(0.11);
            plain->SetPadBottomMargin(0.11);
            plain->SetPadLeftMargin(0.11);
        }
        plain->SetOptTitle(0);

        Int_t icol=0;
        plain->SetFrameBorderMode(icol);
        plain->SetCanvasBorderMode(icol);
        plain->SetPadBorderMode(icol);
        plain->SetPadColor(icol);
        plain->SetCanvasColor(icol);
        plain->SetStatColor(icol);
    }
    plain->SetBarWidth(3);

    //plain->SetPadTopMargin(0.0);	plain->SetPadRightMargin(0.0);	plain->SetPadBottomMargin(0.0);	plain->SetPadLeftMargin(0.00);
}



//convert colors to RGB
 vector<float> getRGBfromHCL(double h, double c0, double l0)
 {
    vector<float> rgbValues;
    float r,g,b;

    float c = c0*l0;
    float m = l0-c;

    float h0 = h/60.;
    int h0i = floor(h0+0.5);

    float x;

    if (h0 < 2) 		x = c*(1-fabs(h0 - 1));
    else
    {
        if	(h0 < 4)	x = c*(1-fabs(h0-2 - 1));
        else		x = c*(1-fabs(h0-4 - 1));
    }

    if (h0<1)
    {
        r = c; g = x; b = 0;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }
    if (h0<2)
    {
        r = x; g = c; b = 0;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }
    if (h0<3)
    {
        r = 0; g = c; b = x;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }
    if (h0<4)
    {
        r = 0; g = x; b = c;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }
    if (h0<5)
    {
        r = x; g = 0; b = c;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }
    if (h0<6)
    {
        r = c; g = 0; b = x;
        rgbValues.push_back(r+m); rgbValues.push_back(g+m); rgbValues.push_back(b+m);
        return rgbValues;
    }

    rgbValues.push_back(0); rgbValues.push_back(0); rgbValues.push_back(0);
    return rgbValues;
 }

// linear function used to generate color codes
 float getLinearC(double factor, double cmin, double cmax)
 {
     return (cmax-(factor*(cmax-cmin)));
 }

 // power function used to generate color codes
 float getScaledColorValue(double factor, double pow, double cmin, double cmax)
 {
     return (cmax-(TMath::Power(factor,pow)*(cmax-cmin)));
 }

 // set up a linear color table with a number of NCont colors and a number of NRGBs scaled color values in between
 // the scale values are set by the a_i e (0,1)
 // values of C and L, bothe e (0,1) are set to default values
 // typically only the color h e (0,360) has to be varied
void  set_plot_styleSingleGradient(float h, Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    float cmin = 1;
    float cmax = 0.3;
    float lmax = 0.9;
    float lmin = 0.3;
    float pow = 1.5;

    vector<float> stop0 = getRGBfromHCL(h,getScaledColorValue(a0,pow,cmin,cmax),getScaledColorValue(a0,pow,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(h,getScaledColorValue(a1,pow,cmin,cmax),getScaledColorValue(a1,pow,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(h,getScaledColorValue(a2,pow,cmin,cmax),getScaledColorValue(a2,pow,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(h,getScaledColorValue(a3,pow,cmin,cmax),getScaledColorValue(a3,pow,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(h,getScaledColorValue(a4,pow,cmin,cmax),getScaledColorValue(a4,pow,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { stop0.at(0), stop1.at(0), stop2.at(0), stop3.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop0.at(1), stop1.at(1), stop2.at(1), stop3.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop0.at(2), stop1.at(2), stop2.at(2), stop3.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

 // set up a two-color table with a number of NCont colors and a number of NRGBs scaled color values in between
 // the scale values are set by the a_i e (0,1)
 // values of C and L, bothe e (0,1) are set to default values
 // typically only the colors h1 and h2 e (0,360) have to be varied
void  set_plot_styleHeatGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    float cmin = 1;
    float cmax = 0.6;
    float lmax = 0.9;
    float lmin = 0.5;
    float pow1 = 0.4;
    float pow2 = 2.2;
    float h1 = 00;
    float h2 = 60;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { stop0.at(0), stop1.at(0), stop2.at(0), stop3.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop0.at(1), stop1.at(1), stop2.at(1), stop3.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop0.at(2), stop1.at(2), stop2.at(2), stop3.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// different gradient, see definition 'set_plot_styleHeatGradient'
void  set_plot_styleHeatGradient2(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    float cmin = 1;
    float cmax = 0.6;
    float lmax = 0.9;
    float lmin = 0.5;
    float pow1 = 0.4;
    float pow2 = 2.2;
    float h1 = 00;
    float h2 = 60;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { stop4.at(0), stop3.at(0), stop2.at(0), stop1.at(0), stop0.at(0) };
    Double_t green[NRGBs] = { stop4.at(1), stop3.at(1), stop2.at(1), stop1.at(1), stop0.at(1) };
    Double_t blue[NRGBs]  = { stop4.at(2), stop3.at(2), stop2.at(2), stop1.at(2), stop0.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// different gradient, see definition 'set_plot_styleHeatGradient'
void  set_plot_styleHeatGradientModified(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    //a0 = 2*a0; a1 = 2*a1;a2 = 2*a2;
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    float cmin = 1;
    float cmax = 0.7;
    float lmax = 0.95;
    float lmin = 0.5;
    float pow1 = 0.8;
    float pow2 = 2.6;
    float h1 = 05;
    float h2 = 60;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { stop4.at(0), stop2.at(0), stop0.at(0), stop2.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop4.at(1), stop2.at(1), stop0.at(1), stop2.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop4.at(2), stop2.at(2), stop0.at(2), stop2.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// different gradient, see definition 'set_plot_styleHeatGradient'
void  set_plot_styleRainbowGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 9;
    const Int_t NCont = 200;

    float cmin = 0.91;
    float cmax = 0.9;
    float lmax = 0.95;
    float lmin = 0.85;
    float pow1 = 0.4; //pow1 = 1.5;
    float pow2 = 2.2; //pow2 = 1.5;
    float h1 = 0;
    float h2 = 250;

    Double_t a01 = (a0+a1)/2.; Double_t a12 = (a1+a2)/2.; Double_t a23 = (a2+a3)/2.; Double_t a34 = (a3+a4)/2.;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop01 = getRGBfromHCL(getLinearC(a01,h1,h2),getScaledColorValue(a01,pow1,cmin,cmax),getScaledColorValue(a01,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop12 = getRGBfromHCL(getLinearC(a12,h1,h2),getScaledColorValue(a12,pow1,cmin,cmax),getScaledColorValue(a12,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop23 = getRGBfromHCL(getLinearC(a23,h1,h2),getScaledColorValue(a23,pow1,cmin,cmax),getScaledColorValue(a23,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop34 = getRGBfromHCL(getLinearC(a34,h1,h2),getScaledColorValue(a34,pow1,cmin,cmax),getScaledColorValue(a34,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a01, a1, a12, a2, a23, a3, a34, a4 };
    Double_t red[NRGBs]   = { stop0.at(0), stop01.at(0), stop1.at(0), stop12.at(0), stop2.at(0), stop23.at(0), stop3.at(0), stop34.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop0.at(1), stop01.at(1), stop1.at(1), stop12.at(1), stop2.at(1), stop23.at(1), stop3.at(1), stop34.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop0.at(2), stop01.at(2), stop1.at(2), stop12.at(2), stop2.at(2), stop23.at(2), stop3.at(2), stop34.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// different gradient, see definition 'set_plot_styleHeatGradient'
void  set_plot_styleAllGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 9;
    const Int_t NCont = 200;

    float cmin = 0.991;
    float cmax = 0.8;
    float lmax = 0.95;
    float lmin = 0.9;
    float pow1 = 1.5; //pow1 = 1.5;
    float pow2 = 1.5; //pow2 = 1.5;
    float h1 = 0;
    float h2 = 359;

    Double_t a01 = (a0+a1)/2.; Double_t a12 = (a1+a2)/2.; Double_t a23 = (a2+a3)/2.; Double_t a34 = (a3+a4)/2.;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop01 = getRGBfromHCL(getLinearC(a01,h1,h2),getScaledColorValue(a01,pow1,cmin,cmax),getScaledColorValue(a01,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop12 = getRGBfromHCL(getLinearC(a12,h1,h2),getScaledColorValue(a12,pow1,cmin,cmax),getScaledColorValue(a12,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop23 = getRGBfromHCL(getLinearC(a23,h1,h2),getScaledColorValue(a23,pow1,cmin,cmax),getScaledColorValue(a23,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop34 = getRGBfromHCL(getLinearC(a34,h1,h2),getScaledColorValue(a34,pow1,cmin,cmax),getScaledColorValue(a34,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a01, a1, a12, a2, a23, a3, a34, a4 };
    Double_t red[NRGBs]   = { stop0.at(0), stop01.at(0), stop1.at(0), stop12.at(0), stop2.at(0), stop23.at(0), stop3.at(0), stop34.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop0.at(1), stop01.at(1), stop1.at(1), stop12.at(1), stop2.at(1), stop23.at(1), stop3.at(1), stop34.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop0.at(2), stop01.at(2), stop1.at(2), stop12.at(2), stop2.at(2), stop23.at(2), stop3.at(2), stop34.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


//formats the colorbar to a non-inhumanic scheme
void  set_plot_style(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4)
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


// different gradient, see definition 'set_plot_styleHeatGradient'
void  set_plot_styleCool()
{
    float a0 = 0, a1 = 0.56, a2= 0.56, a3 = 0.56, a4 = 1;
    const Int_t NRGBs = 5;
    const Int_t NCont = 200;

    float cmin = 1;
    float cmax = 0.3;
    float lmax = 0.99;
    float lmin = 0.6;
    float pow1 = 0.4;
    float pow2 = 2.2;
    float h1 = 230;
    float h2 = 360;

    vector<float> stop0 = getRGBfromHCL(getLinearC(a0,h1,h2),getScaledColorValue(a0,pow1,cmin,cmax),getScaledColorValue(a0,pow2,lmin,lmax));
    vector<float> stop1 = getRGBfromHCL(getLinearC(a1,h1,h2),getScaledColorValue(a1,pow1,cmin,cmax),getScaledColorValue(a1,pow2,lmin,lmax));
    vector<float> stop2 = getRGBfromHCL(getLinearC(a2,h1,h2),getScaledColorValue(a2,pow1,cmin,cmax),getScaledColorValue(a2,pow2,lmin,lmax));
    vector<float> stop3 = getRGBfromHCL(getLinearC(a3,h1,h2),getScaledColorValue(a3,pow1,cmin,cmax),getScaledColorValue(a3,pow2,lmin,lmax));
    vector<float> stop4 = getRGBfromHCL(getLinearC(a4,h1,h2),getScaledColorValue(a4,pow1,cmin,cmax),getScaledColorValue(a4,pow2,lmin,lmax));

    Double_t stops[NRGBs] = { a0, a1, a2, a3, a4 };
    Double_t red[NRGBs]   = { stop0.at(0), stop1.at(0), stop2.at(0), stop3.at(0), stop4.at(0) };
    Double_t green[NRGBs] = { stop0.at(1), stop1.at(1), stop2.at(1), stop3.at(1), stop4.at(1) };
    Double_t blue[NRGBs]  = { stop0.at(2), stop1.at(2), stop2.at(2), stop3.at(2), stop4.at(2) };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

// modes 0: no luminance scaling, 1: dark at max scale 2: dark at min scale
int getScaledColor(float min, float max, float scale, int mode)
{
    vector<float> rgbValues;
    if (mode==0) rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.4,0.91,0.9),getScaledColorValue(scale,2.2,0.8,0.9));
    else if (mode==1) rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.2,0.99,0.5),getScaledColorValue(scale,1.5,0.1,0.99));
    else rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.2,0.99,0.1),getScaledColorValue(scale,0.9,0.99,0.2));

    return TColor::GetColor(rgbValues.at(0),rgbValues.at(1),rgbValues.at(2));
}

// modes 0: no luminance scaling, 1: dark at max scale 2: dark at min scale
vector<float> getScaledColorRGB(float min, float max, float scale, int mode)
{
    vector<float> rgbValues;
    if (mode==0) rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.4,0.91,0.9),getScaledColorValue(scale,2.2,0.8,0.9));
    else if (mode==1) rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.2,0.99,0.5),getScaledColorValue(scale,1.5,0.1,0.99));
    else rgbValues = getRGBfromHCL(getLinearC(scale,min,max),getScaledColorValue(scale,0.2,0.99,0.1),getScaledColorValue(scale,0.9,0.99,0.2));

    return rgbValues;
}


// fashions up the canvas
void  CanvasFashion(TCanvas* c)
{
    c->SetBorderMode(0);
    c->SetFillColor(0);
}

// fashions up TGraphs
// set axis labels and label plot size to large or small using 'setSize'
void TGraphFashion(TGraph* graph, TString xLabel, TString yLabel, Bool_t setSize)
{
    graph->SetFillColor(19);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1);
    graph->GetXaxis()->SetTitle(xLabel);
    graph->GetYaxis()->SetTitle(yLabel);
    graph->GetXaxis()->SetNoExponent(kTRUE);
    graph->GetYaxis()->SetNoExponent(kTRUE);
    if (setSize)
    {
    graph->GetXaxis()->SetLabelSize(0.027);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetYaxis()->SetLabelSize(0.027);
    graph->GetYaxis()->SetTitleSize(0.03);
    }
}

// fashions up TGraphErrors
// set axis labels and label plot size to large or small using 'setSize'
void TGraphErrorFashion(TGraphErrors* graph, TString xLabel, TString yLabel, Bool_t setSize)
{
    graph->SetFillColor(19);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.0);
    graph->GetXaxis()->SetTitle(xLabel);
    graph->GetYaxis()->SetTitle(yLabel);
    graph->GetXaxis()->SetNoExponent(kTRUE);
    graph->GetYaxis()->SetNoExponent(kTRUE);
    if (setSize)
    {
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.03);
    }
}

// fashions up TGraph2Ds
// set axis labels and label plot size to large or small using 'setSize'
void TGraph2DFashion(TGraph2D* graph, TString xLabel, TString yLabel, TString zLabel, Bool_t setSize)
{
    graph->SetFillColor(19);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1);
    graph->GetXaxis()->SetTitle(xLabel);
    graph->GetYaxis()->SetTitle(yLabel);
    graph->GetZaxis()->SetTitle(zLabel);
    graph->GetXaxis()->SetNoExponent(kTRUE);
    graph->GetZaxis()->SetNoExponent(kTRUE);
    graph->GetYaxis()->SetNoExponent(kTRUE);
    if (setSize)
    {
    graph->GetXaxis()->SetLabelSize(0.027);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetYaxis()->SetLabelSize(0.027);
    graph->GetYaxis()->SetTitleSize(0.03);
    graph->GetZaxis()->SetLabelSize(0.027);
    graph->GetZaxis()->SetTitleSize(0.03);
    }
}

// fashions up TMultiGraphs
// set axis labels and label plot size to large or small using 'setSize'
void TMultiGraphFashion(TMultiGraph* graph, TString xLabel, TString yLabel, Bool_t setSize)
{
    graph->GetXaxis()->SetTitle(xLabel);
    graph->GetYaxis()->SetTitle(yLabel);
    graph->GetXaxis()->SetNoExponent(kTRUE);
    graph->GetYaxis()->SetNoExponent(kTRUE);
    if (setSize)
    {
    graph->GetXaxis()->SetLabelSize(0.027);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetYaxis()->SetLabelSize(0.027);
    graph->GetYaxis()->SetTitleSize(0.03);
    }
}

// fashions up TF1 functions
// set axis labels and label plot size to large or small using 'setSize'
void  TF1Fashion(TF1* hist, TString xLabel, TString yLabel, Bool_t setSize)
{
    if (setSize)
    {
    hist->GetXaxis()->SetLabelSize(0.027);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.027);
    hist->GetYaxis()->SetTitleSize(0.03);
    }
    hist->GetYaxis()->SetTitle(xLabel);
    hist->GetXaxis()->SetTitle(yLabel);
    hist->GetXaxis()->SetNoExponent(kTRUE);
    hist->GetYaxis()->SetNoExponent(kTRUE);
    //hist->SetBarOffset(-0.30);
    //hist->SetBarWidth(0.5);
}

// fashions up TH1 histograms
// set axis labels and label plot size to large or small using 'setSize'
void  TH1Fashion(TH1* hist, TString xLabel, TString yLabel, Bool_t setSize)
{
    if (setSize)
    {
    hist->GetXaxis()->SetLabelSize(0.027);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.027);
    hist->GetYaxis()->SetTitleSize(0.03);
    }
    hist->GetYaxis()->SetTitle(xLabel);
    hist->GetXaxis()->SetTitle(yLabel);
    hist->GetXaxis()->SetNoExponent(kTRUE);
    hist->GetYaxis()->SetNoExponent(kTRUE);
    #ifdef ROOTOLDVERSION
    hist->SetBarOffset(-0.30);
    #endif
    //hist->SetBarWidth(0.5);
}

// fashions up TProfile histograms
// set axis labels and label plot size to large or small using 'setSize'
void  TProfileFashion(TProfile* hist, TString xLabel, TString yLabel, Bool_t setSize)
{
    if (setSize)
    {
    hist->GetXaxis()->SetLabelSize(0.027);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.027);
    hist->GetYaxis()->SetTitleSize(0.03);
    }
    hist->GetYaxis()->SetTitle(xLabel);
    hist->GetXaxis()->SetTitle(yLabel);
    hist->GetXaxis()->SetNoExponent(kTRUE);
    hist->GetYaxis()->SetNoExponent(kTRUE);
    #ifdef ROOTOLDVERSION
    hist->SetBarOffset(-0.30);
    #endif
    //hist->SetBarWidth(0.5);
}


// fashions up TH2 histograms for the pixelwise plots
// set axis labels and label plot size to large or small using 'setSize'
void  TH2Fashion(TH2* hist, Bool_t setSize)
{
    //hist->GetXaxis()->SetRange(0,255);
    //hist->GetYaxis()->SetRange(0,255);
    if (setSize)
    {
        hist->GetXaxis()->SetLabelSize(0.02);
        hist->GetXaxis()->SetTitleSize(0.03);
        hist->GetYaxis()->SetLabelSize(0.02);
        hist->GetYaxis()->SetTitleSize(0.03);
        hist->GetZaxis()->SetLabelSize(0.015);
    }
    if (setSize)
    {
        hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.037);
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetTitleSize(0.037);
        hist->GetYaxis()->SetTitleOffset(1.);
    }
    hist->GetXaxis()->SetTitle("Pixel");
    hist->GetYaxis()->SetTitle("Pixel");
    hist->GetXaxis()->SetNoExponent(kTRUE);
    hist->GetYaxis()->SetNoExponent(kTRUE);
    //gStyle->SetOptStat(0);
}

// fashions up TH1 histograms
// set axis labels and label plot size to large or small using 'setSize'
void  TH2Fashion(TH2* hist, TString xLabel, TString yLabel, Bool_t setSize)
{
    //hist->GetXaxis()->SetRange(0,255);
    //hist->GetYaxis()->SetRange(0,255);
    if (setSize)
    {
        hist->GetXaxis()->SetLabelSize(0.02);
        hist->GetXaxis()->SetTitleSize(0.03);
        hist->GetYaxis()->SetLabelSize(0.02);
        hist->GetYaxis()->SetTitleSize(0.03);
        hist->GetZaxis()->SetLabelSize(0.015);
    }
    else
    {
        hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.045);
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetTitleSize(0.045);
        hist->GetZaxis()->SetLabelSize(0.03);
    }
    if (false)
    {
        hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.037);
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetTitleSize(0.037);
        hist->GetYaxis()->SetTitleOffset(1.);
    }
    hist->GetXaxis()->SetTitle(xLabel);
    hist->GetYaxis()->SetTitle(yLabel);
    hist->GetXaxis()->SetNoExponent(kTRUE);
    hist->GetYaxis()->SetNoExponent(kTRUE);
    //gStyle->SetOptStat(0);
}

// fashions up stacked histograms
// set axis labels and label plot size to large or small using 'setSize'
void  THStackFashion(THStack* stack, TString xLabel, TString yLabel, Bool_t setSize)
{
    //hist->GetXaxis()->SetRange(0,255);
    //hist->GetYaxis()->SetRange(0,255);
    if (setSize)
    {
        stack->GetXaxis()->SetLabelSize(0.03);//originally "2"
        stack->GetXaxis()->SetTitleSize(0.03);
        stack->GetYaxis()->SetLabelSize(0.03);//originally "2"
        stack->GetYaxis()->SetTitleSize(0.03);
    }
    else
    {
        stack->GetXaxis()->SetLabelSize(0.035);
        stack->GetXaxis()->SetTitleSize(0.045);
        stack->GetYaxis()->SetLabelSize(0.035);
        stack->GetYaxis()->SetTitleSize(0.045);
    }
    stack->GetXaxis()->SetTitle(xLabel);
    stack->GetYaxis()->SetTitle(yLabel);
    stack->GetXaxis()->SetNoExponent(kTRUE);
    stack->GetYaxis()->SetNoExponent(kTRUE);
    //gStyle->SetOptStat(0);
}



//--------------------------------------------------------------------------
// ++++++++++++++++++++++ Output/Input +++++++++++++++++++++++++++++++++++++


// plots a TMatrix with z values from zmin to zmax to the specified folder and filename
TH2F* printHisto(Double_t zmin, Double_t zmax, Int_t size, TMatrixF data, TString OutputFolder, TString dataname, Bool_t logPlot)
{
    //cout<<"Generating histograms";

    Float_t pixel;
    TH2F* pixelmatrix = new TH2F("p"+dataname, dataname ,size,0.,size-1,size,0.,size-1);

    for(Int_t  l = 0; l < size; l++)
    {
        for(Int_t  m = 0; m < size; m++)
        {
            pixel = data(l,m);

            pixelmatrix->Fill(l,m, pixel);
        }
    }
    //cout<<"...done"<<endl;

    cout<<"Mean x: "<<pixelmatrix->GetMean(1)<<"+/-"<<pixelmatrix->GetMeanError(1)<<" y: "<<pixelmatrix->GetMean(2)<<"+/-"<<pixelmatrix->GetMeanError(2)<<endl;

    TCanvas* c = new TCanvas("c"+dataname,dataname,1280,1280);
    CanvasFashion(c);
    if (logPlot) c->SetLogz();

    TH2Fashion(pixelmatrix,1);
    pixelmatrix->GetZaxis()->SetRangeUser(zmin,zmax);

    set_plot_style(0.00, 0.20, 0.50, 0.7, 1.00);
    pixelmatrix->Draw("COLZ");

    c->SaveAs(OutputFolder+"/"+dataname+".png");
    c->Print(OutputFolder+"/"+dataname+".root");

    return pixelmatrix;
}


// plots a TH1 histogram with y values from zmin to zmax to the specified folder and filename
void  printTH1F(Double_t zmin, Double_t zmax, TH1F* histo, TString OutputFolder, TString dataname)
{
    TString separator;
    if (OutputFolder.Length() > 1) separator = "/";
    else separator = "";

    TCanvas* c = new TCanvas("c"+dataname,dataname,1280,1280);
    CanvasFashion(c);

    TH1Fashion(histo,"Bins","counts", 1);
    if ((zmin*zmax)!=1) histo->GetYaxis()->SetRangeUser(zmin,zmax);

    histo->Draw("");

    c->SaveAs(OutputFolder+separator+dataname+".png");
    //
    histo->Print(OutputFolder+separator+dataname+".root");
}

// plots a TH2 histogram with z values from zmin to zmax to the specified folder and filename
void  printTH2F(Double_t zmin, Double_t zmax, TH2F* histo, TString OutputFolder, TString dataname)
{
    TString separator;
    if (OutputFolder.Length() > 1) separator = "/";
    else separator = "";

    TCanvas* c = new TCanvas("c"+dataname,dataname,1280,1280);
    CanvasFashion(c);

    //TH2Fashion(histo);
    //if ((zmin!=0)&&(zmax!=0)) histo->GetZaxis()->SetRangeUser(zmin,zmax);

    if ((zmin!=0)&&(zmax!=0)) c->SetLogz();

    set_plot_style(0.00, 0.20, 0.50, 0.7, 1.00);
    histo->Draw("COLZ");

    c->SaveAs(OutputFolder+separator+dataname+".png");
    histo->Print(OutputFolder+separator+dataname+".root");
}

// exports a TH2 histogram as a tab separated ASCII matrix
void  storeTH2ToFile(TString OutputFolder, TString filename, TH2F* th2)
{
    ofstream *stream_out;

    TString separator;
    if (OutputFolder.Length() > 1) separator = "/";
    else separator = "";

    stream_out = new ofstream(OutputFolder+separator+filename,ofstream::out);

    Int_t nX = th2->GetNbinsX();
    Int_t nY = th2->GetNbinsY();

    for(Int_t  i = 0; i < nY; i++)
    {
        for(Int_t  j = 0; j < nX; j++)
        {
            *stream_out<<castDoubleToString(th2->GetBinContent(j,i))<<"\t";
        }
        *stream_out<<endl;
    }
    stream_out->close();
}

// exports an ASCII string to a file
void  storeToFile(TString OutputFolder, TString filename, TString what_has_to_be_written, bool new_line)
{
    ofstream *stream_out;

    TString separator;
    if (OutputFolder.Length() > 1) separator = "/";
    else separator = "";

    stream_out = new ofstream(OutputFolder+separator+filename,ofstream::app);

    *stream_out<<what_has_to_be_written;

    if (new_line) *stream_out<<endl;

    stream_out->close();
}

// embeds the drawing fonts into the vector graphics, which is not done by default by ROOT
// uses a powershell script invoking ghostscript
// should be generalized without hard coded file paths
 void embedFonts(TString path)
 {
     ofstream *stream_out;
     stream_out = new ofstream("fontEmbedscript.ps1",ofstream::out);

     *stream_out<<"$msbuild = \"C:\\Program Files\\gs\\gs9.18\\bin\\gswin64c.exe\""<<endl;
     *stream_out<<"$arguments = \"-q -dNOPAUSE -dBATCH -dPDFSETTINGS=/printer  -sDEVICE=pdfwrite  -sOutputFile="<<path<<".bak "<<path<<"\""<<endl;
     *stream_out<<"start-process $msbuild $arguments -Wait"<<endl;
     *stream_out<<"remove-item "<<path<<endl;
     *stream_out<<"rename-item "<<path<<".bak -newname "<<path<<endl;
     stream_out->close();

     system("powershell.exe -ExecutionPolicy Bypass -File fontEmbedscript.ps1");

     // $arguments = "-q -dNOPAUSE -dBATCH -dPDFSETTINGS=/prepress  -sDEVICE=pdfwrite  -sOutputFile=G:\CASCADE\Analyze\Heidi\Efficiency\calc\heidiEfficiencyNew3b.pdf G:\CASCADE\Analyze\Heidi\Efficiency\calc\heidiEfficiencyNew3.pdf"
     //start-process $msbuild $arguments
 }


// 16.01.2015
// exports a TH1 histogram to a tab separated ASCII matrix
void  storeTH1ToFile(TString OutputFolder, TString filename,TH1* histo)
{
    ofstream *stream_out;
    stream_out = new ofstream(OutputFolder+"/"+filename,ofstream::out);

    Int_t entries = histo->GetNbinsX();

    for(Int_t  l = 1; l < entries+1; l++)
    {
        *stream_out<<castDoubleToString(histo->GetBinCenter(l))<<"\t";
        *stream_out<<castDoubleToString(histo->GetBinContent(l))<<"\t"<<"0.0"<<"\t"<<castDoubleToString(sqrt(histo->GetBinContent(l)))<<endl;
    }

    stream_out->close();
}

// reads an image from the specified folder location and transforms it to a TMatrix
TMatrixF readMatrixPNG(TString folder, TString filename)
{
    QRgb pixelValue;
    string imageStr = std::string((folder+filename).Data())+".png";
    QImage matrixImage;
    matrixImage.load(QString::fromStdString(imageStr));
    int matrixWidthHere = matrixImage.width();
    int matrixHeightHere = matrixImage.height();

    TMatrixF matr(matrixWidthHere,matrixHeightHere);
    for (int i = 0; i<matrixWidthHere; ++i)
    {
        for (int j = 0; j<matrixHeightHere; ++j)
        {
            pixelValue = matrixImage.pixel(i,matrixHeightHere-1-j);
            matr(i,j) = qGray(pixelValue);   //cout<<matr(i,j);
        }
    }
    return matr;
}

// reads an n x n ASCII matrix from the specified folder
// requires the size n to be known
TMatrixF  readmatrix(TString folder, TString filename, TString filetype, Int_t counter, Int_t size)
{
    TMatrixF matrix(size,size);
    TString line, dname, separator;
    Float_t temp;
    Int_t i = 0;

    if (folder.Length() > 1) separator = "/";
    else separator = "";
    //fill every value from the file into the according matrix
    if (counter<10)	dname = folder+separator+filename+"0"+castIntToString(counter)+"."+filetype;
    else dname = folder+separator+filename+castIntToString(counter)+"."+filetype;

    if (counter == -1) dname = folder+separator+filename+"."+filetype;

    ifstream input_stream(dname,ios::in);
    while ( line.ReadLine(input_stream) )
    {
        istrstream stream(line.Data());

        for(Int_t  j = 0; j < size; j++)
        {            
            stream >> temp;
            matrix(i,j) = temp;            
        }
        i++;
        if(i == size) break;
    }
    input_stream.close();

    return matrix;
}

// counts the files of a specified filetype in the specified folder
Int_t  countFiles(TString folder, TString filename, TString filetype)
{
    cout<<"Total number of files: ";
    Int_t n = 0;
    ifstream f;

    TString separator;
    if (folder.Length() > 1) separator = "/";
    else separator = "";

    f.open(folder+separator+filename+"1"+"."+filetype);

    while (f.good())
    {
        f.close();
        n++;
        //if (n<10) f.open(folder+"/"+filename+"0"+castIntToString(n)+".dat");
        f.open(folder+separator+filename+castIntToString(n)+"."+filetype);
    }
    n--;
    cout<<n<<endl;
    f.close();
    return n;
}

// deprecated
int countEventsInEllipse(TH2F* histo, double radiusX, double x, double radiusY, double y)
{
    int eventcounter = 0;

    //int nX = 127, nY = 127;
    int x0 = x-2*radiusX;
    int y0 = y-2*radiusY;
    int nX = x+2*radiusX;
    int nY = y+2*radiusY;

    if (x0<1) x0 = 1;
    if (y0<1) y0 = 1;
    if (nX>127) nX = 127;
    if (nY>127) nY = 127;

    for(Int_t  i = x0; i <= nX; i++) {
        for(Int_t  j = y0; j <= nY; j++) {
            if ((TMath::Power(((i-x)/radiusX),2)+TMath::Power(((j-y)/radiusY),2))<=1)
            {
                eventcounter+=histo->GetBinContent(i,j);
            }
        }
    }
    return eventcounter;
}


//--------------------------------------------------------------------------
// ++++++++++++++++++++++ Cast functions +++++++++++++++++++++++++++++++++++

// should be converted to template functions

string  castDoubleToString(double number)
{
    TString s;
    s += number;
    return (string)s;
}

string  castDoubleToString(double number,  Int_t characters)
{
    TString s;
    s += number;
    s.Resize(characters);
    return (string)s;
}

string  castFloatToString(Float_t number)
{
    TString s;
    s += number;
    return (string)s;
}

string  castFloatToString(Float_t number, Int_t characters)
{
    TString s;
    s += number;
    s.Resize(characters);
    return (string)s;
}


string  castIntToString(Int_t &number)
{
    TString s;
    s += number;
    return (string)s;
}

string  castLongToString(Long_t &number)
{
    TString s;
    s += number;
    return (string)s;
}

// converts a hex value to base-10 integers
unsigned int heXheX(const char *value)
{
    struct CHexMap
    {
        char chr;
        int value;
    };
    const unsigned int HexMapL = 16;
    CHexMap HexMap[HexMapL] =
    {
        {'0', 0}, {'1', 1},
        {'2', 2}, {'3', 3},
        {'4', 4}, {'5', 5},
        {'6', 6}, {'7', 7},
        {'8', 8}, {'9', 9},
        {'A', 10}, {'B', 11},
        {'C', 12}, {'D', 13},
        {'E', 14}, {'F', 15}
    };
    char *mstr = (strdup(value));
    char *s = mstr;
    unsigned int result = 0;
    if (*s == '0' && *(s + 1) == 'X') s += 2;
    bool firsttime = true;
    while (*s != '\0')
    {
        bool found = false;
        for (int i = 0; i < HexMapL; i++)
        {
            if (toupper(*s) == HexMap[i].chr)
            {
                if (!firsttime) result <<= 4;
                result |= HexMap[i].value;
                found = true;
                break;
            }
        }
        if (!found) break;
        s++;
        firsttime = false;
    }
    free(mstr);
    return result;
}

const std::string intToHex( int i)
{
    std::ostringstream oss;
    oss << std::hex << i;
    return oss.str();
}

// fills zero values in a TMatrix by neighbor extrapolation
// deprecated
void extrapolateZeroValues(TMatrixF &matrix)
{
    Float_t corrSum, corrCounter;

    Int_t ylow = matrix.GetColLwb();
    Int_t ydim = matrix.GetColUpb()+1;
    Int_t xlow = matrix.GetRowLwb();
    Int_t xdim = matrix.GetRowUpb()+1;


    for(Int_t  y = ylow; y < ydim; y++)
    {
        for(Int_t  x = xlow; x < xdim; x++)
        {
            if (matrix(x,y)==0)
            {
                corrCounter = 0;
                corrSum = 0;

                if (x>xlow){
                    if (y>ylow){
                        if (matrix(x-1,y-1)!=0) {corrSum+=matrix(x-1,y-1); corrCounter++;}
                    }
                    if (matrix(x-1,y)!=0) {corrSum+=matrix(x-1,y); corrCounter++;}
                    if (y<ydim-1){
                        if (matrix(x-1,y+1)!=0) {corrSum+=matrix(x-1,y+1); corrCounter++;}
                    }
                }
                if (x<xdim-1){
                    if (y>ylow){
                        if (matrix(x+1,y-1)!=0) {corrSum+=matrix(x+1,y-1); corrCounter++;}
                    }
                    if (matrix(x+1,y)!=0) {corrSum+=matrix(x+1,y); corrCounter++;}
                    if (y<ydim-1){
                        if (matrix(x+1,y+1)!=0) {corrSum+=matrix(x+1,y+1); corrCounter++;}
                    }
                }
                if (y>ylow){
                    if (matrix(x,y-1)!=0) {corrSum+=matrix(x,y-1); corrCounter++;}
                }
                if (y<ydim-1){
                    if (matrix(x,y+1)!=0) {corrSum+=matrix(x,y+1); corrCounter++;}
                }

                matrix(x,y) = TMath::Nint(0.001+corrSum/corrCounter);
            }
        }
    }
}

// scales down an n x n by 1/factor
TMatrixF reduceMatrix(TMatrixF &matrix, Float_t factor)
{
    Int_t ylow = matrix.GetColLwb();
    Int_t ydim = matrix.GetColUpb()+1;
    Int_t xlow = matrix.GetRowLwb();
    Int_t xdim = matrix.GetRowUpb()+1;

    Float_t temp;
    Float_t max=0;

    TMatrixF redMatrix(xdim,ydim); redMatrix=0;

    for(Int_t  y = ylow; y < ydim; y++)
    {
        for(Int_t  x = xlow; x < xdim; x++)
        {
            if (matrix(x,y)>max) max = matrix(x,y);
        }
    }

    if (max!=0)
    {
        for(Int_t  y = ylow; y < ydim; y++)
        {
            for(Int_t  x = xlow; x < xdim; x++)
            {
                temp = matrix(x,y)-factor*max;
                if (temp>0) redMatrix(x,y) = temp;
            }
        }

    }

    return redMatrix;
}

// calculates the gradient matrix of a TMatrix
TMatrixF getGradientMatrix(TMatrixF &matrix)
{
    Int_t ylow = matrix.GetColLwb();
    Int_t ydim = matrix.GetColUpb()+1;
    Int_t xlow = matrix.GetRowLwb();
    Int_t xdim = matrix.GetRowUpb()+1;

    TMatrixF gradMatrix(xdim,ydim); gradMatrix=0;

    for(Int_t  y = ylow+1; y < ydim-1; y++)
    {
        for(Int_t  x = xlow+1; x < xdim-1; x++)
        {
            gradMatrix(x,y) = sqrt((TMath::Power(matrix(x-1,y)-matrix(x+1,y),2)+TMath::Power(matrix(x,y-1)-matrix(x,y+1),2))/9.);
        }
    }

    return gradMatrix;
}

// calculates the gradient matrix of a TMatrix and stores it to *newMatrix
void getGradientMatrixFromTH2(TH2F* matrix, TH2F* newMatrix)
{
    Int_t ylow = 1;
    Int_t ydim = matrix->GetNbinsY();
    Int_t xlow = 1;
    Int_t xdim = matrix->GetNbinsX();

    for(int x=0;x<=xdim;x++)
    {
        for(int y=0;y<=ydim;y++)
        {
            if ((x<4)||(y<4)||(x>xdim-4)||(y>ydim-4))newMatrix->SetBinContent(x,y,0);

        }
    }

    float gradientX = 0;
    float gradientY = 0;
    float gradientD1 = 0;
    float gradientD2 = 0;

    //TF1* gradientLine = new TF1("gradientLine","pol1",0,10);
    TF1* gradientLine = new TF1("gradientLine","[0]+[1]*x",-1,5);

    int px, p0,p1,p2,p3,p4,p5;
    TGraphErrors* gradientGraph = new TGraphErrors(6);

    for(Int_t  y = ylow+3; y < ydim-3; y++) {

        for(Int_t  x = xlow+3; x < xdim-3; x++) {

            px = matrix->GetBinContent(x,y);

            if(px==0)
            {
                gradientX = 0;
                gradientY = 0;
            }
            else
            {
                p0 = matrix->GetBinContent(x-3,y);
                p1 = matrix->GetBinContent(x-2,y);
                p2 = matrix->GetBinContent(x-1,y);
                p3 = matrix->GetBinContent(x+1,y);
                p4 = matrix->GetBinContent(x+2,y);
                p5 = matrix->GetBinContent(x+3,y);

                if( (p0+p1+p2+p3) > 0)
                {
                    gradientGraph->SetPoint(0,0,p0);  gradientGraph->SetPointError(0,0,TMath::Sqrt(p0));
                    gradientGraph->SetPoint(1,1,p1);  gradientGraph->SetPointError(1,0,TMath::Sqrt(p1));
                    gradientGraph->SetPoint(2,2,p2);  gradientGraph->SetPointError(2,0,TMath::Sqrt(p2));
                    gradientGraph->SetPoint(3,3,p3);  gradientGraph->SetPointError(3,0,TMath::Sqrt(p3));
                    gradientGraph->SetPoint(4,4,p4);  gradientGraph->SetPointError(4,0,TMath::Sqrt(p4));
                    gradientGraph->SetPoint(5,5,p5);  gradientGraph->SetPointError(5,0,TMath::Sqrt(p5));

                    gradientLine->SetParameters(1,0);
                    gradientGraph->Fit("gradientLine","RQ");
                    gradientX = gradientLine->GetParameter(1);
                }
                else gradientX = 0;

                p0 = matrix->GetBinContent(x,y-3);
                p1 = matrix->GetBinContent(x,y-2);
                p2 = matrix->GetBinContent(x,y-1);
                p3 = matrix->GetBinContent(x,y+1);
                p4 = matrix->GetBinContent(x,y+2);
                p5 = matrix->GetBinContent(x,y+3);

                if( (p0+p1+p2+p3) > 0)
                {
                    gradientGraph->SetPoint(0,0,p0);  gradientGraph->SetPointError(0,0,TMath::Sqrt(p0));
                    gradientGraph->SetPoint(1,1,p1);  gradientGraph->SetPointError(1,0,TMath::Sqrt(p1));
                    gradientGraph->SetPoint(2,2,p2);  gradientGraph->SetPointError(2,0,TMath::Sqrt(p2));
                    gradientGraph->SetPoint(3,3,p3);  gradientGraph->SetPointError(3,0,TMath::Sqrt(p3));
                    gradientGraph->SetPoint(4,4,p4);  gradientGraph->SetPointError(4,0,TMath::Sqrt(p4));
                    gradientGraph->SetPoint(5,5,p5);  gradientGraph->SetPointError(5,0,TMath::Sqrt(p5));

                    gradientLine->SetParameters(1,0);
                    gradientGraph->Fit("gradientLine","RQ");
                    gradientY = gradientLine->GetParameter(1);
                }
                else gradientY = 0;

                p0 = matrix->GetBinContent(x-2,y-2);
                p1 = matrix->GetBinContent(x-1,y-1);
                p2 = matrix->GetBinContent(x+1,y+1);
                p3 = matrix->GetBinContent(x+2,y+2);

                if( (p0+p1+p2+p3) < 0)
                {
                    gradientGraph->SetPoint(0,0,p0);  gradientGraph->SetPointError(0,0,TMath::Sqrt(p0));
                    gradientGraph->SetPoint(1,1,p1);  gradientGraph->SetPointError(1,0,TMath::Sqrt(p1));
                    gradientGraph->SetPoint(2,2,p2);  gradientGraph->SetPointError(2,0,TMath::Sqrt(p2));
                    gradientGraph->SetPoint(3,3,p3);  gradientGraph->SetPointError(3,0,TMath::Sqrt(p3));

                    gradientLine->SetParameters(1,0);
                    gradientGraph->Fit("gradientLine","RQ");
                    gradientD1 = gradientLine->GetParameter(1);
                }
                else gradientD1 = 0;

                p0 = matrix->GetBinContent(x+2,y-2);
                p1 = matrix->GetBinContent(x+1,y-1);
                p2 = matrix->GetBinContent(x-1,y+1);
                p3 = matrix->GetBinContent(x-2,y+2);

                if( (p0+p1+p2+p3) < 0)
                {
                    gradientGraph->SetPoint(0,0,p0);  gradientGraph->SetPointError(0,0,TMath::Sqrt(p0));
                    gradientGraph->SetPoint(1,1,p1);  gradientGraph->SetPointError(1,0,TMath::Sqrt(p1));
                    gradientGraph->SetPoint(2,2,p2);  gradientGraph->SetPointError(2,0,TMath::Sqrt(p2));
                    gradientGraph->SetPoint(3,3,p3);  gradientGraph->SetPointError(3,0,TMath::Sqrt(p3));

                    gradientLine->SetParameters(1,0);
                    gradientGraph->Fit("gradientLine","RQ");
                    gradientD2 = gradientLine->GetParameter(1);
                }
                else gradientD2 = 0;
            }
            //newMatrix->SetBinContent(x,y, sqrt((TMath::Power(matrix->GetBinContent(x-1,y)-matrix->GetBinContent(x+1,y),2)+TMath::Power(matrix->GetBinContent(x,y-1)-matrix->GetBinContent(x,y+1),2))/9.) );

            if (gradientX>50) gradientX = 0;
            if (gradientY>50) gradientY = 0;
            if (gradientD1>50) gradientD1 = 0;
            if (gradientD2>50) gradientD2 = 0;
            if (gradientX<-50) gradientX = 0;
            if (gradientY<-50) gradientY = 0;
            if (gradientD1<-50) gradientD1 = 0;
            if (gradientD2<-50) gradientD2 = 0;

            newMatrix->SetBinContent(x,y, sqrt(TMath::Power(gradientX,2)+TMath::Power(gradientY,2)+TMath::Power(gradientD1,2)+TMath::Power(gradientD2,2)) );
        }
    }
    delete gradientLine;
    delete gradientGraph;
}


// converts an ASCII matrix to a TH2 histogram
TH2F* getTH2fromMatrix(TMatrixF &matrix, TString thName, TString thTitle)
{
    Int_t ylow = matrix.GetColLwb();
    Int_t ydim = matrix.GetColUpb()+1;
    Int_t xlow = matrix.GetRowLwb();
    Int_t xdim = matrix.GetRowUpb()+1;

    TH2F* matrixTH = new TH2F(thName, thTitle , xdim-xlow, xlow-0.5 , xdim-1+0.5, ydim-ylow, ylow-0.5 , ydim-1+0.5);

    for(Int_t  y = ylow; y < ydim; y++)
    {
        for(Int_t  x = xlow; x < xdim; x++)
        {
            matrixTH->Fill(x,y,matrix(x,y));
        }
    }

    return matrixTH;
}

// deprecated
TH1F* convoluteGaussian(TH1F* hist, Double_t gaussRelW)
{
    Float_t mDiv = TMath::Abs(hist->GetBinCenter(2) - hist->GetBinCenter(1));
    Int_t nMaxBins =  hist->GetNbinsX();

    Double_t binSum;

    TH1F* convolvedHist = new TH1F(hist->GetName(), hist->GetTitle(), nMaxBins , hist->GetBinCenter(0), hist->GetBinCenter(hist->GetNbinsX()-1));

    for(Int_t  n = 1; n < nMaxBins; n++)
    {
        //if (hist->GetBinContent(n)==0) continue;

        binSum = 0;

        for(Int_t  m = -100; m < 101; m++)
        {
            if ((n-m > 0)&&(n-m < nMaxBins))
            {
                binSum += TMath::Gaus(m*mDiv, 0, gaussRelW/2.35482*hist->GetBinCenter(n), true) * hist->GetBinContent(n-m);
            }
        }
        convolvedHist->SetBinContent(n, binSum);
    }
    return convolvedHist;
}

// removes all histograms from the pointer vector from the memory.
void histoClearUp(vector< TObject* >* vectorTH)
{
    TString name;
    vector<TString> names;

    for(int k=0; k<vectorTH->size();k++)
    {
        //cout<<vectorTH->at(k)->GetName()<<endl;
        if (vectorTH->at(k)->GetName()) names.push_back(vectorTH->at(k)->GetName());
        //cout<<names.at(k)<<" ";
    }

    int length = names.size();

    for(int k=0; k < length;k++)
    {
        //cout<<names.at(k)<<" ";
        if (gDirectory->FindObject(names.at(k))) {delete gDirectory->FindObject(names.at(k));}
    }
    vectorTH->clear();
    cout<<"deleting histos..."<<endl;
}




/**
 *  calculates the effective probability for an interaction using
 *  material density 'density' in g/cm3, atomic weight 'atWeight' in u
 *  single interaction cross section 'cs' in barns, and
 *  the wavelength lambda in Angstroms, which is only relevant for absorption, otherwise set to 0
 *
 *
 * @param density
 * @param atWeight
 * @param cs
 * @param lambda
 * @return wwProb
 */
double getWWProb(double density, double atWeight, double cs, double lambda)
{
    double wwProb = 0;

    if (lambda == 0)	wwProb = 0.06022045 * density / atWeight * cs;
    else				wwProb = 0.06022045 * density / atWeight * cs * lambda / 1.8;

    return wwProb;
}

/**
 * calculates the wavelength lambda in Angstroms for a given energy in MeV
 *
 * @param energy
 * @return the wavelength lambda.
 */
double getLfromE(double energy)
{
    return sqrt(81.81 / energy / 1e9);
}


/**
 * calculates the energy in MeV for a given wavelength lambda in Angstroms
 *
 * @param lambda -the wavelength.
 * @return the energy .
 */
double getEfromL(double lambda)
{
    return 81.81 / TMath::Power(lambda, 2) / 1e9;
}


/**
 * Create a Spline Function for an energy (log) response function.
 *
 * @return spline energy.
 */
TSpline3* getSplinedDetectorEnergyModel()
{
    const int nData = 40;

    float xVal[nData] = { 1.02E-08, 1.82E-08, 3.24E-08, 5.75E-08, 1.02E-07,   1.82E-07,   3.24E-07,   5.75E-07,   1.02E-06,   1.82E-06,   3.24E-06,   5.75E-06,   1.02E-05,   1.82E-05,   3.24E-05,   5.75E-05,   0.000102329,   0.00018197,   0.000323594,   0.00057544,   0.001023293,   0.001819701,   0.003235937,   0.0057544,   0.01023293,   0.018197009,   0.032359365,   0.057543993,   0.102329299,   0.18197009,   0.323593646,   0.57543993,   1.023293018,   1.819700837,   3.235936642,   5.7543993,   10.23293018,   18.19700813,  32.35936737,   57.5439949 };
    float yVal[nData] = { 6.9625 ,7.351 ,7.744 ,8.3515 ,10.429 ,13.8135 ,16.605 ,19.336 ,20.826 ,21.8715 ,22.334 ,22.633 ,22.4995 ,22.208 ,21.594 ,20.939 ,20.1395 ,19.3815 ,18.54 ,17.579 ,16.933 ,16.04 ,15.0365 ,14.2945 ,13.6655 ,12.743 ,11.9555 ,11.3265 ,10.237 ,9.27 ,8.075 ,6.5935 ,5.193 ,3.765 ,2.5505 ,1.6125 ,1.031 ,0.422 ,0.4885 ,0.3255 };

    for (int i = 0; i < nData; i++)
    {
        xVal[i] = TMath::Log10(xVal[i]); //probably easier to spline?
        yVal[i] = yVal[i] * 4.44 / 100.; //scale to 1 maximum value
    }

    TGraph* graphEnergy = new TGraph(nData, xVal, yVal);

    TSpline3* splineEnergy = new TSpline3("splineEnergyModel", graphEnergy);

    delete graphEnergy;

    return splineEnergy;
}

/**
 * Get the number of lines in a File.
 *
 * @param dname - file name
 * @return lengthFile - length of the file..
 */
int getNumberOfFileEntries(TString dname)
{
    int lengthFile = 0;
    TString curLine;

    ifstream input_streamCheck(dname, ios::in);

    while (curLine.ReadLine(input_streamCheck)) { lengthFile++; }

    input_streamCheck.close();

    return lengthFile;
}

// Legendre Polynomials
// inline for computational speed
inline double P0(double x) { return 1.0; }

// n = 1
inline double P1(double x) { return x; }
// n = 2
inline double P2(double x) { return ((3. * x * x) - 1.) * 0.5; }

inline double P7(double x) { return ((429. * x * x * x * x * x * x * x - 693. * x * x * x * x * x + 315. * x * x * x - 35. * x) * 0.0625); }

inline double P8(double x) { return ((6435. * x * x * x * x * x * x * x * x - 12012. * x * x * x * x * x * x + 6930. * x * x * x * x - 1260. * x * x + 35.) * 0.0078125); }

inline double P9(double x) { return ((12155. * x * x * x * x * x * x * x * x * x - 25740. * x * x * x * x * x * x * x + 18018. * x * x * x * x * x - 4620. * x * x * x + 315. * x) * 0.0078125); }

inline double P10(double x) { return ((46189. * x * x * x * x * x * x * x * x * x * x - 109395. * x * x * x * x * x * x * x * x + 90090. * x * x * x * x * x * x - 30030. * x * x * x * x + 3465. * x * x - 63.) * 0.00390625); }

inline double P11(double x) { return ((88179. * x * x * x * x * x * x * x * x * x * x * x - 230945. * x * x * x * x * x * x * x * x * x + 218790. * x * x * x * x * x * x * x - 90090. * x * x * x * x * x + 15015. * x * x * x - 693. * x) * 0.00390625); }

inline double P12(double x) { return ((676039. * x * x * x * x * x * x * x * x * x * x * x * x - 1939938. * x * x * x * x * x * x * x * x * x * x + 2078505. * x * x * x * x * x * x * x * x - 1021020. * x * x * x * x * x * x + 225225. * x * x * x * x - 18018. * x * x + 231.) * 0.0009765625); }

inline double P13(double x) { return ((1300075. * x * x * x * x * x * x * x * x * x * x * x * x * x - 4056234. * x * x * x * x * x * x * x * x * x * x * x + 4849845. * x * x * x * x * x * x * x * x * x - 2771340. * x * x * x * x * x * x * x + 765765. * x * x * x * x * x - 90090. * x * x * x + 3003. * x) * 0.0009765625); }

inline double P14(double x) { return ((5014575. * x * x * x * x * x * x * x * x * x * x * x * x * x * x - 16900975. * x * x * x * x * x * x * x * x * x * x * x * x + 22309287. * x * x * x * x * x * x * x * x * x * x - 14549535. * x * x * x * x * x * x * x * x + 4849845. * x * x * x * x * x * x - 765765. * x * x * x * x + 45045. * x * x - 429.) * 0.00048828125); }

inline double P15(double x) { return ((9694845. * x * x * x * x * x * x * x * x * x * x * x * x * x * x * x - 35102025. * x * x * x * x * x * x * x * x * x * x * x * x * x + 50702925. * x * x * x * x * x * x * x * x * x * x * x - 37182145. * x * x * x * x * x * x * x * x * x + 14549535. * x * x * x * x * x * x * x - 2909907. * x * x * x * x * x + 255255. * x * x * x - 6435. * x) * 0.00048828125); }

inline double P16(double x) { return ((300540195. * x * x * x * x * x * x * x * x * x * x * x * x * x * x * x * x - 1163381400. * x * x * x * x * x * x * x * x * x * x * x * x * x * x + 1825305300. * x * x * x * x * x * x * x * x * x * x * x * x - 1487285800. * x * x * x * x * x * x * x * x * x * x + 669278610. * x * x * x * x * x * x * x * x - 162954792. * x * x * x * x * x * x + 19399380. * x * x * x * x - 875160. * x * x + 6435.) * 0.000030517578125); }

inline double P17(double x)  { return ((583401555.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-2404321560.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x +4071834900.*x*x*x*x*x*x*x*x*x*x*x*x*x-3650610600.*x*x*x*x*x*x*x*x*x*x*x+1859107250.*x*x*x*x*x*x*x*x*x-535422888.*x*x*x*x*x*x*x + 81477396.*x*x*x*x*x-5542680.*x*x*x+109395.*x) *0.000030517578125) ; }

inline double P18(double x)  { return ((2268783825.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-9917826435.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x +18032411700.*x*x*x*x*x*x*x*x*x*x*x*x*x*x-17644617900.*x*x*x*x*x*x*x*x*x*x*x*x+10039179150.*x*x*x*x*x*x*x*x*x*x-3346393050.*x*x*x*x*x*x*x*x + 624660036.*x*x*x*x*x*x-58198140.*x*x*x*x+2078505.*x*x-12155.) *0.000015258789) ; }

inline double P19(double x)  { return ((4418157975.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-20419054425.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x +39671305740.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-42075627300.*x*x*x*x*x*x*x*x*x*x*x*x*x+26466926850.*x*x*x*x*x*x*x*x*x*x*x-10039179150.*x*x*x*x*x*x*x*x*x + 2230928700.*x*x*x*x*x*x*x-267711444.*x*x*x*x*x+14549535.*x*x*x-230945.*x) *0.000015258789) ; }

inline double P20(double x)  { return ((34461632205.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-167890003050.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x +347123925225.*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x-396713057400.*x*x*x*x*x*x*x*x*x*x*x*x*x*x+273491577450.*x*x*x*x*x*x*x*x*x*x*x*x-116454478140.*x*x*x*x*x*x*x*x*x*x + 30117537450.*x*x*x*x*x*x*x*x-4461857400.*x*x*x*x*x*x+334639305.*x*x*x*x-9699690.*x*x+46189.) *0.00000381469727) ; }


/**
 * Defining an legrende polynomial function
 * @param n - degree
 * @param x - x value
 * @return evaluated result of the polynomial
 */
inline double legendre_Pl(int n, double x)
{
    switch (n)
    {
    case 0: return P0(x);
        break;
    case 1: return P1(x);
        break;
    case 2: return P2(x);
        break;
    case 3:  return 0.5 * (5. * x * x * x - 3. * x);
        break;
    case 4: return 0.125 * (35. * x * x * x * x - 30. * x * x + 3.);
        break;
    case 5: return 0.125 * (63. * x * x * x * x * x - 70. * x * x * x + 15. * x);
        break;
    case 6: return 0.0625 * (231. * x * x * x * x * x * x - 315. * x * x * x * x + 105. * x * x + 5.);
        break;
    case 7: return P7(x);
        break;
    case 8: return P8(x);
        break;
    case 9: return P9(x);
        break;
    case 10: return P10(x);
        break;
    case 11: return P11(x);
        break;
    case 12: return P12(x);
        break;
    case 13: return P13(x);
        break;
    case 14: return P14(x);
        break;
    case 15: return P15(x);
        break;
    case 16: return P16(x);
        break;
    case 17: return P17(x);
        break;
    case 18: return P18(x);
        break;
    case 19: return P19(x);
        break;
    case 20: return P20(x);
        break;
    }

    if (x == 1.0)
    {
        return 1.0;
    }

    if (x == -1.0)
    {
        return ((n % 2 == 0) ? 1.0 : -1.0);
    }

    if ((x == 0.0) && (n % 2))
    {
        return 0.0;
    }

    double pnm1(P20(x)) ;
    double pnm2(P19(x)) ;
    double pn(pnm1);

    for (int l = 21 ; l <= n ; l++)
    {
        pn = (((2.0 * (double)l) - 1.0) * x * pnm1 -
            (((double)l - 1.0) * pnm2)) / (double)l;
        pnm2 = pnm1;
        pnm1 = pn;
    }

    return pn;
}



/**
 * replaces the x-axis with a logarithmic x-axis of 100 bins
 * root has no logarithmic binning, so it requires generating a bin array
 * @param h
 *
 */
void logXaxis(TH2* h)
{
    TAxis* axis = h->GetXaxis();
    const Int_t nbins = 100;
    Double_t xmin = axis->GetXmin();
    Double_t xmax = axis->GetXmax();
    Double_t logXmin = TMath::Log10(xmin);
    Double_t logXmax= TMath::Log10(xmax);
    Double_t binwidth = (logXmax - logXmin) / nbins;
    Double_t xBins[nbins + 1];
    xBins[0] = xmin;
    for (Int_t i = 1; i <= nbins; i++)
    {
        xBins[i] = xmin + TMath::Power(10, logXmin + i * binwidth);
    }

    axis->Set(nbins, xBins);
}


/**
 * replaces the x-axis with a logarithmic x-axis of 1000 bins
 * root has no logarithmic binning, so it requires generating a bin array
 * @param h
 *
 */
void logaxis(TH1* h)
{
    TAxis* axis = h->GetXaxis();
    const Int_t nBins = 1000;
    Double_t xMin = axis->GetXmin();
    Double_t xMax = axis->GetXmax();
    Double_t logXmin = TMath::Log10(xMin);
    Double_t logXmax= TMath::Log10(xMax);
    Double_t binWidth = (logXmax- logXmin) / nBins;
    Double_t xBins[nBins + 1];
    xBins[0] = xMin;
    for (Int_t binsItr = 1; binsItr <= nBins; binsItr++)
    {
        xBins[binsItr] = xMin + TMath::Power(10, logXmin + binsItr * binWidth);
    }

    axis->Set(nBins, xBins);
}


/**
 * rebins a logarithmic x-axis from 1000 bins to 100
 * @param h
 *
 */
void rebinX(TH1* h)
{
    const Int_t nBinsPrev = 1000;
    const Int_t nBins = 100;
    TAxis* axis = h->GetXaxis();
    Double_t xMin = axis->GetXmin();
    Double_t xMax = axis->GetXmax();

    Double_t logXmin = TMath::Log10(xMin);
    Double_t logXmax= TMath::Log10(xMax);
    Double_t binwidth = (logXmax- logXmin) / nBins;
    Double_t xBins[nBins + 1];
    xBins[0] = xMin;

    for (Int_t i = 1; i <= nBins; i++)
    {
        xBins[i] = xMin + TMath::Power(10, logXmin + i * binwidth);
    }
    Double_t binContent[nBins + 1];
    for (Int_t i = 1; i <= nBins; i++) { binContent[i] = 0; }

    int binNo;

    for (Int_t i = 1; i <= nBinsPrev + 1; i++)
    {
        binNo = i / 10;
        binContent[binNo] = binContent[binNo] + h->GetBinContent(i);
    }

    for (Int_t i = 1; i <= nBinsPrev + 1; i++)
    {
        h->SetBinContent(i, 0);
    }

    //h->Clear();

    axis->Set(nBins, xBins);

    for (Int_t i = 1; i <= nBins; i++)
    {
        h->Fill(h->GetBinCenter(i), binContent[i]);
    }
}


/**
 * replaces the y-axis with a logarithmic y-axis of 100 bins
 * root has no logarithmic binning, so it requires generating a bin array
 * @param h
 *
 */
void logYaxis(TH2* h)
{
    TAxis* axis = h->GetYaxis();
    const Int_t nbins = 100;
    Double_t xmin = axis->GetXmin();
    Double_t xmax = axis->GetXmax();
    Double_t logXmin = TMath::Log10(xmin);
    Double_t logXmax= TMath::Log10(xmax);
    Double_t binwidth = (logXmax- logXmin) / nbins;
    Double_t xBins[nbins + 1];
    xBins[0] = xmin;
    for (Int_t i = 1; i <= nbins; i++)
    {
        xBins[i] = xmin + TMath::Power(10, logXmin + i * binwidth);
    }

    axis->Set(nbins, xBins);
}


/**
 * calculates neutron transport time in mus
 * z0 and z0Alt in mm, energy in MeV
 *
 * @param z0.
 * @param z0Alt.
 * @param energy
 * @param cosTheta
 * @return the Neutron transport Time .
 */
double calcNeutronDiffTime(double z0, double z0Alt, double energy, double cosTheta)
{
    //if (energy<20)	nSpeed = 3.9560339/sqrt(81.81/energy/1e9) * 1000.*1000.;
    //else 	nSpeed = 3.9560339/sqrt(81.81/20./1e9) * 1000.*1000.;
    //if (energy < nSpeedEnergyCutoff) nSpeed = 3.9560339/sqrt(81.81/nSpeedEnergyCutoff/1e9) * 1000.*1000.;

    if (z0Alt == z0) return 0;

    double nSpeed = 3.9560339 / sqrt(81.81 / energy / 1e9);

    //timeTemp = wwRange/nSpeed;
    return fabs((z0Alt - z0) / cosTheta / nSpeed);
}


/**
 * generates a random number in MeV from an evaporation spectrum with central energy theta
 * pointer to an already seeded Random generator is needed
 * @param theta
 * @param r
 * @return abszRnd  .
 *
 */
double getEvaporationEnergy(double theta, TRandom* r)
{
    bool gotIt;
    double abszRnd, ordRnd;

    TF1* spectrumFuncEvaporation = new TF1("spectrumFuncEvaporation", "x*TMath::Exp(-x/[0])", 0.0000, 20);
    spectrumFuncEvaporation->SetParameter(0, theta / 1e6);
    double maximum = spectrumFuncEvaporation->GetMaximum(0.2, 10);

    gotIt = false;
    while (!gotIt)
    {
        abszRnd = r->Rndm() * 19.;
        ordRnd = r->Rndm() * maximum;
        //gotIt = true;
        if (spectrumFuncEvaporation->Eval(abszRnd) > ordRnd) gotIt = true;
    }

    delete spectrumFuncEvaporation;

    return abszRnd;
}




/**
 * generates a random number in MeV from a moderated Californium source spectrum
 * pointer to an already seeded Random generator is needed
 * @param r
 * @return abszRnd
 */
double getModeratedCfEnergy(TRandom* r)
{
    bool gotIt;
    double abszRnd, ordRnd;

    //this function is in eV, not MeV like the others
    TF1* spectrumMaxwellPhiTempModModFission = new TF1("spectrumMaxwellPhiTempModModFission", "[0]/(TMath::Power([1]/11604.5,2))*TMath::Power(x,2)*TMath::Exp(-x/[1]*11604.5)+[2]*TMath::Power([1]/11604.5,1)/(1+TMath::Power(([3]/x),7)) * 1./(1-[4]/(1+TMath::Power(x/[5],5))) + [6]*x*TMath::SinH(sqrt(2.*(x-[7])/1e6))*TMath::Exp(-(x-[7])/1e6)", 1e-3, 1e7);
    spectrumMaxwellPhiTempModModFission->SetParameters(6868, 243, 9600, 0.7171, 0.6023, 0.2103, 0.003933, -3.394e6);

    //spectrumMaxwellPhiTempModModFission->SetParNames("Intensity th","Temperature", "Intensity ratio","scaling 1","scaling 2 [eV]","scaling 3 [eV]", "FissionPar1", "FissionPar2");
    double maximum = 3800;

    gotIt = false;
    while (!gotIt)
    {
        abszRnd = TMath::Power(10, r->Rndm() * 9. - 2.4);
        ordRnd = r->Rndm() * maximum;
        //gotIt = true;
        if (spectrumMaxwellPhiTempModModFission->Eval(abszRnd) > ordRnd) gotIt = true;
    }

    delete spectrumMaxwellPhiTempModModFission;

    return (1e-6 * abszRnd);
}


/**
 *  generates a random number in MeV from a fission spectrum
 * pointer to an already seeded Random generator is needed
 * @param r
 * @return abszRnd
 */
double getFissionEnergy(TRandom* r)
{
    bool gotIt;
    double abszRnd, ordRnd;

    TF1* spectrumFuncFission = new TF1("spectrumFuncFission", "0.4865*TMath::SinH(sqrt(2*x))*TMath::Exp(-x)", 0.0000, 10);

    double maximum = spectrumFuncFission->GetMaximum(0.2, 3);

    gotIt = false;
    while (!gotIt)
    {
        abszRnd = r->Rndm() * 10.;
        ordRnd = r->Rndm() * maximum;
        //gotIt = true;
        if (spectrumFuncFission->Eval(abszRnd) > ordRnd) gotIt = true;
    }
    delete spectrumFuncFission;

    return abszRnd;
}

/**
 *  generates a random number in MeV from a fission spectrum
 *  pointer to an already seeded Random generator is needed
 * @param r
 * @return abszRnd
 */
double getFissionEnergy2(TRandom* r)
{
    bool gotIt;
    double abszRnd, ordRnd;

    TF1* spectrumFuncFission = new TF1("spectrumFuncFission", "0.4865*TMath::SinH(sqrt(1.03419*x))*TMath::Exp(-x/1.18)", 0.0000, 10); // X-5 MONTE CARLO TEAM,LA-UR-03-1987,Los Alamos National  Laboratory

    double maximum = spectrumFuncFission->GetMaximum(0.2, 3);

    gotIt = false;
    while (!gotIt)
    {
        abszRnd = r->Rndm() * 10.;
        ordRnd = r->Rndm() * maximum;
        //gotIt = true;
        if (spectrumFuncFission->Eval(abszRnd) > ordRnd) gotIt = true;
    }

    delete spectrumFuncFission;

    return abszRnd;
}


/**
 *  generates a random Number representing energy in eV according to a thermal Spectrum
 *  collision partner of given mass massElm and temperature
 *  pointer to an already seeded Random generator is needed
 * @param r
 * @return vNeutron
 */
vector<float> getThermalPDF(const double nEnergy, const float massElm, const float temperature, TRandom* r)
{
    float xRnd, rnd1 = 0, rnd2 = 0, rnd3 = 0, abszRnd, ordRnd;
    vector<float> vNeutron;
    bool gotIt = false;
    bool takeTwo = false;

    double kB = 1.38065e-23;
    //double mN = 1.6749e-27;
    double mN = 939.565 * 1.1111111e-11;
    double sqrtPi = sqrt(TMath::Pi());
    double beta = sqrt((1.66e-27) * massElm * 0.5 / kB / temperature);

    double velR = 0.;
    double velT = 0.;
    double velN = sqrt(2. * nEnergy * 1e6 / mN);

    double y = velN * beta;

    TF1* spectrumFuncTwo = new TF1("spectrumFuncTwo", "2.256758353*TMath::Power(x,2)*TMath::Exp(-x*x)", 0.0000000, 100000);
    TF1* spectrumFuncThree = new TF1("spectrumFuncThree", "2.*TMath::Power(x,3)*TMath::Exp(-x*x)", 0.0000000, 100000);
    //Plot Range x: [0,3.2] y:[0,0.85]

    TF1* spectrumFunc;

    //calc which distribution to take: low energy or high energy
    rnd1 = r->Rndm();
    if (rnd1 < 2. / (sqrtPi * y + 2.))
    {
        takeTwo = false;
        spectrumFunc = spectrumFuncThree;
    }
    else
    {
        takeTwo = true;
        spectrumFunc = spectrumFuncTwo;
    }

    //get a thermal velocity from the Spectrum
    gotIt = false;
    while (!gotIt)
    {
        abszRnd = r->Rndm() * 3.2;
        ordRnd = r->Rndm() * 0.85;        
        if (spectrumFunc->Eval(abszRnd) > ordRnd) gotIt = true;
    }

    velT = abszRnd / beta;

    gotIt = false;
    while (!gotIt)
    {
        rnd2 = r->Rndm() * 2. - 1.;
        rnd3 = r->Rndm();

        velR = sqrt(fabs(TMath::Power(velN, 2) + TMath::Power(velT, 2) - 2. * velN * velT * rnd2));

        if (rnd3 < velR / (velN + velT)) gotIt = true;
    }

    // Energy(speed) and cosine
    vNeutron.push_back(velR * velR * mN / 2.);
    vNeutron.push_back(rnd2);

    delete spectrumFuncTwo; delete spectrumFuncThree;

    return vNeutron;
}

/**
 *  generates a random energy from a logarithmic thermal neutron density function
 *  in eV
 * @param spectrumFunc
 * @param r
 * @return xRnd
 */
double getThermalEnergyLog(const TF1* spectrumFunc, TRandom* r)
{
    double xRnd, yRnd;
    bool gotIt = false;

    while (!gotIt)
    {

        xRnd = TMath::Power(10, r->Rndm() * 2.7 - 8.9);
        yRnd = r->Rndm() * 12.;

        if (spectrumFunc->Eval(xRnd) > yRnd)
        {
            gotIt = true;
        }
    }
    return xRnd;
}


/**
 * generates a random energy from a linear thermal neutron density function
 * @param spectrumFunc
 * @param r
 *
 */
double getThermalEnergy(TF1* spectrumFunc, TRandom* r)
{
    double xRnd, yRnd;
    bool gotIt = false;

    while (!gotIt)
    {
        xRnd = r->Rndm() * 0.19e-6 + 0.1e-9;
        yRnd = r->Rndm();

        if (spectrumFunc->Eval(xRnd) > yRnd)
        {
            gotIt = true;
        }
    }

    return xRnd;
}


/**
 * generates a random energy from a logarithmic thermal neutron density function like the CF moderated spectrum
 * @param spectrumFunc
 * @param r
 *
 */
double getThermalEnergyFromSource(const TF1* spectrumFunc, TRandom* r)
{
    double xRnd, yRnd;
    bool gotIt = false;

    while (!gotIt)
    {
        xRnd = TMath::Power(10, r->Rndm() * 2.2 - 2.5);
        yRnd = r->Rndm() * 3800.;

        if (spectrumFunc->Eval(xRnd) > yRnd)
        {
            gotIt = true;
        }
    }
    return xRnd;
}

/**
 * calculates a scalar product
 * @param rx1
 * @param ry1
 * @param rz1
 * @param rx2
 * @param ry2
 * @param rz2
 * @returns - evaluated scalar product
 */
double scalarProduct(double rx1, double ry1, double rz1, double rx2, double ry2, double rz2)
{
    return(rx1 * rx2 + ry1 * ry2 + rz1 * rz2);
}

/**
 * finding the intersection of a line and a cylinder mantle
 * @param stVx
 * @param stVy
 * @param stVz
 * @param phi
 * @param theta
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @returns - intersection
 */
double intersectCylinderMantle(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad)
{
    double tValue1, tValue2, retValue;

    double sinThetacosPhi = (sin(theta)) * cos(phi);
    double sinThetasinPhi = (sin(theta)) * sin(phi);
    //double cosTheta = (cos(theta));

    double c = sinThetacosPhi * (stVx - px) + sinThetasinPhi * (stVy - py);
    double c2 = pow(px - stVx, 2) + pow(py - stVy, 2);

    tValue1 = -c + sqrt(c * c - (c2 - dRad * dRad));
    tValue2 = -c - sqrt(c * c - (c2 - dRad * dRad));

    double xIs1 = stVx + tValue1 * sinThetacosPhi;
    double yIs1 = stVy + tValue1 * sinThetasinPhi;

    if (scalarProduct(sinThetacosPhi, sinThetasinPhi, 0, (xIs1 - px), (yIs1 - py), 0) >= 0)
    {
        retValue = tValue2;
    }
    else
    {
        retValue = tValue1;
    }

    return retValue;
}


/**
 * finding the intersection of sphere with a line
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @returns - intersection
 */
double intersectSphere(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad)
{
    double tValue1, tValue2, retValue;

    double sinThetacosPhi = (sin(theta)) * cos(phi);
    double sinThetasinPhi = (sin(theta)) * sin(phi);
    double cosTheta = cos(theta);

    double c = sinThetacosPhi * (stVx - px) + sinThetasinPhi * (stVy - py) + cosTheta * (stVz - pz);
    double c2 = pow(px - stVx, 2) + pow(py - stVy, 2) + pow(pz - stVz, 2);

    tValue1 = -c + sqrt(c * c - (c2 - dRad * dRad));
    tValue2 = -c - sqrt(c * c - (c2 - dRad * dRad));

    double xIs1 = stVx + tValue1 * sinThetacosPhi;
    double yIs1 = stVy + tValue1 * sinThetasinPhi;
    /*
    double zIs1 = stVz+tValue1*cos(theta);
    double xIs2 = stVx+tValue2*sin(theta)*cos(phi);
    double yIs2 = stVy+tValue2*sin(theta)*sin(phi);
    double zIs2 = stVz+tValue2*cos(theta);
    */

    bool secondPoint = false;
    if (scalarProduct(sinThetacosPhi, sinThetasinPhi, 0, (xIs1 - px), (yIs1 - py), 0) >= 0)
    {
        secondPoint = true; //take no2
        retValue = tValue2;
    }
    else
    {
        secondPoint = false;
        retValue = tValue1;
    }

    return retValue;
}

/**
 * Deprecated
 * Calculate the Cylindrical hit distance
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @param xref
 * @param yref
 * @returns - hit distance
 */
double calcCylindricalHitDist2(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad, double xRef, double yRef)
{
    double tValue, retValue;

    tValue = intersectCylinderMantle(stVx, stVy, stVz, theta, phi, px, py, pz, dRad);

    double xIs2 = stVx + tValue * sin(theta) * cos(phi);
    double yIs2 = stVy + tValue * sin(theta) * sin(phi);
    //double zIs2 = stVz+tValue*cos(theta);

    retValue = sqrt(pow(xRef - xIs2, 2) + pow(yRef - yIs2, 2));

    return retValue;
}


/**
 * Calculate the Cylindrical hit distance
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @param xref
 * @param yref
 * @returns - hit distance
 */
double calcCylindricalHitDist(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad, double xRef, double yRef)
{
    double tValue, retValue;

    tValue = intersectSphere(stVx, stVy, stVz, theta, phi, px, py, pz, dRad);

    double xIs2 = stVx + tValue * sin(theta) * cos(phi);
    double yIs2 = stVy + tValue * sin(theta) * sin(phi);
    //double zIs2 = stVz+tValue*cos(theta);

    retValue = sqrt(pow(xRef - xIs2, 2) + pow(yRef - yIs2, 2));

    return retValue;
}

/**
 * Calculate the cylinder intersection
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @returns - cylinder intersection
 */
double intersectCylinder(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad)
{
    //takes the projection of the cylinder in the x-y plane, selects between the two intersection points
    double tValue1, tValue2;

    double a = pow((sin(theta) * cos(phi)), 2) + pow((sin(theta) * sin(phi)), 2);
    double p = sin(theta) * cos(phi) * (stVx - px) + sin(theta) * sin(phi) * (stVy - py) / a;
    double q = pow((stVx - px), 2) + pow((stVy - py), 2) - dRad * dRad;

    tValue1 = -0.5 * p + sqrt(pow(0.5 * p, 2) - q);
    tValue2 = -0.5 * p - sqrt(pow(0.5 * p, 2) - q);

    double xIs1 = stVx + tValue1 * sin(theta) * cos(phi);
    double yIs1 = stVy + tValue1 * sin(theta) * sin(phi);
    double zIs1 = stVz + tValue1 * cos(theta);
    double xIs2 = stVx + tValue2 * sin(theta) * cos(phi);
    double yIs2 = stVy + tValue2 * sin(theta) * sin(phi);
    double zIs2 = stVz + tValue2 * cos(theta);

    bool secondPoint = false;
    if (scalarProduct(sin(theta) * cos(phi), sin(theta) * sin(phi), 0, (xIs1 - px), (yIs1 - py), 0) > 0) secondPoint = true; //take no2
    else secondPoint = false;

    return 0; //unfinished
}


/**
 * calculates 3D track ('stVx', 'stVy', 'stVz', 'phi', 'theta') to a line ('px','py','pz') distance
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @returns - distance to line
 */
double getDistanceToLine(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz)
{
    double distance;

    theta = TMath::Pi() - theta;

    double a1 = sin(theta) * sin(phi);
    double a2 = -sin(theta) * cos(phi);

    distance = abs(((px - stVx) * a1 + (py - stVy) * a2) / sqrt(a1 * a1 + a2 * a2));

    return distance;
}


/**
 * calculates 3D track ('stVx', 'stVy', 'stVz', 'phi', 'theta') to point  ('px','py','pz') distance
 * @param stVx
 * @param stVy
 * @param stVz
 * @param theta
 * @param phi
 * @param px
 * @param py
 * @param pz
 * @param dRad
 * @returns - distance to point
 */
double getDistanceToPoint(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz)
{
    double distance;

    theta = TMath::Pi() - theta; //changed 29.11.2018

    double gx = sin(theta) * cos(phi);
    double gy = sin(theta) * sin(phi);
    double gz = -cos(theta); //changed 30.11.2018

    double krx = gy * (stVz - pz) - gz * (stVy - py);
    double kry = gz * (stVx - px) - gx * (stVz - pz);
    double krz = gx * (stVy - py) - gy * (stVx - px);

    distance = sqrt(krx * krx + kry * kry + krz * krz) / sqrt(gx * gx + gy * gy + gz * gz);

    return distance;
}


//DEPRECATED
double getDistanceToPointOLD(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz)
{
    double distance;

    double tx = cos(phi) * (stVz - pz) - tan(theta) * (stVy - py);
    double ty = tan(theta) * (stVx - px) - sin(phi) * (stVz - pz);
    double tz = sin(phi) * (stVy - py) - cos(phi) * (stVx - px);

    distance = sqrt(1. / (1. + tan(theta) * tan(theta))) * sqrt(tx * tx + ty * ty + tz * tz);

    return distance;
}

/**
 * returns a e polynom of order par[1] at x, given scale factor par[0]
 * @param x
 * @param par
 * @returns - legrendian
 */
Double_t legendrian(double* x, Double_t* par)
{
    return par[0] * legendre_Pl(par[1], x[0]);
}

/**
 * returns a normalized sum of Legendre polynoms of order 0 to 10 at x
 * @param x
 * @param par
 * @returns - legrendian 10 fold
 */
Double_t legendrian10fold(double* x, Double_t* par)
{
    x[0] = TMath::Cos(x[0]);
    return 0.159155 * (0.5 * 1. * legendre_Pl(0, x[0]) + 1.5 * par[0] * legendre_Pl(1, x[0]) + 2.5 * par[1] * legendre_Pl(2, x[0]) + 3.5 * par[2] * legendre_Pl(3, x[0]) + 4.5 * par[3] * legendre_Pl(4, x[0]) + 5.5 * par[4] * legendre_Pl(5, x[0]) + 6.5 * par[5] * legendre_Pl(6, x[0]) + 7.5 * par[6] * legendre_Pl(7, x[0]) + 8.5 * par[7] * legendre_Pl(8, x[0]) + 9.5 * par[8] * legendre_Pl(9, x[0]) + 10.5 * par[9] * legendre_Pl(10, x[0]));
}


/**
 * returns a normalized sum of Legendre polynoms of order 0 to 10 at x as a STRING
 * @param par
 * @returns - legrendian 10 fold as string
 */
string legendrian10folded(Float_t* par)
{
    return "0.159155*(0.5*1*legendre_Pl(0, TMath::Cos(x)) + 1.5*" + castDoubleToString(par[0]) + "*legendre_Pl(1, TMath::Cos(x)) + 2.5*" + castDoubleToString(par[1]) + "*legendre_Pl(2, TMath::Cos(x)) + 3.5*" + castDoubleToString(par[2]) + "*legendre_Pl(3,TMath::Cos(x)) + 4.5*" + castDoubleToString(par[3]) + "*legendre_Pl(4,TMath::Cos(x)) + 5.5*" + castDoubleToString(par[4]) + "*legendre_Pl(5, TMath::Cos(x)) + 6.5*" + castDoubleToString(par[5]) + "*legendre_Pl(6, TMath::Cos(x))";
}

/**
 * returns the size of an array
 * @param array
 * @returns - array size
 */
static int ArraySize(double array[])
{
    int i = 0;
    while (array[i] != NULL) i++;
    return i;
}

int legendreElements = 0;

/**
 * returns a normalized sum of Legendre polynoms of order 0 to N (given by the size of the parameter array, which sets the relative scaling factors) at x,
 * @param x
 * @param par
 * @returns - legrendian N fold
 */
static double legendrianNfold(double* x, Double_t* par)
{
    double result = 0.5;

    int elements = legendreElements;

    x[0] = TMath::Cos(x[0]);

    for (int l = 1; l < elements; l++)
    {
        result = result + (2. * (l)+1.) * 0.5 * par[l - 1] * legendre_Pl(l, x[0]);
    }

    return 0.159155 * result;
}


/**
 * converts an endf formatted string to a float number
 * @param str
 * @returns - float of a string
 */
float endfNumberConv(string str)
{
    float pos;
    int length = str.length();

    if (length < 5) return 0;
    if (str.substr(5, 1) == " ") return 0;
    string sign = str.substr(length - 2, 1);
    if ((sign == "+") || (sign == "-")) pos = 2;
    sign = str.substr(length - 3, 1);
    if ((sign == "+") || (sign == "-")) pos = 3;

    string part1 = str.substr(0, length - pos);
    string part2 = str.substr(length - pos, pos);
    float number = atof(part1.c_str()) * TMath::Power(10, atof(part2.c_str()));

    return number;
}


/**
 * searches a value in a matrix by an interval algorithm
 * searches the cumulated probability distribution (at row line) for a given probability (value) and extrapolates linearly between upper and lower bound
 * additionally log interval search can be activated
 * @param matrix
 * @param line
 * @param value
 * @param doLogSearch
 * @returns - index of the horizontal position
 */
float getIndexHorizontalPosition(const TMatrixF& matrix, int line, double value, bool doLogSearch)
{
    float result;

    //int index;

    int counter = 0;

    int cols = matrix.GetNcols();

    int pos;

    //energy is stored in column lower bound
    int left = matrix.GetColLwb() + 1;
    int right = matrix.GetColUpb();

    if (value > matrix(line, right)) return right - 1;
    if (value < matrix(line, left)) return left;

    while (1)
    {
        if (doLogSearch)
        {
            //pos = left+ (log(value)-log(matrix(left,0)))/(log(matrix(right,0))-log(matrix(left,0))) * (right-left);
        }
        else
        {
            //pos = left+ (value-matrix(left,0))/(matrix(right,0)-matrix(left,0)) * (right-left);
        }

        pos = left + (right - left) * 0.5;

        if (pos >= cols - 2) pos = cols - 2;
        if (pos <= 0) pos = 0;

        if (matrix(line, pos) < value)
        {
            if ((pos < cols - 1) && (matrix(line, pos + 1) > value)) break;
            left = pos + 1;
        }
        else
        {
            if ((pos > 0) && (matrix(line, pos - 1) < value)) { pos = pos - 1; break; }
            right = pos - 1;
        }

        counter++;
        if ((counter > 300) && (!doLogSearch)) doLogSearch = true;
        if ((counter > 100) && (doLogSearch)) doLogSearch = false;

        if (counter > 500)
        {
            cout << "ERROR in Matrix Search" << endl;
            if (doLogSearch) cout << "Log  "; cout << "Value: " << value << " failed at position: " << pos << " in " << cols << " columns at line " << line << endl;
            return 0;
        }
    }

    result = pos + (value - matrix(line, pos)) / (matrix(line, pos + 1) - matrix(line, pos));

    if (result >= pos + 1) result = pos + 0.999;

    return result;
}


/**
 * searches a value in a matrix by an interval algorithm
 * searches the index/energy (first column) and extrapolates linearly between upper and lower bound
 * additionally log interval search can be activated
 * @param matrix
 * @param line
 * @param value
 * @param doLogSearch
 * @returns - index of the position
 */
float getIndexPosition(const TMatrixF& matrix, double value, bool doLogSearch)
{
    float result;

    int counter = 0;

    int rows = matrix.GetNrows();

    int pos;

    int left = matrix.GetRowLwb();
    int right = matrix.GetRowUpb();

    if (value > matrix(right, 0)) return right - 1;
    if (value < matrix(left, 0)) return left;

    while (1)
    {
        if (doLogSearch)
        {
            //pos = left+ (log(value)-log(matrix(left,0)))/(log(matrix(right,0))-log(matrix(left,0))) * (right-left);
        }
        else
        {
            //pos = left+ (value-matrix(left,0))/(matrix(right,0)-matrix(left,0)) * (right-left);
        }

        pos = left + (right - left) * 0.5;

        if (pos >= rows - 2) pos = rows - 2;
        if (pos <= 0) pos = 0;

        if (matrix(pos, 0) < value)
        {
            if ((pos < rows - 1) && (matrix(pos + 1, 0) > value)) break;
            left = pos + 1;
        }
        else
        {
            if ((pos > 0) && (matrix(pos - 1, 0) < value)) { pos = pos - 1; break; }
            right = pos - 1;
        }

        counter++;
        if ((counter > 300) && (!doLogSearch)) doLogSearch = true;
        if ((counter > 100) && (doLogSearch)) doLogSearch = false;

        if (counter > 500)
        {
            cout << "ERROR in Matrix Search" << endl;
            if (doLogSearch) cout << "Log  "; cout << "Value: " << value << " failed at position: " << pos << " in " << rows << " rows" << endl;
            return 0;
        }
    }

    result = pos + (value - matrix(pos, 0)) / (matrix(pos + 1, 0) - matrix(pos, 0));

    if (result >= pos + 1) result = pos + 0.999;

    return result;
}




/**
 * searches a matrix for an energy value (first column) and returns a linearly extrapolated cross section value (second column)
 * @param matrix
 * @param energy
 * @returns - cross section value
 */
double calcMeanCS(const TMatrixF& matrix, double energy)
{
    float position = getIndexPosition(matrix, energy, true);
    int index = position;
    float frac = position - index;

    return  (1. - frac) * matrix(index, 1) + frac * matrix(index + 1, 1);
}

/**
 * returns cos(thetha) for high energy angle distributions
 * for a given energy the index in the the angular distribution probability matrix is searched
 * then a linearly interpolated cos(theta) from a random number prob is returned
 * @param angleMatrix
 * @param cumulatedProbMatrix
 * @param energy
 * @param prob
 * @returns - High Energy Cos Theta
 */
double getHighEnergyCosTheta(const TMatrixF& angleMatrix, const TMatrixF& cumulatedProbMatrix, double energy, double prob)
{
    float energyPosition = getIndexPosition(angleMatrix, energy, false);
    int energyIndex = energyPosition;
    float energyFrac = energyPosition - energyIndex;

    float anglePosition1 = getIndexHorizontalPosition(cumulatedProbMatrix, energyIndex, prob, false);
    int angleIndex1 = anglePosition1;

    float anglePosition2 = getIndexHorizontalPosition(cumulatedProbMatrix, energyIndex + 1, prob, false);
    int angleIndex2 = anglePosition2;

    float angleFrac1 = anglePosition1 - angleIndex1; //creates the 0.x fracions
    float angleFrac2 = anglePosition2 - angleIndex2;

    double cosThetaLw = (1. - angleFrac1) * angleMatrix(energyIndex, angleIndex1) + angleFrac1 * angleMatrix(energyIndex, angleIndex1 + 1);
    double cosThetaHigh = (1. - angleFrac2) * angleMatrix(energyIndex + 1, angleIndex2) + angleFrac2 * angleMatrix(energyIndex + 1, angleIndex2 + 1);

    double cosTheta = (1. - energyFrac) * cosThetaLw + energyFrac * cosThetaHigh;

    return cosTheta;
}

/**
 * generates a random number (from a cumulative dsitribution function) in the range min to max
 * pointer to an already seeded Random generator is needed
 * @param spectrumFunc
 * @param min
 * @param max
 * @param r
 * @returns - angle from cumulative function
 */
double getAngleFromCumulativeFunction(const TF1* spectrumFunc, float min, float max, TRandom* r)
{
    const int steps = 51;

    float result = 0;

    // (cumulative Function | xValues)
    TMatrixF cumulativeValueMatrix(steps, 2);

    cumulativeValueMatrix(0, 0) = 0; cumulativeValueMatrix(0, 1) = 0;

    //float stepwidth = fabs(max-min/(steps-1));
    double stepwidth = fabs((max - min) / ((steps - 1) * 1.));

    for (int i = 1; i < steps; i++)
    {
        cumulativeValueMatrix(i, 1) = min + (i - 1) * stepwidth;
        cumulativeValueMatrix(i, 0) = cumulativeValueMatrix(i - 1, 0) + spectrumFunc->Eval(min + (i - 1) * stepwidth);
    }

    result = getIndexPosition(cumulativeValueMatrix, cumulativeValueMatrix(steps - 1, 0) * r->Rndm(), false);
    int bin = result;
    float frac = result - bin;

    if (bin >= (steps - 1)) bin--;
    return (1. - frac) * cumulativeValueMatrix(bin, 1) + frac * cumulativeValueMatrix(bin + 1, 1);
}


/**
 * generates a legendre function according to the polynomial coefficients which are linearly extrapolated from the matrix according to the given energy
 * @param matrix
 * @param energy
 * @returns - mean Angular Distribution
 */
TF1* calcMeanAngularDistribution(const TMatrixF& matrix, double energy)
{
    Float_t cn1[21];
    Float_t cn2[21];

    cn1[0] = 0; cn1[1] = 0;
    cn2[0] = 0; cn2[1] = 0;

    float position = getIndexPosition(matrix, energy, false);

    int index = position;
    float frac = position - index;
    int  l;

    for (l = 1; l < 21; l++)
    {
        cn1[l - 1] = matrix(index, l);
        cn2[l - 1] = matrix(index + 1, l);
        if ((cn1[l - 1] == 0) && (cn2[l - 1] == 0)) break;
        if ((frac == 0) && (cn1[l - 1] == 0)) break;
        if (l == 21) break;

    }

    if (l < 3) { l = 3; frac = 0; }
    const int lastcn = l - 1;

    legendreElements = lastcn;
    TF1* legendrian = new TF1("legendrian", legendrianNfold, 0, TMath::Pi(), lastcn);

    for (int l = 0; l < lastcn; l++)
    {
        legendrian->SetParameter(l, (1 - frac) * cn1[l] + frac * cn2[l]);
    }

    return legendrian;
}

/**
 * transposes a matrix (esp. for input images)
 * @param inputMatrix
 *
 */
void turnInputMatrix(TMatrixF& inputMatrix)
{
    Float_t tmpchg;
    Int_t imageSize = inputMatrix.GetNrows();
    for (Int_t im = 0; im < imageSize / 2.; im++)
    {
        for (Int_t jm = 0; jm < imageSize; jm++)
        {
            tmpchg = inputMatrix(jm, im);
            inputMatrix(jm, im) = inputMatrix(jm, imageSize - 1 - im);
            inputMatrix(jm, imageSize - 1 - im) = tmpchg;
        }
    }
    for (Int_t im = 0; im < imageSize; im++)
    {
        for (Int_t jm = 0; jm < imageSize; jm++)
        {
            if (imageSize - 1 - im < jm)
            {
                tmpchg = inputMatrix(jm, im);
                inputMatrix(jm, im) = inputMatrix(imageSize - 1 - im, imageSize - 1 - jm);
                inputMatrix(imageSize - 1 - im, imageSize - 1 - jm) = tmpchg;
            }
        }
    }
}


/**
 * converts old gray scale identifiers to new ones
 * @param inputMatrix
 * @returns - converted old gray scale
 */
int changeInputMatrixValue(int value)
{
    return 0;
}

