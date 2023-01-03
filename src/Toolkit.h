#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <Riostream.h> 
#include "TFile.h"
#include "TText.h"
//#include "TTree.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h" 
#include "TH3.h" 
#include "TF1.h" 
#include "TF2.h" 
#include "TProfile.h"
#include "Fit/FitResult.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TPad.h"
//#include "TVirtualPad.h"
//#include "TVirtualFitter.h"
#include "TNamed.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TVectorT.h"
#include "TMatrixF.h"
#include "TMatrixD.h" 
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLatex.h"
#include "TColor.h"
#include <string>
#include "RConfig.h"
#include "Rtypes.h" 
#include "TLegend.h"
#include "TSpline.h"
#include "TObject.h" 
#include <strstream>
#include <sstream>
//#include <vector>
#include <float.h>

#include <malloc.h>

#include <QImage>

#include "TRandom3.h"

using namespace std;


	string  castDoubleToString(double number);

	string  castDoubleToString(double number, Int_t intacters);

	string  castFloatToString(Float_t number);

	string  castFloatToString(Float_t number, Int_t intacters);

	string  castIntToString(Int_t &number);

    string  castLongToString(Long_t &number);

	Double_t  lorentzianPeak(Double_t* x, Double_t* par);

	Double_t  gaussoffset(Double_t* x, Double_t* par);

	Double_t errf( Double_t *x, Double_t *par);

	void  rootlogon();

    vector<float> getRGBfromHCL(double h, double c0, double l0);

    float getLinearC(double factor, double cmin, double cmax);

    float getScaledColorValue(double factor, double pow, double cmin, double cmax);

	void  set_plot_styleSingleGradient(float h, Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_styleHeatGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_styleRainbowGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_styleAllGradient(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_styleHeatGradientModified(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_styleHeatGradient2(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	void  set_plot_style(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

	int getScaledColor(float min, float max, float scale, int mode);

	vector<float> getScaledColorRGB(float min, float max, float scale, int mode);

	void set_plot_styleCool();

	void  CanvasFashion(TCanvas* c);

	void  TGraphFashion(TGraph* graph, TString xLabel, TString yLabel, Bool_t setSize);
	
	void  TGraphErrorFashion(TGraphErrors* graph, TString xLabel, TString yLabel, Bool_t setSize);

	void  TGraph2DFashion(TGraph2D* graph, TString xLabel, TString yLabel, TString zLabel, Bool_t setSize);

	void  TMultiGraphFashion(TMultiGraph* graph, TString xLabel, TString yLabel, Bool_t setSize);

	void  TF1Fashion(TF1* hist, TString xLabel, TString yLabel, Bool_t setSize);

	void  TH1Fashion(TH1* hist, TString xLabel, TString yLabel, Bool_t setSize);

	void  TProfileFashion(TProfile* hist, TString xLabel, TString yLabel, Bool_t setSize);

	void  TH2Fashion(TH2* hist, Bool_t setSize);

	void  TH2Fashion(TH2* hist, TString xLabel, TString yLabel, Bool_t setSize);

	void  THStackFashion(THStack* stack, TString xLabel, TString yLabel, Bool_t setSize);	

	TH2F* printHisto(Double_t zmin, Double_t zmax, Int_t size, TMatrixF data, TString OutputFolder, TString dataname, Bool_t logPlot);

	void  printTH1F(Double_t zmin, Double_t zmax, TH1F* histo, TString OutputFolder, TString dataname);

	void  printTH2F(Double_t zmin, Double_t zmax, TH2F* histo, TString OutputFolder, TString dataname);

	void  storeTH2ToFile(TString OutputFolder, TString filename, TH2F* th2);

	void  storeToFile(TString OutputFolder, TString filename, TString what_has_to_be_written, bool new_line);

	void  embedFonts(TString path);

	void  storeTH1ToFile(TString OutputFolder, TString filename,TH1* histo);    

    TMatrixF readMatrixPNG(TString folder, TString filename);

	TMatrixF  readmatrix(TString folder, TString filename, TString filetype, Int_t counter, Int_t size);

	Int_t  countFiles(TString folder, TString filename, TString filetype);

	//void deleteStatsTH(vector<TH1*> *allTHs);

    unsigned int heXheX(const char *value);

	const std::string intToHex( int i);

	void extrapolateZeroValues(TMatrixF &matrix);	

	TMatrixF reduceMatrix(TMatrixF &matrix, Float_t factor);

	TMatrixF getGradientMatrix(TMatrixF &matrix);

    void getGradientMatrixFromTH2(TH2F* matrix, TH2F* newMatrix);

	TH2F* getTH2fromMatrix(TMatrixF &matrix, TString thName, TString thTitle);

	TH1F* convoluteGaussian(TH1F* hist, Double_t gaussRelW);

	int countEventsInEllipse(TH2F* histo, double radiusX, double x, double radiusY, double y);

	void histoClearUp(vector<TObject*>* vectorTH);

    double getWWProb(double density, double atWeight, double cs, double lambda);
    double getLfromE(double energy);
    double getEfromL(double lambda);
    TSpline3* getSplinedDetectorEnergyModel();
    int getNumberOfFileEntries(TString dname);
    inline double legendre_Pl(int n, double x);
    void logXaxis(TH2* h);
    void rebinX(TH1* h);
    void logaxis(TH1* h);
    void logXaxis(TH2* h);
    void logYaxis(TH2* h);
    double getRLuftWasser(float temperature);
    double calcNeutronDiffTime(double z0, double z0Alt, double energy, double cosTheta);
    double getEvaporationEnergy(double theta, TRandom* r);
    double getModeratedCfEnergy(TRandom* r);
    double getFissionEnergy(TRandom* r);
    double getFissionEnergy2(TRandom* r);
    vector<float> getThermalPDF(const double nEnergy, const float massElm, const float temperature, TRandom* r);
    double getThermalEnergyLog(const TF1* spectrumFunc, TRandom* r);
    double getThermalEnergy(TF1* spectrumFunc, TRandom* r);
    double getThermalEnergyFromSource(const TF1* spectrumFunc, TRandom* r);
    double scalarProduct(double rx1, double ry1, double rz1, double rx2, double ry2, double rz2);
    double intersectCylinderMantle(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad);
    double intersectSphere(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad);
    double calcCylindricalHitDist2(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad, double xRef, double yRef);
    double calcCylindricalHitDist(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad, double xRef, double yRef);
    double intersectCylinder(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz, double dRad);
    double getDistanceToLine(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz);
    double getDistanceToPoint(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz);
    double getDistanceToPointOLD(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz);
    Double_t legendrian(double* x, Double_t* par);
    Double_t legendrian10fold(double* x, Double_t* par);
    string legendrian10folded(Float_t* par);
    static int ArraySize(double array[]);
    static double legendrianNfold(double* x, Double_t* par);
    float endfNumberConv(string str);
    float getIndexHorizontalPosition(const TMatrixF& matrix, int line, double value, bool doLogSearch);
    float getIndexPosition(const TMatrixF& matrix, double value, bool doLogSearch);
    double calcMeanCS(const TMatrixF& matrix, double energy);
    double getHighEnergyCosTheta(const TMatrixF& angleMatrix, const TMatrixF& cumulatedProbMatrix, double energy, double prob);
    double getAngleFromCumulativeFunction(const TF1* spectrumFunc, float min, float max, TRandom* r);
    TF1* calcMeanAngularDistribution(const TMatrixF& matrix, double energy);
    void turnInputMatrix(TMatrixF& inputMatrix);
    int changeInputMatrixValue(int value);






