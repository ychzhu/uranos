/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "Toolkit.h"

#include "TRandom3.h"
#include "time.h"


#include <vector>
#include <QImage>
#include <QMessageBox>
#include <QDebug>

#include <algorithm>

//wchar_t BOM = 0xFEFF;


#ifdef _WIN32
    #include <io.h>
    #include <omp.h>
    #define access _access
#elif __linux_
#endif

#ifdef _WIN32
#pragma warning (disable: 4018 4100 4101 4189 4305)
#elif __linux_
#endif

//Initialization of Objects for data transfer
//not yet in order

TRandom3 r; // general random generator in ROOT

VisualizationEnlarge* visualization;  // the 1000 x 1000 pixel window for the bird's eyes view
VisualizationEnlarge2* visualization2;

// objects used for plotting data in QCustomplot (mainly vectors for x and y axis data)
vector< QVector<double>* > allDisplayData;

QVector<double> x(101), y(101);
QVector<double> plotGUIxBinsIncFullSpectrum(10000), plotGUIyBinsIncFullSpectrum(10000);
QVector<double> plotGUIxBinsFootprFunc(200), plotGUIyBinsFootprFunc(200);
QVector<double> plotGUIxBinsIncOnlySpectrum(8000), plotGUIyBinsIncOnlySpectrum(8000);
QVector<double> plotGUIxBinsScatteredSurfaceSpec(8000), plotGUIyBinsScatteredSurfaceSpec(8000);
QVector<double> plotGUIxBinsScatteredSurfaceSpecAlbedo(8000), plotGUIyBinsScatteredSurfaceSpecAlbedo(8000);
QVector<double> plotGUIxBinsCutView(500), plotGUIyBinsCutView(500);
QVector<double> plotGUIxBinsDistanceDet(4000), plotGUIyBinsDistanceDet(4000);
QVector<double> plotGUIxBinsDistanceAlbedoDet(4000), plotGUIyBinsDistanceAlbedoDet(4000);
QVector<double> plotGUIxBinsDistanceDetLayer(4000), plotGUIyBinsDistanceDetLayer(4000);
QVector<double> plotGUIxBinsDistanceAlbedoDetLayer(4000), plotGUIyBinsDistanceAlbedoDetLayer(4000);
QVector<double> xRS(4000), yRS(4000);
QVector<double> plotGUIxBinsScatDepth(2500), plotGUIyBinsScatDepth(2500);
QVector<double> plotGUIxBinsScatDepthDet(2500), plotGUIyBinsScatDepthDet(2500);
QVector<double> plotGUIxBinsScatDepthDetMax(2500), plotGUIyBinsScatDepthDetMax(2500);

QVector<double> plotGUIxBinsFootprFuncLine(4);
QVector<double> plotGUIyBinsFootprFuncLine(4);

QVector<double> x2(501), y2(501);

vector<TH1*> liveTHs; // list of only visible graphs

bool setSize = true;    // large or small labels in the ROOT output
float squareDim = 100000;      //length of the simulated area (square) [mm]
float beamRadius = squareDim * 0.49; //expansion of the beam/source

bool domainCutoff = false;               //cut  neutrons outside the domain: factor
bool domainCutoff2 = false;              //cut  neutrons outside the domain: meters
int domainCutoffMeters = 500;
double domainCutoffFactor = 1.5;

// all ROOT histograms used in the simulation
// energy categories:
// Thermal, epithermal, fast, high enery and detector selection (which is called albedo here)
TH2F* densityTrackMap;
TH2F* densityIntermediateTrackMap;
TH2F* densityFastTrackMap;
TH2F* densityAlbedoTrackMap;
TH2F* densityHighEnergyTrackMap;
TH2F* densityEnergyTrackMap;

TH2F* densityThermalTrackMap;

TH2F* densityTrackMapHighRes;
TH2F* densityIntermediateTrackMapHighRes;
TH2F* densityFastTrackMapHighRes;
TH2F* densityAlbedoTrackMapHighRes;
TH2F* densityHighEnergyTrackMapHighRes;
TH2F* densityEnergyTrackMapHighRes;

TH2F* densityThermalTrackMapHighRes;

TH2F* densityTrackMapHighRes15x;
TH2F* densityIntermediateTrackMapHighRes15x;
TH2F* densityFastTrackMapHighRes15x;
TH2F* densityAlbedoTrackMapHighRes15x;
TH2F* densityHighEnergyTrackMapHighRes15x;
TH2F* densityEnergyTrackMapHighRes15x;

TH2F* densityThermalTrackMapHighRes15x;

TH2F* densityFastTrackMapHighRes2x;

TH2F* densityMapThermal;

TH2F* densityMap;
TH2F* densityMapIntermediate;
TH2F* densityMapFast;
TH2F* densityMapAlbedo;
TH2F* densityMapHighEnergy;

TH2F* densityTrackMapSide;
TH2F* densityTrackMapSideAlbedo;
TH2F* densityTrackMapSideDetector;
TH2F* densityTrackMapSideThermal;

TH1F* cosmicSpectrum = new TH1F("cosmicSpectrum", "Cosmic Spectrum", 10000, 1e-9, 10000);
TH1F* scatteredSurfaceSpectrum = new TH1F("scatteredSurfaceSpectrum", "Scattered Neutron at Surface Spectrum", 10000, 1e-9, 10000);
TH1F* scatteredSurfaceSpectrumBack = new TH1F("scatteredSurfaceSpectrumBack", "Back Scattered Neutron at Surface Spectrum", 10000, 1e-9, 10000);
TH1F* scatteredSurfaceSpectrumHelp = new TH1F("scatteredSurfaceSpectrumHelp", "Non-Scattered Evaporated Neutron at Surface Spectrum", 10000, 1e-9, 10000);

TH1F* detectorDistanceBackScattered;

TH1F* detectorLayerDistance;
TH1F* detectorLayerDistanceBackscattered;

TH1F* detectorDistance;
TH1F* scatteredSurfaceDistance;
TH2F* detectorOriginMap;

TH1F* scatDepth = new TH1F("scatDepth", "Scattered Depth", 2501, -.5, 2500.5);
TH1F* scatteredSurfaceDepth = new TH1F("scatteredSurfaceDepth", "Scattered Depth Neutron on Surface", 2501, -.5, 2500.5);
TH1F* scatteredSurfaceMaxDepth = new TH1F("scatteredSurfaceMaxDepth", "Scattered Maximum Depth Neutron on Surface", 2501, -.5, 2500.5);

TString outputFolder = "";

TString workFolder = "";

string configFilePath = "";
bool configFilePathConfigured = false;

TString inputSpectrumFile = "";

TString detectorResponseFunctionFile = "", detectorResponseFunctionFile2 = "";

// folder for cross sections and angular distributions
TString endfFolder = "";

bool simulationRunning = false;     // GUI
bool silderColorMoved = false;      // GUI
bool silderDetectorMoved = false;   // GUI
bool silderDetectorColorMoved = false; // GUI

bool dontFillunecessaryPlots = false;// disables Scoring of internally used histograms
bool activateFP = false;            // GUI
bool trackAllLayers = false;        // activates tracking for all layers instead of just the detector layer
bool plotTopViewLog = false;        // GUI
bool noGUIMode = false;             // deactivates the GUI for command line run
bool silentMode = false;            // deactivates the message output
bool layerMapsImport = false;       // GUI
bool highResCalc = false;           // track resolution mode with 1000 px
bool highhighResCalc = false;       // extra resolution mode with 1500 px

float densitySideTrackingEnergyCutoff = 20000.; //upper cutoff for sideways tracking... not so important any more

bool useManualColors = false;       // GUI
int manualColor = 0;                // GUI
int manualColorZero = 0;            // GUI

bool useDetectorSensitiveMaterial = false;  // a detector layer override function which makes a normal voxel to a detector
int detectorSensitiveMaterial = 11;

bool doSkyEvaporation = false;      // configures for manual source placement instead of the total source layer
bool noThermalRegime = true;        // thermal cutoff
bool noMultipleScatteringRecording = false;
bool doBatchRun = false;            // batch run with parameters set later
bool doBatchRunDensity = false;      // batch run for soil moisture and density variation
bool doBatchRun2D = false;           // batch run for sm and hum
bool doDetectorBatchRun = false;    // batch run for detector analysis with parameters set later
bool doDetectorBatchRun2 = false;    // batch run for detector analysis with parameters set later
bool doDetectorAngleBatchRun = false;    // batch run for detector analysis with parameters set later
bool regionalthermalcutoff = false;  // kills thermal neutrons which scattered more far away from 0,0 than the radius given
double regionalthermalcutoffradius = 50e3;
bool useExtraCounter = false;       // activates the counter in the live view to display extra counts of for example absorbed thermal neutrons
bool recordSubsurfaceScatterings = true;  // records all scatterings inside the ground

bool doFusion = false;              // 14 MeV neutrons
bool doFission = false;             // Cf Fission spectrum
bool doAmBe = false;                // AmBe Iso spectrum
bool doMonoEnergetic = false;       // monoenergetic neutrons of given energy
bool doThermalSource = false;       // thermal spectrum
bool doNoSource = false;            // cosmic spectrum
bool doModeratedCf = false;         // Cf Spectrum from Heidelberg source

bool useCylindricalDetector = false; // detector shape for the specific detectors
bool useSphericalDetector = true;
bool useySheetDetector = false;
bool usexSheetDetector = false;

bool useRealisticModel = false;     // uses the detector response function instead of upper and lower thl
bool useRealisticModelDetector = false; // uses the realistic (response function) model for the virtual detector
bool useRealisticModelLayer = false; // uses the realistic (response function) model for the detector layer
bool useAdditionalDetectorModel = false; // uses a second response function for the realistic detector model
bool useRoverModel = false;
bool useVolumeSource = false;       // source model which extends the source beyound the source layer until the ground layer by an exponential
bool useHECascadeModel = true;      // model which continues the high energetic part in case of absorption in order to reproduce the correct attenuation length
bool godzillaMode = false;          // allows positioning of the source by the mouse cursor in the birds eye view
bool warnUndefinedMaterial = false; // warns if the material code is not found
bool calcNeutronTime = true;        // adds the (non-relativistic) neutron travel time calculation

bool outputCSfilenames = false;      // command line output of cs files

TH1D* cutView = new TH1D();
vector< TObject* > allTHs;         // these are the containers for all the histograms which are used allTHs is for the ones in the GUI, allTHs2 are the internal ones
vector< TObject* > allTHs2;

float trackMetricFactorModifier = 0.6;  // modifiy the step width of the resolution for different material checks for the voxels

int rangeM = 1;
long nTotal = 0;                    // total number of neutrons to be calculated
float sysTemperature = 293.;        // moderator temperature in K

int maxScatterings = 700;           // maximum number of scatterings of a neutron until killed

float xPosSource = 0;               // coordinates and spatial extension for manual positioning of the source
float yPosSource = 0;
float zPosSource = 0;
float xSizeSource = 0;
float ySizeSource = 0;
float zSizeSource = 0;
int zLayerSource = 0;
float radiusSource = 0;
float sourceEnergy = 1;
int sourceDirection = 0;

double domainLowEdge = 0;           // spatial extension of the domain drawn by sideways tracking
double domainUpperEdge = 0;
double domainZLowEdge = 0;
double domainZUpperEdge = 0;
double domainLowDrawEdge = 0;
double domainUpperDrawEdge = 0;
double domainZLowDrawEdge = 0;
double domainZUpperDrawEdge = 0;

bool reflectiveBoundaries = false;
bool periodicBoundaries = false;
bool differentMaterialHit = false;

bool useRadialBeam = false;         // radial or uniform density distribution from center
bool useRectShape = true;           // rectangular (1) or circular shape (0) of the source

float detRad = 9000;                // radius of the detector [mm]
float detPosX[1] = {0};             // coordinates of the detector
float detPosY[1] = {0};
float detLength = 9000;             // length of the detector sheet [mm]

bool newData = false;               // GUI

bool noTrackRecording = false;      // disables track coordinate scoring
bool clearEveryXNeutrons = false;   // sets data to zero every x neutrons
int clearEveryXNeutronsNumber = 100;  // corresponding number of neutrons
bool setAutoRefreshRate = true;     // sets refresh rate automatically to neutron calculation time
bool setAutoRefreshRateClearing = false;  // clears on every display refresh

bool logScaleR = false;             // GUI
bool logScaleD = false;
bool logScaleRS = false;
bool pausehere = false;
bool stopMode = false;

bool islandSetup = false;           // activates the presets of the GUI
bool riverSetup = false;
bool lakeSetup = false;
bool coastSetup = false;

float islandDiameter = 10000;       // geometry definitions for the presets [mm]
float lakeDiameter = 10000;
float riverDiameter = 10000;
float coastPosition = 0;

bool showDensityTrackMap = false;   // GUI: selection of the radio button for the histogram type to plot
bool showDensityIntermediateTrackMap = false;
bool showDensityFastTrackMap = false;
bool showDensityAlbedoTrackMap = false;

bool showDensityEnergyTrackMap = false;

bool showDensityMapIntermediate = false;
bool showDensityMapFast = false;
bool showDensityThermalTrackMap = false;
bool showDensityMapThermal = false;

bool showDetectorOriginMap = true;
bool updateEnlargedView = false;    // GUI
bool updateEnlargedView2 = false;    // GUI

//activates sideways projection (slows down)
bool showDensityTrackMapSide = false;   // activates sideways tracking (computationally expensive)
float densityTrackMapSideCutOutValue = squareDim * 0.5;
//%%%%%%%%%%%%%%%%%%%%%%%%%

bool showDensityMap = false;        // GUI
bool showDensityMapAlbedo = true;

bool detectorAbsorbing = false;     // sets the detector to either being transparent (0) or absorbing for every neutron passed (1)

bool detFileOutput = false;         // activates the neutron coordinate output for the detector
bool detTrackFileOutput = false;    // activate full neutron track output for detector
bool allTrackFileOutput = false;    // activate full neutron track output
bool otherAllTrackFileOutput = false; //export as json file
bool detLayerFileOutput = false;    // activates the neutron coordinate output for the detector layer
bool detLayerTrackFileOutput = true;    // activates the neutron coordinate output for the detector layer
bool saveEveryXNeutrons = false;    // exports everything in the given interval
int saveEveryXNeutronsPower = 5;    // exports all URANOS files every 10^x neutrons

bool stopRunning = false;           // stops the entire calculation

bool useBasicSpectrum = false;

bool alreadyStarted = false;        // GUI

bool integralSliderMoved = false;   // GUI

bool exportSpectrum = false;        // also exports the precalulated spectrum
bool uranosRootOutput = true;       // exports the internal histograms as well

bool exportEpithermalData = false;  // export options for the neutron density data (check boxes in the export tab)
bool exportFastData = false;
bool exportSelectedData = false;
bool exportIntermediateData = false;
bool exportThermalData = false;
bool exportTrackData = false;
bool exportHighResTrackData = false;

float energylowTHL = 0.000001;      // energy threshold (lower/upper) in MeV for the detector and the detector layer
float energyhighTHL = 0.01;

float downwardScotomaAngle = 0;     // angular cutoff values for the detector layer
float downwardAcceptanceAngle = 6.2832; // scotoma = reject neutrons incoing from a specific downward angle, accpetance = accept only neutrons from that angle
bool createSeparateFolderEachExport = false;    // creates for each new export a new folder
bool exportTemporary = false;

int numberPrecalcNeutrons = 6;      // amount of neutrons which go into the precalculated spectrum as a power of 10^x
int refreshCycle = 100;             // GUI update rate
float refreshTime = 1.;             // GUI update time

Float_t relHumidityAir = 0.400000; 	// humidity air, relative to NN and 20°C    // deprecated
Float_t relHumidityAirVar;                                                      // deprecated
Float_t absHumidityAir = 1.;        // humidity air in g/m^3
Float_t absHumidityAirVar;
float pressureFac = 0.9999;         // air pressure relative to NN and 20°C

//defintion of the soil
Float_t soilSolidFrac = 0.50001; Float_t soilSiFrac = 0.75; Float_t soilAlFrac = 0.25;
Float_t soilSolidFracVar = 0.5;

Float_t soilAirFrac = 0.00000;      //0.25;
Float_t soilWaterFrac = 0.10001;    // volumetric soil moisture for soil
Float_t soilWaterFracVar = 0.1003;  // volumetric soil moisture for soil from general settings (slider)

Float_t rPlants = 0.003;            // plant gas density
double rBoronInSoil = 0.;           // boron concentration in soil [ppm mol/cm3]
double rGdInSoil = 0.;
double rFeInSoil = 0., rMnInSoil = 0., rCInSoil = 0., rNaInSoil = 0., rKInSoil = 0., rTiInSoil = 0., rNInSoil = 0.;
double rCelluloseWaterFrac = 1., rCelluloseFrac = 1.;

long neutrons = 500000;              // number of neutrons
long totalActualNeutrons = 0;        // number of total actual neutrons in loop

float atmDensity = 1013.;            // atmospheric density in cm^2/g
float atmPressure = 101300.;         // atmospheric pressure in Pa
float rigidity = 10.;                // cutoff rigidity in GeV

int startingLayer = 0;               // numbers of the starting layer etc according to the geometry stack (starts at 0 here, in GUI at 1)
int groundLayer = 0;
int detectorLayer = 0;
int detectorID = 0;
int detectorOverrideID = 0;
//bool detectorPixelFound = false;

float detectorHeight = 0;

double evaporationFactor = 0.99;     // amount of neutrons evaporated for each absorbed high energy neutron

float sourceLayerHeight = 0;         // only for specific definition of the source
float xCustomPos = 0, yCustomPos = 0; // for godzilla mode

int nDetectedNeutrons = 0;

float paramMin = 0;                   // parameter counters for the batch mode
float paramMax = 121;                 //40 for detector energy analysis

const int maxLayersAllowed = 300;     // maximum amount of layers, GUI can display only ca 22

QPalette* paletteR = new QPalette();  // GUI
QPalette* paletteB = new QPalette();
QPalette* paletteGray = new QPalette();


vector< float* > neutronCoordinates;  // neutron coordinate vector
vector< double* > geometries;         // vector which stores all the geometry layers

bool haveDifferentSoilMoistures = false; // activates the layer maps
bool useImage = false;                // GUI

int inputMatrixPixels = 1;            // size of the matrix for the matrix layer maps definition
//int inputMatrixMaxPixels = 5000;    //deprecated

// inputPicVector stores material definitions inputPicVector2 stores denstity multiplicators and inputPicVector3 porosity definitions
int inputPics[maxLayersAllowed] = { 0 };
int inputPics2[maxLayersAllowed] = { 0 };
int inputPics3[maxLayersAllowed] = { 0 };
int inputPicSizes[maxLayersAllowed] = { 0 };
vector<TMatrixF> inputPicVector;
vector<TMatrixF> inputPicVector2;
vector<TMatrixF> inputPicVector3;

int inputMaterials[maxLayersAllowed] = { 0 };
vector< vector<float> > materialVector;

// stores the positions of further detector layers
int additionalDetectorLayers[maxLayersAllowed] = { 0 };
int totalAdditionalDetectorLayers = 0;
vector< TH2F* > addDetLayerVec;

int matrixX, matrixY;               // x,y position in meters
double  matrixStartX = -50, matrixStartY = -50; //x,y position in meters of the frame of the matrix
float matrixMetricFactor = 0.5;     // pixel conversion factor equals in meter
float mapMetricFactor = 0;          // same as above

vector< TMatrixF > inputMatrixVector; // stores all the layer maps

float posNo = 0;                    // current layer number

vector< double* > geometryStack;    // vector which holds the layer for the global geometry

// exemplary definition of a geometry
double layer0[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer1[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer2[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer3[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer4[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer5[8] = { -squareDim,-squareDim,squareDim,squareDim,0,0 , 20, posNo++ };
double layer6[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer7[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer8[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer9[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer10[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer11[8] = { -squareDim,-squareDim,squareDim,squareDim,0,0 , 20, posNo++ };
double layer12[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer13[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer14[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer15[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer16[8] = { -squareDim,-squareDim,squareDim,squareDim,-0,0 , 11, posNo++ };
double layer17[8] = { -squareDim,-squareDim,squareDim,squareDim,0,0 , 20, posNo++ };
double layer18[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };
double layer19[8] = { -squareDim,-squareDim,squareDim,squareDim,-0.,0. , 11., posNo++ };


//samples... deprecated
double airUp[8] = { -squareDim,-squareDim,squareDim,squareDim,-1000000.,920000. , 11., posNo++ };
double air0[8] = { -squareDim,-squareDim,squareDim,squareDim,-80000.,30000. , 11., posNo++ };
double air1[8] = { -squareDim,-squareDim,squareDim,squareDim,-50000,47500 , 11, posNo++ };
double air2[8] = { -squareDim,-squareDim,squareDim,squareDim,-2500,500 , 11, posNo++ };
double air3[8] = { -squareDim,-squareDim,squareDim,squareDim,-2000,2000 , 11, posNo++ };
double ground[8] = { -squareDim,-squareDim,squareDim,squareDim,0,1600 , 20, posNo++ };

// source defintion for the AmBe source
float spectrumAmBeBins[53] = { 0., 0.11, 0.33, 0.54, 0.75, 0.97, 1.18, 1.4, 1.61, 1.82, 2.04, 2.25, 2.47, 2.68, 2.9, 3.11, 3.32, 3.54, 3.75, 3.97, 4.18, 4.39, 4.61, 4.82, 5.04, 5.25, 5.47, 5.68, 5.98, 6.11, 6.32, 6.54, 6.75, 6.96, 7.18, 7.39, 7.61, 7.82, 8.03, 8.25, 8.46, 8.68, 8.89, 9.11, 9.32, 9.53, 9.75, 9.96, 10.2, 10.4, 10.6, 10.8, 11. };
float spectrumAmBeValues[53] = { 1.44, 3.34, 3.13, 2.81, 2.5, 2.14, 1.98, 1.75, 1.92, 2.23, 2.15, 2.25, 2.28, 2.95, 3.56, 3.69, 3.46, 3.07, 3, 2.69, 2.86, 3.18, 3.07, 3.33, 3.04, 2.74, 2.33, 2.06, 1.82, 1.77, 2.04, 1.83, 1.63, 1.68, 1.68, 1.88, 1.84, 1.69, 1.44, 0.968, 0.652, 0.426, 0.367, 0.381, 0.506, 0.625, 0.552, 0.468, 0.37, 0.278, 0.151, 0.0363, 0 };

// generates the table in the main tab for the geometry
QStandardItemModel* model;
int row = -1;
int column = -1;
int type = 0;
QWidget* boxWindowPointer;
QPushButton* buttonPointer;
bool inputMatrixDefsShown = false;

// variables used throughout the calculation routine and which are made global
double oldDiffTime = 0;
long nTotalOld = 0;
bool newDataComes = false;
string timeRemainString, timeRemainHoursString, timeRemainMinutesString, timeRemainSecondsString;

float fpSoilMoist = 0.1;
float fpHum = 10;
int killedThermal = 0;
double totalRealTime = 0, neutronRealScalingFactor = 0;

int densityMapButtonID = 4;

float scaleFactor = 0;

// for batch run
int paramInt;
float param, paramEnergy, paramDensity;

// color definitions
QColor lightBlue = QColor(215, 237, 255);
QColor someBlue = QColor(0, 88, 156);
QColor someGreen = QColor(86, 200, 90);
QColor someViolet = QColor(156, 39, 176);
QColor stdGray = QColor(240, 240, 240);

/**
 * QT Main window constructor
 *
 * @param parent - object of the QWidget class which is the base class of all user interface objects.

 */
MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    TCanvas* myDummyCanvas = new TCanvas();
    ui->setupUi(this);

    ui->customPlot2->installEventFilter(this);
    //ui->qcustomPlot2->installEventFilter(this);
    //setWindowIcon(QIcon(":/uranos.ico"));

    setupGraph(0);
    ui->groupBox_Evaporation->setHidden(true);
    //ui->checkBoxThermalData->setHidden(true);
    //ui->checkBoxThermalMap->setHidden(true);
    //ui->radioButton_mapThermal->setHidden(true);
    //ui->radioButton_mapTrackThermal->setHidden(true);

    //setConnectionsEnlarge();

    if (doSkyEvaporation)
    {
        ui->groupBox_Evaporation->setHidden(false);
        on_radioButton_fission_clicked();
    }
    if (!noThermalRegime)
    {
    }

    delete myDummyCanvas;
}

void MainWindow::setConnectionsEnlarge()
{
    connect(visualization, SIGNAL(on_VisualizationEnlargeCursorX(float)), this, SLOT(setSourcePosX(float)));
    connect(visualization, SIGNAL(on_VisualizationEnlargeCursorY(float)), this, SLOT(setSourcePosY(float)));

    connect(visualization, SIGNAL(on_VisualizationEnlargePressed(bool)), this, SLOT(setStopMode(bool)));
    QObject::connect(visualization, SIGNAL(on_VisualizationEnlargePressed(bool)), this, SLOT(setContinue(bool)));
    QObject::connect(visualization, SIGNAL(on_VisualizationEnlargePressed(bool)), this, SLOT(setGodzillaMode(bool)));

    QObject::connect(visualization, SIGNAL(on_VisualizationEnlargeReleased(bool)), this, SLOT(setPause(bool)));
    QObject::connect(visualization, SIGNAL(on_VisualizationEnlargeState(bool)), this, SLOT(setGodzillaMode(bool)));
}

/**
 * Generate all the geometry layers for the minimum config
 */
void MainWindow::setGeometry()
{
    geometries.clear();

    geometries.push_back(airUp);
    geometries.push_back(air0);
    geometries.push_back(air1);
    geometries.push_back(air2);
    geometries.push_back(air3);
    geometries.push_back(ground);
}

/**
 * Set text for the status labels in the top bar.
 *
 * @param  numberLabel The number of the labels.
 * @param  msg The text to be set to the label and displayed.
 */
void MainWindow::setStatus(int numberLabel, string msg)
{
    if (!silentMode)
    {
        if (noGUIMode)
        {
            cout << msg << endl;
        }
        else
        {
            if (numberLabel == 1)  ui->label_Status->setText(QString::fromStdString((string)msg));
            if (numberLabel == 2)  ui->label_Status2->setText(QString::fromStdString((string)msg));
        }
    }
}

/**
 * Set up the Geometry in the main tab, also defines the matrixMetricFactor, the resolution of the image.
 *
 */
void MainWindow::setupGeometry()
{
    geometries.clear();

    for (int i = 0; i < model->rowCount(); i++)
    {
        geometryStack.at(i)[4] = 1000. * (model->data(model->index(i, 0, QModelIndex()), Qt::EditRole).toFloat());
        geometryStack.at(i)[5] = 1000. * (model->data(model->index(i, 1, QModelIndex()), Qt::EditRole).toFloat());
        geometryStack.at(i)[6] = (model->data(model->index(i, 2, QModelIndex()), Qt::EditRole).toFloat());

        geometries.push_back(geometryStack.at(i));
    }

    startingLayer = ui->spinBox_StartingLayer->value() - 1;
    detectorLayer = ui->spinBox_DetectorLayer->value() - 1;
    groundLayer = ui->spinBox_GroundLayer->value() - 1;

    if (useImage)
    {
        for (int k = 0; k < model->rowCount(); k++) { if (inputMatrixPixels < inputPicSizes[k]) inputMatrixPixels = inputPicSizes[k]; }
        matrixMetricFactor = squareDim / 1000. / (inputMatrixPixels * 1.); //meter per pixel
    }

    int mat;

    for (int i = 0; i < maxLayersAllowed; i++)
    {
        additionalDetectorLayers[i] = -1;
    }

    totalAdditionalDetectorLayers = 0;

    // checking whether additional detector layers were defined
    for (int i = 0; i < geometries.size() - 1; i++)
    {
        mat = (int)geometries.at(i)[6];

        if ((mat > 500) && (mat < 600))
        {
            additionalDetectorLayers[i] = totalAdditionalDetectorLayers;

            totalAdditionalDetectorLayers++;

            mat = mat - 500;
            geometries.at(i)[6] = mat;
        }
    }
}

// twodimensional functions which store the footprint paramaters according to Schroen 2017
TF2* p01;
TF2* p11;
TF2* p21;
TF2* p31;

TF2* p0;
TF2* p1;
TF2* p2;
TF2* p3;

bool fpReadSuccess = false;

/**
 * Set up the footprints.
 *
 * @return successful - true if footprints were set, false o/w.
 */
bool setupFootprintFunction()
{
    bool successful = false;

    // the names of the root files and the name of the functions within them are hardcoded here
    TFile readParFtPrnt1(endfFolder + "FootprintParams1.root", "READ");
    if (readParFtPrnt1.IsOpen())
    {
        p01 = (TF2*)readParFtPrnt1.Get("grfunction2");
        p11 = (TF2*)readParFtPrnt1.Get("grfunction4");
        p21 = (TF2*)readParFtPrnt1.Get("grfunction33");
        p31 = (TF2*)readParFtPrnt1.Get("grfunction5");
        successful = true;
    }

    TFile readParFtPrnt2(endfFolder + "FootprintParams2.root", "READ");
    if (readParFtPrnt2.IsOpen())
    {
        p0 = (TF2*)readParFtPrnt2.Get("grfunction2");
        p1 = (TF2*)readParFtPrnt2.Get("grfunction4");
        p2 = (TF2*)readParFtPrnt2.Get("grfunction3");
        p3 = (TF2*)readParFtPrnt2.Get("grfunction5");
    }
    if (!successful) cout << "FootprintParams root files not found in the ENDF folder (only required for graphical representation)" << endl;
    return successful;
}

/**
 * Set up the import dialog.
 * importSettings() imports the uranos.cfg and this method configures the GUI afterwards
 *
 */
void MainWindow::setupImport()
{
    if (importSettings())
    {
        ui->lineEdit_InputSpectrumFolder->setText(QString::fromStdString((string)inputSpectrumFile));
        ui->lineEdit_WorkFolder->setText(QString::fromStdString((string)workFolder));
        if (workFolder == "")  ui->lineEdit_WorkFolder->setText(QString::fromStdString((string)"default"));
        ui->lineEdit_DetectorFile->setText(QString::fromStdString((string)detectorResponseFunctionFile));
        if (detectorResponseFunctionFile == "")  ui->lineEdit_DetectorFile->setText(QString::fromStdString((string)"N/A"));
        ui->lineEdit_OutputFolder->setText(QString::fromStdString((string)outputFolder));
        ui->lineEdit_CrosssectionFolder->setText(QString::fromStdString((string)endfFolder));

        ui->lineEditNeutrinsTotal->setText(QString::fromStdString((string)castLongToString(neutrons)));
        ui->lineEditBeamRad->setText(QString::number(beamRadius / 1000.));
        ui->lineEditSquareDim->setText(QString::number(squareDim / 1000.));

        //ui->lineEditBeamRad->setText(QString::fromStdString((string)castFloatToString(beamRadius/1000.)));
        ui->lineEditRefresh->setText(QString::fromStdString((string)castIntToString(refreshCycle)));

        ui->lineEditTHLlow->setText(QString::number(energylowTHL));
        ui->lineEditTHLhigh->setText(QString::number(energyhighTHL));

        ui->lineEditDetRad->setText(QString::number(detRad / 1000.));
        ui->lineEditDetX->setText(QString::number(detPosX[0] / 1000.));
        ui->lineEditDety->setText(QString::number(detPosY[0] / 1000.));
        ui->lineEditDetLength->setText(QString::number(detLength / 1000.));

        if (uranosRootOutput) ui->checkBoxFileOutputPDF_2->setChecked(true);   else  ui->checkBoxFileOutputPDF_2->setChecked(false);
        if (createSeparateFolderEachExport) ui->checkBoxCreateFolder->setChecked(true);     else  ui->checkBoxCreateFolder->setChecked(false);

        string soilSolidFracTextVar = castFloatToString(((1. - soilSolidFracVar) * 100.), 4) + " %";
        int soilSolidFracIntVar = (1 - soilSolidFracVar) * 100.;
        ui->labelSP->setText(QString::fromStdString((string)soilSolidFracTextVar));
        ui->sliderSoilPorosity->setValue(soilSolidFracIntVar);

        soilSolidFracTextVar = castFloatToString((soilWaterFracVar * 100.), 4) + " %";
        soilSolidFracIntVar = soilWaterFracVar * 200.;
        ui->labelSM->setText(QString::fromStdString((string)soilSolidFracTextVar));
        ui->sliderSoilMoisture->setValue(soilSolidFracIntVar);

        //soilSolidFracTextVar = castFloatToString(10./0.6*(relHumidityAir*1.),4)+" g/m<sup>3</sup>";
        soilSolidFracTextVar = castFloatToString(absHumidityAir, 4) + " g/m<sup>3</sup>";
        //soilSolidFracIntVar = relHumidityAir*50.;
        soilSolidFracIntVar = absHumidityAir * 3.;
        ui->labelHum->setText(QString::fromStdString((string)soilSolidFracTextVar));
        ui->sliderAirHum->setValue(soilSolidFracIntVar);

        ui->sliderAtm1->setValue((int)atmDensity);

        soilSolidFracIntVar = rigidity;
        soilSolidFracTextVar = castFloatToString(soilSolidFracIntVar, 4);
        ui->labelrigidity->setText(QString::fromStdString((string)soilSolidFracTextVar));
        ui->sliderRigidity->setValue(soilSolidFracIntVar * 10);

        if (useRadialBeam) { ui->beamRound->setChecked(true);  ui->beamSquare->setChecked(false); ui->label_21->setText("Radius"); }
        if (useRectShape) { ui->beamSquare->setChecked(true);  ui->beamRound->setChecked(false); ui->label_21->setText("1/2 Side Length"); }

        if (useVolumeSource) ui->checkBox_VolumeSource->setChecked(true); else ui->checkBox_VolumeSource->setChecked(false);
        if (useHECascadeModel) ui->checkBox_HEModel->setChecked(true); else ui->checkBox_HEModel->setChecked(false);

        ui->lineEditScotoma->setText(QString::number(downwardScotomaAngle * 2. / 2. / TMath::Pi() * 360.));
        ui->lineEditAntiScotoma->setText(QString::number(downwardAcceptanceAngle * 2. / 2. / TMath::Pi() * 360.));

        if (useDetectorSensitiveMaterial) ui->checkBox->setChecked(true); else ui->checkBox->setChecked(false);
        ui->lineEditScotoma_2->setText(QString::number(detectorSensitiveMaterial));

        if (noMultipleScatteringRecording) ui->checkBox_NoMultipleScattering->setChecked(true); else ui->checkBox_NoMultipleScattering->setChecked(false);
        if (trackAllLayers) ui->checkBox_TrackAllLayers->setChecked(true); else ui->checkBox_TrackAllLayers->setChecked(false);

        if (useRealisticModelLayer) { ui->radioButton_detectorLayerRealistic->setChecked(true);  ui->radioButton_detectorLayerEnergyBand->setChecked(false); }
        else { ui->radioButton_detectorLayerRealistic->setChecked(false);  ui->radioButton_detectorLayerEnergyBand->setChecked(true); }
        if (useRealisticModelDetector) { ui->radioButton_detectorRealistic->setChecked(true);  ui->radioButton_detectorEnergyBand->setChecked(false); }
        else { ui->radioButton_detectorRealistic->setChecked(false);  ui->radioButton_detectorEnergyBand->setChecked(true); }

        if (ui->checkBoxDetectorOriginMap->isChecked()) { showDetectorOriginMap = true; }
        else { showDetectorOriginMap = false; }

        if (useCylindricalDetector) { ui->radioButton_Cylinder->setChecked(true);  ui->radioButton_Sphere->setChecked(false);  ui->radioButton_ySheet->setChecked(false);  ui->radioButton_xSheet->setChecked(false); }
        if (useCylindricalDetector) { ui->lineEditDetLength->setEnabled(false);  ui->lineEditDetRad->setEnabled(true); }
        //else {ui->radioButton_Cylinder->setChecked(false);  ui->radioButton_Sphere->setChecked(true);}

        if (useSphericalDetector) { ui->radioButton_Cylinder->setChecked(false);  ui->radioButton_Sphere->setChecked(true);   ui->radioButton_ySheet->setChecked(false);  ui->radioButton_xSheet->setChecked(false); }
        if (useSphericalDetector) { ui->lineEditDetLength->setEnabled(false);  ui->lineEditDetRad->setEnabled(true); }
        //else {ui->radioButton_Cylinder->setChecked(true);  ui->radioButton_Sphere->setChecked(false);}

        if (useySheetDetector) { ui->radioButton_Cylinder->setChecked(false);  ui->radioButton_Sphere->setChecked(false);   ui->radioButton_ySheet->setChecked(true);  ui->radioButton_xSheet->setChecked(false); }
        if (useySheetDetector) { ui->lineEditDetLength->setEnabled(true);  ui->lineEditDetRad->setEnabled(false); }

        if (usexSheetDetector) { ui->radioButton_Cylinder->setChecked(false);  ui->radioButton_Sphere->setChecked(false);   ui->radioButton_ySheet->setChecked(false);  ui->radioButton_xSheet->setChecked(true); }
        if (usexSheetDetector) { ui->lineEditDetLength->setEnabled(true);  ui->lineEditDetRad->setEnabled(false); }

        if (layerMapsImport)
        {
            on_pushButton_LoadGeometry_clicked();
            on_checkBox_useImage_clicked();
        }

        if (clearEveryXNeutrons) { ui->checkBoxSaveEvery_2->setChecked(true); }
        else { ui->checkBoxSaveEvery_2->setChecked(false); }
        ui->lineEditClearEveryXNeutrons->setText(QString::number(clearEveryXNeutronsNumber));

        ui->lineEdit_AutoUpdate->setText(QString::number(refreshTime));
        ui->checkBoxAutoRefreshRate->setChecked(setAutoRefreshRate);
        ui->checkBoxClearEveryDisplayRefresh->setChecked(setAutoRefreshRateClearing);

        if (true)
        {
            ui->lineEdit_xPos->setText(QString::number(xPosSource / 1000.));
            ui->lineEdit_yPos->setText(QString::number(yPosSource / 1000.));
            ui->lineEdit_zPos->setText(QString::number(zPosSource / 1000.));
            ui->lineEdit_xSize->setText(QString::number(xSizeSource / 1000.));
            ui->lineEdit_ySize->setText(QString::number(ySizeSource / 1000.));
            ui->lineEdit_zPos_2->setText(QString::number(radiusSource / 1000.));
            ui->lineEdit_SourceEnergy->setText(QString::number(sourceEnergy));
        }
        if (exportTrackData) { ui->checkBoxTrackingData->setChecked(true); }
        else { ui->checkBoxTrackingData->setChecked(false); }
        if (exportHighResTrackData) { ui->checkBoxHighResTrackingData->setChecked(true); highResCalc = true;}
        else { ui->checkBoxHighResTrackingData->setChecked(false);  highResCalc = false; }

        if (domainCutoff) ui->checkBox_DomainCutoff->setChecked(true); else ui->checkBox_DomainCutoff->setChecked(false);
        if (domainCutoff2) ui->checkBox_DomainCutoffMeters->setChecked(true); else ui->checkBox_DomainCutoffMeters->setChecked(false);
        ui->lineEditDomainFactor->setText(QString::number(domainCutoffFactor));
        ui->lineEditDomainMeters->setText(QString::number(domainCutoffMeters));

        switch (densityMapButtonID)
        {
        case 1: on_radioButton_map_clicked(); ui->radioButton_map->click(); break;
        case 2: on_radioButton_mapInter_clicked(); ui->radioButton_mapInter->click(); break;
        case 3: on_radioButton_mapFast_clicked(); ui->radioButton_mapFast->click(); break;
        case 4: on_radioButton_mapAlbedo_clicked(); ui->radioButton_mapAlbedo->click(); break;
        case 5: on_radioButton_mapThermal_clicked(); ui->radioButton_mapThermal->click(); break;
        case 10:on_radioButton_mapTrackEnergy_clicked(); ui->radioButton_mapTrackEnergy->click(); break;
        case 11:on_radioButton_mapTrack_clicked(); ui->radioButton_mapTrack->click(); break;
        case 12:on_radioButton_mapTrackInter_clicked(); ui->radioButton_mapTrackInter->click(); break;
        case 13:on_radioButton_mapTrackFast_clicked(); ui->radioButton_mapTrackFast->click(); break;
        case 14:on_radioButton_mapTrackAlbedo_clicked(); ui->radioButton_mapTrackAlbedo->click(); break;
        case 15:on_radioButton_mapTrackThermal_clicked(); ui->radioButton_mapTrackThermal->click(); break;
        }
    }
}


/**
 * Sets up the graphs in the main window and prepares the geometry.
 *
 */
void MainWindow::setupGraph(int index)
{
    setupRunSpectraGraph(ui->customPlot);
    setupRunBirdsEyeViewGraph(ui->customPlot2);
    setupRunHorizSliceXGraph(ui->customPlot3);
    setupRunNpassDetGraph(ui->customPlot4);
    setupRunDepOfIntGraph(ui->customPlot5);
    setupRunNpassDetLyrGraph(ui->customPlot6);
    setupRunRelIntVsEGraph(ui->WidgetSpectrum);
    setupRunOrginOfNGraph(ui->customPlot10);

    setupTable(ui->tableViewLayer);
    setupRangeFunction(ui->customPlotFP);

    paletteR->setColor(QPalette::Text, Qt::red);
    paletteB->setColor(QPalette::Text, Qt::black);
    paletteGray->setColor(QPalette::Text, Qt::gray);

    geometryStack.push_back(layer0);
    geometryStack.push_back(layer1); geometryStack.push_back(layer2); geometryStack.push_back(layer3); geometryStack.push_back(layer4); geometryStack.push_back(layer5);
    geometryStack.push_back(layer6); geometryStack.push_back(layer7); geometryStack.push_back(layer8); geometryStack.push_back(layer9); geometryStack.push_back(layer10);
    geometryStack.push_back(layer11);  geometryStack.push_back(layer12);  geometryStack.push_back(layer13);  geometryStack.push_back(layer14);  geometryStack.push_back(layer15);
    geometryStack.push_back(layer16);  geometryStack.push_back(layer17);  geometryStack.push_back(layer18);  geometryStack.push_back(layer19);

    for (int i = 20; i < maxLayersAllowed; i++)
    {
        geometryStack.push_back(new double[8]);
        geometryStack.at(i)[0] = -squareDim;
        geometryStack.at(i)[1] = -squareDim;
        geometryStack.at(i)[2] = squareDim;
        geometryStack.at(i)[3] = squareDim;
        geometryStack.at(i)[4] = 0;
        geometryStack.at(i)[5] = 0;
        geometryStack.at(i)[6] = 11;
        geometryStack.at(i)[7] = i;
    }

    //initialize everything with zero
    for (int i = 0; i < (::x.size()); i++) { ::x[i] = 0; ::y[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsIncFullSpectrum.size()); i++) { ::plotGUIxBinsIncFullSpectrum[i] = 0; ::plotGUIyBinsIncFullSpectrum[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsIncOnlySpectrum.size()); i++) { ::plotGUIxBinsIncOnlySpectrum[i] = 0; ::plotGUIyBinsIncOnlySpectrum[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsScatteredSurfaceSpec.size()); i++) { ::plotGUIxBinsScatteredSurfaceSpec[i] = 0; ::plotGUIyBinsScatteredSurfaceSpec[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsScatteredSurfaceSpecAlbedo.size()); i++) { ::plotGUIxBinsScatteredSurfaceSpecAlbedo[i] = 0; ::plotGUIyBinsScatteredSurfaceSpecAlbedo[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsCutView.size()); i++) { ::plotGUIxBinsCutView[i] = 0; ::plotGUIyBinsCutView[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsDistanceDet.size()); i++) { ::plotGUIxBinsDistanceDet[i] = 0; ::plotGUIyBinsDistanceDet[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsDistanceAlbedoDet.size()); i++) { ::plotGUIxBinsDistanceAlbedoDet[i] = 0; ::plotGUIyBinsDistanceAlbedoDet[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsDistanceDetLayer.size()); i++) { ::plotGUIxBinsDistanceDetLayer[i] = 0; ::plotGUIyBinsDistanceDetLayer[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsDistanceAlbedoDetLayer.size()); i++) { ::plotGUIxBinsDistanceAlbedoDetLayer[i] = 0; ::plotGUIyBinsDistanceAlbedoDetLayer[i] = 0; }
    for (int i = 0; i < (::xRS.size()); i++) { ::xRS[i] = 0; ::yRS[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsScatDepth.size()); i++) { ::plotGUIxBinsScatDepth[i] = 0; ::plotGUIyBinsScatDepth[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsScatDepthDet.size()); i++) { ::plotGUIxBinsScatDepthDet[i] = 0; ::plotGUIyBinsScatDepthDet[i] = 0; }
    for (int i = 0; i < (::plotGUIxBinsScatDepthDetMax.size()); i++) { ::plotGUIxBinsScatDepthDetMax[i] = 0; ::plotGUIyBinsScatDepthDetMax[i] = 0; }

    //setupImport();
    //ui->customPlot->replot();
}

void MainWindow::getGraphConfig()
{
    setupImport();
    ui->customPlot->replot();
}

/**
 * Set the output format for vector graphic export (pdf).
 *
 */
void MainWindow::formatForVectorGraphics()
{
    gStyle->SetLineWidth(0.1);
    gStyle->SetHistLineWidth(0.1);
    gStyle->SetTextSize(1.);

    //gStyle->SetFuncWidth(0.1);
    //gStyle->SetLabelSize(0.04,"xy");
    //gStyle->SetTitleSize(0.042,"xy");
    //gStyle->SetTitleOffset(1.2,"x");
    //gStyle->SetLabelOffset(0.004,"xy");
}

/**
 * Footprint function according to Schroen 2017.
 *
 * @param x - x value of the function [m]
 * @param par - par[1] = volumetric soil moisture [Vol%*0.01], par[0] = air humidity [g/m3]
 * @return expValue - evaluated function at x.
 */
double rangeFunctionCombined(double* x, double* par)
{
    double expValue;

    TF1* expFuncD1 = new TF1("expFuncD1", "[0]*(exp(-[1]*x))+[3]*exp(-[2]*x)", 1.5e3, 51e3);
    TF1* expFuncD2 = new TF1("expFuncD2", "[0]*(exp(-[1]*x))+[3]*exp(-[2]*x)", 51.e3, 600e3);

    expFuncD2->SetParameters(p0->Eval(par[0], par[1]), p1->Eval(par[0], par[1]), p2->Eval(par[0], par[1]), p3->Eval(par[0], par[1]));
    expFuncD1->SetParameters(p01->Eval(par[0], par[1]), p11->Eval(par[0], par[1]), p21->Eval(par[0], par[1]), p31->Eval(par[0], par[1]));

    if (x[0] * 1000. < 51e3)
    {
        expValue = expFuncD1->Eval(x[0] * 1000.);
    }
    else
    {
        expValue = expFuncD2->Eval(x[0] * 1000.);
    }
    delete expFuncD1;
    delete expFuncD2;

    return expValue / 1000.;
}


/**
 * Set up the footprint function plot.
 *
 */
void MainWindow::setupRangeFunction(QCustomPlot* customPlot)
{
    customPlot->setBackground(stdGray);
    //customPlot->legend->setBrush(QBrush(stdGray));
    customPlot->xAxis->setLabel("Distance [m]");
    customPlot->yAxis->setLabel("Relative Intensity");
}


/**
 * Formatting plots to specific GUI colors.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::formatPlotToColor(QCustomPlot* customPlot)
{
    customPlot->setBackground(someBlue);
    //customPlot->axisRect()->setBackground(Qt::white);
    customPlot->xAxis->setTickLabelColor(lightBlue);
    customPlot->yAxis->setTickLabelColor(lightBlue);
    customPlot->xAxis->setLabelColor(lightBlue);
    customPlot->yAxis->setLabelColor(lightBlue);
    customPlot->xAxis->setTickPen(QPen(lightBlue));
    customPlot->xAxis->setBasePen(QPen(lightBlue));
    customPlot->xAxis->setSubTickPen(QPen(lightBlue));
    customPlot->yAxis->setTickPen(QPen(lightBlue));
    customPlot->yAxis->setBasePen(QPen(lightBlue));
    customPlot->yAxis->setSubTickPen(QPen(lightBlue));

    customPlot->legend->setTextColor(lightBlue);
    customPlot->legend->setBorderPen(QPen(lightBlue));
    customPlot->legend->setBrush(QBrush(someBlue));
}

/**
 * Set up the Energy vs Relative Intensity graph.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunRelIntVsEGraph(QCustomPlot* customPlot)
{
    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
    customPlot->setBackground(stdGray);
    customPlot->legend->setBrush(QBrush(stdGray));
    customPlot->xAxis->setLabel("Energy [MeV]");
    customPlot->yAxis->setLabel("Relative Intensity");
}


/**
 * Set up the Energy vs no of neutrons graph - spectra graph.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunSpectraGraph(QCustomPlot* customPlot)
{
    formatPlotToColor(customPlot);

    // create graph and assign data to it:
    customPlot->addGraph();
    customPlot->addGraph();
    customPlot->addGraph();

    customPlot->graph(0)->setName("Incoming Spectrum");
    customPlot->graph(1)->setName("Surface Spectrum");
    customPlot->graph(2)->setName("Backscattered Spectrum");

    customPlot->graph(0)->setData(::x, ::y); //allDisplayData.push_back(&::x); allDisplayData.push_back(&::y);

    customPlot->graph(0)->setPen(QPen(lightBlue));
    customPlot->graph(1)->setBrush(QBrush(QColor(215, 237, 255, 90)));

    customPlot->graph(1)->setPen(QPen(someGreen));
    customPlot->graph(2)->setPen(QPen(someViolet));

    // give the axes some labels:
    customPlot->xAxis->setLabel("Energy [MeV]");
    customPlot->yAxis->setLabel("n");
    // set axes ranges, so we see all data:
    customPlot->xAxis->setRange(-1, 1);
    customPlot->yAxis->setRange(0, 1);

    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
}


/**
 * Set up the x [m] vs n graph - horizontal slice along the x axis.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunHorizSliceXGraph(QCustomPlot* customPlot)
{
    customPlot->addGraph();

    formatPlotToColor(customPlot);

    customPlot->graph(0)->setData(::plotGUIxBinsCutView, ::plotGUIyBinsCutView);

    customPlot->graph(0)->setPen(QPen(lightBlue));
    customPlot->graph(0)->setBrush(QBrush(QColor(215, 237, 255, 90)));

    // give the axes some labels:
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("n");
}

/**
 * Set up the x [m] vs n graph - Range Distribution of Neutrons Passing the Detector.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunNpassDetGraph(QCustomPlot* customPlot)
{
    customPlot->addGraph();
    customPlot->addGraph();

    customPlot->graph(1)->setName("Albedo Neutrons");
    customPlot->graph(0)->setName("All Neutrons");

    formatPlotToColor(customPlot);

    customPlot->graph(1)->setData(::plotGUIxBinsDistanceAlbedoDet, ::plotGUIyBinsDistanceAlbedoDet);

    customPlot->graph(0)->setData(::plotGUIxBinsDistanceDet, ::plotGUIyBinsDistanceDet);

    customPlot->graph(0)->setBrush(QBrush(QColor(215, 237, 255, 90)));
    customPlot->graph(0)->setPen(QPen(lightBlue));

    customPlot->graph(1)->setPen(QPen(someGreen));

    // give the axes some labels:
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("n");
    // set axes ranges, so we see all data:
    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignRight);
}

/**
 * Set up the depth vs n graph - Depth Of Interactions.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunDepOfIntGraph(QCustomPlot* customPlot)
{
    customPlot->addGraph();
    customPlot->addGraph();
    customPlot->addGraph();

    formatPlotToColor(customPlot);

    customPlot->graph(0)->setName("All Scatter Centers (normalized)");
    customPlot->graph(1)->setName("Albedo: Last Probe Depth");
    customPlot->graph(2)->setName("Albedo: Maximum Probe Depth");

    customPlot->graph(0)->setData(::plotGUIxBinsScatDepth, ::plotGUIyBinsScatDepth);
    customPlot->graph(0)->setBrush(QBrush(QColor(128, 128, 128, 40))); // first graph will be filled with translucent blue

    customPlot->graph(1)->setData(::plotGUIxBinsScatDepthDet, ::plotGUIyBinsScatDepthDet);

    customPlot->graph(2)->setData(::plotGUIxBinsScatDepthDetMax, ::plotGUIyBinsScatDepthDetMax);

    customPlot->graph(0)->setBrush(QBrush(QColor(215, 237, 255, 90)));
    customPlot->graph(0)->setPen(QPen(lightBlue));

    customPlot->graph(1)->setPen(QPen(someGreen));

    customPlot->graph(2)->setPen(QPen(someViolet));

    // give the axes some labels:
    customPlot->xAxis->setLabel("Depth [mm]");
    customPlot->yAxis->setLabel("n");
    // set axes ranges, so we see all data:

    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignRight);
}

/**
 * Set up the x [m] vs n graph  - Range Distribution of Neutrons Passing the Detector Layer.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunNpassDetLyrGraph(QCustomPlot* customPlot)
{
    customPlot->addGraph();
    customPlot->addGraph();

    formatPlotToColor(customPlot);

    customPlot->graph(0)->setName("All Neutrons");
    customPlot->graph(1)->setName("Albedo Neutrons");

    customPlot->graph(0)->setData(::plotGUIxBinsDistanceDetLayer, ::plotGUIyBinsDistanceDetLayer);

    customPlot->graph(1)->setData(::plotGUIxBinsDistanceAlbedoDetLayer, ::plotGUIyBinsDistanceAlbedoDetLayer);

    customPlot->graph(0)->setBrush(QBrush(QColor(215, 237, 255, 90)));
    customPlot->graph(0)->setPen(QPen(lightBlue));

    customPlot->graph(1)->setPen(QPen(someGreen));

    // give the axes some labels:
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("n");

    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignRight);
}

/**
 * Set up the x [m] vs  y [m] graph - Birds Eye View Graph.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunBirdsEyeViewGraph(QCustomPlot* customPlot)
{
    QColor lightBlue = QColor(215, 237, 255);
    QColor someBlue = QColor(0, 88, 156);
    customPlot->setBackground(someBlue);
    customPlot->axisRect()->setBackground(Qt::white);
    customPlot->xAxis->setTickLabelColor(lightBlue);
    customPlot->yAxis->setTickLabelColor(lightBlue);
    customPlot->xAxis->setLabelColor(lightBlue);
    customPlot->yAxis->setLabelColor(lightBlue);
    customPlot->xAxis->setBasePen(QPen(lightBlue));
    customPlot->xAxis->setSubTickPen(QPen(lightBlue));
    customPlot->yAxis->setTickPen(QPen(lightBlue));
    customPlot->yAxis->setBasePen(QPen(lightBlue));
    customPlot->yAxis->setSubTickPen(QPen(lightBlue));

    // initialization pattern (just optics)
    for (int i = 0; i < 501; ++i)
    {
        ::x2[i] = i / 50.0 - 1;
        ::y2[i] = i / 50.0 - 1;
    }
    QCPColorMap* colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("y [m]");

    colorMap->data()->setSize(250, 250);
    colorMap->data()->setRange(QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.), QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.));

    for (int x = 0; x < 250; ++x)
    {
        for (int y = 0; y < 250; ++y)
        {
            colorMap->data()->setCell(x, y, r.Rndm() * TMath::Abs((x - 125) * (y - 125)));
        }
    }

    colorMap->setGradient(QCPColorGradient::gpJet);
    colorMap->rescaleDataRange(true);
    QCPRange dataR = colorMap->dataRange(); dataR *= 2;
    colorMap->setDataRange(dataR);
    customPlot->rescaleAxes();
    customPlot->replot();
}


/**
 * Set up the x [m] vs  y [m] graph - Orgins Of Neutrons measured by the detector.
 *
 * @param customPlot - An object of QCustomPlot, which is a Qt C++ widget for plotting and data visualization.
 */
void MainWindow::setupRunOrginOfNGraph(QCustomPlot* customPlot)
{
    QColor lightBlue = QColor(215, 237, 255);
    QColor someBlue = QColor(0, 88, 156);
    customPlot->setBackground(someBlue);
    customPlot->axisRect()->setBackground(Qt::white);
    customPlot->xAxis->setTickLabelColor(lightBlue);
    customPlot->yAxis->setTickLabelColor(lightBlue);
    customPlot->xAxis->setLabelColor(lightBlue);
    customPlot->yAxis->setLabelColor(lightBlue);
    customPlot->xAxis->setBasePen(QPen(lightBlue));
    customPlot->xAxis->setSubTickPen(QPen(lightBlue));
    customPlot->yAxis->setTickPen(QPen(lightBlue));
    customPlot->yAxis->setBasePen(QPen(lightBlue));
    customPlot->yAxis->setSubTickPen(QPen(lightBlue));

    for (int i = 0; i < 501; ++i)
    {
        //::x2[i] = i/50.0 - 1;
        //::y2[i] = i/50.0 - 1;
    }

    QCPColorMap* colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("y [m]");

    colorMap->data()->setSize(250, 250);
    colorMap->data()->setRange(QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.), QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.));

    for (int x = 0; x < 250; ++x)
    {
        for (int y = 0; y < 250; ++y)
        {
            colorMap->data()->setCell(x, y, r.Rndm() * TMath::Abs((x - 125) * (y - 125)));
        }
    }
    colorMap->setGradient(QCPColorGradient::gpCold);
    colorMap->rescaleDataRange(true);
    QCPRange dataR = colorMap->dataRange(); dataR *= 2;
    colorMap->setDataRange(dataR);

    QCPLayoutGrid* subLayout = new QCPLayoutGrid;
    ui->customPlot10->plotLayout()->addElement(1, 0, subLayout);
    subLayout->setMargins(QMargins(53, 0, 16, 5));

    QCPColorScale* colorScale = new QCPColorScale(customPlot);
    subLayout->addElement(0, 0, colorScale);
    colorScale->setType(QCPAxis::atBottom);
    colorScale->axis()->setSubTickPen(QPen(lightBlue));
    colorScale->axis()->setTickPen(QPen(lightBlue));
    colorScale->axis()->setTickLabelColor(lightBlue);
    colorScale->axis()->setLabelColor(lightBlue);

    customPlot->rescaleAxes();
    customPlot->replot();
}


/**
 * Set the Integer Type 1
 * currently deprecated
 */
void MainWindow::setIntType1()
{
    type = 1;
    boxWindowPointer->setWindowTitle("_");
    boxWindowPointer->close();
    buttonPointer->setText("Air Layer");
}

/**
 * Set the Integer Type 2
 * currently deprecated
 */
void MainWindow::setIntType2()
{
    type = 2;
    boxWindowPointer->setWindowTitle("_");
    boxWindowPointer->close();
    buttonPointer->setText("Detector Layer");
}

/**
 * Set the Integer Type 3
 * currently deprecated
 */
void MainWindow::setIntType3()
{
    type = 3;
    //boxWindowPointer = widget;
    boxWindowPointer->setWindowTitle("_");
    boxWindowPointer->close();
    buttonPointer->setText("Ground Layer");
}

QDialogButtonBox* buttonBoxPointer;
vector< QPushButton* > buttonVector;

/**
 * MainWindow inherited function for button type configuration for the simulation layers
 * currently deprecated
 * @param -idx
 */
QPushButton* MainWindow::typebutton(QModelIndex idx)
{
    QPushButton* button = new QPushButton();
    buttonPointer = button;

    QDialog* subDialog = new QDialog;
    subDialog->setWindowTitle("Layer Type");

    QWidget* boxWindow = new QWidget;
    //boxWindow->setWindowFlags( Qt::CustomizeWindowHint );
    QVBoxLayout* boxLayout = new QVBoxLayout;
    boxWindowPointer = boxWindow;

    QPushButton* buttonA1 = new QPushButton();
    buttonA1->setText("Air Layer");
    QPushButton* buttonA2 = new QPushButton();
    buttonA2->setText("Detector Layer");
    QPushButton* buttonA3 = new QPushButton();
    buttonA3->setText("Ground Layer");

    QDialogButtonBox* buttonBox = new QDialogButtonBox(Qt::Vertical);
    buttonBoxPointer = buttonBox;
    boxWindow->setGeometry(640, 280, 100, 120);
    buttonBox->addButton(buttonA1, QDialogButtonBox::AcceptRole);
    buttonBox->addButton(buttonA2, QDialogButtonBox::AcceptRole);
    buttonBox->addButton(buttonA3, QDialogButtonBox::AcceptRole);

    boxLayout->addWidget(buttonBox);
    boxWindow->setLayout(boxLayout);

    QObject::connect(button, SIGNAL(clicked()), boxWindow, SLOT(show()));
    // QObject::connect(button, SIGNAL(clicked()), this, SLOT(buttonPointer=boxWindow; ));

    QObject::connect(buttonA1, SIGNAL(clicked()), this, SLOT(setIntType1()));
    QObject::connect(buttonA2, SIGNAL(clicked()), this, SLOT(setIntType2()));
    QObject::connect(buttonA3, SIGNAL(clicked()), this, SLOT(setIntType3()));

    //button->setCheckable(true);
    //QObject::connect(button, SIGNAL(pressed()), ui->tableViewLayer, SLOT(update() ));

    QObject::connect(boxWindow, SIGNAL(windowTitleChanged(QString)), button, SLOT(show()));

    buttonVector.push_back(button);

    return button;
}

/**
 * Set Focus of the mouse for the setupTable function (field which is clicked)
 * @param idx
 */
void MainWindow::setFocus(const QModelIndex& idx)
{
    row = idx.row();
    column = idx.column();
}

/**
 * Set Focus of the mouse for the setupTable function (field which is clicked)
 * @param idx
 */
void MainWindow::setRow(int rowNumber)
{
    row = rowNumber;
}

/**
 * Set up the geometry table for the layer configuration
 * @param table
 */
void MainWindow::setupTable(QTableView* table)
{
    model = new QStandardItemModel(5, 1, ui->tabWidget);

    model->setHorizontalHeaderItem(0, new QStandardItem(QString("Position")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QString("Height")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QString("Material")));
    //model->setHorizontalHeaderItem(3, new QStandardItem(QString("Type")));
    model->setHorizontalHeaderItem(3, new QStandardItem(QString("Matrix")));

    QHeaderView* verticalHeader = table->verticalHeader();
    verticalHeader->sectionResizeMode(QHeaderView::Fixed);
    verticalHeader->setDefaultSectionSize(22);

    table->setModel(model);

    table->setColumnWidth(0, 90);
    table->setColumnWidth(1, 80);
    table->setColumnWidth(2, 60);
    table->setColumnWidth(3, 120);


    QModelIndex id;
    id = model->index(0, 0, QModelIndex());

    connect(table, SIGNAL(clicked(const QModelIndex&)), this, SLOT(setFocus(const QModelIndex&)));
    //connect(table, SIGNAL(currentChanged(const QModelIndex&)), this, SLOT(setFocus(const QModelIndex&)));
    connect(verticalHeader, SIGNAL(sectionClicked(int)), this, SLOT(setRow(int)));

    QSortFilterProxyModel* proxy1 = new QSortFilterProxyModel();
    proxy1->setSourceModel(model);

    for (int rowItr = 0; rowItr < proxy1->rowCount(); ++rowItr)
    {
        for (int colItr = 0; colItr < proxy1->columnCount(); ++colItr)
        {
            //QStandardItem *item= new QStandardItem();
            //item->setText(proxy1->index(z,y).data().toString());
            QStandardItem* item = model->itemFromIndex(model->index(rowItr, colItr, QModelIndex()));

            item->setTextAlignment(Qt::AlignCenter);
            //model->setItem(z,y,item);
        }
    }

    table->setAlternatingRowColors(true);
    //table->resizeRowsToContents();
    //table->resizeColumnsToContents();

    int totalRows = model->rowCount();
    for (int modelRowItr = 1; modelRowItr < totalRows; modelRowItr++)
    {
        model->removeRow(1, QModelIndex());
    }
}

/**
 * Set up the LiveView (Root) Histograms (Histograms which are shown in the GUI)
 *
 */
void setupLiveTHs()
{
    // The histograms need to be defined globally. However, these are just pointers.
    // In order to reinitialize the histograms, the pointers/objects need to be deleted first and then recreated
    delete detectorDistance;
    delete scatteredSurfaceDistance;
    delete detectorDistanceBackScattered;
    delete detectorLayerDistance;
    delete detectorLayerDistanceBackscattered;
    delete detectorOriginMap;

    delete densityTrackMapHighRes;
    delete densityIntermediateTrackMapHighRes;
    delete densityFastTrackMapHighRes;
    delete densityAlbedoTrackMapHighRes;
    delete densityHighEnergyTrackMapHighRes;
    delete densityEnergyTrackMapHighRes;

    delete densityTrackMapHighRes15x;
    delete densityIntermediateTrackMapHighRes15x;
    delete densityFastTrackMapHighRes15x;
    delete densityAlbedoTrackMapHighRes15x;
    delete densityHighEnergyTrackMapHighRes15x;
    delete densityEnergyTrackMapHighRes15x;

    delete densityFastTrackMapHighRes2x;

    delete densityThermalTrackMapHighRes;
    delete densityThermalTrackMapHighRes15x;

    delete densityThermalTrackMap;
    delete densityMapThermal;

    delete densityTrackMap;
    delete densityIntermediateTrackMap;
    delete densityFastTrackMap;
    delete densityAlbedoTrackMap;
    delete densityMap;
    delete densityMapIntermediate;
    delete densityMapFast;
    delete densityMapAlbedo;
    delete densityHighEnergyTrackMap;
    delete densityMapHighEnergy;
    delete densityTrackMapSide;
    delete densityTrackMapSideAlbedo;
    delete densityTrackMapSideDetector;
    delete densityTrackMapSideThermal;
    delete densityEnergyTrackMap;

    domainLowEdge = -squareDim * 0.5 / 1000.;
    domainUpperEdge = squareDim * 0.5 / 1000.;

    densityTrackMapHighRes = new TH2F("densityTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);
    densityIntermediateTrackMapHighRes = new TH2F("densityIntermediateTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);
    densityFastTrackMapHighRes = new TH2F("densityFastTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);
    densityAlbedoTrackMapHighRes = new TH2F("densityAlbedoTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);
    densityHighEnergyTrackMapHighRes = new TH2F("densityHighEnergyTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);

    densityEnergyTrackMapHighRes = new TH2F("densityEnergyTrackMapHighRes", "Detector Layer Neutron Energy Weighted Tracks", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);

    densityTrackMapHighRes15x = new TH2F("densityTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);
    densityIntermediateTrackMapHighRes15x = new TH2F("densityIntermediateTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);
    densityFastTrackMapHighRes15x = new TH2F("densityFastTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);
    densityAlbedoTrackMapHighRes15x = new TH2F("densityAlbedoTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);
    densityHighEnergyTrackMapHighRes15x = new TH2F("densityHighEnergyTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);

    densityEnergyTrackMapHighRes15x = new TH2F("densityEnergyTrackMapHighRes15x", "Detector Layer Neutron Energy Weighted Tracks", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);

    densityFastTrackMapHighRes2x = new TH2F("densityFastTrackMapHighRes2x", "Detector Layer Neutron Track Density", 2000, domainLowEdge, domainUpperEdge, 2000, domainLowEdge, domainUpperEdge);

    densityThermalTrackMapHighRes = new TH2F("densityThermalTrackMapHighRes", "Detector Layer Neutron Track Density", 1000, domainLowEdge, domainUpperEdge, 1000, domainLowEdge, domainUpperEdge);
    densityThermalTrackMapHighRes15x = new TH2F("densityThermalTrackMapHighRes15x", "Detector Layer Neutron Track Density", 1500, domainLowEdge, domainUpperEdge, 1500, domainLowEdge, domainUpperEdge);

    densityTrackMap = new TH2F("densityTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);
    densityIntermediateTrackMap = new TH2F("densityIntermediateTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);
    densityFastTrackMap = new TH2F("densityFastTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);
    densityAlbedoTrackMap = new TH2F("densityAlbedoTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);
    densityHighEnergyTrackMap = new TH2F("densityHighEnergyTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);

    densityEnergyTrackMap = new TH2F("densityEnergyTrackMap", "Detector Layer Neutron Energy Weighted Tracks", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);

    densityThermalTrackMap = new TH2F("densityThermalTrackMap", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);
    densityMapThermal = new TH2F("densityMapThermal", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);

    densityMap = new TH2F("densityMap", "Detector Layer Neutron Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);				//TH2Fashion(densityMap, "x", "y", setSize);
    densityMapIntermediate = new TH2F("densityMapIntermediate", "Intermediate Energy Neutron Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);			//	TH2Fashion(densityMapIntermediate, "x", "y", setSize);
    densityMapFast = new TH2F("densityMapFast", "Fast Neutron Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);				//TH2Fashion(densityMapFast, "x", "y", setSize);
    densityMapAlbedo = new TH2F("densityMapAlbedo", "Detector Layer Albedo Neutron Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);			//	TH2Fashion(densityMapAlbedo, "x", "y", setSize);
    densityMapHighEnergy = new TH2F("densityMapHighEnergy", "Detector Layer Neutron Track Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);

    detectorDistance = new TH1F("detectorDistance", "Distance Distribution for Detector", 4000, -100, domainUpperEdge * 1000.);
    scatteredSurfaceDistance = new TH1F("scatteredSurfaceDistance", "Distance Distribution on Surface", 2000, -100, domainUpperEdge * 1000.);
    detectorOriginMap = new TH2F("detectorOriginMap", "Detector Neutron Origins", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge);

    detectorDistanceBackScattered = new TH1F("detectorDistanceBackScattered", "Distance Distribution Backscattered for Detector", 4000, -100, domainUpperEdge * 1000.);	//TH1Fashion(detectorDistanceBackScattered, "n", "Distance [mm]", setSize);

    detectorLayerDistance = new TH1F("detectorLayerDistance", "Distance Distribution for Detector Layer", 4000, -100, domainUpperEdge * 1000.);
    detectorLayerDistanceBackscattered = new TH1F("detectorLayerDistanceBackscattered", "Albedo Distance Distribution for Detector Layer", 4000, -100, domainUpperEdge * 1000.);

    //domainZLowEdge = -(geometries.at(geometries.size()-1)[4]+geometries.at(geometries.size()-1)[5])/1000.;
    domainZLowEdge = -(geometries.at(geometries.size() - 1)[4] + 0 * geometries.at(geometries.size() - 1)[5]) / 1000.;
    //domainZUpperEdge = -(geometries.at(startingLayer)[4]+geometries.at(startingLayer)[5])/1000.;
    domainZUpperEdge = -(geometries.at(0)[4] + geometries.at(0)[5]) / 1000.;

    densityTrackMapSide = new TH2F("densityTrackMapSide", "Side Neutron Track Density", 2000, domainLowEdge, domainUpperEdge, 2000, domainZLowEdge, domainZUpperEdge);
    densityTrackMapSideAlbedo = new TH2F("densityTrackMapSideAlbedo", "Side Neutron Albedo Track Density", 2000, domainLowEdge, domainUpperEdge, 2000, domainZLowEdge, domainZUpperEdge);
    densityTrackMapSideDetector = new TH2F("densityTrackMapSideDetector", "Side Neutron Detector Track Density", 2000, domainLowEdge, domainUpperEdge, 2000, domainZLowEdge, domainZUpperEdge);
    densityTrackMapSideThermal = new TH2F("densityTrackMapSideThermal", "Side Neutron Thermal Track Density", 2000, domainLowEdge, domainUpperEdge, 2000, domainZLowEdge, domainZUpperEdge);

    for(int i = 0; i < addDetLayerVec.size(); i++)
    {
       delete  addDetLayerVec.at(i);
    }

    addDetLayerVec.clear();

    for(int i = 0; i < totalAdditionalDetectorLayers; i++)
    {
        TString addLayerName = "densityMapAlbedo_"+castIntToString(i);
        addDetLayerVec.push_back(new TH2F(addLayerName, "Detector Layer Neutron Density", 500, domainLowEdge, domainUpperEdge, 500, domainLowEdge, domainUpperEdge));
    }

    liveTHs.clear();

    liveTHs.push_back(scatDepth); liveTHs.push_back(scatteredSurfaceDistance); liveTHs.push_back(detectorDistance); liveTHs.push_back(detectorDistanceBackScattered);
    liveTHs.push_back(detectorLayerDistance); liveTHs.push_back(detectorLayerDistanceBackscattered);
    liveTHs.push_back(scatteredSurfaceSpectrumBack); liveTHs.push_back(scatteredSurfaceSpectrum); liveTHs.push_back(cosmicSpectrum);
    liveTHs.push_back(densityMapAlbedo); liveTHs.push_back(densityMapFast); liveTHs.push_back(densityMapIntermediate); liveTHs.push_back(densityMap); liveTHs.push_back(densityMapHighEnergy);
    liveTHs.push_back(densityTrackMap); liveTHs.push_back(densityIntermediateTrackMap); liveTHs.push_back(densityFastTrackMap); liveTHs.push_back(densityAlbedoTrackMap); liveTHs.push_back(densityHighEnergyTrackMap);
    liveTHs.push_back(densityEnergyTrackMap);
    liveTHs.push_back(densityMapThermal); liveTHs.push_back(densityThermalTrackMap); liveTHs.push_back(densityThermalTrackMapHighRes); liveTHs.push_back(densityThermalTrackMapHighRes15x);
    liveTHs.push_back(densityTrackMapHighRes); liveTHs.push_back(densityIntermediateTrackMapHighRes); liveTHs.push_back(densityFastTrackMapHighRes); liveTHs.push_back(densityFastTrackMapHighRes2x); liveTHs.push_back(densityAlbedoTrackMapHighRes); liveTHs.push_back(densityHighEnergyTrackMapHighRes);
    liveTHs.push_back(densityEnergyTrackMapHighRes);
    liveTHs.push_back(densityTrackMapHighRes15x); liveTHs.push_back(densityIntermediateTrackMapHighRes15x); liveTHs.push_back(densityFastTrackMapHighRes15x); liveTHs.push_back(densityAlbedoTrackMapHighRes15x); liveTHs.push_back(densityHighEnergyTrackMapHighRes15x);
    liveTHs.push_back(densityEnergyTrackMapHighRes15x);

    liveTHs.push_back(scatteredSurfaceDepth); liveTHs.push_back(scatteredSurfaceMaxDepth);
    liveTHs.push_back(detectorOriginMap);

    if (showDensityTrackMapSide)
    {
        liveTHs.push_back(densityTrackMapSide);
        liveTHs.push_back(densityTrackMapSideAlbedo);
        liveTHs.push_back(densityTrackMapSideDetector);
        liveTHs.push_back(densityTrackMapSideThermal);
    }
}

/**
 * waiting routine.
 *
 * @param millisecondsToWait - No of milli seconds to wait.
 */
void delay(int millisecondsToWait)
{
    QTime dieTime = QTime::currentTime().addMSecs(millisecondsToWait);
    while (QTime::currentTime() < dieTime)
    {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
    }
}

/**
 * redraws the high resolution view of the Bird's Eyes density plot
 *
 */
void MainWindow::redrawEnlargedView()
{
    visualization->setColorScheme(0);
    if (ui->radioButton_NeutronNight->isChecked()) { visualization->setColorScheme(1); }
    if (ui->radioButton_NeutronCold->isChecked()) { visualization->setColorScheme(2); }
    if (ui->radioButton_NeutronPolar->isChecked()) { visualization->setColorScheme(3); }
    if (ui->radioButton_NeutronHot->isChecked()) { visualization->setColorScheme(4); }
    if (ui->radioButton_NeutronThermal->isChecked()) { visualization->setColorScheme(5); }
    if (ui->radioButton_NeutronGrayScale->isChecked()) { visualization->setColorScheme(6); }

    visualization->setmanualColorZero(manualColorZero);
    visualization->setmanualColor(manualColor);
    visualization->setuseManualColor(useManualColors);
    visualization->setsilderColorMoved(silderColorMoved);
    visualization->sethorizontalSliderColorZeroValue(ui->horizontalSliderColorZero->value());
    visualization->sethorizontalSliderValue(ui->horizontalSliderColor->value());
    visualization->setplotTopViewLog(ui->checkBoxLogarithmic->isChecked());

    if (showDensityMapThermal) { visualization->plotGraph(densityMapThermal, 500, squareDim); }

    if (showDensityMap) { visualization->plotGraph(densityMap, 500, squareDim); }
    if (showDensityMapIntermediate) { visualization->plotGraph(densityMapIntermediate, 500, squareDim); }
    if (showDensityMapFast) { visualization->plotGraph(densityMapFast, 500, squareDim); }
    if (showDensityMapAlbedo) { visualization->plotGraph(densityMapAlbedo, 500, squareDim); }

    if (showDensityTrackMap) { visualization->plotGraph(densityTrackMapHighRes, 1000, squareDim); }
    if (showDensityIntermediateTrackMap) { visualization->plotGraph(densityIntermediateTrackMapHighRes, 1000, squareDim); }
    if (showDensityFastTrackMap) { visualization->plotGraph(densityFastTrackMapHighRes, 1000, squareDim); }
    if (showDensityAlbedoTrackMap) { visualization->plotGraph(densityAlbedoTrackMapHighRes, 1000, squareDim); }

    if (showDensityThermalTrackMap) { visualization->plotGraph(densityThermalTrackMapHighRes, 1000, squareDim); }

    if (showDensityEnergyTrackMap) { visualization->plotGraph(densityEnergyTrackMapHighRes, 1000, squareDim); }
}

void MainWindow::redrawEnlargedView2()
{
    visualization2->setColorScheme2(0);
    if (ui->radioButton_NeutronNight->isChecked()) { visualization2->setColorScheme2(1); }
    if (ui->radioButton_NeutronCold->isChecked()) { visualization2->setColorScheme2(2); }
    if (ui->radioButton_NeutronPolar->isChecked()) { visualization2->setColorScheme2(3); }
    if (ui->radioButton_NeutronHot->isChecked()) { visualization2->setColorScheme2(4); }
    if (ui->radioButton_NeutronThermal->isChecked()) { visualization2->setColorScheme2(5); }
    if (ui->radioButton_NeutronGrayScale->isChecked()) { visualization2->setColorScheme2(6); }

    visualization2->setmanualColorZero2(manualColorZero);
    visualization2->setmanualColor2(manualColor);
    visualization2->setuseManualColor2(useManualColors);
    visualization2->setsilderColorMoved2(silderColorMoved);
    visualization2->sethorizontalSliderColorZeroValue2(ui->horizontalSliderColorZero->value());
    visualization2->sethorizontalSliderValue2(ui->horizontalSliderColor->value());
    visualization2->setplotTopViewLog2(ui->checkBoxLogarithmic->isChecked());

    if (showDensityMapThermal) { visualization2->plotGraph2(densityMapThermal, 500, squareDim); }

    if (showDensityMap) { visualization2->plotGraph2(densityMap, 500, squareDim); }
    if (showDensityMapIntermediate) { visualization2->plotGraph2(densityMapIntermediate, 500, squareDim); }
    if (showDensityMapFast) { visualization2->plotGraph2(densityMapFast, 500, squareDim); }
    if (showDensityMapAlbedo) { visualization2->plotGraph2(densityMapAlbedo, 500, squareDim); }

    if (showDensityTrackMap) { visualization2->plotGraph2(densityTrackMapHighRes15x, 1500, squareDim); }
    if (showDensityIntermediateTrackMap) { visualization2->plotGraph2(densityIntermediateTrackMapHighRes15x, 1500, squareDim); }
    if (showDensityFastTrackMap) { visualization2->plotGraph2(densityFastTrackMapHighRes15x, 1500, squareDim); }
    if (showDensityAlbedoTrackMap) { visualization2->plotGraph2(densityAlbedoTrackMapHighRes15x, 1500, squareDim); }

    if (showDensityThermalTrackMap) { visualization2->plotGraph2(densityThermalTrackMapHighRes15x, 1500, squareDim); }

    if (showDensityEnergyTrackMap) { visualization2->plotGraph2(densityEnergyTrackMapHighRes15x, 1500, squareDim); }
}


double entryWeight, maximumWeight;


/**
 * redraws all graphs in the GUI.
 *
 * @param difftime - time since last refresh (in order to calculate neutrons per second).
 */
void  MainWindow::redrawNeutronMap(double difftime)
{
    TCanvas* myDummyCanvas = new TCanvas();

    // detector tab
    if ((ui->tabWidget_live->currentIndex() == 3) || (!simulationRunning))
    {
        int elementCount = ui->customPlot10->plotLayout()->elementCount();

        for (int elementItr = 0; elementItr < elementCount; elementItr++)
        {
            if (qobject_cast<QCPLayoutGrid*>(ui->customPlot10->plotLayout()->elementAt(elementItr)))  ui->customPlot10->plotLayout()->removeAt(elementItr);
        }
        ui->customPlot10->plotLayout()->simplify();

        QCPLayoutGrid* subLayout = new QCPLayoutGrid;
        ui->customPlot10->plotLayout()->addElement(1, 0, subLayout);
        subLayout->setMargins(QMargins(53, 0, 16, 5));
        QCPColorScale* colorScale = new QCPColorScale(ui->customPlot10);
        subLayout->addElement(0, 0, colorScale);
        colorScale->setType(QCPAxis::atBottom);

        QColor lightBlue = QColor(215, 237, 255);
        colorScale->axis()->setSubTickPen(QPen(lightBlue));
        colorScale->axis()->setTickPen(QPen(lightBlue));
        colorScale->axis()->setTickLabelColor(lightBlue);
        colorScale->axis()->setLabelColor(lightBlue);

        ui->customPlot10->clearPlottables();

        float sliderDetectorValue = ui->horizontalSliderDetector->value();
        float rangeValue = 2. + sliderDetectorValue / 200. * 500.; //in m

        const int unitCells = 30;

        if (ui->horizontalSliderDetector->value() <= 1)
        {
            QCPColorMap* colorMapDetector = new QCPColorMap(ui->customPlot10->xAxis, ui->customPlot10->yAxis);

            colorMapDetector->data()->setSize(250, 250);
            colorMapDetector->data()->setRange(QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001), QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001));

            for (int keyIdxItr = 0; keyIdxItr < 500; keyIdxItr += 2)
            {
                for (int valIdxItr = 0; valIdxItr < 500; valIdxItr += 2)
                {
                    colorMapDetector->data()->setCell((int)keyIdxItr / 2, (int)valIdxItr / 2, ::detectorOriginMap->GetBinContent(keyIdxItr, valIdxItr) + ::detectorOriginMap->GetBinContent(keyIdxItr + 1, valIdxItr) + ::detectorOriginMap->GetBinContent(keyIdxItr, valIdxItr + 1) + ::detectorOriginMap->GetBinContent(keyIdxItr + 1, valIdxItr + 1));
                }
            }

            colorMapDetector->setDataRange(QCPRange(0, 2.2));
            colorMapDetector->setGradient(QCPColorGradient::gpCold);

            string detectorLabelMaxText = "-";
            string detectorLabelText = "-";

            ui->label_DetectorRangeMaximum->setText(QString::fromStdString(detectorLabelMaxText));
            ui->label_DetectorRange->setText(QString::fromStdString(detectorLabelText));

            colorScale->setDataRange(QCPRange(0, 2.2));
            colorScale->setGradient(QCPColorGradient::gpCold);
        }
        else
        {
            QCPColorMap* colorMapDetectorRelativeSquare = new QCPColorMap(ui->customPlot10->xAxis, ui->customPlot10->yAxis);

            colorMapDetectorRelativeSquare->data()->setSize(250, 250);
            colorMapDetectorRelativeSquare->data()->setRange(QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001), QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001));

            float allDetectorNeutrons = detectorOriginMap->GetEntries();
            float detectorOriginMapValue;
            float detectorOriginMapDrawValue, detectorOriginMapDrawValueMax = 0;
            bool mapValueFound = false;

            float pixX, pixY, sumMatrixX, sumMatrixY, sumMatrixXp1, sumMatrixYp1;

            TMatrixF sumMatrix(unitCells, unitCells);

            for (int x = 0; x < 500; x += 1)
            {
                pixX = (-squareDim * 0.5 + squareDim * 0.002 * x - detPosX[0]) * 0.001; //relative to Det in m
                for (int y = 0; y < 500; y += 1)
                {
                    detectorOriginMapValue = detectorOriginMap->GetBinContent(x, y);
                    if (detectorOriginMapValue < 0.1) continue;

                    pixY = (-squareDim * 0.5 + squareDim * 0.002 * y - detPosY[0]) * 0.001;

                    for (int m = 0; m < unitCells; m += 1)
                    {
                        mapValueFound = false;
                        sumMatrixX = -rangeValue + 2. * rangeValue * m / (unitCells * 1.);
                        sumMatrixXp1 = sumMatrixX + 2. * rangeValue / (unitCells * 1.);
                        if ((pixX >= sumMatrixX) && (pixX < sumMatrixXp1))
                        {
                            for (int n = 0; n < unitCells; n += 1)
                            {
                                sumMatrixY = -rangeValue + 2. * rangeValue * n / (unitCells * 1.);
                                sumMatrixYp1 = sumMatrixY + 2. * rangeValue / (unitCells * 1.);
                                if ((pixY >= sumMatrixY) && (pixY < sumMatrixYp1))
                                {
                                    sumMatrix(m, n) = sumMatrix(m, n) + detectorOriginMapValue;
                                    mapValueFound = true;
                                    break;
                                }
                            }
                        }
                        if (mapValueFound) break;
                    }
                }
            }

            for (int x = 0; x < 250; x += 1)
            {
                pixX = (-squareDim * 0.5 + squareDim * 0.004 * x - detPosX[0]) * 0.001; //relative to Det in m
                for (int y = 0; y < 250; y += 1)
                {
                    detectorOriginMapDrawValue = 0;

                    pixY = (-squareDim * 0.5 + squareDim * 0.004 * y - detPosY[0]) * 0.001;

                    for (int m = 0; m < unitCells; m += 1)
                    {
                        mapValueFound = false;
                        sumMatrixX = -rangeValue + 2. * rangeValue * m / (unitCells * 1.);
                        sumMatrixXp1 = sumMatrixX + 2. * rangeValue / (unitCells * 1.);
                        if ((pixX >= sumMatrixX) && (pixX < sumMatrixXp1))
                        {
                            for (int n = 0; n < unitCells; n += 1)
                            {
                                sumMatrixY = -rangeValue + 2. * rangeValue * n / (unitCells * 1.);
                                sumMatrixYp1 = sumMatrixY + 2. * rangeValue / (unitCells * 1.);
                                if ((pixY >= sumMatrixY) && (pixY < sumMatrixYp1))
                                {
                                    detectorOriginMapDrawValue = sumMatrix(m, n);
                                    if (detectorOriginMapDrawValue > detectorOriginMapDrawValueMax) detectorOriginMapDrawValueMax = detectorOriginMapDrawValue;
                                    mapValueFound = true;
                                    break;
                                }
                            }
                        }
                        if (mapValueFound) break;
                    }
                    colorMapDetectorRelativeSquare->data()->setCell((int)x, (int)y, detectorOriginMapDrawValue / allDetectorNeutrons);
                }
            }
            float detectorDataScaleRange = 1.05 * (0.01 + ui->horizontalSliderDetectorColor->value() / 100. * detectorOriginMapDrawValueMax / allDetectorNeutrons);
            colorMapDetectorRelativeSquare->setDataRange(QCPRange(0, detectorDataScaleRange));
            colorMapDetectorRelativeSquare->setGradient(QCPColorGradient::gpJet);

            int detectorLabelMaxVal = detectorDataScaleRange * 100;
            int detectorLabelTextVal = rangeValue;
            string detectorLabelMaxText = (string)castIntToString(detectorLabelMaxVal) + " %";
            string detectorLabelText = (string)castIntToString(detectorLabelTextVal) + " m";

            ui->label_DetectorRangeMaximum->setText(QString::fromStdString(detectorLabelMaxText));
            ui->label_DetectorRange->setText(QString::fromStdString(detectorLabelText));

            colorScale->setDataRange(QCPRange(0, 100 * detectorDataScaleRange));
            colorScale->setGradient(QCPColorGradient::gpJet);
        }

        ui->customPlot10->rescaleAxes();
        ui->customPlot10->update();
        ui->customPlot10->replot();
    }

    // bird's eye view tab
    if ((ui->tabWidget_live->currentIndex() == 0) || (!simulationRunning))
    {
        ui->customPlot2->clearPlottables();

        QCPColorMap* colorMap = new QCPColorMap(ui->customPlot2->xAxis, ui->customPlot2->yAxis);

        colorMap->data()->setSize(250, 250);
        colorMap->data()->setRange(QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.), QCPRange(-::squareDim / 2. / 1000., ::squareDim / 2. / 1000.));

        if (plotTopViewLog) { colorMap->setDataScaleType(QCPAxis::stLogarithmic); }

        double allEntries = 0;
        long allEntriesInt = 0;
        double allEntriesMaximum = 0;
        float actualEntry = 0;
        float minZ = 0;

        if (plotTopViewLog) minZ = 1;

        if (true)
        {
            if (ui->checkBoxGradient2->isChecked())
                //if (false)
            {
                TH2F* densityMapCopy = new TH2F();
                TH2F* densityMapRed = new TH2F();

                if (showDensityMap) { densityMapRed = (TH2F*)densityMap->Clone(""); densityMapCopy = (TH2F*)densityMap->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }
                if (showDensityMapIntermediate) { densityMapRed = (TH2F*)densityMapIntermediate->Clone(""); densityMapCopy = (TH2F*)densityMapIntermediate->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }
                if (showDensityMapAlbedo) { densityMapRed = (TH2F*)densityMapAlbedo->Clone(""); densityMapCopy = (TH2F*)densityMapAlbedo->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }

                for (int x = 0; x < 250; x += 1)
                {
                    for (int y = 0; y < 250; y += 1)
                    {
                        actualEntry = densityMapCopy->GetBinContent(x, y);
                        colorMap->data()->setCell((int)x, (int)y, actualEntry);
                        allEntries += actualEntry;
                        if (actualEntry > allEntriesMaximum) allEntriesMaximum = actualEntry;
                    }
                }

                entryWeight = (allEntries) / 250. / 250. * 2.; maximumWeight = 2. * allEntriesMaximum;

                delete densityMapCopy;
                delete densityMapRed;
            }
            else
            {
                for (int x = 0; x < 500; x += 2)
                {
                    for (int y = 0; y < 500; y += 2)
                    {
                        if (showDensityTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityTrackMap->GetBinContent(x, y) + ::densityTrackMap->GetBinContent(x + 1, y) + ::densityTrackMap->GetBinContent(x, y + 1) + ::densityTrackMap->GetBinContent(x + 1, y + 1));
                        if (showDensityIntermediateTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityIntermediateTrackMap->GetBinContent(x, y) + ::densityIntermediateTrackMap->GetBinContent(x + 1, y) + ::densityIntermediateTrackMap->GetBinContent(x, y + 1) + ::densityIntermediateTrackMap->GetBinContent(x + 1, y + 1));
                        if (showDensityFastTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityFastTrackMap->GetBinContent(x, y) + ::densityFastTrackMap->GetBinContent(x + 1, y) + ::densityFastTrackMap->GetBinContent(x, y + 1) + ::densityFastTrackMap->GetBinContent(x + 1, y + 1));
                        if (showDensityAlbedoTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityAlbedoTrackMap->GetBinContent(x, y) + ::densityAlbedoTrackMap->GetBinContent(x + 1, y) + ::densityAlbedoTrackMap->GetBinContent(x, y + 1) + ::densityAlbedoTrackMap->GetBinContent(x + 1, y + 1));

                        if (showDensityEnergyTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityEnergyTrackMap->GetBinContent(x, y) + ::densityEnergyTrackMap->GetBinContent(x + 1, y) + ::densityEnergyTrackMap->GetBinContent(x, y + 1) + ::densityEnergyTrackMap->GetBinContent(x + 1, y + 1));

                        if (showDensityThermalTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityThermalTrackMap->GetBinContent(x, y) + ::densityThermalTrackMap->GetBinContent(x + 1, y) + ::densityThermalTrackMap->GetBinContent(x, y + 1) + ::densityThermalTrackMap->GetBinContent(x + 1, y + 1));
                        if (showDensityMapThermal)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapThermal->GetBinContent(x, y) + ::densityMapThermal->GetBinContent(x + 1, y) + ::densityMapThermal->GetBinContent(x, y + 1) + ::densityMapThermal->GetBinContent(x + 1, y + 1));


                        if (showDensityMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMap->GetBinContent(x, y) + ::densityMap->GetBinContent(x + 1, y) + ::densityMap->GetBinContent(x, y + 1) + ::densityMap->GetBinContent(x + 1, y + 1));
                        if (showDensityMapIntermediate)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapIntermediate->GetBinContent(x, y) + ::densityMapIntermediate->GetBinContent(x + 1, y) + ::densityMapIntermediate->GetBinContent(x, y + 1) + ::densityMapIntermediate->GetBinContent(x + 1, y + 1));
                        if (showDensityMapFast)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapFast->GetBinContent(x, y) + ::densityMapFast->GetBinContent(x + 1, y) + ::densityMapFast->GetBinContent(x, y + 1) + ::densityMapFast->GetBinContent(x + 1, y + 1));
                        if (showDensityMapAlbedo)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapAlbedo->GetBinContent(x, y) + ::densityMapAlbedo->GetBinContent(x + 1, y) + ::densityMapAlbedo->GetBinContent(x, y + 1) + ::densityMapAlbedo->GetBinContent(x + 1, y + 1));

                        if (showDensityMapAlbedo)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapAlbedo->GetBinContent(x, y) + ::densityMapAlbedo->GetBinContent(x + 1, y) + ::densityMapAlbedo->GetBinContent(x, y + 1) + ::densityMapAlbedo->GetBinContent(x + 1, y + 1));
                    }
                }

                if (showDensityTrackMap) { entryWeight = (densityTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityTrackMap->GetBinContent(densityTrackMap->GetMaximumBin()); }
                if (showDensityIntermediateTrackMap) { entryWeight = (densityIntermediateTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityIntermediateTrackMap->GetBinContent(densityIntermediateTrackMap->GetMaximumBin()); }
                if (showDensityFastTrackMap) { entryWeight = (densityFastTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityFastTrackMap->GetBinContent(densityFastTrackMap->GetMaximumBin()); }
                if (showDensityAlbedoTrackMap) { entryWeight = (densityAlbedoTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityAlbedoTrackMap->GetBinContent(densityAlbedoTrackMap->GetMaximumBin()); }

                if (showDensityEnergyTrackMap) { entryWeight = (densityEnergyTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityEnergyTrackMap->GetBinContent(densityEnergyTrackMap->GetMaximumBin()); }

                if (showDensityThermalTrackMap) { entryWeight = (densityThermalTrackMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityThermalTrackMap->GetBinContent(densityThermalTrackMap->GetMaximumBin()); }
                if (showDensityMapThermal) { entryWeight = (densityMapThermal->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityMapThermal->GetBinContent(densityMapThermal->GetMaximumBin()); }

                if (showDensityMap) { entryWeight = (densityMap->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityMap->GetBinContent(densityMap->GetMaximumBin()); }
                if (showDensityMapIntermediate) { entryWeight = (densityMapIntermediate->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityMapIntermediate->GetBinContent(densityMapIntermediate->GetMaximumBin()); }
                if (showDensityMapFast) { entryWeight = (densityMapFast->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityMapFast->GetBinContent(densityMapFast->GetMaximumBin()); }
                if (showDensityMapAlbedo) { entryWeight = (densityMapAlbedo->GetEntries()) / (250. * 250.) * 4.; maximumWeight = 4. * densityMapAlbedo->GetBinContent(densityMapAlbedo->GetMaximumBin()); }

                /*
                allEntriesInt = 0;
                for (int x=120; x<380; x+=1)
                {
                    for (int y=120; y<380; y+=1)
                    {
                        allEntriesInt+=  ::densityMapAlbedo->GetBinContent(x,y);
                    }
                }
                */
                allEntriesInt = entryWeight * 250. * 250. * 4.;

                totalRealTime = (allEntriesInt*1.)/neutronRealScalingFactor;
                string numberofDetectorLayerNs;
                if (totalRealTime < 1)  numberofDetectorLayerNs = castDoubleToString(totalRealTime*1000.,6)+" ms";
                else  numberofDetectorLayerNs = castDoubleToString(totalRealTime,6)+" s";
                ui->label_detectorLayerNs->setText(QString::fromStdString(numberofDetectorLayerNs));

                allEntriesInt = densityMapHighEnergy->GetEntries();
                if (useExtraCounter) ui->label_detectorLayerNs2->setText(QString::fromStdString((string)castLongToString(allEntriesInt)));
            }
        }

        float colorScaleMax = (ui->horizontalSliderColor->value() * 1.) / 200. * maximumWeight;
        if (colorScaleMax < 1) colorScaleMax = 1;

        colorMap->setGradient(QCPColorGradient::gpJet);

        if (ui->radioButton_NeutronNight->isChecked()) { colorMap->setGradient(QCPColorGradient::gpNight); }
        if (ui->radioButton_NeutronCold->isChecked()) { colorMap->setGradient(QCPColorGradient::gpCold); }
        if (ui->radioButton_NeutronPolar->isChecked()) { colorMap->setGradient(QCPColorGradient::gpThermal); }
        if (ui->radioButton_NeutronHot->isChecked()) { colorMap->setGradient(QCPColorGradient::gpHot); }
        if (ui->radioButton_NeutronThermal->isChecked()) { colorMap->setGradient(QCPColorGradient::gpPolar); }
        if (ui->radioButton_NeutronGrayScale->isChecked()) { colorMap->setGradient(QCPColorGradient::gpGrayscale); }

        if ((!silderColorMoved) && ((!ui->checkBoxClearEveryDisplayRefresh->isChecked()) || (!ui->checkBoxSaveEvery_2->isChecked())))
        {
            if (entryWeight * 1.1 < 5) colorMap->setDataRange(QCPRange(minZ, 5));
            else colorMap->setDataRange(QCPRange(minZ, entryWeight * 1.1));
        }
        else
        {
            if (ui->horizontalSliderColorZero->value() == 0) colorMap->setDataRange(QCPRange(minZ, colorScaleMax));
            else colorMap->setDataRange(QCPRange((ui->horizontalSliderColorZero->value() * 1.) / 200. * colorScaleMax, colorScaleMax));
        }

        if (manualColorZero < minZ) manualColorZero = minZ;

        int lowColorValue = (ui->horizontalSliderColorZero->value() * 1.) / 200. * colorScaleMax;
        int highColorValue = colorScaleMax;

        if (ui->checkBoxManual->isChecked())
        {
            useManualColors = true;
            colorMap->setDataRange(QCPRange(manualColorZero, manualColor));
        }
        else
        {
            useManualColors = false;
            ui->lineEditManualColorZero->setText(QString::fromStdString((string)castIntToString(lowColorValue)));
            ui->lineEditManualColor->setText(QString::fromStdString((string)castIntToString(highColorValue)));
        }

        string colorLabelText = "["+(string)castIntToString(lowColorValue)+"-"+(string)castIntToString(highColorValue)+"]";
        ui->label_ColorRange->setText(QString::fromStdString(colorLabelText));

        //colorMap->rescaleDataRange(true);
        ui->customPlot2->rescaleAxes(); // colorMap
        ui->customPlot2->replot();

        if (visualization != 0x0)
        {
            if (visualization->isVisible()) updateEnlargedView = true;
            else  updateEnlargedView = false;
        }

        if (visualization2 != 0x0)
        {
            if (visualization2->isVisible()) updateEnlargedView2 = true;
            else  updateEnlargedView2 = false;
        }

        if (updateEnlargedView)
        {
            redrawEnlargedView();
        }

        if (updateEnlargedView2)
        {
            redrawEnlargedView2();
        }

        // spectrum plots
        ui->customPlot->graph(0)->data()->clear();

        bool rebingraph = false;
        int j = 0;
        for (int i = 100; i < 999; ++i)
        {
            if (!rebingraph)
            {
                ::plotGUIxBinsIncOnlySpectrum[i - 100] = cosmicSpectrum->GetBinLowEdge(i + 1);
                ::plotGUIyBinsIncOnlySpectrum[i - 100] = cosmicSpectrum->GetBinContent(i + 1);
                ::plotGUIxBinsScatteredSurfaceSpec[i - 100] = scatteredSurfaceSpectrum->GetBinLowEdge(i + 1);
                ::plotGUIyBinsScatteredSurfaceSpec[i - 100] = scatteredSurfaceSpectrum->GetBinContent(i + 1);
                ::plotGUIxBinsScatteredSurfaceSpecAlbedo[i - 100] = scatteredSurfaceSpectrumBack->GetBinLowEdge(i + 1);
                ::plotGUIyBinsScatteredSurfaceSpecAlbedo[i - 100] = scatteredSurfaceSpectrumBack->GetBinContent(i + 1);
            }
            else
            {
                ::plotGUIxBinsIncOnlySpectrum[j] = cosmicSpectrum->GetBinLowEdge(i + 1);
                ::plotGUIyBinsIncOnlySpectrum[j] = cosmicSpectrum->GetBinContent(i + 1) + cosmicSpectrum->GetBinContent(i + 2);
                ::plotGUIxBinsScatteredSurfaceSpec[j] = scatteredSurfaceSpectrum->GetBinLowEdge(i + 1);
                ::plotGUIyBinsScatteredSurfaceSpec[j] = scatteredSurfaceSpectrum->GetBinContent(i + 1) + scatteredSurfaceSpectrum->GetBinContent(i + 2);
                ::plotGUIxBinsScatteredSurfaceSpecAlbedo[j] = scatteredSurfaceSpectrumBack->GetBinLowEdge(i + 1);
                ::plotGUIyBinsScatteredSurfaceSpecAlbedo[j] = scatteredSurfaceSpectrumBack->GetBinContent(i + 1) + scatteredSurfaceSpectrumBack->GetBinContent(i + 2);
                j++;
                i++;
            }
        }
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        ui->customPlot->xAxis->setTicker(logTicker);
        ui->customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot->graph(0)->setData(::plotGUIxBinsIncOnlySpectrum, ::plotGUIyBinsIncOnlySpectrum);
        ui->customPlot->graph(1)->setData(::plotGUIxBinsScatteredSurfaceSpec, ::plotGUIyBinsScatteredSurfaceSpec);
        ui->customPlot->graph(2)->setData(::plotGUIxBinsScatteredSurfaceSpecAlbedo, ::plotGUIyBinsScatteredSurfaceSpecAlbedo);
        ui->customPlot->update();
        ui->customPlot->rescaleAxes();
        ui->customPlot->replot();
    }

    if ((ui->tabWidget_live->currentIndex() == 2) || (!simulationRunning))
    {
        ui->customPlot3->graph(0)->data()->clear();  // cut view

        if (showDensityTrackMap)               cutView = densityTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityIntermediateTrackMap)   cutView = densityIntermediateTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityFastTrackMap)           cutView = densityFastTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityAlbedoTrackMap)         cutView = densityAlbedoTrackMap->ProjectionX("proj", 240, 260);

        if (showDensityEnergyTrackMap)         cutView = densityEnergyTrackMap->ProjectionX("proj", 240, 260);

        if (showDensityThermalTrackMap)        cutView = densityThermalTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityMapThermal)             cutView = densityMapThermal->ProjectionX("proj", 240, 260);

        if (showDensityMap)                    cutView = densityMap->ProjectionX("proj", 240, 260);
        if (showDensityMapIntermediate)        cutView = densityMapIntermediate->ProjectionX("proj", 240, 260);
        if (showDensityMapFast)                cutView = densityMapFast->ProjectionX("proj", 240, 260);
        if (showDensityMapAlbedo)              cutView = densityMapAlbedo->ProjectionX("proj", 240, 260);


        for (int i = 0; (i < cutView->GetNbinsX()) && (i < 500); ++i)
        {
            ::plotGUIxBinsCutView[i] = cutView->GetBinCenter(i + 1);
            ::plotGUIyBinsCutView[i] = cutView->GetBinContent(i + 1);
        }

        ui->customPlot3->graph(0)->setData(::plotGUIxBinsCutView, ::plotGUIyBinsCutView);
        ui->customPlot3->update();
        ui->customPlot3->rescaleAxes();
        ui->customPlot3->replot();
    }

    int allDetectorNeutrons = 0;
    int allDetectorAlbedoNeutrons = 0;

    if ((ui->tabWidget_live->currentIndex() == 1) || (!simulationRunning))
    {
        ui->customPlot4->graph(0)->data()->clear(); // detector range distribution
        ui->customPlot4->graph(1)->data()->clear();

        for (int i = 0; i < 4000 - 1; ++i)
        {
            ::plotGUIxBinsDistanceAlbedoDet[i] = detectorDistanceBackScattered->GetBinLowEdge(i + 1) * 0.001;
            ::plotGUIyBinsDistanceAlbedoDet[i] = detectorDistanceBackScattered->GetBinContent(i + 1);
            ::plotGUIxBinsDistanceDet[i] = detectorDistance->GetBinLowEdge(i + 1) * 0.001;
            ::plotGUIyBinsDistanceDet[i] = detectorDistance->GetBinContent(i + 1);
            allDetectorNeutrons += plotGUIyBinsDistanceDet[i];
            allDetectorAlbedoNeutrons += plotGUIyBinsDistanceAlbedoDet[i];
        }

        ui->customPlot4->graph(1)->setData(::plotGUIxBinsDistanceAlbedoDet, ::plotGUIyBinsDistanceAlbedoDet);
        ui->customPlot4->graph(0)->setData(::plotGUIxBinsDistanceDet, ::plotGUIyBinsDistanceDet);

        float neutronDetectionFraction = 100. * allDetectorAlbedoNeutrons / (1. * allDetectorNeutrons);
        ui->detectorRatioBar->setValue(neutronDetectionFraction);

        if (logScaleR)
        {
            ui->customPlot4->yAxis->setScaleType(QCPAxis::stLogarithmic);
            QSharedPointer<QCPAxisTickerLog> logTicker4(new QCPAxisTickerLog);
            ui->customPlot4->yAxis->setTicker(logTicker4);
        }
        else
        {
            ui->customPlot4->yAxis->setScaleType(QCPAxis::stLinear);

            QSharedPointer<QCPAxisTickerFixed> linTicker4(new QCPAxisTickerFixed);
            linTicker4->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
            linTicker4->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

            ui->customPlot4->yAxis->setTicker(linTicker4);
        }

        ui->customPlot4->update();
        ui->customPlot4->rescaleAxes();
        ui->customPlot4->replot();
    }

    if ((ui->tabWidget_live->currentIndex() == 2) || (!simulationRunning))
    {
        ui->customPlot5->graph(0)->data()->clear();  // depth of interaction
        ui->customPlot5->graph(1)->data()->clear();

        float maxDepthCount = 0;
        float maxDepthSurfaceCount = 0;
        for (int i = 0; i < 20; ++i)
        {
            if (scatDepth->GetBinContent(i + 1) > maxDepthCount) maxDepthCount = scatDepth->GetBinContent(i + 1);
            if (scatteredSurfaceDepth->GetBinContent(i + 1) > maxDepthSurfaceCount) maxDepthSurfaceCount = scatteredSurfaceDepth->GetBinContent(i + 1);
        }

        for (int i = 0; i < 1500 - 1; ++i)
        {
            ::plotGUIxBinsScatDepth[i] = scatDepth->GetBinLowEdge(i + 1);
            ::plotGUIyBinsScatDepth[i] = scatDepth->GetBinContent(i + 1) / maxDepthCount * maxDepthSurfaceCount * 0.75;
            ::plotGUIxBinsScatDepthDet[i] = scatteredSurfaceDepth->GetBinLowEdge(i + 1);
            ::plotGUIyBinsScatDepthDet[i] = scatteredSurfaceDepth->GetBinContent(i + 1);
            ::plotGUIxBinsScatDepthDetMax[i] = scatteredSurfaceMaxDepth->GetBinLowEdge(i + 1);
            ::plotGUIyBinsScatDepthDetMax[i] = scatteredSurfaceMaxDepth->GetBinContent(i + 1);
        }

        ui->customPlot5->graph(0)->setData(::plotGUIxBinsScatDepth, ::plotGUIyBinsScatDepth);

        ui->customPlot5->graph(1)->setData(::plotGUIxBinsScatDepthDet, ::plotGUIyBinsScatDepthDet);

        ui->customPlot5->graph(2)->setData(::plotGUIxBinsScatDepthDetMax, ::plotGUIyBinsScatDepthDetMax);

        if (logScaleD)
        {
            QSharedPointer<QCPAxisTickerLog> logTicker5(new QCPAxisTickerLog);
            ui->customPlot5->yAxis->setTicker(logTicker5);
            ui->customPlot5->yAxis->setScaleType(QCPAxis::stLogarithmic);
        }
        else
        {
            ui->customPlot5->yAxis->setScaleType(QCPAxis::stLinear);

            QSharedPointer<QCPAxisTickerFixed> linTicker5(new QCPAxisTickerFixed);
            linTicker5->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
            linTicker5->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

            ui->customPlot5->yAxis->setTicker(linTicker5);
        }

        ui->customPlot5->update();
        ui->customPlot5->rescaleAxes();
        ui->customPlot5->replot();
    }

    if ((ui->tabWidget_live->currentIndex() == 1) || (!simulationRunning))
    {
        ui->customPlot6->graph(0)->data()->clear();  //detector layer distance

        int allDetectorLayerNeutrons = 0;
        int allDetectorLayerAlbedoNeutrons = 0;
        for (int i = 0; i < 4000 - 1; ++i)
        {
            ::plotGUIxBinsDistanceDetLayer[i] = detectorLayerDistance->GetBinLowEdge(i + 1) * 0.001;
            ::plotGUIyBinsDistanceDetLayer[i] = detectorLayerDistance->GetBinContent(i + 1);
            ::plotGUIxBinsDistanceAlbedoDetLayer[i] = detectorLayerDistanceBackscattered->GetBinLowEdge(i + 1) * 0.001;
            ::plotGUIyBinsDistanceAlbedoDetLayer[i] = detectorLayerDistanceBackscattered->GetBinContent(i + 1);
            allDetectorLayerNeutrons += plotGUIyBinsDistanceDetLayer[i];
            allDetectorLayerAlbedoNeutrons += plotGUIyBinsDistanceAlbedoDetLayer[i];
        }

        ui->customPlot6->graph(0)->setData(::plotGUIxBinsDistanceDetLayer, ::plotGUIyBinsDistanceDetLayer);
        ui->customPlot6->graph(1)->setData(::plotGUIxBinsDistanceAlbedoDetLayer, ::plotGUIyBinsDistanceAlbedoDetLayer);
        if (logScaleRS)
        {
            QSharedPointer<QCPAxisTickerLog> logTicker6(new QCPAxisTickerLog);
            ui->customPlot6->yAxis->setTicker(logTicker6);
            ui->customPlot6->yAxis->setScaleType(QCPAxis::stLogarithmic);
        }
        else
        {
            ui->customPlot6->yAxis->setScaleType(QCPAxis::stLinear);

            QSharedPointer<QCPAxisTickerFixed> linTicker6(new QCPAxisTickerFixed);
            linTicker6->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
            linTicker6->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

            ui->customPlot6->yAxis->setTicker(linTicker6);
        }

        float neutronLayerDetectionFraction = 100. * allDetectorLayerAlbedoNeutrons / (1. * allDetectorLayerNeutrons);
        ui->detectorLayerRatioBar->setValue(neutronLayerDetectionFraction);

        if (true)
        {
            int detectorFootprintSum = 0, detectorLayerFootprintSum = 0;
            int detectorFootprint63 = 0, detectorFootprint86 = 0;
            int detectorLayerFootprint63 = 0, detectorLayerFootprint86 = 0;
            int nearField = 0, nearFieldLayer = 0;

            for (int i = 0; i < 4000; i++)
            {
                detectorFootprintSum += plotGUIyBinsDistanceAlbedoDet[i];
                detectorLayerFootprintSum += plotGUIyBinsDistanceAlbedoDetLayer[i];
                if ((detectorFootprintSum > 0.6321 * allDetectorAlbedoNeutrons) && (detectorFootprint63 < 0.1)) { detectorFootprint63 = 0.001 * detectorDistanceBackScattered->GetBinCenter(i + 1); }
                else {}
                if ((detectorFootprintSum > 0.8650 * allDetectorAlbedoNeutrons) && (detectorFootprint86 < 0.1)) { detectorFootprint86 = 0.001 * detectorDistanceBackScattered->GetBinCenter(i + 1); }
                else {}
                if ((detectorLayerFootprintSum > 0.6321 * allDetectorLayerAlbedoNeutrons) && (detectorLayerFootprint63 < 0.1)) { detectorLayerFootprint63 = 0.001 * detectorLayerDistanceBackscattered->GetBinCenter(i + 1); }
                else {}
                if ((detectorLayerFootprintSum > 0.8650 * allDetectorLayerAlbedoNeutrons) && (detectorLayerFootprint86 < 0.1)) { detectorLayerFootprint86 = 0.001 * detectorLayerDistanceBackscattered->GetBinCenter(i + 1); }
                else {}
                if ((plotGUIxBinsDistanceAlbedoDet[i] > 15) && (nearField < 0.1)) { nearField = 100. * detectorFootprintSum / (1. * allDetectorAlbedoNeutrons); }
                if ((plotGUIxBinsDistanceAlbedoDetLayer[i] > 15) && (nearFieldLayer < 0.1)) { nearFieldLayer = 100. * detectorLayerFootprintSum / (1. * allDetectorLayerAlbedoNeutrons); }

                if ((detectorFootprint63 > 1) && (detectorFootprint86 > 1) && (detectorLayerFootprint63 > 1) && (detectorLayerFootprint86 > 1)) break;
            }
            string detectorFootprint63s = (string)castIntToString(detectorFootprint63) + " m";
            string detectorFootprint86s = (string)castIntToString(detectorFootprint86) + " m";
            string detectorLayerFootprint63s = (string)castIntToString(detectorLayerFootprint63) + " m";
            string detectorLayerFootprint86s = (string)castIntToString(detectorLayerFootprint86) + " m";
            string nearFields = (string)castIntToString(nearField) + " %";
            string nearFieldLayers = (string)castIntToString(nearFieldLayer) + " %";

            ui->label_Footprint63D->setText(QString::fromStdString(detectorFootprint63s));
            ui->label_Footprint86D->setText(QString::fromStdString(detectorFootprint86s));
            ui->label_Footprint63DL->setText(QString::fromStdString(detectorLayerFootprint63s));
            ui->label_Footprint86DL->setText(QString::fromStdString(detectorLayerFootprint86s));
            ui->label_FootprintNFL->setText(QString::fromStdString(nearFieldLayers));
            ui->label_FootprintNF->setText(QString::fromStdString(nearFields));
        }

        ui->customPlot6->update();
        ui->customPlot6->rescaleAxes();
        ui->customPlot6->replot();
    }

    if (difftime > 0)
    {
        int npers = (nTotal * 1. / difftime);
        //int npers = ((nTotal-nTotalOld)*1./difftime);
        string numberofNperS = "(" + (string)castIntToString(npers) + "/s)";
        ui->label_npers->setText(QString::fromStdString(numberofNperS));

        if (ui->checkBoxAutoRefreshRate->isChecked())
        {
            if (((refreshTime * npers) / refreshCycle > 1.1) || ((refreshTime * npers) / refreshCycle < 0.9))
            {
                refreshCycle = refreshTime * npers; if (refreshCycle < 10) refreshCycle = 10;
                ui->lineEditRefresh->setText(QString::fromStdString((string)castIntToString(refreshCycle)));
            }
        }

        int timeRemainHours = ((neutrons - nTotal) / npers) / 3600;
        int timeRemainMinutes = (((neutrons - nTotal) / npers) - timeRemainHours * 3600) / 60;
        int timeRemainSeconds = ((neutrons - nTotal) / npers) - timeRemainHours * 3600 - timeRemainMinutes * 60;

        if (timeRemainHours < 10) timeRemainHoursString = "0" + (string)castIntToString(timeRemainHours);
        else  timeRemainHoursString = (string)castIntToString(timeRemainHours);
        if (timeRemainMinutes < 10) timeRemainMinutesString = "0" + (string)castIntToString(timeRemainMinutes);
        else  timeRemainMinutesString = (string)castIntToString(timeRemainMinutes);
        if (timeRemainSeconds < 10) timeRemainSecondsString = "0" + (string)castIntToString(timeRemainSeconds);
        else  timeRemainSecondsString = (string)castIntToString(timeRemainSeconds);

        timeRemainString = "-" + timeRemainHoursString + ":" + timeRemainMinutesString + ":" + timeRemainSecondsString;
        ui->label_timeRemain->setText(QString::fromStdString(timeRemainString));

        oldDiffTime = difftime;
        nTotalOld = nTotal;
    }

    string numberofDetectedNs = (string)castIntToString(nDetectedNeutrons);
    ui->label_detectorNs->setText(QString::fromStdString(numberofDetectedNs));

    if (simulationRunning)
    {
        //string numberofN = "Neutrons: "+castIntToString(nTotal);
        string numberofN = (string)castLongToString(nTotal);
        ui->neutronCountView->setText(QString::fromStdString(numberofN));

        int progress = 100. * (nTotal * 1.) / (neutrons * 1.);
        ui->progressBar->setValue(progress);

    }
    ui->progressBar->repaint();

    newDataComes = false;

    delete myDummyCanvas;

    delay(10);
}

/**
 * import settings from config file and saves the settings to the respective variables.
 *
 * @return - true if the file is not empty , false o/w.
 */
bool MainWindow::importSettings()
{
    TString curLine;
    string fileName = "Uranos.cfg";
    int lineCounter = 0, tempInt = 0;

    if ((noGUIMode) || (configFilePathConfigured)) fileName = configFilePath;

    ifstream input_stream(fileName, ios::in);
    while (curLine.ReadLine(input_stream))
    {
        istrstream stream(curLine.Data());

        if (lineCounter == 0) stream >> outputFolder;
        if (lineCounter == 1) { stream >> workFolder; if (workFolder == "default") workFolder = ""; }
        if (lineCounter == 2) stream >> inputSpectrumFile;
        if (lineCounter == 3) { stream >> detectorResponseFunctionFile; if ((detectorResponseFunctionFile == "default") || (detectorResponseFunctionFile == "N/A") || (detectorResponseFunctionFile == "n/a")) detectorResponseFunctionFile = ""; }
        if (lineCounter == 4) stream >> endfFolder;

        if (lineCounter == 5) stream >> neutrons;
        if (lineCounter == 6) stream >> squareDim;
        if (lineCounter == 7) stream >> beamRadius;

        if (lineCounter == 8) stream >> refreshCycle;

        if (lineCounter == 9) stream >> energylowTHL;
        if (lineCounter == 10) stream >> energyhighTHL;

        if (lineCounter == 11) stream >> detRad;
        if (lineCounter == 12) stream >> detPosX[0];
        if (lineCounter == 13) stream >> detPosY[0];

        if (lineCounter == 14) { stream >> tempInt; if (tempInt == 1) { detectorAbsorbing = true; ui->checkBoxTransparent->setChecked(false); } else { detectorAbsorbing = false; ui->checkBoxTransparent->setChecked(true); } }
        if (lineCounter == 15) { stream >> tempInt; if (tempInt == 1) { detFileOutput = true; ui->checkBoxFileOutput->setChecked(true); } else { detFileOutput = false; ui->checkBoxFileOutput->setChecked(false); } }

        if (lineCounter == 16) stream >> soilWaterFracVar;
        //if (lineCounter==16) stream >> relHumidityAir ;
        if (lineCounter == 17) stream >> absHumidityAir;
        if (lineCounter == 18) stream >> soilSolidFracVar; soilSolidFrac = soilSolidFracVar;
        if (lineCounter == 19) stream >> atmDensity;
        if (lineCounter == 20) stream >> rigidity;

        if (lineCounter == 21) { stream >> tempInt; numberPrecalcNeutrons = tempInt; }

        if (lineCounter == 22) { stream >> tempInt; if (tempInt == 1) { useRadialBeam = true; useRectShape = false; } else { useRadialBeam = false; useRectShape = true; } }
        if (lineCounter == 23) { stream >> tempInt; if (tempInt == 1) useVolumeSource = true; else useVolumeSource = false; }
        if (lineCounter == 24) { stream >> tempInt; if (tempInt == 1) useHECascadeModel = true; else useHECascadeModel = false; }

        float angleFloat;
        if (lineCounter == 25) { stream >> angleFloat; downwardScotomaAngle = fabs(0.5 * (angleFloat / 360. * 2. * TMath::Pi())); }
        if (lineCounter == 26) { stream >> angleFloat; downwardAcceptanceAngle = fabs(0.5 * (angleFloat / 360. * 2. * TMath::Pi())); }

        if (lineCounter == 27) { stream >> tempInt; if (tempInt == 1) useDetectorSensitiveMaterial = true; else useDetectorSensitiveMaterial = false; }
        if (lineCounter == 28) stream >> detectorSensitiveMaterial;

        if (lineCounter == 29) { stream >> tempInt; if (tempInt == 1) noMultipleScatteringRecording = true; else noMultipleScatteringRecording = false; }
        if (lineCounter == 30) { stream >> tempInt; if (tempInt == 1) trackAllLayers = true; else trackAllLayers = false; }
        if (lineCounter == 31) { stream >> tempInt; if (tempInt == 1) useRealisticModelLayer = true; else useRealisticModelLayer = false; }
        if (lineCounter == 32) { stream >> tempInt; if (tempInt == 1) useRealisticModelDetector = true; else useRealisticModelDetector = false; }

        if (lineCounter == 33) { stream >> tempInt; if (tempInt == 1) { useCylindricalDetector = true; useSphericalDetector = false; useySheetDetector = false; usexSheetDetector = false; } }
        if (lineCounter == 34) { stream >> tempInt; if (tempInt == 1) { useSphericalDetector = true; useCylindricalDetector = false; useySheetDetector = false; usexSheetDetector = false; } }

        if (lineCounter == 35) { stream >> tempInt; if (tempInt == 1) uranosRootOutput = true; else uranosRootOutput = false; }
        if (lineCounter == 36) { stream >> tempInt; if (tempInt == 1) createSeparateFolderEachExport = true; else createSeparateFolderEachExport = false; }

        if (true)
        {
            if (lineCounter == 37) { stream >> tempInt; if (tempInt == 1) ui->checkBoxEpithermalMap->setChecked(true); else ui->checkBoxEpithermalMap->setChecked(false); }
            if (lineCounter == 38) { stream >> tempInt; if (tempInt == 1) { ui->checkBoxEpithermalData->setChecked(true); exportEpithermalData = true; } else { ui->checkBoxEpithermalData->setChecked(false); exportEpithermalData = false; } }
            if (lineCounter == 39) { stream >> tempInt; if (tempInt == 1) ui->checkBoxIntermediateMap->setChecked(true); else ui->checkBoxIntermediateMap->setChecked(false); }
            if (lineCounter == 40) { stream >> tempInt; if (tempInt == 1) { ui->checkBoxIntermediateData->setChecked(true); exportIntermediateData = true; } else { ui->checkBoxIntermediateData->setChecked(false); exportIntermediateData = false; } }
            if (lineCounter == 41) { stream >> tempInt; if (tempInt == 1) ui->checkBoxFastMap->setChecked(true); else ui->checkBoxFastMap->setChecked(false); }
            if (lineCounter == 42) { stream >> tempInt; if (tempInt == 1) { ui->checkBoxFastData->setChecked(true); exportFastData = true; } else { ui->checkBoxFastData->setChecked(false); exportFastData = false; } }
            if (lineCounter == 43) { stream >> tempInt; if (tempInt == 1) ui->checkBoxSelectedMap->setChecked(true); else ui->checkBoxSelectedMap->setChecked(false); }
            if (lineCounter == 44) { stream >> tempInt; if (tempInt == 1) { ui->checkBoxSelectedData->setChecked(true); exportSelectedData = true; } else { ui->checkBoxSelectedData->setChecked(false); exportSelectedData = false; } }
            if (lineCounter == 45) { stream >> tempInt; if (tempInt == 1) ui->checkBoxDetectorOriginMap->setChecked(true); else ui->checkBoxDetectorOriginMap->setChecked(false); }
            if (lineCounter == 46) { stream >> tempInt; if (tempInt == 1) ui->checkBoxDetectorOriginData->setChecked(true); else ui->checkBoxDetectorOriginData->setChecked(false); }
            if (lineCounter == 47) { stream >> tempInt; if (tempInt == 1) ui->checkBoxDetectorDistanceData->setChecked(true); else ui->checkBoxDetectorDistanceData->setChecked(false); }
            if (lineCounter == 48) { stream >> tempInt; if (tempInt == 1) ui->checkBoxDetectorLayerDistanceData->setChecked(true); else ui->checkBoxDetectorLayerDistanceData->setChecked(false); }
            if (lineCounter == 49) { stream >> tempInt; if (tempInt == 1) ui->checkBoxThermalMap->setChecked(true); else ui->checkBoxThermalMap->setChecked(false); }
            if (lineCounter == 50) { stream >> tempInt; if (tempInt == 1) { ui->checkBoxThermalData->setChecked(true); exportThermalData = true; } else { ui->checkBoxThermalData->setChecked(false); exportThermalData = false; } }
            if (lineCounter == 51) { stream >> tempInt; if (tempInt == 1) ui->checkBoxTravelDistGraph->setChecked(true); else ui->checkBoxTravelDistGraph->setChecked(false); }
            if (lineCounter == 52) { stream >> tempInt; if (tempInt == 1) { detLayerFileOutput = true; ui->checkBoxFileOutput3->setChecked(true); } else { detLayerFileOutput = false; ui->checkBoxFileOutput3->setChecked(false); } }
            if (lineCounter == 53) { stream >> tempInt; if (tempInt == 1) ui->checkBox_activateThermal->setChecked(true); else ui->checkBox_activateThermal->setChecked(false); }
            if (lineCounter == 54) { stream >> tempInt; if (tempInt == 1) { ui->checkBox_useImage->setChecked(true); layerMapsImport = true; useImage = true;} else { ui->checkBox_useImage->setChecked(false); layerMapsImport = false; useImage = false;} }
            if (lineCounter == 55) { stream >> tempInt; if (tempInt == 1) showDensityTrackMapSide = true; else showDensityTrackMapSide = false; }
            if (lineCounter == 56) stream >> densityTrackMapSideCutOutValue;
        }

        if (lineCounter == 57) { stream >> tempInt; if (tempInt == 1) clearEveryXNeutrons = true; else clearEveryXNeutrons = false; }
        if (lineCounter == 58) stream >> clearEveryXNeutronsNumber;
        if (lineCounter == 59) { stream >> tempInt; if (tempInt == 1) setAutoRefreshRate = true; else setAutoRefreshRate = false; }
        if (lineCounter == 60) stream >> refreshTime;

        if (true)
        {
            if (lineCounter == 61) stream >> xPosSource;
            if (lineCounter == 62) stream >> yPosSource;
            if (lineCounter == 63) stream >> zPosSource;
            if (lineCounter == 64) stream >> xSizeSource;
            if (lineCounter == 65) stream >> ySizeSource;
            if (lineCounter == 66) stream >> radiusSource;
            if (lineCounter == 67) stream >> sourceEnergy;
        }

        if (lineCounter == 68) { stream >> tempInt; if (tempInt == 1) setAutoRefreshRateClearing = true; else setAutoRefreshRateClearing = false; }
        if (lineCounter == 69) { stream >> tempInt; if (tempInt == 1) exportTrackData = true; else exportTrackData = false; }
        if (lineCounter == 70) { stream >> tempInt; if (tempInt == 1) exportHighResTrackData = true; else exportHighResTrackData = false; }
        if (lineCounter == 71) { stream >> tempInt; if (tempInt == 1) { useySheetDetector = true; useCylindricalDetector = false; useSphericalDetector = false; usexSheetDetector = false; } }
        if (lineCounter == 72) { stream >> tempInt; if (tempInt == 1) { usexSheetDetector = true; useCylindricalDetector = false; useySheetDetector = false; useSphericalDetector = false; } }
        if (lineCounter == 73) stream >> detLength;

        if (lineCounter == 74) { stream >> tempInt; if (tempInt == 1) domainCutoff = true; else domainCutoff = false; }
        if (lineCounter == 75) stream >> domainCutoffFactor;
        if (lineCounter == 76) { stream >> tempInt; if (tempInt == 1) domainCutoff2 = true; else domainCutoff2 = false; }
        if (lineCounter == 77) stream >> domainCutoffMeters;

        if (lineCounter == 78) { stream >> tempInt; if (tempInt == 1) { detTrackFileOutput = true; ui->checkBoxFileOutput2->setChecked(true); } else { detTrackFileOutput = false; ui->checkBoxFileOutput2->setChecked(false); } }
        if (lineCounter == 79) { stream >> tempInt; if (tempInt == 1) { allTrackFileOutput = true; ui->checkBoxExportAllTracks->setChecked(true); } else { allTrackFileOutput = false; ui->checkBoxExportAllTracks->setChecked(false); } }

        if (lineCounter == 80) { stream >> tempInt; if ((((tempInt == 5) || (tempInt == 15)) && (!noThermalRegime)) || (noThermalRegime)) densityMapButtonID = tempInt; }

        if (lineCounter == 81)  {stream >> rPlants ; rPlants = rPlants*0.001; stream >> rCelluloseFrac; if ((rCelluloseFrac < 0)||(rCelluloseFrac > 3)) {rCelluloseFrac = 1;} else {rCelluloseFrac = rCelluloseFrac/1.5; stream >> rCelluloseWaterFrac;}  }

        if (lineCounter == 82) { stream >> tempInt; if (tempInt == 1) { reflectiveBoundaries = true; } else { reflectiveBoundaries = false; } ui->checkBox_ReflectiveBoundaries->setChecked(reflectiveBoundaries); }
        if (lineCounter == 83) { stream >> tempInt; if (tempInt == 1) { periodicBoundaries = true; } else { periodicBoundaries = false; } ui->checkBox_PeriodicBoundaries->setChecked(periodicBoundaries); }
        if (lineCounter == 84) { stream >> rBoronInSoil; stream >> rCInSoil; stream >> rNInSoil; stream >> rNaInSoil; stream >> rKInSoil; stream >> rTiInSoil; stream >> rMnInSoil; stream >> rFeInSoil; stream >> rGdInSoil; }

        lineCounter++;
    }

    fpReadSuccess = setupFootprintFunction();

    fpSoilMoist = soilWaterFracVar;
    soilWaterFrac = soilWaterFracVar;
    atmPressure = atmDensity * 100.;
    fpHum = 10. / 0.6 * (relHumidityAir * 1.) / 50.;
    pressureFac = (atmDensity * 1.) / 1020.;

    string wrkFldrName = (string)workFolder;
    visualization->setWorkFolder(wrkFldrName);
    visualization2->setWorkFolder2(wrkFldrName);

    if (lineCounter > 1) return true;

    return false;
}



/**
 * export GUI settings to the config file.
 *
 */
void MainWindow::exportSettings(string str)
{
    ofstream* stream_out;
    if (str == "")    stream_out = new ofstream("Uranos.cfg", ofstream::out);
    else  stream_out = new ofstream(str, ofstream::out);

    *stream_out << outputFolder << endl;
    if (workFolder == "")  *stream_out << "default" << endl;
    else *stream_out << workFolder << endl;
    *stream_out << inputSpectrumFile << endl;
    if (detectorResponseFunctionFile == "")  *stream_out << "default" << endl;
    else *stream_out << detectorResponseFunctionFile << endl;
    *stream_out << endfFolder << endl;

    *stream_out << std::fixed << std::setprecision(3);

    *stream_out << neutrons << "\t Number of Neutrons" << endl;
    *stream_out << squareDim << "\t Dimension [mm]" << endl;
    *stream_out << beamRadius << "\t Beam Radius [mm]" << endl;

    *stream_out << refreshCycle << "\t #Neutrons Refresh Rate" << endl;
    *stream_out << std::fixed << std::setprecision(9);
    *stream_out << energylowTHL << "\t Lower THL" << endl;
    *stream_out << energyhighTHL << "\t Higher THL" << endl;
    *stream_out << std::fixed << std::setprecision(3);

    *stream_out << detRad << "\t Detector Radius [mm]" << endl;
    *stream_out << detPosX[0] << "\t Detector Pos X [mm]" << endl;
    *stream_out << detPosY[0] << "\t Detector Pos Y [mm]" << endl;

    *stream_out << detectorAbsorbing << "\t Detector Absorbing [bool]" << endl;
    *stream_out << detFileOutput << "\t Detector File Output [bool]" << endl;

    *stream_out << soilWaterFracVar << "\t Soil Volumetric Water Fraction" << endl;
    //*stream_out<<relHumidityAir<<"\t Air Humidity [internal, relative to NTP]"<<endl;
    *stream_out << absHumidityAir << "\t Air Humidity [g/m3]" << endl;
    *stream_out << soilSolidFracVar << "\t Soil Solidity (1-porosity)" << endl;
    *stream_out << atmDensity << "\t Atmospheric Depth [g/cm2]" << endl;
    *stream_out << rigidity << "\t Cutoff Rigidity [GeV]" << endl;

    *stream_out << numberPrecalcNeutrons << "\t #Neutrons Precalulated Spectrum [Power]" << endl;

    *stream_out << useRadialBeam << "\t Radial Source [bool]" << endl;
    *stream_out << useVolumeSource << "\t Volume Source [bool]" << endl;
    *stream_out << useHECascadeModel << "\t HE Cascade Model [bool]" << endl;

    float angleFloat = fabs(2. * (downwardScotomaAngle * 360. / 2. / TMath::Pi())); if (angleFloat > 361) angleFloat = 0;
    *stream_out << angleFloat << "\t Downward Scotoma Angle [0-360]" << endl;
    angleFloat = fabs(2. * (downwardAcceptanceAngle * 360. / 2. / TMath::Pi())); if (angleFloat > 361) angleFloat = 0;
    *stream_out << angleFloat << "\t Downward Acceptance Angle [0-360]" << endl;

    *stream_out << useDetectorSensitiveMaterial << "\t Only Record in Material [bool]" << endl;
    *stream_out << detectorSensitiveMaterial << "\t Record Material No [bool]" << endl;

    *stream_out << noMultipleScatteringRecording << "\t Exclude Multiple Scattering Recording [bool]" << endl;
    *stream_out << trackAllLayers << "\t Track all Layers [bool]" << endl;

    *stream_out << useRealisticModelLayer << "\t Layer Physics Model [bool]" << endl;
    *stream_out << useRealisticModelDetector << "\t Detector Physics Model [bool]" << endl;

    *stream_out << useCylindricalDetector << "\t Cylindrical Detector [bool]" << endl;
    *stream_out << useSphericalDetector << "\t Spherical Detector [bool]" << endl;

    *stream_out << uranosRootOutput << "\t ROOT Output [bool]" << endl;
    *stream_out << createSeparateFolderEachExport << "\t Separate Folder Each Export [bool]" << endl;

    *stream_out << ui->checkBoxEpithermalMap->isChecked() << "\t Export Epithermal Map [bool]" << endl;
    *stream_out << ui->checkBoxEpithermalData->isChecked() << "\t Export Epithermal Data [bool]" << endl;
    *stream_out << ui->checkBoxIntermediateMap->isChecked() << "\t Export Intermediat Map [bool]" << endl;
    *stream_out << ui->checkBoxIntermediateData->isChecked() << "\t Export Intermediat Data [bool]" << endl;
    *stream_out << ui->checkBoxFastMap->isChecked() << "\t Export Fast Map [bool]" << endl;
    *stream_out << ui->checkBoxFastData->isChecked() << "\t Export Fast Data [bool]" << endl;
    *stream_out << ui->checkBoxSelectedMap->isChecked() << "\t Export Detector Energy Map [bool]" << endl;
    *stream_out << ui->checkBoxSelectedData->isChecked() << "\t Export Detector Energy Data [bool]" << endl;

    *stream_out << ui->checkBoxDetectorOriginMap->isChecked() << "\t Export Detector Origins Map [bool]" << endl;
    *stream_out << ui->checkBoxDetectorOriginData->isChecked() << "\t Export Detector Origins Data [bool]" << endl;
    *stream_out << ui->checkBoxDetectorDistanceData->isChecked() << "\t Export Detector Distance Data [bool]" << endl;
    *stream_out << ui->checkBoxDetectorLayerDistanceData->isChecked() << "\t Export Detector Layer Distance Data [bool]" << endl;

    *stream_out << ui->checkBoxThermalMap->isChecked() << "\t Export x Map [bool]" << endl;
    *stream_out << ui->checkBoxThermalData->isChecked() << "\t Export x Data [bool]" << endl;

    *stream_out << ui->checkBoxTravelDistGraph->isChecked() << "\t Export Travel Distance Graph [bool]" << endl;
    *stream_out << ui->checkBoxFileOutput3->isChecked() << "\t Export Detector Layer File Output [bool]" << endl;
    *stream_out << ui->checkBox_activateThermal->isChecked() << "\t Export Detector Physics Ext. Mod. [bool]" << endl;
    *stream_out << ui->checkBox_useImage->isChecked() << "\t Use Layer Maps [bool]" << endl;
    *stream_out << showDensityTrackMapSide << "\t Sideways Tracking [bool]" << endl;
    *stream_out << densityTrackMapSideCutOutValue << "\t Sideways Tracking y +/-Cutout [mm]" << endl;

    *stream_out << ui->checkBoxSaveEvery_2->isChecked() << "\t Clear Every x Neutrons [bool]" << endl;
    *stream_out << clearEveryXNeutronsNumber << "\t Clear Every x Neutrons Number" << endl;
    *stream_out << ui->checkBoxAutoRefreshRate->isChecked() << "\t Refresh Rate Auto Update [bool]" << endl;
    *stream_out << refreshTime << "\t  Refresh Rate Auto Update Time [s]" << endl;

    if (true)
    {
        if (fabs(xPosSource) < 0.0001) *stream_out << "0" << endl; else  *stream_out << xPosSource << endl;
        if (fabs(yPosSource) < 0.0001) *stream_out << "0" << endl; else  *stream_out << yPosSource << endl;
        if (fabs(zPosSource) < 0.0001) *stream_out << "0" << endl; else  *stream_out << zPosSource << endl;
        if (xSizeSource < 0.0001) *stream_out << "0" << endl; else  *stream_out << xSizeSource << endl;
        if (ySizeSource < 0.0001) *stream_out << "0" << endl; else  *stream_out << ySizeSource << endl;
        if (radiusSource < 0.0001) *stream_out << "0" << endl; else  *stream_out << radiusSource << endl;
        *stream_out << sourceEnergy << endl;
    }
    *stream_out << ui->checkBoxClearEveryDisplayRefresh->isChecked() << "\t Set Clearing to Refresh Rate [bool]" << endl;
    *stream_out << ui->checkBoxTrackingData->isChecked() << "\t Export Neutron Track Data [bool]" << endl;
    *stream_out << ui->checkBoxHighResTrackingData->isChecked() << "\t Export High Resolution Neutron Track Data [bool]" << endl;
    *stream_out << useySheetDetector << "\t Detector Sheet along y-Axis [bool]" << endl;
    *stream_out << usexSheetDetector << "\t Detector Sheet along y-Axis [bool]" << endl;
    *stream_out << detLength << "\t Detector Sheet Length [mm]" << endl;
    *stream_out << domainCutoff << "\t use Domain Cutoff Factor [bool]" << endl;
    *stream_out << domainCutoffFactor << "\t Domain Cutoff Factor [float]" << endl;
    *stream_out << domainCutoff2 << "\t use Domain Cutoff Distance [bool]" << endl;
    *stream_out << domainCutoffMeters << "\t Domain Cutoff Distance [m]" << endl;
    *stream_out << detTrackFileOutput << "\t Detector Neutron Track File Output [bool]" << endl;
    *stream_out << allTrackFileOutput << "\t All Neutron Track File Output [bool]" << endl;
    *stream_out << densityMapButtonID << "\t Energy Display Range for Birds-Eye View [int]" << endl;
    *stream_out<<(rPlants*1000.)<<"\t"<<rCelluloseFrac*1.5<<"\t"<<rCelluloseWaterFrac<<"\t Plant Gas Density Multiplicator [kg/m3]\t Plant dry density [g/cm3] \t Plant water density [g/cm3] [float]"<<endl;
    *stream_out << reflectiveBoundaries << "\t Reflective Boundary Conditions [bool]" << endl;
    *stream_out << periodicBoundaries << "\t Periodic Boundary Conditions [bool]" << endl;
    //*stream_out<<rBoronInSoil<<"\t Boron Density in Soil #19 [1e-6 g/cm3]"<<endl;
    *stream_out << rBoronInSoil << "\t" << rCInSoil << "\t" << rNInSoil << "\t" << rNaInSoil << "\t" << rKInSoil << "\t" << rTiInSoil << "\t" << rMnInSoil << "\t" << rFeInSoil << "\t" << rGdInSoil << " Element Density in Soil #19 [1e-6 g/cm3]" << endl;
}


/**
 * destructor of the main window with export GUI settings to be saved.
 *
 */
MainWindow::~MainWindow()
{
    if (!noGUIMode)
    {
        if (configFilePathConfigured)
        {
            exportSettings(configFilePath);
        }
        else
        {
            exportSettings("");
        }
    }
    delete ui;
}

/**
 * export data to be saved in files, structured according to the settings of the GUI.
 *
 */
void MainWindow::exportToSave()
{
    setStatus(1, "Exporting Data");
    if (!noGUIMode) {setStatus(2, "");    delay(5);}

    gROOT->ForceStyle();

    int digitsTotalActualN = log10((double)totalActualNeutrons);
    int digitsTotalN = log10((double)neutrons);
    string zeros = "";
    for (int i = 0; i < (digitsTotalN - digitsTotalN); i++)
    {
        zeros = zeros + "0";
    }

    TDatime* actualTime = new TDatime();
    actualTime->Set();

    int actualYear = actualTime->GetYear();
    int actualMonth = actualTime->GetMonth();
    int actualDay = actualTime->GetDay();
    int actualHour = actualTime->GetHour();
    int actualMinute = actualTime->GetMinute();

    string actualMonthString, actualDayString, actualHourString, actualMinuteString;
    if (actualMonth > 9) actualMonthString = (string)castIntToString(actualMonth); else actualMonthString = "0" + (string)castIntToString(actualMonth);
    if (actualDay > 9) actualDayString = (string)castIntToString(actualDay); else actualDayString = "0" + (string)castIntToString(actualDay);
    if (actualHour > 9) actualHourString = (string)castIntToString(actualHour); else actualHourString = "0" + (string)castIntToString(actualHour);
    if (actualMinute > 9) actualMinuteString = (string)castIntToString(actualMinute); else actualMinuteString = "0" + (string)castIntToString(actualMinute);

    //string datString = (string)castIntToString(actualYear)+"-"+(string)castIntToString(actualMonth)+"-"+(string)castIntToString(actualDay)+"_"+(string)castIntToString(actualHour)+"-"+(string)castIntToString(actualMinute);
    string datString = (string)castIntToString(actualYear) + "" + actualMonthString + "" + actualDayString + "-" + actualHourString + "" + actualMinuteString + "_N" + zeros + castLongToString(totalActualNeutrons);

    string folderString = (string)outputFolder;
    std::replace(folderString.begin(), folderString.end(), '/', '\\');
    string foldercommandstring = "mkdir " + folderString + datString;
    const char* foldercommandstring2 = foldercommandstring.c_str();

    string folderMod = "";

    if (createSeparateFolderEachExport)
    {
        //mkdir(outputFolder+datString);
        system(foldercommandstring2);
        folderMod = datString + "/";
        datString = "";
    }
    else
    {
        if (exportTemporary)
        {
            string datStringTemp = datString + "Temp";
            string foldercommandstringTemp = "mkdir " + folderString + datStringTemp;
            const char* foldercommandstringTemp2 = foldercommandstringTemp.c_str();
            system(foldercommandstringTemp2);
            folderMod = datStringTemp + "/";
            datStringTemp = "";
        }
    }

    //if (doBatchRun) {datString+=castDoubleToString(soilWaterFrac,5)+"_"+castDoubleToString(paramDensity,5)+"_"+castIntToString(paramInt);}
    if ((doBatchRun) || (doDetectorBatchRun) || (doDetectorBatchRun2) || (doDetectorAngleBatchRun) || (doBatchRunDensity) || (doBatchRun2D)) { datString = castIntToString(paramInt); }
    //if (doBatchRunDensity) {datString+=castDoubleToString(soilWaterFrac,5)+"_"+castDoubleToString(paramDensity,5)+"_"+castIntToString(paramInt);}

    int rebinLevel = 1;
    if (ui->spinBoxRebinning->value() > 1) rebinLevel = ui->spinBoxRebinning->value();

    if (true)
    {
        string expSetF = (string)outputFolder + (string)folderMod + "Uranos_" + (string)datString + ".cfg";
        exportSettings(expSetF);
    }

    ofstream* stream_out;

    Int_t nX, nY;

    if (exportSelectedData)
    {
        stream_out = new ofstream(outputFolder + folderMod + "densityMapSelected_" + datString + ".csv", ofstream::out);

        nX = densityMapAlbedo->GetNbinsX() + 1;
        nY = densityMapAlbedo->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(densityMapAlbedo->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();

        if (exportTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapSelected_" + datString + ".csv", ofstream::out);

            nX = densityAlbedoTrackMap->GetNbinsX() + 1;
            nY = densityAlbedoTrackMap->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityAlbedoTrackMap->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (exportHighResTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapSelectedHighRes_" + datString + ".csv", ofstream::out);

            nX = densityAlbedoTrackMapHighRes->GetNbinsX() + 1;
            nY = densityAlbedoTrackMapHighRes->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityAlbedoTrackMapHighRes->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (totalAdditionalDetectorLayers > 0)
        {
            for (int k = 0; k < maxLayersAllowed; k++)
            {
                int vecPos = additionalDetectorLayers[k];
                int tp = k + 1;

                if (vecPos < 0) continue;

                stream_out = new ofstream(outputFolder + folderMod + "densityMapSelected_L"+castIntToString(tp)+"_" + datString + ".csv", ofstream::out);

                nX = addDetLayerVec.at(vecPos)->GetNbinsX() + 1;
                nY = addDetLayerVec.at(vecPos)->GetNbinsY() + 1;

                for (Int_t i = 1; i < nY; i++)
                {
                    for (Int_t j = 1; j < nX; j++)
                    {
                        *stream_out << castDoubleToString(addDetLayerVec.at(vecPos)->GetBinContent(j, i)) << "\t";
                    }
                    *stream_out << endl;
                }

                stream_out->close();
            }

        }
    }

    if (exportEpithermalData)
    {
        stream_out = new ofstream(outputFolder + folderMod + "densityMapEpithermal_" + datString + ".csv", ofstream::out);

        nX = densityMap->GetNbinsX() + 1;
        nY = densityMap->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(densityMap->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();

        if (exportTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapEpithermal_" + datString + ".csv", ofstream::out);

            nX = densityTrackMap->GetNbinsX() + 1;
            nY = densityTrackMap->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityTrackMap->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (exportHighResTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapEpithermalHighRes_" + datString + ".csv", ofstream::out);

            nX = densityTrackMapHighRes->GetNbinsX() + 1;
            nY = densityTrackMapHighRes->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityTrackMapHighRes->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }
    }

    if (exportIntermediateData)
    {
        stream_out = new ofstream(outputFolder + folderMod + "densityMapIntermediateEnergy_" + datString + ".csv", ofstream::out);

        nX = densityMapIntermediate->GetNbinsX() + 1;
        nY = densityMapIntermediate->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(densityMapIntermediate->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();

        if (exportTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapIntermediateEnergy_" + datString + ".csv", ofstream::out);

            nX = densityIntermediateTrackMap->GetNbinsX() + 1;
            nY = densityIntermediateTrackMap->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityIntermediateTrackMap->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (exportHighResTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapIntermediateEnergyHighRes_" + datString + ".csv", ofstream::out);

            nX = densityIntermediateTrackMapHighRes->GetNbinsX() + 1;
            nY = densityIntermediateTrackMapHighRes->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityIntermediateTrackMapHighRes->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }
    }

    if (exportFastData)
    {
        stream_out = new ofstream(outputFolder + folderMod + "densityMapFastNeutron_" + datString + ".csv", ofstream::out);

        nX = densityMapFast->GetNbinsX() + 1;
        nY = densityMapFast->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(densityMapFast->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();

        if (exportTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapFastNeutron_" + datString + ".csv", ofstream::out);

            nX = densityFastTrackMap->GetNbinsX() + 1;
            nY = densityFastTrackMap->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityFastTrackMap->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (exportHighResTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapFastNeutronHighRes_" + datString + ".csv", ofstream::out);

            nX = densityFastTrackMapHighRes->GetNbinsX() + 1;
            nY = densityFastTrackMapHighRes->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityFastTrackMapHighRes->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }
    }

    if (exportThermalData)
    {
        stream_out = new ofstream(outputFolder + folderMod + "densityMapThermalNeutron_" + datString + ".csv", ofstream::out);

        nX = densityMapThermal->GetNbinsX() + 1;
        nY = densityMapThermal->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(densityMapThermal->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();

        if (exportTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapThermalNeutron_" + datString + ".csv", ofstream::out);

            nX = densityThermalTrackMap->GetNbinsX() + 1;
            nY = densityThermalTrackMap->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityThermalTrackMap->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }

        if (exportHighResTrackData)
        {
            stream_out = new ofstream(outputFolder + folderMod + "densityTrackMapThermalNeutronHighRes_" + datString + ".csv", ofstream::out);

            nX = densityThermalTrackMapHighRes->GetNbinsX() + 1;
            nY = densityThermalTrackMapHighRes->GetNbinsY() + 1;

            for (Int_t i = 1; i < nY; i++)
            {
                for (Int_t j = 1; j < nX; j++)
                {
                    *stream_out << castDoubleToString(densityThermalTrackMapHighRes->GetBinContent(j, i)) << "\t";
                }
                *stream_out << endl;
            }
            stream_out->close();
        }
    }

    if (!(ui->checkBoxDetectorOriginData->isChecked()));
    else
    {
        stream_out = new ofstream(outputFolder + folderMod + "detectorOrigins" + datString + ".csv", ofstream::out);

        nX = detectorOriginMap->GetNbinsX() + 1;
        nY = detectorOriginMap->GetNbinsY() + 1;

        for (Int_t i = 1; i < nY; i++)
        {
            for (Int_t j = 1; j < nX; j++)
            {
                *stream_out << castDoubleToString(detectorOriginMap->GetBinContent(j, i)) << "\t";
            }
            *stream_out << endl;
        }
        stream_out->close();
    }

    TH2Fashion(densityTrackMap, "x [m]", "y [m]", setSize);
    TH2Fashion(densityIntermediateTrackMap, "x [m]", "y [m]", setSize);
    TH2Fashion(densityFastTrackMap, "x [m]", "y [m]", setSize);
    TH2Fashion(densityAlbedoTrackMap, "x [m]", "y [m]", setSize);

    TH2Fashion(densityEnergyTrackMap, "x [m]", "y [m]", setSize);

    TH2Fashion(densityThermalTrackMap, "x [m]", "y [m]", setSize);
    TH2Fashion(densityMapThermal, "x [m]", "y [m]", setSize);

    TH2Fashion(densityMapAlbedo, "x [m]", "y [m]", setSize);
    TH2Fashion(densityMap, "x [m]", "y [m]", setSize);
    TH2Fashion(densityMapFast, "x [m]", "y [m]", setSize);
    TH2Fashion(densityMapIntermediate, "x [m]", "y [m]", setSize);

    TH2Fashion(densityTrackMapSide, "x [m]", "z [m]", setSize);
    TH2Fashion(densityTrackMapSideAlbedo, "x [m]", "z [m]", setSize);
    TH2Fashion(densityTrackMapSideDetector, "x [m]", "z [m]", setSize);
    TH2Fashion(densityTrackMapSideThermal, "x [m]", "z [m]", setSize);

    TH2Fashion(detectorOriginMap, "x [m]", "y [m]", setSize);

    /*
    float SetTitleOffsetValue = 1.;
    densityTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityIntermediateTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityFastTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityAlbedoTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);

    densityEnergyTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);

    densityThermalTrackMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityMapThermal->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);

    densityMapAlbedo->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityMapFast->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityMapIntermediate->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);

    densityTrackMapSide->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityTrackMapSideAlbedo->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityTrackMapSideDetector->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    densityTrackMapSideThermal->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);

    detectorOriginMap->GetYaxis()->SetTitleOffset(SetTitleOffsetValue);
    */

    TH1Fashion(detectorDistance, "n", "Distance [mm]", setSize);
    TH1Fashion(detectorDistanceBackScattered, "n", "Albedo Neutron Distance to Detector [mm]", setSize);
    detectorDistanceBackScattered->GetYaxis()->SetTitleOffset(1.3);

    TH1Fashion(detectorLayerDistance, "n", "Distance [mm]", setSize);
    TH1Fashion(detectorLayerDistanceBackscattered, "n", "Albedo Neutron Distance in the Detector Layer [mm]", setSize);
    detectorLayerDistanceBackscattered->GetYaxis()->SetTitleOffset(1.3);

    TH1Fashion(scatteredSurfaceSpectrumHelp, "n", "Energy [MeV]", setSize);

    if (exportSpectrum)
    {
        TH1Fashion(scatteredSurfaceSpectrum, "n", "Energy [MeV]", setSize);
        TH1Fashion(scatteredSurfaceSpectrumBack, "n", "Energy [MeV]", setSize);

        TFile fOut(outputFolder + folderMod + "spectrum_" + datString + ".root", "RECREATE");
        scatteredSurfaceSpectrum->Write();
        scatteredSurfaceSpectrumBack->Write();
        fOut.Close();
    }

    scatteredSurfaceDistance->GetYaxis()->SetTitleOffset(1.3);

    TH1Fashion(scatDepth, "n", "Depth [mm]", setSize);
    TH1Fashion(scatteredSurfaceDepth, "n", "Depth [mm]", setSize);
    TH1Fashion(scatteredSurfaceMaxDepth, "n", "Depth [mm]", setSize);

    set_plot_styleHeatGradient(0.00, 0.25, 0.50, 0.75, 1.00);

    if (ui->radioButton_expGray->isChecked()) set_plot_styleSingleGradient(205, 0.00, 0.25, 0.50, 0.75, 0.99500);;
    if (ui->radioButton_expCool->isChecked()) set_plot_styleCool();
    if (ui->radioButton_expHeatInv->isChecked()) set_plot_styleHeatGradient2(0.00, 0.25, 0.50, 0.75, 0.99500);
    if (ui->radioButton_expDeviation->isChecked()) set_plot_styleHeatGradientModified(0.00, 0.25, 0.50, 0.75, 0.99500);
    if (ui->radioButton_expRainbow->isChecked()) set_plot_styleRainbowGradient(0.00, 0.25, 0.50, 0.75, 0.99500);

    if (true)
    {
        if (!(ui->checkBoxSelectedMap->isChecked()));
        else
        {
            TCanvas* cDensity2 = new TCanvas("cDensity2", "cDensity2", 1600, 1600);
            CanvasFashion(cDensity2);

            cDensity2->SetGrid(0, 0);

            TH2F* densityMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityMapCopy = (TH2F*)densityMapAlbedo->Clone("");
                densityMapCopy->RebinX(rebinLevel); densityMapCopy->RebinY(rebinLevel);

                densityMapCopy->Draw("COLZ"); densityMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityMapAlbedo->Draw("COLZ"); densityMapAlbedo->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity2->SaveAs(outputFolder + folderMod + "densityMapSelectedEnergy_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity2->SaveAs(outputFolder + folderMod + "densityMapSelectedEnergy_" + datString + ".pdf");

            delete densityMapCopy;
            delete cDensity2;

            TCanvas* cDensity2b = new TCanvas("cDensity2b", "cDensity2b", 1600, 1600);
            CanvasFashion(cDensity2b);

            TH2F* densityTrackMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityTrackMapCopy = (TH2F*)densityAlbedoTrackMap->Clone("");
                densityTrackMapCopy->RebinX(rebinLevel); densityTrackMapCopy->RebinY(rebinLevel);

                densityTrackMapCopy->Draw("COLZ"); densityTrackMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityAlbedoTrackMap->Draw("COLZ"); densityAlbedoTrackMap->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity2b->SaveAs(outputFolder + folderMod + "densityTrackMapSelectedEnergy_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity2b->SaveAs(outputFolder + folderMod + "densityTrackMapSelectedEnergy_" + datString + ".pdf");

            delete densityTrackMapCopy;
            delete cDensity2b;
        }

        if (!(ui->checkBoxFastMap->isChecked()));
        else
        {
            TCanvas* cDensity3 = new TCanvas("cDensity3", "cDensity3", 1600, 1600);
            CanvasFashion(cDensity3);

            cDensity3->SetGrid(0, 0);

            TH2F* densityMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityMapCopy = (TH2F*)densityMapFast->Clone("");
                densityMapCopy->RebinX(rebinLevel); densityMapCopy->RebinY(rebinLevel);

                densityMapCopy->Draw("COLZ"); densityMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityMapFast->Draw("COLZ"); densityMapFast->GetZaxis()->SetLabelSize(0.025);
            }
            cDensity3->SaveAs(outputFolder + folderMod + "densityMapFast_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity3->SaveAs(outputFolder + folderMod + "densityMapFast_" + datString + ".pdf");

            delete densityMapCopy;
            delete cDensity3;

            TCanvas* cDensity3b = new TCanvas("cDensity3b", "cDensity3b", 1600, 1600);
            CanvasFashion(cDensity3b);

            TH2F* densityTrackMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityTrackMapCopy = (TH2F*)densityFastTrackMap->Clone("");
                densityTrackMapCopy->RebinX(rebinLevel); densityTrackMapCopy->RebinY(rebinLevel);

                densityTrackMapCopy->Draw("COLZ"); densityTrackMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityFastTrackMap->Draw("COLZ"); densityFastTrackMap->GetZaxis()->SetLabelSize(0.025);
            }
            cDensity3b->SaveAs(outputFolder + folderMod + "densityTrackMapFast_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity3b->SaveAs(outputFolder + folderMod + "densityTrackMapFast_" + datString + ".pdf");

            delete densityTrackMapCopy;
            delete cDensity3b;
        }
        if (!(ui->checkBoxIntermediateMap->isChecked()));
        else
        {
            TCanvas* cDensity4 = new TCanvas("cDensity4", "cDensity4", 1600, 1600);
            CanvasFashion(cDensity4);

            cDensity4->SetGrid(0, 0);

            TH2F* densityMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityMapCopy = (TH2F*)densityMapIntermediate->Clone("");
                densityMapCopy->RebinX(rebinLevel); densityMapCopy->RebinY(rebinLevel);

                densityMapCopy->Draw("COLZ"); densityMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityMapIntermediate->Draw("COLZ"); densityMapIntermediate->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity4->SaveAs(outputFolder + folderMod + "densityMapIntermediate_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity4->SaveAs(outputFolder + folderMod + "densityMapIntermediate_" + datString + ".pdf");

            delete densityMapCopy;
            delete cDensity4;


            TCanvas* cDensity4b = new TCanvas("cDensity4b", "cDensity4b", 1600, 1600);
            CanvasFashion(cDensity4b);

            cDensity4b->SetGrid(0, 0);

            TH2F* densityTrackMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityTrackMapCopy = (TH2F*)densityIntermediateTrackMap->Clone("");
                densityTrackMapCopy->RebinX(rebinLevel); densityTrackMapCopy->RebinY(rebinLevel);

                densityTrackMapCopy->Draw("COLZ"); densityTrackMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityIntermediateTrackMap->Draw("COLZ"); densityIntermediateTrackMap->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity4b->SaveAs(outputFolder + folderMod + "densityTrackMapIntermediate_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity4b->SaveAs(outputFolder + folderMod + "densityTrackMapIntermediate_" + datString + ".pdf");

            delete densityTrackMapCopy;
            delete cDensity4b;
        }

        if (!(ui->checkBoxEpithermalMap->isChecked()));
        else
        {
            TCanvas* cDensity5 = new TCanvas("cDensity5", "cDensity5", 1600, 1600);
            CanvasFashion(cDensity5);

            cDensity5->SetGrid(0, 0);

            TH2F* densityMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityMapCopy = (TH2F*)densityMap->Clone("");
                densityMapCopy->RebinX(rebinLevel); densityMapCopy->RebinY(rebinLevel);

                densityMapCopy->Draw("COLZ"); densityMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityMap->Draw("COLZ"); densityMap->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity5->SaveAs(outputFolder + folderMod + "densityMapEpithermal_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked())   cDensity5->SaveAs(outputFolder + folderMod + "densityMapEpithermal_" + datString + ".pdf");

            delete densityMapCopy;
            delete cDensity5;

            TCanvas* cDensity5b = new TCanvas("cDensity5b", "cDensity5b", 1600, 1600);
            CanvasFashion(cDensity5b);

            cDensity5b->SetGrid(0, 0);

            TH2F* densityTrackMapCopy = new TH2F();

            if (rebinLevel > 1)
            {
                densityTrackMapCopy = (TH2F*)densityTrackMap->Clone("");
                densityTrackMapCopy->RebinX(rebinLevel); densityTrackMapCopy->RebinY(rebinLevel);

                densityTrackMapCopy->Draw("COLZ"); densityTrackMapCopy->GetZaxis()->SetLabelSize(0.025);
            }
            else
            {
                densityTrackMap->Draw("COLZ"); densityTrackMap->GetZaxis()->SetLabelSize(0.025);
            }

            cDensity5b->SaveAs(outputFolder + folderMod + "densityTrackMapEpithermal_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked())   cDensity5b->SaveAs(outputFolder + folderMod + "densityTrackMapEpithermal_" + datString + ".pdf");

            delete densityTrackMapCopy;
            delete cDensity5b;
        }
        if (!(ui->checkBoxThermalMap->isChecked()));
        else
        {
            if (!noThermalRegime)
            {
                TCanvas* cDensity6 = new TCanvas("cDensity6", "cDensity6", 1600, 1600);
                CanvasFashion(cDensity6);

                cDensity6->SetGrid(0, 0);

                TH2F* densityMapCopy = new TH2F();

                if (rebinLevel > 1)
                {
                    densityMapCopy = (TH2F*)densityMapThermal->Clone("");
                    densityMapCopy->RebinX(rebinLevel); densityMapCopy->RebinY(rebinLevel);

                    densityMapCopy->Draw("COLZ"); densityMapCopy->GetZaxis()->SetLabelSize(0.025);
                }
                else
                {
                    densityMapThermal->Draw("COLZ"); densityMapThermal->GetZaxis()->SetLabelSize(0.025);
                }

                cDensity6->SaveAs(outputFolder + folderMod + "densityMapThermal_" + datString + ".png");

                if (ui->checkBoxFileOutputPDF->isChecked()) cDensity6->SaveAs(outputFolder + folderMod + "densityMapThermal_" + datString + ".pdf");

                delete densityMapCopy;
                delete cDensity6;

                TCanvas* cDensity6b = new TCanvas("cDensity6b", "cDensity6b", 1600, 1600);
                CanvasFashion(cDensity6b);

                cDensity6b->SetGrid(0, 0);

                TH2F* densityTrackMapCopy = new TH2F();

                if (rebinLevel > 1)
                {
                    densityTrackMapCopy = (TH2F*)densityThermalTrackMap->Clone("");
                    densityTrackMapCopy->RebinX(rebinLevel); densityTrackMapCopy->RebinY(rebinLevel);

                    densityTrackMapCopy->Draw("COLZ"); densityTrackMapCopy->GetZaxis()->SetLabelSize(0.025);
                }
                else
                {
                    densityThermalTrackMap->Draw("COLZ"); densityThermalTrackMap->GetZaxis()->SetLabelSize(0.025);
                }

                cDensity6b->SaveAs(outputFolder + folderMod + "densityTrackMapThermal_" + datString + ".png");

                if (ui->checkBoxFileOutputPDF->isChecked()) cDensity6b->SaveAs(outputFolder + folderMod + "densityTrackMapThermal_" + datString + ".pdf");

                delete densityTrackMapCopy;
                delete cDensity6b;
            }
        }

        if (!(ui->checkBoxDetectorOriginMap->isChecked()));
        else
        {
            TCanvas* cDensity2k = new TCanvas("cDensity2k", "cDensity2k", 1600, 1600);
            CanvasFashion(cDensity2k);

            cDensity2k->SetGrid(0, 0);

            detectorOriginMap->Draw("COLZ"); detectorOriginMap->GetZaxis()->SetLabelSize(0.025);

            cDensity2k->SaveAs(outputFolder + folderMod + "detectorOriginsMap_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked()) cDensity2k->SaveAs(outputFolder + folderMod + "detectorOriginsMap_" + datString + ".pdf");
        }


        if (!(ui->checkBoxTravelDistGraph->isChecked()));
        else
        {
            TCanvas* cDetector = new TCanvas("cDetector", "cDetector", 3240, 1080);
            CanvasFashion(cDetector);

            cDetector->SetGrid(0, 0);

            cDetector->Divide(3, 1);

            if (ui->checkBoxLogTravelDistance->isChecked())
            {
                TPad* p1 = (TPad*)cDetector->cd(1); p1->SetLogy(); TPad* p2 = (TPad*)cDetector->cd(2); p2->SetLogy(); TPad* p3 = (TPad*)cDetector->cd(3); p3->SetLogy();
            }

            cDetector->cd(1);
            scatteredSurfaceDistance->Draw("");

            //detectorDistance->GetYaxis()->SetRangeUser(0,15);
            cDetector->cd(2);
            detectorLayerDistanceBackscattered->Draw("");

            cDetector->cd(3);
            detectorDistanceBackScattered->Draw("");

            cDetector->SaveAs(outputFolder + folderMod + "NeutronDistances_" + datString + ".png");

            if (ui->checkBoxFileOutputPDF->isChecked())
            {
                //formatForVectorGraphics();
                //gStyle->SetLineWidth(0.1);
                //gStyle->SetHistLineWidth(0.1);
                cDetector->SaveAs(outputFolder + folderMod + "NeutronDistances_" + datString + ".pdf");
            }

            delete cDetector;
        }

        if (!(ui->checkBoxDetectorLayerDistanceData->isChecked()));
        else
        {
            //stream_out = new ofstream(outputFolder+folderMod+"AlbedoNeutronLayerDistances_"+datString+".dat",ofstream::out);

            stream_out = new ofstream(outputFolder + folderMod + "AlbedoNeutronLayerDistances_" + datString + ".csv", ofstream::out);

            *stream_out << "Distance [mm]" << "\t" << "Events" << endl;

            nX = detectorLayerDistanceBackscattered->GetNbinsX();

            for (Int_t i = 1; i < nX; i++)
            {
                *stream_out << castDoubleToString(detectorLayerDistanceBackscattered->GetBinCenter(i)) << "\t" << castDoubleToString(detectorLayerDistanceBackscattered->GetBinContent(i)) << endl;
            }
            stream_out->close();
        }

        if (!(ui->checkBoxDetectorDistanceData->isChecked()));
        else
        {
            stream_out = new ofstream(outputFolder + folderMod + "AlbedoNeutronDetectorDistances_" + datString + ".csv", ofstream::out);

            *stream_out << "Distance [mm]" << "\t" << "Events" << endl;

            nX = detectorDistanceBackScattered->GetNbinsX();

            for (Int_t i = 1; i < nX; i++)
            {
                *stream_out << castDoubleToString(detectorDistanceBackScattered->GetBinCenter(i)) << "\t" << castDoubleToString(detectorDistanceBackScattered->GetBinContent(i)) << endl;
            }
            stream_out->close();
        }
    }

    if (uranosRootOutput)
    {
        TFile f(outputFolder + folderMod + "liveViewGraphs_" + datString + ".root", "RECREATE");

        for (Int_t k = 0; k < liveTHs.size(); k++)
        {
            liveTHs.at(k)->Write();
        }
        f.Close();
    }

    if (uranosRootOutput)
    {
        TFile fileUROutp(outputFolder + folderMod + "uranosRawHistos_" + datString + ".root", "RECREATE");

        for (Int_t k = 0; k < allTHs.size(); k++)
        {
            allTHs.at(k)->Write();
        }
        fileUROutp.Close();
    }
    delay(100);

    if (!noGUIMode) setStatus(1, "");
}

/**
 * button click event function for export to save.
 *
 */
void MainWindow::on_pushButton_clicked()
{
    if (alreadyStarted)
    {
        exportToSave();
    }

    delay(10);
}

double rangeIntegral = -1, rangeUpToIntegral, integralResult;

/**
 * Replot the Footprint graphs.
 *
 */
void MainWindow::replotFootprint()
{
    if (!fpReadSuccess) return;

    float lowEdge = 1.5;
    float highEdge = 1. * squareDim / 1000.; if (highEdge > 600) highEdge = 600.;

    ::plotGUIyBinsFootprFuncLine[0] = 50.;
    ::plotGUIyBinsFootprFuncLine[1] = 50.;
    ::plotGUIyBinsFootprFuncLine[2] = 50.;
    ::plotGUIyBinsFootprFuncLine[3] = 0.;

    TF1* rangeFunction = new TF1("rangeFunction", rangeFunctionCombined, 0, 600, 2);
    rangeFunction->SetNpx(400);
    //rangeFunction->SetLineWidth(2);
    rangeFunction->SetParameters(fpHum, fpSoilMoist);

    string integralResultText;

    if (integralSliderMoved)
    {
        if (rangeIntegral < 0) rangeIntegral = rangeFunction->Integral(1, 500);
        rangeUpToIntegral = rangeFunction->Integral(1, plotGUIxBinsFootprFuncLine[1]);

        integralResult = 100. * rangeUpToIntegral / rangeIntegral; if (integralResult > 100) integralResult = 100;

        integralResultText = castDoubleToString(integralResult, 5) + " %";
        ui->labelIntegral->setText(QString::fromStdString(integralResultText));
    }

    integralSliderMoved = false;

    if (ui->customPlotFP->plottableCount() < 3) { ui->customPlotFP->addGraph(); ui->customPlotFP->addGraph(); }
    ui->customPlotFP->graph(0)->data()->clear();
    ui->customPlotFP->graph(0)->setName("Footprint");
    ui->customPlotFP->graph(0)->setPen(QPen(someBlue));
    ui->customPlotFP->graph(0)->setBrush(QBrush(QColor(215, 237, 255, 90)));


    ui->customPlotFP->graph(1)->data()->clear();
    ui->customPlotFP->graph(1)->setPen(QPen(QColor(120, 120, 120)));
    ui->customPlotFP->graph(1)->setBrush(QBrush(QColor(120, 120, 120, 30)));
    ui->customPlotFP->graph(1)->setData(::plotGUIxBinsFootprFuncLine, ::plotGUIyBinsFootprFuncLine);

    for (int i = 0; i < 200; ++i)
    {
        ::plotGUIxBinsFootprFunc[i] = (highEdge - lowEdge) * (i + 1) / (200.);
        ::plotGUIyBinsFootprFunc[i] = rangeFunction->Eval(::plotGUIxBinsFootprFunc[i]);
    }

    ui->customPlotFP->graph(0)->setData(::plotGUIxBinsFootprFunc, ::plotGUIyBinsFootprFunc);
    if (ui->checkBoxFPLog->isChecked())
    {
        QSharedPointer<QCPAxisTickerLog> logTickerFP(new QCPAxisTickerLog);
        ui->customPlotFP->yAxis->setTicker(logTickerFP);
        ui->customPlotFP->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else
    {
        ui->customPlotFP->yAxis->setScaleType(QCPAxis::stLinear);

        QSharedPointer<QCPAxisTickerFixed> linTickerFP(new QCPAxisTickerFixed);
        linTickerFP->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
        linTickerFP->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

        ui->customPlotFP->yAxis->setTicker(linTickerFP);
    }
    ui->customPlotFP->update();
    ui->customPlotFP->rescaleAxes();
    ui->customPlotFP->yAxis->setRange(0.1, 23);

    ui->customPlotFP->replot();
}

/**
 * Function to be called On clicking View Spectrum button in Showcase tab.
 * This is the old Sato 2008 function, which is not used any more in the simulation itself (but only for graphical representation)
 */
void MainWindow::on_pushButton_2_clicked()
{
    TH1F* precalculatedSpectrumC = new TH1F("precalculatedSpectrumC", "precalculated Spectrum", 10000, 1e-9, 10000);
    logaxis(precalculatedSpectrumC);

    //TString inputSpectrumFile = "G:/Analyse/Simulation/Cosmics/results19/allHistos0.99_0.00.root";

    generateNormalizedSpectrum(precalculatedSpectrumC, TMath::Power(10, ::numberPrecalcNeutrons), inputSpectrumFile);

    TString phiBmean = "0.229*TMath::Power(x/2.31,0.721)*exp(-x/2.31)+0.0516*exp(-TMath::Power(log10(x)-log10(126),2)/(2*TMath::Power(log10(2.17),2)) ) + 0.00108*log10(x/3.33*TMath::Power(10.,12.)) * (1+tanh(1.62*log10(x/9.59*TMath::Power(10.,8.)))) *(1-tanh(1.48*log10(x/299.)))";

    TString a1 = "(12.9 + 15.7 / (1+exp( ([1]-5.62) /1.79) )  )";			// parameter: cut-off rigidity r_c
    TString a2 = "(0.00706 + 0.00057 / (1+exp( ([1]-5.99) /1.94) )  )";		// parameter: cut-off rigidity r_c
    TString a3 = "(0.975 - 0.210 / (1+exp( ([1]-0.99) /2.24) )  )";			// parameter: cut-off rigidity r_c
    TString a4 = "(0.0084 + 0.00441 / (1+exp( ([1]-2.24) /2.66) )  )";		// parameter: cut-off rigidity r_c

    TString a5 = "(-0.00701 + 0.0258 / (1+exp( ([1]-10.9) /2.38) )  )";		// parameter: cut-off rigidity r_c
    TString a9 = "(642 - 189 / (1+exp( ([1]-2.32) /0.897) )  )";			// parameter: cut-off rigidity r_c
    TString a10 = "(0.00112 + 0.000181 / (1+exp( ([1]-8.84) /0.587) )  )";	// parameter: cut-off rigidity r_c
    TString a11 = "(1.26 - 0.958 / (1+exp( ([1]-3.18) /1.47) )  )";			// parameter: cut-off rigidity r_c

    TString c4 = "(" + a5 + "+0.000171*[0]/(1+0.53*exp(0.00136*[0])) )";               // parameter: atmospheric density
    TString c12 = "(" + a9 + "*(exp(-" + a10 + "*[0]))+" + a11 + "*exp(-0.0133*[0]) )";        // parameter: atmospheric density

    TString phiL = "(" + a1 + "*(exp(-" + a2 + "*[0]) - " + a3 + "*exp(-" + a4 + "*[0])) )";		// parameter: atmospheric density

    TString g3 = "(-25.2+2.73/([2]+0.0715))";									//parameter: weight fraction of water w
    TString g5 = "(0.348+3.35*[2]-1.57*TMath::Power([2],2) )";					//parameter: weight fraction of water w

    TString phiB = "0.229*TMath::Power(x/2.31,0.721)*exp(-x/2.31)+" + c4 + "*exp(-TMath::Power(log10(x)-log10(126),2)/(2*TMath::Power(log10(2.17),2)) ) + 0.00108*log10(x/3.33*TMath::Power(10.,12.)) * (1+tanh(1.62*log10(x/9.59*TMath::Power(10.,8.)))) *(1-tanh(1.48*log10(x/" + c12 + ")))";

    TString fG = "TMath::Power(10.,-0.0235-0.0129*(log10(x)-" + g3 + ")*(1-tanh(0.969*log10(x/" + g5 + ") ) )  )";

    TF1* phiAllFunc = new TF1("phiAllFunc", phiL + "*(" + phiB + ")*(" + fG + ")", 0.000000001, 10000);

    TString phiFunctionString;
    if (useBasicSpectrum) phiFunctionString = phiL + "*(" + phiBmean + ")*(" + fG + ")";
    else phiFunctionString = phiL + "*(" + phiB + ")*(" + fG + ")";

    if (useBasicSpectrum) phiAllFunc = new TF1("phiAllFunc", phiL + "*(" + phiBmean + ")*(" + fG + ")", 0.000000001, 10000);

    phiAllFunc->SetParameter(0, atmDensity);		//Atm Density
    phiAllFunc->SetParameter(1, rigidity);		// cut-off rigidity
    phiAllFunc->SetParameter(2, 0.99);		// water fraction

    float funcMax = phiAllFunc->GetMaximum(0.01, 900);

    if (ui->WidgetSpectrum->plottableCount() < 2) ui->WidgetSpectrum->addGraph();
    ui->WidgetSpectrum->graph(0)->data()->clear();
    ui->WidgetSpectrum->graph(0)->setPen(QPen(someBlue));
    ui->WidgetSpectrum->graph(0)->setName("Randomized Incoming Spectrum");

    for (int i = 100; i < 1001; ++i)
    {
        ::plotGUIxBinsIncOnlySpectrum[i - 100] = precalculatedSpectrumC->GetBinLowEdge(i + 1);
        ::plotGUIyBinsIncOnlySpectrum[i - 100] = precalculatedSpectrumC->GetBinContent(i + 1);
    }
    QSharedPointer<QCPAxisTickerLog> logTickerWS(new QCPAxisTickerLog);
    ui->WidgetSpectrum->xAxis->setTicker(logTickerWS);
    ui->WidgetSpectrum->xAxis->setScaleType(QCPAxis::stLogarithmic);
    ui->WidgetSpectrum->graph(0)->setData(::plotGUIxBinsIncOnlySpectrum, ::plotGUIyBinsIncOnlySpectrum);
    ui->WidgetSpectrum->update();
    ui->WidgetSpectrum->rescaleAxes();

    float specMax = ui->WidgetSpectrum->yAxis->range().upper;

    if (ui->WidgetSpectrum->plottableCount() < 2) ui->WidgetSpectrum->addGraph();
    ui->WidgetSpectrum->graph(1)->data()->clear();
    ui->WidgetSpectrum->graph(1)->setName("Full Analytical Spectrum");
    ui->WidgetSpectrum->graph(1)->setPen(QPen(Qt::black));

    for (int i = 100; i < 1001; ++i)
    {
        ::plotGUIxBinsIncFullSpectrum[i - 100] = precalculatedSpectrumC->GetBinLowEdge(i + 1);
        ::plotGUIyBinsIncFullSpectrum[i - 100] = 0.95 * specMax / funcMax * phiAllFunc->Eval(::plotGUIxBinsIncFullSpectrum[i - 100]);
    }

    ui->WidgetSpectrum->graph(1)->setData(::plotGUIxBinsIncFullSpectrum, ::plotGUIyBinsIncFullSpectrum);
    ui->WidgetSpectrum->update();
    //ui->WidgetSpectrum->rescaleAxes();

    ui->WidgetSpectrum->replot();

    delete phiAllFunc;
    delete precalculatedSpectrumC;
}

/**
 * Gets the splined detector response function energy model from a matrix.
 *
 * @param matrix - matrix name
 * @return splineEnergy - spline function with logarithmic energy and detection probability normalized to 1.
 */
TSpline3* MainWindow::getSplinedEnergyModelFromMatrix(TMatrixF* matrix, bool logXValues, bool normalize)
{
    const int nData = matrix->GetNrows();

    double xMin = 0, xMax = 0;

    //double xVal[nData] = {0};
    //double yVal[nData] = {0};

    double* xVal = new double[nData+6];
    double* yVal = new double[nData+6];

    for (int i = 0; i < nData; i++)
    {
        xVal[i+3] = (*matrix)(i,0);
        yVal[i+3] = (*matrix)(i,1);

        if (i == 0) xMin = xVal[i+3];
        xMax = xVal[i+3];

        if (logXValues)
        {
            xVal[i+3] = TMath::Log10(xVal[i+3]);        //probably easier to spline?
            if (xVal[i+3] < - 15) xVal[i+3] = -15;
            if (xVal[i+3] > 1000) xVal[i+3] = 1000;
        }
        /*
        if (normalize)
        {
            if ((ymax > 1) || (ymax < 0.01)) yVal[i] = yVal[i] / ymax;  //scale to 1 maximum value //only normalize if values are larger than 1
        }
        */
    }

    if (logXValues)
    {
        if (xMin > -4.8)
        {
            xVal[0] = -5.5;
            xVal[1] = xVal[3]-0.1;
            xVal[2] = xVal[3]-0.05;
            yVal[0] = 0;
            yVal[1] = 0;
            yVal[2] = 0;
        }
        else
        {
            xVal[0] = -5.3;
            xVal[1] = -5.2;
            xVal[2] = -5.1;
            yVal[0] = yVal[3];
            yVal[1] = yVal[3];
            yVal[2] = yVal[3];
        }
        if (xMax < 9)
        {
            xVal[nData+3] = xVal[nData+2] + 0.1;
            xVal[nData+4] = xVal[nData+2] + 0.2;
            xVal[nData+5] = xVal[nData+2] + 0.3;
            yVal[nData+3] = 0;
            yVal[nData+4] = 0;
            yVal[nData+5] = 0;
        }
        else
        {
            xVal[nData+3] = 9.1;
            xVal[nData+4] = 9.2;
            xVal[nData+5] = 9.3;
            yVal[nData+3] = yVal[nData+2];
            yVal[nData+4] = yVal[nData+2];
            yVal[nData+5] = yVal[nData+2];
        }
    }

    vector<double> xVal2, yVal2;
    double xDiff = 0, yDiff = 0, yRel = 0;

    xVal2.push_back(xVal[0]);
    yVal2.push_back(yVal[0]);

    for (int i = 0; i < nData+5; i++)
    {
        xDiff = xVal[i+1] - xVal[i];
        yRel = yVal[i+1] / yVal[i];
        if (xDiff > 0.3)
        {
            for (float j = 0.3; j < xDiff; j += 0.3)
            {
                xVal2.push_back(xVal[i] + j);
                yVal2.push_back(yVal[i] + TMath::Power(10,j - xDiff) * (yVal[i+1] - yVal[i]));
                //cout<<TMath::Log10(j - xDiff)<<" "<<xDiff<<" "<<xVal2.back()<<" "<<yVal2.back()<<endl; ;
            }
        }
        else
        {
            //if ((yRel > 1.25) || (yRel < 0.75)) // doesn't work for now, there's a mistake somewhere
            if (false)
            {
                if ((yRel > 10) || (yRel < 0.1))
                //if (false)
                {

                    if ((yRel > 1) && (abs(yVal[i]) > 0.001))
                    {
                        for (float j = 0.25; j < 1; j += 0.25)
                        {
                            xVal2.push_back(TMath::Log10(j)*xDiff + xVal[i+1]);
                            yVal2.push_back(j * yVal[i+1]);
                        }
                    }
                    else
                    {
                        for (float j = 0.25; j < 1; j += 0.25)
                        {
                            xVal2.push_back(TMath::Log10(1-j)*xDiff + xVal[i+1]);
                            yVal2.push_back((1-j) * yVal[i]);
                        }
                    }
                }
                else
                {
                    if (yRel > 1)
                    {
                        for (float j = 1.25; j < yRel; j += 0.25)
                        {
                            xVal2.push_back(TMath::Log10((j-1) / (yRel-1))*xDiff + xVal[i+1]);
                            yVal2.push_back(j * yVal[i]);
                        }
                    }
                    else
                    {
                        yRel = 1. / yRel;

                        for (float j = 1.25; j < yRel; j += 0.25)
                        {
                            xVal2.push_back(TMath::Log10((j-1) / (yRel-1))*xDiff + xVal[i+1]);
                            yVal2.push_back(yVal[i] / j);
                        }
                    }
                }
                //cout<<" "<<xVal2.back()<<" "<<yVal2.back()<<endl;
            }
            else
            {
                xVal2.push_back(xVal[i+1]);
                yVal2.push_back(yVal[i+1]);
            }
        }
    }

    const int nDataMod = xVal2.size();

    double* xValMod = new double[nDataMod];
    double* yValMod = new double[nDataMod];

    for (int i = 0; i < nDataMod; i++)
    {
        xValMod[i] = xVal2.at(i);
        yValMod[i] = yVal2.at(i);
    }

    TGraph* graphEnergyFromFile = new TGraph(nData, xValMod, yValMod);

    TString name = "Spline" + castFloatToString(yVal[nData/2],6);

    TSpline3* splineEnergy = new TSpline3(name, graphEnergyFromFile, "");

    delete graphEnergyFromFile;

    return splineEnergy;
}

/**
 * Gets the splined detector response function energy model from a matrix.
 *
 * @param matrix - matrix name
 * @return splineEnergy - spline function with logarithmic energy and detection probability normalized to 1.
 */
TSpline5* MainWindow::getSplinedEnergyModelFromMatrix2(TMatrixF* matrix, bool logXValues, bool normalize)
{
    const int nData = matrix->GetNrows();

    double xMin = 0, xMax = 0;

    //double xVal[nData] = {0};
    //double yVal[nData] = {0};

    double* xVal = new double[nData+6];
    double* yVal = new double[nData+6];

    for (int i = 0; i < nData; i++)
    {
        xVal[i+3] = (*matrix)(i,0);
        yVal[i+3] = (*matrix)(i,1);

        if (i == 0) xMin = xVal[i+3];
        xMax = xVal[i+3];

        if (logXValues)
        {
            xVal[i+3] = TMath::Log10(xVal[i+3]);        //probably easier to spline?
            if (xVal[i+3] < - 15) xVal[i+3] = -15;
            if (xVal[i+3] > 1000) xVal[i+3] = 1000;
        }
        /*
        if (normalize)
        {
            if ((ymax > 1) || (ymax < 0.01)) yVal[i] = yVal[i] / ymax;  //scale to 1 maximum value //only normalize if values are larger than 1
        }
        */
    }

    if (logXValues)
    {
        if (xMin > -4.8)
        {
            xVal[0] = -5.5;
            xVal[1] = xVal[3]-0.1;
            xVal[2] = xVal[3]-0.05;
            yVal[0] = 0;
            yVal[1] = 0;
            yVal[2] = 0;
        }
        else
        {
            xVal[0] = -5.3;
            xVal[1] = -5.2;
            xVal[2] = -5.1;
            yVal[0] = yVal[3];
            yVal[1] = yVal[3];
            yVal[2] = yVal[3];
        }
        if (xMax < 9)
        {
            xVal[nData+3] = xVal[nData+2] + 0.1;
            xVal[nData+4] = xVal[nData+2] + 0.2;
            xVal[nData+5] = xVal[nData+2] + 0.3;
            yVal[nData+3] = 0;
            yVal[nData+4] = 0;
            yVal[nData+5] = 0;
        }
        else
        {
            xVal[nData+3] = 9.1;
            xVal[nData+4] = 9.2;
            xVal[nData+5] = 9.3;
            yVal[nData+3] = yVal[nData+2];
            yVal[nData+4] = yVal[nData+2];
            yVal[nData+5] = yVal[nData+2];
        }
    }

    TGraph* graphEnergyFromFile = new TGraph(nData, xVal, yVal);

    TString name = "Spline" + castFloatToString(yVal[nData/2],6);

    TSpline5* splineEnergy = new TSpline5(name, graphEnergyFromFile, "");

    delete graphEnergyFromFile;

    return splineEnergy;
}

/**
 * Gets the splined detector response function energy model from a file.
 *
 * @param dname - file name
 * @return splineEnergy - spline function with logarithmic energy and detection probability normalized to 1.
 */
TSpline3* MainWindow::getSplinedDetectorEnergyModelFromFile(TString dname, int linesToSkip, bool logXValues, bool normalize)
{
    const int nData = getNumberOfFileEntries(dname) - linesToSkip;

    TString curLine;

    //double xVal[nData] = {0};
    //double yVal[nData] = {0};

    double* xVal = new double[nData];
    double* yVal = new double[nData];

    int lineCounter = 0;
    double temp;
    double ymax = 0;

    ifstream input_stream(dname, ios::in);

    while (curLine.ReadLine(input_stream))
    {
        istrstream stream(curLine.Data());

        if (lineCounter >= linesToSkip)
        {
            stream >> temp;
            xVal[lineCounter - linesToSkip] = temp;
            if (temp < 0) { if (!noGUIMode) { setStatus(2, "Det. Energy File: Neg. x-Value"); delay(1500); } temp = 0; }
            stream >> temp;
            yVal[lineCounter - linesToSkip] = temp;
            if (temp > ymax) ymax = temp;
            if (temp < 0) { if (!noGUIMode) { setStatus(2, "Det. Energy File: Neg. y-Value"); delay(1500); } temp = 0; }
        }

        lineCounter++;
    }
    input_stream.close();

    if (logXValues || normalize)
    {
        for (int i = 0; i < nData; i++)
        {
            if (logXValues)
            {
                xVal[i] = TMath::Log10(xVal[i]);        //probably easier to spline?
                if (xVal[i] < - 15) xVal[i] = -15;
                if (xVal[i] > 1000) xVal[i] = 1000;
            }
            if (normalize)
            {
                if ((ymax > 1) || (ymax < 0.01)) yVal[i] = yVal[i] / ymax;  //scale to 1 maximum value //only normalize if values are larger than 1
            }
        }
    }

    TGraph* graphEnergyFromFile = new TGraph(nData, xVal, yVal);

    TSpline3* splineEnergy = new TSpline3("splineEnergy" + dname, graphEnergyFromFile);

    delete graphEnergyFromFile;

    return splineEnergy;
}

/**
 * Function to be called on change of value of Sampling magnitude textbox in Showcase tab
 * @param arg1
 */
void MainWindow::on_spinBox_2_valueChanged(int arg1)
{
    ::numberPrecalcNeutrons = ui->spinBox_2->value();
}

/**
 * to set new data flag
 *
 */
void declareNewData()
{
    if (newData) newData = false;
    else     newData = true;
}


/**
 * Clear up the liveView data histograms.
 *
 * @param vectorTH
 */
void histoLiveClearUp(vector< TH1* >* vectorTH)
{
    for (int k = 0; k < vectorTH->size(); k++)
    {
        for (int l = 0; l < vectorTH->at(k)->GetNbinsX(); l++)
        {
            vectorTH->at(k)->SetBinContent(l, 0);
            for (int m = 0; m < vectorTH->at(k)->GetNbinsY(); m++)
            {
                vectorTH->at(k)->SetBinContent(l, m, 0);
            }
            vectorTH->at(k)->Reset("");
            //vectorTH->at(k)->ResetStats();
        }
    }
}


/**
 * deletes the statistics boxes in the TH1 because they are drawn automatically
 * @param allTHs
 *
 */
void deleteStatsTH(TH1* allTHs)
{
    TPaveStats* stA = (TPaveStats*)allTHs->FindObject("stats");
    if (stA != 0x0) { delete stA; stA = 0; }
}


float oldTemperature = 0;
double rLuftDensity = 0;
double rLuftWaterDensity = 0;
double rLuftWater = 0, rLuft = 0;

/**
 *  Formula: ps = exp(-6094.4642/T + 21.1249952 - 0.027245552*T + 1.6853396*10-5*T^2 + 2.4575506*ln(T))
 *  with T in K, ps in Pa and temperature in K
 *
 * @param temperature temperature of the water body.
 * @return steam density in g/cm^3.
 */
double getRLuftWasser(float temperature)
{
    if (fabs(temperature - oldTemperature) > 1)
    {
        double psVap = exp(-6094.4642 / temperature + 21.1249952 - 0.027245552 * temperature + 1.6853396 * 1e-5 * temperature * temperature + 2.4575506 * log(temperature)); //Saetttigungsdampfdruck

        rLuftWaterDensity = psVap / (461.5 * temperature) * 1e-6;

        oldTemperature = temperature;
    }

    return rLuftWaterDensity;
}

//deprecated
TMatrixF* csMatrixHP;
TMatrixF* csMatrixOP;
TMatrixF* csMatrixNP;
TMatrixF* csMatrixHabsP;
TMatrixF* csMatrixOabsP;
TMatrixF* csMatrixNabsP;

double csMatrixHPpreviousEnergy, csMatrixHPpreviousCS;
double csMatrixOPpreviousEnergy, csMatrixOPpreviousCS;
double csMatrixNPpreviousEnergy, csMatrixNPpreviousCS;
double csMatrixOabsPpreviousEnergy, csMatrixOabsPpreviousCS;
double csMatrixHabsPpreviousEnergy, csMatrixHabsPpreviousCS;
double csMatrixNabsPpreviousEnergy, csMatrixNabsPpreviousCS;

/**
 * draws a function, mainly used for legendre polynomials
 * @param function.
 * @param outputFolder.
 * @param filename
 *
 */
void printTF1(TF1* function, TString outputFolder, TString filename)
{
    TCanvas* cfunc = new TCanvas("cfunc", "cfunc Spectrum", 1280, 1280);
    CanvasFashion(cfunc);

    //cfunc->SetLogy();
    //function->GetYaxis()->SetRangeUser(0.5e-2,0.8);
    function->Draw("");

    cfunc->SaveAs(outputFolder + "/" + filename + ".png");
}

// deprecated
bool checkDetectorHit(double x, double y, double z, double deltaZ, double phi, double theta)
{
    float	xt = cos(phi) * TMath::Abs(tan(theta) * deltaZ) + x;
    float	yt = sin(phi) * TMath::Abs(tan(theta) * deltaZ) + y;

    return true;
}


/**
 * generates a random number in MeV from an AmBe spectrum
 * pointer to an already seeded Random generator is needed
 * @param seeded random generator r
 * @return abszRnd
 */
double getAmBeEnergy(TRandom* r)
{
    bool gotIt;
    double abszRnd, ordRnd;

    gotIt = false;
    while (!gotIt)
    {
        abszRnd = r->Rndm() * 11.;
        ordRnd = r->Rndm() * 3.7;
        for (int x = 0; x < 53; x++)
        {
            if (abszRnd < spectrumAmBeBins[x])
            {
                if (ordRnd < spectrumAmBeValues[x])  gotIt = true;
                break;
            }
        }
    }
    return abszRnd;
}


/**
 * reads a textfile with two columns and returns a matrix with two columns
 * @param folder
 * @param filename
 * @returns - sigmaRedMatrix
 */
TMatrixF readSigmaEnergy(TString folder, TString filename)
{
    TString line;
    float temp;
    int linecounter = 1, lengthFile = 0;
    //string dname = (string)folder + "/" + (string)filename;
    string dname = (string)folder + (string)filename;

    float minDeviation = 1.005;
    bool useMatrixReduction = true;

    if (outputCSfilenames) cout << "Reading Energy Cross Sections from " << dname;

    ifstream input_streamCheck(dname, ios::in);

    while (line.ReadLine(input_streamCheck)) { lengthFile++; }

    input_streamCheck.close();

    if (lengthFile < 8) { cout << "File not available " << dname << endl; TMatrixF sigmaMatrixFail(1, 2); return sigmaMatrixFail; }

    const int matrixSize = lengthFile - 8;
    TMatrixF sigmaMatrix(matrixSize, 2);

    ifstream input_stream(dname, ios::in);
    while (line.ReadLine(input_stream))
    {
        istrstream stream(line.Data());

        if (linecounter > 8)
        {
            stream >> temp;
            sigmaMatrix(linecounter - 9, 0) = temp;
            stream >> temp;
            sigmaMatrix(linecounter - 9, 1) = temp;
        }
        linecounter++;
    }
    input_stream.close();

    if (!useMatrixReduction) return sigmaMatrix;

    int lineCounter2 = 1;
    float lastValue, actualValue;
    TMatrixF sigmaCopyMatrix(matrixSize, 2);

    lastValue = sigmaMatrix(0, 0);
    sigmaCopyMatrix(0, 0) = sigmaMatrix(0, 0);
    sigmaCopyMatrix(0, 1) = sigmaMatrix(0, 1);

    for (int l = 1; l < matrixSize; l++)
    {
        actualValue = sigmaMatrix(l, 0);
        if ((actualValue > minDeviation * lastValue) || (l == matrixSize - 1) ||  (l == matrixSize - 2) ||  (l == matrixSize - 3))
        {
            for (int m = 0; m < 2; m++)
            {
                sigmaCopyMatrix(lineCounter2, m) = sigmaMatrix(l, m);
            }
            lineCounter2++;
            lastValue = sigmaMatrix(l, 0);
        }
        else {}
    }
    lineCounter2 = lineCounter2 - 1;

    int killedElements = matrixSize - lineCounter2;
    const int newSize = lineCounter2;

    if (outputCSfilenames) cout << " (killed " << killedElements << " of " << matrixSize << ", " << lineCounter2 << " left)" << endl;

    TMatrixF sigmaRedMatrix(newSize, 2);

    for (int l = 0; l < lineCounter2; l++)
    {
        for (int m = 0; m < 2; m++)
        {
            sigmaRedMatrix(l, m) = sigmaCopyMatrix(l, m);
        }
    }

    return sigmaRedMatrix;
}

/**
 * changes the cross sections of a matrix by adding a number, which is given in standard deviations, both for the low energy and high energy regime
 * @param sigmaMatrix
 * @param factorLowE
 * @param factorHighE
 */
void modifyCSmatrix(TMatrixF* sigmaMatrix, float factorLowE, float factorHighE)
{
    float energy;
    double sigma;

    for (int l = 0; l < sigmaMatrix->GetNrows(); l++)
    {
        energy = (*sigmaMatrix)(l, 0);
        sigma = (*sigmaMatrix)(l, 1);
        if (energy < 1E6) (*sigmaMatrix)(l, 1) = sigma + sigma * factorLowE;
        if (energy >= 1E6) (*sigmaMatrix)(l, 1) = sigma + sigma * factorHighE;
    }
}


/**
 * reads ENDF tabulated angular coefficients
 * exports vector<TMatrixF> consisting of:
 * TMatrixF for the angles
 * TMatrixF for the cumulated probability distribution
 * int coefficients = maximum coefficients +1
 * @param folder
 * @param filename
 * @param coefficients
 * @returns - vector
 */
vector<TMatrixF> readAngularTabulatedCoefficients(TString folder, TString filename, const int coefficients)
{
    TString line;
    string lineStr, part, temp;
    string tempFl;
    int lineCounter = -1, lengthFile = 0, matrixElm, lineNumber = 0;
    int minusThere = 0;
    //int templength, columns = 0;
    bool setSkipLines;
    //string dname = (string)folder + "/" + (string)filename;
    string dname = (string)folder + (string)filename;

    vector<TMatrixF> dataVec;

    int matrCounter = -1;

    if (outputCSfilenames) cout << "Reading Angular Tabulated Distributions from " << dname;

    ifstream input_streamCheck(dname, ios::in);

    while (line.ReadLine(input_streamCheck)) { lineStr = line;  if ((lineStr.substr(0, 3) == "0.0") && (lineStr.substr(31, 1) == "0"))lengthFile++; }
    input_streamCheck.close();

    if (outputCSfilenames) cout << " (" << lengthFile << " items)" << endl;

    if (lengthFile < 1) { cout << "File not available " << dname << endl; TMatrixF sigmaMatrixFail(1, 2); dataVec.push_back(sigmaMatrixFail); dataVec.push_back(sigmaMatrixFail); return dataVec; }

    const int matrixSize = lengthFile;
    TMatrixF angleMatrix(matrixSize, coefficients);
    TMatrixF cumulatedProbMatrix(matrixSize, coefficients);

    for (int l = 0; l < matrixSize; l++)
    {
        for (int m = 0; m < coefficients; m++)
        {
            angleMatrix(l, m) = 0;
            cumulatedProbMatrix(l, m) = 0;
        }
    }

    matrixElm = 0;

    // a file with angular coefficients consists of 6+x columns of 11 characters for a number. Columns are not separated by a delimiter.

    ifstream input_stream(dname, ios::in);

    while (line.ReadLine(input_stream))
    {
        setSkipLines = false;

        lineStr = line;

        if ((lineStr.substr(0, 3) == "0.0") && (lineStr.substr(31, 1) == "0"))
        {
            lineCounter = 0;
            matrCounter++;
            //matrCounter = 0;
            angleMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(10, 11));
            cumulatedProbMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(10, 11));
            lineCounter++;
            setSkipLines = true;
        }

        if (setSkipLines)
        {
            line.ReadLine(input_stream);
        }
        else
        {
            if (lineStr.substr(0, 1) == "-") { minusThere = 1; }
            else { minusThere = 0; }
            angleMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(0, 10 + minusThere));
            if (lineCounter > 1) { cumulatedProbMatrix(matrCounter, lineCounter) = cumulatedProbMatrix(matrCounter, lineCounter - 1) + endfNumberConv(lineStr.substr(10 + minusThere, 11)); }
            else cumulatedProbMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(10 + minusThere, 11));
            lineCounter++;

            if (lineCounter > coefficients - 1);
            else
            {
                angleMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(21 + minusThere, 11));
                if (lineCounter > 1) { cumulatedProbMatrix(matrCounter, lineCounter) = cumulatedProbMatrix(matrCounter, lineCounter - 1) + endfNumberConv(lineStr.substr(32 + minusThere, 11)); }
                else cumulatedProbMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(32 + minusThere, 11));
                lineCounter++;

                angleMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(43 + minusThere, 11));
                if (lineCounter > 1) { cumulatedProbMatrix(matrCounter, lineCounter) = cumulatedProbMatrix(matrCounter, lineCounter - 1) + endfNumberConv(lineStr.substr(54 + minusThere, 11)); }
                else cumulatedProbMatrix(matrCounter, lineCounter) = endfNumberConv(lineStr.substr(54 + minusThere, 11));
                lineCounter++;
            }
        }

        lineNumber++;
    }

    int maxIndex = 0;
    for (int l = 0; l < matrixSize; l++)
    {
        maxIndex = 0;
        for (int m = 0; m < coefficients; m++)
        {
            if ((maxIndex == 0) && (angleMatrix(l, m) == 0) && (cumulatedProbMatrix(l, m) == 0)) { maxIndex = m - 1; }

            if (maxIndex != 0)
            {
                angleMatrix(l, m) = angleMatrix(l, maxIndex);
                cumulatedProbMatrix(l, m) = cumulatedProbMatrix(l, maxIndex);
            }
        }
    }

    double maxValue;
    for (int l = 0; l < matrixSize; l++)
    {
        maxValue = cumulatedProbMatrix(l, coefficients - 1);
        for (int m = 1; m < coefficients; m++)
        {
            cumulatedProbMatrix(l, m) = cumulatedProbMatrix(l, m) / maxValue;
            if (cumulatedProbMatrix(l, m) == 1) angleMatrix(l, m) = 1;
        }
    }

    if (false)
    {
        for (int l = 1; l < 4; l++)
        {
            for (int m = 0; m < coefficients; m++)
            {
                cout << angleMatrix(l, m) << " ";
                cout << cumulatedProbMatrix(l, m) << endl;
            }
        }
    }

    dataVec.push_back(angleMatrix); dataVec.push_back(cumulatedProbMatrix);
    return dataVec;
}


/**
 * reads ENDF legendre polynomial tabulated angular coefficients
 * exports as a matrix
 * @param folder
 * @param filename
 * @returns - sigmaRedMatrix
 */
TMatrixF readAngularCoefficients(TString folder, TString filename)
{
    TString line;
    string lineStr, part;
    string tempFl;
    int linecounter = -1, lengthFile = 0, matrixElm;
    int minusThere = 0;
    int columns = 0;
    //string dname = (string)folder + "/" + (string)filename;
    string dname = (string)folder + (string)filename;

    if (outputCSfilenames) cout << "Reading Angular Distributions from " << dname;

    ifstream input_streamCheck(dname, ios::in);

    while (line.ReadLine(input_streamCheck)) { lengthFile++; }
    input_streamCheck.close();

    if (lengthFile < 1) { cout << "File not available " << dname << endl; TMatrixF sigmaMatrixFail(1, 2); return sigmaMatrixFail; }

    const int matrixSize = lengthFile;
    TMatrixF sigmaMatrix(matrixSize, 21);

    for (int l = 0; l < matrixSize; l++)
    {
        for (int m = 0; m < 21; m++)
        {
            sigmaMatrix(l, m) = 0;
        }
    }

    matrixElm = 0;

    // a file with angular coefficients consists of 6+1 columns of 11 characters for a number. Columns are not separated by a delimiter.

    ifstream input_stream(dname, ios::in);

    while (line.ReadLine(input_stream))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if (!(stream.good())) break;

        if ((linecounter > 0) && (endfNumberConv(lineStr.substr(0, 11)) < sigmaMatrix(linecounter - 1, 0))) { columns = 1; }
        else
        {
            stream >> tempFl;  if (!(stream.good())) break;
            linecounter++;
            sigmaMatrix(linecounter, 0) = endfNumberConv(tempFl);
            matrixElm = 0;
            columns = 0;
        }

        for (int l = 0; l < 6; l++)
        {
            part = "0";
            matrixElm++; if (matrixElm > 20) break;
            if (columns == 0) part = lineStr.substr(11 + l * 11, 11);
            if (columns == 1)
            {
                if (l == 0) { if (lineStr.substr(0, 1) == "-") minusThere = 1; else minusThere = 0; }
                if (minusThere == 0)
                {
                    if (l == 0) part = lineStr.substr(0, 10);
                    else part = lineStr.substr(10 + (l - 1) * 11, 11);
                }
                else
                {
                    part = lineStr.substr(0 + l * 11, 11);
                }
            }
            //part = lineStr.substr(11-(11*columns)+l*11,11);

            sigmaMatrix(linecounter, matrixElm) = endfNumberConv(part);
            //if ( endfNumberConv(part)==0) break;
        }
    }

    input_stream.close();

    const int matrixSizeRed = linecounter + 1;
    TMatrixF sigmaRedMatrix(matrixSizeRed, 21);

    for (int l = 0; l < matrixSizeRed; l++)
    {
        for (int m = 0; m < 21; m++)
        {
            sigmaRedMatrix(l, m) = sigmaMatrix(l, m);

        }
    }
    if (outputCSfilenames) cout << " (" << matrixSizeRed << " elements)" << endl;

    return sigmaRedMatrix;
}


//containers for the Sato 2016 parammeters and file names of the containers
float gji[10][8];
string paramGjiFile = "param_gji.dat";

float acorr1[10][28];
string acorr1File = "param_acorr1.dat";

float acorr2[10][28];
string acorr2File = "param_acorr2.dat";

float acorr3[10][28];
string acorr3File = "param_acorr3.dat";

float acorr6[10][28];
string acorr6File = "param_acorr6.dat";

float acorr7[10][28];
string acorr7File = "param_acorr7.dat";

float piAng[10][1008];
string parampiAngularFile = "param_piAngular.dat";

/**
 * load the Sato 2016 parameter data from the file folder path
 * @param fileFolder
 * @returns - error
 */
bool MainWindow::loadParamData(string fileFolder)
{
    TString line;
    string dname;
    string lineStr, part;
    string tempFl;
    //remove unused vars
    int linecounter = -1, lengthFile = 0;
    int i = 0, j = 0;
    float temp;
    bool error = false;

    cout << "Reading Sato 2016 Parameters from " << fileFolder << endl;

    i = 0; j = 0;
    dname = fileFolder + paramGjiFile;

    ifstream input_stream(dname, ios::in);
    while (line.ReadLine(input_stream))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            gji[l][j] = temp;
        }
        j++;
    }
    input_stream.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + paramGjiFile);
    }

    i = 0; j = 0;
    dname = fileFolder + acorr1File;

    ifstream input_stream1(dname, ios::in);
    while (line.ReadLine(input_stream1))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            acorr1[l][j] = temp;
        }
        j++;
    }
    input_stream1.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + acorr1File);
    }

    i = 0; j = 0;
    dname = fileFolder + acorr2File;

    ifstream input_stream2(dname, ios::in);
    while (line.ReadLine(input_stream2))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            acorr2[l][j] = temp;
        }
        j++;
    }
    input_stream2.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + acorr2File);
    }

    i = 0; j = 0;
    dname = fileFolder + acorr3File;

    ifstream input_stream3(dname, ios::in);
    while (line.ReadLine(input_stream3))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            acorr3[l][j] = temp;
        }
        j++;
    }
    input_stream3.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + acorr3File);
    }

    i = 0; j = 0;
    dname = fileFolder + acorr6File;

    ifstream input_stream6(dname, ios::in);
    while (line.ReadLine(input_stream6))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            acorr6[l][j] = temp;
        }
        j++;
    }
    input_stream6.close();

    if (j < 1)
    {
        error = true;
        if (!noGUIMode) setStatus(2, "Failed to load " + acorr6File);
    }

    i = 0; j = 0;
    dname = fileFolder + acorr7File;

    ifstream input_stream7(dname, ios::in);
    while (line.ReadLine(input_stream7))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            acorr7[l][j] = temp;
        }
        j++;
    }
    input_stream7.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + acorr7File);
    }

    i = 0; j = 0;
    dname = fileFolder + parampiAngularFile;

    ifstream input_stream8(dname, ios::in);
    while (line.ReadLine(input_stream8))
    {
        istrstream stream(line.Data());

        lineStr = line;
        if ((!(stream.good()))) break;

        for (int l = 0; l < 10; l++)
        {
            stream >> temp;
            piAng[l][j] = temp;
        }
        j++;
    }
    input_stream8.close();

    if (j < 1)
    {
        error = true;
        setStatus(2, "Failed to load " + parampiAngularFile);
    }

    return error;
}

//the specific values the a parameters in Sato 2016 were sampled for
const float atmIndex[28] = { 0, 0.15,0.46,1.24,2.9,6.7,13.2,24.1,40.6,58.8,79.9,109.05,148.5,202.5,299.5,424.5,558,666,718,755.5,794.5,834.5,876.5,920,965,1000,1024,10000 };
float retValue[11]; //container


/**
 * returns the index for the set of parameters for altitude corrected fluxes
 * @param atm
 * @param column
 * @param matrixNumber
 * @returns - the index
 */
float* getPar(float atm, int column, int matrixNumber)
{
    if (atm < 0) return 0;
    if ((column < 0) || (column > 9)) return 0;
    if ((matrixNumber < 0) || (matrixNumber > 5)) return 0;

    int line;
    for (int i = 27; i >= 0; i--)
    {
        if (atmIndex[i] <= atm)
        {
            line = i;
            break;
        }
    }

    if (matrixNumber == 1) { retValue[1] = acorr1[column][line]; retValue[2] = acorr1[column + 1][line];  retValue[3] = acorr1[column + 2][line];  retValue[4] = acorr1[column + 3][line]; retValue[5] = acorr1[column + 4][line];  retValue[6] = acorr1[column][line + 1]; retValue[7] = acorr1[column + 1][line + 1];  retValue[8] = acorr1[column + 2][line + 1];  retValue[9] = acorr1[column + 3][line + 1]; retValue[10] = acorr1[column + 4][line + 1]; }
    if (matrixNumber == 2) { retValue[1] = acorr2[column][line]; retValue[2] = acorr2[column + 1][line];  retValue[3] = acorr2[column + 2][line];  retValue[4] = acorr2[column + 3][line]; retValue[5] = acorr2[column + 4][line];  retValue[6] = acorr2[column][line + 1]; retValue[7] = acorr2[column + 1][line + 1];  retValue[8] = acorr2[column + 2][line + 1];  retValue[9] = acorr2[column + 3][line + 1]; retValue[10] = acorr2[column + 4][line + 1]; }
    if (matrixNumber == 3) { retValue[1] = acorr3[column][line]; retValue[2] = acorr3[column + 1][line];  retValue[3] = acorr3[column + 2][line];  retValue[4] = acorr3[column + 3][line]; retValue[5] = acorr3[column + 4][line];  retValue[6] = acorr3[column][line + 1]; retValue[7] = acorr3[column + 1][line + 1];  retValue[8] = acorr3[column + 2][line + 1];  retValue[9] = acorr3[column + 3][line + 1]; retValue[10] = acorr3[column + 4][line + 1]; }
    if (matrixNumber == 4) { retValue[1] = acorr6[column][line]; retValue[2] = acorr6[column + 1][line];  retValue[3] = acorr6[column + 2][line];  retValue[4] = acorr6[column + 3][line]; retValue[5] = acorr6[column + 4][line];  retValue[6] = acorr6[column][line + 1]; retValue[7] = acorr6[column + 1][line + 1];  retValue[8] = acorr6[column + 2][line + 1];  retValue[9] = acorr6[column + 3][line + 1]; retValue[10] = acorr6[column + 4][line + 1]; }
    if (matrixNumber == 5) { retValue[1] = acorr7[column][line]; retValue[2] = acorr7[column + 1][line];  retValue[3] = acorr7[column + 2][line];  retValue[4] = acorr7[column + 3][line]; retValue[5] = acorr7[column + 4][line];  retValue[6] = acorr7[column][line + 1]; retValue[7] = acorr7[column + 1][line + 1];  retValue[8] = acorr7[column + 2][line + 1];  retValue[9] = acorr7[column + 3][line + 1]; retValue[10] = acorr7[column + 4][line + 1]; }

    retValue[0] = line;

    return retValue;
}

//the specific values the b parameters in Sato 2016 were sampled for
const float atmAngIndex[18] = { 0.765, 1.901, 4.091, 9.503, 17.103, 31.263, 49.983, 67.763, 92.163, 125.563, 171.163, 233.563, 364.963, 483.563, 631.563, 774.683, 897.803, 1036.361 };
const float rigAngIndex[7] = { 1, 3, 5, 7, 10, 15, 20 };

//containers
double retValueAngulardLowrLow[10];
double retValueAngulardLowrHigh[10];
double retValueAngulardHighrLow[10];
double retValueAngulardHighrHigh[10];

//definitions
const float alLow = -0.05;
const float alHigh = 0.05;
const float phiMin = 0.001;

//runtime variables
double ratioAtm = 0;
double ratioRigidity = 0;
float parMin = 0.01;
int rigID = -1;
int atmID = -1;

double atmOld = -1;
double rigidityOld = -1;
bool newParSet = false;


/**
 * parameters atmDepth(g/cm^2), cutoff rigidity(GV), energy (MeV), cosTheta
 * @param atm
 * @param rigidityHere
 * @param parameterNo
 * @returns - Angular Par
 */
double* getAngularPar(float atm, float rigidityHere, int parameterNo)
{
    int parameter = parameterNo - 1;
    int line;

    if ((rigidityHere == rigidityOld) && (atmOld == atm))
    {
        newParSet = false;
    }
    else
    {
        newParSet = true;
        rigidityOld = rigidityHere;
        atmOld = atm;
    }

    //make sure that parameters are only calculated once
    if (newParSet)
    {
        if (atm < 0.765) atmID = 1;
        if (atm > 1036.36) atmID = 17;
        if (rigidityHere < 1) rigID = 1;
        if (rigidityHere > 20) rigID = 6;

        if (atmID < 0)
        {
            for (int i = 17; i >= 1; i--)
            {
                if (atmAngIndex[i] <= atm)
                {
                    atmID = i + 1;
                    break;
                }
            }
        }

        if (rigID < 0)
        {
            for (int i = 6; i >= 1; i--)
            {
                if (rigAngIndex[i] <= rigidityHere)
                {
                    rigID = i + 1;
                    break;
                }
            }
        }
        ratioAtm = (atm - atmAngIndex[atmID - 1]) / (atmAngIndex[atmID] - atmAngIndex[atmID - 1]);
        ratioRigidity = (rigidityHere - rigAngIndex[rigID - 1]) / (rigAngIndex[rigID] - rigAngIndex[rigID - 1]);
    }

    for (int j = 0; j < 10; j++)
    {
        line = parameter * 126 + (atmID - 1) * 7 + (rigID - 1);
        retValueAngulardLowrLow[j] = piAng[j][line];
        line = parameter * 126 + (atmID - 1) * 7 + (rigID);
        retValueAngulardLowrHigh[j] = piAng[j][line];
        line = parameter * 126 + (atmID) * 7 + (rigID - 1);
        retValueAngulardHighrLow[j] = piAng[j][line];
        line = parameter * 126 + (atmID) * 7 + (rigID);
        retValueAngulardHighrHigh[j] = piAng[j][line];

        if ((j == 3) || (j == 6) || (j == 9))
        {
            if (retValueAngulardLowrLow[j] < parMin) retValueAngulardLowrLow[j] = parMin;
            if (retValueAngulardLowrHigh[j] < parMin) retValueAngulardLowrHigh[j] = parMin;
            if (retValueAngulardHighrLow[j] < parMin) retValueAngulardHighrLow[j] = parMin;
            if (retValueAngulardHighrHigh[j] < parMin) retValueAngulardHighrHigh[j] = parMin;
        }
    }

    return retValueAngulardLowrLow;
}



/**
 * parameters atmDepth(g/cm^2), cutoff rigidity(GV), energy (MeV), cosTheta
 * calculates the angular dependance function of the CR spectrum
 * @param atm
 * @param rigidityHere
 * @param energy
 * @param cosTheta
 * @returns - angular dependance
 */
float calcParaAdep(float atm, float rigidityHere, double energy, double cosTheta)
{
    double p1Ang[32];
    double sumValues[16];
    double sumTotal[4];
    double tmp;
    double rightofCosTheta = (cosTheta - alLow) / (alHigh - alLow);
    double bAng[4];
    double cAng[2];

    //missing here: adding of fluxes without blackhole

    //the following functions call getAngularPar, which fills the retValueAngular(dLow/High/rLow/High) containers
    getAngularPar(atm, rigidityHere, 1);
    p1Ang[0] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((+retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((+retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[1] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((+retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((+retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[2] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((+retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((+retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[3] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((+retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((+retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (p1Ang[0] < 0) p1Ang[0] = 0;
    if (p1Ang[1] < 0) p1Ang[1] = 0;
    if (p1Ang[2] < 0) p1Ang[2] = 0;
    if (p1Ang[3] < 0) p1Ang[3] = 0;

    getAngularPar(atm, rigidityHere, 2);
    p1Ang[4] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((+retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((+retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[5] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((+retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((+retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[6] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((+retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((+retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[7] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((+retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((+retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (phiMin < p1Ang[0]) tmp = phiMin; else tmp = p1Ang[0]; if (tmp - p1Ang[0] > p1Ang[4]) p1Ang[4] = tmp - p1Ang[0];
    if (phiMin < p1Ang[1]) tmp = phiMin; else tmp = p1Ang[1]; if (tmp - p1Ang[1] > p1Ang[5]) p1Ang[5] = tmp - p1Ang[1];
    if (phiMin < p1Ang[2]) tmp = phiMin; else tmp = p1Ang[2]; if (tmp - p1Ang[2] > p1Ang[6]) p1Ang[6] = tmp - p1Ang[2];
    if (phiMin < p1Ang[3]) tmp = phiMin; else tmp = p1Ang[3]; if (tmp - p1Ang[3] > p1Ang[7]) p1Ang[7] = tmp - p1Ang[3];

    getAngularPar(atm, rigidityHere, 3);
    p1Ang[8] = TMath::Power(10, retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9])));
    p1Ang[9] = TMath::Power(10, retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9])));
    p1Ang[10] = TMath::Power(10, retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9])));
    p1Ang[11] = TMath::Power(10, retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9])));

    getAngularPar(atm, rigidityHere, 4);
    p1Ang[12] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((+retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[13] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((+retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[14] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((+retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[15] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((+retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (p1Ang[12] < 0) p1Ang[12] = 0;
    if (p1Ang[13] < 0) p1Ang[13] = 0;
    if (p1Ang[14] < 0) p1Ang[14] = 0;
    if (p1Ang[15] < 0) p1Ang[15] = 0;

    getAngularPar(atm, rigidityHere, 5);
    p1Ang[16] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[17] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[18] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[19] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (phiMin < p1Ang[12]) tmp = phiMin; else tmp = p1Ang[12]; if (tmp - p1Ang[12] > p1Ang[16]) p1Ang[16] = tmp - p1Ang[12];
    if (phiMin < p1Ang[13]) tmp = phiMin; else tmp = p1Ang[13]; if (tmp - p1Ang[13] > p1Ang[17]) p1Ang[17] = tmp - p1Ang[13];
    if (phiMin < p1Ang[14]) tmp = phiMin; else tmp = p1Ang[14]; if (tmp - p1Ang[14] > p1Ang[18]) p1Ang[18] = tmp - p1Ang[14];
    if (phiMin < p1Ang[15]) tmp = phiMin; else tmp = p1Ang[15]; if (tmp - p1Ang[15] > p1Ang[19]) p1Ang[19] = tmp - p1Ang[15];

    getAngularPar(atm, rigidityHere, 6);
    p1Ang[20] = TMath::Power(10, retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9])));
    p1Ang[21] = TMath::Power(10, retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9])));
    p1Ang[22] = TMath::Power(10, retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9])));
    p1Ang[23] = TMath::Power(10, retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9])));

    getAngularPar(atm, rigidityHere, 7);
    p1Ang[24] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[25] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[26] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[27] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (p1Ang[24] < 1) p1Ang[24] = 1;
    if (p1Ang[25] < 1) p1Ang[25] = 1;
    if (p1Ang[26] < 1) p1Ang[26] = 1;
    if (p1Ang[27] < 1) p1Ang[27] = 1;

    getAngularPar(atm, rigidityHere, 8);
    p1Ang[28] = retValueAngulardLowrLow[0] + retValueAngulardLowrLow[1] / (1. + exp((retValueAngulardLowrLow[2] - log10(energy)) / retValueAngulardLowrLow[3])) + retValueAngulardLowrLow[4] / (1. + exp((retValueAngulardLowrLow[5] - log10(energy)) / retValueAngulardLowrLow[6])) + retValueAngulardLowrLow[7] / (1. + exp((retValueAngulardLowrLow[8] - log10(energy)) / retValueAngulardLowrLow[9]));
    p1Ang[29] = retValueAngulardLowrHigh[0] + retValueAngulardLowrHigh[1] / (1. + exp((retValueAngulardLowrHigh[2] - log10(energy)) / retValueAngulardLowrHigh[3])) + retValueAngulardLowrHigh[4] / (1. + exp((retValueAngulardLowrHigh[5] - log10(energy)) / retValueAngulardLowrHigh[6])) + retValueAngulardLowrHigh[7] / (1. + exp((retValueAngulardLowrHigh[8] - log10(energy)) / retValueAngulardLowrHigh[9]));
    p1Ang[30] = retValueAngulardHighrLow[0] + retValueAngulardHighrLow[1] / (1. + exp((retValueAngulardHighrLow[2] - log10(energy)) / retValueAngulardHighrLow[3])) + retValueAngulardHighrLow[4] / (1. + exp((retValueAngulardHighrLow[5] - log10(energy)) / retValueAngulardHighrLow[6])) + retValueAngulardHighrLow[7] / (1. + exp((retValueAngulardHighrLow[8] - log10(energy)) / retValueAngulardHighrLow[9]));
    p1Ang[31] = retValueAngulardHighrHigh[0] + retValueAngulardHighrHigh[1] / (1. + exp((retValueAngulardHighrHigh[2] - log10(energy)) / retValueAngulardHighrHigh[3])) + retValueAngulardHighrHigh[4] / (1. + exp((retValueAngulardHighrHigh[5] - log10(energy)) / retValueAngulardHighrHigh[6])) + retValueAngulardHighrHigh[7] / (1. + exp((retValueAngulardHighrHigh[8] - log10(energy)) / retValueAngulardHighrHigh[9]));
    if (p1Ang[28] < 0) p1Ang[28] = 0;
    if (p1Ang[29] < 0) p1Ang[29] = 0;
    if (p1Ang[30] < 0) p1Ang[30] = 0;
    if (p1Ang[31] < 0) p1Ang[31] = 0;

    for (int i = 0; i < 4; i++)
    {
        sumValues[i] = p1Ang[i] * (1. - fabs(alLow)) + p1Ang[i + 4] * (1. - TMath::Power(fabs(alLow), (1 + p1Ang[i + 8])) / (1 + p1Ang[i + 8]));
        if (sumValues[i] < 0) sumValues[i] = 0;
        sumValues[i + 4] = (p1Ang[i] + p1Ang[i + 4] * TMath::Power(fabs(alLow), p1Ang[i + 8]) + p1Ang[i + 12] + p1Ang[i + 16] * TMath::Power(alHigh, p1Ang[i + 20])) * 0.5 * (alHigh - alLow);
        if (sumValues[i + 4] < 0) sumValues[i + 4] = 0;
        sumValues[i + 8] = p1Ang[i + 12] * (p1Ang[i + 24] - alHigh) + p1Ang[i + 16] * (TMath::Power(p1Ang[i + 24], p1Ang[i + 20] + 1.) - TMath::Power(alHigh, p1Ang[i + 20] + 1)) / (p1Ang[i + 20] + 1);
        if (sumValues[i + 8] < 0) sumValues[i + 8] = 0;
        if (p1Ang[i + 24] < 1) { sumValues[i + 12] = (p1Ang[i + 12] + p1Ang[i + 16] * TMath::Power(p1Ang[i + 24], p1Ang[i + 20]) + p1Ang[i + 28]) * 0.5 * (1. - p1Ang[i + 24]);  if (sumValues[i + 12] < 0) sumValues[i + 12] = 0; }
        else sumValues[i + 12] = 0;
        sumTotal[i] = (sumValues[i] + sumValues[i + 4] + sumValues[i + 8] + sumValues[i + 12]) * 2. * TMath::Pi();
        if (sumTotal[i] < 1e-20) sumTotal[i] = 1e-20;
    }

    for (int i = 0; i < 4; i++)
    {
        if (cosTheta <= alLow) bAng[i] = p1Ang[i] + p1Ang[i + 4] * TMath::Power(fabs(cosTheta), p1Ang[i + 8]);
        else
        {
            if (cosTheta <= alHigh) bAng[i] = (p1Ang[i] + p1Ang[i + 4] * TMath::Power(fabs(alLow), p1Ang[i + 8])) * (1. - rightofCosTheta) + (p1Ang[i + 12] + p1Ang[i + 16] * TMath::Power(alHigh, p1Ang[i + 20])) * rightofCosTheta;
            else
            {
                if (cosTheta <= p1Ang[i + 24]) bAng[i] = p1Ang[i + 12] + p1Ang[i + 16] * TMath::Power((cosTheta), p1Ang[i + 20]);
                else  bAng[i] = (p1Ang[i + 12] + p1Ang[i + 16] * TMath::Power(p1Ang[i + 24], p1Ang[i + 20])) * (1. - (cosTheta - p1Ang[i + 24]) / (1. - p1Ang[i + 24])) + p1Ang[i + 28] * (cosTheta - p1Ang[i + 24]) / (1. - p1Ang[i + 24]);
            }
        }
    }

    cAng[0] = bAng[0] / sumTotal[0] * (1. - ratioRigidity) + bAng[1] / sumTotal[1] * ratioRigidity;
    cAng[1] = bAng[2] / sumTotal[2] * (1. - ratioRigidity) + bAng[3] / sumTotal[3] * ratioRigidity;

    double specAng = cAng[0] * (1. - ratioAtm) + cAng[1] * ratioAtm;
    if (specAng < 0) specAng = 0;

    // black hole correction
    double a9bh = 0.052682 + (-0.000002) / (1. + exp((0.51974 - log10(energy)) / 0.8592)) + 0.94732 / (1. + exp((0.53833 - log10(energy)) / 0.83756));
    double a10bh = 0.1397 + (0.41101) / (1. + exp((-0.95252 - log10(energy)) / 0.66911)) + 0.44929 / (1. + exp((0.1879 - log10(energy)) / 0.284));
    double a11bh = 1.4369 + (-1.1579) / (1. + exp((0.25629 - log10(energy)) / 0.32432));

    double bhFactor = a9bh + (a10bh - a9bh) * TMath::Power(cosTheta, a11bh);

    return specAng * bhFactor;
}


/**
 * variable: energy (MeV)
 * parameters atmDepth(g/cm^2), cutoff rigidity(GV), cosTheta
 * returns the angular dependance function of the CR spectrum
 * @param x
 * @param par
 * @returns - angular dependance
 */
Double_t  angularFunc(Double_t* x, Double_t* par)
{
    return calcParaAdep(par[0], par[1], x[0], par[2]);
}

//variables for angularAllFunc of the Sato 2016 parameter data set
double b1CorrSmin = -1, b2CorrSmin = -1, b3CorrSmin = -1, b4CorrSmin = -1, b5CorrSmin = -1, b6CorrSmin = -1, b7CorrSmin = -1, b8CorrSmin = -1, b9CorrSmin = -1;
double b1CorrSmax = -1, b2CorrSmax = -1, b3CorrSmax = -1, b4CorrSmax = -1, b5CorrSmax = -1, b6CorrSmax = -1, b7CorrSmax = -1, b8CorrSmax = -1, b9CorrSmax = -1;
double FFPSmax = -1, FFPSmin = -1;
double PhiiS = -1, c4n = -1, c12 = -1, powHere = -1;
bool alreadyCalculated = false;
double atmDensityAlready = -1, wValueAlready = -1, rigidityAlready = -1;


/**
 * variable: energy (MeV)
 * parameters atmDepth(g/cm^2), cutoff rigidity(GV), cosTheta
 * returns the value for the mean basic spectrum (phiBmean) times the solar cylcle correction (Phiis) times the upper/lower rigidity correction (correction)
 * times the angular dependance including the black hole correction (angleCorrection)
 * @param x
 * @param par
 * @returns - angular dependance including angleCorrection
 */
Double_t angularAllFunc(Double_t* x, Double_t* par)
{
    double atmDensityHere = par[0];
    double rigidityHere = par[1];

    //double cosTheta = par[2];

    double wValue = 20.; // solar activity

    double FFP = 370. + 0.3 * TMath::Power(wValue, 1.45);

    //calculate angular dependance
    double angleCorrection = calcParaAdep(par[0], par[1], x[0], par[2]);

    //the 'alreadyCalculated' just makes sure, that the set of parameters, which is static will not be loaded more than once
    if (!alreadyCalculated)
    {
        FFPSmax = 370. + 0.3 * TMath::Power(150, 1.45);
        FFPSmin = 370. + 0.3 * TMath::Power(0, 1.45);

        double c1 = (gji[0][0]) + (gji[1][0]) * rigidityHere + (gji[2][0]) / (1. + exp((rigidityHere - gji[3][0]) / (gji[4][0])));			// parameter: cut-off rigidity r_c
        double c2 = (gji[0][1]) + (gji[1][1]) * rigidityHere + (gji[2][1]) / (1. + exp((rigidityHere - gji[3][1]) / (gji[4][1])));
        double c3 = (gji[0][2]) + (gji[1][2]) * rigidityHere + (gji[2][2]) / (1. + exp((rigidityHere - gji[3][2]) / (gji[4][2])));
        double c4 = (gji[0][3]) + (gji[1][3]) * rigidityHere + (gji[2][3]) / (1. + exp((rigidityHere - gji[3][3]) / (gji[4][3])));

        double c5 = (gji[0][4]) + (gji[1][4]) * rigidityHere + (gji[2][4]) / (1. + exp((rigidityHere - gji[3][4]) / (gji[4][4])));
        double c9 = (gji[0][5]) + (gji[1][5]) * rigidityHere + (gji[2][5]) / (1. + exp((rigidityHere - gji[3][5]) / (gji[4][5])));
        double c10 = (gji[0][6]) + (gji[1][6]) * rigidityHere + (gji[2][6]) / (1. + exp((rigidityHere - gji[3][6]) / (gji[4][6])));
        double c11 = (gji[0][7]) + (gji[1][7]) * rigidityHere + (gji[2][7]) / (1. + exp((rigidityHere - gji[3][7]) / (gji[4][7])));

        double Phii = c1 * exp(-c2 * atmDensityHere) - c3 * exp(-c4 * atmDensityHere);

        // high altitude correction
        double c1H = (gji[5][0]) + (gji[6][0]) * rigidityHere + (gji[7][0]) / (1. + exp((rigidityHere - gji[8][0]) / (gji[9][0])));			// parameter: cut-off rigidity r_c
        double c2H = (gji[5][1]) + (gji[6][1]) * rigidityHere + (gji[7][1]) / (1. + exp((rigidityHere - gji[8][1]) / (gji[9][1])));
        double c3H = (gji[5][2]) + (gji[6][2]) * rigidityHere + (gji[7][2]) / (1. + exp((rigidityHere - gji[8][2]) / (gji[9][2])));
        double c4H = (gji[5][3]) + (gji[6][3]) * rigidityHere + (gji[7][3]) / (1. + exp((rigidityHere - gji[8][3]) / (gji[9][3])));

        // high altitude correction
        double PhiiH = c1H * exp(-c2H * atmDensityHere) - c3H * exp(-c4H * atmDensityHere);

        double b1 = (0.16 + 0.00842 * rigidityHere + (-1.43) / (1. + exp((rigidityHere - 2.86) / 1.62)));
        double b2 = (0.0119 + (-0.000352) * rigidityHere + (-0.0117) / (1. + exp((rigidityHere - 9.83) / 6.68)));

        powHere = b1 + b2 * atmDensityHere;

        double sa2 = (Phii - PhiiH) / (TMath::Power(FFPSmin, powHere) - TMath::Power(FFPSmax, powHere));
        double sa1 = Phii - sa2 * TMath::Power(FFPSmin, powHere);

        PhiiS = sa1 + sa2 * TMath::Power(FFP, powHere);

        c4n = c5 + 0.000198 * atmDensityHere / (1. + 0.567 * exp(0.00115 * atmDensityHere));
        c12 = c9 * (exp(-c10 * atmDensityHere)) + c11 * exp(-0.328 * atmDensityHere);

        double lowAtm, highAtm;

        //loads the set of parameters for the rigidity as well as for atmospheric depth minimum and maximum condition
        getPar(atmDensityHere, 0, 1); lowAtm = atmIndex[(int)retValue[0]]; highAtm = atmIndex[1 + (int)retValue[0]];
        double b1CorrSminLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b1CorrSminUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b1CorrSmin = b1CorrSminLowHeight + (atmDensityHere - lowAtm) * (b1CorrSminUppHeight - b1CorrSminLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 0, 2);
        double b2CorrSminLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b2CorrSminUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b2CorrSmin = b2CorrSminLowHeight + (atmDensityHere - lowAtm) * (b2CorrSminUppHeight - b2CorrSminLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 0, 3);
        double b3CorrSminLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b3CorrSminUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b3CorrSmin = b3CorrSminLowHeight + (atmDensityHere - lowAtm) * (b3CorrSminUppHeight - b3CorrSminLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 0, 4);
        double b6CorrSminLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b6CorrSminUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b6CorrSmin = b6CorrSminLowHeight + (atmDensityHere - lowAtm) * (b6CorrSminUppHeight - b6CorrSminLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 0, 5);
        double b7CorrSminLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b7CorrSminUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b7CorrSmin = b7CorrSminLowHeight + (atmDensityHere - lowAtm) * (b7CorrSminUppHeight - b7CorrSminLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 5, 1);
        double b1CorrSmaxLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b1CorrSmaxUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b1CorrSmax = b1CorrSmaxLowHeight + (atmDensityHere - lowAtm) * (b1CorrSmaxUppHeight - b1CorrSmaxLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 5, 2);
        double b2CorrSmaxLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b2CorrSmaxUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b2CorrSmax = b2CorrSmaxLowHeight + (atmDensityHere - lowAtm) * (b2CorrSmaxUppHeight - b2CorrSmaxLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 5, 3);
        double b3CorrSmaxLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b3CorrSmaxUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b3CorrSmax = b3CorrSmaxLowHeight + (atmDensityHere - lowAtm) * (b3CorrSmaxUppHeight - b3CorrSmaxLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 5, 4);
        double b6CorrSmaxLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b6CorrSmaxUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b6CorrSmax = b6CorrSmaxLowHeight + (atmDensityHere - lowAtm) * (b6CorrSmaxUppHeight - b6CorrSmaxLowHeight) / (highAtm - lowAtm);

        getPar(atmDensityHere, 5, 5);
        double b7CorrSmaxLowHeight = retValue[1] + retValue[2] * rigidityHere + retValue[3] / (1. + exp((rigidityHere - retValue[4]) / retValue[5]));
        double b7CorrSmaxUppHeight = retValue[6] + retValue[7] * rigidityHere + retValue[8] / (1. + exp((rigidityHere - retValue[9]) / retValue[10]));
        b7CorrSmax = b7CorrSmaxLowHeight + (atmDensityHere - lowAtm) * (b7CorrSmaxUppHeight - b7CorrSmaxLowHeight) / (highAtm - lowAtm);

        b4CorrSmin = 0.8341 + (-0.008648) * rigidityHere + (0.03596) / (1. + exp((rigidityHere - 7.603) / 0.4569));
        b4CorrSmax = 0.7672 + (-0.00863) * rigidityHere + (0.07556) / (1. + exp((rigidityHere - 8.185) / 0.4569));
        b5CorrSmin = 1.073 + (0.006603) * rigidityHere + (-0.5479) / (1. + exp((rigidityHere - 17.11) / 7.046));
        b5CorrSmax = 1.173 + (-0.009488) * rigidityHere + (-0.5682) / (1. + exp((rigidityHere - 9.954) / 2.61));

        b8CorrSmin = 3.573 + (-0.09026) * rigidityHere + (0.3519) / (1. + exp((rigidityHere - 8.899) / 0.5996));
        b9CorrSmin = 168.6 + (18.85) * rigidityHere + (-88.42) / (1. + exp((rigidityHere - 5.831) / 0.231));
        b8CorrSmax = 0.2814 + (0.06762) * rigidityHere + (2.756) / (1. + exp((rigidityHere - 9.428) / 2.317));
        b9CorrSmax = 1623. + (1.041) * rigidityHere + (-1550.) / (1. + exp((rigidityHere - 16.17) / 5.049));
    }
    if ((atmDensityHere > 0.98 * atmDensityAlready) && (atmDensityHere < 1.02 * atmDensityAlready) && (wValue > 0.98 * wValueAlready) && (wValue < 1.02 * wValueAlready) && (rigidityHere > 0.98 * rigidityAlready) && (rigidityHere < 1.02 * rigidityAlready))
    {
        alreadyCalculated = true;
    }
    else
    {
        atmDensityAlready = atmDensityHere;
        wValueAlready = wValue;
        rigidityAlready = rigidityHere;
        alreadyCalculated = false;
    }

    double phiBmean = 0.209 * TMath::Power(x[0] / 2.08, 0.675) * exp(-x[0] / 2.08) + c4n * exp(-TMath::Power(log10(x[0]) - log10(123.), 2) / (2. * TMath::Power(log10(2.01), 2))) + 0.000873 * log10(x[0] / 2.36 * TMath::Power(10., 14.)) * (1. + tanh(1.34 * log10(x[0] / 1.23 * TMath::Power(10., 7.)))) * (1. - tanh(1.19 * log10(x[0] / c12)));

    double sMinCorrection = TMath::Power(10, b1CorrSmin) + (b2CorrSmin * log10(x[0]) + b3CorrSmin) * (1. - tanh(b4CorrSmin * log10(x[0] / b5CorrSmin))) + (b6CorrSmin * log10(x[0]) + b7CorrSmin) * (1. + tanh(b8CorrSmin * log10(x[0] / b9CorrSmin)));
    double sMaxCorrection = TMath::Power(10, b1CorrSmax) + (b2CorrSmax * log10(x[0]) + b3CorrSmax) * (1. - tanh(b4CorrSmax * log10(x[0] / b5CorrSmax))) + (b6CorrSmax * log10(x[0]) + b7CorrSmax) * (1. + tanh(b8CorrSmax * log10(x[0] / b9CorrSmax)));

    double A2 = (sMinCorrection - sMaxCorrection) / (TMath::Power(FFPSmin, powHere) - TMath::Power(FFPSmax, powHere));
    double A1 = (sMinCorrection)-A2 * TMath::Power(FFPSmin, powHere);
    double correction = A1 + A2 * TMath::Power(FFP, powHere);

    return  PhiiS * phiBmean * correction * angleCorrection;
}



/**
 * generates normalized input spectrum from a modified cosmic ray spectrum with a number (entries) of random numbers
 * the cosmic ray spectrum is generated according to the formula given by Sato
 * then the reflected part is sperated from the incoming part
 * deprecated Sato 2008 version
 * @param precalculatedSpectrum
 * @param entries
 * @param file1
 */
void generateNormalizedSpectrum(TH1F* precalculatedSpectrum, int entries, TString file1)
{
    bool gotIt;
    int nBins;
    float xRnd, yRnd;
    float number;

    float backscatteredFactor = 1.05 * (1. - sourceLayerHeight / 1000. / 200. * 0.7);

    if (useVolumeSource) backscatteredFactor = 0;

    if (useHECascadeModel) backscatteredFactor = 1.25;

    TRandom3 r;
    time_t timeMsec = time(NULL);
    r.SetSeed(timeMsec);

    TString histoName1 = "scatteredSurfaceSpectrumRel";

    TString outputFolder = "";

    TFile f1(file1, "READ");

    TH1F* scatteredSurfaceSpectrumRel = (TH1F*)f1.Get(histoName1);

    TString phiBmean = "0.229*TMath::Power(x/2.31,0.721)*exp(-x/2.31)+0.0516*exp(-TMath::Power(log10(x)-log10(126),2)/(2*TMath::Power(log10(2.17),2)) ) + 0.00108*log10(x/3.33*TMath::Power(10.,12.)) * (1+tanh(1.62*log10(x/9.59*TMath::Power(10.,8.)))) *(1-tanh(1.48*log10(x/299.)))";
    //phiB = "1";

    //TF1* phiBFunc = new TF1("phiBFunc",phiB,0.0000000001,10000);

    TString a1 = "(12.9 + 15.7 / (1+exp( ([1]-5.62) /1.79) )  )";			// parameter: cut-off rigidity r_c
    TString a2 = "(0.00706 + 0.00057 / (1+exp( ([1]-5.99) /1.94) )  )";		// parameter: cut-off rigidity r_c
    TString a3 = "(0.975 - 0.210 / (1+exp( ([1]-0.99) /2.24) )  )";			// parameter: cut-off rigidity r_c
    TString a4 = "(0.0084 + 0.00441 / (1+exp( ([1]-2.24) /2.66) )  )";		// parameter: cut-off rigidity r_c

    TString a5 = "(-0.00701 + 0.0258 / (1+exp( ([1]-10.9) /2.38) )  )";		// parameter: cut-off rigidity r_c
    TString a9 = "(642 - 189 / (1+exp( ([1]-2.32) /0.897) )  )";			// parameter: cut-off rigidity r_c
    TString a10 = "(0.00112 + 0.000181 / (1+exp( ([1]-8.84) /0.587) )  )";	// parameter: cut-off rigidity r_c
    TString a11 = "(1.26 - 0.958 / (1+exp( ([1]-3.18) /1.47) )  )";			// parameter: cut-off rigidity r_c

    TString c4 = "(" + a5 + "+0.000171*[0]/(1+0.53*exp(0.00136*[0])) )";               // parameter: atmospheric density
    TString c12 = "(" + a9 + "*(exp(-" + a10 + "*[0]))+" + a11 + "*exp(-0.0133*[0]) )";        // parameter: atmospheric density

    TString phiL = "(" + a1 + "*(exp(-" + a2 + "*[0]) - " + a3 + "*exp(-" + a4 + "*[0])) )";		// parameter: atmospheric density

    TString g3 = "(-25.2+2.73/([2]+0.0715))";									//parameter: weight fraction of water w
    TString g5 = "(0.348+3.35*[2]-1.57*TMath::Power([2],2) )";					//parameter: weight fraction of water w

    TString phiB = "0.229*TMath::Power(x/2.31,0.721)*exp(-x/2.31)+" + c4 + "*exp(-TMath::Power(log10(x)-log10(126),2)/(2*TMath::Power(log10(2.17),2)) ) + 0.00108*log10(x/3.33*TMath::Power(10.,12.)) * (1+tanh(1.62*log10(x/9.59*TMath::Power(10.,8.)))) *(1-tanh(1.48*log10(x/" + c12 + ")))";

    TString fG = "TMath::Power(10.,-0.0235-0.0129*(log10(x)-" + g3 + ")*(1-tanh(0.969*log10(x/" + g5 + ") ) )  )";
    //TF1* fgFunc = new TF1("fgFunc", "[0]*0+[1]*0+"+fG,0.0000000001,10000); fgFunc->SetParameter(2,0.9);

    TF1* phiAllFunc = new TF1("phiAllFunc", phiL + "*(" + phiB + ")*(" + fG + ")", 0.000000001, 10000);

    TString phiFunctionString;
    if (useBasicSpectrum) phiFunctionString = phiL + "*(" + phiBmean + ")*(" + fG + ")";
    else phiFunctionString = phiL + "*(" + phiB + ")*(" + fG + ")";

    //TString phiT = "1./1.3E6 * x * TMath::Power(1./[0],2)*TMath::Exp(-x/[0])";
    // number density

    //TString phiT = "1/2.E6*2./sqrt(TMath::Pi() )*sqrt(x) * TMath::Power(1./[3],1.5)*exp(-x/[3])";
    TString phiT = "(0.118+0.144*TMath::Exp(-3.87*[2]))/(1.+0.653*TMath::Exp(-42.8*[2]))* TMath::Power(1./[3],2)*TMath::Exp(-x/[3])"; //should be E/[3]^2
    phiT = "0";

    //TF1* phiTFunc = new TF1("phiTFunc", phiT,0.0000000001,10000);
    //phiTFunc->SetParameter(0,0.025e-6);  //thermal rEnergy 25 meV

    //TF1* phiAllFunc = new TF1("phiAllFunc",phiL+"*("+phiB+")",0.0000000001,10000);
    if (useBasicSpectrum)
    {
        if (noThermalRegime)   phiAllFunc = new TF1("phiAllFunc", phiL + "*(" + phiBmean + ")*(" + fG + ")", 0.000000001, 10000);
        else
        {
            phiAllFunc = new TF1("phiAllFunc", phiL + "*((" + phiB + ")*(" + fG + ")+" + phiT + ")", 0.000000001, 10000);
            phiAllFunc->SetParameter(3, 0.025e-6);  //thermal rEnergy 25 meV
        }
    }
    else
    {
        if (!noThermalRegime)
        {
            phiAllFunc = new TF1("phiAllFunc", phiL + "*((" + phiB + ")*(" + fG + ")+" + phiT + ")", 0.000000001, 10000);
            phiAllFunc->SetParameter(3, 0.025e-6);  //thermal rEnergy 25 meV
        }
    }


    phiAllFunc->SetParameter(0, atmDensity);		//Atm Density
    phiAllFunc->SetParameter(1, rigidity);		// cut-off rigidity
    phiAllFunc->SetParameter(2, 0.99);		// water fraction

    //phiAllFunc->SaveAs(outputFolder+"/CosmicSpectrumFunction.root");

    TH1F* cosmicSpectrumHere = new TH1F("cosmicSpectrumHere", "Cosmic Spectrum", 10000, 1e-9, 10000);					TH1Fashion(cosmicSpectrumHere, "n", "Energy [MeV]", true);
    logaxis(cosmicSpectrumHere);
    TH1F* cosmicSpectrum2 = new TH1F("cosmicSpectrum2", "Cosmic Spectrum 2", 10000, 1e-9, 10000);						TH1Fashion(cosmicSpectrum2, "n", "Energy [MeV]", true);
    logaxis(cosmicSpectrum2);
    TH1F* cosmicSpectrum3 = new TH1F("cosmicSpectrum3", "Cosmic Spectrum 3", 10000, 1e-9, 10000);						TH1Fashion(cosmicSpectrum3, "n", "Energy [MeV]", true);
    logaxis(cosmicSpectrum3);

    nBins = cosmicSpectrumHere->GetNbinsX();

    float funcMax = phiAllFunc->GetMaximum(0.05, 800);

    for (Int_t k = 0; k < entries; k++)
    {
        gotIt = false;
        while (!gotIt)
        {
            if (noThermalRegime)  xRnd = TMath::Power(10, r.Rndm() * 12. - 7.7);
            else xRnd = TMath::Power(10, r.Rndm() * 12.3 - 8.);
            yRnd = r.Rndm() * funcMax;
            //if (xRnd < 2e-8) continue;
            if (phiAllFunc->Eval(xRnd) > yRnd)
            {
                gotIt = true;
            }
            //gotIt = true;
        }
        cosmicSpectrumHere->Fill(xRnd);
    }

    // for some reason for ROOT 5.34.23+ there is an undefined error here as soon as there is two times the operation GetBinContent on two different histograms. Reason ist absolutely unclear.
    for (Int_t k = 1; k < (nBins - 1); k++)
    {
        number = 0;

        if (cosmicSpectrum2->GetBinCenter(k) < (1e-7))
        {
            number = (1. - backscatteredFactor * 1.) * cosmicSpectrumHere->GetBinContent(k);
        }
        else
        {
            number = (1. - backscatteredFactor * scatteredSurfaceSpectrumRel->GetBinContent(k)) * cosmicSpectrumHere->GetBinContent(k);
        }
        if ((number > 1e-9) && (number < 1e9)) { ; }
        else  number = 0.;

        //number = 1-(1-number)/5.;
        cosmicSpectrum2->SetBinContent(k, number);
    }

    float allEntries = 0;
    float actualBin = 0;
    float largestBin = 0;
    for (Int_t k = 1; k < nBins; k++)
    {
        actualBin = cosmicSpectrum2->GetBinContent(k);
        allEntries += actualBin;
        if (actualBin > largestBin) largestBin = actualBin;
    }

    for (Int_t k = 1; k < precalculatedSpectrum->GetNbinsX(); k++)
    {
        actualBin = cosmicSpectrum2->GetBinContent(k);
        //cosmicSpectrum3->SetBinContent(k,cosmicSpectrum2->GetBinContent(k)/allEntries);
        if (actualBin > 0) cosmicSpectrum3->SetBinContent(k, actualBin / largestBin);
        else cosmicSpectrum3->SetBinContent(k, 0);
        if (actualBin > 0) precalculatedSpectrum->SetBinContent(k, actualBin / largestBin);
        else precalculatedSpectrum->SetBinContent(k, 0);
    }

    delete cosmicSpectrumHere;
    delete cosmicSpectrum2;
    delete cosmicSpectrum3;
    delete phiAllFunc;
    delete scatteredSurfaceSpectrumRel;
}




// ********************
// ********************
// MAIN ROUTINE
// ********************
// ********************
bool cosmicNSimulator(MainWindow* uiM)
{
    try {
        //using namespace Ui ;

        if (!noGUIMode) uiM->setStatus(1, "Activate Live Mode");
        if (!noGUIMode) {uiM->setStatus(2, "");        delay(10);}

        setupLiveTHs();

        // Parameter batch run
        //float paramMin = 0;
        //float paramMax = 60; //40 for detector energy analysis

        if ((doDetectorBatchRun) || ((doDetectorBatchRun2))) paramMax = 40;
        if (doDetectorAngleBatchRun) paramMax = 20;

        //float matrixMetricFactor = 0; // pixel conversion factor equals in meter
        bool drawDensityMap = true;

        //number density conversion to volumetric water content in soil
        const float soilStrechFactor = 1.37; // 1.321; changed to 1.37 22.07.2019

        //generate Geometry

        // [x0,y0,xEnd,yEnd,z0,height,material,layernumber]
        cout << "No. of layers: " << geometries.size() << endl;
        cout << "No. of neutrons: " << neutrons << endl;

        // ********************
        // Variables

        // timer variables
        time_t timeMsec = time(NULL);
        time_t start, end, diffmean, actualTime, oldTime, pause1, pause2;
        time_t startTemp, endTemp;

        // these are all variables which are used in one way or another, mainly in the calculation routine, some are deprecated (from the time when this simulation was called CosmicNeutronSimulator and was stripped down to URANOS)
        double xRnd, yRnd;
        double length, lengthZ, distToDet, distToCyl, energy, energy1e6, energy1e6Log, energyInitial, energyAtInterface, absDepth, energyNew, energyOld, maxgenerationHeight, generationHeight;
        double theta, thetaScat, thetaCMS, thetaFromTab, thetaAlt, thetaCalc, thetaOld, phiAlt, phi, phiScat, phiOld, nSpeed, timeTemp, trackTemp, timeTrans, timeNtr, pauseTime;
        double dcosTheta, theta12;
        double lambda, wwRangeAbs, wwRange, wwRangeSi, wwRangeIn, probSi, xTrack, xHere, yHere, zTrack, z0Here, yTrack, thetaTr, phiTr, energyTr, batchVar, trackMetricFactor, trackMetricFactorTemp, aspectRatioSide;
        double sumCs = 0, cs = 0, csAdd = 0, csH = 0, csO = 0, csSi = 0, csN = 0, csW = 0, csAl = 0, csFe = 0, csC = 0, csAr = 0, csB = 0, csB10 = 0, csB11 = 0, csHe3 = 0, csGd155 = 0, csGd157 = 0, csPb206 = 0, csPb207 = 0, csPb208 = 0, csPb = 0;
        double csF = 0, csNi = 0, csCr52 = 0, csCr53 = 0, csMn55 = 0, csCl = 0, csNa = 0, csTi48 = 0, csK = 0, csSolids = 0, csPlants = 0, csLuft = 0;
        double asAdd = 0, asAl = 0, asSi = 0, asH = 0, asO = 0, asN = 0, asC = 0, asB = 0, asB10 = 0, asB11 = 0, asF = 0, asGd155 = 0, asGd157 = 0, asAr = 0, asNi = 0, asCr52 = 0, asCr53 = 0;
        double asMn55 = 0, asPb206 = 0, asPb207 = 0, asPb208 = 0, asPb = 0, asCl = 0, asNa = 0, asTi48 = 0, asK = 0, asBoden = 0, asPlants = 0, asHe3 = 0, asAll = 0;
        //double inEVAl = 0, inEVSi = 0, inEVO = 0, inEVN = 0, inEVC = 0, inEVAr = 0, inEVBoden = 0,  inEVPlants = 0, csInH = 0, csInO = 0, csInSi = 0, csInN = 0, csInW = 0, csInAl = 0, csInAr = 0;
        double csIn = 0;
        double sumCsArray[3];
        float temp;
        float tmpEvapo, scale, attenuationLength, funcMax, funcMin, moistureAtInterface, weight, weightThermal, alpha;
        double tempDist, tempDist2, tempDist3, tempDist4, temp2, tempz, tempzEnd, tempRnd;
        double previousTheta, previousPhi, previousX, previousY, previousZ0, previousCosPhi, previousSinPhi, previousTanTheta, previousCosTheta, previousSoilX, previousSoilY, previousSoilZ0, thermalizedX, thermalizedY, thermalizedZ0;
        double r1, phi1, x, y, xt, yt, xtEnd, ytEnd, xStart, yStart, xAtInterface, yAtInterface, zAtInterface, xAlt, yAlt, xLastScattered, yLastScattered, zLastScattered, xySqrt, z0, z0Alt, gAlt, z0max, zPlane, yPlane, xPlane, vectorDirFactor, previousTrackEnergy, xIs, yIs;
        double cosPhiTr, sinPhiTr, tanThetaTr, cosPhi, sinPhi, cosTheta, tanTheta, inelasticEnergyLoss, tValue;
        double subsurfaceScatteringMean;
        float currentlayer, scatterLayer, scatterings = 0;
        bool stillFirstLayer, absorbBreak, continuedThisLayer, hasBeenInSoil, slowedDown, newLayer, hasbeenEvaporized, isEvaporationNeutron, measuredOnce, started, sAbsorbed, thermalized, useCumulativeCalc, hasBeenReorded, hasPassedSurface;
        bool detectorHit, gotIt, intoThermalization, check, reverseDir, reverseDirAlt, scattered, scatteredThisLayer, leaveLayer, scatteredInelastic, scatteredElastic, scatteredEvaporating, materialNotFound;
        int progressN, thermalizedLayer, inputMatrixValue, inputMatrixSelector, binNo, graphCounter, graphCounterMG, graphN, maxlayer, element, elementID, material, tempInt, tt, escapes = 0, hits = 0;
        int glayerCounter, ztrCounter, sideNumber, lCounter, gotItCounter, evaporationCounter = 0;
        int actualMaterial, selectedMaterial, startMaterial, actualMaterial2, selectedMaterial2, startMaterial2, actualMaterial3, selectedMaterial3, startMaterial3;
        int absMaterialNo, currentG, startingLayerHere;
        long loopNumber = -1, neutronTrackNeutrons = 0;
        bool domainWallHit, differentMaterial1, differentMaterial2, differentMaterial3, oneLastCheck, neutronAbsorbedbyDetector;
        bool doNormalCalc;
        bool foundSomething, enteredDetectorLayer = false, detectorLayerOverride = false, detectorOverride = false;

        const double pi = TMath::Pi();
        const double piHalf = 0.5 * TMath::Pi();

        // for drawing the neutron graphs
        vector<float> rgbValues;
        vector<TGraph*> graphs;
        vector<float> vNeutron;
        TGraph* neutronPath = new TGraph();
        TMultiGraph* multigraph = new TMultiGraph();

        // neutron track coordinate vectors for interpolating between points
        vector< double > neutronTrackCoordinates;
        vector< double > neutronTrackCoordinatesFullSet;
        vector< double > neutronTrackCoordinates2;

        // neutron coordinate vectors for simulation
        vector< vector<float*> > allNeutronCoordinates;
        vector<float> subsurfaceScatterings; subsurfaceScatterings.clear();
        vector<double> xCoord, yCoord, zCoord, eCoord;

        //x,y,z0,theta,phi,g, energy
        vector<double> nEvapoVector;


        //these vectors store values used in the cross section/inelastic variables calculation step;
        TMatrixF* angularProb;
        vector<float> allInelastics;
        vector<float> allInelastics11, inelasticEnergyLossVec11, allInelastics20, inelasticEnergyLossVec20;
        vector<TMatrixF*> allInelasticAngulars11, allInelasticAngulars20;
        vector<int> allInelasticElements11, allInelasticElements20;
        vector<float> inelasticEnergyLossVec;
        vector<int> allMaterialsElements, allPlantsElements, allInelasticElements, allBodenElements, allAbsorptionMaterialsElements;
        vector<float> allBoden, allPlants, allMaterials, allAbsorptionMaterials;
        vector<TMatrixF*> allInelasticAngulars;

        int debugOutput = 0;                    // activate debug output... mostly deprectated

        float thetaBeam = 0;                    // nadir angle in rad
        float beamX = 0;                        //squareDim*0.5; //x and y position of the beam in mm
        float beamY = 0;                        //squareDim*0.5;

        bool useUniformBeam = true;             //radial distribution function

        bool logSpectrum = false;
        bool usePrecalculatedSpectrum = false;  // use spectrum calculated by Sato
        bool usePrecalculatedAsciiSpectrum = false;  // use spectrum from source file
        bool useSato2016 = true;                // use spectrum from Sato 2016
        bool usePrecalculatedSato2016Spectrum = true;   // use spectrum from Sato 2016
        bool useHomogenousSpectrum = false;     //equal distribution for all spectrum bins
        //bool noThermalRegime = true;          // thermal cutoff
        bool noHighEnergyRegime = false;        // high energy cutoff
        // these options set up specific cases
        bool doTheWaterThing = false;
        bool doTheCoastalTransect = false;
        bool doTheZredaThing = false;
        bool doSingleDetector = true;


        evaporationFactor = 0.42;  //0.75 for Paul ? 1+0.1 for sigma  0.88 for tau to omega 0.97 for tau to omega10
        float heEvaporationFactor = 0.77;
        //for HE Cascade Model //0.33 //0.32 set for first calculations with 570m attenuation length //0.77 for 1450m attenuation length //0.91 for sigma // 0.72 for tau?
        // changed to 0.81 for 1250m attenuation length 12.12.2018
        // changed back to 0.77 for 1250m attenuation length 19.12.2018 //0.6 leads to an attenuation length of 830m //0.77 leads to an attenuation length of 1350 m
        int evaporationAdditions = 1;

        const int satoArrayLength = 200;
        double satoEnergyArray[satoArrayLength] = {0};
        double satoAngleArray[satoArrayLength] = {0};
        int satoCounter = satoArrayLength;

        //for drawing
        bool drawSingleNeutronGraphs = false;   // draw neutron paths
        TString neutronPicSubPath = "/NeutronPics";
        bool drawSingleNeutronPropagation = false;
        const float timeFrame = 2.25e-7;        //sec //2e-7eV is good for cosmics
        const int trackIterationCutoff = 400;
        const float nSpeedEnergyCutoff = 0.05;
        const int numberOfFrames = 1000;
        nDetectedNeutrons = 0;
        int timeRunRemainHours, timeRunRemainMinutes, timeRunRemainSeconds, nRunPerS;
        double nRPerS;
        int nspers;

        attenuationLength = 148;                // for artificial source in the material
        killedThermal = 0;

        bool useHumAtmProfile = true;           // exponential humidity profile

        int neutronTrackNeutronsToExport = 100;

        float DetectorCylinderRadius = 1000;    // in mm
        detectorHeight = geometries.at(detectorLayer)[4] + 0.5 * geometries.at(detectorLayer)[5]; // sets the detector height to the vertical center of the detector layer

        // for error estimations: very cross sections by number of standard deviations
        /*
        bool varyCrosssections = false;
        float sigmaAbsorbFactor = -2.;
        float sigmaElasticFactor = -2.;
        float sigmaInelasticFactor = -2.;
        float sigmaHeFactor = -2.;
        float sigmaHeEstimator = .1;
        */

        // for single detector counts
        int det1 = 0, det2 = 0, det3 = 0, det4 = 0, det5 = 0, det6 = 0, det7 = 0, det8 = 0, det9 = 0, det10 = 0, det11 = 0, det0 = 0, deta = 0, detb = 0, detc = 0, detd = 0, dete = 0;
        int	det12 = 0, det13 = 0, det14 = 0, det15 = 0, det16 = 0;

        // single detector spatial settings
        float xtD, ytD;

        bool detHitted, detHitSelected, recordedThisLayer;
        bool detectorRealisticallyHitted, layerRealisticallyHitted;

        if ((soilSolidFracVar > 0.49) && (soilSolidFracVar < 0.51)); else soilSolidFrac = soilSolidFracVar;

        bool useFastCutOff = false;     // cuts of high energy calculation

        // output streams
        ofstream detOutputFile;
        ofstream detLayerOutputFile;
        ofstream detTrackOutputFile;
        ofstream detLayerTrackOutputFile;
        ofstream allTrackOutputFile;

        TH1F* precalculatedSpectrum = new TH1F("precalculatedSpectrum", "precalculated Spectrum", 10000, 1e-9, 10000);	 	TH1Fashion(precalculatedSpectrum, "n", "Energy [MeV]", true);
        logaxis(precalculatedSpectrum);
        allTHs.push_back(precalculatedSpectrum);

        sourceLayerHeight = TMath::Abs((geometries.at(startingLayer)[4] + 0.5 * geometries.at(startingLayer)[5]) - geometries.at(groundLayer)[4]);

        // Sato 2008 function
        if (usePrecalculatedSpectrum)
        {
            if (!noGUIMode) {uiM->setStatus(1, "Generating Input Spectrum");            delay(5);}
            generateNormalizedSpectrum(precalculatedSpectrum, 500000, inputSpectrumFile);
        }

        // Sato 2016 function
        TF1* angAllFunction = new TF1("angAllFunction", angularAllFunc, 0.000000001, 10000, 3); angAllFunction->SetNpx(3000);
        angAllFunction->SetParameters(atmDensity, rigidity, TMath::Cos(0));

        vector<double> tAngleSplVec;
        vector<TGraph*> tGrSplVec;
        vector<TSpline3*> tSplVec;
        int supportPts = 300;   // number of points to calculate the cosmic neutron spectrum
        int supportAngles = 50; // number of angles to calculate the cosmic neutron spectrum

        if (usePrecalculatedSato2016Spectrum)
        {
            if (!noGUIMode) {uiM->setStatus(1, "Generating Sato 2016 Input Spectrum");            delay(5);}

            double value, energyvalue, result;
            double intSum, intSumAll = 0, xSpl,ySpl;

            for (int i = 0; i < supportAngles; i++)
            {
                tGrSplVec.push_back(new TGraph());
            }

            for (Int_t k = 0; k < supportAngles; k++)
            {
                value = k * TMath::Pi() / 2. / (supportAngles * 1.);
                angAllFunction->SetParameters(atmDensity, rigidity, TMath::Cos(value));

                tAngleSplVec.push_back(value);

                intSum = 0;

                //ofstream dOutputFile;
                //if (k==0)    dOutputFile.open (outputFolder+"spectrum.dat", ios::out | ios::app);

                for (Int_t l = 0; l < supportPts; l++)
                {
                    energyvalue = (l * 1.) / (supportPts * 1.) * 12. - 7.7;
                    result = angAllFunction->Eval(TMath::Power(10, energyvalue));
                    intSum += result;

                    //if (k==0)   dOutputFile<<energyvalue<<"\t"<<result<<endl;

                    tGrSplVec.at(k)->SetPoint(l, energyvalue, intSum);
                }
                intSumAll += intSum;
                // if (k==0)    dOutputFile.close();

                for (int m = 0; m < supportPts; m++)
                {
                    tGrSplVec.at(k)->GetPoint(m, xSpl, ySpl);
                    tGrSplVec.at(k)->SetPoint(m, ySpl / intSum, xSpl);
                }
            }

            TString stringName;
            for (int i = 0; i < supportAngles; i++)
            {
                stringName = "angSp" + castIntToString(i);
                tSplVec.push_back(new TSpline3(stringName, tGrSplVec.at(i))); tSplVec.at(i)->SetNpx(1000);
            }

            neutronRealScalingFactor = intSumAll*2.*pi*100000./(supportPts*supportAngles*1.) ; //n/m^2/s but with 4pi or 2pi?

            if (useRadialBeam) neutronRealScalingFactor = neutronRealScalingFactor * (beamRadius*beamRadius*pi); //n/s
            else neutronRealScalingFactor = neutronRealScalingFactor * (beamRadius*beamRadius*4.);
        }

        // function to generate a spectrum
        // deprecated at the moment mostly
        // nice to keep here, will be killed when there's time to tidy up
        TString phiB = "0.229*TMath::Power(x/2.31,0.721)*exp(-x/2.31)+0.0516*exp(-TMath::Power(log10(x)-log10(126),2)/(2*TMath::Power(log10(2.17),2)) ) + 0.00108*log10(x/3.33*TMath::Power(10,12)) * (1+tanh(1.62*log10(x/9.59*TMath::Power(10,8)))) *(1-tanh(1.48*log10(x/299.)))";
        //phiB = "1";
        TF1* phiBFunc = new TF1("phiBFunc", phiB, 0.0000000001, 10000);

        //from Sato 2008 paper
        TString phiT = "1./0.05 * TMath::Power(x/[0],2)*TMath::Exp(-x/[0])"; //changed to this one on 18.01.2018
        // thermal flux
        // TString phiT = "1./1.3E6 * x * TMath::Power(1./[0],2)*TMath::Exp(-x/[0])";
        // TString phiT = "1./1.3E6 * TMath::Power(x/[0],2)*TMath::Exp(-x/[0])";
        // number density
        //TString phiT = "1/2.E6*2./sqrt(TMath::Pi() )*sqrt(x) * TMath::Power(1./[0],1.5)*TMath::Exp(-x/[0])";

        TF1* phiTFunc = new TF1("phiTFunc", phiT, 0.0000000001, 10000);
        phiTFunc->SetParameter(0, 0.0253e-6);  //thermal Energy 25 meV

        //Thermal Neutron Spectrum
        TF1* spectrumMaxwellPhi = new TF1("spectrumMaxwellPhi", "[0]/(TMath::Power([1],2))*x*TMath::Exp(-x/[1])", 0.1e-9, 0.4e-6); spectrumMaxwellPhi->SetNpx(2000);
        spectrumMaxwellPhi->SetParameters(0.67e-7, 0.0253e-6);

        //Thermal Neutron Spectrum with fixed parameters, normalized to max(y) = 1
        TF1* spectrumMaxwellPhiLinMod = new TF1("spectrumMaxwellPhiLinMod", "105000000.*x*TMath::Exp(-x/0.0253e-6)", 0.1e-9, 0.4e-6); spectrumMaxwellPhiLinMod->SetNpx(10000);

        //Spectrum of the CF Source
        TF1* spectrumMaxwellPhiTempModMod = new TF1("spectrumMaxwellPhiTempModMod", "[0]/(TMath::Power([1]/11604.5,2))*TMath::Power(x,2)*TMath::Exp(-x/[1]*11604.5)+[2]*TMath::Power([1]/11604.5,1)/(1+TMath::Power(([3]/x),7)) * 1./(1-[4]/(1+TMath::Power(x/[5],5)))", 1e-3, 1e1);
        spectrumMaxwellPhiTempModMod->SetNpx(10000);
        spectrumMaxwellPhiTempModMod->SetParameters(6868, 243, 9938, 0.07228, 0.5978, 0.2027);


        TSpline3* spectrumModel;
        if (inputSpectrumFile.Length() > 8)
        {
            TFile testFile(inputSpectrumFile, "READ");
            if (testFile.IsZombie()) //TODO: change to ASCII detection
            {
                spectrumModel = uiM->getSplinedDetectorEnergyModelFromFile(inputSpectrumFile, 0, true, true);
                cout << "Using spectrum model from: " << inputSpectrumFile << endl;
                usePrecalculatedAsciiSpectrum = true;       // use spectrum from source file
                useSato2016 = false;                        // use spectrum from Sato 2016
                usePrecalculatedSato2016Spectrum = false;   // use spectrum from Sato 2016
                usePrecalculatedSpectrum = false;
            }
            else
            {
                testFile.Close();
            }
        }

        if (!noGUIMode) {uiM->setStatus(1, "Generating Detector Model");        delay(5);}

        // detector response function
        // spline for He3 Rover scaled to full efficiency of 1 on LOG10 (energy)!
        TSpline3* detectorEnergyModel;
        if (detectorResponseFunctionFile.Length() < 8)
        {
            detectorEnergyModel = getSplinedDetectorEnergyModel();
        }
        else
        {
            detectorEnergyModel = uiM->getSplinedDetectorEnergyModelFromFile(detectorResponseFunctionFile, 0, true, true);
            cout << "Using energy model from: " << detectorResponseFunctionFile << endl;
        }

        TSpline3* detectorEnergyModel2;
        if (detectorResponseFunctionFile2.Length() > 6)
        {
            if (useAdditionalDetectorModel)
            {
                detectorEnergyModel2 = uiM->getSplinedDetectorEnergyModelFromFile(detectorResponseFunctionFile2, 0, true, true);
                cout << "Using additional energy model from: " << detectorResponseFunctionFile2 << endl;
                //useAdditionalDetectorModel = true;
            }
        }

        TF1* detectorAngleModel = new TF1("detectorAngleModel", "[0]+[1]*TMath::Exp(-x/[2])", 0, 1.6);
        detectorAngleModel->SetParameters(1.32171, -0.2542, -0.92);
        detectorAngleModel->SetParNames("Offset", "Base", "Tau");

        if (!noGUIMode) {uiM->setStatus(1, "Generating Histograms");        delay(5);}

        // Histograms to be filled during runtime
        TH1F* scatteredSurfaceSpectrumRel = new TH1F("scatteredSurfaceSpectrumRel", "Relative Scattered Neutron at Surface Spectrum", 10000, 1e-9, 10000);  TH1Fashion(scatteredSurfaceSpectrumRel, "n", "Energy [MeV]", setSize);
        TH1F* scatteredSurfaceSpectrumOnce = new TH1F("scatteredSurfaceSpectrumOnce", "Scattered Neutron at Surface Spectrum Once", 10000, 1e-9, 10000);  TH1Fashion(scatteredSurfaceSpectrumOnce, "n", "Energy [MeV]", setSize);
        TH1F* scatteredSurfaceSpectrumMultiplicity = new TH1F("scatteredSurfaceSpectrumMultiplicity", "Scattered Neutron at Surface Spectrum Multiplicity", 10000, 1e-9, 10000);  TH1Fashion(scatteredSurfaceSpectrumMultiplicity, "n", "Energy [MeV]", setSize);
        TH1F* surfaceSpectrumLikewise = new TH1F("surfaceSpectrumLikewise", "Neutrons at Surface Spectrum like Paper", 10000, 1e-9, 10000);  TH1Fashion(surfaceSpectrumLikewise, "n", "Energy [MeV]", setSize);

        //BinLogX(cosmicSpectrum);
        // creates logarithmically binned histograms from equally spaced bins
        logaxis(cosmicSpectrum);
        logaxis(scatteredSurfaceSpectrum);	logaxis(scatteredSurfaceSpectrumBack);	logaxis(scatteredSurfaceSpectrumRel); logaxis(scatteredSurfaceSpectrumOnce); logaxis(scatteredSurfaceSpectrumMultiplicity); logaxis(surfaceSpectrumLikewise);
        logaxis(scatteredSurfaceSpectrumHelp);
        // all root objects are added to this vector to store them after the run
        //allTHs.push_back(cosmicSpectrum); allTHs.push_back(scatteredSurfaceSpectrum); allTHs.push_back(scatteredSurfaceSpectrumBack);
        allTHs.push_back(scatteredSurfaceSpectrumOnce);  allTHs.push_back(scatteredSurfaceSpectrumRel); allTHs.push_back(scatteredSurfaceSpectrumMultiplicity); allTHs.push_back(surfaceSpectrumLikewise);
        //allTHs.push_back(scatteredSurfaceSpectrumHelp);

        TH1F* thetas1 = new TH1F("thetas1", "theta distribution neutrons leaving the ground", 1000, 0, 1.1 * pi);			TH1Fashion(thetas1, "n", "Angle", setSize);
        TH1F* thetas2 = new TH1F("thetas2", "1/sin(theta) distribution neutrons leaving the ground", 1000, 0, 1.1 * pi);    TH1Fashion(thetas2, "n", "Angle", setSize);
        //TH2F* thetas3 = new TH2F("thetas3", "theta Distribution 3" ,1000,0,1.1*pi, 1000,0,2);                           TH2Fashion(thetas3, "Angle", "range", setSize);
        TH2F* thetaEnergy = new TH2F("thetaEnergy", "Energy vs theta distribution neutrons leaving the ground", 1000, 0, 1.1 * pi, 1000, 1e-9, 10000);				TH2Fashion(thetaEnergy, "Angle", "Energy [MeV]", setSize);
        logYaxis(thetaEnergy);
        allTHs.push_back(thetas1); allTHs.push_back(thetas2); //allTHs.push_back(thetas3);
        allTHs.push_back(thetaEnergy);

        //TH1F* slowDownDepth = new TH1F("slowDownDepth", "slow Down Depth" ,5000,0,5000);

        TH1F* absDist = new TH1F("absDist", "Absorbed Layer", 135, -.5, 134.5);                                           TH1Fashion(absDist, "n", "Layer", setSize);
        TH1F* scatDist = new TH1F("scatDist", "Scattered Layer", 135, -.5, 134.5);                                        TH1Fashion(scatDist, "n", "Layer", setSize);
        TH1F* scatteringNs = new TH1F("scatteringNs", "Number of Scatterings per Neutron", 301, -.5, 300.5);              TH1Fashion(scatteringNs, "n", "Number of scatterings", setSize);
        TH1F* scatLengths = new TH1F("scatLengths", "Scattering Lengths ", 1000, -1., 150.);                              TH1Fashion(scatLengths, "n", "length [m]", setSize);
        //allTHs.push_back(slowDownDepth);
        //allTHs.push_back(scatDepth);
        allTHs.push_back(scatteringNs);

        //deprecated
        allTHs2.push_back(scatLengths);
        allTHs2.push_back(scatDist); allTHs2.push_back(absDist);

        TH1F* absEnergy = new TH1F("absEnergy", "Absorbed Neutron Energy", 10000, 1e-9, 10000);                           TH1Fashion(absEnergy, "n", "Energy [MeV]", setSize);
        logaxis(absEnergy);
        TH1D* absMaterial = new TH1D("absMaterial", "Absorbed Material", 33, 6.5, 39.5);                                  TH1Fashion(absMaterial, "n", "Material No", setSize);
        TH1D* scatMaterial = new TH1D("scatMaterial", "Scattered Material", 33, 6.5, 39.5);                               TH1Fashion(scatMaterial, "n", "Material No", setSize);
        TH1D* scatElement = new TH1D("scatElement", "Scattered Element", 60, -.5, 59.5);                                  TH1Fashion(scatElement, "n", "Element", setSize);
        TH1D* absElement = new TH1D("absElement", "Absorbed Element", 60, -.5, 59.5);                                     TH1Fashion(absElement, "n", "Element", setSize);
        TH1D* absElementL1 = new TH1D("absElementL1", "Absorbed Element", 60, -.5, 59.5);									TH1Fashion(absElementL1, "n", "Element", setSize);
        TH1D* absElementL2 = new TH1D("absElementL2", "Absorbed Element", 60, -.5, 59.5);									TH1Fashion(absElementL2, "n", "Element", setSize);
        TH1D* absElementL3 = new TH1D("absElementL3", "Absorbed Element", 60, -.5, 59.5);									TH1Fashion(absElementL3, "n", "Element", setSize);

        TH1F* energyLoss = new TH1F("energyLoss", "Scattered Neutron Energy Loss", 10000, 1e-9, 10000);                   TH1Fashion(energyLoss, "n", "Energy [MeV]", setSize);
        logaxis(energyLoss);
        allTHs.push_back(absEnergy); allTHs.push_back(absMaterial); allTHs.push_back(scatMaterial);  allTHs.push_back(absElement); allTHs.push_back(scatElement);
        allTHs2.push_back(absElementL1); allTHs2.push_back(absElementL2); allTHs2.push_back(absElementL3);
        allTHs.push_back(energyLoss);

        TH1F* detectorSpectrum = new TH1F("detectorSpectrum", "Detector Energy Spectrum", 1000, 1e-9, 10000);             TH1Fashion(detectorSpectrum, "n", "Energy [MeV]", setSize);
        logaxis(detectorSpectrum);
        TH1F* detectorSpectrumNear = new TH1F("detectorSpectrumNear", "Detector Energy Spectrum all Neutrons", 1000, 1e-9, 10000);			TH1Fashion(detectorSpectrumNear, "n", "Energy [MeV]", setSize);
        logaxis(detectorSpectrumNear);
        TH1F* detectorSpectrumOrigin = new TH1F("detectorSpectrumOrigin", "Detector Energy Neutron Origin Spectrum", 1000, 1e-9, 10000);			TH1Fashion(detectorSpectrumOrigin, "n", "Energy [MeV]", setSize);
        logaxis(detectorSpectrumOrigin);
        TH1F* detectorSpectrumOriginMoisture = new TH1F("detectorSpectrumOriginMoisture", "Detector Neutron Origin Moisture", 201, -0.0025, 1.0025);			TH1Fashion(detectorSpectrumOriginMoisture, "n", "Soil Moisture", setSize);

        TH1F* detectorLayerTheta = new TH1F("detectorLayerTheta", "Angle Distribution incoming Neutrons for Detector Layer", 100, 0., pi * 1.25);	TH1Fashion(detectorLayerTheta, "n", "Angle Theta", setSize);
        TH1F* detectorLayerTheta2 = new TH1F("detectorLayerTheta2", "Angle Distribution incoming Neutrons for Detector Layer Normalized to theta", 100, 0., pi * 1.25);	TH1Fashion(detectorLayerTheta2, "n", "Angle Theta", setSize);

        TH1F* detectorDistanceWeighted = new TH1F("detectorDistanceWeighted", "Weighted Distance Distribution for Detector", 4000, -1000, squareDim * 0.5);	TH1Fashion(detectorDistanceWeighted, "n (weighted)", "Distance [mm]", setSize);
        TH1F* detectorTheta = new TH1F("detectorTheta", "Angle Distribution incoming Neutrons for Detector", 100, 0., pi * 1.25);	TH1Fashion(detectorTheta, "n", "Angle Theta", setSize);
        TH1F* detectorTheta2 = new TH1F("detectorTheta2", "Angle Distribution incoming Neutrons for Detector Normalized to theta", 100, 0., pi * 1.25);	TH1Fashion(detectorTheta2, "n", "Angle Theta", setSize);
        TH1F* detectorPhi = new TH1F("detectorPhi", "Angle Distribution incoming Neutrons for Detector", 100, 0., 2. * pi * 1.25);	TH1Fashion(detectorPhi, "n", "Angle Phi", setSize);
        TH1F* detectorPhi2 = new TH1F("detectorPhi2", "Angle Distribution incoming Neutrons for Detector with Soil Contact", 100, -pi, 1. * pi * 1.25);	TH1Fashion(detectorPhi2, "n", "Angle Phi", setSize);
        //TH2F* detectorPhiEnergy = new TH2F("detectorPhiEnergy", "Energy vs Phi Distribution for Detector" ,100,0.,2.*pi*1.25, 1000,1e-9,10000);				TH2Fashion(detectorPhiEnergy, "Angle Phi", "Energy [MeV]", setSize);
        //logYaxis(detectorPhiEnergy);
        //TH2F* detectorThetaEnergy = new TH2F("detectorThetaEnergy", "Energy vs Theta Distribution for Detector" ,100,0.,pi*1.25, 1000,1e-9,10000);				TH2Fashion(detectorThetaEnergy, "Angle Theta", "Energy [MeV]", setSize);
        //logYaxis(detectorPhiEnergy);
        TH2F* detectorDistanceDepth = new TH2F("detectorDistanceDepth", "Distance vs Last Depth Distribution for Detector", 1500, -1000, squareDim * 0.59, 1500, 0.01, 1500);		TH2Fashion(detectorDistanceDepth, "Distance [mm]", "Depth [mm]", setSize);
        TH2F* detectorDistanceDepth2 = new TH2F("detectorDistanceDepth2", "Distance vs Depth Distribution for Detector", 1500, -1000, squareDim * 0.59, 1500, 0.01, 1500);		TH2Fashion(detectorDistanceDepth2, "Distance [mm]", "Depth [mm]", setSize);
        TH2F* detectorDistanceMaxDepth = new TH2F("detectorDistanceMaxDepth", "Distance vs Maximum Depth Distribution for Detector", 1500, -1000, squareDim * 0.59, 1500, 0.01, 1500);		TH2Fashion(detectorDistanceMaxDepth, "Distance [mm]", "Depth [mm]", setSize);
        TH2F* detectorDistanceEnergy = new TH2F("detectorDistanceEnergy", "Energy vs Distance Distribution", 4000, -1000, squareDim * 0.59, 1000, 1e-9, 10000);				TH2Fashion(detectorDistanceEnergy, "Distance", "Energy [MeV]", setSize);
        logYaxis(detectorDistanceEnergy);
        TH1F* detectorTimeTrans = new TH1F("detectorTimeTrans", "Detector Neutron Transport Time", 900000, 0, 900000);                   TH1Fashion(detectorTimeTrans, "n", "Time [mus]", setSize);
        TH1F* detectorLayerTimeTrans = new TH1F("detectorLayerTimeTrans", "Detector Layer Neutron Transport Time", 900000, 0, 900000);                   TH1Fashion(detectorLayerTimeTrans, "n", "Time [mus]", setSize);

        //TH2F* depthEnergyScattered = new TH2F("depthEnergyScattered", "Energy vs Depth Distribution for Penetration" ,125,-0.5,124.5, 1000,1e-9,10000);				TH2Fashion(depthEnergyScattered, "Layer", "Energy [MeV]", setSize);
        //logYaxis(depthEnergyScattered);
        TH2F* depthEnergyScattered2 = new TH2F("depthEnergyScattered2", "Energy vs Depth Distribution for Penetration", 2500, 0.1, 2500.5, 1000, 1e-9, 10000);				TH2Fashion(depthEnergyScattered2, "Depth [mm]", "Energy [MeV]", setSize);
        logYaxis(depthEnergyScattered2);
        TH2F* depthEnergyAbsorbed2 = new TH2F("depthEnergyAbsorbed2", "Energy vs Depth Distribution for Absorption", 2500, 0.001, 2500.5, 1000, 1e-9, 10000);				TH2Fashion(depthEnergyAbsorbed2, "Depth [mm]", "Energy [MeV]", setSize);
        logYaxis(depthEnergyAbsorbed2);
        //TH2F* depthEnergyAbsorbed = new TH2F("depthEnergyAbsorbed", "Energy vs Depth Distribution for Absorption" ,125,-0.5,124.5, 1000,1e-9,10000);				TH2Fashion(depthEnergyAbsorbed, "Layer", "Energy [MeV]", setSize);
        //logYaxis(depthEnergyAbsorbed);
        TH2F* depthThermalized = new TH2F("depthThermalized", "Energy vs Depth of Thermalization reaching surface", 2500, 0.001, 2500.5, 1000, 1e-9, 10000);				TH2Fashion(depthThermalized, "Depth [mm]", "Energy [MeV]", setSize);
        logYaxis(depthThermalized);

        TH2F* energyAbsorbedCosmic = new TH2F("energyAbsorbedCosmic", "Energy Original vs Energy Absorbed", 1000, 1e-9, 10000, 1000, 1e-9, 10000);						TH2Fashion(energyAbsorbedCosmic, "Energy [MeV]", "Energy [MeV]", setSize);
        logXaxis(energyAbsorbedCosmic); logYaxis(energyAbsorbedCosmic);

        TH1F* scatteredSurfaceDetectionDistance = new TH1F("scatteredSurfaceDetectionDistance", "Distance Distribution on Surface for Detection", 1000, 0, squareDim * 0.5);			TH1Fashion(scatteredSurfaceDetectionDistance, "n", "Distance [mm]", setSize);
        TH1F* scatteredSurfaceLayer = new TH1F("scatteredSurfaceLayer", "Scattered Layer Neutron on Surface", 135, -.5, 134.5);					TH1Fashion(scatteredSurfaceLayer, "n", "Layer", setSize);
        TH1F* scatteredSurfaceMaxLayer = new TH1F("scatteredSurfaceMaxLayer", "Scattered Deepest Layer Neutron on Surface", 135, -.5, 134.5);		TH1Fashion(scatteredSurfaceMaxLayer, "n", "Layer", setSize);

        TH1F* scatteredLayer = new TH1F("scatteredLayer", "Scattered Layer Neutron", 22, -.5, 21.5);					TH1Fashion(scatteredLayer, "n", "Layer", setSize);
        allTHs2.push_back(scatteredLayer);

        allTHs.push_back(detectorTimeTrans);  allTHs.push_back(detectorLayerTimeTrans);
        //logaxis(scatteredSurfaceDepth);
        allTHs.push_back(detectorSpectrum); allTHs.push_back(detectorSpectrumNear); allTHs.push_back(detectorSpectrumOrigin);  allTHs.push_back(detectorSpectrumOriginMoisture);
        //allTHs.push_back(detectorDistance);
        allTHs.push_back(scatteredSurfaceDetectionDistance); allTHs.push_back(detectorTheta);  allTHs.push_back(detectorPhi); allTHs.push_back(detectorTheta2);  allTHs.push_back(detectorLayerTheta);  allTHs.push_back(detectorLayerTheta2);// allTHs.push_back(detectorPhiEnergy); allTHs.push_back(detectorThetaEnergy);
        allTHs.push_back(detectorPhi2);
        //allTHs.push_back(scatteredSurfaceDistance);

        //allTHs.push_back(scatteredSurfaceDepth);
        allTHs.push_back(depthEnergyAbsorbed2);  allTHs.push_back(depthEnergyScattered2); allTHs.push_back(detectorDistanceDepth); allTHs.push_back(detectorDistanceDepth2); allTHs.push_back(detectorDistanceMaxDepth); allTHs.push_back(detectorDistanceEnergy);
        //allTHs.push_back(depthEnergyScattered);
        //allTHs.push_back(depthEnergyAbsorbed);
        //allTHs.push_back(scatteredSurfaceMaxDepth);
        allTHs.push_back(energyAbsorbedCosmic); allTHs.push_back(depthThermalized);
        //allTHs.push_back(detectorDistanceBackScattered);

        //Deprecated
        allTHs2.push_back(scatteredSurfaceLayer);
        allTHs2.push_back(detectorDistanceWeighted);
        allTHs2.push_back(scatteredSurfaceMaxLayer);

        TGraphErrors* quantile = new TGraphErrors();		TGraphErrorFashion(quantile, "Distance [mm]", "Quantil", true);	quantile->SetTitle("Quantiles in Distance");
        quantile->SetMarkerSize(2);
        allTHs.push_back(quantile);

        materialVector.clear();

        // reads input matrix definitions from files
        // can be a symmetric ASCII matrix or grayscale images
        // examples: a flower or a jumping goat
        if (useImage)
        {
            inputPicVector.clear(); inputPicVector2.clear(); inputPicVector3.clear();

            for (int ic = 0; ic < maxLayersAllowed; ic++)
            {
                TMatrixF* imageMatrix = new TMatrixF(inputPicSizes[ic], inputPicSizes[ic]);
                TMatrixF* imageMatrix2 = new TMatrixF(inputPicSizes[ic], inputPicSizes[ic]);
                TMatrixF* imageMatrix3 = new TMatrixF(inputPicSizes[ic], inputPicSizes[ic]);
                inputPicVector.push_back(*imageMatrix);
                inputPicVector2.push_back(*imageMatrix2);
                inputPicVector3.push_back(*imageMatrix3);
            }

            if (!noGUIMode) {uiM->setStatus(1, "Reading Matrices");   delay(5);}

            TString add = ""; // the additional letter after the layer number

            int picNo = 0;

            {uiM->setStatus(1, "Reading ASCII and png maps...");            delay(5);}

            for (int i = 0; i < model->rowCount(); i++)
            {
                picNo = i + 1;

                if (inputPics[i] == 1) add = "";
                if (inputPics[i] == 2) add = "";

                if ((inputPics[i] == 1) || (inputPics[i] == 3) || (inputPics[i] == 5))
                {
                    //transferMatrix  = readmatrix(workFolder, castIntToString(picNo), "dat", -1, inputMatrixPixels); turnInputMatrix(inputPicVector.at(i));
                    inputPicVector.at(i) = readmatrix(workFolder, castIntToString(picNo) + add, "dat", -1, inputPicSizes[i]); turnInputMatrix(inputPicVector.at(i));
                     //status bar message
                    string tpp = "#" + castIntToString(picNo);
                    if (!noGUIMode) {uiM->setStatus(2, tpp);   delay(2);}
                }

                if ((inputPics[i] == 2) || (inputPics[i] == 4) || (inputPics[i] == 6))
                {

                    inputPicVector.at(i) = readMatrixPNG(workFolder, castIntToString(picNo) + add); //TMatrixF mm =  readMatrixPNG(workFolder, castIntToString(picNo));
                    string tpp = "#" + castIntToString(picNo);
                    if (!noGUIMode) { uiM->setStatus(2, tpp);   delay(2);}

                }

                if (inputPics2[i] == 1) add = "d";
                if (inputPics2[i] == 2) add = "d";

                if (inputPics2[i] == 1)
                {
                    inputPicVector2.at(i) = readmatrix(workFolder, castIntToString(picNo) + add, "dat", -1, inputPicSizes[i]); turnInputMatrix(inputPicVector2.at(i));
                    string tpp = "#" + (string)castIntToString(picNo);
                    if (!noGUIMode) {uiM->setStatus(2, tpp);   delay(2);}
                }

                if (inputPics2[i] == 2)
                {
                    inputPicVector2.at(i) = readMatrixPNG(workFolder, castIntToString(picNo) + add);
                    string tpp = "#" + (string)castIntToString(picNo);
                    if (!noGUIMode) {uiM->setStatus(2, tpp);   delay(2);}
                }

                if (inputPics3[i] == 1) add = "p";
                if (inputPics3[i] == 2) add = "p";

                if (inputPics3[i] == 1)
                {
                    inputPicVector3.at(i) = readmatrix(workFolder, castIntToString(picNo) + add, "dat", -1, inputPicSizes[i]); turnInputMatrix(inputPicVector3.at(i));
                    string tpp = "#" + (string)castIntToString(picNo);
                    if (!noGUIMode) { uiM->setStatus(2, tpp);   delay(2);}
                }

                if (inputPics3[i] == 2)
                {
                    inputPicVector3.at(i) = readMatrixPNG(workFolder, castIntToString(picNo) + add);
                    string tpp = "#" + (string)castIntToString(picNo);
                    if (!noGUIMode) {uiM->setStatus(2, tpp);   delay(2); }
                }

            }

            matrixMetricFactor = squareDim / 1000. / (inputMatrixPixels * 1.);
        }
        if (!noGUIMode) uiM->setStatus(2, "");

        // initialize random generator
        r.SetSeed(time(NULL));

        soilWaterFrac = soilWaterFracVar;


        //Cross sections
        // not used when reading from matrix, however, some variables are used throughot the code when calculating temporary results
        //absorption
        //Float_t asAlu = 0.23, asCu = 3.8, asKapton = 7.8, asBoron = 3835., asLiF = 940., asPoly = 11.3, asEStahl = 3.1, asArCO2 = 0.541, asZirkon = 0.185, asBronze = 3.6, asBoronNat = 767., asHDPE = .66;
        Float_t asQuarz = 0.344, asAl2O3 = 0.64, asFe = 0.,asLuft = 3.13, asWater = 0.66,  asSolids = 0.;

        //element identifiers
        const int nH1 = 11, nHe3 = 23, nLi6 = 36, nB10 = 510, nB11 = 511, nC12 = 612, nN14 = 714, nO16 = 816, nF19 = 919, nNa23 = 1123, nAl27 = 1327, nSi28 = 1428, nS32 = 1632, nCl35 = 1735, nK39 = 1939, nTi48 = 2248;
        const int nAr40 = 1840, nCr52 = 2452, nCr53 = 2453, nMn55 = 2555, nFe56 = 2656, nNi58 = 2858, nCu63 = 2963, nCu65 = 2965, nGd155 = 64155, nGd157 = 64157, nPb208 = 82208, nPb206 = 82206, nPb207 = 82207;

        //scattering incoherent
        //Float_t siAlu = 0.0082, siCu = 0.55, siKapton = 804, siBoron = 3, siLiF = 1798.1, siPoly = 1766.6, siEStahl = 1.2, siArCO2 = 0.176, siLuft = 0.0022, siWater = 160.6, siZirkon = 0.2 , siBronze = 0.51;
        //Float_t siAlu = 1.503, siCu = 8.03, siKapton = 968.1, siBoron = 3.1, siLiF = 4.84, siPoly = 1775.6, siEStahl = 10.85, siArCO2 = 2.68, siLuft = 20.1, siWater = 168.12, siZirkon = 6.46, siBronze = 0.51, siBoronNat = 5.24;
        //Float_t siQuarz = 10.63, siAl2O3 = 15.7, siFe = 0, siHDPE = 0, siHe3 = 0.;
        //Float_t siHydrogen = 82.02, siOxygen = 4.232, siNitrogen = 11.51, siSilicon = 2.167;

        //densities rho
        //const double rAlu = 2.7, rCu = 8.94, rKapton = 1.43, rBoron = 2.34, rLiF = 2.6, rPoly = 1.14, rArCO2 = 0.0018, rCr = 7.14, rZirkon = 6.5, rBronze = 8.9, rSalt = 2.16, rBoronNat = 2.46, rB4CNat = 2.51, rB4C = 2.42, rTi = 4.5, rKalium = 0.856 ;
        const double  rQuarz = 2.5, rGraphit = 2.2, rAl2O3 = 3.94, rFe = 7.87, rAlMg = 2.66, rEStahl = 8.03, rPb = 11.342, rTNT = 1.654;
        //Float_t rAsphalt = 1.0;
        //double rHDPE2;
        double rWater = 0.99, rGd2O3 = 7.41, rMethane = 0.000656, rHDPE = 0.95, rPVC = 1.45, rPVCc = 1.56 ; //rLuft = pressureFac*0.0012, rLuftWater = 1.*2.*0.0000177*relHumidityAir;
        double rHe3 = 0.000125, rBF3 = 0.00276, rDiesel = 0.83, rSaltWater = 1., rCellulose = 1.5, rKieferTr = 0.44, rBucheTr = 0.59, rGeneral = 1.;


        rLuft = atmPressure / 287.058 / sysTemperature / 1e3;
        rLuftWater = absHumidityAir / 1e6;
        //rLuftWater = getRLuftWasser(sysTemperature);
        //rLuft = getRLuft(sysTemperature,atmPressure);

        // ps = exp(-6094.4642/T + 21.1249952 - 0.027245552*T + 1.6853396*10-5*T^2 + 2.4575506*ln(T)) with T in K and ps in Pa
        // for 10^5 Pa and 290 K: Sättigungsdruck ps: 1936.3 Pa is rLuftWater= relHum * ps / (RD + T) = 0.0000144678 for RD = 461.5 Nm/(kgK)

        //float deadMaterialFactor = 0.565;
        //float deadMaterialFactor = 0.355; //for 1.5 mum (23% absorbed in boron, 10 % cut away from the remaining, plus dead material)
        float deadMaterialFactor = 0.31; //for 1.2 mum (18% absorbed in boron, 8 % cut away from the remaining, plus dead material)
        //float deadMaterialFactor = 0.288; //for 1.0 mum (15% absorbed in boron, 7 % cut away from the remaining, plus dead material)
        //float deadMaterialFactor = 0.282; //for 0.8 mum (12.2% absorbed in boron, 7 % cut away from the remaining, plus dead material)

        rHe3 = 1. * rHe3 * 1.5 * pi / 4.; //rHe3 = rHe3*1.72;
        rBF3 = rBF3 * 0.7;//(0.342+0.02+0.0109);
        //rBF3 = rBF3*1.*(0.427+0.22+0.02+0.0109); //0.9*pi/4./5.; 0.5*rBF3*pi/4.;// //0.55;// * 0.9*5.; 1 mum 2layer B4C BF3 equals 0.91 cm BF3 or maybe rather 0.37 bar BF3 equals a cylindrical tube of 50mm (but made recatangular) with 1.3 mum B4C coating
        //rGd2O3 = rGd2O3/1000.;
        rHDPE = rHDPE * 1.00;
        rMethane = 8. * rMethane;

        Float_t rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFrac * rWater;
        Float_t rSolids = soilSiFrac * rQuarz + soilAlFrac * rAl2O3;

        Float_t rCelluloseMix = rCelluloseFrac*rCellulose+rCelluloseWaterFrac*rWater;

        Float_t rMaterial = 1.0;
        if (drawSingleNeutronPropagation) rBoden = rBoden / 5.;
        if (drawSingleNeutronPropagation) rWater = rWater / 5.;


        //atomic weight
        //Float_t wAlu = 27., wCu = 63.5, wKapton = 382., wBoron = 10.013, wLiF = 86.42, wPoly = 226., wEStahl = 55.7, wArCO2 = 34.9, wLuft = 28., wWater = 18., wZirkon = 91.2, wBronze = 67.5;
        const Float_t wAlu = 26.98,  wB10 = 10., wB11 = 11., wBoron = 10.013, wEStahl = 55.7,  wLuft = 28.965;
        const Float_t wWater = 18.015;
        //const Float_t wCu = 63.5, wKapton = 382., wZirkon = 91.2, wBronze = 67.5, wPoly = 226.32, wArCO2 = 40.6, wLiF = 25.94, wBoronNat = 10.8, wB4CNat = 13.8, wB4C = 13.;
        const Float_t wQuarz = 60.1, wAl2O3 = 101.96, wNi = 58.7, wNa = 23., wCl = 35.45, wFe = 55.85, wCr52 = 52., wCr53 = 53., wMn55 = 55., wAr = 39.95, wSulfur = 32., wHe3 = 3., wBF3 = 67.82;
        const Float_t wGd2O3 = 362.49, wMethane = 16.04, wSalt = 58.5, wGraphit = 12., wDiesel = 167., wPb = 207.2, wPb206 = 206., wPb207 = 207., wPb208 = 208., wTi = 47.867, wTi48 = 48., wKalium = 39.1, wTNT = 11.03, wPVC = 62.5, wPVCc = 367.15, wCellulose = 162.14;

        //rCelluloseFrac = 0.333;
        //rCelluloseWaterFrac = 0*0.25;
        //other atomic weights
        const Float_t wHydrogen = 1., wCarbon = 12., wOxygen = 16., wNitrogen = 14., wSilicon = 28.09, wAl = 26.98, wF = 19., wGd155 = 155., wGd157 = 157.;
        Float_t wSaltWater = 18.;
        Float_t wBoden = soilSolidFracVar * (soilSiFrac * wQuarz + soilAlFrac * wAl2O3) + soilWaterFrac * wWater * soilStrechFactor;
        Float_t wSolids = soilSiFrac * wQuarz + soilAlFrac * wAl2O3;
        Float_t wPlants = wWater + 0.2 * wCarbon;
        // http://www.tll.de/ainfo/pdf/biog1205.pdf
        Float_t wAsphalt = 0.11 * wCarbon + 0.14 * wHydrogen + 0.25 * wSilicon + 0.5 * wOxygen;
        Float_t wCatLitter = 0.44 * wHydrogen + 0.12 * wSilicon + 0.44 * wOxygen;
        const Float_t wHDPE = 14.025;
        double nCellulose = rCelluloseFrac*rCellulose/wCellulose;
        double nCelluloseWater = rCelluloseWaterFrac*rWater/wWater;
        double nCelluloseMix = nCellulose + nCelluloseWater;
        Float_t wCelluloseMix = rCelluloseMix /(nCelluloseMix); if (wCelluloseMix == 0) wCelluloseMix = 1.;


        double asHLast = 0, asOLast = 0, asNLast = 0, asArLast = 0, asSiLast = 0, asAlLast = 0, asB10Last = 0, asCLast = 0;
        double csHLast = 0, csOLast = 0, csNLast = 0, csArLast = 0, csSiLast = 0, csAlLast = 0, csB10Last = 0, csCLast = 0;
        double lastEnergy = 0, lastEnergy11 = 0, lastEnergy19 = 0, lastEnergy20 = 0, lastEnergy25 = 0;

        float nSalt = 1.;
        float saltConcentration = 3.5; // %

        rSaltWater = 0.99 + saltConcentration / 100. * 0.8;
        nSalt = (saltConcentration * wWater) / ((100. - saltConcentration) * wSalt);
        wSaltWater = wWater + nSalt * wSalt;

        //Carbon 11 %    Hydrogen 14 %    Silicon 25 %     Oxygen 50%

        //Geometries
        //scheme: (x1,y1,x2,y2 (all (mm)),z (mm),thickness (mm),type,layer)
        //(types: B=1, Al=2, Cu=3, Kapton=4, LiF = 5, Polyamid = 6, Zirkon = 7, natural Boron = 8, water = 9, air (dry) = 10, air (wet) = 11, Quarz = 12)

        // ****************************************
        // Reading all cross sections from matrices

        if (!noGUIMode) {uiM->setStatus(1, "Reading Angular Tabulated Cross Sections");        delay(5);}
        else {cout << "Reading Cross Sections" <<endl;}

        vector<TMatrixF> angularHighEnergySi = readAngularTabulatedCoefficients(endfFolder, "angularSi28tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyN = readAngularTabulatedCoefficients(endfFolder, "angularN14tabulated2.txt", 182);
        vector<TMatrixF> angularHighEnergyO = readAngularTabulatedCoefficients(endfFolder, "angularO16tabulated2.txt", 182);
        vector<TMatrixF> angularHighEnergyAl = readAngularTabulatedCoefficients(endfFolder, "angularAl27tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyC = readAngularTabulatedCoefficients(endfFolder, "angularC12tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyAr = readAngularTabulatedCoefficients(endfFolder, "angularAr40tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyF = readAngularTabulatedCoefficients(endfFolder, "angularF19tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyFe = readAngularTabulatedCoefficients(endfFolder, "angularFe56tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyNa = readAngularTabulatedCoefficients(endfFolder, "angularNa23tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyCl = readAngularTabulatedCoefficients(endfFolder, "angularCl35tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyMn = readAngularTabulatedCoefficients(endfFolder, "angularMn55tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyPb206 = readAngularTabulatedCoefficients(endfFolder, "angularPb206tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyPb207 = readAngularTabulatedCoefficients(endfFolder, "angularPb207tabulated.txt", 182);
        vector<TMatrixF> angularHighEnergyPb208 = readAngularTabulatedCoefficients(endfFolder, "angularPb208tabulated.txt", 182);

        if (!noGUIMode) {uiM->setStatus(1, "Reading Angular Cross Sections");        delay(5);}

        static TMatrixF angularSi = readAngularCoefficients(endfFolder, "angularSi28.txt");
        static TMatrixF angularN = readAngularCoefficients(endfFolder, "angularN14.txt");
        static TMatrixF angularO = readAngularCoefficients(endfFolder, "angularO16.txt");
        static TMatrixF angularH = readAngularCoefficients(endfFolder, "angularH1.txt");
        static TMatrixF angularAl = readAngularCoefficients(endfFolder, "angularAl27.txt");
        static TMatrixF angularC = readAngularCoefficients(endfFolder, "angularC12.txt");
        static TMatrixF angularAr = readAngularCoefficients(endfFolder, "angularAr40.txt");
        static TMatrixF angularFe = readAngularCoefficients(endfFolder, "angularFe56.txt");
        static TMatrixF angularS = readAngularCoefficients(endfFolder, "angularS32.txt");
        static TMatrixF angularHe3 = readAngularCoefficients(endfFolder, "angularHe3.txt");
        static TMatrixF angularB10 = readAngularCoefficients(endfFolder, "angularB10.txt");
        static TMatrixF angularB11 = readAngularCoefficients(endfFolder, "angularB11.txt");
        static TMatrixF angularGd155 = readAngularCoefficients(endfFolder, "angularGd155.txt");
        static TMatrixF angularGd157 = readAngularCoefficients(endfFolder, "angularGd157.txt");
        static TMatrixF angularF = readAngularCoefficients(endfFolder, "angularF19.txt");
        static TMatrixF angularNa = readAngularCoefficients(endfFolder, "angularNa23.txt");
        static TMatrixF angularCl35 = readAngularCoefficients(endfFolder, "angularCl35.txt");
        static TMatrixF angularCr52 = readAngularCoefficients(endfFolder, "angularCr52.txt");
        static TMatrixF angularCr53 = readAngularCoefficients(endfFolder, "angularCr53.txt");
        static TMatrixF angularNi58 = readAngularCoefficients(endfFolder, "angularNi58.txt");
        static TMatrixF angularMn55 = readAngularCoefficients(endfFolder, "angularMn55.txt");
        static TMatrixF angularPb206 = readAngularCoefficients(endfFolder, "angularPb206.txt");
        static TMatrixF angularPb207 = readAngularCoefficients(endfFolder, "angularPb207.txt");
        static TMatrixF angularPb208 = readAngularCoefficients(endfFolder, "angularPb208.txt");
        //TMatrixF angularK39,angularTi48;
        static TMatrixF angularK39 = readAngularCoefficients(endfFolder, "angularK39.txt");
        static TMatrixF angularTi48 = readAngularCoefficients(endfFolder, "angularTi48.txt");

        if (!noGUIMode) {uiM->setStatus(1, "Reading Elastic Cross Sections");        delay(5);}

        static TMatrixF sigmaSi = readSigmaEnergy(endfFolder, "elasticSi28.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaSi, sigmaElasticFactor*0.01,sigmaElasticFactor*0.08);
        static TMatrixF sigmaN = readSigmaEnergy(endfFolder, "elasticN14.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaN, sigmaElasticFactor*0.04,sigmaElasticFactor*0.02); //not known in fact, taken from N15
        static TMatrixF sigmaO = readSigmaEnergy(endfFolder, "elasticO16.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaO, sigmaElasticFactor*0.02,sigmaElasticFactor*0.04);
        static TMatrixF sigmaH = readSigmaEnergy(endfFolder, "elasticH1mod.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaH, sigmaElasticFactor*0.003,sigmaElasticFactor*0.024);
        static TMatrixF sigmaAl = readSigmaEnergy(endfFolder, "elasticAl27.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaAl, sigmaElasticFactor*0.018,sigmaElasticFactor*0.13);
        static TMatrixF sigmaC = readSigmaEnergy(endfFolder, "elasticC12.txt");		//if (varyCrosssections) modifyCSmatrix(&sigmaC, sigmaElasticFactor*0.005,sigmaElasticFactor*0.10);
        static TMatrixF sigmaAr = readSigmaEnergy(endfFolder, "elasticAr40.txt");
        static TMatrixF sigmaFe = readSigmaEnergy(endfFolder, "elasticFe56.txt");
        static TMatrixF sigmaS = readSigmaEnergy(endfFolder, "elasticS32.txt");
        static TMatrixF sigmaHe3 = readSigmaEnergy(endfFolder, "elasticHe3.txt");
        static TMatrixF sigmaB10 = readSigmaEnergy(endfFolder, "elasticB10.txt");
        static TMatrixF sigmaB11 = readSigmaEnergy(endfFolder, "elasticB11.txt");
        static TMatrixF sigmaF = readSigmaEnergy(endfFolder, "elasticF19.txt");
        static TMatrixF sigmaGd155 = readSigmaEnergy(endfFolder, "elasticGd155.txt");
        static TMatrixF sigmaGd157 = readSigmaEnergy(endfFolder, "elasticGd157.txt");
        static TMatrixF sigmaNa = readSigmaEnergy(endfFolder, "elasticNa23.txt");
        static TMatrixF sigmaCl35 = readSigmaEnergy(endfFolder, "elasticCl35.txt");
        static TMatrixF sigmaCr52 = readSigmaEnergy(endfFolder, "elasticCr52.txt");
        static TMatrixF sigmaCr53 = readSigmaEnergy(endfFolder, "elasticCr53.txt");
        static TMatrixF sigmaNi58 = readSigmaEnergy(endfFolder, "elasticNi58.txt");
        static TMatrixF sigmaMn55 = readSigmaEnergy(endfFolder, "elasticMn55.txt");
        static TMatrixF sigmaPb206 = readSigmaEnergy(endfFolder, "elasticPb206.txt");
        static TMatrixF sigmaPb207 = readSigmaEnergy(endfFolder, "elasticPb207.txt");
        static TMatrixF sigmaPb208 = readSigmaEnergy(endfFolder, "elasticPb208.txt");
        static TMatrixF sigmaK39 = readSigmaEnergy(endfFolder, "elasticK39.txt");
        static TMatrixF sigmaTi48 = readSigmaEnergy(endfFolder, "elasticTi48.txt");


        TSpline3* sigmaSiSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaSi, true, false);
        TSpline3* sigmaNSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaN, true, false);
        TSpline3* sigmaOSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaO, true, false);
        TSpline3* sigmaHSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaH, true, false);
        TSpline3* sigmaAlSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaAl, true, false);
        TSpline3* sigmaCSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaC, true, false);
        TSpline3* sigmaArSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaAr, true, false);
        TSpline3* sigmaFeSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaFe, true, false);
        TSpline3* sigmaHe3Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaHe3, true, false);
        TSpline3* sigmaSSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaS, true, false);
        TSpline3* sigmaB10Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaB10, true, false);
        TSpline3* sigmaB11Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaB11, true, false);
        TSpline3* sigmaFSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaF, true, false);
        TSpline3* sigmaGd155Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaGd155, true, false);
        TSpline3* sigmaGd157Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaGd157, true, false);
        TSpline3* sigmaNaSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaNa, true, false);
        TSpline3* sigmaClSpline = uiM->getSplinedEnergyModelFromMatrix(&sigmaCl35, true, false);
        TSpline3* sigmaCr52Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaCr52, true, false);
        TSpline3* sigmaCr53Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaCr53, true, false);
        TSpline3* sigmaNi58Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaNi58, true, false);
        TSpline3* sigmaMn55Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaMn55, true, false);
        TSpline3* sigmaPb206Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaPb206, true, false);
        TSpline3* sigmaPb207Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaPb207, true, false);
        TSpline3* sigmaPb208Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaPb208, true, false);
        TSpline3* sigmaK39Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaK39, true, false);
        TSpline3* sigmaTi48Spline = uiM->getSplinedEnergyModelFromMatrix(&sigmaTi48, true, false);

        csMatrixHP = &sigmaH;
        csMatrixOP = &sigmaO;
        csMatrixNP = &sigmaN;

        if (!noGUIMode) {uiM->setStatus(1, "Reading Absorption Cross Sections");        delay(5);}

        static TMatrixF absorbSi = readSigmaEnergy(endfFolder, "absorbSi28.txt");           //if (varyCrosssections) modifyCSmatrix(&absorbSi, sigmaAbsorbFactor*0.20,sigmaAbsorbFactor*0.40);
        static TMatrixF absorbN = readSigmaEnergy(endfFolder, "absorbN14.txt");             //if (varyCrosssections) modifyCSmatrix(&absorbN, sigmaAbsorbFactor*0.30,sigmaAbsorbFactor*0.075);
        //that (b) is the (n,p) reaction n+N -> C+p MT103
        static TMatrixF absorbNb = readSigmaEnergy(endfFolder, "absorbN14b.txt");           //if (varyCrosssections) modifyCSmatrix(&absorbNb, sigmaAbsorbFactor*0.05,sigmaAbsorbFactor*0.10);
        static TMatrixF absorbO = readSigmaEnergy(endfFolder, "absorbO16.txt");             //if (varyCrosssections) modifyCSmatrix(&absorbO, sigmaAbsorbFactor*0.10,sigmaAbsorbFactor*0.10);
        static TMatrixF absorbH = readSigmaEnergy(endfFolder, "absorbH1.txt");              //if (varyCrosssections) modifyCSmatrix(&absorbH, sigmaAbsorbFactor*0.026,sigmaAbsorbFactor*0.25);
        static TMatrixF absorbAl = readSigmaEnergy(endfFolder, "absorbAl27.txt");           //if (varyCrosssections) modifyCSmatrix(&absorbAl, sigmaAbsorbFactor*0.017,sigmaAbsorbFactor*0.60);
        static TMatrixF absorbC = readSigmaEnergy(endfFolder, "absorbC12.txt");             //if (varyCrosssections) modifyCSmatrix(&absorbC, sigmaAbsorbFactor*0.03,sigmaAbsorbFactor*0.20);
        static TMatrixF absorbAr = readSigmaEnergy(endfFolder, "absorbAr40.txt");
        static TMatrixF absorbFe = readSigmaEnergy(endfFolder, "absorbFe56.txt");
        static TMatrixF absorbS = readSigmaEnergy(endfFolder, "absorbS32.txt");
        static TMatrixF absorbHe3 = readSigmaEnergy(endfFolder, "absorbHe3.txt"); //z gamma
        static TMatrixF absorbB10b = readSigmaEnergy(endfFolder, "absorbB10.txt");
        static TMatrixF absorbB10 = readSigmaEnergy(endfFolder, "absorbMt107B10.txt"); //n alpha
        static TMatrixF absorbB11 = readSigmaEnergy(endfFolder, "absorbMt107B11.txt"); //n alpha
        static TMatrixF absorbF = readSigmaEnergy(endfFolder, "absorbF19.txt"); //z gamma
        static TMatrixF absorbGd155 = readSigmaEnergy(endfFolder, "absorbGd155.txt");
        static TMatrixF absorbGd157 = readSigmaEnergy(endfFolder, "absorbGd157.txt");
        static TMatrixF absorbNa = readSigmaEnergy(endfFolder, "absorbNa23.txt");
        static TMatrixF absorbCl35 = readSigmaEnergy(endfFolder, "absorbCl35.txt");
        static TMatrixF absorbCr52 = readSigmaEnergy(endfFolder, "absorbCr52.txt");
        static TMatrixF absorbCr53 = readSigmaEnergy(endfFolder, "absorbCr53.txt");
        static TMatrixF absorbNi58 = readSigmaEnergy(endfFolder, "absorbNi58.txt");
        static TMatrixF absorbMn55 = readSigmaEnergy(endfFolder, "absorbMn55.txt");
        static TMatrixF absorbPb206 = readSigmaEnergy(endfFolder, "absorbPb206.txt");
        static TMatrixF absorbPb207 = readSigmaEnergy(endfFolder, "absorbPb207.txt");
        static TMatrixF absorbPb208 = readSigmaEnergy(endfFolder, "absorbPb208.txt");
        static TMatrixF absorbK39 = readSigmaEnergy(endfFolder, "absorbK39.txt");
        static TMatrixF absorbTi48 = readSigmaEnergy(endfFolder, "absorbTi48.txt");

        TSpline3* absorbSiSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbSi, true, false);
        TSpline3* absorbNSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbN, true, false);
        TSpline3* absorbNbSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbNb, true, false);
        TSpline3* absorbOSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbO, true, false);
        TSpline3* absorbHSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbH, true, false);
        TSpline3* absorbAlSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbAl, true, false);
        TSpline3* absorbCSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbC, true, false);
        TSpline3* absorbArSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbAr, true, false);
        TSpline3* absorbFeSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbFe, true, false);
        TSpline3* absorbSSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbS, true, false);
        TSpline3* absorbHe3Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbHe3, true, false);
        TSpline3* absorbB10bSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbB10b, true, false);
        TSpline3* absorbB10Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbB10, true, false);
        TSpline3* absorbB11Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbB11, true, false);
        TSpline3* absorbFSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbF, true, false);
        TSpline3* absorbGd155Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbGd155, true, false);
        TSpline3* absorbGd157Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbGd157, true, false);
        TSpline3* absorbNaSpline = uiM->getSplinedEnergyModelFromMatrix(&absorbNa, true, false);
        TSpline3* absorbCl35Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbCl35, true, false);
        TSpline3* absorbCr52Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbCr52, true, false);
        TSpline3* absorbCr53Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbCr53, true, false);
        TSpline3* absorbNi58Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbNi58, true, false);
        TSpline3* absorbMn55Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbMn55, true, false);
        TSpline3* absorbPb206Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbPb206, true, false);
        TSpline3* absorbPb207Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbPb207, true, false);
        TSpline3* absorbPb208Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbPb208, true, false);
        TSpline3* absorbK39Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbK39, true, false);
        TSpline3* absorbTi48Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbTi48, true, false);

        /*
        int ptrc = 0;
        TGraph* splineGr = new TGraph();
        TGraph* csGr  = new TGraph();
        TGraph* csDiffGr  = new TGraph();
        int imax = 2e8;
        double diffaR, aR, bR, energyEval, sumDiff = 0;

        for (float i = -3; i < 7; i+=0.005)
        {
            energyEval = TMath::Power(10,i*1.);

            aR = calcMeanCS(absorbNb, energyEval);
            bR = absorbNbSpline->Eval(i*1.);

            diffaR = aR-bR; sumDiff+= abs(diffaR); //if (diffaR > 0.5) diffaR = 0.5; if (diffaR < -0.5) diffaR = -0.5;
            csGr->SetPoint(ptrc,energyEval, aR);
            splineGr->SetPoint(ptrc,energyEval, bR);
            csDiffGr->SetPoint(ptrc,energyEval, diffaR);
            ptrc++;
        }
        cout<<"Sum: "<<diffaR*100.<<endl;

        TCanvas* cRatesPH2 = new TCanvas("cRatesPH2", "cRatesPH2", 3200, 1600);
        CanvasFashion(cRatesPH2);
        csDiffGr->Draw("APL");
        cRatesPH2->SetLogx();
        cRatesPH2->SaveAs(outputFolder + "/cssDallabsNb.png");

        TCanvas* cRatesPH = new TCanvas("cRatesPH", "cRatesPH", 3200, 1600);
        CanvasFashion(cRatesPH);
        splineGr->Draw("APL");
        cRatesPH->SetLogx();
        TGraphFashion(splineGr, "energy", "value", true); splineGr->SetMarkerColor(kBlue); splineGr->SetLineColor(kBlue);
        csGr->Draw("SAMEPL");
        cRatesPH->SaveAs(outputFolder + "/cssallabsNb.png");
        */

        csMatrixHabsP = &absorbH;
        csMatrixOabsP = &absorbO;
        csMatrixNabsP = &absorbNb;

        static TMatrixF absorbMt209H = readSigmaEnergy(endfFolder, "absorbMt209H1.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt209H, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5H = readSigmaEnergy(endfFolder, "absorbMt5H1.txt");		//if (varyCrosssections) modifyCSmatrix(&absorbMt5H, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5N = readSigmaEnergy(endfFolder, "absorbMt5N14.txt");		//if (varyCrosssections) modifyCSmatrix(&absorbMt5N, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103N = readSigmaEnergy(endfFolder, "absorbMt103N14.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt103N, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt104N = readSigmaEnergy(endfFolder, "absorbMt104N14.txt");
        static TMatrixF absorbMt105N = readSigmaEnergy(endfFolder, "absorbMt105N14.txt");
        static TMatrixF absorbMt107N = readSigmaEnergy(endfFolder, "absorbMt107N14.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt107N, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt108N = readSigmaEnergy(endfFolder, "absorbMt108N14.txt");
        static TMatrixF absorbMt209N = readSigmaEnergy(endfFolder, "absorbMt209N14.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt209N, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5C = readSigmaEnergy(endfFolder, "absorbMt5C12.txt");		//if (varyCrosssections) modifyCSmatrix(&absorbMt5C, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103C = readSigmaEnergy(endfFolder, "absorbMt103C12.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt103C, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt107C = readSigmaEnergy(endfFolder, "absorbMt107C12.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt107C, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt209C = readSigmaEnergy(endfFolder, "absorbMt209C12.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt209C, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5O = readSigmaEnergy(endfFolder, "absorbMt5O16.txt");		//if (varyCrosssections) modifyCSmatrix(&absorbMt5O, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103O = readSigmaEnergy(endfFolder, "absorbMt103O16.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt103O, sigmaAbsorbFactor*0.50,sigmaAbsorbFactor*0.50);
        static TMatrixF absorbMt105O = readSigmaEnergy(endfFolder, "absorbMt105O16.txt");
        static TMatrixF absorbMt107O = readSigmaEnergy(endfFolder, "absorbMt107O16.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt107O, sigmaAbsorbFactor*0.10,sigmaAbsorbFactor*0.10);
        static TMatrixF absorbMt209O = readSigmaEnergy(endfFolder, "absorbMt209O16.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt209O, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5Al = readSigmaEnergy(endfFolder, "absorbMt5Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt5Al, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103Al = readSigmaEnergy(endfFolder, "absorbMt103Al27.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt103Al, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt107Al = readSigmaEnergy(endfFolder, "absorbMt107Al27.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt5Al, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt209Al = readSigmaEnergy(endfFolder, "absorbMt209Al27.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt107Al, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5Si = readSigmaEnergy(endfFolder, "absorbMt5Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt209Al, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103Si = readSigmaEnergy(endfFolder, "absorbMt103Si28.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt5Si, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt107Si = readSigmaEnergy(endfFolder, "absorbMt107Si28.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt107Si, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt209Si = readSigmaEnergy(endfFolder, "absorbMt209Si28.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt209Si, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5Ar = readSigmaEnergy(endfFolder, "absorbMt5Ar40.txt");	//if (varyCrosssections) modifyCSmatrix(&absorbMt5Ar, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt103Ar = readSigmaEnergy(endfFolder, "absorbMt103Ar40.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt103Ar, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt107Ar = readSigmaEnergy(endfFolder, "absorbMt107Ar40.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt107Ar, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt209Ar = readSigmaEnergy(endfFolder, "absorbMt209Ar40.txt");//if (varyCrosssections) modifyCSmatrix(&absorbMt209Ar, sigmaHeFactor*sigmaHeEstimator,sigmaHeFactor*sigmaHeEstimator);
        static TMatrixF absorbMt5Fe = readSigmaEnergy(endfFolder, "absorbMt5Fe56.txt");
        static TMatrixF absorbMt103Fe = readSigmaEnergy(endfFolder, "absorbMt103Fe56.txt");
        static TMatrixF absorbMt107Fe = readSigmaEnergy(endfFolder, "absorbMt107Fe56.txt");
        static TMatrixF absorbMt209Fe = readSigmaEnergy(endfFolder, "absorbMt209Fe56.txt");
        static TMatrixF absorbMt5S = readSigmaEnergy(endfFolder, "absorbMt5S32.txt");
        static TMatrixF absorbMt103S = readSigmaEnergy(endfFolder, "absorbMt103S32.txt");
        static TMatrixF absorbMt107S = readSigmaEnergy(endfFolder, "absorbMt107S32.txt");
        static TMatrixF absorbMt103He3 = readSigmaEnergy(endfFolder, "absorbMt103He3.txt");
        static TMatrixF absorbMt104He3 = readSigmaEnergy(endfFolder, "absorbMt104He3.txt");
        static TMatrixF absorbMt103F = readSigmaEnergy(endfFolder, "absorbMt103F19.txt");
        static TMatrixF absorbMt107F = readSigmaEnergy(endfFolder, "absorbMt107F19.txt");
        static TMatrixF absorbMt5Na = readSigmaEnergy(endfFolder, "absorbMt5Na23.txt");
        static TMatrixF absorbMt103Na = readSigmaEnergy(endfFolder, "absorbMt103Na23.txt");
        static TMatrixF absorbMt107Na = readSigmaEnergy(endfFolder, "absorbMt107Na23.txt");
        static TMatrixF absorbMt5Cl35 = readSigmaEnergy(endfFolder, "absorbMt5Cl35.txt");
        static TMatrixF absorbMt103Cl35 = readSigmaEnergy(endfFolder, "absorbMt103Cl35.txt");
        static TMatrixF absorbMt107Cl35 = readSigmaEnergy(endfFolder, "absorbMt107Cl35.txt");
        static TMatrixF absorbMt5Cr52 = readSigmaEnergy(endfFolder, "absorbMt5Cr52.txt");
        static TMatrixF absorbMt103Cr52 = readSigmaEnergy(endfFolder, "absorbMt103Cr52.txt");
        static TMatrixF absorbMt107Cr52 = readSigmaEnergy(endfFolder, "absorbMt107Cr52.txt");
        static TMatrixF absorbMt5Cr53 = readSigmaEnergy(endfFolder, "absorbMt5Cr53.txt");
        static TMatrixF absorbMt103Cr53 = readSigmaEnergy(endfFolder, "absorbMt103Cr53.txt");
        static TMatrixF absorbMt107Cr53 = readSigmaEnergy(endfFolder, "absorbMt107Cr53.txt");
        static TMatrixF absorbMt5Ni58 = readSigmaEnergy(endfFolder, "absorbMt5Ni58.txt");
        static TMatrixF absorbMt103Ni58 = readSigmaEnergy(endfFolder, "absorbMt103Ni58.txt");
        static TMatrixF absorbMt107Ni58 = readSigmaEnergy(endfFolder, "absorbMt107Ni58.txt");
        static TMatrixF absorbMt5Mn55 = readSigmaEnergy(endfFolder, "absorbMt5Mn55.txt");
        static TMatrixF absorbMt103Mn55 = readSigmaEnergy(endfFolder, "absorbMt103Mn55.txt");
        static TMatrixF absorbMt107Mn55 = readSigmaEnergy(endfFolder, "absorbMt107Mn55.txt");
        static TMatrixF absorbMt5Pb206 = readSigmaEnergy(endfFolder, "absorbMt5Pb206.txt");
        static TMatrixF absorbMt103Pb206 = readSigmaEnergy(endfFolder, "absorbMt103Pb206.txt");
        static TMatrixF absorbMt107Pb206 = readSigmaEnergy(endfFolder, "absorbMt107Pb206.txt");
        static TMatrixF absorbMt209Pb206 = readSigmaEnergy(endfFolder, "absorbMt209Pb206.txt");
        static TMatrixF absorbMt5Pb207 = readSigmaEnergy(endfFolder, "absorbMt5Pb207.txt");
        static TMatrixF absorbMt103Pb207 = readSigmaEnergy(endfFolder, "absorbMt103Pb207.txt");
        static TMatrixF absorbMt107Pb207 = readSigmaEnergy(endfFolder, "absorbMt107Pb207.txt");
        static TMatrixF absorbMt209Pb207 = readSigmaEnergy(endfFolder, "absorbMt209Pb207.txt");
        static TMatrixF absorbMt5Pb208 = readSigmaEnergy(endfFolder, "absorbMt5Pb208.txt");
        static TMatrixF absorbMt103Pb208 = readSigmaEnergy(endfFolder, "absorbMt103Pb208.txt");
        static TMatrixF absorbMt107Pb208 = readSigmaEnergy(endfFolder, "absorbMt107Pb208.txt");
        static TMatrixF absorbMt209Pb208 = readSigmaEnergy(endfFolder, "absorbMt209Pb208.txt");
        static TMatrixF absorbMt103K39 = readSigmaEnergy(endfFolder, "absorbMt103K39.txt");
        static TMatrixF absorbMt107K39 = readSigmaEnergy(endfFolder, "absorbMt107K39.txt");
        static TMatrixF absorbMt5Ti48 = readSigmaEnergy(endfFolder, "absorbMt5Ti48.txt");
        static TMatrixF absorbMt103Ti48 = readSigmaEnergy(endfFolder, "absorbMt103Ti48.txt");
        static TMatrixF absorbMt104Ti48 = readSigmaEnergy(endfFolder, "absorbMt104Ti48.txt");
        static TMatrixF absorbMt107Ti48 = readSigmaEnergy(endfFolder, "absorbMt107Ti48.txt");
        static TMatrixF absorbMt207Ti48 = readSigmaEnergy(endfFolder, "absorbMt207Ti48.txt");

        TSpline3* absorbMt103Cl35Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbMt103Cl35, true, false);
        TSpline3* absorbMt107Cl35Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbMt107Cl35, true, false);
        TSpline3* absorbMt103He3Spline = uiM->getSplinedEnergyModelFromMatrix(&absorbMt103He3, true, false);

        if (!noGUIMode) {uiM->setStatus(1, "Reading Inelastic Cross Sections");        delay(5);}

        static TMatrixF sigmaInMt51N = readSigmaEnergy(endfFolder, "inelasticMt51N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt51N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt52N = readSigmaEnergy(endfFolder, "inelasticMt52N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt52N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt53N = readSigmaEnergy(endfFolder, "inelasticMt53N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt53N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt54N = readSigmaEnergy(endfFolder, "inelasticMt54N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt54N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt55N = readSigmaEnergy(endfFolder, "inelasticMt55N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt55N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt56N = readSigmaEnergy(endfFolder, "inelasticMt56N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt56N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt57N = readSigmaEnergy(endfFolder, "inelasticMt57N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt57N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt58N = readSigmaEnergy(endfFolder, "inelasticMt58N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt58N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt59N = readSigmaEnergy(endfFolder, "inelasticMt59N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt59N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt60N = readSigmaEnergy(endfFolder, "inelasticMt60N14.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt60N, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        vector<TMatrixF*> sigmaInNVec;
        sigmaInNVec.push_back(&sigmaInMt51N); sigmaInNVec.push_back(&sigmaInMt52N); sigmaInNVec.push_back(&sigmaInMt53N); sigmaInNVec.push_back(&sigmaInMt54N); sigmaInNVec.push_back(&sigmaInMt55N);
        sigmaInNVec.push_back(&sigmaInMt56N); sigmaInNVec.push_back(&sigmaInMt57N); sigmaInNVec.push_back(&sigmaInMt58N); sigmaInNVec.push_back(&sigmaInMt59N); sigmaInNVec.push_back(&sigmaInMt60N);
        vector<double> inelasticEnergyLossN;
        inelasticEnergyLossN.push_back(2.31); inelasticEnergyLossN.push_back(3.95); inelasticEnergyLossN.push_back(4.91); inelasticEnergyLossN.push_back(5.11); inelasticEnergyLossN.push_back(5.69);
        inelasticEnergyLossN.push_back(5.83); inelasticEnergyLossN.push_back(6.2);  inelasticEnergyLossN.push_back(6.446);  inelasticEnergyLossN.push_back(7.029);  inelasticEnergyLossN.push_back(7.9669);

        static TMatrixF sigmaInMt51O = readSigmaEnergy(endfFolder, "inelasticMt51O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt51O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt52O = readSigmaEnergy(endfFolder, "inelasticMt52O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt52O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt53O = readSigmaEnergy(endfFolder, "inelasticMt53O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt53O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt54O = readSigmaEnergy(endfFolder, "inelasticMt54O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt54O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt55O = readSigmaEnergy(endfFolder, "inelasticMt55O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt55O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt56O = readSigmaEnergy(endfFolder, "inelasticMt56O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt56O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt57O = readSigmaEnergy(endfFolder, "inelasticMt57O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt57O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt58O = readSigmaEnergy(endfFolder, "inelasticMt58O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt58O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt59O = readSigmaEnergy(endfFolder, "inelasticMt59O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt59O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt60O = readSigmaEnergy(endfFolder, "inelasticMt60O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt60O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt61O = readSigmaEnergy(endfFolder, "inelasticMt61O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt61O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt62O = readSigmaEnergy(endfFolder, "inelasticMt62O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt62O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt63O = readSigmaEnergy(endfFolder, "inelasticMt63O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt63O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt64O = readSigmaEnergy(endfFolder, "inelasticMt64O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt64O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt65O = readSigmaEnergy(endfFolder, "inelasticMt65O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt65O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt66O = readSigmaEnergy(endfFolder, "inelasticMt66O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt66O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt67O = readSigmaEnergy(endfFolder, "inelasticMt67O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt67O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt68O = readSigmaEnergy(endfFolder, "inelasticMt68O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt68O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt69O = readSigmaEnergy(endfFolder, "inelasticMt69O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt69O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt70O = readSigmaEnergy(endfFolder, "inelasticMt70O16.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt70O, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        vector<TMatrixF*> sigmaInOVec;
        sigmaInOVec.push_back(&sigmaInMt51O); sigmaInOVec.push_back(&sigmaInMt52O); sigmaInOVec.push_back(&sigmaInMt53O); sigmaInOVec.push_back(&sigmaInMt54O); sigmaInOVec.push_back(&sigmaInMt55O);
        sigmaInOVec.push_back(&sigmaInMt56O); sigmaInOVec.push_back(&sigmaInMt57O); sigmaInOVec.push_back(&sigmaInMt58O); sigmaInOVec.push_back(&sigmaInMt59O); sigmaInOVec.push_back(&sigmaInMt60O);
        sigmaInOVec.push_back(&sigmaInMt61O); sigmaInOVec.push_back(&sigmaInMt62O); sigmaInOVec.push_back(&sigmaInMt63O); sigmaInOVec.push_back(&sigmaInMt64O); sigmaInOVec.push_back(&sigmaInMt65O);
        sigmaInOVec.push_back(&sigmaInMt66O); sigmaInOVec.push_back(&sigmaInMt67O); sigmaInOVec.push_back(&sigmaInMt68O); sigmaInOVec.push_back(&sigmaInMt69O); sigmaInOVec.push_back(&sigmaInMt70O);
        vector<double> inelasticEnergyLossO;
        inelasticEnergyLossO.push_back(6.05); inelasticEnergyLossO.push_back(6.13); inelasticEnergyLossO.push_back(6.91); inelasticEnergyLossO.push_back(7.12); inelasticEnergyLossO.push_back(8.87);
        inelasticEnergyLossO.push_back(9.63);  inelasticEnergyLossO.push_back(9.85);  inelasticEnergyLossO.push_back(10.36);  inelasticEnergyLossO.push_back(10.96);  inelasticEnergyLossO.push_back(11.08);
        inelasticEnergyLossO.push_back(11.1);  inelasticEnergyLossO.push_back(11.52);  inelasticEnergyLossO.push_back(11.6);  inelasticEnergyLossO.push_back(12.05);  inelasticEnergyLossO.push_back(12.44);
        inelasticEnergyLossO.push_back(12.53);  inelasticEnergyLossO.push_back(12.8);   inelasticEnergyLossO.push_back(12.97);  inelasticEnergyLossO.push_back(13.02);  inelasticEnergyLossO.push_back(13.09);

        static TMatrixF sigmaInMt51Si = readSigmaEnergy(endfFolder, "inelasticMt51Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt51Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt52Si = readSigmaEnergy(endfFolder, "inelasticMt52Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt52Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt53Si = readSigmaEnergy(endfFolder, "inelasticMt53Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt53Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt54Si = readSigmaEnergy(endfFolder, "inelasticMt54Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt54Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt55Si = readSigmaEnergy(endfFolder, "inelasticMt55Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt55Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt56Si = readSigmaEnergy(endfFolder, "inelasticMt56Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt56Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt57Si = readSigmaEnergy(endfFolder, "inelasticMt57Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt57Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt58Si = readSigmaEnergy(endfFolder, "inelasticMt58Si28.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt58Si, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        vector<TMatrixF*> sigmaInSiVec;
        sigmaInSiVec.push_back(&sigmaInMt51Si); sigmaInSiVec.push_back(&sigmaInMt52Si); sigmaInSiVec.push_back(&sigmaInMt53Si); sigmaInSiVec.push_back(&sigmaInMt54Si); sigmaInSiVec.push_back(&sigmaInMt55Si);
        sigmaInSiVec.push_back(&sigmaInMt56Si); sigmaInSiVec.push_back(&sigmaInMt57Si); sigmaInSiVec.push_back(&sigmaInMt58Si);
        vector<double> inelasticEnergyLossSi;
        inelasticEnergyLossSi.push_back(1.78); inelasticEnergyLossSi.push_back(4.62); inelasticEnergyLossSi.push_back(4.98); inelasticEnergyLossSi.push_back(6.28); inelasticEnergyLossSi.push_back(6.69);
        inelasticEnergyLossSi.push_back(6.88); inelasticEnergyLossSi.push_back(6.889);  inelasticEnergyLossSi.push_back(7.381);

        static TMatrixF sigmaInMt51C = readSigmaEnergy(endfFolder, "inelasticMt51C12.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt51C, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt52C = readSigmaEnergy(endfFolder, "inelasticMt52C12.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt52C, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt53C = readSigmaEnergy(endfFolder, "inelasticMt53C12.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt53C, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        //54-57 like 53
        static TMatrixF sigmaInMt58C = readSigmaEnergy(endfFolder, "inelasticMt58C12.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt58C, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        vector<TMatrixF*> sigmaInCVec;
        sigmaInCVec.push_back(&sigmaInMt51C); sigmaInCVec.push_back(&sigmaInMt52C); sigmaInCVec.push_back(&sigmaInMt53C);
        sigmaInCVec.push_back(&sigmaInMt53C); sigmaInCVec.push_back(&sigmaInMt53C); sigmaInCVec.push_back(&sigmaInMt53C); sigmaInCVec.push_back(&sigmaInMt53C);
        sigmaInCVec.push_back(&sigmaInMt58C);
        vector<double> inelasticEnergyLossC;
        inelasticEnergyLossC.push_back(4.43); inelasticEnergyLossC.push_back(7.65); inelasticEnergyLossC.push_back(9.64);
        inelasticEnergyLossC.push_back(10.8); inelasticEnergyLossC.push_back(11.8); inelasticEnergyLossC.push_back(12.7); inelasticEnergyLossC.push_back(13.35);
        inelasticEnergyLossC.push_back(14.08);

        static TMatrixF sigmaInMt51Al = readSigmaEnergy(endfFolder, "inelasticMt51Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt51Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt52Al = readSigmaEnergy(endfFolder, "inelasticMt52Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt52Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt53Al = readSigmaEnergy(endfFolder, "inelasticMt53Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt53Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt54Al = readSigmaEnergy(endfFolder, "inelasticMt54Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt54Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt55Al = readSigmaEnergy(endfFolder, "inelasticMt55Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt55Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt56Al = readSigmaEnergy(endfFolder, "inelasticMt56Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt56Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt57Al = readSigmaEnergy(endfFolder, "inelasticMt57Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt57Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        static TMatrixF sigmaInMt58Al = readSigmaEnergy(endfFolder, "inelasticMt58Al27.txt");	//if (varyCrosssections) modifyCSmatrix(&sigmaInMt58Al, sigmaInelasticFactor*0.40,sigmaInelasticFactor*0.40);
        vector<TMatrixF*> sigmaInAlVec;
        sigmaInAlVec.push_back(&sigmaInMt51Al); sigmaInAlVec.push_back(&sigmaInMt52Al); sigmaInAlVec.push_back(&sigmaInMt53Al); sigmaInAlVec.push_back(&sigmaInMt54Al); sigmaInAlVec.push_back(&sigmaInMt55Al); sigmaInAlVec.push_back(&sigmaInMt56Al); sigmaInAlVec.push_back(&sigmaInMt57Al); sigmaInAlVec.push_back(&sigmaInMt58Al);
        vector<double> inelasticEnergyLossAl;
        inelasticEnergyLossAl.push_back(0.84); inelasticEnergyLossAl.push_back(1.01); inelasticEnergyLossAl.push_back(2.21); inelasticEnergyLossAl.push_back(2.73); inelasticEnergyLossAl.push_back(2.98); inelasticEnergyLossAl.push_back(3.0); inelasticEnergyLossAl.push_back(3.67); inelasticEnergyLossAl.push_back(3.96);

        static TMatrixF sigmaInMt51Fe = readSigmaEnergy(endfFolder, "inelasticMt51Fe56.txt");
        static TMatrixF sigmaInMt52Fe = readSigmaEnergy(endfFolder, "inelasticMt52Fe56.txt");
        static TMatrixF sigmaInMt53Fe = readSigmaEnergy(endfFolder, "inelasticMt53Fe56.txt");
        static TMatrixF sigmaInMt54Fe = readSigmaEnergy(endfFolder, "inelasticMt54Fe56.txt");
        static TMatrixF sigmaInMt55Fe = readSigmaEnergy(endfFolder, "inelasticMt55Fe56.txt");
        static TMatrixF sigmaInMt56Fe = readSigmaEnergy(endfFolder, "inelasticMt56Fe56.txt");
        static TMatrixF sigmaInMt57Fe = readSigmaEnergy(endfFolder, "inelasticMt57Fe56.txt");
        static TMatrixF sigmaInMt58Fe = readSigmaEnergy(endfFolder, "inelasticMt58Fe56.txt");
        vector<TMatrixF*> sigmaInFeVec;
        sigmaInFeVec.push_back(&sigmaInMt51Fe); sigmaInFeVec.push_back(&sigmaInMt52Fe); sigmaInFeVec.push_back(&sigmaInMt53Fe); sigmaInFeVec.push_back(&sigmaInMt54Fe); sigmaInFeVec.push_back(&sigmaInMt55Fe); sigmaInFeVec.push_back(&sigmaInMt56Fe); sigmaInFeVec.push_back(&sigmaInMt57Fe); sigmaInFeVec.push_back(&sigmaInMt58Fe);
        vector<double> inelasticEnergyLossFe;
        inelasticEnergyLossFe.push_back(0.847); inelasticEnergyLossFe.push_back(2.085); inelasticEnergyLossFe.push_back(2.658); inelasticEnergyLossFe.push_back(2.942); inelasticEnergyLossFe.push_back(2.96); inelasticEnergyLossFe.push_back(3.12); inelasticEnergyLossFe.push_back(3.123); inelasticEnergyLossFe.push_back(3.37);

        static TMatrixF sigmaInMt51S = readSigmaEnergy(endfFolder, "inelasticMt51Si28.txt");
        static TMatrixF sigmaInMt52S = readSigmaEnergy(endfFolder, "inelasticMt52Si28.txt");
        static TMatrixF sigmaInMt53S = readSigmaEnergy(endfFolder, "inelasticMt53Si28.txt");
        static TMatrixF sigmaInMt54S = readSigmaEnergy(endfFolder, "inelasticMt54Si28.txt");
        static TMatrixF sigmaInMt55S = readSigmaEnergy(endfFolder, "inelasticMt55Si28.txt");
        vector<TMatrixF*> sigmaInSVec;
        sigmaInSVec.push_back(&sigmaInMt51S); sigmaInSVec.push_back(&sigmaInMt52S); sigmaInSVec.push_back(&sigmaInMt53S); sigmaInSVec.push_back(&sigmaInMt54S); sigmaInSVec.push_back(&sigmaInMt55S);
        vector<double> inelasticEnergyLossS;
        inelasticEnergyLossS.push_back(2.23); inelasticEnergyLossS.push_back(3.77899); inelasticEnergyLossS.push_back(4.27999); inelasticEnergyLossS.push_back(4.45899); inelasticEnergyLossS.push_back(4.69499);

        static TMatrixF sigmaInMt51Ar = readSigmaEnergy(endfFolder, "inelasticMt51Ar40.txt");
        static TMatrixF sigmaInMt52Ar = readSigmaEnergy(endfFolder, "inelasticMt52Ar40.txt");
        static TMatrixF sigmaInMt53Ar = readSigmaEnergy(endfFolder, "inelasticMt53Ar40.txt");
        static TMatrixF sigmaInMt54Ar = readSigmaEnergy(endfFolder, "inelasticMt54Ar40.txt");
        static TMatrixF sigmaInMt55Ar = readSigmaEnergy(endfFolder, "inelasticMt55Ar40.txt");
        vector<TMatrixF*> sigmaInArVec;
        sigmaInArVec.push_back(&sigmaInMt51Ar); sigmaInArVec.push_back(&sigmaInMt52Ar); sigmaInArVec.push_back(&sigmaInMt53Ar); sigmaInArVec.push_back(&sigmaInMt54Ar); sigmaInArVec.push_back(&sigmaInMt55Ar);
        vector<double> inelasticEnergyLossAr;
        inelasticEnergyLossAr.push_back(1.4609); inelasticEnergyLossAr.push_back(2.1208); inelasticEnergyLossAr.push_back(2.5241); inelasticEnergyLossAr.push_back(2.893); inelasticEnergyLossAr.push_back(3.2083);

        static TMatrixF sigmaInMt51B10 = readSigmaEnergy(endfFolder, "inelasticMt51B10.txt");
        static TMatrixF sigmaInMt52B10 = readSigmaEnergy(endfFolder, "inelasticMt52B10.txt");
        static TMatrixF sigmaInMt53B10 = readSigmaEnergy(endfFolder, "inelasticMt53B10.txt");
        static TMatrixF sigmaInMt54B10 = readSigmaEnergy(endfFolder, "inelasticMt54B10.txt");
        vector<TMatrixF*> sigmaInB10Vec;
        sigmaInB10Vec.push_back(&sigmaInMt51B10); sigmaInB10Vec.push_back(&sigmaInMt52B10); sigmaInB10Vec.push_back(&sigmaInMt53B10); sigmaInB10Vec.push_back(&sigmaInMt54B10);
        vector<double> inelasticEnergyLossB10;
        inelasticEnergyLossB10.push_back(0.718299); inelasticEnergyLossB10.push_back(1.7402); inelasticEnergyLossB10.push_back(2.154); inelasticEnergyLossB10.push_back(3.58699);

        static TMatrixF sigmaInMt51B11 = readSigmaEnergy(endfFolder, "inelasticMt51B11.txt");
        static TMatrixF sigmaInMt52B11 = readSigmaEnergy(endfFolder, "inelasticMt52B11.txt");
        static TMatrixF sigmaInMt53B11 = readSigmaEnergy(endfFolder, "inelasticMt53B11.txt");
        static TMatrixF sigmaInMt54B11 = readSigmaEnergy(endfFolder, "inelasticMt54B11.txt");
        vector<TMatrixF*> sigmaInB11Vec;
        sigmaInB11Vec.push_back(&sigmaInMt51B11); sigmaInB11Vec.push_back(&sigmaInMt52B11); sigmaInB11Vec.push_back(&sigmaInMt53B11); sigmaInB11Vec.push_back(&sigmaInMt54B11);
        vector<double> inelasticEnergyLossB11;
        inelasticEnergyLossB11.push_back(2.12469); inelasticEnergyLossB11.push_back(4.44489); inelasticEnergyLossB11.push_back(5.02031); inelasticEnergyLossB11.push_back(6.7429);

        static TMatrixF sigmaInMt51F19 = readSigmaEnergy(endfFolder, "inelasticMt51F19.txt");
        static TMatrixF sigmaInMt52F19 = readSigmaEnergy(endfFolder, "inelasticMt52F19.txt");
        static TMatrixF sigmaInMt53F19 = readSigmaEnergy(endfFolder, "inelasticMt53F19.txt");
        static TMatrixF sigmaInMt54F19 = readSigmaEnergy(endfFolder, "inelasticMt54F19.txt");
        vector<TMatrixF*> sigmaInF19Vec;
        sigmaInF19Vec.push_back(&sigmaInMt51F19); sigmaInF19Vec.push_back(&sigmaInMt52F19); sigmaInF19Vec.push_back(&sigmaInMt53F19); sigmaInF19Vec.push_back(&sigmaInMt54F19);
        vector<double> inelasticEnergyLossF19;
        inelasticEnergyLossF19.push_back(0.1099); inelasticEnergyLossF19.push_back(0.197); inelasticEnergyLossF19.push_back(1.346); inelasticEnergyLossF19.push_back(1.459);

        static TMatrixF sigmaInMt51Gd155 = readSigmaEnergy(endfFolder, "inelasticMt51Gd155.txt");
        static TMatrixF sigmaInMt52Gd155 = readSigmaEnergy(endfFolder, "inelasticMt52Gd155.txt");
        static TMatrixF sigmaInMt53Gd155 = readSigmaEnergy(endfFolder, "inelasticMt53Gd155.txt");
        static TMatrixF sigmaInMt54Gd155 = readSigmaEnergy(endfFolder, "inelasticMt54Gd155.txt");
        vector<TMatrixF*> sigmaInGd155Vec;
        sigmaInGd155Vec.push_back(&sigmaInMt51Gd155); sigmaInGd155Vec.push_back(&sigmaInMt52Gd155); sigmaInGd155Vec.push_back(&sigmaInMt53Gd155); sigmaInGd155Vec.push_back(&sigmaInMt54Gd155);
        vector<double> inelasticEnergyLossGd155;
        inelasticEnergyLossGd155.push_back(0.06); inelasticEnergyLossGd155.push_back(0.08655); inelasticEnergyLossGd155.push_back(0.10531); inelasticEnergyLossGd155.push_back(0.10758);

        static TMatrixF sigmaInMt51Gd157 = readSigmaEnergy(endfFolder, "inelasticMt51Gd157.txt");
        static TMatrixF sigmaInMt52Gd157 = readSigmaEnergy(endfFolder, "inelasticMt52Gd157.txt");
        static TMatrixF sigmaInMt53Gd157 = readSigmaEnergy(endfFolder, "inelasticMt53Gd157.txt");
        static TMatrixF sigmaInMt54Gd157 = readSigmaEnergy(endfFolder, "inelasticMt54Gd157.txt");
        vector<TMatrixF*> sigmaInGd157Vec;
        sigmaInGd157Vec.push_back(&sigmaInMt51Gd157); sigmaInGd157Vec.push_back(&sigmaInMt52Gd157); sigmaInGd157Vec.push_back(&sigmaInMt53Gd157); sigmaInGd157Vec.push_back(&sigmaInMt54Gd157);
        vector<double> inelasticEnergyLossGd157;
        inelasticEnergyLossGd157.push_back(0.05453); inelasticEnergyLossGd157.push_back(0.06392); inelasticEnergyLossGd157.push_back(0.11572); inelasticEnergyLossGd157.push_back(0.3145);

        static TMatrixF sigmaInMt51Na = readSigmaEnergy(endfFolder, "inelasticMt51Na23.txt");
        static TMatrixF sigmaInMt52Na = readSigmaEnergy(endfFolder, "inelasticMt52Na23.txt");
        static TMatrixF sigmaInMt53Na = readSigmaEnergy(endfFolder, "inelasticMt53Na23.txt");
        static TMatrixF sigmaInMt54Na = readSigmaEnergy(endfFolder, "inelasticMt54Na23.txt");
        static TMatrixF sigmaInMt55Na = readSigmaEnergy(endfFolder, "inelasticMt55Na23.txt");
        static TMatrixF sigmaInMt56Na = readSigmaEnergy(endfFolder, "inelasticMt56Na23.txt");
        vector<TMatrixF*> sigmaInNaVec;
        sigmaInNaVec.push_back(&sigmaInMt51Na); sigmaInNaVec.push_back(&sigmaInMt52Na); sigmaInNaVec.push_back(&sigmaInMt53Na); sigmaInNaVec.push_back(&sigmaInMt54Na); sigmaInNaVec.push_back(&sigmaInMt55Na); sigmaInNaVec.push_back(&sigmaInMt56Na);
        vector<double> inelasticEnergyLossNa;
        inelasticEnergyLossNa.push_back(0.44); inelasticEnergyLossNa.push_back(2.087); inelasticEnergyLossNa.push_back(2.393); inelasticEnergyLossNa.push_back(2.64); inelasticEnergyLossNa.push_back(2.705); inelasticEnergyLossNa.push_back(2.983);

        static TMatrixF sigmaInMt51Cl35 = readSigmaEnergy(endfFolder, "inelasticMt51Cl35.txt");
        static TMatrixF sigmaInMt52Cl35 = readSigmaEnergy(endfFolder, "inelasticMt52Cl35.txt");
        static TMatrixF sigmaInMt53Cl35 = readSigmaEnergy(endfFolder, "inelasticMt53Cl35.txt");
        static TMatrixF sigmaInMt54Cl35 = readSigmaEnergy(endfFolder, "inelasticMt54Cl35.txt");
        static TMatrixF sigmaInMt55Cl35 = readSigmaEnergy(endfFolder, "inelasticMt55Cl35.txt");
        static TMatrixF sigmaInMt56Cl35 = readSigmaEnergy(endfFolder, "inelasticMt56Cl35.txt");
        vector<TMatrixF*> sigmaInCl35Vec;
        sigmaInCl35Vec.push_back(&sigmaInMt51Cl35); sigmaInCl35Vec.push_back(&sigmaInMt52Cl35); sigmaInCl35Vec.push_back(&sigmaInMt53Cl35); sigmaInCl35Vec.push_back(&sigmaInMt54Cl35); sigmaInCl35Vec.push_back(&sigmaInMt55Cl35); sigmaInCl35Vec.push_back(&sigmaInMt56Cl35);
        vector<double> inelasticEnergyLossCl35;
        inelasticEnergyLossCl35.push_back(1.21944); inelasticEnergyLossCl35.push_back(1.76315); inelasticEnergyLossCl35.push_back(2.6456); inelasticEnergyLossCl35.push_back(2.6936); inelasticEnergyLossCl35.push_back(3.00274); inelasticEnergyLossCl35.push_back(3.168);

        static TMatrixF sigmaInMt51Cr52 = readSigmaEnergy(endfFolder, "inelasticMt51Cr52.txt");
        static TMatrixF sigmaInMt52Cr52 = readSigmaEnergy(endfFolder, "inelasticMt52Cr52.txt");
        static TMatrixF sigmaInMt53Cr52 = readSigmaEnergy(endfFolder, "inelasticMt53Cr52.txt");
        static TMatrixF sigmaInMt54Cr52 = readSigmaEnergy(endfFolder, "inelasticMt54Cr52.txt");
        static TMatrixF sigmaInMt55Cr52 = readSigmaEnergy(endfFolder, "inelasticMt55Cr52.txt");
        vector<TMatrixF*> sigmaInCr52Vec;
        sigmaInCr52Vec.push_back(&sigmaInMt51Cr52); sigmaInCr52Vec.push_back(&sigmaInMt52Cr52); sigmaInCr52Vec.push_back(&sigmaInMt53Cr52); sigmaInCr52Vec.push_back(&sigmaInMt54Cr52); sigmaInCr52Vec.push_back(&sigmaInMt55Cr52);
        vector<double> inelasticEnergyLossCr52;
        inelasticEnergyLossCr52.push_back(1.434); inelasticEnergyLossCr52.push_back(2.37); inelasticEnergyLossCr52.push_back(2.647); inelasticEnergyLossCr52.push_back(2.768); inelasticEnergyLossCr52.push_back(2.965);

        static TMatrixF sigmaInMt51Cr53 = readSigmaEnergy(endfFolder, "inelasticMt51Cr53.txt");
        static TMatrixF sigmaInMt52Cr53 = readSigmaEnergy(endfFolder, "inelasticMt52Cr53.txt");
        static TMatrixF sigmaInMt53Cr53 = readSigmaEnergy(endfFolder, "inelasticMt53Cr53.txt");
        static TMatrixF sigmaInMt54Cr53 = readSigmaEnergy(endfFolder, "inelasticMt54Cr53.txt");
        static TMatrixF sigmaInMt55Cr53 = readSigmaEnergy(endfFolder, "inelasticMt55Cr53.txt");
        vector<TMatrixF*> sigmaInCr53Vec;
        sigmaInCr53Vec.push_back(&sigmaInMt51Cr53); sigmaInCr53Vec.push_back(&sigmaInMt52Cr53); sigmaInCr53Vec.push_back(&sigmaInMt53Cr53); sigmaInCr53Vec.push_back(&sigmaInMt54Cr53); sigmaInCr53Vec.push_back(&sigmaInMt55Cr53);
        vector<double> inelasticEnergyLossCr53;
        inelasticEnergyLossCr53.push_back(0.564); inelasticEnergyLossCr53.push_back(1.006); inelasticEnergyLossCr53.push_back(1.29);
        inelasticEnergyLossCr53.push_back(1.537); inelasticEnergyLossCr53.push_back(1.974);

        static TMatrixF sigmaInMt51Ni58 = readSigmaEnergy(endfFolder, "inelasticMt51Ni58.txt");
        static TMatrixF sigmaInMt52Ni58 = readSigmaEnergy(endfFolder, "inelasticMt52Ni58.txt");
        static TMatrixF sigmaInMt53Ni58 = readSigmaEnergy(endfFolder, "inelasticMt53Ni58.txt");
        static TMatrixF sigmaInMt54Ni58 = readSigmaEnergy(endfFolder, "inelasticMt54Ni58.txt");
        vector<TMatrixF*> sigmaInNi58Vec;
        sigmaInNi58Vec.push_back(&sigmaInMt51Ni58); sigmaInNi58Vec.push_back(&sigmaInMt52Ni58); sigmaInNi58Vec.push_back(&sigmaInMt53Ni58); sigmaInNi58Vec.push_back(&sigmaInMt54Ni58);
        vector<double> inelasticEnergyLossNi58;
        inelasticEnergyLossNi58.push_back(1.454); inelasticEnergyLossNi58.push_back(2.495); inelasticEnergyLossNi58.push_back(2.776); inelasticEnergyLossNi58.push_back(2.903);

        static TMatrixF sigmaInMt51Mn55 = readSigmaEnergy(endfFolder, "inelasticMt51Mn55.txt");
        static TMatrixF sigmaInMt52Mn55 = readSigmaEnergy(endfFolder, "inelasticMt52Mn55.txt");
        static TMatrixF sigmaInMt53Mn55 = readSigmaEnergy(endfFolder, "inelasticMt53Mn55.txt");
        static TMatrixF sigmaInMt54Mn55 = readSigmaEnergy(endfFolder, "inelasticMt54Mn55.txt");
        static TMatrixF sigmaInMt55Mn55 = readSigmaEnergy(endfFolder, "inelasticMt55Mn55.txt");
        static TMatrixF sigmaInMt56Mn55 = readSigmaEnergy(endfFolder, "inelasticMt56Mn55.txt");
        vector<TMatrixF*> sigmaInMn55Vec;
        sigmaInMn55Vec.push_back(&sigmaInMt51Mn55); sigmaInMn55Vec.push_back(&sigmaInMt52Mn55); sigmaInMn55Vec.push_back(&sigmaInMt53Mn55);
        sigmaInMn55Vec.push_back(&sigmaInMt54Mn55); sigmaInMn55Vec.push_back(&sigmaInMt55Mn55); sigmaInMn55Vec.push_back(&sigmaInMt56Mn55);
        vector<double> inelasticEnergyLossMn55;
        inelasticEnergyLossMn55.push_back(0.126); inelasticEnergyLossMn55.push_back(0.984); inelasticEnergyLossMn55.push_back(1.289);
        inelasticEnergyLossMn55.push_back(1.292); inelasticEnergyLossMn55.push_back(1.293); inelasticEnergyLossMn55.push_back(1.528);

        static TMatrixF sigmaInMt51Pb206 = readSigmaEnergy(endfFolder, "inelasticMt51Pb206.txt");
        static TMatrixF sigmaInMt52Pb206 = readSigmaEnergy(endfFolder, "inelasticMt52Pb206.txt");
        static TMatrixF sigmaInMt53Pb206 = readSigmaEnergy(endfFolder, "inelasticMt53Pb206.txt");
        static TMatrixF sigmaInMt54Pb206 = readSigmaEnergy(endfFolder, "inelasticMt54Pb206.txt");
        static TMatrixF sigmaInMt55Pb206 = readSigmaEnergy(endfFolder, "inelasticMt55Pb206.txt");

        vector<TMatrixF*> sigmaInPb206Vec;
        sigmaInPb206Vec.push_back(&sigmaInMt51Pb206); sigmaInPb206Vec.push_back(&sigmaInMt52Pb206); sigmaInPb206Vec.push_back(&sigmaInMt53Pb206);
        sigmaInPb206Vec.push_back(&sigmaInMt54Pb206); sigmaInPb206Vec.push_back(&sigmaInMt55Pb206);
        vector<double> inelasticEnergyLossPb206;
        inelasticEnergyLossPb206.push_back(0.8031); inelasticEnergyLossPb206.push_back(1.166); inelasticEnergyLossPb206.push_back(1.34054);
        inelasticEnergyLossPb206.push_back(1.4668); inelasticEnergyLossPb206.push_back(1.684);

        static TMatrixF sigmaInMt51Pb207 = readSigmaEnergy(endfFolder, "inelasticMt51Pb207.txt");
        static TMatrixF sigmaInMt52Pb207 = readSigmaEnergy(endfFolder, "inelasticMt52Pb207.txt");
        static TMatrixF sigmaInMt53Pb207 = readSigmaEnergy(endfFolder, "inelasticMt53Pb207.txt");
        static TMatrixF sigmaInMt54Pb207 = readSigmaEnergy(endfFolder, "inelasticMt54Pb207.txt");
        static TMatrixF sigmaInMt55Pb207 = readSigmaEnergy(endfFolder, "inelasticMt55Pb207.txt");

        vector<TMatrixF*> sigmaInPb207Vec;
        sigmaInPb207Vec.push_back(&sigmaInMt51Pb207); sigmaInPb207Vec.push_back(&sigmaInMt52Pb207); sigmaInPb207Vec.push_back(&sigmaInMt53Pb207);
        sigmaInPb207Vec.push_back(&sigmaInMt54Pb207); sigmaInPb207Vec.push_back(&sigmaInMt55Pb207);
        vector<double> inelasticEnergyLossPb207;
        inelasticEnergyLossPb207.push_back(0.5697); inelasticEnergyLossPb207.push_back(0.8978); inelasticEnergyLossPb207.push_back(1.633);
        inelasticEnergyLossPb207.push_back(2.3399); inelasticEnergyLossPb207.push_back(2.6235);

        static TMatrixF sigmaInMt51Pb208 = readSigmaEnergy(endfFolder, "inelasticMt51Pb208.txt");
        static TMatrixF sigmaInMt52Pb208 = readSigmaEnergy(endfFolder, "inelasticMt52Pb208.txt");
        static TMatrixF sigmaInMt53Pb208 = readSigmaEnergy(endfFolder, "inelasticMt53Pb208.txt");
        static TMatrixF sigmaInMt54Pb208 = readSigmaEnergy(endfFolder, "inelasticMt54Pb208.txt");
        static TMatrixF sigmaInMt55Pb208 = readSigmaEnergy(endfFolder, "inelasticMt55Pb208.txt");

        vector<TMatrixF*> sigmaInPb208Vec;
        sigmaInPb208Vec.push_back(&sigmaInMt51Pb208); sigmaInPb208Vec.push_back(&sigmaInMt52Pb208); sigmaInPb208Vec.push_back(&sigmaInMt53Pb208);
        sigmaInPb208Vec.push_back(&sigmaInMt54Pb208); sigmaInPb208Vec.push_back(&sigmaInMt55Pb208);
        vector<double> inelasticEnergyLossPb208;
        inelasticEnergyLossPb208.push_back(2.61455); inelasticEnergyLossPb208.push_back(3.1977); inelasticEnergyLossPb208.push_back(3.4751);
        inelasticEnergyLossPb208.push_back(3.7084); inelasticEnergyLossPb208.push_back(3.9198);

        static TMatrixF sigmaInMt51K39 = readSigmaEnergy(endfFolder, "inelasticMt51K39.txt");
        static TMatrixF sigmaInMt52K39 = readSigmaEnergy(endfFolder, "inelasticMt52K39.txt");
        static TMatrixF sigmaInMt53K39 = readSigmaEnergy(endfFolder, "inelasticMt53K39.txt");
        static TMatrixF sigmaInMt54K39 = readSigmaEnergy(endfFolder, "inelasticMt54K39.txt");
        vector<TMatrixF*> sigmaInK39Vec;
        sigmaInK39Vec.push_back(&sigmaInMt51K39); sigmaInK39Vec.push_back(&sigmaInMt52K39); sigmaInK39Vec.push_back(&sigmaInMt53K39); sigmaInK39Vec.push_back(&sigmaInMt54K39);
        vector<double> inelasticEnergyLossK39;
        inelasticEnergyLossK39.push_back(2.523); inelasticEnergyLossK39.push_back(2.825); inelasticEnergyLossK39.push_back(3.019); inelasticEnergyLossK39.push_back(3.5979);

        static TMatrixF sigmaInMt51Ti48 = readSigmaEnergy(endfFolder, "inelasticMt51Ti48.txt");
        static TMatrixF sigmaInMt52Ti48 = readSigmaEnergy(endfFolder, "inelasticMt52Ti48.txt");
        static TMatrixF sigmaInMt53Ti48 = readSigmaEnergy(endfFolder, "inelasticMt53Ti48.txt");
        static TMatrixF sigmaInMt54Ti48 = readSigmaEnergy(endfFolder, "inelasticMt54Ti48.txt");
        vector<TMatrixF*> sigmaInTi48Vec;
        sigmaInTi48Vec.push_back(&sigmaInMt51Ti48); sigmaInTi48Vec.push_back(&sigmaInMt52Ti48); sigmaInTi48Vec.push_back(&sigmaInMt53Ti48); sigmaInTi48Vec.push_back(&sigmaInMt54Ti48);
        vector<double> inelasticEnergyLossTi48;
        inelasticEnergyLossTi48.push_back(0.9835); inelasticEnergyLossTi48.push_back(2.2956); inelasticEnergyLossTi48.push_back(2.421); inelasticEnergyLossTi48.push_back(2.997);

        if (!noGUIMode) {uiM->setStatus(1, "Reading Inelastic Angular Cross Sections");        delay(5);}

        static TMatrixF angularInelasticMt51N = readAngularCoefficients(endfFolder, "angularInelasticMt51N14.txt");
        static TMatrixF angularInelasticMt52N = readAngularCoefficients(endfFolder, "angularInelasticMt52N14.txt");
        static TMatrixF angularInelasticMt53N = readAngularCoefficients(endfFolder, "angularInelasticMt53N14.txt");
        static TMatrixF angularInelasticMt54N = readAngularCoefficients(endfFolder, "angularInelasticMt54N14.txt");
        static TMatrixF angularInelasticMt55N = readAngularCoefficients(endfFolder, "angularInelasticMt55N14.txt");
        static TMatrixF angularInelasticMt56N = readAngularCoefficients(endfFolder, "angularInelasticMt56N14.txt");
        static TMatrixF angularInelasticMt57N = readAngularCoefficients(endfFolder, "angularInelasticMt57N14.txt");
        static TMatrixF angularInelasticMt58N = readAngularCoefficients(endfFolder, "angularInelasticMt58N14.txt");
        static TMatrixF angularInelasticMt59N = readAngularCoefficients(endfFolder, "angularInelasticMt59N14.txt");
        static TMatrixF angularInelasticMt60N = readAngularCoefficients(endfFolder, "angularInelasticMt60N14.txt");
        vector<TMatrixF> sigmaInNAngularVec;
        sigmaInNAngularVec.push_back(angularInelasticMt51N); sigmaInNAngularVec.push_back(angularInelasticMt52N); sigmaInNAngularVec.push_back(angularInelasticMt53N);
        sigmaInNAngularVec.push_back(angularInelasticMt54N); sigmaInNAngularVec.push_back(angularInelasticMt55N); sigmaInNAngularVec.push_back(angularInelasticMt56N);
        sigmaInNAngularVec.push_back(angularInelasticMt57N); sigmaInNAngularVec.push_back(angularInelasticMt58N); sigmaInNAngularVec.push_back(angularInelasticMt59N); sigmaInNAngularVec.push_back(angularInelasticMt60N);

        static TMatrixF angularInelasticMt51O = readAngularCoefficients(endfFolder, "angularInelasticMt51O16.txt");
        static TMatrixF angularInelasticMt52O = readAngularCoefficients(endfFolder, "angularInelasticMt52O16.txt");
        static TMatrixF angularInelasticMt53O = readAngularCoefficients(endfFolder, "angularInelasticMt53O16.txt");
        static TMatrixF angularInelasticMt54O = readAngularCoefficients(endfFolder, "angularInelasticMt54O16.txt");
        static TMatrixF angularInelasticMt55O = readAngularCoefficients(endfFolder, "angularInelasticMt55O16.txt");
        static TMatrixF angularInelasticMt56O = readAngularCoefficients(endfFolder, "angularInelasticMt56O16.txt");
        static TMatrixF angularInelasticMt57O = readAngularCoefficients(endfFolder, "angularInelasticMt57O16.txt"); //only copies of 56!
        static TMatrixF angularInelasticMt58O = readAngularCoefficients(endfFolder, "angularInelasticMt56O16.txt");
        static TMatrixF angularInelasticMt59O = readAngularCoefficients(endfFolder, "angularInelasticMt57O16.txt");
        static TMatrixF angularInelasticMt60O = readAngularCoefficients(endfFolder, "angularInelasticMt56O16.txt");
        static TMatrixF angularInelasticMt61O = readAngularCoefficients(endfFolder, "angularInelasticMt57O16.txt");
        static TMatrixF angularInelasticMt62O = readAngularCoefficients(endfFolder, "angularInelasticMt56O16.txt");
        static TMatrixF angularInelasticMt63O = readAngularCoefficients(endfFolder, "angularInelasticMt57O16.txt");
        static TMatrixF angularInelasticMt64O = readAngularCoefficients(endfFolder, "angularInelasticMt56O16.txt");
        static TMatrixF angularInelasticMt65O = readAngularCoefficients(endfFolder, "angularInelasticMt57O16.txt");
        vector<TMatrixF> sigmaInOAngularVec;
        sigmaInOAngularVec.push_back(angularInelasticMt51O); sigmaInOAngularVec.push_back(angularInelasticMt52O); sigmaInOAngularVec.push_back(angularInelasticMt53O); sigmaInOAngularVec.push_back(angularInelasticMt54O);
        sigmaInOAngularVec.push_back(angularInelasticMt55O); sigmaInOAngularVec.push_back(angularInelasticMt56O);  sigmaInOAngularVec.push_back(angularInelasticMt57O);
        sigmaInOAngularVec.push_back(angularInelasticMt58O); sigmaInOAngularVec.push_back(angularInelasticMt59O);  sigmaInOAngularVec.push_back(angularInelasticMt60O);
        sigmaInOAngularVec.push_back(angularInelasticMt61O); sigmaInOAngularVec.push_back(angularInelasticMt62O); sigmaInOAngularVec.push_back(angularInelasticMt63O);
        sigmaInOAngularVec.push_back(angularInelasticMt64O); sigmaInOAngularVec.push_back(angularInelasticMt65O); sigmaInOAngularVec.push_back(angularInelasticMt65O);
        sigmaInOAngularVec.push_back(angularInelasticMt65O); sigmaInOAngularVec.push_back(angularInelasticMt65O); sigmaInOAngularVec.push_back(angularInelasticMt65O);
        sigmaInOAngularVec.push_back(angularInelasticMt65O);

        static TMatrixF angularInelasticMt51Si = readAngularCoefficients(endfFolder, "angularInelasticMt51Si28.txt");
        static TMatrixF angularInelasticMt52Si = readAngularCoefficients(endfFolder, "angularInelasticMt52Si28.txt");
        static TMatrixF angularInelasticMt53Si = readAngularCoefficients(endfFolder, "angularInelasticMt53Si28.txt");
        static TMatrixF angularInelasticMt54Si = readAngularCoefficients(endfFolder, "angularInelasticMt54Si28.txt");
        static TMatrixF angularInelasticMt55Si = readAngularCoefficients(endfFolder, "angularInelasticMt55Si28.txt");
        static TMatrixF angularInelasticMt56Si = readAngularCoefficients(endfFolder, "angularInelasticMt56Si28.txt");
        static TMatrixF angularInelasticMt57Si = readAngularCoefficients(endfFolder, "angularInelasticMt57Si28.txt");
        static TMatrixF angularInelasticMt58Si = readAngularCoefficients(endfFolder, "angularInelasticMt58Si28.txt");
        vector<TMatrixF> sigmaInSiAngularVec;
        sigmaInSiAngularVec.push_back(angularInelasticMt51Si); sigmaInSiAngularVec.push_back(angularInelasticMt52Si); sigmaInSiAngularVec.push_back(angularInelasticMt53Si);
        sigmaInSiAngularVec.push_back(angularInelasticMt54Si); sigmaInSiAngularVec.push_back(angularInelasticMt55Si); sigmaInSiAngularVec.push_back(angularInelasticMt56Si);
        sigmaInSiAngularVec.push_back(angularInelasticMt57Si);  sigmaInSiAngularVec.push_back(angularInelasticMt58Si);
        //beware! mt52 angular distribution at mt53 because of tabulated form of 53! //deprecated 02/2018

        static TMatrixF angularInelasticMt51C = readAngularCoefficients(endfFolder, "angularInelasticMt51C12.txt");
        static TMatrixF angularInelasticMt52C = readAngularCoefficients(endfFolder, "angularInelasticMt52C12.txt");
        static TMatrixF angularInelasticMt53C = readAngularCoefficients(endfFolder, "angularInelasticMt53C12.txt");
        static TMatrixF angularInelasticMt58C = readAngularCoefficients(endfFolder, "angularInelasticMt58C12.txt");
        vector<TMatrixF> sigmaInCAngularVec;
        sigmaInCAngularVec.push_back(angularInelasticMt51C); sigmaInCAngularVec.push_back(angularInelasticMt52C); sigmaInCAngularVec.push_back(angularInelasticMt53C);
        sigmaInCAngularVec.push_back(angularInelasticMt53C); sigmaInCAngularVec.push_back(angularInelasticMt53C); sigmaInCAngularVec.push_back(angularInelasticMt53C); sigmaInCAngularVec.push_back(angularInelasticMt53C);
        sigmaInCAngularVec.push_back(angularInelasticMt58C);

        static TMatrixF angularInelasticMt51Al = readAngularCoefficients(endfFolder, "angularInelasticMt51Al27.txt");
        static TMatrixF angularInelasticMt52Al = readAngularCoefficients(endfFolder, "angularInelasticMt52Al27.txt");
        static TMatrixF angularInelasticMt53Al = readAngularCoefficients(endfFolder, "angularInelasticMt53Al27.txt");
        static TMatrixF angularInelasticMt54Al = readAngularCoefficients(endfFolder, "angularInelasticMt54Al27.txt");
        static TMatrixF angularInelasticMt55Al = readAngularCoefficients(endfFolder, "angularInelasticMt55Al27.txt");
        static TMatrixF angularInelasticMt56Al = readAngularCoefficients(endfFolder, "angularInelasticMt56Al27.txt");
        static TMatrixF angularInelasticMt57Al = readAngularCoefficients(endfFolder, "angularInelasticMt57Al27.txt");
        static TMatrixF angularInelasticMt58Al = readAngularCoefficients(endfFolder, "angularInelasticMt58Al27.txt");
        vector<TMatrixF> sigmaInAlAngularVec;
        sigmaInAlAngularVec.push_back(angularInelasticMt51Al); sigmaInAlAngularVec.push_back(angularInelasticMt52Al); sigmaInAlAngularVec.push_back(angularInelasticMt53Al);
        sigmaInAlAngularVec.push_back(angularInelasticMt54Al); sigmaInAlAngularVec.push_back(angularInelasticMt55Al); sigmaInAlAngularVec.push_back(angularInelasticMt56Al);
        sigmaInAlAngularVec.push_back(angularInelasticMt57Al);   sigmaInAlAngularVec.push_back(angularInelasticMt58Al);

        static TMatrixF angularInelasticMt51Fe = readAngularCoefficients(endfFolder, "angularInelasticMt51Fe56.txt");
        static TMatrixF angularInelasticMt52Fe = readAngularCoefficients(endfFolder, "angularInelasticMt52Fe56.txt");
        static TMatrixF angularInelasticMt53Fe = readAngularCoefficients(endfFolder, "angularInelasticMt53Fe56.txt");
        static TMatrixF angularInelasticMt54Fe = readAngularCoefficients(endfFolder, "angularInelasticMt54Fe56.txt");
        static TMatrixF angularInelasticMt55Fe = readAngularCoefficients(endfFolder, "angularInelasticMt55Fe56.txt");
        static TMatrixF angularInelasticMt56Fe = readAngularCoefficients(endfFolder, "angularInelasticMt56Fe56.txt");
        static TMatrixF angularInelasticMt57Fe = readAngularCoefficients(endfFolder, "angularInelasticMt57Fe56.txt");
        static TMatrixF angularInelasticMt58Fe = readAngularCoefficients(endfFolder, "angularInelasticMt58Fe56.txt");
        vector<TMatrixF> sigmaInFeAngularVec;
        sigmaInFeAngularVec.push_back(angularInelasticMt51Fe); sigmaInFeAngularVec.push_back(angularInelasticMt52Fe); sigmaInFeAngularVec.push_back(angularInelasticMt53Fe); sigmaInFeAngularVec.push_back(angularInelasticMt54Fe); sigmaInFeAngularVec.push_back(angularInelasticMt55Fe); sigmaInFeAngularVec.push_back(angularInelasticMt56Fe); sigmaInFeAngularVec.push_back(angularInelasticMt57Fe); sigmaInFeAngularVec.push_back(angularInelasticMt58Fe);

        static TMatrixF angularInelasticMt51S = readAngularCoefficients(endfFolder, "angularInelasticMt51S32.txt");
        static TMatrixF angularInelasticMt52S = readAngularCoefficients(endfFolder, "angularInelasticMt52S32.txt");
        static TMatrixF angularInelasticMt53S = readAngularCoefficients(endfFolder, "angularInelasticMt53S32.txt");
        static TMatrixF angularInelasticMt54S = readAngularCoefficients(endfFolder, "angularInelasticMt54S32.txt");
        static TMatrixF angularInelasticMt55S = readAngularCoefficients(endfFolder, "angularInelasticMt55S32.txt");
        vector<TMatrixF> sigmaInSAngularVec;
        sigmaInSAngularVec.push_back(angularInelasticMt51S); sigmaInSAngularVec.push_back(angularInelasticMt52S); sigmaInSAngularVec.push_back(angularInelasticMt53S); sigmaInSAngularVec.push_back(angularInelasticMt54S); sigmaInSAngularVec.push_back(angularInelasticMt55S);

        static TMatrixF angularInelasticMt51Ar = readAngularCoefficients(endfFolder, "angularInelasticMt51Ar40.txt");
        static TMatrixF angularInelasticMt52Ar = readAngularCoefficients(endfFolder, "angularInelasticMt52Ar40.txt");
        static TMatrixF angularInelasticMt53Ar = readAngularCoefficients(endfFolder, "angularInelasticMt53Ar40.txt");
        static TMatrixF angularInelasticMt54Ar = readAngularCoefficients(endfFolder, "angularInelasticMt54Ar40.txt");
        static TMatrixF angularInelasticMt55Ar = readAngularCoefficients(endfFolder, "angularInelasticMt55Ar40.txt");
        vector<TMatrixF> sigmaInArAngularVec;
        sigmaInArAngularVec.push_back(angularInelasticMt51Ar); sigmaInArAngularVec.push_back(angularInelasticMt52Ar); sigmaInArAngularVec.push_back(angularInelasticMt53Ar);
        sigmaInArAngularVec.push_back(angularInelasticMt54Ar); sigmaInArAngularVec.push_back(angularInelasticMt55Ar);

        static TMatrixF angularInelasticMt51B10 = readAngularCoefficients(endfFolder, "angularInelasticMt51B10.txt");
        static TMatrixF angularInelasticMt52B10 = readAngularCoefficients(endfFolder, "angularInelasticMt52B10.txt");
        static TMatrixF angularInelasticMt53B10 = readAngularCoefficients(endfFolder, "angularInelasticMt53B10.txt");
        static TMatrixF angularInelasticMt54B10 = readAngularCoefficients(endfFolder, "angularInelasticMt54B10.txt");
        vector<TMatrixF> sigmaInB10AngularVec;
        sigmaInB10AngularVec.push_back(angularInelasticMt51B10); sigmaInB10AngularVec.push_back(angularInelasticMt52B10); sigmaInB10AngularVec.push_back(angularInelasticMt53B10); sigmaInB10AngularVec.push_back(angularInelasticMt54B10);

        static TMatrixF angularInelasticMt51B11 = readAngularCoefficients(endfFolder, "angularInelasticMt51B11.txt");
        static TMatrixF angularInelasticMt52B11 = readAngularCoefficients(endfFolder, "angularInelasticMt52B11.txt");
        static TMatrixF angularInelasticMt53B11 = readAngularCoefficients(endfFolder, "angularInelasticMt53B11.txt");
        static TMatrixF angularInelasticMt54B11 = readAngularCoefficients(endfFolder, "angularInelasticMt54B11.txt");
        vector<TMatrixF> sigmaInB11AngularVec;
        sigmaInB11AngularVec.push_back(angularInelasticMt51B11); sigmaInB11AngularVec.push_back(angularInelasticMt52B11); sigmaInB11AngularVec.push_back(angularInelasticMt53B11); sigmaInB11AngularVec.push_back(angularInelasticMt54B11);

        static TMatrixF angularInelasticMt51F19 = readAngularCoefficients(endfFolder, "angularInelasticMt51F19.txt");
        static TMatrixF angularInelasticMt52F19 = readAngularCoefficients(endfFolder, "angularInelasticMt52F19.txt");
        static TMatrixF angularInelasticMt53F19 = readAngularCoefficients(endfFolder, "angularInelasticMt53F19.txt");
        static TMatrixF angularInelasticMt54F19 = readAngularCoefficients(endfFolder, "angularInelasticMt54F19.txt");
        vector<TMatrixF> sigmaInF19AngularVec;
        sigmaInF19AngularVec.push_back(angularInelasticMt51F19); sigmaInF19AngularVec.push_back(angularInelasticMt52F19); sigmaInF19AngularVec.push_back(angularInelasticMt53F19); sigmaInF19AngularVec.push_back(angularInelasticMt54F19);

        static TMatrixF angularInelasticMt51Gd155 = readAngularCoefficients(endfFolder, "angularInelasticMt51Gd155.txt");
        static TMatrixF angularInelasticMt52Gd155 = readAngularCoefficients(endfFolder, "angularInelasticMt52Gd155.txt");
        static TMatrixF angularInelasticMt53Gd155 = readAngularCoefficients(endfFolder, "angularInelasticMt53Gd155.txt");
        static TMatrixF angularInelasticMt54Gd155 = readAngularCoefficients(endfFolder, "angularInelasticMt54Gd155.txt");
        vector<TMatrixF> sigmaInGd155AngularVec;
        sigmaInGd155AngularVec.push_back(angularInelasticMt51Gd155); sigmaInGd155AngularVec.push_back(angularInelasticMt52Gd155); sigmaInGd155AngularVec.push_back(angularInelasticMt53Gd155); sigmaInGd155AngularVec.push_back(angularInelasticMt54Gd155);

        static TMatrixF angularInelasticMt51Gd157 = readAngularCoefficients(endfFolder, "angularInelasticMt51Gd157.txt");
        static TMatrixF angularInelasticMt52Gd157 = readAngularCoefficients(endfFolder, "angularInelasticMt52Gd157.txt");
        static TMatrixF angularInelasticMt53Gd157 = readAngularCoefficients(endfFolder, "angularInelasticMt53Gd157.txt");
        static TMatrixF angularInelasticMt54Gd157 = readAngularCoefficients(endfFolder, "angularInelasticMt54Gd157.txt");
        vector<TMatrixF> sigmaInGd157AngularVec;
        sigmaInGd157AngularVec.push_back(angularInelasticMt51Gd157); sigmaInGd157AngularVec.push_back(angularInelasticMt52Gd157); sigmaInGd157AngularVec.push_back(angularInelasticMt53Gd157); sigmaInGd157AngularVec.push_back(angularInelasticMt54Gd157);

        static TMatrixF angularInelasticMt51Na = readAngularCoefficients(endfFolder, "angularInelasticMt51Na23.txt");
        static TMatrixF angularInelasticMt52Na = readAngularCoefficients(endfFolder, "angularInelasticMt52Na23.txt");
        static TMatrixF angularInelasticMt53Na = readAngularCoefficients(endfFolder, "angularInelasticMt53Na23.txt");
        static TMatrixF angularInelasticMt54Na = readAngularCoefficients(endfFolder, "angularInelasticMt54Na23.txt");
        static TMatrixF angularInelasticMt55Na = readAngularCoefficients(endfFolder, "angularInelasticMt55Na23.txt");
        static TMatrixF angularInelasticMt56Na = readAngularCoefficients(endfFolder, "angularInelasticMt56Na23.txt");
        vector<TMatrixF> sigmaInNaAngularVec;
        sigmaInNaAngularVec.push_back(angularInelasticMt51Na); sigmaInNaAngularVec.push_back(angularInelasticMt52Na); sigmaInNaAngularVec.push_back(angularInelasticMt53Na);
        sigmaInNaAngularVec.push_back(angularInelasticMt54Na); sigmaInNaAngularVec.push_back(angularInelasticMt55Na); sigmaInNaAngularVec.push_back(angularInelasticMt56Na);

        static TMatrixF angularInelasticMt51Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt51Cl35.txt");
        static TMatrixF angularInelasticMt52Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt52Cl35.txt");
        static TMatrixF angularInelasticMt53Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt53Cl35.txt");
        static TMatrixF angularInelasticMt54Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt54Cl35.txt");
        static TMatrixF angularInelasticMt55Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt55Cl35.txt");
        static TMatrixF angularInelasticMt56Cl35 = readAngularCoefficients(endfFolder, "angularInelasticMt56Cl35.txt");
        vector<TMatrixF> sigmaInCl35AngularVec;
        sigmaInCl35AngularVec.push_back(angularInelasticMt51Cl35); sigmaInCl35AngularVec.push_back(angularInelasticMt52Cl35); sigmaInCl35AngularVec.push_back(angularInelasticMt53Cl35);
        sigmaInCl35AngularVec.push_back(angularInelasticMt54Cl35); sigmaInCl35AngularVec.push_back(angularInelasticMt55Cl35); sigmaInCl35AngularVec.push_back(angularInelasticMt56Cl35);

        static TMatrixF angularInelasticMt51Cr52 = readAngularCoefficients(endfFolder, "angularInelasticMt51Cr52.txt");
        static TMatrixF angularInelasticMt52Cr52 = readAngularCoefficients(endfFolder, "angularInelasticMt52Cr52.txt");
        static TMatrixF angularInelasticMt53Cr52 = readAngularCoefficients(endfFolder, "angularInelasticMt53Cr52.txt");
        static TMatrixF angularInelasticMt54Cr52 = readAngularCoefficients(endfFolder, "angularInelasticMt54Cr52.txt");
        static TMatrixF angularInelasticMt55Cr52 = readAngularCoefficients(endfFolder, "angularInelasticMt55Cr52.txt");
        vector<TMatrixF> sigmaInCr52AngularVec;
        sigmaInCr52AngularVec.push_back(angularInelasticMt51Cr52); sigmaInCr52AngularVec.push_back(angularInelasticMt52Cr52); sigmaInCr52AngularVec.push_back(angularInelasticMt53Cr52);
        sigmaInCr52AngularVec.push_back(angularInelasticMt54Cr52); sigmaInCr52AngularVec.push_back(angularInelasticMt55Cr52);

        static TMatrixF angularInelasticMt51Cr53 = readAngularCoefficients(endfFolder, "angularInelasticMt51Cr53.txt");
        static TMatrixF angularInelasticMt52Cr53 = readAngularCoefficients(endfFolder, "angularInelasticMt52Cr53.txt");
        static TMatrixF angularInelasticMt53Cr53 = readAngularCoefficients(endfFolder, "angularInelasticMt53Cr53.txt");
        static TMatrixF angularInelasticMt54Cr53 = readAngularCoefficients(endfFolder, "angularInelasticMt54Cr53.txt");
        static TMatrixF angularInelasticMt55Cr53 = readAngularCoefficients(endfFolder, "angularInelasticMt55Cr53.txt");
        vector<TMatrixF> sigmaInCr53AngularVec;
        sigmaInCr53AngularVec.push_back(angularInelasticMt51Cr53); sigmaInCr53AngularVec.push_back(angularInelasticMt52Cr53); sigmaInCr53AngularVec.push_back(angularInelasticMt53Cr53);
        sigmaInCr53AngularVec.push_back(angularInelasticMt54Cr53); sigmaInCr53AngularVec.push_back(angularInelasticMt55Cr53);

        static TMatrixF angularInelasticMt51Ni58 = readAngularCoefficients(endfFolder, "angularInelasticMt51Ni58.txt");
        static TMatrixF angularInelasticMt52Ni58 = readAngularCoefficients(endfFolder, "angularInelasticMt52Ni58.txt");
        static TMatrixF angularInelasticMt53Ni58 = readAngularCoefficients(endfFolder, "angularInelasticMt53Ni58.txt");
        static TMatrixF angularInelasticMt54Ni58 = readAngularCoefficients(endfFolder, "angularInelasticMt54Ni58.txt");
        vector<TMatrixF> sigmaInNi58AngularVec;
        sigmaInNi58AngularVec.push_back(angularInelasticMt51Ni58); sigmaInNi58AngularVec.push_back(angularInelasticMt52Ni58); sigmaInNi58AngularVec.push_back(angularInelasticMt53Ni58);
        sigmaInNi58AngularVec.push_back(angularInelasticMt54Ni58);

        static TMatrixF angularInelasticMt51Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt51Mn55.txt");
        static TMatrixF angularInelasticMt52Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt52Mn55.txt");
        static TMatrixF angularInelasticMt53Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt53Mn55.txt");
        static TMatrixF angularInelasticMt54Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt54Mn55.txt");
        static TMatrixF angularInelasticMt55Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt55Mn55.txt");
        static TMatrixF angularInelasticMt56Mn55 = readAngularCoefficients(endfFolder, "angularInelasticMt56Mn55.txt");
        vector<TMatrixF> sigmaInMn55AngularVec;
        sigmaInMn55AngularVec.push_back(angularInelasticMt51Mn55); sigmaInMn55AngularVec.push_back(angularInelasticMt52Mn55); sigmaInMn55AngularVec.push_back(angularInelasticMt53Mn55);
        sigmaInMn55AngularVec.push_back(angularInelasticMt54Mn55); sigmaInMn55AngularVec.push_back(angularInelasticMt55Mn55); sigmaInMn55AngularVec.push_back(angularInelasticMt56Mn55);

        static TMatrixF angularInelasticMt51Pb206 = readAngularCoefficients(endfFolder, "angularInelasticMt51Pb206.txt");
        static TMatrixF angularInelasticMt52Pb206 = readAngularCoefficients(endfFolder, "angularInelasticMt52Pb206.txt");
        static TMatrixF angularInelasticMt53Pb206 = readAngularCoefficients(endfFolder, "angularInelasticMt53Pb206.txt");
        static TMatrixF angularInelasticMt54Pb206 = readAngularCoefficients(endfFolder, "angularInelasticMt54Pb206.txt");
        static TMatrixF angularInelasticMt55Pb206 = readAngularCoefficients(endfFolder, "angularInelasticMt55Pb206.txt");

        vector<TMatrixF> sigmaInPb206AngularVec;
        sigmaInPb206AngularVec.push_back(angularInelasticMt51Pb206); sigmaInPb206AngularVec.push_back(angularInelasticMt52Pb206); sigmaInPb206AngularVec.push_back(angularInelasticMt53Pb206);
        sigmaInPb206AngularVec.push_back(angularInelasticMt54Pb206); sigmaInPb206AngularVec.push_back(angularInelasticMt55Pb206);

        static TMatrixF angularInelasticMt51Pb207 = readAngularCoefficients(endfFolder, "angularInelasticMt51Pb207.txt");
        static TMatrixF angularInelasticMt52Pb207 = readAngularCoefficients(endfFolder, "angularInelasticMt52Pb207.txt");
        static TMatrixF angularInelasticMt53Pb207 = readAngularCoefficients(endfFolder, "angularInelasticMt53Pb207.txt");
        static TMatrixF angularInelasticMt54Pb207 = readAngularCoefficients(endfFolder, "angularInelasticMt54Pb207.txt");
        static TMatrixF angularInelasticMt55Pb207 = readAngularCoefficients(endfFolder, "angularInelasticMt55Pb207.txt");

        vector<TMatrixF> sigmaInPb207AngularVec;
        sigmaInPb207AngularVec.push_back(angularInelasticMt51Pb207); sigmaInPb207AngularVec.push_back(angularInelasticMt52Pb207); sigmaInPb207AngularVec.push_back(angularInelasticMt53Pb207);
        sigmaInPb207AngularVec.push_back(angularInelasticMt54Pb207); sigmaInPb207AngularVec.push_back(angularInelasticMt55Pb207);

        static TMatrixF angularInelasticMt51Pb208 = readAngularCoefficients(endfFolder, "angularInelasticMt51Pb208.txt");
        static TMatrixF angularInelasticMt52Pb208 = readAngularCoefficients(endfFolder, "angularInelasticMt52Pb208.txt");
        static TMatrixF angularInelasticMt53Pb208 = readAngularCoefficients(endfFolder, "angularInelasticMt53Pb208.txt");
        static TMatrixF angularInelasticMt54Pb208 = readAngularCoefficients(endfFolder, "angularInelasticMt54Pb208.txt");
        static TMatrixF angularInelasticMt55Pb208 = readAngularCoefficients(endfFolder, "angularInelasticMt55Pb208.txt");

        vector<TMatrixF> sigmaInPb208AngularVec;
        sigmaInPb208AngularVec.push_back(angularInelasticMt51Pb208); sigmaInPb208AngularVec.push_back(angularInelasticMt52Pb208); sigmaInPb208AngularVec.push_back(angularInelasticMt53Pb208);
        sigmaInPb208AngularVec.push_back(angularInelasticMt54Pb208); sigmaInPb208AngularVec.push_back(angularInelasticMt55Pb208);

        static TMatrixF angularInelasticMt51K39 = readAngularCoefficients(endfFolder, "angularInelasticMt51K39.txt");
        static TMatrixF angularInelasticMt52K39 = readAngularCoefficients(endfFolder, "angularInelasticMt52K39.txt");
        static TMatrixF angularInelasticMt53K39 = readAngularCoefficients(endfFolder, "angularInelasticMt53K39.txt");
        static TMatrixF angularInelasticMt54K39 = readAngularCoefficients(endfFolder, "angularInelasticMt54K39.txt");
        vector<TMatrixF> sigmaInK39AngularVec;
        sigmaInK39AngularVec.push_back(angularInelasticMt51K39); sigmaInK39AngularVec.push_back(angularInelasticMt52K39); sigmaInK39AngularVec.push_back(angularInelasticMt53K39); sigmaInK39AngularVec.push_back(angularInelasticMt54K39);

        static TMatrixF angularInelasticMt51Ti48 = readAngularCoefficients(endfFolder, "angularInelasticMt51Ti48.txt");
        static TMatrixF angularInelasticMt52Ti48 = readAngularCoefficients(endfFolder, "angularInelasticMt52Ti48.txt");
        static TMatrixF angularInelasticMt53Ti48 = readAngularCoefficients(endfFolder, "angularInelasticMt53Ti48.txt");
        static TMatrixF angularInelasticMt54Ti48 = readAngularCoefficients(endfFolder, "angularInelasticMt54Ti48.txt");
        vector<TMatrixF> sigmaInTi48AngularVec;
        sigmaInTi48AngularVec.push_back(angularInelasticMt51Ti48); sigmaInTi48AngularVec.push_back(angularInelasticMt52Ti48); sigmaInTi48AngularVec.push_back(angularInelasticMt53Ti48); sigmaInTi48AngularVec.push_back(angularInelasticMt54Ti48);


        if (!noGUIMode) cout << endl;
        if (!noGUIMode) {uiM->setStatus(1, "Starting simulation...");        delay(5);}

        // for parameter batch run
        int parInt1, parInt2;
        const int arrayLength = 11;
        //double soilWaterFracs[10] = {0.03,0.05,0.07,0.10,0.15,0.20,0.25,0.30,0.4,0.5};
        //double snowHeights[10] = {0.001,1,2,5,10,15,20,25,30,50};

        //double soilWaterFracs[10] = {0.01,0.02,0.04,0.08,0.12,0.18,0.22,0.27,0.35,0.45};
        //double snowHeights[10] = {0.001,3,4,7,12,18,40,60,80,100};

        //double soilWaterFracs[15] = {0.02001,0.03001,0.04001,0.05001,0.06001,0.08001,0.10001, 0.12001,0.15001,0.20001, 0.25001,0.30001,0.35001,0.45001,0.50};
        //double snowHeights[15] = {0.0001,1,3,5,7,10,15,20,25,30,40,60,80,100,150};

        //double soilWaterFracs[7] = {0.02001,0.05001,0.10001, 0.15001, 0.25001,0.35001,0.50};
        //double snowHeights[7] = {175,200,250,300,400,500,1000};

        //double soilWaterFracs[arrayLength] = {0.02001,0.03001,0.04001,0.05001,0.06001,0.08001,0.10001, 0.12001,0.15001,0.20001, 0.25001,0.30001,0.35001,0.45001,0.50};
        //double snowHeights[arrayLength] = {0.0001,1,5,10,15,20,30,60,80,100,150,250,400,500,1000};

        //double soilWaterFracs[arrayLength] = {0.04001,0.08001,0.10001, 0.12001,0.15001,0.20001, 0.25001,0.30001,0.35001,0.50};
        double snowHeights[arrayLength] = {0.0001,2,4,6,8,10,12,14,16,18 };

        double soilWaterFracs[arrayLength] = { 0.01001, 0.03001,0.06001,0.10001, 0.12001,0.15001,0.20001, 0.25001,0.30001,0.40001,0.50001 };
        //double soilWaterFracs[arrayLength] = {0.70001, 0.90001,10.0001,0.000, 0.0,0.0,0.0, 0.0,0.0,0.0000,0.0};
        double airhums[arrayLength] = {1,2,4,7,10,12,15,18,21,27,35 };

        // this for loop around the whole simulation is the parameter-defined batch run
        for (param = paramMin; param < paramMax; param += 1)
        {
            paramInt = param;
            if (doDetectorAngleBatchRun)  paramEnergy = 1e-4;//
            //if (doBatchRun) paramEnergy = TMath::Power(10,r.Rndm()*10.-8.);

            if ((doDetectorBatchRun) || ((doDetectorBatchRun2))) paramEnergy = TMath::Power(10, (param / (paramMax - paramMin)) * 10. - 7.99);
            if ((!doBatchRun) && (!doDetectorBatchRun) && (!doDetectorBatchRun2) && (!doDetectorAngleBatchRun) && (!doBatchRunDensity) && (!doBatchRun2D)) param = 1e9;

            if (doBatchRun)
            {
                //xr = xPosSource + (r.Rndm()*2.-1.)*0.5*xSizeSource;
                //yr = yPosSource + (r.Rndm()*2.-1.)*0.5*ySizeSource;
               //rHe3 = 0.000125 * pi/4. *param/10.;
            }

            if (doBatchRunDensity)
            {
                parInt1 = param / arrayLength;
                parInt2 = paramInt % arrayLength;
                paramDensity = snowHeights[parInt1];
                soilWaterFrac = soilWaterFracs[parInt2];
                soilWaterFracVar = soilWaterFrac;
                //if (parInt2 > 2) continue;
            }

            if (doBatchRun2D)
            {
                parInt1 = param / arrayLength;
                parInt2 = paramInt % arrayLength;
                soilWaterFrac = soilWaterFracs[parInt2];
                soilWaterFracVar = soilWaterFrac;
                //relHumidityAir = (1*airhums[parInt1])/10.*0.6;
                absHumidityAir = (1. * airhums[parInt1]);
                //if (parInt2 > 2) continue;
            }

            if (detFileOutput)
            {
                detOutputFile.open(outputFolder + "detectorNeutronHitData.dat", ios::out | ios::app);
                //detOutputFile << "Detector ID" <<"\t"<<"Neutron number"<<"\t"<<"Number of scatterings"<<"\t"<< "previous x [m]" <<"\t"<< "previous y [m]" <<"\t"<< "previous depth [m]" <<"\t"<< "Nadir angle" <<"\t"<< "Azimuth angle" <<"\t"<< "Energy [MeV]" <<"\t"<< "Energy at interface [MeV]" <<"\t"<< "Footprint crow-flight distance [m] (deprecated)" << "\t"<< "x at interface [m]" <<"\t"<< "y at interface [m]" << "\t"<<"z at interface [m]" <<"\t"<< "Maximum depth [m]" <<"\t"<<"Moisture at interface"<<"\t"<<"Counted using Physics Model"<<"\t"<<"Soil contact"<<endl;
                detOutputFile << "Detector_ID" << "\t" << "Neutron_Number" << "\t" << "Number_of_Scatterings" << "\t" << "previous_x_[m]" << "\t" << "previous_y_[m]" << "\t" << "previous_Depth_[m]" << "\t" << "Nadir_Angle" << "\t" << "Azimuth_Angle" << "\t" << "Energy_[MeV]" << "\t" << "Energy_at_Interface_[MeV]" << "\t" << "Footprint_Crow-Flight_Distance_[m]_(deprecated)" << "\t" << "x_at_interface_[m]" << "\t" << "y_at_Interface_[m]" << "\t" << "z_at_Interface_[m]" << "\t" << "previous_x_in_Soil_[m]" << "\t" << "previous_y_in_Soil_[m]" << "\t" << "previous_Depth_in_Soil_[m]" << "\t" << "previous_x_thermalized_[m]" << "\t" << "previous_y_thermalized_[m]" << "\t" << "previous_Depth_thermalized_[m]" << "\t" << "maximum_Depth_[m]" << "\t" << "mean_Probe_Depth_[m]"<< "\t" << "Time_[mus]" << "\t" << "moisture_at_Interface" << "\t" << "Recorded_using_Physics_Model" << "\t" << "Soil_Contact" << endl;
                detOutputFile.close();
            }

            if (detLayerFileOutput)
            {
                detLayerOutputFile.open(outputFolder + "detectorLayerNeutronHitData.dat", ios::out | ios::app);
                //detLayerOutputFile << "Neutron number"<<"\t"<<"Number of scatterings"<<"\t"<< "x [m]" <<"\t"<<"y [m]" <<"\t"<<"z [m]" <<"\t"<<"previous x [m]" <<"\t"<< "previous y [m]" <<"\t"<< "previous depth [m]" <<"\t"<< "Nadir angle" <<"\t"<< "Azimuth angle" <<"\t"<< "Energy [MeV]" <<"\t"<< "Energy at interface [MeV]" <<"\t"<< "Footprint crow-flight distance [m] (deprecated)" << "\t"<< "x at interface [m]" <<"\t"<< "y at interface [m]" << "\t"<<"z at interface [m]" <<"\t"<< "Maximum depth [m]" <<"\t"<<"Moisture at interface"<<"\t"<<"Counted using Physics Model"<<"\t"<<"Soil contact"<<endl;
                detLayerOutputFile << "Neutron_Number" << "\t" << "Number_of_Scatterings" << "\t" << "x_[m]" << "\t" << "y_[m]" << "\t" << "z_[m]" << "\t" << "previous_x_[m]" << "\t" << "previous_y_[m]" << "\t" << "previous_Depth_[m]" << "\t" << "Nadir_Angle" << "\t" << "Azimuth_Angle" << "\t" << "Energy_[MeV]" << "\t" << "Energy_at_Interface_[MeV]" << "\t" << "Footprint_Crow-Flight_Distance_[m]_(deprecated)" << "\t" << "x_at_Interface_[m]" << "\t" << "y_at_Interface_[m]" << "\t" << "z_at_Interface_[m]" << "\t" << "previous_x_in_Soil [m]" << "\t" << "previous_y_in_Soil_[m]" << "\t" << "previous_Depth_in_Soil_[m]" << "\t" << "previous_x_thermalized_[m]" << "\t" << "previous_y_thermalized_[m]" << "\t" << "previous_Depth_thermalized_[m]" << "\t" << "maximum_Depth_[m]" << "mean_Probe_Depth_[m]"<< "\t" <<  "\t" << "Time_[mus]" << "\t" << "moisture_at_Interface" << "\t" << "Recorded_using_Physics_Model" << "\t" << "Soil_contact" << endl;
                detLayerOutputFile.close();
            }

            if (detTrackFileOutput)
            {
                detTrackOutputFile.open(outputFolder + "detectorNeutronTrackHitData.dat", ios::out | ios::app);
                detTrackOutputFile << "Detector_ID" << "\t" << "Neutron_Number" << "\t" << "x_[m]" << "\t" << "y_[m]" << "\t" << "z_[m]" << "\t" << "Nadir_Angle" << "\t" << "Azimuth_Angle" << "\t" << "Energy_[MeV]" << "\t" << "Time_[mus]" << "\t" << "Element_ID_[ZZAA]" << "\t" << "Material_ID" << "\t" << "Counted_using_Physics_Model" << endl;
                detTrackOutputFile.close();
            }

            if (detLayerTrackFileOutput)
            {
                detLayerTrackOutputFile.open(outputFolder + "detectorNeutronTrackHitData.dat", ios::out | ios::app);
                detLayerTrackOutputFile << "Detector_ID" << "\t" << "Neutron_Number" << "\t" << "x_[m]" << "\t" << "y_[m]" << "\t" << "z_[m]" << "\t" << "Nadir_Angle" << "\t" << "Azimuth_Angle" << "\t" << "Energy_[MeV]" << "\t" << "Time_[mus]" << "\t" << "Element_ID_[ZZAA]" << "\t" << "Material_ID" << "\t" << "Counted_using_Physics_Model" << endl;
                detLayerTrackOutputFile.close();
            }

            if (allTrackFileOutput)
            {
                allTrackOutputFile.open(outputFolder + "NeutronTrackData.dat", ios::out | ios::app);
                allTrackOutputFile << "Neutron_Number" << "\t" << "x_[m]" << "\t" << "y_[m]" << "\t" << "z_[m]" << "\t" << "Nadir_Angle" << "\t" << "Azimuth_Angle" << "\t" << "Energy_[MeV]" << "\t" << "Time_[mus]" << "\t" << "Element_ID_[ZZAA]" << "\t" << "Material_ID" << endl;
                allTrackOutputFile.close();
            }

            detectorHeight = geometries.at(detectorLayer)[4] + 0.5 * geometries.at(detectorLayer)[5];

            if (uranosRootOutput) { dontFillunecessaryPlots = false; }
            else { dontFillunecessaryPlots = true; }

            if (!noThermalRegime) { maxScatterings = 700; }

            mapMetricFactor = squareDim / 1000. / 500.;
            started = false;

            // start the timer
            time(&start);
            time(&diffmean);
            time(&oldTime);
            pauseTime = 0;

            aspectRatioSide = (domainZUpperEdge - domainZLowEdge) / squareDim * 1000.;
            if (refreshCycle >= neutrons) refreshCycle = neutrons * 0.1;

            float eFoldingDepth = 150. / rBoden * 10.; //specifically used for distributiong a source in the ground
            bool continueTracking = true;
            nEvapoVector.clear();

            relHumidityAirVar = relHumidityAir;
            absHumidityAirVar = absHumidityAir;

            alreadyStarted = true;

            // //////////////////// MAIN LOOP
            long n;
            // start the simulation
            for (n = 0; n < neutrons + 1; n++)
            {
                if (n == 0)
                {
                    if (!noGUIMode) {uiM->setStatus(1, "");                    delay(5);}
                }
                totalActualNeutrons = n;

                // event display output
                if ((n % refreshCycle == 0) && (n > 9))
                {
                    nTotal = n;
                    declareNewData();
                    time(&diffmean);
                    //time (&actualTime);
                    newDataComes = true;

                    if (!noGUIMode) uiM->redrawNeutronMap(difftime(diffmean, start) - pauseTime);
                    //uiM->redrawNeutronMap(difftime(diffmean,oldTime));
                    //time (&oldTime);
                }

                if ((n % (refreshCycle / 5 ) == 0) && (n > 9))
                {
                    if (!noGUIMode) delay(1);
                }

                if ((n % refreshCycle == 0) && (n > 99) && (n < 5100))
                {
                    //uiM->setStatus(1,"");
                    //delay(5);
                }

                if ((saveEveryXNeutrons) && (n % refreshCycle == 0) && (n % ((int)TMath::Power(10, saveEveryXNeutronsPower) * 1) == 0) && (n > 10))
                {
                    exportTemporary = true;
                    uiM->exportToSave();
                    delay(300);
                    exportTemporary = false;
                }

                if ((false) && (n > 2))
                {
                    densityTrackMapSide->SaveAs(outputFolder + "/singleGraphs/" + castLongToString(n) + ".root");
                }

                if (stopRunning)
                {
                    stopRunning = false;
                    n = neutrons + 10;
                    //neutrons = 1;
                    param = 1e9;
                    continue;
                }

                // initialize variables for the neutron
                theta = thetaBeam;
                detectorHit = false;
                maxlayer = 0;
                hasBeenReorded = false;
                hasPassedSurface = false;
                hasBeenInSoil = false;
                sAbsorbed = false;
                measuredOnce = false;
                scatteredEvaporating = false;
                hasbeenEvaporized = false;
                isEvaporationNeutron = false;
                slowedDown = false;
                continuedThisLayer = false;
                scatteredThisLayer = false;
                leaveLayer = false;
                intoThermalization = false;
                //recordedThisLayer = false;

                if ((drawSingleNeutronGraphs) && (!started))
                {
                    graphCounter = 1;
                    graphs.clear();
                    if (neutronPath != 0x0) delete neutronPath;		if (multigraph != 0x0) delete multigraph;
                    neutronPath = new TGraph();
                    multigraph = new TMultiGraph();
                    graphs.push_back(new TGraph());
                    graphN = 0;
                    graphCounterMG = 0;

                    //activate this if all Multigraph paths shall be painted in one graph
                    //started = true;
                }

                subsurfaceScatterings.clear();
                subsurfaceScatteringMean = 0;

                //generate neutron Energy according to cosmic spectrum

                // initialize energy
                // store in xRnd
                xRnd = 1;
                if (!(nEvapoVector.size() > 0))
                {
                    if (usePrecalculatedAsciiSpectrum)
                    {
                        gotIt = false;
                        gotItCounter = 0;
                        while (!gotIt)
                        {
                            if (noHighEnergyRegime) xRnd = TMath::Power(10, r.Rndm() * 9. - 7.7);
                            else
                            {
                                if (noThermalRegime) xRnd = TMath::Power(10, r.Rndm() * 12. - 7.7);
                                else  xRnd = TMath::Power(10, r.Rndm() * 12.3 - 8.);
                            }

                            yRnd = r.Rndm();

                            if (spectrumModel->Eval(TMath::Log10(xRnd)) > yRnd)
                            {
                                gotIt = true;
                            }
                            gotItCounter++;
                            if (gotItCounter > 1000) gotIt = true;
                        }
                    }

                    if (usePrecalculatedSpectrum)
                    {
                        gotIt = false;
                        gotItCounter = 0;
                        while (!gotIt)
                        {
                            if (noHighEnergyRegime) xRnd = TMath::Power(10, r.Rndm() * 9. - 7.7);
                            else
                            {
                                if (noThermalRegime) xRnd = TMath::Power(10, r.Rndm() * 12. - 7.7);
                                else  xRnd = TMath::Power(10, r.Rndm() * 12.3 - 8.);
                            }
                            yRnd = r.Rndm();

                            binNo = precalculatedSpectrum->FindBin(xRnd);

                            if (precalculatedSpectrum->GetBinContent(binNo) > yRnd)
                            {
                                gotIt = true;
                            }
                            if (useHomogenousSpectrum) gotIt = true;
                            gotItCounter++;
                            if (gotItCounter > 1000) gotIt = true;
                        }
                    }
                    else
                    {
                        if (useSato2016)
                        {
                            if (usePrecalculatedSato2016Spectrum)
                            {
                                //calculates array of angles, otherwise takes precalculated
                                if (satoCounter < satoArrayLength)
                                {
                                    theta = satoAngleArray[satoCounter];
                                }
                                else
                                {
                                    for (int sp = 0; sp < satoArrayLength; sp++)
                                    {
                                        gotIt = false;
                                        gotItCounter = 0;
                                        // calculates the angle first
                                        while (!gotIt)
                                        {
                                            xRnd = r.Rndm() * piHalf;
                                            yRnd = r.Rndm();
                                            //angle distribution from Sato (integrated over the spectrum)
                                            // that 1.9 gives the angular distribution from Nesterenok
                                            if (yRnd < exp(-1.9 * (1. - cos(xRnd)))) { gotIt = true; }
                                            gotItCounter++;
                                            if (gotItCounter > 1000) gotIt = true;
                                        }
                                        theta = thetaBeam + xRnd;
                                        satoAngleArray[sp] = theta;
                                        //cout<<satoAngleArray[sp]<<" "<<gotItCounter<<" "<<yRnd<<" "<<exp(-1.9 * (1. - cos(xRnd)))<<endl;
                                    }
                                    theta = satoAngleArray[0];
                                }
                            }
                            else
                            {
                                gotIt = false;
                                gotItCounter = 0;
                                // calculates the angle first
                                while (!gotIt)
                                {
                                    xRnd = r.Rndm() * piHalf;
                                    yRnd = r.Rndm();
                                    //angle distribution from Sato (integrated over the spectrum)
                                    // that 1.9 gives the angular distribution from Nesterenok
                                    if (yRnd < exp(-1.9 * (1. - cos(xRnd)))) { gotIt = true; }
                                    gotItCounter++;
                                    if (gotItCounter > 1000) gotIt = true;
                                }
                                theta = thetaBeam + xRnd;
                            }

                            if (usePrecalculatedSato2016Spectrum)
                            {
                                double evalValue;
                                int aa;

                                //calculates array of energies, otherwise takes precalculated
                                if (satoCounter < satoArrayLength)
                                {
                                    xRnd = satoEnergyArray[satoCounter];
                                    satoCounter++;
                                }
                                else
                                {
                                    //#pragma omp for
                                    for (int sp = 0; sp < satoArrayLength; sp++)
                                    {
                                        for (aa = 0; aa < supportAngles; aa++)
                                        {
                                            if (tAngleSplVec.at(aa) > satoAngleArray[sp]) break;
                                        }
                                        if (aa > 0) aa--;

                                        evalValue = tGrSplVec.at(aa)->Eval(r.Rndm());
                                        xRnd = TMath::Power(10, evalValue);

                                        satoEnergyArray[sp] = xRnd;
                                    }

                                    xRnd = satoEnergyArray[0];
                                    satoCounter = 1;
                                }
                            }
                            else
                            {
                                angAllFunction->SetParameters(atmDensity, rigidity, cos(theta));

                                float funcMaxHere = 1.05 * angAllFunction->Eval(150.);

                                gotIt = false;
                                gotItCounter = 0;
                                while (!gotIt)
                                {
                                    if (noHighEnergyRegime) xRnd = TMath::Power(10, r.Rndm() * 9. - 7.7);
                                    else
                                    {
                                        if (noThermalRegime) xRnd = TMath::Power(10, r.Rndm() * 10. - 6.);
                                        else  xRnd = TMath::Power(10, r.Rndm() * 12.3 - 8.);
                                    }
                                    yRnd = r.Rndm() * funcMaxHere;

                                    if (angAllFunction->Eval(xRnd) > yRnd)
                                    {
                                        gotIt = true;
                                    }
                                    gotItCounter++;
                                    if (gotItCounter > 1000) gotIt = true;
                                }
                            }
                        }
                        else
                        {
                            if (false)
                            {
                                gotIt = false;
                                gotItCounter = 0;
                                while (!gotIt)
                                {
                                    if (noHighEnergyRegime) xRnd = TMath::Power(10, r.Rndm() * 9. - 7.7);
                                    else
                                    {
                                        if (noThermalRegime) xRnd = TMath::Power(10, r.Rndm() * 12. - 7.7);
                                        else  xRnd = TMath::Power(10, r.Rndm() * 12.3 - 8.);
                                    }
                                    yRnd = r.Rndm() * 0.14;

                                    if (phiBFunc->Eval(xRnd) > yRnd)
                                    {
                                        gotIt = true;
                                    }
                                    if (useHomogenousSpectrum) gotIt = true;
                                    gotItCounter++;
                                    if (gotItCounter > 1000) gotIt = true;
                                }
                            }
                            //else xRnd = 0.000001;
                        }
                    }
                }

                if (doTheWaterThing)
                {
                    xRnd = r.Gaus(14.1, 0.1);
                }

                if (doTheZredaThing)
                {
                    gotIt = false;
                    while (!gotIt)
                    {
                        xRnd = TMath::Power(10, r.Rndm() * 12. - 7.7);
                        yRnd = r.Rndm();

                        binNo = precalculatedSpectrum->FindBin(xRnd);

                        if (precalculatedSpectrum->GetBinContent(binNo) > yRnd)
                        {
                            gotIt = true;
                        }
                    }
                    //if ((xRnd>0.3)||(r.Rndm()>0.55))	xRnd = getEvaporationEnergy(2.e6, &r);
                    xRnd = getEvaporationEnergy(2.e6, &r);
                }

                if (doSkyEvaporation)
                {
                    if (!doNoSource)
                    {
                        if (doFusion) xRnd = r.Gaus(14.1, 0.1); //13.6?
                        if (doFission) xRnd = getFissionEnergy(&r);
                        if (doAmBe) xRnd = getAmBeEnergy(&r);
                        if (doThermalSource) xRnd = getThermalEnergy(spectrumMaxwellPhiLinMod, &r); ;
                        if (doMonoEnergetic) xRnd = r.Gaus(sourceEnergy, 0.1 * sourceEnergy);
                        if (doModeratedCf) xRnd = getModeratedCfEnergy(&r);
                    }
                }

                energyInitial = xRnd;

                if ((doDetectorBatchRun) || (doDetectorBatchRun2)) energyInitial = paramEnergy;
                if ((doDetectorAngleBatchRun) && (doMonoEnergetic)) energyInitial = r.Gaus(sourceEnergy, 0.1 * sourceEnergy);

                if (!(nEvapoVector.size() > 0)) cosmicSpectrum->Fill(energyInitial);

                startingLayerHere = startingLayer;

                //generate neutron position
                if (useRadialBeam)
                {
                    //generate spatial position
                    r1 = r.Rndm() * beamRadius;
                    phi1 = r.Rndm() * 2. * pi - pi;
                    //set ray initial values
                    x = r1 * cos(-phi1) + beamX;
                    y = r1 * sin(-phi1) + beamY;
                }

                if (useUniformBeam)
                {
                    check = true;
                    while (check)
                    {
                        xt = r.Rndm() * 2. - 1.;
                        yt = r.Rndm() * 2. - 1.;
                        if (useRectShape) check = false;
                        else if ((TMath::Power(xt, 2) + TMath::Power(yt, 2)) < 1) check = false;
                    }
                    x = xt * beamRadius + beamX;
                    y = yt * beamRadius + beamY;
                }

                if (doTheWaterThing) { x = 0; y = 0; }
                if (doTheZredaThing) { x = 0; y = 0; }

                if (godzillaMode)
                {
                    x = xCustomPos; y = yCustomPos;
                    xPosSource = xCustomPos;
                    yPosSource = yCustomPos;
                }

                if (doSkyEvaporation)
                {
                    if (radiusSource > 0)
                    {
                        check = true;
                        while (check)
                        {
                            xt = r.Rndm() * 2. - 1.;
                            yt = r.Rndm() * 2. - 1.;
                            if ((TMath::Power(xt, 2) + TMath::Power(yt, 2)) < 1) check = false;
                        }
                        x = xt * radiusSource + xPosSource;
                        y = yt * radiusSource + yPosSource;
                    }
                    else
                    {
                        if (xSizeSource > 0)
                        {
                            x = xPosSource + (r.Rndm() * 2. - 1.) * 0.5 * xSizeSource;
                        }
                        else
                        {
                            x = xPosSource;
                        }
                        if (ySizeSource > 0)
                        {
                            y = yPosSource + (r.Rndm() * 2. - 1.) * 0.5 * ySizeSource;
                        }
                        else
                        {
                            y = yPosSource;
                        }
                    }
                }

                cs = 0; length = 0;
                timeNtr = 0;

                xStart = x; yStart = y;
                moistureAtInterface = 0;
                xAtInterface = x;
                yAtInterface = y;
                zAtInterface = z0;
                xLastScattered = x; yLastScattered = y; zLastScattered = z0;

                if (doTheWaterThing) z0 = 0;
                else
                {
                    scatteredThisLayer = true;
                    hasbeenEvaporized = true;
                    z0 = r.Rndm() * geometries.at(startingLayer)[5] + geometries.at(startingLayer)[4];
                }

                if (doTheZredaThing)
                {
                    scatteredThisLayer = true;
                    hasbeenEvaporized = true;
                    if (xRnd > 0.3) { z0 = r.Exp(eFoldingDepth); }
                    else z0 = r.Rndm() * geometries.at(startingLayer)[5] + geometries.at(startingLayer)[4];
                }

                if (doSkyEvaporation) { z0 = 0.5 * geometries.at(startingLayer)[5] + geometries.at(startingLayer)[4]; }
                if (doSkyEvaporation) { z0 = zPosSource; scatteredThisLayer = true; }

                z0max = z0;

                //set initial values
                //angles phi (radial), theta (inclination[resp beam direction]) and z0 (starting depth)
                //reverseDir tells direction of theta (and direction to pass geometry)
                //scattered and currentlayer are flags
                //phi = 0.0;
                phi = r.Rndm() * 2. * pi; //theta = r.Rndm()*pi/3.;

                if (useSato2016)
                {
                }
                else
                {
                    if ((energyInitial > 9) && ((!doSkyEvaporation) || (doNoSource)))
                    {
                        // that gives the angular distribution from Nesterenok
                        gotIt = false;
                        gotItCounter = 0;
                        while (!gotIt)
                        {
                            xRnd = r.Rndm() * piHalf;
                            yRnd = r.Rndm();

                            if (yRnd < exp(-2.4 * (1. - cos(xRnd)))) { gotIt = true; }
                            gotItCounter++;
                            if (gotItCounter > 1000) gotIt = true;
                        }
                        theta = thetaBeam + xRnd;
                        //int thetaInt = theta/pi*180.;
                        //r.Rndm()*pi/(energyInitial/3.);
                    }
                    else theta = thetaBeam + TMath::ACos(r.Rndm());
                }

                if (doTheWaterThing) theta = TMath::ACos(r.Rndm());
                if (doTheZredaThing)
                {
                    if (xRnd > 0.3) theta = 0 * piHalf + r.Rndm() * piHalf + piHalf;
                    else  theta = thetaBeam + r.Rndm() * piHalf;
                }

                if (doSkyEvaporation)
                {
                    if ((doBatchRun) || (doBatchRunDensity) || (doBatchRun2D)) theta = TMath::ACos(r.Rndm());
                    else theta = TMath::ACos(2. * r.Rndm() - 1.);

                    scatteredThisLayer = true;

                    //if (doBatchRun) theta = pi/2.;
                    //if (doBatchRun) theta = param/(paramMax-paramMin)*pi/2.;
                    //if (doBatchRun) theta = 0;
                    if (doDetectorBatchRun) { theta = 0; batchVar = paramEnergy; }
                    if (doDetectorBatchRun2) { theta = TMath::ACos(r.Rndm()); batchVar = paramEnergy; }
                    if (doDetectorAngleBatchRun) {theta = param/(paramMax-paramMin)*pi/2.; batchVar = theta;}

                    //theta = scaleFactor/100. *pi;
                }

                // set just one direction
                // phi = 1.5*pi; theta = pi*0.5;

                scattered = false;

                scatterLayer = 0;

                thermalizedLayer = -1;
                thermalized = false;
                timeTrans = 0;
                stillFirstLayer = true;
                enteredDetectorLayer = false;

                if (useVolumeSource)
                {
                    maxgenerationHeight = (geometries.at(startingLayer)[4] - geometries.at(groundLayer)[4]);
                    generationHeight = (geometries.at(startingLayer)[4] + 10. * attenuationLength / rLuft * TMath::Log(1. - r.Rndm() * (1. - TMath::Exp(-(maxgenerationHeight / 10. * rLuft) / (attenuationLength)))));

                    for (int i = 0; i < geometries.size() - 1; i++)
                    {
                        if ((geometries.at(i)[4] < generationHeight) && ((geometries.at(i)[4] + geometries.at(i)[5]) > generationHeight))
                        {
                            startingLayerHere = i;
                            z0 = generationHeight;

                            break;
                        }
                    }
                }

                if (!(nEvapoVector.size() > 0))
                {
                    if (!dontFillunecessaryPlots)
                    {
                        if (scatterings != 0) scatteringNs->Fill(scatterings);
                    }

                    if (true)
                    {
                        if (((nTotal <= 100000) && (n % 500 == 0)) || ((nTotal > 100000) && (n % 2000 == 0)))
                        {
                            nRPerS = double (difftime(diffmean, start) - pauseTime);

                            if (nRPerS > 0)
                            {
                                nRunPerS = (nTotal * 1.) / (nRPerS);
                                nspers = nRunPerS;
                                timeRunRemainHours = ((neutrons - nTotal) / nRunPerS) / 3600;
                                timeRunRemainMinutes = (((neutrons - nTotal) / nRunPerS) - timeRunRemainHours * 3600) / 60;
                                timeRunRemainSeconds = ((neutrons - nTotal) / nRunPerS) - timeRunRemainHours * 3600 - timeRunRemainMinutes * 60;
                            }

                            cout << "\r" << castFloatToString(100. * (nTotal * 1.) / (neutrons * 1.), 6) << " % completed ";
                            if (nRPerS > 0)
                            {
                                cout<< " (" << timeRunRemainHours << ":"; if (timeRunRemainMinutes < 10) cout<<"0"; cout<<timeRunRemainMinutes << ":"; if (timeRunRemainSeconds < 10) cout<<"0"; cout<<timeRunRemainSeconds <<") ["<<(string)castIntToString(nspers)<< " n/s]" ;
                            }
                            cout.flush();
                        }
                    }
                    else
                    {
                        if (n % 100000 == 0) cout << "|" << endl;
                        else { if (n % 10000 == 0) cout << ":"; else if (n % 2000 == 0) cout << "."; }
                    }

                    if ((doDetectorBatchRun) || (doDetectorBatchRun2) || (doDetectorAngleBatchRun))
                    {
                        if (!noGUIMode)
                        {
                            if (n % 10000 == 0) delay(1);
                            if (n % 50000 == 0) delay(1);
                        }
                    }
                    else
                    {
                        if (!noGUIMode)
                        {
                            if (n % 2000 == 0) delay(1);
                            if (n % 10000 == 0) delay(1);
                        }
                    }

                    if ((clearEveryXNeutronsNumber > 0) || (setAutoRefreshRateClearing))
                    {
                        if (((n % clearEveryXNeutronsNumber == 0) && (clearEveryXNeutrons)) || ((n % refreshCycle == 0) && (setAutoRefreshRateClearing)))
                        {
                            if (!noThermalRegime) { densityThermalTrackMap->Reset();  densityThermalTrackMapHighRes->Reset(); densityThermalTrackMapHighRes15x->Reset(); densityMapThermal->Reset(); }

                            densityTrackMap->Reset(); densityIntermediateTrackMap->Reset(); densityFastTrackMap->Reset(); densityAlbedoTrackMap->Reset();  densityEnergyTrackMap->Reset();
                            densityTrackMapHighRes->Reset(); densityIntermediateTrackMapHighRes->Reset(); densityFastTrackMapHighRes->Reset();   densityAlbedoTrackMapHighRes->Reset(); densityHighEnergyTrackMapHighRes->Reset();  densityFastTrackMapHighRes2x->Reset();
                            densityTrackMapHighRes15x->Reset(); densityIntermediateTrackMapHighRes15x->Reset(); densityFastTrackMapHighRes15x->Reset();   densityAlbedoTrackMapHighRes15x->Reset(); densityHighEnergyTrackMapHighRes15x->Reset();  densityEnergyTrackMapHighRes15x->Reset();
                            densityMapIntermediate->Reset(); densityMapFast->Reset(); densityMapAlbedo->Reset(); densityMap->Reset(); densityMapHighEnergy->Reset();  densityEnergyTrackMapHighRes->Reset();

                            if (showDensityTrackMapSide) { densityTrackMapSide->Reset();  densityTrackMapSideAlbedo->Reset(); densityTrackMapSideDetector->Reset(); densityTrackMapSideThermal->Reset(); }
                        }
                    }
                }

                scatterings = 0;

                if (nEvapoVector.size() > 0)
                {
                    x = nEvapoVector.at(0);
                    y = nEvapoVector.at(1);
                    z0 = nEvapoVector.at(2);
                    xLastScattered = x; yLastScattered = y; zLastScattered = z0;

                    theta = nEvapoVector.at(3);
                    phi = nEvapoVector.at(4);
                    startingLayerHere = nEvapoVector.at(5);
                    energy = nEvapoVector.at(6);
                    timeNtr = nEvapoVector.at(7);

                    if (nEvapoVector.at(8) > 0) { hasBeenInSoil = true; }

                    nEvapoVector.erase(nEvapoVector.begin(), nEvapoVector.begin() + 9);

                    energyInitial = energy;

                    xAtInterface = x;
                    yAtInterface = y;
                    zAtInterface = z0;

                    scatteredThisLayer = true;
                    hasbeenEvaporized = true;
                    isEvaporationNeutron = true;

                    z0max = z0;
                    if (z0 >= geometries.at(groundLayer)[4]) { subsurfaceScatterings.push_back(z0); }

                    n--;
                }

                if (theta > piHalf) reverseDir = true;
                else reverseDir = false;

                energyAtInterface = energyInitial;
                energy = energyInitial;

                elementID = -1;

                if (drawSingleNeutronGraphs)
                {
                    neutronPath->SetPoint(0, x / 1000., -z0);
                    scale = (log10(energyInitial) + 7.7) / 12.;
                    rgbValues = getRGBfromHCL(getLinearC(scale, 0, 240), getScaledColorValue(scale, 0.4, 0.91, 0.9), getScaledColorValue(scale, 2.2, 0.8, 0.9));
                    neutronPath->SetLineColor(TColor::GetColor(rgbValues.at(0), rgbValues.at(1), rgbValues.at(2)));

                    graphs.at(graphN)->SetPoint(0, x / 1000., -z0);
                    graphs.at(graphN)->SetLineColor(TColor::GetColor(rgbValues.at(0), rgbValues.at(1), rgbValues.at(2)));
                }
                //fills the number of scatterings for the last neutron (which probably escaped...)


                if (drawSingleNeutronPropagation)
                {
                    float* coordinates = new float[5];
                    coordinates[0] = x;
                    coordinates[1] = y;
                    coordinates[2] = z0;  coordinates[3] = energy; coordinates[4] = 0;
                    if (neutronCoordinates.size() < numberOfFrames + 1) neutronCoordinates.push_back(coordinates);
                }

                thetaOld = 1e9; phiOld = 1e9;
                xAlt = x; yAlt = y; z0Alt = z0; phiAlt = phi; thetaAlt = theta; gAlt = startingLayerHere;
                previousTheta = theta; previousPhi = phi;
                previousX = x; previousY = y; previousZ0 = z0;
                previousSoilX = 0; previousSoilY = 0; previousSoilZ0 = 0; thermalizedX = 0; thermalizedY = 0; thermalizedZ0 = 0;
                previousCosPhi = cos(phi); previousSinPhi = sin(phi); previousTanTheta = tan(theta), previousCosTheta = cos(theta);
                //
                //initialization complete
                //

                //how many layers have been traversed
                glayerCounter = 0;
                newLayer = true;
                neutronTrackCoordinates.clear();
                neutronTrackCoordinates2.clear();
                neutronTrackCoordinatesFullSet.clear();

                absMaterialNo = 0;
                actualMaterial = -1;
                startMaterial = -1;
                actualMaterial2 = -1;
                startMaterial2 = -1;
                actualMaterial3 = -1;
                startMaterial3 = -1;
                neutronAbsorbedbyDetector = false;

                loopNumber++;

                //begin going through the stack
                for (Int_t g = startingLayerHere; (g < geometries.size()) && (g >= 0); )
                {
                    allInelastics.clear(); allInelasticAngulars.clear(); allInelasticElements.clear();
                    allBoden.clear(); allBodenElements.clear(); allPlants.clear(); allPlantsElements.clear(); allMaterials.clear(); allMaterialsElements.clear(); allAbsorptionMaterials.clear(); allAbsorptionMaterialsElements.clear();
                    inelasticEnergyLossVec.clear();
                    scatteredInelastic = false; sAbsorbed = false; scatteredElastic = false; scatteredEvaporating = false; absorbBreak = false;
                    asAll = 0; csIn = 0; inelasticEnergyLoss = 0; //wwRangeIn = 1e9; wwRangeSi = 1e9;
                    weight = 10; element = 15;

                    if (glayerCounter > 30000) { cout << "-!"; break; }

                    if (pausehere)
                    {
                        time(&pause1);
                        //uiM->redrawTopView();
                        //uiM->redrawNeutronMap(-1);
                        while ((pausehere) && (!godzillaMode))  delay(20);
                        while ((pausehere) && (godzillaMode))  delay(2);
                        time(&pause2);
                        pauseTime += difftime(pause2, pause1);
                    }

                    // getting the actual layer
                    currentlayer = geometries.at(g)[7];

                    tanTheta = tan(theta);
                    cosTheta = cos(theta);
                    cosPhi = cos(phi);
                    sinPhi = sin(phi);

                    // temporal coordinates for the actual layer
                    tempz = geometries.at(g)[4];
                    xt = cosPhi * fabs(tanTheta * (tempz - z0)) + x;
                    yt = sinPhi * fabs(tanTheta * (tempz - z0)) + y;

                    tempzEnd = geometries.at(g)[4] + geometries.at(g)[5];
                    xtEnd = cosPhi * fabs(tanTheta * (tempzEnd - z0)) + x;
                    ytEnd = sinPhi * fabs(tanTheta * (tempzEnd - z0)) + y;
                    //}

                    energy1e6 = energy * 1e6;
                    energy1e6Log = TMath::Log10(energy1e6);

                    if (calcNeutronTime)
                    {
                        timeNtr += calcNeutronDiffTime(z0, z0Alt, energy, cosTheta);
                    }

                    thetaOld = theta;
                    phiOld = phi;

                    if ((g == (groundLayer - 1)) && (isEvaporationNeutron) && (hasBeenInSoil) && (scatterings == 0)) {scatteredSurfaceSpectrumHelp->Fill(energy); }

                    if ((z0 < tempzEnd) && (z0 > tempz) && (!scatteredThisLayer) && (!continuedThisLayer)) scatteredThisLayer = true;

                    if (useImage)
                    {
                        if (inputPicSizes[g] > 0)
                        {
                            inputMatrixPixels = inputPicSizes[g];
                            matrixMetricFactor = squareDim * 0.001 / (inputPicSizes[g] * 1.);
                        }
                    }

                    // start coordinate tracking
                    // this has to be moved downwards after material selection
                    if (!noTrackRecording)
                    {
                        if ((((g == detectorLayer) && (!enteredDetectorLayer)) || (trackAllLayers)) && (!continuedThisLayer) && (newLayer))
                        {
                            neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                            //neutronTrackCoordinates.clear();
                            enteredDetectorLayer = true;
                            if ((hasbeenEvaporized) || (scatteredThisLayer))
                            {
                                neutronTrackCoordinates.push_back(x);
                                neutronTrackCoordinates.push_back(y);
                                neutronTrackCoordinates.push_back(z0);
                            }
                            else
                            {
                                if (theta < piHalf)
                                {
                                    neutronTrackCoordinates.push_back(xt);
                                    neutronTrackCoordinates.push_back(yt);
                                    neutronTrackCoordinates.push_back(tempz);
                                }
                                else
                                {
                                    neutronTrackCoordinates.push_back(xtEnd);
                                    neutronTrackCoordinates.push_back(ytEnd);
                                    neutronTrackCoordinates.push_back(tempzEnd);
                                }
                            }

                            neutronTrackCoordinates.push_back(theta);
                            neutronTrackCoordinates.push_back(phi);
                            neutronTrackCoordinates.push_back(energy);
                            neutronTrackCoordinates.push_back(timeNtr);
                            if (isEvaporationNeutron) { neutronTrackCoordinates.push_back(-2); neutronTrackCoordinates.push_back(-2); }
                            else { neutronTrackCoordinates.push_back(-1); neutronTrackCoordinates.push_back(-1); }
                        }
                        if (g != detectorLayer) enteredDetectorLayer = false;

                        //if ((g!=detectorLayer)&&(showDensityTrackMapSide)&&(energy<densitySideTrackingEnergyCutoff)&&(!continuedThisLayer)&&(newLayer))
                        if ((true) && (showDensityTrackMapSide) && (energy < densitySideTrackingEnergyCutoff) && (!continuedThisLayer) && (newLayer))
                        {
                            neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                            //neutronTrackCoordinates2.clear();
                            if ((hasbeenEvaporized) || (scatteredThisLayer))
                            {
                                neutronTrackCoordinates2.push_back(x);
                                neutronTrackCoordinates2.push_back(y);
                                neutronTrackCoordinates2.push_back(z0);
                            }
                            else
                            {
                                if (theta < piHalf)
                                {
                                    neutronTrackCoordinates2.push_back(xt);
                                    neutronTrackCoordinates2.push_back(yt);
                                    neutronTrackCoordinates2.push_back(tempz);
                                }
                                else
                                {
                                    neutronTrackCoordinates2.push_back(xtEnd);
                                    neutronTrackCoordinates2.push_back(ytEnd);
                                    neutronTrackCoordinates2.push_back(tempzEnd);
                                }
                            }
                            neutronTrackCoordinates2.push_back(theta);
                            neutronTrackCoordinates2.push_back(phi);
                            neutronTrackCoordinates2.push_back(energy);
                        }
                    }

                    if (newLayer) newLayer = false;

                    if ((domainCutoff || domainCutoff2) && (energy < 20))
                    {
                        if (domainCutoff)
                        {
                            if (x < -domainCutoffFactor * squareDim) break;
                            if (x > domainCutoffFactor * squareDim) break;
                            if (y < -domainCutoffFactor * squareDim) break;
                            if (y > domainCutoffFactor * squareDim) break;
                        }

                        if (domainCutoff2)
                        {
                            float dcutoffRange = squareDim * 0.5 + domainCutoffMeters * 1000.;

                            //physics criteria
                            if ((x < -dcutoffRange) || (x < -dcutoffRange)) break;
                            if ((x > dcutoffRange) || (x > dcutoffRange)) break;
                            if ((y < -dcutoffRange) || (y < -dcutoffRange)) break;
                            if ((y > dcutoffRange) || (y > dcutoffRange)) break;
                        }
                    }

                    //MATERIAL SELECTION here
                    material = (int)geometries.at(g)[6];

                    if (material == 9)
                    {
                        rWater = 0.99;
                    }

                    // just for snow
                    if (material == 8)
                    {
                        material = 9;
                        rWater = 0.03;
                        if (doBatchRunDensity) rWater = paramDensity;
                    }
                    else
                    {
                        rWater = 0.99;
                    }

                    if (material == 21)
                    {
                        if (doBatchRunDensity) rPlants = paramDensity / 1000.;
                    }

                    if (g >= groundLayer)
                    {
                        if (islandSetup) { if (TMath::Power(x, 2) + TMath::Power(y, 2) > TMath::Power(lakeDiameter, 2))  material = 9; }
                        if (riverSetup) { if ((x > -riverDiameter * 0.5) && (x < riverDiameter * 0.5))  material = 9; }
                        if (lakeSetup) { if (TMath::Power(x, 2) + TMath::Power(y, 2) < TMath::Power(lakeDiameter, 2))  material = 9; }
                        if (coastSetup) { if ((true) && (x > coastPosition))  material = 9; }
                    }

                    detectorLayerOverride = false;
                    detectorOverride = false;
                    inputMatrixValue = -1;

                    if (haveDifferentSoilMoistures)
                    {
                        if ((g >= 0) && (useImage))
                        {
                            if ((scatteredThisLayer) || (continuedThisLayer))
                            {
                                matrixX = (x * 0.001 - matrixStartX) / matrixMetricFactor;
                                matrixY = (y * 0.001 - matrixStartY) / matrixMetricFactor;
                            }
                            else
                            {
                                if (theta > piHalf)
                                {
                                    matrixX = (xtEnd * 0.001 - matrixStartX) / matrixMetricFactor;
                                    matrixY = (ytEnd * 0.001 - matrixStartY) / matrixMetricFactor;
                                }
                                else
                                {
                                    matrixX = (xt * 0.001 - matrixStartX) / matrixMetricFactor;
                                    matrixY = (yt * 0.001 - matrixStartY) / matrixMetricFactor;
                                }

                            }
                            if (matrixX > inputMatrixPixels - 1) matrixX = inputMatrixPixels - 1; if (matrixY > inputMatrixPixels - 1) matrixY = inputMatrixPixels - 1;
                            if (matrixX < 0) matrixX = 0; if (matrixY < 0) matrixY = 0;

                            currentG = g;

                            if ((inputPics2[currentG] == 1) || (inputPics2[currentG] == 2))
                            {
                                materialNotFound = false;

                                inputMatrixValue = (inputPicVector2.at(currentG))(matrixX, matrixY);
                                actualMaterial2 = inputMatrixValue;

                                if ((inputMatrixValue >= 0) && (inputMatrixValue < 251))
                                {
                                    if (inputMatrixValue == 0) rGeneral = 0.001;
                                    if ((inputMatrixValue > 0) && (inputMatrixValue < 201))
                                    {
                                        rGeneral = (inputMatrixValue * 1.) / 100.;
                                    }
                                    if ((inputMatrixValue > 200) && (inputMatrixValue < 251))
                                    {
                                        rGeneral = inputMatrixValue - 200;  materialNotFound = false;
                                        //rGeneral = 13.;
                                    }
                                }
                                else rGeneral = 1.;
                            }
                            else rGeneral = 1.;

                            actualMaterial = -1;

                            // Material definition
                            if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2))
                            {
                                materialNotFound = false;

                                inputMatrixValue = (inputPicVector.at(currentG))(matrixX, matrixY);

                                actualMaterial = inputMatrixValue;

                                //inputMatrixValue = changeInputMatrixValue(inputMatrixValue);

                                //material = (int) geometries.at(g)[6];

                                switch (inputMatrixValue)
                                {
                                case 0:   material = 11; detectorLayerOverride = true; break; // Air
                                case 1:   material = 11; break; //Air
                                case 201: material = 20; soilWaterFrac = 0.1; rBoden = 2.0; break; //concrete wall
                                case 202: material = 20; soilWaterFrac = 0.03; rBoden = 1.4; break; //stones
                                case 203: material = 20; soilWaterFrac = 0.1; rBoden = 0.1 * 1.5; break; //House //.05x
                                case 204: material = 20; soilWaterFrac = 0.1; rBoden = 2.0; break; //concrete street
                                case 205: material = 20; soilWaterFrac = 0.1; rBoden = 0.03 * 2.0; break; //container
                                case 206: material = 26; break; //Aluminium
                                case 207: material = 10; break; //Copy Air (Building)
                                case 208: material = 11; break; //Copy Air with a some humidity
                                case 209: material = 23; rMaterial = 1. * 1.1; break; //cat litter
                                case 210: material = 24; rMaterial = 2.58; break; //Asphalt
                                case 211: material = 20; soilWaterFrac = 0.1; rBoden = 2.0; break;
                                case 212: material = 12; soilWaterFrac = 0.0001; break;
                                case 213: material = 20; soilWaterFrac = 0.2; rBoden = 1.5; break;
                                case 214: material = 19; soilWaterFrac = soilWaterFracVar; rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFracVar * rWater; break;
                                case 215: material = 21; rPlants = 0.030;   break; //Plant gas (Cellulose (25% by weight) plus water))
                                case 216: material = 21; rPlants = 0.020;   break; //Plant gas
                                case 217: material = 21; rPlants = 0.015;   break; //Plant gas
                                case 218: material = 21; rPlants = 0.011;   break; //Plant gas
                                case 219: material = 21; rPlants = 0.008;   break; //Plant gas
                                case 220: material = 21; rPlants = 0.005;  break;  //Plants
                                case 221: material = 21; rPlants = 0.003;   break;  //Tree 0.1 * 0.003?
                                case 222: material = 21; rPlants = 0.002; break;  //Wooden House
                                case 223: material = 21; rPlants = 0.5; break;  //Wood
                                case 224: material = 29; break; //Gd2O3
                                case 225: material = 25; break; //HDPE
                                case 226: material = 26; break; //Aluminum
                                case 227: material = 27; break; //He-3
                                case 228: material = 28; break; //BF3
                                case 229: material = 15; break; //Iron
                                case 230: material = 30; break; //HDPE with approx 3% natural Boron
                                case 231: material = 31; break; //PVC
                                case 232: material = 32; break; //Steel 304L
                                case 233: material = 33; break; //Methane
                                case 234: material = 34; break; //Diesel
                                case 236: material = 36; break; //Graphite
                                case 237: material = 37; break; //Lead
                                case 238: material = 38; break; //TNT
                                case 239: material = 21; break; //Plants
                                case 240: material = 9; rWater = 0.03; break; //Snow (lowest density, new)
                                case 241: material = 9; rWater = 0.1; break; //Snow (new, wet)
                                case 242: material = 9; rWater = 0.3; break; //Snow (old, dry)
                                case 243: material = 9; rWater = 0.5; break; //Snow (old, wet)
                                case 244: material = 9; rWater = 0.85; break; //Snow (lowest density)
                                case 247: material = 7; saltConcentration = 3.5; break;
                                case 248: material = 41; break;
                                case 251: material = 11; detectorLayerOverride = true; break; // Air
                                case 252: material = 11; detectorOverride = true; detectorOverrideID = matrixX * 10000 + matrixY; break; // Air
                                case 254: material = 9; rWater = 0.99; break;  //Water
                                //case 120: if (currentlayer<groundlayer){material = 11;} else {material = 20; soilWaterFrac = 0.3; rBoden = 1.5;} break;
                                case 100: material = 20; soilWaterFrac = soilWaterFracVar; rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFracVar * rWater; break; //Ground
                                default: materialNotFound = true;
                                }

                                if (detectorOverride) { detectorHeight = geometries.at(currentG)[4] + 0.5 * geometries.at(currentG)[5]; }

                                if ((inputMatrixValue > 1) && (inputMatrixValue < 100))
                                {
                                    material = 20; soilWaterFrac = (inputMatrixValue * 1.) / 200.;
                                    rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFrac * rWater;
                                    materialNotFound = false;
                                } //Ground

                                if ((inputMatrixValue > 100) && (inputMatrixValue < 171))
                                {
                                    material = 20; soilWaterFrac = (inputMatrixValue * 1.) / 200.;
                                    rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFrac * rWater;
                                    materialNotFound = false;
                                }

                                if ((inputMatrixValue > 239) && (inputMatrixValue < 240)) // deprecated for now
                                {
                                    material = inputMatrixValue; materialNotFound = false;
                                }

                                if (actualMaterial == 208) { absHumidityAir = 1e-2; relHumidityAir = 0.05; }
                                else { absHumidityAir = absHumidityAirVar; relHumidityAir = relHumidityAirVar; }

                                if (materialNotFound)
                                {
                                    if (warnUndefinedMaterial)  cout << "Material Number not found: " << inputMatrixValue << " at (" << matrixX << "," << matrixY << ") in layer " << currentG + 1 << endl;
                                    material = 11;
                                }
                                if ((warnUndefinedMaterial) && (rGeneral > 1.5)) { cout << "Density Multiplicator: " << inputMatrixValue << " at (" << matrixX << "," << matrixY << ") in layer " << currentG + 1 << endl; }
                            }
                            else
                            {
                                rGeneral = 1.;
                                rMaterial = 1.;
                                saltConcentration = 3.5;
                                if ((int)geometries.at(g)[6] != 8) rWater = 0.99; //not Snow
                                rPlants = 0.003;
                                soilWaterFrac = soilWaterFracVar;
                                absHumidityAir = absHumidityAirVar;
                                relHumidityAir = relHumidityAirVar;
                                rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFrac * rWater;
                                wBoden = soilSolidFracVar * (soilSiFrac * wQuarz + soilAlFrac * wAl2O3) + soilWaterFrac * wWater * soilStrechFactor;
                            }

                            if ((inputPics3[currentG] == 1) || (inputPics3[currentG] == 2))
                            {
                                materialNotFound = false;

                                inputMatrixValue = (inputPicVector3.at(currentG))(matrixX, matrixY);
                                actualMaterial3 = inputMatrixValue;

                                if ((inputMatrixValue > 0) && (inputMatrixValue < 101))
                                {
                                    soilSolidFracVar = 1. - (inputMatrixValue * 1.) / 100.;
                                    //rBoden = soilSolidFracVar*(soilSiFrac*rQuarz+soilAlFrac*rAl2O3)+soilWaterFracVar*rWater;
                                    //wBoden = soilSolidFracVar*(soilSiFrac*wQuarz+soilAlFrac*wAl2O3)+soilWaterFrac*wWater*soilStrechFactor;
                                }
                                else
                                {
                                    soilSolidFracVar = soilSolidFrac;
                                    materialNotFound = true;
                                }
                            }
                        }
                    }

                    if ((doBatchRun2D) && (material == 11))
                    {
                        if (useHumAtmProfile)
                        {
                            if (fabs((geometries.at(g)[4] + geometries.at(g)[5] - geometries.at(groundLayer)[4])) < 200000) absHumidityAir = 1. * airhums[parInt1];
                            else  absHumidityAir = (1. * airhums[parInt1]) * exp(-fabs((geometries.at(g)[4] + geometries.at(g)[5] - (geometries.at(groundLayer)[4] - 200000)) / 2300000.));
                        }
                    }

                    if ((material == 18) || (material == 19) || (material == 20))
                    {
                        rBoden = soilSolidFracVar * (soilSiFrac * rQuarz + soilAlFrac * rAl2O3) + soilWaterFracVar * rWater;
                        wBoden = soilSolidFracVar * (soilSiFrac * wQuarz + soilAlFrac * wAl2O3) + soilWaterFrac * wWater * soilStrechFactor;
                    }
                    if (material == 7) {rSaltWater = 0.99 + saltConcentration / 100. * 0.8; nSalt = (saltConcentration * wWater) / ((100. - saltConcentration) * wSalt); wSaltWater = wWater + nSalt * wSalt; }

                    if (!hasPassedSurface)
                    {
                        if ((g > startingLayer) && (material != 10) && (material != 11))
                        {
                            if (scatteredThisLayer)
                            {
                                energyAtInterface = energy;
                                xAtInterface = x;
                                yAtInterface = y;

                                if ((scatteredThisLayer) || (continuedThisLayer)) zAtInterface = z0;
                                else
                                {
                                    if (theta < piHalf)
                                    {
                                        zAtInterface = tempz;
                                    }
                                    else
                                    {
                                        zAtInterface = tempzEnd;
                                    }
                                }
                                if (material == 18) moistureAtInterface = soilWaterFrac;
                                if (material == 19) moistureAtInterface = soilWaterFrac;
                                if (material == 20) moistureAtInterface = soilWaterFrac;
                                if (material == 41) moistureAtInterface = rCelluloseWaterFrac;
                                if (material == 9) moistureAtInterface = 1;
                                if (material > 20) moistureAtInterface = 0.;

                                hasPassedSurface = true;
                            }
                        }
                    }

                    if ((material != 10) && (material != 11) && (hasPassedSurface) && (!hasBeenInSoil)) hasBeenInSoil = true;

                    //calculate interaction length for absorbtion (wwRangeAbs) and probability for incoherent scattering (probSi <- no spatial information due to small WW probabilities)
                    //materials: water = 9, air (dry) = 10, air (wet) = 11, Quarz = 12, Al2O3 = 13; Soil(wc%) = 20, Plants = 21;
                    //Elements: Hydrogen = 1; Oxygen = 16; Nitrogen = 14; Silicon = 28;
                    //21 = Plants is basically Water with Biomass, which is 4x H2O plus 1xH2O plus 1 Carbon

                    switch (material)
                    {
                    case 7: csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log); csNa = sigmaNaSpline->Eval(energy1e6Log);  csCl = sigmaClSpline->Eval(energy1e6Log);
                        if (csNa < 0) csNa = 0; if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csCl < 0) csCl = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); csNa = calcMeanCS(sigmaNa, energy * 1e6); csCl = calcMeanCS(sigmaCl35, energy * 1e6);
                        csW = csO + 2. * csH + nSalt * (csNa + csCl);
                        cs = csW;
                        break;
                    case 9:  csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csO < 0) csO = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        csW = csO + 2. * csH;
                        cs = csW;
                        break;
                    case 10: csN = sigmaNSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log);  csAr = sigmaArSpline->Eval(energy1e6Log);
                        if (csN < 0) csN = 0; if (csO < 0) csO = 0; if (csAr < 0) csAr = 0;
                         //csN = calcMeanCS(sigmaN, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); csAr = calcMeanCS(sigmaAr, energy * 1e6);
                        cs = 2. * (0.78 * csN + 0.21 * csO) + 0.0093 * csAr;
                        break;
                    case 11: if (energy == lastEnergy11) { csH = csHLast; csO = csOLast; csN = csNLast; csAr = csArLast; }
                        else
                        {
                           csN = sigmaNSpline->Eval(energy1e6Log); csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log); csAr = sigmaArSpline->Eval(energy1e6Log);
                           if (csN < 0) csN = 0; if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csAr < 0) csAr = 0;
                           //csN = calcMeanCS(sigmaN, energy * 1e6);  csH = calcMeanCS(sigmaH, energy * 1e6);  csO = calcMeanCS(sigmaO, energy * 1e6); csAr = calcMeanCS(sigmaAr, energy * 1e6);
                           csHLast = csH; csOLast = csO; csNLast = csN; csArLast = csAr;
                        }
                        csLuft = 2. * (0.78 * csN + 0.21 * csO) + 0.0093 * csAr;
                        csW = csO + 2. * csH;
                        rLuftWater = absHumidityAir / 1e6;
                        cs = csLuft * rLuft / wLuft + csW * rLuftWater / wWater;
                        break;
                    case 12:  csSi = sigmaSiSpline->Eval(energy1e6Log); csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csSi < 0) csSi = 0; if (csO < 0) csO = 0;
                        //csSi = calcMeanCS(sigmaSi, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        cs = 2. * csO + csSi;
                        break;
                    case 13:  csAl = sigmaAlSpline->Eval(energy1e6Log); csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csAl < 0) csAl = 0; if (csO < 0) csO = 0;
                        //csAl = calcMeanCS(sigmaAl, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        cs = 3. * csO + 2. * csAl;
                        break;
                    case 15: csFe = sigmaFeSpline->Eval(energy1e6Log);
                        if (csFe < 0) csFe = 0;
                        //csFe = calcMeanCS(sigmaFe, energy * 1e6);
                        cs = csFe;
                        break;
                    case 18: csAl = calcMeanCS(sigmaAl, energy * 1e6);  csSi = calcMeanCS(sigmaSi, energy * 1e6); csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        csB10 = calcMeanCS(sigmaB10, energy * 1e6); csGd155 = calcMeanCS(sigmaGd155, energy * 1e6); csGd157 = calcMeanCS(sigmaGd157, energy * 1e6);
                        csMn55 = calcMeanCS(sigmaMn55, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6); csFe = calcMeanCS(sigmaFe, energy * 1e6); csN = calcMeanCS(sigmaN, energy * 1e6);  csNa = calcMeanCS(sigmaNa, energy * 1e6);
                        csTi48 = calcMeanCS(sigmaTi48, energy * 1e6); csK = calcMeanCS(sigmaK39, energy * 1e6);
                        csAdd = rBoronInSoil / wBoron * 0.2 * csB10 + rGdInSoil * (0.148 * csGd155 / wGd155 + 0.1565 * csGd157 / wGd157) + rNInSoil / wNitrogen * csN + rFeInSoil / wFe * csFe + rMnInSoil / wMn55 * csMn55 + rCInSoil / wCarbon * csC + rNaInSoil / wNa * csNa + rKInSoil / wKalium * csK + 0.85 * rTiInSoil / wTi48 * csTi48;
                        csSolids = (soilSiFrac * csSi + 2. * soilAlFrac * csAl) * soilSolidFracVar + ((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar;
                        csW = (2. * csH + csO) * soilWaterFrac;
                        cs = csSolids + csW * soilStrechFactor + wBoden / rBoden * 1e-6 * csAdd;
                        break;
                    case 19: //if (energy==lastEnergy19){csH = csHLast; csO = csOLast; csAl = csAlLast; csSi = csSiLast;}
                             //else{
                        csAl = calcMeanCS(sigmaAl, energy * 1e6);  csSi = calcMeanCS(sigmaSi, energy * 1e6); csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); //csB10 = calcMeanCS(sigmaB10,energy*1e6);
                        if (csAl < 0) csAl = 0; if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csSi < 0) csSi = 0;
                        csHLast = csH; csOLast = csO; csAlLast = csAl; csSiLast = csSi; //csB10Last = csB10;
                        //}
                        //cs = (soilSiFrac*csSi+2.*soilAlFrac*csAl)*soilSolidFracVar+csO*(2.25*soilSolidFracVar+soilWaterFrac)+2.*csH*soilWaterFrac;
                        csSolids = (soilSiFrac * csSi + 2. * soilAlFrac * csAl) * soilSolidFracVar + ((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar;
                        csW = (2. * csH + csO) * soilWaterFrac;
                        cs = csSolids + csW * soilStrechFactor;
                        break;
                    case 20: if (energy == lastEnergy20){csH = csHLast; csO = csOLast; csAl = csAlLast; csSi = csSiLast;}
                        else
                        {
                            csAl = sigmaAlSpline->Eval(energy1e6Log);  csSi = sigmaSiSpline->Eval(energy1e6Log); csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log);
                            if (csAl < 0) csAl = 0; if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csSi < 0) csSi = 0;
                            //csAl = calcMeanCS(sigmaAl, energy * 1e6); csSi = calcMeanCS(sigmaSi, energy * 1e6); csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        }
                        csHLast = csH; csOLast = csO; csAlLast = csAl; csSiLast = csSi;
                        //cs = (soilSiFrac*csSi+2.*soilAlFrac*csAl)*soilSolidFracVar+csO*(2.25*soilSolidFracVar+soilWaterFrac)+2.*csH*soilWaterFrac;
                        csSolids = (soilSiFrac * csSi + 2. * soilAlFrac * csAl) * soilSolidFracVar + ((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar;
                        csW = (2. * csH + csO) * soilWaterFrac;
                        //cs = csSolids + csW*soilWaterFrac*soilStrechFactor; break; //error here
                        cs = csSolids + csW * soilStrechFactor;
                        break;
                    case 21: csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log); csN = sigmaNSpline->Eval(energy1e6Log);
                        if (csC < 0) csC = 0; if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csN < 0) csN = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6); csN = calcMeanCS(sigmaN, energy * 1e6);
                        csLuft = 2. * (0.78 * csN + 0.22 * csO);
                        csPlants = csO + 2. * csH + 0.2 * csC;
                        csW = csO + 2. * csH;
                        //rLuftWater = 1.*2.*0.0000177*relHumidityAir;
                        rLuftWater = absHumidityAir / 1e6;
                        cs = csLuft * rLuft / wLuft + csW * rLuftWater / wWater + csPlants * rPlants / wPlants;
                        break;
                    case 23: csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log); csSi = sigmaSiSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csSi < 0) csSi = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); csSi = calcMeanCS(sigmaSi, energy * 1e6);
                        cs = 0.44 * csH + 0.44 * csO + 0.12 * csSi;
                        break;
                    case 24: csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log); csSi = sigmaSiSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csO < 0) csO = 0; if (csSi < 0) csSi = 0; if (csC < 0) csC = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6); csSi = calcMeanCS(sigmaSi, energy * 1e6);
                        cs = 0.14 * csH + 0.5 * csO + 0.11 * csC + 0.25 * csSi;
                        break;
                    case 25:
                        if (energy == lastEnergy25){csH = csHLast; csC = csCLast;}
                        else
                        {
                            csH = sigmaHSpline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log);
                            if (csH < 0) csH = 0; if (csO < 0) csO = 0;
                            //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);
                        }
                        csHLast = csH; csCLast = csC;
                        cs = 2. * csH + csC;
                        break;
                    case 26:  csAl = sigmaAlSpline->Eval(energy1e6Log);
                        if (csAl < 0) csAl = 0;
                        //csAl = calcMeanCS(sigmaAl, energy * 1e6);
                        cs = csAl;  break;
                    case 27: csHe3 = sigmaHe3Spline->Eval(energy1e6Log); //csHe3 = calcMeanCS(sigmaHe3, energy * 1e6);
                        if (csHe3 < 0) csHe3 = 0;
                        cs = csHe3;  break;
                    case 28: csB10 = sigmaB10Spline->Eval(energy1e6Log);  csF = sigmaFSpline->Eval(energy1e6Log);
                        if (csB10 < 0) csB10 = 0; if (csF < 0) csF = 0;
                        //csB10 = calcMeanCS(sigmaB10, energy * 1e6); csF = calcMeanCS(sigmaF, energy * 1e6);
                        cs = csB10 + 3. * csF;
                        break;
                    case 29:  csGd155 = sigmaGd155Spline->Eval(energy1e6Log);  csGd157 = sigmaGd157Spline->Eval(energy1e6Log);  csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csGd155 < 0) csGd155 = 0; if (csGd157 < 0) csGd157 = 0;  if (csO < 0) csO = 0;
                        //csGd155 = calcMeanCS(sigmaGd155, energy * 1e6);  csGd157 = calcMeanCS(sigmaGd157, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        cs = 3. * csO + 2. * 0.148 * csGd155 + 2. * 0.1565 * csGd157;
                        break;
                    case 30: csB10 = sigmaB10Spline->Eval(energy1e6Log); csB11 = sigmaB11Spline->Eval(energy1e6Log);
                        if (csB10 < 0) csB10 = 0; if (csB11 < 0) csB11 = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);
                        csB10 = calcMeanCS(sigmaB10, energy * 1e6); csB11 = calcMeanCS(sigmaB11, energy * 1e6);
                        cs = 2. * csH + csC + 0.04 * (0.2 * csB10 + 0.8 * csB11);
                        break;
                    case 31: csH = sigmaHSpline->Eval(energy1e6Log);  csC = sigmaCSpline->Eval(energy1e6Log);  csCl = sigmaClSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csC < 0) csC = 0; if (csCl < 0) csCl = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);  csCl = calcMeanCS(sigmaCl35, energy * 1e6);
                        cs = 3. * csH + 2. * csC + csCl;
                        break;
                    case 32: csFe = sigmaFeSpline->Eval(energy1e6Log); csCr52 = sigmaCr52Spline->Eval(energy1e6Log); csCr53 = sigmaCr53Spline->Eval(energy1e6Log);
                        csNi = sigmaNi58Spline->Eval(energy1e6Log); csMn55 = sigmaMn55Spline->Eval(energy1e6Log);  csSi = sigmaSiSpline->Eval(energy1e6Log);
                        if (csNi < 0) csNi = 0; if (csMn55 < 0) csMn55 = 0; if (csSi < 0) csSi = 0; if (csFe < 0) csFe = 0;  if (csCr52 < 0) csCr52 = 0; if (csCr53 < 0) csCr53 = 0;
                        //csFe = calcMeanCS(sigmaFe, energy * 1e6); csCr52 = calcMeanCS(sigmaCr52, energy * 1e6); csCr53 = calcMeanCS(sigmaCr53, energy * 1e6);
                        //csNi = calcMeanCS(sigmaNi58, energy * 1e6); csMn55 = calcMeanCS(sigmaMn55, energy * 1e6); csSi = calcMeanCS(sigmaSi, energy * 1e6);
                        cs = 0.68 * csFe + 0.19 * (0.86 * csCr52 + 0.14 * csCr53) + 0.09 * csNi + 0.02 * csSi + 0.02 * csMn55;
                        break;
                    case 33: csH = sigmaHSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csC < 0) csC = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);
                        cs = 4. * csH + csC;
                        break;
                    case 34:  csH = sigmaHSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csC < 0) csC = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);
                        cs = 23. * csH + 12. * csC;
                        break;
                    case 36:  csC = sigmaCSpline->Eval(energy1e6Log);
                        if (csC < 0) csC = 0;
                        //csC = calcMeanCS(sigmaC, energy * 1e6);
                        cs = csAl;  break;
                    case 37: csPb206 = sigmaPb206Spline->Eval(energy1e6Log); csPb207 = sigmaPb207Spline->Eval(energy1e6Log); csPb208 = sigmaPb208Spline->Eval(energy1e6Log);
                        if (csPb206 < 0) csPb206 = 0;  if (csPb207 < 0) csPb207 = 0;  if (csPb208 < 0) csPb208 = 0;
                        //csPb206 = calcMeanCS(sigmaPb206, energy * 1e6); csPb207 = calcMeanCS(sigmaPb207, energy * 1e6); csPb208 = calcMeanCS(sigmaPb208, energy * 1e6);
                        cs = 0.241 * csPb206 + 0.221 * csPb207 + 0.524 * csPb208;
                        break;
                    case 38:  csH = sigmaHSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log);
                        csN = sigmaNSpline->Eval(energy1e6Log); csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csC < 0) csC = 0;  if (csN < 0) csN = 0; if (csO < 0) csO = 0;
                        //csH = calcMeanCS(sigmaH, energy * 1e6); csC = calcMeanCS(sigmaC, energy * 1e6);
                        //csN = calcMeanCS(sigmaN, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        cs = 0.35 * csC + 0.238 * csH + 0.144 * csN + 0.286 * csO;
                        break;
                    case 41: csH = sigmaHSpline->Eval(energy1e6Log); csC = sigmaCSpline->Eval(energy1e6Log); csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csH < 0) csH = 0; if (csC < 0) csC = 0; if (csO < 0) csO = 0;
                        //csH = calcMeanCS(sigmaH,energy*1e6); csC = calcMeanCS(sigmaC,energy*1e6); csO = calcMeanCS(sigmaO,energy*1e6);
                        csW = csO + 2.*csH;
                        cs = ((6. * csC + 10. * csH + 5. * csO) * nCellulose + csW * nCelluloseWater) / nCelluloseMix;
                        break;
                    case 100: csN = sigmaNSpline->Eval(energy1e6Log); csO = sigmaOSpline->Eval(energy1e6Log);
                        if (csN < 0) csN = 0; if (csO < 0) csO = 0;
                        // csN = calcMeanCS(sigmaN, energy * 1e6); csO = calcMeanCS(sigmaO, energy * 1e6);
                        cs = 2. * (0.78 * csN + 0.22 * csO);
                        break;
                    }

                    switch (material)
                    {
                    case 7: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e67) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3 * calcMeanCS(absorbMt209O, energy * 1e6);
                        if (energy1e6 < 1e6) {asCl = absorbCl35Spline->Eval(energy1e6Log) + absorbMt103Cl35Spline->Eval(energy1e6Log) + absorbMt107Cl35Spline->Eval(energy1e6Log); } else {asCl = calcMeanCS(absorbCl35, energy1e6);} if (energy1e6 > 2e7) {asCl += calcMeanCS(absorbMt5Cl35, energy1e6);}
                        if (energy1e6 < 1e6) {asNa = absorbNaSpline->Eval(energy1e6Log);} else {asNa = calcMeanCS(absorbNa, energy1e6);} if (energy1e6 > 3e6) { asNa += calcMeanCS(absorbMt5Na, energy * 1e6) + calcMeanCS(absorbMt103Na, energy * 1e6) + calcMeanCS(absorbMt107Na, energy1e6);}
                        if (asNa < 0) asNa = 0; if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asCl < 0) asCl = 0;
                        //asNa = calcMeanCS(absorbNa, energy * 1e6) + calcMeanCS(absorbMt5Na, energy * 1e6) + calcMeanCS(absorbMt103Na, energy * 1e6) + calcMeanCS(absorbMt107Na, energy * 1e6);
                        //asCl = calcMeanCS(absorbCl35, energy * 1e6) + calcMeanCS(absorbMt5Cl35, energy * 1e6) + calcMeanCS(absorbMt103Cl35, energy * 1e6) + calcMeanCS(absorbMt107Cl35, energy * 1e6);
                        asWater = 2. * asH + asO;
                        asAll = asWater + nSalt * (asNa + asCl);
                        break;
                    case 9:  asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (asH < 0) asH = 0; if (asO < 0) asO = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3 * calcMeanCS(absorbMt209O, energy * 1e6);
                        asWater = 2. * asH + asO;
                        asAll = asWater;
                        break;
                    case 10: if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asN = absorbNSpline->Eval(energy1e6Log) + absorbNbSpline->Eval(energy1e6Log);} else {asN = calcMeanCS(absorbN, energy1e6) + calcMeanCS(absorbNb, energy1e6);} if (energy1e6 > 1e6) {asN += calcMeanCS(absorbMt5N, energy1e6) + calcMeanCS(absorbMt104N,energy1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N,energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N,energy1e6);}
                        if (energy1e6 < 1e5) {asAr = absorbArSpline->Eval(energy1e6Log);} else {asAr = calcMeanCS(absorbAr, energy1e6);} if (energy1e6 > 3e6) {asAr += calcMeanCS(absorbMt5Ar, energy1e6) + calcMeanCS(absorbMt103Ar, energy1e6) + calcMeanCS(absorbMt107Ar, energy1e6) + 3. * calcMeanCS(absorbMt209Ar,energy1e6);}
                        //asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy * 1e6) + calcMeanCS(absorbMt107N, energy * 1e6) + calcMeanCS(absorbMt108N, energy * 1e6) + 3. * calcMeanCS(absorbMt209N, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //asAr = calcMeanCS(absorbAr, energy * 1e6) + calcMeanCS(absorbMt5Ar, energy * 1e6) + calcMeanCS(absorbMt103Ar, energy * 1e6) + calcMeanCS(absorbMt107Ar, energy * 1e6) + 3. * calcMeanCS(absorbMt209Ar, energy * 1e6);
                        if (asN < 0) asN = 0; if (asO < 0) asO = 0; if (asAr < 0) asAr = 0;
                        asLuft = 0.78 * 2. * asN + 0.21 * 2. * asO + 0.0093 * asAr;
                        asAll = asLuft;
                        break;
                        /*case 11: asLuft = 0.78*2.*(calcMeanCS(absorbN,energy*1e6)+calcMeanCS(absorbMt5N,energy*1e6)+calcMeanCS(absorbMt103N,energy*1e6)+calcMeanCS(absorbMt107N,energy*1e6)+3.*calcMeanCS(absorbMt209N,energy*1e6)+calcMeanCS(absorbNb,energy*1e6) )
                                 + 0.21*2*(calcMeanCS(absorbO,energy*1e6)+calcMeanCS(absorbMt5O,energy*1e6)+calcMeanCS(absorbMt103O,energy*1e6)+calcMeanCS(absorbMt107O,energy*1e6)+3.*calcMeanCS(absorbMt209O,energy*1e6))
                                 + 0.01*(calcMeanCS(absorbAr,energy*1e6)+calcMeanCS(absorbMt5Ar,energy*1e6)+calcMeanCS(absorbMt103Ar,energy*1e6)+calcMeanCS(absorbMt107Ar,energy*1e6)+3.*calcMeanCS(absorbMt209Ar,energy*1e6));
                                 asWater = 2.* (calcMeanCS(absorbH,energy*1e6)+calcMeanCS(absorbMt5H,energy*1e6)+3.*calcMeanCS(absorbMt209H,energy*1e6)) + calcMeanCS(absorbO,energy*1e6)+calcMeanCS(absorbMt5O,energy*1e6)+calcMeanCS(absorbMt103O,energy*1e6)+calcMeanCS(absorbMt107O,energy*1e6)+3.*calcMeanCS(absorbMt209O,energy*1e6); asAll = asLuft;  break;
                                 */
                    case 11: if (energy == lastEnergy11) { asH = asHLast; asO = asOLast; asN = asNLast; asAr = asArLast; }
                           else
                           {
                            //if (energy > 0.1) cout<<energy<<" "<<endl;
                            asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                            if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) {asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                            if (energy1e6 < 1e6) {asN = absorbNSpline->Eval(energy1e6Log) + absorbNbSpline->Eval(energy1e6Log);} else {asN = calcMeanCS(absorbN, energy1e6) + calcMeanCS(absorbNb, energy1e6);} if (energy1e6 > 1e6) {asN += calcMeanCS(absorbMt5N, energy1e6) + calcMeanCS(absorbMt104N,energy1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N,energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N,energy1e6);}
                            if (energy1e6 < 1e5) {asAr = absorbArSpline->Eval(energy1e6Log);} else {asAr = calcMeanCS(absorbAr, energy1e6);} if (energy1e6 > 3e6) {asAr += calcMeanCS(absorbMt5Ar, energy1e6) + calcMeanCS(absorbMt103Ar, energy1e6) + calcMeanCS(absorbMt107Ar, energy1e6) + 3. * calcMeanCS(absorbMt209Ar,energy1e6);}
                            //if (energy > 0.1) cout<<asH<<" "<<asO<<" "<<asN<<" "<<asAr<<endl;
                            if (asN < 0) asN = 0; if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asAr < 0) asAr = 0;
                            //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                            //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                            //asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy * 1e6) + calcMeanCS(absorbMt107N, energy * 1e6) + calcMeanCS(absorbMt108N, energy * 1e6) + 3. * calcMeanCS(absorbMt209N, energy * 1e6);
                            //asAr = calcMeanCS(absorbAr, energy * 1e6) + calcMeanCS(absorbMt5Ar, energy * 1e6) + calcMeanCS(absorbMt103Ar, energy * 1e6) + calcMeanCS(absorbMt107Ar, energy * 1e6) + 3. * calcMeanCS(absorbMt209Ar, energy * 1e6);
                            //if (energy > 0.1) cout<<asH<<" "<<asO<<" "<<asN<<" "<<asAr<<endl;
                            asHLast = asH; asOLast = asO; asNLast = asN; asArLast = asAr;
                           }
                           asLuft = 0.78 * 2. * asN + 0.21 * 2. * asO + 0.0093 * asAr;
                           asWater = asO + 2. * asH;
                           asAll = asLuft * rLuft / wLuft + asWater * rLuftWater / wWater;
                           break;
                    case 12: if (energy1e6 < 1e6) {asSi = absorbSiSpline->Eval(energy1e6Log);} else {asSi = calcMeanCS(absorbSi, energy1e6);} if (energy1e6 > 2e6) { asSi += calcMeanCS(absorbMt5Si, energy1e6) + calcMeanCS(absorbMt103Si, energy1e6) + calcMeanCS(absorbMt107Si, energy1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (asSi < 0) asSi = 0; if (asO < 0) asO = 0;
                        //asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        asQuarz = 2. * asO + asSi;
                        asAll = asQuarz;
                        break;
                    case 13:  if (energy1e6 < 1e6) {asAl = absorbAlSpline->Eval(energy1e6Log);} else {asAl = calcMeanCS(absorbAl, energy1e6);} if (energy1e6 > 2e6) { asAl += calcMeanCS(absorbMt5Al, energy1e6) + calcMeanCS(absorbMt103Al, energy1e6) + calcMeanCS(absorbMt107Al, energy1e6) + 3. * calcMeanCS(absorbMt209Al, energy1e6); }
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        //asAl = calcMeanCS(absorbAl, energy * 1e6) + calcMeanCS(absorbMt5Al, energy * 1e6) + calcMeanCS(absorbMt103Al, energy * 1e6) + calcMeanCS(absorbMt107Al, energy * 1e6) + 3. * calcMeanCS(absorbMt209Al, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        if (asAl < 0) asAl = 0; if (asO < 0) asO = 0;
                        asAl2O3 = 3. * asO + 2. * asAl;
                        asAll = asAl2O3;
                        break;
                    case 15: if (energy1e6 < 1e6) {asFe = absorbFeSpline->Eval(energy1e6Log);} else {asFe = calcMeanCS(absorbFe, energy1e6);} if (energy1e6 > 1e6) { asFe += calcMeanCS(absorbMt5Fe, energy1e6) + calcMeanCS(absorbMt103Fe, energy1e6) + calcMeanCS(absorbMt107Fe, energy1e6) + 3. * calcMeanCS(absorbMt209Fe, energy1e6);}
                        //asFe = 1. * (calcMeanCS(absorbFe, energy * 1e6) + calcMeanCS(absorbMt5Fe, energy * 1e6) + calcMeanCS(absorbMt103Fe, energy * 1e6) + calcMeanCS(absorbMt107Fe, energy * 1e6) + 3. * calcMeanCS(absorbMt209Fe, energy * 1e6));
                        if (asFe < 0) asFe = 0;
                        asAll = asFe;
                        break;
                    case 18: asAl = calcMeanCS(absorbAl, energy * 1e6) + calcMeanCS(absorbMt5Al, energy * 1e6) + calcMeanCS(absorbMt103Al, energy * 1e6) + calcMeanCS(absorbMt107Al, energy * 1e6) + 3. * calcMeanCS(absorbMt209Al, energy1e6);
                        asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);
                        asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        asB10 = calcMeanCS(absorbB10, energy * 1e6);
                        asMn55 = calcMeanCS(absorbMn55, energy * 1e6) + calcMeanCS(absorbMt5Mn55, energy * 1e6) + calcMeanCS(absorbMt103Mn55, energy * 1e6) + calcMeanCS(absorbMt107Mn55, energy * 1e6);
                        asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N, energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N, energy1e6);
                        asFe = 1. * (calcMeanCS(absorbFe, energy * 1e6) + calcMeanCS(absorbMt5Fe, energy * 1e6) + calcMeanCS(absorbMt103Fe, energy * 1e6) + calcMeanCS(absorbMt107Fe, energy * 1e6) + 3. * calcMeanCS(absorbMt209Fe, energy * 1e6));
                        asGd155 = calcMeanCS(absorbGd155, energy * 1e6); asGd157 = calcMeanCS(absorbGd157, energy * 1e6);
                        asNa = calcMeanCS(absorbNa, energy * 1e6) + calcMeanCS(absorbMt5Na, energy * 1e6) + calcMeanCS(absorbMt103Na, energy * 1e6) + calcMeanCS(absorbMt107Na, energy * 1e6);
                        asK = calcMeanCS(absorbK39, energy * 1e6) + calcMeanCS(absorbMt103K39, energy * 1e6) + calcMeanCS(absorbMt107K39, energy * 1e6);
                        asTi48 = calcMeanCS(absorbTi48, energy * 1e6) + calcMeanCS(absorbMt5Ti48, energy * 1e6) + calcMeanCS(absorbMt103Ti48, energy * 1e6) + calcMeanCS(absorbMt104Ti48, energy * 1e6) + calcMeanCS(absorbMt107Ti48, energy * 1e6);

                        asAdd = rBoronInSoil / wBoron * 0.2 * asB10 + rGdInSoil * (0.148 * asGd155 / wGd155 + 0.1565 * asGd157 / wGd157) + rNInSoil / wNitrogen * asN + rFeInSoil / wFe * asFe + rMnInSoil / wMn55 * asMn55 + rCInSoil / wCarbon * asC + rNaInSoil / wNa * asNa + rKInSoil / wKalium * asK + 0.85 * rTiInSoil / wTi48 * asTi48;
                        asBoden = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * (2.25 * soilSolidFracVar + soilWaterFrac) + 2. * asH * soilWaterFrac; //asAll = asBoden;
                        asWater = asH + 2. * asO;
                        asSolids = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * ((soilSiFrac * 2. + 3. * soilAlFrac) * asO * soilSolidFracVar);
                        asAll = asSolids + asWater * soilWaterFrac * soilStrechFactor + wBoden / rBoden * 1e-6 * asAdd;
                        break;
                    case 19: if (energy == lastEnergy19) {asH = asHLast; asO = asOLast; asAl = asAlLast; asSi = asSiLast; asB10Last = asB10; }
                        else
                        {
                        asAl = calcMeanCS(absorbAl, energy * 1e6) + calcMeanCS(absorbMt5Al, energy * 1e6) + calcMeanCS(absorbMt103Al, energy * 1e6) + calcMeanCS(absorbMt107Al, energy * 1e6) + 3. * calcMeanCS(absorbMt209Al, energy * 1e6);
                        asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        asB10 = calcMeanCS(absorbB10, energy * 1e6);
                        asHLast = asH; asOLast = asO; asAlLast = asAl; asSiLast = asSi; asB10Last = asB10;
                        }
                           asBoden = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * (2.25 * soilSolidFracVar + soilWaterFrac) + 2. * asH * soilWaterFrac; //asAll = asBoden;
                           asWater = asH + 2. * asO;
                           asSolids = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * ((soilSiFrac * 2. + 3. * soilAlFrac) * asO * soilSolidFracVar);
                           asAll = asSolids + asWater * soilWaterFrac * soilStrechFactor + asB10 * rBoronInSoil * 1e-6 * wBoden / rBoden / wBoron;
                           break;
                    case 20: if (energy == lastEnergy20) {asH = asHLast; asO = asOLast; asAl = asAlLast; asSi = asSiLast; }
                        else
                        {
                        //  if (energy > 0.0001) cout<<energy<<" ";
                        asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asAl = absorbAlSpline->Eval(energy1e6Log);} else {asAl = calcMeanCS(absorbAl, energy1e6);} if (energy1e6 > 2e6) { asAl += calcMeanCS(absorbMt5Al, energy1e6) + calcMeanCS(absorbMt103Al, energy1e6) + calcMeanCS(absorbMt107Al, energy1e6) + 3. * calcMeanCS(absorbMt209Al, energy1e6); }
                        if (energy1e6 < 1e6) {asSi = absorbSiSpline->Eval(energy1e6Log);} else {asSi = calcMeanCS(absorbSi, energy1e6);} if (energy1e6 > 2e6) { asSi += calcMeanCS(absorbMt5Si, energy1e6) + calcMeanCS(absorbMt103Si, energy1e6) + calcMeanCS(absorbMt107Si, energy1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);}
                        //asBoden = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * (2.25 * soilSolidFracVar + soilWaterFrac) + 2. * asH * soilWaterFrac;
                        if (asSi < 0) asSi = 0; if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asAl < 0) asAl = 0;
                        //if (energy > 0.0001) cout<<asBoden<<" ";
                        //if (asBoden < 0) cout<<energy<<" "<<asH<<" "<<asO<<" "<<asAl<<" "<<asSi<<" "<<endl;
                        //asAl = calcMeanCS(absorbAl, energy * 1e6) + calcMeanCS(absorbMt5Al, energy * 1e6) + calcMeanCS(absorbMt103Al, energy * 1e6) + calcMeanCS(absorbMt107Al, energy * 1e6) + 3. * calcMeanCS(absorbMt209Al, energy * 1e6);
                        //asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //if (asBoden < 0) cout<<energy<<" "<<asO<<" "<<asAl<<" "<<asSi<<" "<<endl;
                        //asBoden = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * (2.25 * soilSolidFracVar + soilWaterFrac) + 2. * asH * soilWaterFrac;
                        //if (energy > 0.0001) cout<<asBoden<<endl;
                        asHLast = asH; asOLast = asO; asAlLast = asAl; asSiLast = asSi;
                        }
                           asBoden = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * (2.25 * soilSolidFracVar + soilWaterFrac) + 2. * asH * soilWaterFrac; //asAll = asBoden;
                           asWater = asH + 2. * asO;
                           asSolids = (soilSiFrac * asSi + 2. * soilAlFrac * asAl) * soilSolidFracVar + asO * ((soilSiFrac * 2. + 3. * soilAlFrac) * asO * soilSolidFracVar);
                           asAll = asSolids + asWater * soilWaterFrac * soilStrechFactor;
                           break;
                    case 21: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asN = absorbNSpline->Eval(energy1e6Log) + absorbNbSpline->Eval(energy1e6Log);} else {asN = calcMeanCS(absorbN, energy1e6) + calcMeanCS(absorbNb, energy1e6);} if (energy1e6 > 1e6) {asN += calcMeanCS(absorbMt5N, energy1e6) + calcMeanCS(absorbMt104N,energy1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N,energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N,energy1e6);}
                        if (energy1e6 < 1e5) {asAr = absorbArSpline->Eval(energy1e6Log);} else {asAr = calcMeanCS(absorbAr, energy1e6);} if (energy1e6 > 3e6) {asAr += calcMeanCS(absorbMt5Ar, energy1e6) + calcMeanCS(absorbMt103Ar, energy1e6) + calcMeanCS(absorbMt107Ar, energy1e6) + 3. * calcMeanCS(absorbMt209Ar,energy1e6);}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asN < 0) asN = 0; if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asAr < 0) asAr = 0; if (asC < 0) asC = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        //asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy * 1e6) + calcMeanCS(absorbMt107N, energy * 1e6) + calcMeanCS(absorbMt108N, energy * 1e6) + 3. * calcMeanCS(absorbMt209N, energy * 1e6);
                        //asAr = calcMeanCS(absorbAr, energy * 1e6) + calcMeanCS(absorbMt5Ar, energy * 1e6) + calcMeanCS(absorbMt103Ar, energy * 1e6) + calcMeanCS(absorbMt107Ar, energy * 1e6) + 3. * calcMeanCS(absorbMt209Ar, energy * 1e6);
                        asLuft = 0.78 * 2. * asN + 0.21 * 2. * asO + 0.0093 * asAr;
                        asWater = asO + 2. * asH;
                        asPlants = asO + 2. * asH + 0.2 * asC;
                        asAll = asLuft * rLuft / wLuft + asWater * rLuftWater / wWater + asPlants * rPlants / wPlants;
                        break;
                    case 23: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asSi = absorbSiSpline->Eval(energy1e6Log);} else {asSi = calcMeanCS(absorbSi, energy1e6);} if (energy1e6 > 2e6) { asSi += calcMeanCS(absorbMt5Si, energy1e6) + calcMeanCS(absorbMt103Si, energy1e6) + calcMeanCS(absorbMt107Si, energy1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);}
                        if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asSi < 0) asSi = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        asAll = 0.44 * asO + 0.44 * asH + 0.12 * asSi;
                        break;
                    case 24: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asSi = absorbSiSpline->Eval(energy1e6Log);} else {asSi = calcMeanCS(absorbSi, energy1e6);} if (energy1e6 > 2e6) { asSi += calcMeanCS(absorbMt5Si, energy1e6) + calcMeanCS(absorbMt103Si, energy1e6) + calcMeanCS(absorbMt107Si, energy1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asC < 0) asC = 0; if (asSi < 0) asSi = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        //asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        asAll = 0.5 * asO + 0.14 * asH + 0.11 * asC + 0.25 * asSi;
                        break;
                    case 25: if (energy == lastEnergy25) { asH = asHLast; asC = asCLast; }
                        else
                        {
                            //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                            //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                            asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);} // absorbMt5HSpline->Eval(energy1e6Log) + 3. * absorbMt209HSpline->Eval(energy1e6Log);
                            if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);} //absorbMt5OSpline->Eval(energy1e6Log) + absorbMt103OSpline->Eval(energy1e6Log) + absorbMt105OSpline->Eval(energy1e6Log) + absorbMt107OSpline->Eval(energy1e6Log) + 3. * absorbMt209OSpline->Eval(energy1e6Log);
                            if (asH < 0) asH = 0; if (asO < 0) asO = 0;
                        }
                        asHLast = asH; asCLast = asC;
                        asAll = 2. * asH + 1. * asC;
                        break;
                    case 26: if (energy1e6 < 1e6) {asAl = absorbAlSpline->Eval(energy1e6Log);} else {asAl = calcMeanCS(absorbAl, energy1e6);} if (energy1e6 > 2e6) { asAl += calcMeanCS(absorbMt5Al, energy1e6) + calcMeanCS(absorbMt103Al, energy1e6) + calcMeanCS(absorbMt107Al, energy1e6) + 3. * calcMeanCS(absorbMt209Al, energy1e6); }
                         if (asAl < 0) asAl = 0;
                        //asAl = calcMeanCS(absorbAl, energy * 1e6) + calcMeanCS(absorbMt5Al, energy * 1e6) + calcMeanCS(absorbMt103Al, energy * 1e6) + calcMeanCS(absorbMt107Al, energy * 1e6) + 3. * calcMeanCS(absorbMt209Al, energy * 1e6);
                        asAll = asAl;
                        break;
                    case 27: if (energy1e6 < 1e6) {asHe3 = absorbMt103He3Spline->Eval(energy1e6Log);} else {asHe3 = calcMeanCS(absorbMt103He3, energy1e6);} if (energy1e6 > 4e6) { asHe3 += calcMeanCS(absorbMt104He3, energy1e6) ; }
                         if (asHe3 < 0) asHe3 = 0;
                        //asHe3 = + calcMeanCS(absorbMt103He3, energy * 1e6);
                        asAll = asHe3;
                        break;
                    case 28: if (energy1e6 < 1e6) {asF = absorbFSpline->Eval(energy1e6Log);} else {asF = calcMeanCS(absorbF, energy1e6);} if (energy1e6 > 3e6) { asF += calcMeanCS(absorbMt103F, energy1e6) + calcMeanCS(absorbMt107F, energy1e6); }
                        if (energy1e6 < 1e6) {asB10 = absorbB10Spline->Eval(energy1e6Log);} else {asB10 = 0;}
                        if (asF < 0) asF = 0;  if (asB10 < 0) asB10 = 0;
                        //asF = calcMeanCS(absorbF, energy * 1e6) +
                        //asB10 = calcMeanCS(absorbB10, energy * 1e6);
                        asAll = asB10 + 3. * asF;
                        break;
                    case 29: if (energy1e6 < 1e6) {asGd155 = absorbGd155Spline->Eval(energy1e6Log);} else {asGd155 = 0;}
                        if (energy1e6 < 1e6) {asGd157 = absorbGd157Spline->Eval(energy1e6Log);} else {asGd157 = 0;}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (asO < 0) asO = 0;  if (asGd155 < 0) asGd155 = 0; if (asGd157 < 0) asGd157 = 0;
                        //asGd155 = calcMeanCS(absorbGd155, energy * 1e6); asGd157 = calcMeanCS(absorbGd157, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        asAll = 3. * asO + 2. * 0.148 * asGd155 + 2. * 0.1565 * asGd157;
                        break;
                    case 30: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asB10 = absorbB10Spline->Eval(energy1e6Log);} else {asB10 = 0;}
                        if (energy1e6 < 1e6) {asB11 = absorbB10Spline->Eval(energy1e6Log);} else {asB11 = 0;}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asH < 0) asH = 0;  if (asB10 < 0) asB10 = 0; if (asC < 0) asC = 0;  if (asB11 < 0) asB11 = 0;
                        //asB10 = calcMeanCS(absorbB10, energy * 1e6); asB11 = calcMeanCS(absorbB11, energy * 1e6);
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        asAll = 2. * asH + 1. * asC + 0.04 * (0.2 * asB10 + 0.8 * asB11);
                        break;
                    case 31: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (energy1e6 < 1e6) {asCl = absorbCl35Spline->Eval(energy1e6Log) + absorbMt103Cl35Spline->Eval(energy1e6Log) + absorbMt107Cl35Spline->Eval(energy1e6Log); } else {asCl = calcMeanCS(absorbCl35, energy1e6);} if (energy1e6 > 2e7) {asCl += calcMeanCS(absorbMt5Cl35, energy1e6);}
                        if (asH < 0) asH = 0;  if (asC < 0) asC = 0; if (asCl < 0) asCl = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        //asCl = calcMeanCS(absorbCl35, energy * 1e6) + calcMeanCS(absorbMt5Cl35, energy * 1e6) + calcMeanCS(absorbMt103Cl35, energy * 1e6) + calcMeanCS(absorbMt107Cl35, energy * 1e6);
                        asAll = 3. * asH + 2. * asC + asCl;
                        break;
                    case 32: if (energy1e6 < 1e6) {asFe = absorbFeSpline->Eval(energy1e6Log);} else {asFe = calcMeanCS(absorbFe, energy1e6);} if (energy1e6 > 1e6) { asFe += calcMeanCS(absorbMt5Fe, energy1e6) + calcMeanCS(absorbMt103Fe, energy1e6) + calcMeanCS(absorbMt107Fe, energy1e6) + 3. * calcMeanCS(absorbMt209Fe, energy1e6);}
                        if (energy1e6 < 1e6) {asNi = absorbNi58Spline->Eval(energy1e6Log);} else {asNi = calcMeanCS(absorbNi58, energy1e6);} if (energy1e6 > 1e6) { asNi += calcMeanCS(absorbMt5Ni58, energy1e6) + calcMeanCS(absorbMt103Ni58, energy1e6) + calcMeanCS(absorbMt107Ni58, energy1e6);}
                        if (energy1e6 < 1e6) {asCr52 = absorbCr52Spline->Eval(energy1e6Log);} else {asCr52 = calcMeanCS(absorbCr52, energy1e6);} if (energy1e6 > 1e6) { asCr52 += calcMeanCS(absorbMt5Cr52, energy1e6) + calcMeanCS(absorbMt103Cr52, energy1e6) + calcMeanCS(absorbMt107Cr52, energy1e6);}
                        if (energy1e6 < 1e6) {asCr53 = absorbCr53Spline->Eval(energy1e6Log);} else {asCr53 = calcMeanCS(absorbCr53, energy1e6);} if (energy1e6 > 1e6) { asCr53 += calcMeanCS(absorbMt5Cr53, energy1e6) + calcMeanCS(absorbMt103Cr53, energy1e6) + calcMeanCS(absorbMt107Cr53, energy1e6);}
                        if (energy1e6 < 1e6) {asMn55 = absorbMn55Spline->Eval(energy1e6Log);} else {asMn55 = calcMeanCS(absorbMn55, energy1e6);} if (energy1e6 > 3e6) { asMn55 += calcMeanCS(absorbMt5Mn55, energy1e6) + calcMeanCS(absorbMt103Mn55, energy1e6) + calcMeanCS(absorbMt107Mn55, energy1e6);}
                        if (energy1e6 < 1e6) {asSi = absorbSiSpline->Eval(energy1e6Log);} else {asSi = calcMeanCS(absorbSi, energy1e6);} if (energy1e6 > 2e6) { asSi += calcMeanCS(absorbMt5Si, energy1e6) + calcMeanCS(absorbMt103Si, energy1e6) + calcMeanCS(absorbMt107Si, energy1e6) + 3. * calcMeanCS(absorbMt209Si, energy1e6);}
                        if (asFe < 0) asFe = 0; if (asNi < 0) asNi = 0; if (asCr52 < 0) asCr52 = 0; if (asCr53 < 0) asCr53 = 0; if (asMn55 < 0) asMn55 = 0; if (asSi < 0) asSi = 0;
                        //asFe = calcMeanCS(absorbFe, energy * 1e6) + calcMeanCS(absorbMt5Fe, energy * 1e6) + calcMeanCS(absorbMt103Fe, energy * 1e6) + calcMeanCS(absorbMt107Fe, energy * 1e6);
                        //asNi = calcMeanCS(absorbNi58, energy * 1e6) + calcMeanCS(absorbMt5Ni58, energy * 1e6) + calcMeanCS(absorbMt103Ni58, energy * 1e6) + calcMeanCS(absorbMt107Ni58, energy * 1e6);
                        //asCr52 = calcMeanCS(absorbCr52, energy * 1e6) + calcMeanCS(absorbMt5Cr52, energy * 1e6) + calcMeanCS(absorbMt103Cr52, energy * 1e6) + calcMeanCS(absorbMt107Cr52, energy * 1e6);
                        //asCr53 = calcMeanCS(absorbCr53, energy * 1e6) + calcMeanCS(absorbMt5Cr53, energy * 1e6) + calcMeanCS(absorbMt103Cr53, energy * 1e6) + calcMeanCS(absorbMt107Cr53, energy * 1e6);
                        //asMn55 = calcMeanCS(absorbMn55, energy * 1e6) + calcMeanCS(absorbMt5Mn55, energy * 1e6) + calcMeanCS(absorbMt103Mn55, energy * 1e6) + calcMeanCS(absorbMt107Mn55, energy * 1e6);
                        //asSi = calcMeanCS(absorbSi, energy * 1e6) + calcMeanCS(absorbMt5Si, energy * 1e6) + calcMeanCS(absorbMt103Si, energy * 1e6) + calcMeanCS(absorbMt107Si, energy * 1e6) + 3. * calcMeanCS(absorbMt209Si, energy * 1e6);
                        asAll = 0.68 * asFe + 0.19 * (0.86 * asCr52 + 0.14 * asCr53) + 0.09 * asNi + 0.02 * asSi + 0.02 * asMn55;
                        break;
                    case 33: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asH < 0) asH = 0; if (asC < 0) asC = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        asAll = 4. * asH + 1. * asC;
                        break;
                    case 34: asH = absorbHSpline->Eval(energy1e6Log); if (energy1e6 > 2e8) { asH += calcMeanCS(absorbMt5H, energy1e6) + 3. * calcMeanCS(absorbMt209H, energy1e6);}
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asH < 0) asH = 0; if (asC < 0) asC = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        asAll = 23. * asH + 12. * asC;
                        break;
                    case 36: if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (asC < 0) asC = 0;
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        asAll = asC;
                        break;
                    case 37: if (energy1e6 < 1e6) {asPb206 = absorbPb206Spline->Eval(energy1e6Log);} else {asPb206 = calcMeanCS(absorbPb206, energy1e6);} if (energy1e6 > 3e6) { asPb206 += calcMeanCS(absorbMt5Pb206, energy1e6) + calcMeanCS(absorbMt103Pb206, energy1e6) + calcMeanCS(absorbMt107Pb206, energy1e6) + 3. * calcMeanCS(absorbMt209Pb206, energy1e6); }
                        if (energy1e6 < 1e6) {asPb207 = absorbPb207Spline->Eval(energy1e6Log);} else {asPb207 = calcMeanCS(absorbPb207, energy1e6);} if (energy1e6 > 3e6) { asPb207 += calcMeanCS(absorbMt5Pb207, energy1e6) + calcMeanCS(absorbMt103Pb207, energy1e6) + calcMeanCS(absorbMt107Pb207, energy1e6) + 3. * calcMeanCS(absorbMt209Pb207, energy1e6); }
                        if (energy1e6 < 1e6) {asPb208 = absorbPb208Spline->Eval(energy1e6Log);} else {asPb208 = calcMeanCS(absorbPb208, energy1e6);} if (energy1e6 > 3e6) { asPb208 += calcMeanCS(absorbMt5Pb208, energy1e6) + calcMeanCS(absorbMt103Pb208, energy1e6) + calcMeanCS(absorbMt107Pb208, energy1e6) + 3. * calcMeanCS(absorbMt209Pb208, energy1e6); }
                        if (asPb206 < 0) asPb206 = 0; if (asPb207 < 0) asPb207 = 0; if (asPb208 < 0) asPb208 = 0;
                        //asPb206 = calcMeanCS(absorbPb206, energy * 1e6) + calcMeanCS(absorbMt5Pb206, energy * 1e6) + calcMeanCS(absorbMt103Pb206, energy * 1e6) + calcMeanCS(absorbMt107Pb206, energy * 1e6) + 3. * calcMeanCS(absorbMt209Pb206, energy * 1e6);
                        //asPb207 = calcMeanCS(absorbPb207, energy * 1e6) + calcMeanCS(absorbMt5Pb207, energy * 1e6) + calcMeanCS(absorbMt103Pb207, energy * 1e6) + calcMeanCS(absorbMt107Pb207, energy * 1e6) + 3. * calcMeanCS(absorbMt209Pb207, energy * 1e6);
                        //asPb206 = calcMeanCS(absorbPb208, energy * 1e6) + calcMeanCS(absorbMt5Pb208, energy * 1e6) + calcMeanCS(absorbMt103Pb208, energy * 1e6) + calcMeanCS(absorbMt107Pb208, energy * 1e6) + 3. * calcMeanCS(absorbMt209Pb208, energy * 1e6);
                        asAll = 0.241 * asPb206 + 0.221 * asPb207 + 0.524 * asPb208;
                        break;
                    case 38: asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 1e6) {asN = absorbNSpline->Eval(energy1e6Log) + absorbNbSpline->Eval(energy1e6Log);} else {asN = calcMeanCS(absorbN, energy1e6) + calcMeanCS(absorbNb, energy1e6);} if (energy1e6 > 1e6) {asN += calcMeanCS(absorbMt5N, energy1e6) + calcMeanCS(absorbMt104N,energy1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N,energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N,energy1e6);}
                        if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asN < 0) asN = 0; if (asC < 0) asC = 0;
                        //asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        //asC = calcMeanCS(absorbC, energy * 1e6) + calcMeanCS(absorbMt5C, energy * 1e6) + calcMeanCS(absorbMt103C, energy * 1e6) + calcMeanCS(absorbMt107C, energy * 1e6) + 3. * calcMeanCS(absorbMt209C, energy * 1e6);
                        //asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy * 1e6) + calcMeanCS(absorbMt107N, energy * 1e6) + calcMeanCS(absorbMt108N, energy * 1e6) + 3. * calcMeanCS(absorbMt209N, energy * 1e6);
                        asAll = 0.35 * asC + 0.238 * asH + 0.144 * asN + 0.286 * asO;
                        break;
                    case 41: asH = calcMeanCS(absorbH, energy * 1e6) + calcMeanCS(absorbMt5H, energy * 1e6) + 3. * calcMeanCS(absorbMt209H, energy * 1e6);
                        if (energy1e6 < 1e6) {asC = absorbCSpline->Eval(energy1e6Log);} else {asC = calcMeanCS(absorbC, energy1e6);} if (energy1e6 > 6e6) {asC += calcMeanCS(absorbMt5C, energy1e6) + calcMeanCS(absorbMt103C, energy1e6) + calcMeanCS(absorbMt107C, energy1e6) + 3. * calcMeanCS(absorbMt209C, energy1e6);}
                        if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (asH < 0) asH = 0; if (asO < 0) asO = 0; if (asC < 0) asC = 0;
                        //asH = calcMeanCS(absorbH,energy*1e6)+calcMeanCS(absorbMt5H,energy*1e6)+3.*calcMeanCS(absorbMt209H,energy*1e6);
                        //asO = calcMeanCS(absorbO,energy*1e6)+calcMeanCS(absorbMt5O,energy*1e6)+calcMeanCS(absorbMt103O,energy*1e6)+calcMeanCS(absorbMt105O,energy*1e6)+calcMeanCS(absorbMt107O,energy*1e6)+3.*calcMeanCS(absorbMt209O,energy*1e6);
                        //asC = calcMeanCS(absorbC,energy*1e6)+calcMeanCS(absorbMt5C,energy*1e6)+calcMeanCS(absorbMt103C,energy*1e6)+calcMeanCS(absorbMt107C,energy*1e6)+3.*calcMeanCS(absorbMt209C,energy*1e6);
                        asWater = asH + 2. * asO;
                        asAll = ((6. * asC + 10. * asH + 5. * csO) * nCellulose + asWater * nCelluloseWater) / nCelluloseMix;
                        break;
                    case 100: if (energy1e6 < 1e6) {asO = absorbOSpline->Eval(energy1e6Log);} else {asO = calcMeanCS(absorbO, energy1e6);} if (energy1e6 > 2e6) { asO += calcMeanCS(absorbMt5O, energy1e6) + calcMeanCS(absorbMt103O, energy1e6) + calcMeanCS(absorbMt105O, energy1e6) + calcMeanCS(absorbMt107O, energy1e6) + 3. * calcMeanCS(absorbMt209O, energy1e6);}
                        if (energy1e6 < 2e7) {asN = absorbNSpline->Eval(energy1e6Log) + absorbNbSpline->Eval(energy1e6Log);} else {asN = calcMeanCS(absorbN, energy1e6) + calcMeanCS(absorbNb, energy1e6);} if (energy1e6 > 1e6) {asN += calcMeanCS(absorbMt5N, energy1e6) + calcMeanCS(absorbMt104N,energy1e6) + calcMeanCS(absorbMt105N, energy1e6) + calcMeanCS(absorbMt107N,energy1e6) + calcMeanCS(absorbMt108N, energy1e6) + 3. * calcMeanCS(absorbMt209N,energy1e6);}
                        if (asO < 0) asO = 0; if (asN < 0) asN = 0;
                        //asN = calcMeanCS(absorbN, energy * 1e6) + calcMeanCS(absorbNb, energy * 1e6) + calcMeanCS(absorbMt5N, energy * 1e6) + calcMeanCS(absorbMt104N, energy * 1e6) + calcMeanCS(absorbMt105N, energy * 1e6) + calcMeanCS(absorbMt107N, energy * 1e6) + calcMeanCS(absorbMt108N, energy * 1e6) + 3. * calcMeanCS(absorbMt209N, energy * 1e6);
                        //asO = calcMeanCS(absorbO, energy * 1e6) + calcMeanCS(absorbMt5O, energy * 1e6) + calcMeanCS(absorbMt103O, energy * 1e6) + calcMeanCS(absorbMt105O, energy * 1e6) + calcMeanCS(absorbMt107O, energy * 1e6) + 3. * calcMeanCS(absorbMt209O, energy * 1e6);
                        asLuft = 0.78 * 2. * asN + 0.22 * 2. * asO;
                        asAll = asLuft;
                        break;
                    }

                    if (material == 7)
                    {
                        allMaterials.push_back(csO); allMaterials.push_back(2. * csH); allMaterials.push_back(nSalt * csNa); allMaterials.push_back(nSalt * csCl);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(1); allMaterialsElements.push_back(23); allMaterialsElements.push_back(35);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(asO); allAbsorptionMaterials.push_back(2. * asH); allAbsorptionMaterials.push_back(nSalt * asNa); allAbsorptionMaterials.push_back(nSalt * asCl);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(23); allAbsorptionMaterialsElements.push_back(35);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 9)
                    {
                        allMaterials.push_back(csO); allMaterials.push_back(2. * csH);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(1);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(asO); allAbsorptionMaterials.push_back(2. * asH);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 10)
                    {
                        allMaterials.push_back(0.21 * 2. * csO); allMaterials.push_back(0.78 * 2. * csN); allMaterials.push_back(0.0093 * csAr);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(14); allMaterialsElements.push_back(40);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.21 * 2. * asO); allAbsorptionMaterials.push_back(0.78 * 2. * asN); allAbsorptionMaterials.push_back(0.0093 * asAr);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(14);  allAbsorptionMaterialsElements.push_back(40);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 11)
                    {
                        allMaterials.push_back(0.21 * 2. * csO * rLuft / wLuft + csO * rLuftWater / wWater); allMaterials.push_back(2. * csH * rLuftWater / wWater); allMaterials.push_back(0.78 * 2. * csN * rLuft / wLuft); allMaterials.push_back(0.0093 * csAr * rLuft / wLuft);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(1); allMaterialsElements.push_back(14);  allMaterialsElements.push_back(40);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.21 * 2. * asO * rLuft / wLuft + asO * rLuftWater / wWater); allAbsorptionMaterials.push_back(2. * asH * rLuftWater / wWater); allAbsorptionMaterials.push_back(0.78 * 2. * asN * rLuft / wLuft); allAbsorptionMaterials.push_back(0.0093 * asAr * rLuft / wLuft);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(14);  allAbsorptionMaterialsElements.push_back(40);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 12)
                    {
                        allMaterials.push_back(csSi); allMaterials.push_back(2. * csO);
                        allMaterialsElements.push_back(28); allMaterialsElements.push_back(16);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(asSi); allAbsorptionMaterials.push_back(3. * asO);
                        allAbsorptionMaterialsElements.push_back(28); allAbsorptionMaterialsElements.push_back(16);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 13)
                    {
                        allMaterials.push_back(2. * csAl); allMaterials.push_back(3. * csO);
                        allMaterialsElements.push_back(27); allMaterialsElements.push_back(16);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(2. * asAl); allAbsorptionMaterials.push_back(3. * asO);
                        allAbsorptionMaterialsElements.push_back(27); allAbsorptionMaterialsElements.push_back(16);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 15)
                    {
                        allMaterials.push_back(csFe);
                        allMaterialsElements.push_back(56);
                        allMaterials.at(0) = allMaterials.at(0) / cs;

                        allAbsorptionMaterials.push_back(asFe);
                        allAbsorptionMaterialsElements.push_back(56);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                    }
                    //make a cumulative distribution function
                    if (material == 18)
                    {
                        allBoden.push_back(soilSiFrac * csSi * soilSolidFracVar); allBoden.push_back(2. * soilAlFrac * csAl * soilSolidFracVar); allBoden.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar); allBoden.push_back(2. * csH * soilWaterFrac * soilStrechFactor);
                        allBoden.push_back(wBoden / rBoden * 1e-6 * rBoronInSoil / wBoron * 0.2 * csB10); allBoden.push_back(wBoden / rBoden * 1e-6 * rGdInSoil / wGd155 * 0.148 * csGd155);  allBoden.push_back(wBoden / rBoden * 1e-6 * rGdInSoil / wGd157 * 0.1565 * csGd157);
                        allBoden.push_back(wBoden / rBoden * 1e-6 * rNInSoil / wNitrogen * csN); allBoden.push_back(wBoden / rBoden * 1e-6 * rCInSoil / wCarbon * csC); allBoden.push_back(wBoden / rBoden * 1e-6 * rNaInSoil / wNa * csNa);
                        allBoden.push_back(wBoden / rBoden * 1e-6 * rKInSoil / wKalium * csK); allBoden.push_back(wBoden / rBoden * 1e-6 * rTiInSoil / wTi48 * csTi48 * 0.85);
                        allBoden.push_back(wBoden / rBoden * 1e-6 * rMnInSoil / wMn55 * csMn55); allBoden.push_back(wBoden / rBoden * 1e-6 * rFeInSoil / wFe * csFe);
                        allBodenElements.push_back(28); allBodenElements.push_back(27); allBodenElements.push_back(16);	allBodenElements.push_back(1);
                        allBodenElements.push_back(10); allBodenElements.push_back(155); allBodenElements.push_back(157);
                        allBodenElements.push_back(14); allBodenElements.push_back(12);  allBodenElements.push_back(23); allBodenElements.push_back(39); allBodenElements.push_back(48); allBodenElements.push_back(55); allBodenElements.push_back(56);
                        allBoden.at(0) = allBoden.at(0) / cs;
                        for (int k = 1; k < allBoden.size(); k++) { allBoden.at(k) = allBoden.at(k - 1) + allBoden.at(k) / cs; }

                        allAbsorptionMaterials.push_back(soilSiFrac * asSi * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * soilAlFrac * asAl * soilSolidFracVar); allAbsorptionMaterials.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * asO) * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * asH * soilWaterFrac * soilStrechFactor);
                        allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rBoronInSoil / wBoron * 0.2 * asB10); allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rGdInSoil / wGd155 * 0.148 * asGd155);  allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rGdInSoil / wGd157 * 0.1565 * asGd157);
                        allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rNInSoil / wNitrogen * asN); allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rCInSoil / wCarbon * asC); allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rNaInSoil / wNa * asNa);
                        allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rKInSoil / wKalium * asK); allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rTiInSoil / wTi48 * asTi48 * 0.85);
                        allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rMnInSoil / wMn55 * asMn55); allAbsorptionMaterials.push_back(wBoden / rBoden * 1e-6 * rFeInSoil / wFe * asFe);

                        allAbsorptionMaterialsElements.push_back(28); allAbsorptionMaterialsElements.push_back(27); allAbsorptionMaterialsElements.push_back(16);	allAbsorptionMaterialsElements.push_back(1);
                        allAbsorptionMaterialsElements.push_back(10); allAbsorptionMaterialsElements.push_back(155); allAbsorptionMaterialsElements.push_back(157);
                        allAbsorptionMaterialsElements.push_back(14); allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(23);
                        allAbsorptionMaterialsElements.push_back(39); allAbsorptionMaterialsElements.push_back(48); allAbsorptionMaterialsElements.push_back(55); allAbsorptionMaterialsElements.push_back(56);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 19)
                    {
                        allBoden.push_back(soilSiFrac * csSi * soilSolidFracVar); allBoden.push_back(2. * soilAlFrac * csAl * soilSolidFracVar); allBoden.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar); allBoden.push_back(2. * csH * soilWaterFrac * soilStrechFactor);
                        allBodenElements.push_back(28); allBodenElements.push_back(27); allBodenElements.push_back(16);	allBodenElements.push_back(1);
                        allBoden.at(0) = allBoden.at(0) / cs;
                        for (int k = 1; k < allBoden.size(); k++) { allBoden.at(k) = allBoden.at(k - 1) + allBoden.at(k) / cs; }

                        allAbsorptionMaterials.push_back(soilSiFrac * asSi * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * soilAlFrac * asAl * soilSolidFracVar); allAbsorptionMaterials.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * asO) * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * asH * soilWaterFrac * soilStrechFactor); allAbsorptionMaterials.push_back(asB10 * rBoronInSoil * wBoden / rBoden);
                        allAbsorptionMaterialsElements.push_back(28); allAbsorptionMaterialsElements.push_back(27); allAbsorptionMaterialsElements.push_back(16);	allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(10);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 20)
                    {
                        allBoden.push_back(soilSiFrac * csSi * soilSolidFracVar); allBoden.push_back(2. * soilAlFrac * csAl * soilSolidFracVar); allBoden.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * csO) * soilSolidFracVar); allBoden.push_back(2. * csH * soilWaterFrac * soilStrechFactor);
                        allBodenElements.push_back(28); allBodenElements.push_back(27); allBodenElements.push_back(16);	allBodenElements.push_back(1);
                        allBoden.at(0) = allBoden.at(0) / cs;
                        for (int k = 1; k < allBoden.size(); k++) { allBoden.at(k) = allBoden.at(k - 1) + allBoden.at(k) / cs; }

                        allAbsorptionMaterials.push_back(soilSiFrac * asSi * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * soilAlFrac * asAl * soilSolidFracVar); allAbsorptionMaterials.push_back(((soilSiFrac * 2. + 3. * soilAlFrac) * asO) * soilSolidFracVar); allAbsorptionMaterials.push_back(2. * asH * soilWaterFrac * soilStrechFactor);
                        allAbsorptionMaterialsElements.push_back(28); allAbsorptionMaterialsElements.push_back(27); allAbsorptionMaterialsElements.push_back(16);	allAbsorptionMaterialsElements.push_back(1);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 21)
                    {
                        allPlants.push_back(0.21 * 2. * csO * rLuft / wLuft + csO * rLuftWater / wWater + csO * rPlants / wPlants); allPlants.push_back(2. * csH * rLuftWater / wWater + 2. * csH * rPlants / wPlants); allPlants.push_back(0.2 * csC * rPlants / wPlants);  allPlants.push_back(0.78 * 2. * csN * rLuft / wPlants);
                        allPlantsElements.push_back(16); allPlantsElements.push_back(1); allPlantsElements.push_back(12); allPlantsElements.push_back(14);
                        allPlants.at(0) = allPlants.at(0) / cs;
                        for (int k = 1; k < allPlants.size(); k++) { allPlants.at(k) = allPlants.at(k - 1) + allPlants.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.21 * 2. * asO * rLuft / wLuft + asO * rLuftWater / wWater + asO * rPlants / wPlants); allAbsorptionMaterials.push_back(2. * asH * rLuftWater / wWater + 2. * asH * rPlants / wPlants); allAbsorptionMaterials.push_back(0.2 * asC * rPlants / wPlants);  allAbsorptionMaterials.push_back(0.78 * 2. * asN * rLuft / wPlants);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(14);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 23)
                    {
                        allMaterials.push_back(0.44 * csO); allMaterials.push_back(0.44 * csH); allMaterials.push_back(0.12 * csSi);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(1); allMaterialsElements.push_back(28);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.44 * asO); allAbsorptionMaterials.push_back(0.44 * asH); allAbsorptionMaterials.push_back(0.12 * asSi);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(28);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 24)
                    {
                        allMaterials.push_back(0.5 * csO); allMaterials.push_back(0.14 * csH); allMaterials.push_back(0.11 * csC); allMaterials.push_back(0.25 * csSi);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(1); allMaterialsElements.push_back(12); allMaterialsElements.push_back(28);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.5 * asO); allAbsorptionMaterials.push_back(0.14 * asH); allAbsorptionMaterials.push_back(0.11 * asC); allAbsorptionMaterials.push_back(0.25 * asSi);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(28);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 25)
                    {
                        allMaterials.push_back(2. * csH); allMaterials.push_back(1. * csC);
                        allMaterialsElements.push_back(1); allMaterialsElements.push_back(12);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(2. * asH); allAbsorptionMaterials.push_back(1. * asC);
                        allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 26)
                    {
                        allMaterials.push_back(csAl);
                        allMaterialsElements.push_back(27);
                        allMaterials.at(0) = allMaterials.at(0) / cs;

                        allAbsorptionMaterials.push_back(asAl);
                        allAbsorptionMaterialsElements.push_back(27);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                    }
                    if (material == 27)
                    {
                        allMaterials.push_back(csHe3);
                        allMaterialsElements.push_back(3);
                        allMaterials.at(0) = allMaterials.at(0) / cs;

                        allAbsorptionMaterials.push_back(asHe3);
                        allAbsorptionMaterialsElements.push_back(3);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                    }
                    if (material == 28)
                    {
                        allMaterials.push_back(csB10); allMaterials.push_back(3. * csF);
                        allMaterialsElements.push_back(10); allMaterialsElements.push_back(19);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(asB10); allAbsorptionMaterials.push_back(3. * asF);
                        allAbsorptionMaterialsElements.push_back(10); allAbsorptionMaterialsElements.push_back(19);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 29)
                    {
                        allMaterials.push_back(2. * 0.148 * csGd155); allMaterials.push_back(2. * 0.1565 * csGd157); allMaterials.push_back(3. * csO);
                        allMaterialsElements.push_back(155); allMaterialsElements.push_back(157); allMaterialsElements.push_back(16);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(2. * 0.148 * asGd155); allAbsorptionMaterials.push_back(2. * 0.1565 * asGd157); allAbsorptionMaterials.push_back(3. * asO);
                        allAbsorptionMaterialsElements.push_back(155); allAbsorptionMaterialsElements.push_back(157); allAbsorptionMaterialsElements.push_back(16);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 30)
                    {
                        allMaterials.push_back(2. * csH); allMaterials.push_back(1. * csC); allMaterials.push_back(0.04 * 0.2 * csB10); allMaterials.push_back(0.04 * 0.8 * csB11);
                        allMaterialsElements.push_back(1); allMaterialsElements.push_back(12); allMaterialsElements.push_back(10); allMaterialsElements.push_back(11);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(2. * asH); allAbsorptionMaterials.push_back(1. * asC); allAbsorptionMaterials.push_back(0.04 * 0.2 * asB10); allAbsorptionMaterials.push_back(0.04 * 0.8 * asB11);
                        allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(10); allAbsorptionMaterialsElements.push_back(11);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 31)
                    {
                        allMaterials.push_back(3. * csH); allMaterials.push_back(2. * csC);  allMaterials.push_back(csCl);
                        allMaterialsElements.push_back(1); allMaterialsElements.push_back(12); allMaterialsElements.push_back(35);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(3. * asH); allAbsorptionMaterials.push_back(2. * asC);  allAbsorptionMaterials.push_back(1. * asCl);
                        allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(35);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 32)
                    {
                        allMaterials.push_back(0.68 * csFe); allMaterials.push_back(0.19 * 0.86 * csCr52); allMaterials.push_back(0.19 * 0.14 * csCr53);  allMaterials.push_back(0.09 * csNi); allMaterials.push_back(0.02 * csSi);  allMaterials.push_back(0.02 * csMn55);
                        allMaterialsElements.push_back(56); allMaterialsElements.push_back(52);  allMaterialsElements.push_back(53);  allMaterialsElements.push_back(58); allMaterialsElements.push_back(28);  allMaterialsElements.push_back(55);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.72 * asFe); allAbsorptionMaterials.push_back(0.19 * 0.86 * asCr52); allAbsorptionMaterials.push_back(0.19 * 0.14 * asCr53);  allAbsorptionMaterials.push_back(0.09 * asNi);
                        allAbsorptionMaterialsElements.push_back(56); allAbsorptionMaterialsElements.push_back(52);  allAbsorptionMaterialsElements.push_back(53);  allAbsorptionMaterialsElements.push_back(58);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 33)
                    {
                        allMaterials.push_back(4. * csH); allMaterials.push_back(1. * csC);
                        allMaterialsElements.push_back(1); allMaterialsElements.push_back(12);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(4. * asH); allAbsorptionMaterials.push_back(1. * asC);
                        allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 34)
                    {
                        allMaterials.push_back(23. * csH); allMaterials.push_back(12. * csC);
                        allMaterialsElements.push_back(1); allMaterialsElements.push_back(12);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(23. * asH); allAbsorptionMaterials.push_back(12. * asC);
                        allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(12);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 36)
                    {
                        allMaterials.push_back(csC);
                        allMaterialsElements.push_back(12);
                        allMaterials.at(0) = allMaterials.at(0) / cs;

                        allAbsorptionMaterials.push_back(asC);
                        allAbsorptionMaterialsElements.push_back(12);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                    }
                    if (material == 37)
                    {
                        allMaterials.push_back(0.241 * csPb206); allMaterials.push_back(0.221 * csPb207); allMaterials.push_back(0.524 * csPb208);
                        allMaterialsElements.push_back(206); allMaterialsElements.push_back(207); allMaterialsElements.push_back(208);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.241 * asPb206); allAbsorptionMaterials.push_back(0.221 * asPb207); allAbsorptionMaterials.push_back(0.524 * asPb208);
                        allAbsorptionMaterialsElements.push_back(206); allAbsorptionMaterialsElements.push_back(207); allAbsorptionMaterialsElements.push_back(208);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 38)
                    {
                        allMaterials.push_back(0.35 * csC); allMaterials.push_back(0.238 * csH); allMaterials.push_back(0.144 * csN); allMaterials.push_back(0.286 * csO);
                        allMaterialsElements.push_back(12); allMaterialsElements.push_back(1); allMaterialsElements.push_back(14); allMaterialsElements.push_back(16);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.35 * asC); allAbsorptionMaterials.push_back(0.238 * asH); allAbsorptionMaterials.push_back(0.144 * asN);  allAbsorptionMaterials.push_back(0.286 * asO);
                        allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(14);  allAbsorptionMaterialsElements.push_back(16);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }
                    if (material == 41)
                    {
                        allMaterials.push_back(6.*csC*nCellulose/nCelluloseMix); allMaterials.push_back(csH*(10.*nCellulose/nCelluloseMix + 1.*nCelluloseWater/nCelluloseMix)); allMaterials.push_back(csO*(5.*nCellulose/nCelluloseMix+2.*nCelluloseWater/nCelluloseMix));
                        allMaterialsElements.push_back(12); allMaterialsElements.push_back(1); allMaterialsElements.push_back(16);
                        allMaterials.at(0) = allMaterials.at(0)/cs;
                        for(int k = 1; k < allMaterials.size(); k++) {allMaterials.at(k) = allMaterials.at(k-1)+allMaterials.at(k)/cs;}

                         allAbsorptionMaterials.push_back(6.*asC*nCellulose/nCelluloseMix); allAbsorptionMaterials.push_back(asH*(10.*nCellulose/nCelluloseMix + 1.*nCelluloseWater/nCelluloseMix)); allAbsorptionMaterials.push_back(asO*(5.*nCellulose/nCelluloseMix+2.*nCelluloseWater/nCelluloseMix));
                         allAbsorptionMaterialsElements.push_back(12); allAbsorptionMaterialsElements.push_back(1); allAbsorptionMaterialsElements.push_back(16);
                         allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0)/asAll;
                         for(int k = 1; k < allAbsorptionMaterials.size(); k++) {allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k-1)+allAbsorptionMaterials.at(k)/asAll;}
                     }
                    if (material == 100)
                    {
                        allMaterials.push_back(0.22 * 2. * csO); allMaterials.push_back(0.78 * 2. * csN);
                        allMaterialsElements.push_back(16); allMaterialsElements.push_back(14);
                        allMaterials.at(0) = allMaterials.at(0) / cs;
                        for (int k = 1; k < allMaterials.size(); k++) { allMaterials.at(k) = allMaterials.at(k - 1) + allMaterials.at(k) / cs; }

                        allAbsorptionMaterials.push_back(0.22 * 2. * asO); allAbsorptionMaterials.push_back(0.78 * 2. * asN);
                        allAbsorptionMaterialsElements.push_back(16); allAbsorptionMaterialsElements.push_back(14);
                        allAbsorptionMaterials.at(0) = allAbsorptionMaterials.at(0) / asAll;
                        for (int k = 1; k < allAbsorptionMaterials.size(); k++) { allAbsorptionMaterials.at(k) = allAbsorptionMaterials.at(k - 1) + allAbsorptionMaterials.at(k) / asAll; }
                    }


                    // H, C, O, N, Ar, Al, Si

                    // inelastic cross sections
                    if ((energy > 0.55) && (energy < 50))
                    {
                        switch (material)
                        {
                        case 7:  for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                              for (int k = 0; k < sigmaInNaVec.size(); k++) { allInelastics.push_back(nSalt * calcMeanCS(*(sigmaInNaVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNaAngularVec.at(k))); allInelasticElements.push_back(23); inelasticEnergyLossVec.push_back(inelasticEnergyLossNa.at(k)); };
                              for (int k = 0; k < sigmaInCl35Vec.size(); k++) { allInelastics.push_back(nSalt * calcMeanCS(*(sigmaInCl35Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCl35AngularVec.at(k))); allInelasticElements.push_back(35); inelasticEnergyLossVec.push_back(inelasticEnergyLossCl35.at(k)); } break;
                        case 9: {for (int  k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }; } break;
                        case 10: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(0.21 * 2. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back(0.78 * 2. * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); };
                               for (int k = 0; k < sigmaInArVec.size(); k++) { allInelastics.push_back(0.0093 * calcMeanCS(*(sigmaInArVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInArAngularVec.at(k))); allInelasticElements.push_back(40); inelasticEnergyLossVec.push_back(inelasticEnergyLossAr.at(k)); };  break;
                        case 11: if (energy == lastEnergy11) { allInelastics = allInelastics11; allInelasticAngulars = allInelasticAngulars11; allInelasticElements = allInelasticElements11; inelasticEnergyLossVec = inelasticEnergyLossVec11; }
                               else {
                               for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back((0.21 * 2. * rLuft / wLuft + rLuftWater / wWater) * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back(0.78 * 2. * rLuft / wLuft * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); };
                               for (int k = 0; k < sigmaInArVec.size(); k++) { allInelastics.push_back(0.0093 * rLuft / wLuft * calcMeanCS(*(sigmaInArVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInArAngularVec.at(k))); allInelasticElements.push_back(40); inelasticEnergyLossVec.push_back(inelasticEnergyLossAr.at(k)); };
                               allInelastics11 = allInelastics; allInelasticAngulars11 = allInelasticAngulars; allInelasticElements11 = allInelasticElements; inelasticEnergyLossVec11 = inelasticEnergyLossVec;
                               }
                               break;
                        case 12: { for (int  k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(2. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }; }
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); } break;
                        case 13: { for (int  k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(3. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }; }
                                { for (int  k = 0; k < sigmaInAlVec.size(); k++) { allInelastics.push_back(2. * calcMeanCS(*(sigmaInAlVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInAlAngularVec.at(k))); allInelasticElements.push_back(27); inelasticEnergyLossVec.push_back(inelasticEnergyLossAl.at(k)); }; } break;
                        case 15: { for (int  k = 0; k < sigmaInFeVec.size(); k++) { allInelastics.push_back(calcMeanCS(*(sigmaInFeVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInFeAngularVec.at(k))); allInelasticElements.push_back(56); inelasticEnergyLossVec.push_back(inelasticEnergyLossFe.at(k)); }; } break;
                        case 18: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back((2.25 * soilSolidFracVar + soilStrechFactor * soilWaterFrac) * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16);  inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(soilSiFrac * soilSolidFracVar * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); }
                               for (int k = 0; k < sigmaInAlVec.size(); k++) { allInelastics.push_back(2. * soilAlFrac * soilSolidFracVar * calcMeanCS(*(sigmaInAlVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInAlAngularVec.at(k))); allInelasticElements.push_back(27); inelasticEnergyLossVec.push_back(inelasticEnergyLossAl.at(k)); }
                               for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rNInSoil / wNitrogen * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); }
                               for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rCInSoil / wCarbon * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }
                               for (int k = 0; k < sigmaInNaVec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rNaInSoil / wNa * calcMeanCS(*(sigmaInNaVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNaAngularVec.at(k))); allInelasticElements.push_back(23); inelasticEnergyLossVec.push_back(inelasticEnergyLossNa.at(k)); }
                               for (int k = 0; k < sigmaInK39Vec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rKInSoil / wKalium * calcMeanCS(*(sigmaInK39Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInK39AngularVec.at(k))); allInelasticElements.push_back(39); inelasticEnergyLossVec.push_back(inelasticEnergyLossK39.at(k)); }
                               for (int k = 0; k < sigmaInTi48Vec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rTiInSoil / wTi48 * calcMeanCS(*(sigmaInTi48Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInTi48AngularVec.at(k))); allInelasticElements.push_back(48); inelasticEnergyLossVec.push_back(inelasticEnergyLossTi48.at(k)); }
                               for (int k = 0; k < sigmaInMn55Vec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rMnInSoil / wMn55 * calcMeanCS(*(sigmaInMn55Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInMn55AngularVec.at(k))); allInelasticElements.push_back(55); inelasticEnergyLossVec.push_back(inelasticEnergyLossMn55.at(k)); }
                               for (int k = 0; k < sigmaInFeVec.size(); k++) { allInelastics.push_back(wBoden / rBoden * 1e-6 * rFeInSoil / wFe * calcMeanCS(*(sigmaInFeVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInFeAngularVec.at(k))); allInelasticElements.push_back(56); inelasticEnergyLossVec.push_back(inelasticEnergyLossFe.at(k)); }
                               break;
                        case 19: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back((2.25 * soilSolidFracVar + soilStrechFactor * soilWaterFrac) * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16);  inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(soilSiFrac * soilSolidFracVar * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); }
                               for (int k = 0; k < sigmaInAlVec.size(); k++) { allInelastics.push_back(2. * soilAlFrac * soilSolidFracVar * calcMeanCS(*(sigmaInAlVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInAlAngularVec.at(k))); allInelasticElements.push_back(27); inelasticEnergyLossVec.push_back(inelasticEnergyLossAl.at(k)); }
                               break;
                        case 20: if (energy == lastEnergy20){allInelastics = allInelastics20; allInelasticAngulars = allInelasticAngulars20; allInelasticElements = allInelasticElements20; inelasticEnergyLossVec = inelasticEnergyLossVec20;}
                                else {
                                for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back((2.25 * soilSolidFracVar + soilStrechFactor * soilWaterFrac) * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16);  inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }
                                for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(soilSiFrac * soilSolidFracVar * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); }
                                for (int k = 0; k < sigmaInAlVec.size(); k++) { allInelastics.push_back(2. * soilAlFrac * soilSolidFracVar * calcMeanCS(*(sigmaInAlVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInAlAngularVec.at(k))); allInelasticElements.push_back(27); inelasticEnergyLossVec.push_back(inelasticEnergyLossAl.at(k)); }
                                allInelastics20 = allInelastics; allInelasticAngulars20 = allInelasticAngulars; allInelasticElements20 = allInelasticElements; inelasticEnergyLossVec20 = inelasticEnergyLossVec;
                                }
                            break;
                        case 21: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back((0.21 * 2. * rLuft / wLuft + rLuftWater / wWater + rPlants / wPlants) * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16);  inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(0.2 * rPlants / wPlants * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }
                               for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back((0.78 * 2. * rLuft / wLuft) * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); }    break;
                        case 23: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(0.44 * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(0.12 * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); } break;
                        case 24: for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(0.5 * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(0.11 * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); };
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(0.25 * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); } break;
                        case 25: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  break;
                        case 26: for (int k = 0; k < sigmaInAlVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInAlVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInAlAngularVec.at(k))); allInelasticElements.push_back(27); inelasticEnergyLossVec.push_back(inelasticEnergyLossAl.at(k)); } break;
                        case 28: for (int k = 0; k < sigmaInB10Vec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInB10Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInB10AngularVec.at(k))); allInelasticElements.push_back(10); inelasticEnergyLossVec.push_back(inelasticEnergyLossB10.at(k)); }
                               for (int k = 0; k < sigmaInF19Vec.size(); k++) { allInelastics.push_back(3. * calcMeanCS(*(sigmaInF19Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInF19AngularVec.at(k))); allInelasticElements.push_back(19); inelasticEnergyLossVec.push_back(inelasticEnergyLossF19.at(k)); } break;
                        case 29: for (int k = 0; k < sigmaInGd155Vec.size(); k++) { allInelastics.push_back(2. * 0.148 * calcMeanCS(*(sigmaInGd155Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInGd155AngularVec.at(k))); allInelasticElements.push_back(155); inelasticEnergyLossVec.push_back(inelasticEnergyLossGd155.at(k)); }
                               for (int k = 0; k < sigmaInGd157Vec.size(); k++) { allInelastics.push_back(2. * 0.1565 * calcMeanCS(*(sigmaInGd157Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInGd157AngularVec.at(k))); allInelasticElements.push_back(157); inelasticEnergyLossVec.push_back(inelasticEnergyLossGd157.at(k)); }
                               for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(3. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); }
                        case 30: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  break;
                                for (int k = 0; k < sigmaInB10Vec.size(); k++) { allInelastics.push_back(0.02 * (0.2 * calcMeanCS(*(sigmaInB10Vec.at(k)), energy * 1e6))); allInelasticAngulars.push_back(&(sigmaInB10AngularVec.at(k))); allInelasticElements.push_back(10); inelasticEnergyLossVec.push_back(inelasticEnergyLossB10.at(k)); }
                                for (int k = 0; k < sigmaInB11Vec.size(); k++) { allInelastics.push_back(0.02 * (0.8 * calcMeanCS(*(sigmaInB11Vec.at(k)), energy * 1e6))); allInelasticAngulars.push_back(&(sigmaInB11AngularVec.at(k))); allInelasticElements.push_back(11); inelasticEnergyLossVec.push_back(inelasticEnergyLossB11.at(k)); }  break;
                        case 31: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  ;
                                for (int k = 0; k < sigmaInCl35Vec.size(); k++) { allInelastics.push_back(nSalt * calcMeanCS(*(sigmaInCl35Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCl35AngularVec.at(k))); allInelasticElements.push_back(35); inelasticEnergyLossVec.push_back(inelasticEnergyLossCl35.at(k)); } break;
                        case 32: for (int k = 0; k < sigmaInFeVec.size(); k++) { allInelastics.push_back(0.68 * calcMeanCS(*(sigmaInFeVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInFeAngularVec.at(k))); allInelasticElements.push_back(56); inelasticEnergyLossVec.push_back(inelasticEnergyLossFe.at(k)); };
                               for (int k = 0; k < sigmaInNi58Vec.size(); k++) { allInelastics.push_back(0.09 * calcMeanCS(*(sigmaInNi58Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNi58AngularVec.at(k))); allInelasticElements.push_back(58); inelasticEnergyLossVec.push_back(inelasticEnergyLossNi58.at(k)); };
                               for (int k = 0; k < sigmaInCr52Vec.size(); k++) { allInelastics.push_back(0.19 * 0.86 * calcMeanCS(*(sigmaInCr52Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCr52AngularVec.at(k))); allInelasticElements.push_back(52); inelasticEnergyLossVec.push_back(inelasticEnergyLossCr52.at(k)); };
                               for (int k = 0; k < sigmaInCr53Vec.size(); k++) { allInelastics.push_back(0.19 * 0.14 * calcMeanCS(*(sigmaInCr53Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCr53AngularVec.at(k))); allInelasticElements.push_back(53); inelasticEnergyLossVec.push_back(inelasticEnergyLossCr53.at(k)); }
                               for (int k = 0; k < sigmaInSiVec.size(); k++) { allInelastics.push_back(0.02 * calcMeanCS(*(sigmaInSiVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInSiAngularVec.at(k))); allInelasticElements.push_back(28); inelasticEnergyLossVec.push_back(inelasticEnergyLossSi.at(k)); }
                               for (int k = 0; k < sigmaInMn55Vec.size(); k++) { allInelastics.push_back(0.02 * calcMeanCS(*(sigmaInMn55Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInMn55AngularVec.at(k))); allInelasticElements.push_back(55); inelasticEnergyLossVec.push_back(inelasticEnergyLossMn55.at(k)); } break;
                        case 33: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  break;
                        case 34: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  break;
                        case 36: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(1. * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); }  break;
                        case 37: for (int k = 0; k < sigmaInPb206Vec.size(); k++) { allInelastics.push_back(0.241 * calcMeanCS(*(sigmaInPb206Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInPb206AngularVec.at(k))); allInelasticElements.push_back(206); inelasticEnergyLossVec.push_back(inelasticEnergyLossPb206.at(k)); };
                               for (int k = 0; k < sigmaInPb207Vec.size(); k++) { allInelastics.push_back(0.221 * calcMeanCS(*(sigmaInPb207Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInPb207AngularVec.at(k))); allInelasticElements.push_back(207); inelasticEnergyLossVec.push_back(inelasticEnergyLossPb207.at(k)); };
                               for (int k = 0; k < sigmaInPb208Vec.size(); k++) { allInelastics.push_back(0.524 * calcMeanCS(*(sigmaInPb208Vec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInPb208AngularVec.at(k))); allInelasticElements.push_back(208); inelasticEnergyLossVec.push_back(inelasticEnergyLossPb208.at(k)); } break;
                        case 38: for (int k = 0; k < sigmaInCVec.size(); k++) { allInelastics.push_back(0.350 * calcMeanCS(*(sigmaInCVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k)); };
                               for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(0.286 * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                               for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back(0.144 * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); }; break;
                        case 41: for(int  k = 0; k < sigmaInCVec.size(); k++) {allInelastics.push_back(6.*nCellulose/nCelluloseMix *calcMeanCS(*(sigmaInCVec.at(k)),energy*1e6)); allInelasticAngulars.push_back(&(sigmaInCAngularVec.at(k))); allInelasticElements.push_back(12); inelasticEnergyLossVec.push_back(inelasticEnergyLossC.at(k));};
                                 for(int  k = 0; k < sigmaInOVec.size(); k++) {allInelastics.push_back((10.*nCellulose/nCelluloseMix + 1.*nCelluloseWater/nCelluloseMix)*calcMeanCS(*(sigmaInOVec.at(k)),energy*1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k));}; break;
                        case 100: {for (int k = 0; k < sigmaInOVec.size(); k++) { allInelastics.push_back(0.22 * 2. * calcMeanCS(*(sigmaInOVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInOAngularVec.at(k))); allInelasticElements.push_back(16); inelasticEnergyLossVec.push_back(inelasticEnergyLossO.at(k)); };
                                for (int k = 0; k < sigmaInNVec.size(); k++) { allInelastics.push_back(0.78 * 2. * calcMeanCS(*(sigmaInNVec.at(k)), energy * 1e6)); allInelasticAngulars.push_back(&(sigmaInNAngularVec.at(k))); allInelasticElements.push_back(14); inelasticEnergyLossVec.push_back(inelasticEnergyLossN.at(k)); }; } break;
                        }

                        //make a cumulative distribution function
                        if (allInelastics.size() > 0)
                        {
                            for (int k = 0; k < allInelastics.size(); k++) { csIn += allInelastics.at(k); }
                            allInelastics.at(0) = allInelastics.at(0) / csIn;
                            for (int k = 1; k < allInelastics.size(); k++) { allInelastics.at(k) = allInelastics.at(k - 1) + allInelastics.at(k) / csIn; }
                        }
                    }
                    // in the case of those materials the last energy and cross section assigned to it are stored in order to not calculate it twice when traversing through layers
                    lastEnergy11 = 0; lastEnergy19 = 0; lastEnergy20 = 0; lastEnergy25 = 0;
                    if (material == 11) lastEnergy11 = energy;
                    if (material == 19) lastEnergy19 = energy;
                    if (material == 20) lastEnergy20 = energy;
                    if (material == 25) lastEnergy25 = energy;

                    lastEnergy = energy;

                    // CALCULATION OF TOTAL CROSSSECTION AND SELECTION OF PROCESSES
                    sumCs = cs + asAll + csIn;
                    //Elastic + Absorption + Ienalstic
                    sumCsArray[0] = cs / sumCs; sumCsArray[1] = sumCsArray[0] + asAll / sumCs;  sumCsArray[2] = sumCsArray[1] + csIn / sumCs;


                    // DECISION FOR THE TYPE OF SCATTERING PROCESS
                    tempRnd = r.Rndm();
                    if (tempRnd < sumCsArray[0]) { scatteredElastic = true; scatteredInelastic = false; sAbsorbed = false; }
                    else
                        if (tempRnd < sumCsArray[1]) { sAbsorbed = true; scatteredElastic = false; scatteredInelastic = false; }
                        else { scatteredInelastic = true; scatteredElastic = false; sAbsorbed = false; }


                    if ((material == 21) || (material == 11))
                    {
                        // cout<<energy<<" "<<material<<" "<<sumCs<<" "<<cs<<" "<<asAll<<" "<<csIn<<endl;
                    }

                    sumCs = sumCs * rGeneral;

                    // calculate range
                    switch (material)
                    {
                    case 7:  wwRange = -TMath::Log(r.Rndm()) / getWWProb(rSaltWater, wSaltWater, sumCs, 0);	break;
                    case 9:  wwRange = -TMath::Log(r.Rndm()) / getWWProb(rWater, wWater, sumCs, 0);	break;
                    case 10: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rLuft, wLuft, sumCs, 0);	break;
                        //case 11: wwRange = -TMath::Log(r.Rndm())/(getWWProb(rLuft, wLuft, sumCs,0) + getWWProb(rLuftWater, wWater, csW, 0)) ; break;
                    case 11: wwRange = -TMath::Log(r.Rndm()) / getWWProb(1.0001, 1.0001, sumCs, 0); break; //the density is already included in sumCs
                    case 12: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rQuarz, wQuarz, sumCs, 0);  break;
                    case 13: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rAl2O3, wAl2O3, sumCs, 0);   break;
                    case 15: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rFe, wFe, sumCs, 0);   break;
                    case 18: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rBoden, wBoden, sumCs, 0);   break;
                    case 19: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rBoden, wBoden, sumCs, 0);   break;
                    case 20: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rBoden, wBoden, sumCs, 0);   break;
                    case 21: wwRange = -TMath::Log(r.Rndm()) / getWWProb(1.0001, 1.0001, sumCs, 0);   break; //the density is already included in sumCs
                    case 23: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rMaterial, wCatLitter, sumCs, 0);   break;
                    case 24: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rMaterial, wAsphalt, sumCs, 0);   break;
                    case 25: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rHDPE, wHDPE, sumCs, 0);   break;
                    case 26: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rAlMg, wAlu, sumCs, 0);   break;
                    case 27: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rHe3, wHe3, sumCs, 0);   break;
                    case 28: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rBF3, wBF3, sumCs, 0);   break;
                    case 29: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rGd2O3, wGd2O3, sumCs, 0);   break;
                    case 30: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rHDPE, wHDPE, sumCs, 0);   break;
                    case 31: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rPVC, wPVC, sumCs, 0);   break;
                    case 32: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rEStahl, wEStahl, sumCs, 0);   break;
                    case 33: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rMethane, wMethane, sumCs, 0);   break;
                    case 34: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rDiesel, wDiesel, sumCs, 0);   break;
                    case 36: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rGraphit, wGraphit, sumCs, 0);   break;
                    case 37: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rPb, wPb, sumCs, 0);   break;
                    case 38: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rTNT, wTNT, sumCs, 0);   break;
                    case 41: wwRange = -TMath::Log(r.Rndm())/getWWProb(rCelluloseMix, wCelluloseMix, sumCs,0);   break;
                    case 100: wwRange = -TMath::Log(r.Rndm()) / getWWProb(rLuft, wLuft, sumCs, 0);	break;
                    }


                    //deprecated
                    switch (material)
                    {
                    case 7:   weightThermal = wSaltWater;	break;
                    case 9:   weightThermal = wWater;	break;
                    case 10:  weightThermal = wLuft;	break;
                    case 11:  weightThermal = wLuft; break;
                    case 12:  weightThermal = wQuarz;  break;
                    case 13:  weightThermal = wAl2O3;   break;
                    case 15:  weightThermal = wFe;   break;
                    case 18:  weightThermal = wBoden;  break;
                    case 19:  weightThermal = wBoden;  break;
                    case 20:  weightThermal = wBoden;  break;
                    case 21:  weightThermal = wWater;  break;
                    case 23:  weightThermal = wCatLitter;   break;
                    case 24:  weightThermal = wAsphalt;   break;
                    case 25:  weightThermal = wHDPE;   break; //wHDPE
                    case 26:  weightThermal = wAlu;   break;
                    case 27:  weightThermal = wHe3;   break;
                    case 28:  weightThermal = wBF3;   break;
                    case 29:  weightThermal = wGd2O3;   break;
                    case 30:  weightThermal = wHDPE;   break;
                    case 31:  weightThermal = wPVC;   break;
                    case 32:  weightThermal = wEStahl;   break;
                    case 33:  weightThermal = wMethane;   break;
                    case 34:  weightThermal = wHDPE;   break; //that is Diesel
                    case 36:  weightThermal = wGraphit;  break;
                    case 37:  weightThermal = wPb;  break;
                    case 38:  weightThermal = wTNT;  break;
                    case 41:  weightThermal = wCellulose;  break;
                    case 100:  weightThermal = wLuft;	break;
                    }

                    if ((!scatteredThisLayer) && (!continuedThisLayer))
                    {
                        length = fabs(geometries.at(g)[5] / cosTheta);
                    }
                    else
                    {
                        if (!reverseDir) { length = fabs((geometries.at(g)[4] + geometries.at(g)[5] - z0) / cosTheta); }
                        else { length = fabs((z0 - geometries.at(g)[4]) / cosTheta); }
                    }

                    detectorRealisticallyHitted = false;
                    layerRealisticallyHitted = false;

                    if ((useRealisticModelDetector) || (useRealisticModelLayer)) useRealisticModel = true; else useRealisticModel = false;

                    if ((useRealisticModel) && ((g == detectorLayer) || ((detectorLayerOverride) || (detectorOverride)) || (additionalDetectorLayers[g] >= 0)))
                    {
                        //transform the angle of 0....pi to 0...pi/2.
                        if (!useAdditionalDetectorModel)
                        {
                            if (detectorAngleModel->Eval(piHalf - fabs(theta - piHalf)) > r.Rndm())
                            {
                                detectorRealisticallyHitted = true;
                                layerRealisticallyHitted = true;
                            }
                        }
                        else { detectorRealisticallyHitted = true; layerRealisticallyHitted = true; }

                        if ((detectorRealisticallyHitted) || (layerRealisticallyHitted))
                        {
                            //hits the top or just uses one response function
                            if ((!useAdditionalDetectorModel) || (!reverseDir))
                            {
                                if (detectorEnergyModel->Eval(TMath::Log10(energy)) > r.Rndm())
                                {
                                    detectorRealisticallyHitted = true; layerRealisticallyHitted = true;
                                }
                                else
                                {
                                    detectorRealisticallyHitted = false; layerRealisticallyHitted = false;
                                }
                            }
                            else
                            {
                                // probability of hitting the bottom is cos(theta)*A_bottom / (cos(theta)*A_bottom + sin(theta)*A_side)
                                if (r.Rndm() < (110. / (110. + fabs(tanTheta) * 37.)))
                                {//bottom

                                    if (detectorEnergyModel2->Eval(TMath::Log10(energy)) > r.Rndm())
                                    {
                                        detectorRealisticallyHitted = true; layerRealisticallyHitted = true;
                                    }
                                    else
                                    {
                                        detectorRealisticallyHitted = false; layerRealisticallyHitted = false;
                                    }
                                }
                                else
                                {//side
                                    if (detectorEnergyModel->Eval(TMath::Log10(energy)) > r.Rndm())
                                    {
                                        detectorRealisticallyHitted = true; layerRealisticallyHitted = true;
                                    }
                                    else
                                    {
                                        detectorRealisticallyHitted = false; layerRealisticallyHitted = false;
                                    }
                                }
                            }
                        }
                    }

                   if (((drawDensityMap) && (g == detectorLayer) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) || (detectorLayerOverride))
                   {
                        if (((useDetectorSensitiveMaterial) && (material == detectorSensitiveMaterial)) || (!useDetectorSensitiveMaterial))
                        {
                            if ((theta > downwardScotomaAngle) && (theta < downwardAcceptanceAngle))
                            {
                                if ((energy < 2e-7) && (energy > 1e-9)) {if (scatteredThisLayer) densityMapThermal->Fill(x * 0.001, y * 0.001); else {if (!reverseDir) densityMapThermal->Fill(xt * 0.001, yt * 0.001); else densityMapThermal->Fill(xtEnd * 0.001, yt * 0.001); } }

                                if ((energy < 0.001) && (energy > 0.000001))
                                {
                                    if (scatteredThisLayer) densityMap->Fill(x * 0.001, y * 0.001);
                                    else
                                    {
                                        if (!reverseDir) densityMap->Fill(xt * 0.001, yt * 0.001); else densityMap->Fill(xtEnd * 0.001, yt * 0.001);
                                    }
                                }
                                if ((energy < 0.5) && (energy > 0.001)) { if (scatteredThisLayer) densityMapIntermediate->Fill(x * 0.001, y * 0.001);  else { if (!reverseDir) densityMapIntermediate->Fill(xt * 0.001, yt * 0.001);  else densityMapIntermediate->Fill(xtEnd * 0.001, yt * 0.001); } }
                                if ((energy < 10) && (energy > 0.5)) { if (scatteredThisLayer) densityMapFast->Fill(x * 0.001, y * 0.001);  else { if (!reverseDir) densityMapFast->Fill(xt * 0.001, yt * 0.001);  else densityMapFast->Fill(xtEnd * 0.001, yt * 0.001); } }
                                if (((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) { if (scatteredThisLayer) densityMapAlbedo->Fill(x * 0.001, y * 0.001);  else { if (!reverseDir) densityMapAlbedo->Fill(xt * 0.001, yt * 0.001);  else densityMapAlbedo->Fill(xtEnd * 0.001, yt * 0.001); } }
                                if ((energy < 100000) && (energy > 20)) { if (scatteredThisLayer) densityMapHighEnergy->Fill(x * 0.001, y * 0.001);  else { if (!reverseDir) densityMapHighEnergy->Fill(xt * 0.001, yt * 0.001); else densityMapHighEnergy->Fill(xtEnd * 0.001, yt * 0.001); } }
                            }
                        }
                    }

                    if (totalAdditionalDetectorLayers > 0)
                    {
                        if (((additionalDetectorLayers[g] >= 0) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) )
                        {
                            if (((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) { if (scatteredThisLayer) addDetLayerVec.at(additionalDetectorLayers[g])->Fill(x * 0.001, y * 0.001);  else { if (!reverseDir) addDetLayerVec.at(additionalDetectorLayers[g])->Fill(xt * 0.001, yt * 0.001);  else addDetLayerVec.at(additionalDetectorLayers[g])->Fill(xtEnd * 0.001, yt * 0.001); } }
                        }
                    }

                    if (((g == detectorLayer) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) || (detectorLayerOverride))
                    {
                        if (((useDetectorSensitiveMaterial) && (material == detectorSensitiveMaterial)) || (!useDetectorSensitiveMaterial))
                        {
                            //if ((((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) && (theta > downwardScotomaAngle) && (theta < downwardAcceptanceAngle)) detectorSpectrumNear->Fill(energy);
                        }
                    }

                    detHitted = false;
                    detHitSelected = false;
                    distToDet = 1e12;
                    distToCyl = 1e12;

                    if (((g == detectorLayer) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) || (detectorLayerOverride))
                    {
                        subsurfaceScatteringMean = 0;

                        if (((useDetectorSensitiveMaterial) && (material == detectorSensitiveMaterial)) || (!useDetectorSensitiveMaterial))
                        {
                            if ((((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) && (theta > downwardScotomaAngle) && (theta < downwardAcceptanceAngle))
                            {
                                //detectorLayerTheta2->Fill(theta, fabs(1./sin(theta)) );
                                detectorLayerTheta->Fill(theta);
                            }
                        }

                        if ((((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) && (theta > downwardScotomaAngle) && (theta < downwardAcceptanceAngle))
                        {
                            if ((scatteredThisLayer) || (continuedThisLayer))  xySqrt = sqrt(TMath::Power(xAtInterface - x, 2) + TMath::Power(yAtInterface - y, 2));
                            else  if (theta < piHalf) xySqrt = sqrt(TMath::Power(xAtInterface - xt, 2) + TMath::Power(yAtInterface - yt, 2));
                            else xySqrt = sqrt(TMath::Power(xAtInterface - xtEnd, 2) + TMath::Power(yAtInterface - ytEnd, 2));

                            detectorLayerDistance->Fill(xySqrt);
                            if (hasPassedSurface)  detectorLayerDistanceBackscattered->Fill(xySqrt);
                            if (hasPassedSurface) { if (recordSubsurfaceScatterings) { for (int sElm = 0; sElm < subsurfaceScatterings.size(); sElm++) { detectorDistanceDepth2->Fill(xySqrt, subsurfaceScatterings.at(sElm)); } } }
                            //if (hasPassedSurface) { if (recordSubsurfaceScatterings) { subsurfaceScatteringMean = 0; for (int sElm = 0; sElm < subsurfaceScatterings.size(); sElm++) {subsurfaceScatteringMean += subsurfaceScatterings.at(sElm); } if (subsurfaceScatterings.size() > 0) { subsurfaceScatteringMean = subsurfaceScatteringMean/( subsurfaceScatterings.size()*1.);} }}
                            if ((hasPassedSurface) && (recordSubsurfaceScatterings) && (subsurfaceScatterings.size() > 0))
                            {
                                if (subsurfaceScatterings.size() <= 9)
                                {
                                    subsurfaceScatteringMean = subsurfaceScatterings.back();
                                }
                                else
                                {
                                    subsurfaceScatteringMean = 0;
                                    sort(subsurfaceScatterings.begin(), subsurfaceScatterings.end());

                                    for (int sElm = 0; sElm < subsurfaceScatterings.size(); sElm++)
                                    {
                                        if (sElm > 0.864 * subsurfaceScatterings.size()) {subsurfaceScatteringMean = subsurfaceScatterings.at(sElm); break;} else {}
                                    }
                                }
                            }

                            if (calcNeutronTime)
                            {
                                if ((scatteredThisLayer) || (continuedThisLayer))  timeNtr += calcNeutronDiffTime(z0, z0Alt, energy, cosTheta);
                                else  if (theta < piHalf) timeNtr += calcNeutronDiffTime(tempz, z0, energy, cosTheta);
                                else timeNtr += calcNeutronDiffTime(tempzEnd, z0, energy, cosTheta);

                                if (timeNtr > 1)  detectorLayerTimeTrans->Fill(timeNtr);
                            }
                        }
                    }

                    //this one is the new solution
                    if (((g == detectorLayer) || (detectorOverride)) && (doSingleDetector) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer)))
                    {
                        if (detectorOverride)
                        {
                            detHitted = true;
                        }
                        else
                        {
                            for (lCounter = 0; lCounter < 1; lCounter++)
                            {
                                if (useCylindricalDetector)
                                {
                                    distToCyl = getDistanceToLine(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight);

                                    if (distToCyl < detRad)
                                    {
                                        if (getDistanceToPoint(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight) < sqrt(pow(0.5 * detectorHeight, 2) + pow(detRad, 2)))
                                        {
                                            distToDet = distToCyl;
                                        }
                                    }
                                    else distToDet = 1e9;
                                }

                                else
                                {
                                    distToDet = getDistanceToPoint(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight);
                                }

                                if (distToDet <= detRad)
                                {
                                    tempDist = sqrt(TMath::Power((x - detPosX[lCounter]), 2) + TMath::Power((y - detPosY[lCounter]), 2) + TMath::Power((z0 - detectorHeight), 2));

                                    if (true)
                                    {
                                        detHitSelected = true;

                                        //check if the neutron anyway intersects with the layer
                                        if (useCylindricalDetector)
                                        {
                                            tempDist3 = sqrt(TMath::Power(xt - detPosX[lCounter], 2) + TMath::Power(yt - detPosY[lCounter], 2));

                                            if (tempDist3 < detRad) detHitted = true;

                                            tempDist4 = sqrt(TMath::Power(xtEnd - detPosX[lCounter], 2) + TMath::Power(ytEnd - detPosY[lCounter], 2));

                                            if (tempDist4 < detRad) detHitted = true;
                                        }
                                        if (useSphericalDetector)
                                        {
                                            if (tempDist < detRad)
                                            {
                                                detHitted = true;
                                            }
                                        }
                                    }
                                }

                                //check if the neutron anyway intersects with the layer
                                if (useySheetDetector)
                                {
                                    vectorDirFactor = (detPosY[lCounter] - y) / (sinPhi * tanTheta);
                                    zPlane = z0 + vectorDirFactor;
                                    xPlane = x + vectorDirFactor * cosPhi * tanTheta;

                                    if ((fabs(xPlane - detPosX[lCounter]) < 0.5 * detLength) && (fabs(zPlane - detectorHeight) < 0.5 * geometries.at(detectorLayer)[5]))
                                    {
                                        if (fabs(vectorDirFactor) < wwRange) detHitSelected = true;
                                    }
                                }

                                if (usexSheetDetector)
                                {
                                    vectorDirFactor = (detPosX[lCounter] - x) / (cosPhi * tanTheta);
                                    zPlane = z0 + vectorDirFactor;
                                    yPlane = y + vectorDirFactor * sinPhi * tanTheta;

                                    if ((fabs(yPlane - detPosY[lCounter]) < 0.5 * detLength) && (fabs(zPlane - detectorHeight) < 0.5 * geometries.at(detectorLayer)[5]))
                                    {
                                        if (fabs(vectorDirFactor) < wwRange) detHitSelected = true;
                                        //tempDist3 = sqrt(TMath::Power(x-detPosX[lCounter],2)+TMath::Power(y-detPosY[lCounter],2)+TMath::Power(z0-detectorHeight,2));
                                    }
                                }

                                if ((detHitSelected) && (!detHitted)) //that is detHitSelected from cylindrical and spherical
                                {
                                    if (scalarProduct(detPosX[lCounter] - x, detPosY[lCounter] - y, detectorHeight - z0, cosPhi * fabs(tanTheta * cosTheta), sinPhi * fabs(tanTheta * cosTheta), cosTheta) >= 0)
                                    {
                                        detHitted = true;
                                    }
                                }
                                if (detHitted)  detectorID = lCounter;
                            }
                        }


                        // xySqrt calculation for xSheet and ySheet detector misssing!

                        if ((detHitted) && (useCylindricalDetector) && (!detectorOverride))
                        {
                            if (TMath::Power(xtEnd - detPosX[lCounter], 2) + TMath::Power(ytEnd - detPosY[lCounter], 2) < detRad * detRad)
                            {
                                xySqrt = sqrt(TMath::Power(xtEnd - xAtInterface, 2) + TMath::Power(ytEnd - yAtInterface, 2));
                            }
                            else
                            {
                                if (TMath::Power(xt - detPosX[lCounter], 2) + TMath::Power(yt - detPosY[lCounter], 2) < detRad * detRad)
                                {
                                    xySqrt = sqrt(TMath::Power(xt - xAtInterface, 2) + TMath::Power(yt - yAtInterface, 2));
                                }
                                else
                                {
                                    tValue = intersectCylinderMantle(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight, detRad);

                                    if (((z0 + tValue * cosTheta) < tempzEnd) && ((z0 + tValue * cosTheta) > tempz))
                                    {
                                        xIs = x + tValue * sin(theta) * cosPhi;
                                        yIs = y + tValue * sin(theta) * sinPhi;

                                        xySqrt = sqrt(pow(xAtInterface - xIs, 2) + pow(yAtInterface - yIs, 2));
                                    }
                                    else //too high or too low
                                    {
                                        detHitted = false;
                                    }
                                }
                            }
                        }

                        if ((detHitted) && (useSphericalDetector) && (!detectorOverride))
                        {
                            if (2. * detRad > geometries.at(detectorLayer)[5])
                            {
                                //cut of the sphere at the detector layer
                                tempDist2 = sqrt(TMath::Power(detRad, 2) - TMath::Power(0.5 * geometries.at(detectorLayer)[5], 2));

                                tempDist3 = sqrt(TMath::Power(xt - detPosX[lCounter], 2) + TMath::Power(yt - detPosY[lCounter], 2));

                                if (tempDist3 < tempDist2) xySqrt = sqrt(TMath::Power(xt - xAtInterface, 2) + TMath::Power(yt - yAtInterface, 2));     //top sphere cut
                                else
                                {
                                    tempDist4 = sqrt(TMath::Power(xtEnd - detPosX[lCounter], 2) + TMath::Power(ytEnd - detPosY[lCounter], 2));
                                    if (tempDist4 < tempDist2)  xySqrt = sqrt(TMath::Power(xtEnd - xAtInterface, 2) + TMath::Power(ytEnd - yAtInterface, 2)); //bottom sphere cut
                                    else
                                    {
                                        xySqrt = calcCylindricalHitDist(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight, detRad, xAtInterface, yAtInterface);
                                    }
                                }
                            }
                            else
                            {
                                xySqrt = calcCylindricalHitDist(x, y, z0, theta, phi, detPosX[lCounter], detPosY[lCounter], detectorHeight, detRad, xAtInterface, yAtInterface);
                            }
                        }

                        if ((continuedThisLayer) || (detectorOverride))
                        {
                            xySqrt = sqrt(TMath::Power(xAtInterface - x, 2) + TMath::Power(yAtInterface - y, 2));
                            detectorID = detectorOverrideID;
                        }
                    }


                    // this is the virtual detector
                    if ((detHitted) && (doSingleDetector))
                    {
                        if ((((energy < energyhighTHL) && (energy > energylowTHL) && (!useRealisticModelDetector)) || (detectorRealisticallyHitted)) && (theta > downwardScotomaAngle) && (theta < downwardAcceptanceAngle))
                        {
                            detectorDistance->Fill(xySqrt);
                            if (hasPassedSurface)
                            {
                                detectorSpectrumNear->Fill(energy);
                                detectorDistanceBackScattered->Fill(xySqrt);

                                if (true)
                                {
                                    detectorOriginMap->Fill(xAtInterface / 1000., yAtInterface / 1000.);
                                }
                            }
                            nDetectedNeutrons++;

                            if (!dontFillunecessaryPlots)
                            {
                                detectorTheta->Fill(theta);
                                detectorTheta2->Fill(theta, fabs(1. / sin(theta)));
                                //if (hasPassedSurface) detectorTheta2->Fill(theta);
                                detectorPhi->Fill(fabs(fmod(phi, 2. * pi)));
                            }

                            if (calcNeutronTime)
                            {
                                if (timeNtr > 1) detectorTimeTrans->Fill(timeNtr + calcNeutronDiffTime(detectorHeight, z0, energy, cosTheta));
                            }
                        }

                        if (!dontFillunecessaryPlots)
                        {
                            detectorSpectrum->Fill(energy);
                            detectorSpectrumOrigin->Fill(energyAtInterface);

                            detectorDistanceDepth->Fill(xySqrt, z0);
                            detectorDistanceMaxDepth->Fill(xySqrt, z0max);

                            if (hasPassedSurface) detectorDistanceEnergy->Fill(xySqrt, energy);
                            if (hasPassedSurface) detectorSpectrumOriginMoisture->Fill(moistureAtInterface);
                        }

                        if (detFileOutput)
                        {
                            if (doBatchRun2D) detOutputFile.open(outputFolder + "detectorNeutronHitData_" + castIntToString(paramInt) + ".dat", ios::out | ios::app);
                            else detOutputFile.open(outputFolder + "detectorNeutronHitData.dat", ios::out | ios::app);

                            detOutputFile << detectorID << "\t" << n << "\t" << scatterings << "\t" << xLastScattered * 0.001 << "\t" << yLastScattered * 0.001 << "\t" << zLastScattered * 0.001 << "\t" << theta << "\t" << phi << "\t" << energy << "\t" << energyAtInterface << "\t" << xySqrt * 0.001 << "\t" << xAtInterface * 0.001 << "\t" << yAtInterface * 0.001 << "\t" << zAtInterface * 0.001 << "\t" << previousSoilX * 0.001 << "\t" << previousSoilY * 0.001 << "\t" << previousSoilZ0 * 0.001 << "\t" << thermalizedX * 0.001 << "\t" << thermalizedY * 0.001 << "\t" << thermalizedZ0 * 0.001 << "\t" << z0max * 0.001 << "\t" << subsurfaceScatteringMean * 0.001 << "\t" << timeNtr << "\t" << moistureAtInterface << "\t" << detectorRealisticallyHitted << "\t" << hasPassedSurface << "\n";

                            detOutputFile.close();
                        }

                        if (detTrackFileOutput)
                        {
                            if (neutronTrackCoordinatesFullSet.size() > 0)
                            {
                                if (doBatchRun2D) detTrackOutputFile.open(outputFolder + "detectorNeutronTrackHitData_" + castIntToString(paramInt) + ".dat", ios::out | ios::app);
                                else detTrackOutputFile.open(outputFolder + "detectorNeutronTrackHitData.dat", ios::out | ios::app);

                                for (int nt = 0; nt <= (neutronTrackCoordinatesFullSet.size() - 9); nt += 9)
                                {
                                    if ((nt > 0) && ((neutronTrackCoordinatesFullSet[nt + 5] != neutronTrackCoordinatesFullSet[nt - 4]) || (neutronTrackCoordinatesFullSet[nt + 7] == 0) || (neutronTrackCoordinatesFullSet[nt + 7] == 1)))
                                    {
                                        detTrackOutputFile << detectorID << "\t" << n << "\t" << neutronTrackCoordinatesFullSet[nt + 0] * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << "\t" << detectorRealisticallyHitted << endl;
                                    }
                                    else
                                    {
                                        if (nt == 0) detTrackOutputFile << detectorID << "\t" << n << "\t" << neutronTrackCoordinatesFullSet.at(nt + 0) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << "\t" << detectorRealisticallyHitted << endl;
                                    }
                                }
                                detTrackOutputFile.close();
                            }
                        }

                        if (detectorAbsorbing)
                        {
                            absorbBreak = true;
                            neutronAbsorbedbyDetector = true;
                        }
                        //breaks if detector absorbs
                        //break;
                    }

                    if (detLayerFileOutput)
                    {
                        if (((g == detectorLayer) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) || (detectorLayerOverride))
                        {
                            if (doBatchRun2D) detLayerOutputFile.open(outputFolder + "detectorLayerNeutronHitData_" + castIntToString(paramInt) + ".dat", ios::out | ios::app);
                            else detLayerOutputFile.open(outputFolder + "detectorLayerNeutronHitData.dat", ios::out | ios::app);

                            detLayerOutputFile << n << "\t" << scatterings << "\t";

                            if (scatteredThisLayer)  detLayerOutputFile << x * 0.001 << "\t" << y * 0.001 << "\t" << z0 * 0.001 << "\t";
                            else
                            {
                                if (theta < piHalf) detLayerOutputFile << xt * 0.001 << "\t" << yt * 0.001 << "\t" << tempz * 0.001 << "\t";
                                else detLayerOutputFile << xtEnd * 0.001 << "\t" << ytEnd * 0.001 << "\t" << tempzEnd * 0.001 << "\t";
                            }

                            //detLayerOutputFile << xLastScattered*0.001  <<"\t"<< yLastScattered*0.001  <<"\t"<< zLastScattered*0.001  <<"\t"<< theta <<"\t"<< phi <<"\t"<< energy <<"\t"<< energyAtInterface <<"\t"<< xySqrt*0.001 << "\t"<< xAtInterface*0.001 <<"\t"<< yAtInterface*0.001 <<"\t"<< zAtInterface*0.001 <<"\t"<< z0max*0.001 <<"\t"<<moistureAtInterface<<"\t"<<layerRealisticallyHitted<<"\t"<<hasPassedSurface<<"\n";
                            detLayerOutputFile << xLastScattered * 0.001 << "\t" << yLastScattered * 0.001 << "\t" << zLastScattered * 0.001 << "\t" << theta << "\t" << phi << "\t" << energy << "\t" << energyAtInterface << "\t" << xySqrt * 0.001 << "\t" << xAtInterface * 0.001 << "\t" << yAtInterface * 0.001 << "\t" << zAtInterface * 0.001 << "\t" << previousSoilX * 0.001 << "\t" << previousSoilY * 0.001 << "\t" << previousSoilZ0 * 0.001 << "\t" << thermalizedX * 0.001 << "\t" << thermalizedY * 0.001 << "\t" << thermalizedZ0 * 0.001 << "\t" << z0max * 0.001 << "\t" << subsurfaceScatteringMean * 0.001 << "\t" << timeNtr << "\t" <<  moistureAtInterface << "\t" << layerRealisticallyHitted << "\t" << hasPassedSurface << "\n";
                            detLayerOutputFile.close();
                        }
                    }

                    if (detLayerTrackFileOutput)
                    {
                        if (((g == detectorLayer) && (!continuedThisLayer) && ((!noMultipleScatteringRecording) || (!scatteredThisLayer))) || (detectorLayerOverride))
                        {
                            if (neutronTrackCoordinatesFullSet.size() > 0)
                            {
                                if (doBatchRun2D) detLayerTrackOutputFile.open(outputFolder + "detectorLayerNeutronTrackHitData_" + castIntToString(paramInt) + ".dat", ios::out | ios::app);
                                else detLayerTrackOutputFile.open(outputFolder + "detectorLayerNeutronTrackHitData.dat", ios::out | ios::app);

                                for (int nt = 0; nt <= (neutronTrackCoordinatesFullSet.size() - 9); nt += 9)
                                {
                                    if ((nt > 0) && ((neutronTrackCoordinatesFullSet[nt + 5] != neutronTrackCoordinatesFullSet[nt - 4]) || (neutronTrackCoordinatesFullSet[nt + 7] == 0) || (neutronTrackCoordinatesFullSet[nt + 7] == 1)))
                                    {
                                        detLayerTrackOutputFile << scatterings << "\t" << n << "\t" << neutronTrackCoordinatesFullSet[nt + 0] * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << "\t" << layerRealisticallyHitted << endl;
                                    }
                                    else
                                    {
                                        if (nt == 0) detLayerTrackOutputFile << scatterings << "\t" << n << "\t" << neutronTrackCoordinatesFullSet.at(nt + 0) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << "\t" << layerRealisticallyHitted << endl;
                                    }
                                }
                                detLayerTrackOutputFile.close();
                            }
                        }
                    }

                    //this limit is mainly meant to cancel calculations of very long tracks in environments with poor absorbers.
                    //Typical limits of 300 scatterings are well above the typical number of 50
                    if (scatterings > maxScatterings)
                    {
                        absorbBreak = true;
                        killedThermal++;
                        //break;
                    }


                    //here the calculation of the scattering process starts
                    if ((wwRange < length) || (reflectiveBoundaries || periodicBoundaries))
                    {
                        foundSomething = false;
                        domainWallHit = false;
                        differentMaterialHit = false;

                        lengthZ = fabs(wwRange * cosTheta);
                        currentG = g;

                        if ((scatteredThisLayer) || (continuedThisLayer))
                        {
                            z0Here = z0;
                            xHere = x;
                            yHere = y;
                        }
                        else
                        {
                            if (theta < piHalf)
                            {
                                z0Here = tempz;
                                xHere = xt;
                                yHere = yt;
                            }
                            else
                            {
                                z0Here = tempzEnd;
                                xHere = xtEnd;
                                yHere = ytEnd;
                            }
                        }

                        if (useImage)
                        {
                            if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2) || (inputPics2[currentG] == 1) || (inputPics2[currentG] == 2) || (inputPics3[currentG] == 1) || (inputPics3[currentG] == 2))
                            {
                                trackMetricFactor = trackMetricFactorModifier * fabs(1000. * matrixMetricFactor / tanTheta);
                                if (trackMetricFactor < 0.0005 * matrixMetricFactor) trackMetricFactor = 0.0005 * matrixMetricFactor;
                                if (trackMetricFactor > 1e9 * matrixMetricFactor) trackMetricFactor = 1e9 * matrixMetricFactor;

                                ztrCounter = 0;

                                continueTracking = true;
                                differentMaterial1 = false;
                                differentMaterial2 = false;
                                differentMaterial3 = false;
                                oneLastCheck = true;

                                for (double ztr = z0Here; continueTracking;)
                                {
                                    if (ztrCounter > 10000) continueTracking = false;

                                    if (ztr == z0Here)
                                    {
                                        if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2))     startMaterial = actualMaterial;
                                        if ((inputPics2[currentG] == 1) || (inputPics2[currentG] == 2))   startMaterial2 = actualMaterial2;
                                        if ((inputPics3[currentG] == 1) || (inputPics3[currentG] == 2))   startMaterial3 = actualMaterial3;
                                    }

                                    if (theta < piHalf)
                                    {
                                        ztr = ztr + trackMetricFactor;
                                        if (ztr > tempzEnd)
                                        {
                                            continueTracking = false;
                                            ztr = tempzEnd;
                                        }
                                    }
                                    else
                                    {
                                        ztr = ztr - trackMetricFactor;
                                        if (ztr < tempz)
                                        {
                                            continueTracking = false;
                                            ztr = tempz;
                                        }
                                    }

                                    if (fabs(z0Here - ztr) > lengthZ)
                                    {
                                        if ((trackMetricFactor > lengthZ) && (oneLastCheck))
                                        {
                                            if (theta < piHalf) ztr = ztr - trackMetricFactor + 0.8 * lengthZ;
                                            else ztr = ztr + trackMetricFactor - 0.8 * lengthZ;

                                            oneLastCheck = false;
                                        }
                                        else
                                        {
                                            continueTracking = false;
                                        }
                                    }

                                    if ((scatteredThisLayer) || (continuedThisLayer))
                                    {
                                        xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + x;
                                        yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + y;
                                    }
                                    else
                                    {
                                        if (theta > piHalf)
                                        {
                                            xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + xtEnd;
                                            yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + ytEnd;
                                        }
                                        else
                                        {
                                            xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + xt;
                                            yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + yt;
                                        }
                                    }

                                    if (continueTracking)
                                    {
                                        matrixX = (xTrack * 0.001 - matrixStartX) / matrixMetricFactor;
                                        matrixY = (yTrack * 0.001 - matrixStartY) / matrixMetricFactor;

                                        if (matrixX > inputMatrixPixels - 1) matrixX = inputMatrixPixels - 1; if (matrixY > inputMatrixPixels - 1) matrixY = inputMatrixPixels - 1;
                                        if (matrixX < 0) matrixX = 0; if (matrixY < 0) matrixY = 0;

                                        if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2)) { selectedMaterial = (inputPicVector.at(currentG))(matrixX, matrixY); if ((selectedMaterial != actualMaterial) || (selectedMaterial != startMaterial)) differentMaterial1 = true; else differentMaterial1 = false; }
                                        if ((inputPics2[currentG] == 1) || (inputPics2[currentG] == 2)) { selectedMaterial2 = (inputPicVector2.at(currentG))(matrixX, matrixY); if ((selectedMaterial2 != actualMaterial2) || (selectedMaterial2 != startMaterial2)) differentMaterial2 = true; else differentMaterial2 = false; }
                                        if ((inputPics3[currentG] == 1) || (inputPics3[currentG] == 2)) { selectedMaterial3 = (inputPicVector3.at(currentG))(matrixX, matrixY); if ((selectedMaterial3 != actualMaterial3) || (selectedMaterial3 != startMaterial3)) differentMaterial3 = true; else differentMaterial3 = false; }
                                    }

                                    if (((differentMaterial1) || (differentMaterial2) || (differentMaterial3)) && (continueTracking))
                                    {
                                        //scatteredThisLayer = true;
                                        //continuedThisLayer = true;

                                        continueTracking = false;
                                        foundSomething = true;
                                        differentMaterialHit = true;
                                        zTrack = ztr;
                                    }

                                    ztrCounter++;
                                }
                            }
                        }

                        if (reflectiveBoundaries || periodicBoundaries)
                        {
                            sideNumber = 0;

                            //if (!((inputPics[currentG] == 1)||(inputPics[currentG] == 2)))
                            if (true)
                            {
                                if (!foundSomething)
                                {
                                    vectorDirFactor = (-0.5 * squareDim - xHere) / (cosPhi * tanTheta);
                                    zPlane = z0Here + vectorDirFactor;
                                    yPlane = yHere + (zPlane - z0Here) * sinPhi * tanTheta;
                                    if ((fabs(yPlane) < 0.5 * squareDim) && (cosPhi < 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                    {
                                        domainWallHit = true; foundSomething = true; sideNumber = 3;
                                        if (reflectiveBoundaries) xPlane = -squareDim * 0.5;
                                        if (periodicBoundaries) xPlane = squareDim * 0.5;
                                    }
                                }
                                if (!foundSomething)
                                {
                                    vectorDirFactor = (0.5 * squareDim - xHere) / (cosPhi * tanTheta);
                                    zPlane = z0Here + vectorDirFactor;
                                    yPlane = yHere + (zPlane - z0Here) * sinPhi * tanTheta;
                                    if ((fabs(yPlane) < 0.5 * squareDim) && (cosPhi > 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                    {
                                        domainWallHit = true; foundSomething = true; sideNumber = 1;
                                        if (reflectiveBoundaries) xPlane = squareDim * 0.5;
                                        if (periodicBoundaries) xPlane = -squareDim * 0.5;
                                    }
                                }
                                if (!foundSomething)
                                {
                                    vectorDirFactor = (0.5 * squareDim - yHere) / (sinPhi * tanTheta);
                                    zPlane = z0Here + vectorDirFactor;
                                    xPlane = xHere + (zPlane - z0Here) * cosPhi * tanTheta;
                                    if ((fabs(xPlane) < 0.5 * squareDim) && (sinPhi > 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                    {
                                        domainWallHit = true;  foundSomething = true; sideNumber = 2;
                                        if (reflectiveBoundaries) yPlane = squareDim * 0.5;
                                        if (periodicBoundaries) yPlane = -squareDim * 0.5;
                                    }
                                }
                                if (!foundSomething)
                                {
                                    vectorDirFactor = (-0.5 * squareDim - yHere) / (sinPhi * tanTheta);
                                    zPlane = z0Here + vectorDirFactor;
                                    xPlane = xHere + (zPlane - z0Here) * cosPhi * tanTheta;
                                    if ((fabs(xPlane) < 0.5 * squareDim) && (sinPhi < 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                    {
                                        domainWallHit = true;  foundSomething = true; sideNumber = 4;
                                        if (reflectiveBoundaries) yPlane = -squareDim * 0.5;
                                        if (periodicBoundaries) yPlane = squareDim * 0.5;
                                    }
                                }
                                if (fabs(zPlane - z0Here) > lengthZ)
                                {
                                    if (domainWallHit) { domainWallHit = false;  foundSomething = false; }
                                }
                            }
                        }

                        if (foundSomething)
                        {
                            continuedThisLayer = true;

                            if (domainWallHit)
                            {
                                z0 = zPlane;
                                x = xPlane;
                                y = yPlane;

                                if (reflectiveBoundaries)
                                {
                                    if ((sideNumber == 1) || (sideNumber == 3)) phi = pi - phi;
                                    else phi = 2. * pi - phi;
                                }
                                else
                                {
                                }
                            }

                            if (differentMaterialHit) //useimage
                            {
                                x = xTrack;
                                y = yTrack;
                                z0 = zTrack;
                            }

                            if (!noTrackRecording)
                            {
                                if ((scatteredThisLayer) || (continuedThisLayer))
                                {
                                    if ((g == detectorLayer) || (trackAllLayers))
                                    {
                                        neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                        neutronTrackCoordinates.push_back(x);
                                        neutronTrackCoordinates.push_back(y);
                                        neutronTrackCoordinates.push_back(z0);
                                        neutronTrackCoordinates.push_back(theta);
                                        neutronTrackCoordinates.push_back(phi);
                                        neutronTrackCoordinates.push_back(energy);
                                        neutronTrackCoordinates.push_back(timeNtr);
                                        neutronTrackCoordinates.push_back(-10);
                                        neutronTrackCoordinates.push_back(-10);
                                    }
                                    if (showDensityTrackMapSide)
                                    {
                                        neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                        neutronTrackCoordinates2.push_back(x);
                                        neutronTrackCoordinates2.push_back(y);
                                        neutronTrackCoordinates2.push_back(z0);
                                        neutronTrackCoordinates2.push_back(theta);
                                        neutronTrackCoordinates2.push_back(phi);
                                        neutronTrackCoordinates2.push_back(energy);
                                    }
                                }
                            }
                        }
                        //else  {continuedThisLayer = false; }
                    }


                    // neutron is absorbed or scattered
                    if (wwRange < length)
                    {
                        if (!foundSomething)
                        {
                            z0Alt = z0;	energyOld = energy; phiAlt = phi; thetaAlt = theta;	xAlt = x; yAlt = y;
                            if (reverseDir) reverseDirAlt = true; else reverseDirAlt = false;

                            if ((!scatteredThisLayer) && (!continuedThisLayer))
                            {
                                if (theta < piHalf) { z0 = geometries.at(g)[4] + fabs(wwRange * cosTheta); }
                                else { z0 = geometries.at(g)[4] + geometries.at(g)[5] - fabs(wwRange * cosTheta); }
                            }
                            else
                            {
                                if (theta < piHalf) { z0 = z0 + fabs(wwRange * cosTheta); }
                                else { z0 = z0 - fabs(wwRange * cosTheta); }
                            }

                            if (z0max < z0) z0max = z0; //maximum penetration depth

                            //set new coordiantes for trajectory vector
                            x = cosPhi * fabs(tanTheta * (z0Alt - z0)) + x;
                            y = sinPhi * fabs(tanTheta * (z0Alt - z0)) + y;

                            xLastScattered = x; yLastScattered = y; zLastScattered = z0;

                            if (sAbsorbed)
                            {
                                absMaterialNo = material;

                                element = -1;

                                if (element < 0)
                                {
                                    tempInt = allAbsorptionMaterials.size() - 1;
                                    temp = r.Rndm();
                                    for (int k33 = 0; k33 < allAbsorptionMaterials.size(); k33++) { if (temp < allAbsorptionMaterials.at(k33)) { tempInt = k33;   k33 = allAbsorptionMaterials.size(); } }
                                    element = allAbsorptionMaterialsElements.at(tempInt);
                                }

                                switch (element)
                                {
                                case 1: elementID = nH1; break;
                                case 3: elementID = nHe3; break;
                                case 10: elementID = nB10; break;
                                case 11: elementID = nB11; break;
                                case 12: elementID = nC12; break;
                                case 14: elementID = nN14; break;
                                case 16: elementID = nO16; break;
                                case 19: elementID = nF19; break;
                                case 23: elementID = nNa23; break;
                                case 27: elementID = nAl27; break;
                                case 28: elementID = nSi28; break;
                                case 32: elementID = nS32; break;
                                case 35: elementID = nCl35; break;
                                case 39: elementID = nK39; break;
                                case 40: elementID = nAr40; break;
                                case 48: elementID = nTi48; break;
                                case 52: elementID = nCr52; break;
                                case 53: elementID = nCr53; break;
                                case 55: elementID = nMn55; break;
                                case 56: elementID = nFe56; break;
                                case 58: elementID = nNi58; break;
                                case 155: elementID = nGd155; break;
                                case 157: elementID = nGd157; break;
                                case 206: elementID = nPb206; break;
                                case 207: elementID = nPb207; break;
                                case 208: elementID = nPb208; break;
                                default: elementID = 1;
                                }

                                if (!dontFillunecessaryPlots)
                                {
                                    absDist->Fill(currentlayer);
                                    if ((material == 27) || (material == 28))
                                    {
                                        absEnergy->Fill(energy);
                                    }
                                    absMaterial->Fill(material);
                                    //depthEnergyAbsorbed->Fill(currentlayer,energyInitial);

                                    if (material == 28)
                                    {
                                        //temp = r.Rndm();
                                        //if (temp*asAll > 3.*asF) element = 10; else element = 19;
                                        if (r.Rndm() > deadMaterialFactor)
                                        { //account for absorption with the boron itself
                                            absElement->Fill(element);
                                        }
                                        if (element == 10) absElementL1->Fill(currentlayer + 1);
                                    }
                                    else
                                    {
                                        if (element == 1) absElement->Fill(element);
                                        else absElement->Fill(element);
                                    }
                                    if (absMaterialNo == 27) absElement->Fill(3);
                                    if (absMaterialNo == 25) absElementL2->Fill(currentlayer); //Moderator
                                    if (absMaterialNo == 29) absElementL3->Fill(currentlayer); //absorber
                                }

                                absDepth = z0;

                                if (!dontFillunecessaryPlots)
                                {
                                    depthEnergyAbsorbed2->Fill(absDepth, energyAtInterface);

                                    if ((thermalizedLayer >= 0)) depthThermalized->Fill(absDepth, energyAtInterface);
                                }

                                if (drawSingleNeutronPropagation)
                                {
                                    if (energy < 20)	nSpeed = 3.9560339 / sqrt(81.81 / energy / 1e9) * 1000. * 1000.;
                                    else 	nSpeed = 3.9560339 / sqrt(81.81 / 20. / 1e9) * 1000. * 1000.;
                                    if (energy < nSpeedEnergyCutoff) nSpeed = 3.9560339 / sqrt(81.81 / nSpeedEnergyCutoff / 1e9) * 1000. * 1000.;

                                    //timeTemp = wwRange/nSpeed;
                                    timeTemp = fabs((z0Alt - z0) / cos(thetaAlt) / nSpeed);
                                    tt = 0;
                                    while (timeFrame * tt < timeTemp - timeTrans)
                                    {
                                        trackTemp = nSpeed * (timeFrame * tt + timeTrans);

                                        if ((!scatteredThisLayer) && (false))
                                        {
                                            if (!reverseDir) tempz = geometries.at(g)[4] + trackTemp * cos(theta);
                                            else   tempz = geometries.at(g)[4] + geometries.at(g)[5] - fabs(trackTemp * cos(theta));
                                        }
                                        else
                                        {
                                            if (!reverseDir)  tempz = z0Alt + fabs(trackTemp * cos(theta));
                                            else  tempz = z0Alt - fabs(trackTemp * cos(theta));
                                        }

                                        float* coordinates = new float[5];
                                        coordinates[0] = cos(phi) * TMath::Abs(tan(theta) * (tempz - z0Alt)) + xAlt;
                                        coordinates[1] = sin(phi) * TMath::Abs(tan(theta) * (tempz - z0Alt)) + yAlt;
                                        coordinates[2] = tempz;  coordinates[3] = energy; coordinates[4] = 0;
                                        if (neutronCoordinates.size() < numberOfFrames + 1)
                                        {
                                            neutronCoordinates.push_back(coordinates);
                                        }
                                        tt++;
                                        if (tt > trackIterationCutoff) break;
                                    }
                                }

                                if (drawSingleNeutronGraphs)
                                {
                                    if (!reverseDir) tempz = z0 + fabs(wwRange * cos(theta));
                                    else tempz = z0 - fabs(wwRange * cos(theta));
                                    x = cos(phi) * TMath::Abs(tan(theta) * (tempz - z0)) + x;
                                    neutronPath->SetPoint(graphCounter, x/1000.,-tempz); graphCounter++; neutronPath->SetTitle(("Neutron No. "+castLongToString(n)+" absorbed").c_str());
                                    graphCounterMG++;
                                    graphs.at(graphN)->SetPoint(graphCounterMG,x/1000., -tempz); graphs.at(graphN)->SetTitle(("Neutron No. "+castLongToString(n)+" absorbed").c_str());
                                }

                                if (!noTrackRecording)
                                {
                                    if ((g == detectorLayer) || (trackAllLayers))
                                    {
                                        neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                        neutronTrackCoordinates.push_back(x);
                                        neutronTrackCoordinates.push_back(y);
                                        neutronTrackCoordinates.push_back(z0);
                                        neutronTrackCoordinates.push_back(theta);
                                        neutronTrackCoordinates.push_back(phi);
                                        neutronTrackCoordinates.push_back(energy);
                                        neutronTrackCoordinates.push_back(timeNtr);
                                        neutronTrackCoordinates.push_back(elementID);
                                        neutronTrackCoordinates.push_back(material);
                                    }
                                    if (showDensityTrackMapSide)
                                    {
                                        neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                        neutronTrackCoordinates2.push_back(x);
                                        neutronTrackCoordinates2.push_back(y);
                                        neutronTrackCoordinates2.push_back(z0);
                                        neutronTrackCoordinates2.push_back(theta);
                                        neutronTrackCoordinates2.push_back(phi);
                                        neutronTrackCoordinates2.push_back(energy);
                                    }
                                }
                                absorbBreak = true;
                            }

                            if ((sAbsorbed) && (element > 1) && (energy > 15.5))
                            {
                                if (evaporationFactor < 1.0)
                                {
                                    tmpEvapo = r.Rndm();

                                    if (tmpEvapo < evaporationFactor)
                                    {
                                        scatteredEvaporating = true; // sAbsorbed = false;
                                    }
                                    else scatteredEvaporating = false;
                                }

                                if (evaporationAdditions > 0)
                                {
                                    scatteredEvaporating = true;
                                    //sAbsorbed = false;
                                }

                                if ((energy > 15.5) && (useHECascadeModel))
                                {
                                    //if((scatteredEvaporating||scatteredElastic)&&(r.Rndm()<.4))
                                    if (r.Rndm() < heEvaporationFactor) //
                                    {
                                        nEvapoVector.reserve(nEvapoVector.size() + 9);

                                        //nEvapoVector.clear();
                                        nEvapoVector.push_back(x); //put in new coordinates of the particle
                                        nEvapoVector.push_back(y);
                                        nEvapoVector.push_back(z0);
                                        nEvapoVector.push_back(theta + ((2. * r.Rndm() - 1.) * 0.02)); // angular distribution for the HE cascade model
                                        nEvapoVector.push_back(phi);
                                        nEvapoVector.push_back(g);
                                        nEvapoVector.push_back((.95 + r.Rndm() * 0.05) * energy); // energy distribution for the HE cascade model
                                        nEvapoVector.push_back(timeNtr);
                                        if (hasBeenInSoil) nEvapoVector.push_back(1); else nEvapoVector.push_back(0);
                                        //n++;
                                    }
                                }
                            }

                            if (scatteredEvaporating)
                            {
                                if (evaporationAdditions > 0)
                                {
                                    //nEvapoVector.clear();
                                    for (int ee = 0; ee < evaporationAdditions; ee++)
                                    {
                                        nEvapoVector.reserve(nEvapoVector.size() + 9);

                                        nEvapoVector.push_back(x);
                                        nEvapoVector.push_back(y);
                                        nEvapoVector.push_back(z0);
                                        //nEvapoVector.push_back(r.Rndm()*pi);  //changed from TMath::ACos(2.*r.Rndm()-1.) 28.11.2018
                                        nEvapoVector.push_back(TMath::ACos(2. * r.Rndm() - 1.)); // changed back to TMath::ACos(2.*r.Rndm()-1.) 13.12.2018
                                        nEvapoVector.push_back((2. * r.Rndm() - 1.) * pi);
                                        nEvapoVector.push_back(g);
                                        nEvapoVector.push_back(getEvaporationEnergy(2.e6, &r)); //changed from 1.5 21.10.2017
                                        nEvapoVector.push_back(timeNtr);
                                        if (hasBeenInSoil) nEvapoVector.push_back(1); else nEvapoVector.push_back(0);
                                        evaporationCounter++;
                                    }
                                }

                                if ((tmpEvapo < evaporationFactor) && (true))
                                {
                                    nEvapoVector.reserve(nEvapoVector.size() + 9);

                                    nEvapoVector.push_back(x);
                                    nEvapoVector.push_back(y);
                                    nEvapoVector.push_back(z0);
                                    //nEvapoVector.push_back(r.Rndm()*pi);  //changed from TMath::ACos(2.*r.Rndm()-1.) 28.11.2018
                                    nEvapoVector.push_back(TMath::ACos(2. * r.Rndm() - 1.)); // changed back to TMath::ACos(2.*r.Rndm()-1.) 13.12.2018
                                    nEvapoVector.push_back((2. * r.Rndm() - 1.) * pi);
                                    nEvapoVector.push_back(g);
                                    nEvapoVector.push_back(getEvaporationEnergy(2.e6, &r)); //changed from 1.5 21.10.2017
                                    nEvapoVector.push_back(timeNtr);
                                    if (hasBeenInSoil) nEvapoVector.push_back(1); else nEvapoVector.push_back(0);
                                    evaporationCounter++;
                                }

                                scatteredEvaporating = false;
                            }

                            //check if scattered incoherently
                            if (scatteredInelastic || scatteredElastic)
                            {
                                scatteredThisLayer = true;

                                if (drawSingleNeutronGraphs)
                                {
                                    neutronPath->SetPoint(graphCounter, x * 0.001, -z0); graphCounter++;

                                    graphCounterMG++;
                                    graphs.at(graphN)->SetPoint(graphCounterMG, x * 0.001, -z0);
                                }

                                if (!dontFillunecessaryPlots) scatMaterial->Fill(material);

                                if (!sAbsorbed)
                                {
                                    if (!scatteredInelastic)  // that means scatteredElastic
                                    {
                                        temp = r.Rndm();
                                        //temp2 = r.Rndm();
                                        element = -1;

                                        switch (material)
                                        {
                                        case 9:  if (temp * csW > 2. * csH) element = 16; else element = 1; break;
                                            //case 10: if (temp*cs > 2.*csN*0.78) element = 16; else element = 14; break;
                                            //case 11: if (temp2*rLuft/14.44*cs > rLuftWater/18.*csW) {if (temp*cs > 2.*csN*0.78) element = 16; else element = 14;} else {if (temp*csW > 2.*csH) element = 16; else element = 1;} break;
                                        case 12: if (temp * cs > csSi) element = 16; else element = 28; break;
                                        case 13: if (temp * cs > 2. * csAl) element = 16; else element = 27; break;
                                        case 25: if (temp * cs > 2. * csH) element = 12; else element = 1; break;
                                        case 26: element = 27; break;
                                        case 27: element = 3; break;
                                            //case 31: if (temp*cs > 4.*csH) element = 12; else element = 1; break;
                                            //case 34: if (temp*cs > 23.*csH) element = 12; else element = 1; break;
                                        }

                                        if (material == 18)
                                        {
                                            tempInt = allBoden.size() - 1;
                                            temp = r.Rndm();
                                            for (int k33 = 0; k33 < allBoden.size(); k33++)
                                            {
                                                if (temp < allBoden.at(k33)) { tempInt = k33;   k33 = allBoden.size(); }
                                            }
                                            element = allBodenElements.at(tempInt);
                                        }
                                        if (material == 19)
                                        {
                                            tempInt = allBoden.size() - 1;
                                            temp = r.Rndm();
                                            for (int k33 = 0; k33 < allBoden.size(); k33++)
                                            {
                                                if (temp < allBoden.at(k33)) { tempInt = k33;   k33 = allBoden.size(); }
                                            }
                                            element = allBodenElements.at(tempInt);
                                        }
                                        if (material == 20)
                                        {
                                            tempInt = allBoden.size() - 1;
                                            temp = r.Rndm();
                                            for (int k33 = 0; k33 < allBoden.size(); k33++)
                                            {
                                                if (temp < allBoden.at(k33)) { tempInt = k33;   k33 = allBoden.size(); }
                                            }
                                            element = allBodenElements.at(tempInt);
                                        }
                                        if (material == 21)
                                        {
                                            tempInt = allPlants.size() - 1;
                                            temp = r.Rndm();
                                            for (int k33 = 0; k33 < allPlants.size(); k33++)
                                            {
                                                if (temp < allPlants.at(k33)) { tempInt = k33;   k33 = allPlants.size(); }
                                            }
                                            element = allPlantsElements.at(tempInt);
                                        }

                                        if (element < 0)
                                        {
                                            tempInt = allMaterials.size() - 1;
                                            temp = r.Rndm();
                                            for (int k33 = 0; k33 < allMaterials.size(); k33++) { if (temp < allMaterials.at(k33)) { tempInt = k33;   k33 = allMaterials.size(); } }
                                            element = allMaterialsElements.at(tempInt);
                                        }

                                        switch (element)
                                        {
                                        case 1: weight = wHydrogen; angularProb = &angularH;  elementID = nH1; break;
                                        case 3: weight = wHe3; angularProb = &angularHe3; elementID = nHe3; break;
                                        case 10: weight = wB10; angularProb = &angularB10; elementID = nB10; break;
                                        case 11: weight = wB11; angularProb = &angularB11; elementID = nB11; break;
                                        case 12: weight = wCarbon; angularProb = &angularC;  elementID = nC12; break;
                                        case 14: weight = wNitrogen; angularProb = &angularN;  elementID = nN14; break;
                                        case 16: weight = wOxygen;  angularProb = &angularO; elementID = nO16; break;
                                        case 19: weight = wF;  angularProb = &angularF; elementID = nF19; break;
                                        case 23: weight = wNa;  angularProb = &angularNa; elementID = nNa23; break;
                                        case 27: weight = wAl;  angularProb = &angularAl; elementID = nAl27; break;
                                        case 28: weight = wSilicon;  angularProb = &angularSi; elementID = nSi28; break;
                                        case 32: weight = wSulfur;  angularProb = &angularS; elementID = nS32; break;
                                        case 35: weight = wCl;  angularProb = &angularCl35; elementID = nCl35; break;
                                        case 39: weight = wKalium;  angularProb = &angularK39; elementID = nK39; break;
                                        case 40: weight = wAr;  angularProb = &angularAr; elementID = nAr40; break;
                                        case 48: weight = wTi;  angularProb = &angularTi48; elementID = nTi48; break;
                                        case 52: weight = wCr52;  angularProb = &angularCr52; elementID = nCr52; break;
                                        case 53: weight = wCr53;  angularProb = &angularCr53; elementID = nCr53; break;
                                        case 55: weight = wMn55;  angularProb = &angularMn55; elementID = nMn55; break;
                                        case 56: weight = wFe;  angularProb = &angularFe; elementID = nFe56; break;
                                        case 58: weight = wNi;  angularProb = &angularNi58;  elementID = nNi58; break;
                                        case 155: weight = wGd155;  angularProb = &angularGd155; elementID = nGd155; break;
                                        case 157: weight = wGd157;  angularProb = &angularGd157; elementID = nGd157; break;
                                        case 206: weight = wPb206;  angularProb = &angularPb206; elementID = nPb206; break;
                                        case 207: weight = wPb207;  angularProb = &angularPb207; elementID = nPb207; break;
                                        case 208: weight = wPb208;  angularProb = &angularPb208; elementID = nPb208; break;
                                        default: weight = 15; angularProb = &angularSi; elementID = 1;
                                        }
                                    }
                                    else // scatteredInelastic
                                    {
                                        tempInt = allInelastics.size() - 1;
                                        temp = r.Rndm();
                                        for (int k33 = 0; k33 < allInelastics.size(); k33++)
                                        {
                                            if (temp < allInelastics.at(k33)) { tempInt = k33;   k33 = allInelastics.size(); }
                                        }
                                        element = allInelasticElements.at(tempInt);

                                        switch (element)
                                        {
                                        case 1: weight = wHydrogen; elementID = nH1; break;
                                        case 3: weight = wHe3; elementID = nHe3; break;
                                        case 10: weight = wB10; elementID = nB10; break;
                                        case 11: weight = wB11; elementID = nB11; break;
                                        case 12: weight = wCarbon; elementID = nC12; break;
                                        case 14: weight = wNitrogen; elementID = nN14; break;
                                        case 16: weight = wOxygen; elementID = nO16; break;
                                        case 19: weight = wF; elementID = nF19; break;
                                        case 23: weight = wNa; elementID = nNa23; break;
                                        case 27: weight = wAl; elementID = nAl27; break;
                                        case 28: weight = wSilicon; elementID = nSi28; break;
                                        case 32: weight = wSulfur; elementID = nS32; break;
                                        case 35: weight = wCl; elementID = nCl35; break;
                                        case 39: weight = wKalium; elementID = nK39; break;
                                        case 40: weight = wAr; elementID = nAr40; break;
                                        case 48: weight = wTi; elementID = nTi48; break;
                                        case 52: weight = wCr52; elementID = nCr52; break;
                                        case 53: weight = wCr53; elementID = nCr53; break;
                                        case 55: weight = wMn55; elementID = nMn55; break;
                                        case 56: weight = wFe; elementID = nFe56; break;
                                        case 58: weight = wNi; elementID = nNi58; break;
                                        case 155: weight = wGd155; elementID = nGd155; break;
                                        case 157: weight = wGd157; elementID = nGd157; break;
                                        case 206: weight = wPb206; elementID = nPb206; break;
                                        case 207: weight = wPb207; elementID = nPb207; break;
                                        case 208: weight = wPb208; elementID = nPb208; break;
                                        default: weight = 15; elementID = 1;
                                        }

                                        inelasticEnergyLoss = inelasticEnergyLossVec.at(tempInt);
                                        angularProb = allInelasticAngulars.at(tempInt);
                                    }
                                }


                                //////////////////////////////////////////////Set new angles
                                //scattering without thermal energy of the target

                                if ((scatteredElastic) || (scatteredInelastic))
                                //elastic scattering
                                {
                                    if (!dontFillunecessaryPlots) scatElement->Fill(element);

                                    if (!scatteredInelastic) // that means scatteredElastic
                                    {
                                        if (((energy >= 20.5) && (element > 1)) && (!(((element == 16) || (element == 140) || (element == 39) || (element == 48) || (element == 10) || (element == 155) || (element == 157)) && (energy < 30.5))))
                                        { // elements and energies for which tabulated energy distributions exist
                                            switch (element)
                                            {
                                            case 12: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyC.at(0), angularHighEnergyC.at(1), energy * 1e6, r.Rndm()); break;
                                            case 14: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyN.at(0), angularHighEnergyN.at(1), energy * 1e6, r.Rndm()); break;
                                            case 16: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyO.at(0), angularHighEnergyO.at(1), energy * 1e6, r.Rndm()); break;
                                            case 19: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyF.at(0), angularHighEnergyF.at(1), energy * 1e6, r.Rndm()); break;
                                            case 23: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyNa.at(0), angularHighEnergyNa.at(1), energy * 1e6, r.Rndm()); break;
                                            case 27: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyAl.at(0), angularHighEnergyAl.at(1), energy * 1e6, r.Rndm()); break;
                                            case 28: thetaFromTab = getHighEnergyCosTheta(angularHighEnergySi.at(0), angularHighEnergySi.at(1), energy * 1e6, r.Rndm()); break;
                                            case 35: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyCl.at(0), angularHighEnergyCl.at(1), energy * 1e6, r.Rndm()); break;
                                            case 40: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyAr.at(0), angularHighEnergyAr.at(1), energy * 1e6, r.Rndm()); break;
                                            case 55: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyMn.at(0), angularHighEnergyMn.at(1), energy * 1e6, r.Rndm()); break;
                                            case 56: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyFe.at(0), angularHighEnergyFe.at(1), energy * 1e6, r.Rndm()); break;
                                            case 206: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyPb206.at(0), angularHighEnergyPb206.at(1), energy * 1e6, r.Rndm()); break;
                                            case 207: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyPb207.at(0), angularHighEnergyPb207.at(1), energy * 1e6, r.Rndm()); break;
                                            case 208: thetaFromTab = getHighEnergyCosTheta(angularHighEnergyPb208.at(0), angularHighEnergyPb208.at(1), energy * 1e6, r.Rndm()); break;
                                            default:  thetaFromTab = getHighEnergyCosTheta(angularHighEnergySi.at(0), angularHighEnergySi.at(1), energy * 1e6, r.Rndm()); break;
                                            }

                                            thetaCMS = TMath::ACos(thetaFromTab);
                                            alpha = TMath::Power((1. * weight - 1.) / (1. * weight + 1.), 2);
                                            energyNew = 0.5 * (1. + alpha + (1. - alpha) * cos(thetaCMS)) * energyOld;
                                        }
                                        else
                                        {
                                            if (((energy > 1e-2) && (element > 2)) || ((energy > 1) && (element == 1)))
                                            {
                                                TF1* angularProbFunc = calcMeanAngularDistribution(*angularProb, energy * 1e6);

                                                funcMax = angularProbFunc->Eval(0) * 1.2; //changed 10.02.2018 // that is actually an approximation, there are rare cases where this would underestimate the function.

                                                //funcMax = 1;
                                                //if (funcMax>20*funcMin) useCumulativeCalc = true; else useCumulativeCalc = false;
                                                useCumulativeCalc = false;

                                                if (!useCumulativeCalc)
                                                {
                                                    gotIt = false;
                                                    gotItCounter = 0;
                                                    while (!gotIt)
                                                    {
                                                        xRnd = TMath::ACos(2. * r.Rndm() - 1.);
                                                        yRnd = r.Rndm() * funcMax;

                                                        if (angularProbFunc->Eval(xRnd) > yRnd)
                                                        {
                                                            gotIt = true;
                                                        }
                                                        gotItCounter++;
                                                        if (gotItCounter > 1000) gotIt = true;
                                                    }
                                                    thetaCMS = xRnd;
                                                }
                                                else
                                                {
                                                    thetaCMS = pi / 180. * getAngleFromCumulativeFunction(angularProbFunc, 0, pi, &r);
                                                }
                                                if (angularProbFunc != 0x0) delete angularProbFunc;

                                                alpha = TMath::Power((1. * weight - 1.) / (1. * weight + 1.), 2);
                                                energyNew = 0.5 * (1. + alpha + (1. - alpha) * cos(thetaCMS)) * energyOld;
                                            }
                                            //scattering with thermal energy of the target
                                            else
                                            {
                                                //only computational speed issues
                                                float thermalcutoff = 0.17115e-6; //changed on 15.01.2018
                                                if ((energy > thermalcutoff) || (noThermalRegime))
                                                {
                                                    thetaCMS = TMath::ACos(2. * r.Rndm() - 1.);
                                                    alpha = TMath::Power((weight - 1.) / (weight + 1.), 2);
                                                    energyNew = 0.5 * (1. + alpha + (1. - alpha) * cos(thetaCMS)) * energyOld;
                                                }
                                                else
                                                {
                                                    //these are gases
                                                    if ((material == 10) || (material == 11) || (material == 27) || (material == 28) || (material == 33))
                                                    {
                                                        //weightThermal = weight;
                                                        vNeutron = getThermalPDF(energyOld, weightThermal, 300., &r);

                                                        thetaCMS = TMath::ACos(vNeutron.at(1));
                                                        //energyNew = vNeutron.at(0)/1e6;
                                                        alpha = TMath::Power((weightThermal - 1.) / (weightThermal + 1.), 2);
                                                        energyNew = 0.5 * (1. + alpha + (1. - alpha) * cos(thetaCMS)) * vNeutron.at(0) / 1e6;

                                                        //if ((energyNew > thermalcutoff)&&(!noThermalRegime)) energyNew = 0.5*(1.+alpha+(1.-alpha)*TMath::Cos(thetaCMS))*vNeutron.at(0)/1e6;

                                                        intoThermalization = true;
                                                        if (energyNew > thermalcutoff) intoThermalization = false;
                                                        if (energyNew < 1e-11) energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                    }
                                                    else
                                                    {
                                                        thetaCMS = TMath::ACos(2. * r.Rndm() - 1.);
                                                        if (false)
                                                        {
                                                            alpha = TMath::Power((weightThermal - 1.) / (weightThermal + 1.), 2);
                                                            energyNew = 0.5 * (1. + alpha + (1. - alpha) * cos(thetaCMS)) * energyOld;
                                                        }
                                                        else
                                                        {
                                                            //energyNew = getThermalEnergyFromSource(spectrumMaxwellPhiTempModMod, &r)/1e6*1.0; // deprecated
                                                            energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                            if (energyNew > thermalcutoff) energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                            if (energyNew > thermalcutoff) energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                            if (energyNew > thermalcutoff) energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                            if (energyNew < 1e-11) energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                                        }

                                                        intoThermalization = true;
                                                    }
                                                    // liquid velocity distribution?  m^3v^5/(kT)^3 * EXP [-mv^2/kT]
                                                    // Reference https://www.physicsforums.com/threads/liquid-molecular-velocity-distribution.231145/
                                                    // energy distribution 4E^2/(kT)^3 * EXP [-2E/kT]
                                                }
                                            }
                                        }
                                    }
                                    //inelastic
                                    else
                                    {
                                        TF1* angularProbFunc = calcMeanAngularDistribution(*angularProb, energy * 1e6);

                                        funcMax = angularProbFunc->Eval(0) * 1.5; //changed 10.02.2018 // that is actually an approximation, there are rare cases where this would underestimate the function.

                                        gotIt = false;
                                        gotItCounter = 0;
                                        while (!gotIt)
                                        {
                                            xRnd = TMath::ACos(2. * r.Rndm() - 1.);
                                            yRnd = r.Rndm() * funcMax;
                                            if (angularProbFunc->Eval(xRnd) > yRnd)
                                            {
                                                gotIt = true;
                                            }

                                            gotItCounter++;
                                            if (gotItCounter > 1000) gotIt = true;
                                        }
                                        thetaCMS = xRnd;

                                        energyNew = energyOld - inelasticEnergyLoss;
                                        if (energyNew < 1e-11)
                                        {
                                            energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);  //mostly due to duplicated inelastic cross sections
                                        }

                                        if (angularProbFunc != 0x0) delete angularProbFunc;
                                    }

                                    //both calculations are equivalent
                                    //thetaScat = thetaCMS/2. + TMath::ATan(sqrt(alpha)*TMath::Tan(thetaCMS/2.));

                                    if (energy < 0.1e-6)
                                    {
                                        thetaScat = TMath::ACos((1. + weightThermal * cos(thetaCMS)) / (sqrt(weightThermal * weightThermal + 1. + 2. * weightThermal * cos(thetaCMS))));
                                    }
                                    else
                                    {
                                        thetaScat = TMath::ACos((1. + weight * cos(thetaCMS)) / (sqrt(weight * weight + 1. + 2. * weight * cos(thetaCMS))));
                                    }

                                    phiScat = r.Rndm() * 2. * pi - pi;

                                    doNormalCalc = true;

                                    // just checks here for simple cases
                                    if (theta == 0)
                                    {
                                        theta = thetaScat;
                                        phi = phi + phiScat;
                                        doNormalCalc = false;
                                    }
                                    else
                                    {
                                        if (thetaScat == 0)
                                        {
                                            theta = theta;
                                            phi = phi + phiScat;
                                            doNormalCalc = false;
                                        }
                                        else
                                        {
                                            if (phiScat == 0)
                                            {
                                                theta = theta + thetaScat;
                                                phi = phi;
                                                doNormalCalc = false;
                                            }
                                        }
                                    }

                                    // thetas1->Fill(theta);

                                    if (doNormalCalc)
                                    {
                                        thetaCalc = theta;
                                        dcosTheta = cos(thetaScat);
                                        theta12 = cos(thetaCalc) * cos(thetaScat) + sin(thetaCalc) * sin(thetaScat) * cos(pi + 1 * phiScat);
                                        theta = TMath::ACos(theta12);

                                        if (phiScat > 0) phi = phi + TMath::ACos((dcosTheta - theta12 * cos(thetaCalc)) / (sin(theta12) * sin(thetaCalc)));
                                        else phi = phi - TMath::ACos((dcosTheta - theta12 * cos(thetaCalc)) / (sin(theta12) * sin(thetaCalc)));
                                    }

                                    if (theta > pi) { theta = 2. * pi - theta;  phi = phi + 1. * pi; }
                                    if (theta < 0) { theta = fabs(theta);  phi = phi + 1. * pi; }

                                    if (phi < 0) phi = phi + 2. * pi;
                                    if (phi > 2. * pi) phi = phi - 2. * pi;

                                    //thetas2->Fill(thetaScat,  1./sin(thetaScat));
                                    //thetas1->Fill(theta,  1./sin(theta));
                                } //end of else scatteredevaporating

                                //detectorPhi2->Fill(fabs(fmod(phi,2.*pi)));

                                //if (energyNew < 1.5e-7) energyNew = getThermalEnergy(phiTFunc, &r);
                                //if (((energyOld < 1.5e-7)&&(energyNew < 9.e-8))||(energyNew < 0) )

                                if (energyNew < 0) // that sould not happen (and normally does not)
                                {
                                    cout << "E";
                                    energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                    if (energyNew > 1.75e-7)  energyNew = getThermalEnergy(spectrumMaxwellPhiLinMod, &r);
                                }

                                energy = energyNew;

                                if (!dontFillunecessaryPlots)
                                {
                                    energyLoss->Fill(energyOld - energyNew);
                                }

                                //lambda = getLfromE(energyNew);

                                if ((thermalizedLayer < 0) && (energyNew < 5e-7)) { thermalizedLayer = currentlayer; thermalized = true; }

                                //if ((energyNew < 1.46e-6)&&(!slowedDown))	//Indium resonance
                                //if ((energyNew < 1.8e-6)&&(energyNew > 1.1e-6) ) //Indium resonance

                                if ((energyNew < 80e-9) && (energyNew > 1e-8))
                                {
                                    slowedDown = true;
                                    if (!dontFillunecessaryPlots)
                                    {
                                        //slowDownDepth->Fill(sqrt((z0-0.)*(z0-0.)+(x-xStart)*(x-xstart)+(y-yStart)*(y-yStart)));
                                    }
                                }

                                if (!dontFillunecessaryPlots) scatDist->Fill(currentlayer);

                                scatDepth->Fill(z0);
                                scatterLayer = currentlayer;
                                if (currentlayer > maxlayer) maxlayer = currentlayer;
                                if ((recordSubsurfaceScatterings) && (g >= groundLayer)) { subsurfaceScatterings.push_back(z0); }

                                if (theta > piHalf)
                                {
                                    if (debugOutput > 0) cout << "reverse ";
                                    reverseDir = true;
                                }
                                else { reverseDir = false; }

                                scatteredThisLayer = true;
                                scattered = true;
                                scatterings++;

                                if (drawSingleNeutronGraphs)
                                {
                                    scale = (log10(energyNew) + 7.7) / 12.; if (scale < 0) scale = 0.;
                                    rgbValues = getRGBfromHCL(getLinearC(scale, 0, 240), getScaledColorValue(scale, 0.4, 0.91, 0.9), getScaledColorValue(scale, 2.2, 0.8, 0.9));
                                    graphs.push_back(new TGraph());
                                    graphN++; graphCounterMG = 0;
                                    graphs.at(graphN)->SetMarkerStyle(20);
                                    if (element == 1) graphs.at(graphN)->SetMarkerColor(2);
                                    if (element == 14) graphs.at(graphN)->SetMarkerColor(8);
                                    if (element == 16) graphs.at(graphN)->SetMarkerColor(4);
                                    if (element == 27) graphs.at(graphN)->SetMarkerColor(16);
                                    if (element == 28) graphs.at(graphN)->SetMarkerColor(43);
                                    graphs.at(graphN)->SetLineColor(TColor::GetColor(rgbValues.at(0), rgbValues.at(1), rgbValues.at(2)));
                                    graphs.at(graphN)->SetPoint(graphCounterMG, x / 1000., -z0);
                                }

                                if (drawSingleNeutronPropagation)
                                {
                                    if (energyOld < 20)	nSpeed = 3.9560339 / sqrt(81.81 / energyOld / 1e9) * 1000. * 1000.;
                                    else 	nSpeed = 3.9560339 / sqrt(81.81 / 20. / 1e9) * 1000. * 1000.;
                                    if (energyOld < nSpeedEnergyCutoff) nSpeed = 3.9560339 / sqrt(81.81 / nSpeedEnergyCutoff / 1e9) * 1000. * 1000.;

                                    timeTemp = fabs((z0Alt - z0) / cos(thetaAlt) / nSpeed);

                                    tt = 0; trackTemp = 0;

                                    while (timeTrans + timeFrame * tt < timeTemp)
                                    {
                                        trackTemp = nSpeed * (timeFrame * tt + timeTrans);

                                        if ((!scatteredThisLayer) && (!stillFirstLayer) && (false))
                                        {
                                            if (!reverseDirAlt) tempz = geometries.at(g)[4] + trackTemp * cos(thetaAlt);
                                            else   tempz = geometries.at(g)[4] + geometries.at(g)[5] - fabs(trackTemp * cos(thetaAlt));
                                        }
                                        else
                                        {
                                            if (!reverseDirAlt)  tempz = z0Alt + fabs(trackTemp * cos(thetaAlt));
                                            else { tempz = z0Alt - fabs(trackTemp * cos(thetaAlt)); }
                                        }

                                        float* coordinates = new float[5];
                                        coordinates[0] = cos(phiAlt) * TMath::Abs(tan(thetaAlt) * (tempz - z0Alt)) + xAlt;
                                        coordinates[1] = sin(phiAlt) * TMath::Abs(tan(thetaAlt) * (tempz - z0Alt)) + yAlt;
                                        coordinates[2] = tempz;  coordinates[3] = energyOld; coordinates[4] = 0;

                                        if (neutronCoordinates.size() < numberOfFrames + 1)
                                        {
                                            neutronCoordinates.push_back(coordinates);
                                        }
                                        tt++;
                                        if (tt > trackIterationCutoff) break;
                                    }

                                    if (tt > 0) timeTrans = timeFrame - (timeTemp - trackTemp / nSpeed);
                                    else timeTrans = timeTrans - timeTemp;
                                }
                            }

                            if (!noTrackRecording)
                            {
                                if ((g == detectorLayer) || (trackAllLayers))
                                {
                                    neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                    neutronTrackCoordinates.push_back(x);
                                    neutronTrackCoordinates.push_back(y);
                                    neutronTrackCoordinates.push_back(z0);
                                    neutronTrackCoordinates.push_back(theta);
                                    neutronTrackCoordinates.push_back(phi);
                                    neutronTrackCoordinates.push_back(energy);
                                    neutronTrackCoordinates.push_back(timeNtr);
                                    neutronTrackCoordinates.push_back(elementID);
                                    neutronTrackCoordinates.push_back(material);
                                }
                                if (showDensityTrackMapSide)
                                {
                                    neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                    neutronTrackCoordinates2.push_back(x);
                                    neutronTrackCoordinates2.push_back(y);
                                    neutronTrackCoordinates2.push_back(z0);
                                    neutronTrackCoordinates2.push_back(theta);
                                    neutronTrackCoordinates2.push_back(phi);
                                    neutronTrackCoordinates2.push_back(energy);
                                }
                            }
                        } //if (!foundSomething)
                    } //if wwRange<length
                    else
                    {
                        foundSomething = false;
                        differentMaterialHit = false;
                        domainWallHit = false;
                        currentG = g;

                        if ((scatteredThisLayer) || (continuedThisLayer))
                        {
                            z0Here = z0;
                            xHere = x;
                            yHere = y;
                        }
                        else
                        {
                            if (theta < piHalf)
                            {
                                z0Here = tempz;
                                xHere = xt;
                                yHere = yt;
                            }
                            else
                            {
                                z0Here = tempzEnd;
                                xHere = xtEnd;
                                yHere = ytEnd;
                            }
                        }

                        if (useImage)
                        {
                            if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2) || (inputPics2[currentG] == 1) || (inputPics2[currentG] == 2) || (inputPics3[currentG] == 1) || (inputPics3[currentG] == 2))
                            {
                                trackMetricFactor = trackMetricFactorModifier * fabs(1000. * matrixMetricFactor / tanTheta);
                                if (trackMetricFactor < 0.0005 * matrixMetricFactor) trackMetricFactor = 0.0005 * matrixMetricFactor;
                                if (trackMetricFactor > 1e9 * matrixMetricFactor) trackMetricFactor = 1e9 * matrixMetricFactor;

                                ztrCounter = 0;

                                if ((scatteredThisLayer) || (continuedThisLayer)) z0Here = z0;
                                else
                                {
                                    if (theta < piHalf) z0Here = tempz;
                                    else  z0Here = tempzEnd;
                                }

                                continueTracking = true;
                                differentMaterial1 = false;
                                differentMaterial2 = false;
                                differentMaterial3 = false;

                                for (double ztr = z0Here; continueTracking;)
                                {
                                    if (ztrCounter > 10000) continueTracking = false;

                                    if (ztr == z0Here)
                                    {
                                        if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2))  startMaterial = actualMaterial;
                                        if ((inputPics2[currentG] == 1) || (inputPics2[currentG] == 2))  startMaterial2 = actualMaterial2;
                                        if ((inputPics3[currentG] == 1) || (inputPics3[currentG] == 2))  startMaterial3 = actualMaterial3;
                                    }

                                    if (theta < piHalf)
                                    {
                                        ztr = ztr + trackMetricFactor;
                                        if (ztr > tempzEnd)
                                        {
                                            continueTracking = false;
                                            ztr = tempzEnd;
                                        }
                                    }
                                    else
                                    {
                                        ztr = ztr - trackMetricFactor;
                                        if (ztr < tempz)
                                        {
                                            continueTracking = false;
                                            ztr = tempz;
                                        }
                                    }

                                    if ((scatteredThisLayer) || (continuedThisLayer))
                                    {
                                        xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + x;
                                        yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + y;
                                    }
                                    else
                                    {
                                        if (theta > piHalf)
                                        {
                                            xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + xtEnd;
                                            yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + ytEnd;
                                        }
                                        else
                                        {
                                            xTrack = cosPhi * fabs(tanTheta * (ztr - z0Here)) + xt;
                                            yTrack = sinPhi * fabs(tanTheta * (ztr - z0Here)) + yt;
                                        }
                                    }

                                    if (continueTracking)
                                    {
                                        matrixX = (xTrack * 0.001 - matrixStartX) / matrixMetricFactor;
                                        matrixY = (yTrack * 0.001 - matrixStartY) / matrixMetricFactor;

                                        if (matrixX > inputMatrixPixels - 1) matrixX = inputMatrixPixels - 1; if (matrixY > inputMatrixPixels - 1) matrixY = inputMatrixPixels - 1;
                                        if (matrixX < 0) matrixX = 0; if (matrixY < 0) matrixY = 0;

                                        if ((inputPics[currentG] == 1) || (inputPics[currentG] == 2)) { selectedMaterial = (inputPicVector.at(currentG))(matrixX, matrixY); if ((selectedMaterial != actualMaterial) || (selectedMaterial != startMaterial)) differentMaterial1 = true; else differentMaterial1 = false; }
                                        if ((inputPics2[currentG] == 1) || (inputPics2[currentG] == 2)) { selectedMaterial2 = (inputPicVector2.at(currentG))(matrixX, matrixY); if ((selectedMaterial2 != actualMaterial2) || (selectedMaterial2 != startMaterial2)) differentMaterial2 = true; else differentMaterial2 = false; }
                                        if ((inputPics3[currentG] == 1) || (inputPics3[currentG] == 2)) { selectedMaterial3 = (inputPicVector3.at(currentG))(matrixX, matrixY); if ((selectedMaterial3 != actualMaterial3) || (selectedMaterial3 != startMaterial3)) differentMaterial3 = true; else differentMaterial3 = false; }
                                    }

                                    if (((differentMaterial1) || (differentMaterial2) || (differentMaterial3)) && (continueTracking))
                                    {
                                        continueTracking = false;
                                        foundSomething = true;
                                        differentMaterialHit = true;
                                        zTrack = ztr;
                                    }
                                    ztrCounter++;
                                }
                            }
                        }

                        if (reflectiveBoundaries || periodicBoundaries)
                        {
                            sideNumber = 0;

                            if (!foundSomething)
                            {
                                vectorDirFactor = (-0.5 * squareDim - xHere) / (cosPhi * tanTheta);
                                zPlane = z0Here + vectorDirFactor;
                                yPlane = yHere + (zPlane - z0Here) * sinPhi * tanTheta;
                                if ((fabs(yPlane) < 0.5 * squareDim) && (cosPhi < 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                {
                                    domainWallHit = true; foundSomething = true; sideNumber = 3;
                                    if (reflectiveBoundaries) xPlane = -squareDim * 0.5;
                                    if (periodicBoundaries) xPlane = squareDim * 0.5;
                                }
                            }
                            if (!foundSomething)
                            {
                                vectorDirFactor = (0.5 * squareDim - xHere) / (cosPhi * tanTheta);
                                zPlane = z0Here + vectorDirFactor;
                                yPlane = yHere + (zPlane - z0Here) * sinPhi * tanTheta;
                                if ((fabs(yPlane) < 0.5 * squareDim) && (cosPhi > 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                {
                                    domainWallHit = true; foundSomething = true; sideNumber = 1;
                                    if (reflectiveBoundaries) xPlane = squareDim * 0.5;
                                    if (periodicBoundaries) xPlane = -squareDim * 0.5;
                                }
                            }
                            if (!foundSomething)
                            {
                                vectorDirFactor = (0.5 * squareDim - yHere) / (sinPhi * tanTheta);
                                zPlane = z0Here + vectorDirFactor;
                                xPlane = xHere + (zPlane - z0Here) * cosPhi * tanTheta;
                                if ((fabs(xPlane) < 0.5 * squareDim) && (sinPhi > 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                {
                                    domainWallHit = true;  foundSomething = true; sideNumber = 2;
                                    if (reflectiveBoundaries) yPlane = squareDim * 0.5;
                                    if (periodicBoundaries) yPlane = -squareDim * 0.5;
                                }
                            }
                            if (!foundSomething)
                            {
                                vectorDirFactor = (-0.5 * squareDim - yHere) / (sinPhi * tanTheta);
                                zPlane = z0Here + vectorDirFactor;
                                xPlane = xHere + (zPlane - z0Here) * cosPhi * tanTheta;
                                if ((fabs(xPlane) < 0.5 * squareDim) && (sinPhi < 0) && (zPlane > tempz) && (zPlane < tempzEnd))
                                {
                                    domainWallHit = true;  foundSomething = true; sideNumber = 4;
                                    if (reflectiveBoundaries) yPlane = -squareDim * 0.5;
                                    if (periodicBoundaries) yPlane = squareDim * 0.5;
                                }
                            }
                        }

                        if (foundSomething)
                        {
                            continuedThisLayer = true;

                            if (domainWallHit)
                            {
                                z0 = zPlane;
                                x = xPlane;
                                y = yPlane;

                                if (reflectiveBoundaries)
                                {
                                    //theta = pi-theta;
                                    if ((sideNumber == 1) || (sideNumber == 3)) phi = pi - phi;
                                    else phi = 2. * pi - phi;
                                }
                                else
                                {
                                }
                            }
                            if (differentMaterialHit) //useimage
                            {
                                x = xTrack;
                                y = yTrack;
                                z0 = zTrack;
                            }

                            if (!noTrackRecording)
                            {
                                if ((scatteredThisLayer) || (continuedThisLayer))
                                {
                                    if ((g == detectorLayer) || (trackAllLayers))
                                    {
                                        neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                        neutronTrackCoordinates.push_back(x);
                                        neutronTrackCoordinates.push_back(y);
                                        neutronTrackCoordinates.push_back(z0);
                                        neutronTrackCoordinates.push_back(theta);
                                        neutronTrackCoordinates.push_back(phi);
                                        neutronTrackCoordinates.push_back(energy);
                                        neutronTrackCoordinates.push_back(timeNtr);
                                        neutronTrackCoordinates.push_back(0);
                                        neutronTrackCoordinates.push_back(material);
                                    }
                                    if (showDensityTrackMapSide)
                                    {
                                        neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                        neutronTrackCoordinates2.push_back(x);
                                        neutronTrackCoordinates2.push_back(y);
                                        neutronTrackCoordinates2.push_back(z0);
                                        neutronTrackCoordinates2.push_back(theta);
                                        neutronTrackCoordinates2.push_back(phi);
                                        neutronTrackCoordinates2.push_back(energy);
                                    }
                                }
                            }
                        }
                        else { continuedThisLayer = false; }

                        if (!foundSomething) { leaveLayer = true; }
                    }

                    if ((g == groundLayer - 1) && (!scatteredThisLayer) && (!continuedThisLayer))
                    {
                        scatteredSurfaceSpectrum->Fill(energy);

                        if (false)
                        {
                            if (reverseDir) thetas1->Fill(theta);
                            if (reverseDir) { temp2 = fabs(1. / sin(theta)); if (temp2 > 100) temp2 = 100; thetas2->Fill(theta, temp2); }
                            if (reverseDir) thetaEnergy->Fill(theta, energy);
                        }

                        if ((hasPassedSurface) && true) scatteredSurfaceSpectrumBack->Fill(energy);

                        if (hasPassedSurface)
                        {
                            if (!hasBeenReorded) scatteredSurfaceSpectrumOnce->Fill(energy);
                             hasBeenReorded = true;
                        }

                        if (drawSingleNeutronGraphs)
                        {
                            neutronPath->SetPoint(graphCounter, x * 0.001, 4); graphCounter++; neutronPath->SetTitle("Surface");
                            graphCounterMG++;
                            graphs.at(graphN)->SetPoint(graphCounterMG, x * 0.001, 4); graphs.at(graphN)->SetTitle("Surface");
                        }
                    }

                    if ((g == groundLayer - 1) && (!scatteredThisLayer) && (reverseDir) && (!continuedThisLayer))
                    {


                        xySqrt = sqrt(TMath::Power(xAtInterface - xt, 2) + TMath::Power(yAtInterface - yt, 2));

                        if (hasPassedSurface)
                        {
                            if (xySqrt > 0) scatteredSurfaceDistance->Fill(xySqrt);
                            if ((xySqrt > 45) && ((energy > energylowTHL) && (energy < energyhighTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) scatteredSurfaceDetectionDistance->Fill(xySqrt);
                            scatteredSurfaceDepth->Fill(z0);
                            scatteredSurfaceMaxDepth->Fill(z0max);
                            if (!dontFillunecessaryPlots)
                            {
                                 depthEnergyScattered2->Fill(z0max, energyAtInterface);
                            }
                        }

                        if (!dontFillunecessaryPlots)
                        {
                            energyAbsorbedCosmic->Fill(energyInitial, energy);
                        }
                    }

                    if ((noThermalRegime) && (energy < 5e-8))
                    {
                        absorbBreak = true; killedThermal++;
                    }

                    if ((!noThermalRegime) && (energy < 5e-8))
                    {
                        if ((scatteredThisLayer) && (regionalthermalcutoff))
                        {
                            if ((true) && (TMath::Power(x, 2) + TMath::Power(y, 2) > TMath::Power(regionalthermalcutoffradius, 2)))
                            {
                                absorbBreak = true; killedThermal++;
                            }
                            if (z0 - geometries.at(groundLayer)[4] > 6e2)
                            {
                                absorbBreak = true; killedThermal++;
                            }
                        }
                    }

                    if (((((g == detectorLayer) || (detectorLayerOverride)) || (trackAllLayers)) && ((leaveLayer) || (absorbBreak)) && ((drawDensityMap) && (!noTrackRecording))))
                    {
                        if (((useDetectorSensitiveMaterial) && (material == detectorSensitiveMaterial)) || (!useDetectorSensitiveMaterial))
                        {
                            if (absorbBreak)
                            {
                                if (detHitted)
                                {
                                    neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                    neutronTrackCoordinates.push_back(x);
                                    neutronTrackCoordinates.push_back(y);
                                    neutronTrackCoordinates.push_back(z0);
                                    neutronTrackCoordinates.push_back(theta);
                                    neutronTrackCoordinates.push_back(phi);
                                    neutronTrackCoordinates.push_back(energy);
                                    neutronTrackCoordinates.push_back(timeNtr);
                                    neutronTrackCoordinates.push_back(1);
                                    neutronTrackCoordinates.push_back(0);
                                }
                            }
                            else
                            {
                                //endpoints of the track
                                if (theta > piHalf)
                                {
                                    neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                    neutronTrackCoordinates.push_back(xt);
                                    neutronTrackCoordinates.push_back(yt);
                                    neutronTrackCoordinates.push_back(tempz);
                                }
                                else
                                {
                                    neutronTrackCoordinates.reserve(neutronTrackCoordinates.size() + 9);
                                    neutronTrackCoordinates.push_back(xtEnd);
                                    neutronTrackCoordinates.push_back(ytEnd);
                                    neutronTrackCoordinates.push_back(tempzEnd);
                                }
                                neutronTrackCoordinates.push_back(theta);
                                neutronTrackCoordinates.push_back(phi);
                                neutronTrackCoordinates.push_back(energy);
                                neutronTrackCoordinates.push_back(timeNtr);
                                neutronTrackCoordinates.push_back(0);
                                neutronTrackCoordinates.push_back(0);
                            }

                            try
                            {
                                //#pragma omp for
                                for (int nt = 0; nt < (neutronTrackCoordinates.size() - 9); nt += 9)
                                {
                                    thetaTr = neutronTrackCoordinates.at(3 + nt);
                                    phiTr = neutronTrackCoordinates.at(4 + nt);
                                    energyTr = neutronTrackCoordinates.at(5 + nt);

                                    cosPhiTr = cos(phiTr);
                                    sinPhiTr = sin(phiTr);
                                    tanThetaTr = tan(thetaTr);

                                    trackMetricFactor = fabs(mapMetricFactor * 1000. / tanThetaTr)*0.4;

                                    continueTracking = true;
                                    ztrCounter = 0;

                                    for (double ztr = neutronTrackCoordinates.at(2 + nt); continueTracking;)
                                    {
                                        if (ztrCounter > 10000) continueTracking = false;

                                        if ((ztrCounter == 0) && (nt % 18 == 0))
                                        {
                                        }
                                        else
                                        {
                                            if (thetaTr < piHalf)
                                            {
                                                ztr += trackMetricFactor;
                                                if (ztr > neutronTrackCoordinates.at(2 + 9 + nt)) continueTracking = false;
                                            }
                                            else
                                            {
                                                ztr -= trackMetricFactor;
                                                if (ztr < neutronTrackCoordinates.at(2 + 9 + nt)) continueTracking = false;
                                            }

                                            if (!continueTracking) break;
                                        }

                                        ztrCounter++;

                                        xTrack = cosPhiTr * fabs(tanThetaTr * (ztr - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(0 + nt);
                                        yTrack = sinPhiTr * fabs(tanThetaTr * (ztr - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(1 + nt);

                                        if ((xTrack < squareDim * 0.5) && (xTrack > -squareDim * 0.5) && (yTrack < squareDim * 0.5) && (yTrack > -squareDim * 0.5))
                                        {
                                            if ((energyTr < 2e-7) && (energyTr > 1e-9))  densityThermalTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);
                                            if ((energyTr < 0.001) && (energyTr > 0.000001))  densityTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);
                                            if ((energyTr < 0.5) && (energyTr > 0.001)) densityIntermediateTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);
                                            if ((energyTr < 10) && (energyTr > 0.5)) densityFastTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);
                                            if (((energyTr < energyhighTHL) && (energyTr > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) densityAlbedoTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);
                                            if ((energyTr < 100000) && (energyTr > 20)) densityHighEnergyTrackMap->Fill(xTrack * 0.001, yTrack * 0.001);

                                            // dothewaterthing
                                            //if ((energyTr<1.7e-6)&&(energyTr>1.2e-6)) thetas2->Fill( 1./1000.*sqrt(TMath::Power(xAtInterface-xTrack,2)+TMath::Power(yAtInterface-yTrack,2)+TMath::Power(zPosSource-ztr,2))  );
                                            //if ((energyTr<0.7e-7)&&(energyTr>1.2e-9)) thetas1->Fill( 1./1000.*sqrt(TMath::Power(xAtInterface-xTrack,2)+TMath::Power(yAtInterface-yTrack,2)+TMath::Power(zPosSource-ztr,2)) );

                                            densityEnergyTrackMap->SetBinContent(densityEnergyTrackMap->GetXaxis()->FindBin(xTrack * 0.001), densityEnergyTrackMap->GetYaxis()->FindBin(yTrack * 0.001), log(energyTr) + 20.);

                                            if (highResCalc)
                                            {
                                                const short hrmax = 2; // oversampling ratio
                                                for (int hr = 0; hr < hrmax; hr++)
                                                {
                                                    //that is because mapMetricFactor is defined for a 500 px map and so is trackMetricFactor and therefore an intermediate value for 1000 px has to be inserted
                                                    if (hr != 0)
                                                    {
                                                        if (thetaTr < piHalf)
                                                        {
                                                            xTrack = cosPhiTr * fabs(tanThetaTr * ((ztr + 1./(hrmax*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(0 + nt);
                                                            yTrack = sinPhiTr * fabs(tanThetaTr * ((ztr + 1./(hrmax*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(1 + nt);
                                                        }
                                                        else
                                                        {
                                                            xTrack = cosPhiTr * fabs(tanThetaTr * ((ztr - 1./(hrmax*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(0 + nt);
                                                            yTrack = sinPhiTr * fabs(tanThetaTr * ((ztr - 1./(hrmax*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(1 + nt);
                                                        }
                                                    }

                                                    if ((energyTr < 2e-7) && (energyTr > 1e-9))  densityThermalTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 0.001) && (energyTr > 0.000001))  densityTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 0.5) && (energyTr > 0.001)) densityIntermediateTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 10) && (energyTr > 0.5)) densityFastTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if (((energyTr < energyhighTHL) && (energyTr > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) densityAlbedoTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 100000) && (energyTr > 20)) densityHighEnergyTrackMapHighRes->Fill(xTrack * 0.001, yTrack * 0.001);

                                                    densityEnergyTrackMapHighRes->SetBinContent(densityEnergyTrackMapHighRes->GetXaxis()->FindBin(xTrack * 0.001), densityEnergyTrackMapHighRes->GetYaxis()->FindBin(yTrack * 0.001), log(energyTr) + 20.);
                                                }
                                            }

                                            if ((!noGUIMode) && (highhighResCalc))
                                            {
                                                const short hrmax2 = 6;
                                                for (int hr = 0; hr < hrmax2; hr++)
                                                {
                                                    //that is because mapMetricFactor is defined for a 500 px map and so is trackMetricFactor and therefore an intermediate value for 1000 px has to be inserted
                                                    if (hr != 0)
                                                    {
                                                        if (thetaTr < piHalf)
                                                        {
                                                            xTrack = cosPhiTr * fabs(tanThetaTr * ((ztr + 1./(hrmax2*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(0 + nt);
                                                            yTrack = sinPhiTr * fabs(tanThetaTr * ((ztr + 1./(hrmax2*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(1 + nt);
                                                        }
                                                        else
                                                        {
                                                            xTrack = cosPhiTr * fabs(tanThetaTr * ((ztr - 1./(hrmax2*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(0 + nt);
                                                            yTrack = sinPhiTr * fabs(tanThetaTr * ((ztr - 1./(hrmax2*1.) * (hr*1.) * trackMetricFactor) - neutronTrackCoordinates.at(2 + nt))) + neutronTrackCoordinates.at(1 + nt);
                                                        }
                                                    }

                                                    if ((energyTr < 2e-7) && (energyTr > 1e-9))  densityThermalTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 0.001) && (energyTr > 0.000001))  densityTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 0.5) && (energyTr > 0.001)) densityIntermediateTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 10) && (energyTr > 0.5)) densityFastTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if (((energyTr < energyhighTHL) && (energyTr > energylowTHL) && (!useRealisticModelLayer)) || (layerRealisticallyHitted)) densityAlbedoTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);
                                                    if ((energyTr < 100000) && (energyTr > 20)) densityHighEnergyTrackMapHighRes15x->Fill(xTrack * 0.001, yTrack * 0.001);

                                                    densityEnergyTrackMapHighRes15x->SetBinContent(densityEnergyTrackMapHighRes15x->GetXaxis()->FindBin(xTrack * 0.001), densityEnergyTrackMapHighRes15x->GetYaxis()->FindBin(yTrack * 0.001), log(energyTr) + 20.);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            catch (std::out_of_range& e)
                            {
                                //cout<<"Out of Range Exception in Tracking";
                                //cout<<n<<endl<<g<<endl<<energy<<endl<<phi<<endl<<theta<<endl<<x<<endl<<y<<endl<<z0<<endl;
                            }

                            if ((detTrackFileOutput) || (allTrackFileOutput))
                            {
                                neutronTrackCoordinatesFullSet.reserve(neutronTrackCoordinatesFullSet.size() + neutronTrackCoordinates.size());
                                neutronTrackCoordinatesFullSet.insert(neutronTrackCoordinatesFullSet.end(), neutronTrackCoordinates.begin(), neutronTrackCoordinates.end());
                            }

                            neutronTrackCoordinates.clear();
                        }
                    }

                    if (((showDensityTrackMapSide) && (!noTrackRecording) && (energy < densitySideTrackingEnergyCutoff) && ((leaveLayer) || (absorbBreak))) && (true))
                    {
                        if (absorbBreak);
                        else
                        {
                            //endpoints of the track
                            if (theta > piHalf)
                            {
                                neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                neutronTrackCoordinates2.push_back(xt);
                                neutronTrackCoordinates2.push_back(yt);
                                neutronTrackCoordinates2.push_back(tempz);
                                neutronTrackCoordinates2.push_back(theta);
                                neutronTrackCoordinates2.push_back(phi);
                                neutronTrackCoordinates2.push_back(energy);
                            }
                            else
                            {
                                neutronTrackCoordinates2.reserve(neutronTrackCoordinates2.size() + 6);
                                neutronTrackCoordinates2.push_back(xtEnd);
                                neutronTrackCoordinates2.push_back(ytEnd);
                                neutronTrackCoordinates2.push_back(tempzEnd);
                                neutronTrackCoordinates2.push_back(theta);
                                neutronTrackCoordinates2.push_back(phi);
                                neutronTrackCoordinates2.push_back(energy);
                            }
                        }

                        trackMetricFactorTemp = 1.0005 * fabs((domainZUpperEdge - domainZLowEdge) * 0.5);

                        for (int nt = 0; nt < (neutronTrackCoordinates2.size() - 6); nt += 6)
                        {
                            thetaTr = neutronTrackCoordinates2[3 + nt];
                            phiTr = neutronTrackCoordinates2[4 + nt];
                            energyTr = neutronTrackCoordinates2[5 + nt];

                            cosPhiTr = cos(phiTr);
                            sinPhiTr = sin(phiTr);
                            tanThetaTr = tan(thetaTr);

                            trackMetricFactor = fabs(trackMetricFactorTemp * cos(TMath::ATan(tanThetaTr * aspectRatioSide)));

                            if (trackMetricFactor < 0.001 * trackMetricFactorTemp) trackMetricFactor = 0.001 * trackMetricFactorTemp;
                            if (trackMetricFactor > 1e5 * trackMetricFactorTemp) trackMetricFactor = 1e5 * trackMetricFactorTemp;

                            continueTracking = true;
                            ztrCounter = 0;

                            for (double ztr = neutronTrackCoordinates2.at(2 + nt); continueTracking;)
                            {
                                if (ztrCounter > 40000) continueTracking = false;

                                if (ztrCounter != 0)
                                {
                                    if (thetaTr < piHalf)
                                    {
                                        ztr += trackMetricFactor;
                                        if (ztr > neutronTrackCoordinates2.at(2 + 6 + nt)) continueTracking = false;
                                    }
                                    else
                                    {
                                        ztr -= trackMetricFactor;
                                        if (ztr < neutronTrackCoordinates2.at(2 + 6 + nt)) continueTracking = false;
                                    }

                                    if (!continueTracking) break;
                                }
                                ztrCounter++;

                                xTrack = cosPhiTr * fabs(tanThetaTr * (ztr - neutronTrackCoordinates2.at(2 + nt))) + neutronTrackCoordinates2.at(0 + nt);
                                yTrack = sinPhiTr * fabs(tanThetaTr * (ztr - neutronTrackCoordinates2.at(2 + nt))) + neutronTrackCoordinates2.at(1 + nt);

                                if ((xTrack < squareDim * 0.5) && (xTrack > -squareDim * 0.5) && (yTrack < densityTrackMapSideCutOutValue) && (yTrack > -densityTrackMapSideCutOutValue))
                                {
                                    if ((-ztr > domainZLowEdge * 1000.) && (-ztr < domainZUpperEdge * 1000.) && (true))
                                    {
                                        densityTrackMapSide->Fill(xTrack * 0.001, -0.001 * ztr);
                                        if (hasBeenInSoil) densityTrackMapSideAlbedo->Fill(xTrack * 0.001, -0.001 * ztr);
                                        if ((energyTr < energyhighTHL) && (energyTr > energylowTHL)) densityTrackMapSideDetector->Fill(xTrack * 0.001, -0.001 * ztr);
                                        if ((!noThermalRegime) && (energyTr < 2e-7)) densityTrackMapSideThermal->Fill(xTrack * 0.001, -0.001 * ztr);
                                    }
                                }
                            }
                        }
                        neutronTrackCoordinates2.clear();
                    }

                    previousTheta = theta; previousPhi = phi;
                    previousX = x; previousY = y, previousZ0 = z0;
                    previousCosPhi = cosPhi; previousSinPhi = sinPhi; previousTanTheta = tanTheta, previousCosTheta = cosTheta;

                    if ((material == 18) || (material == 19) || (material == 20) || (material == 9) || (material == 21)) { previousSoilX = x; previousSoilY = y; previousSoilZ0 = z0; }
                    if (thermalized) { if (thermalizedX == thermalizedY) { thermalizedX = x; thermalizedY = y; thermalizedZ0 = z0; } }

                    if (theta > piHalf) reverseDir = true;   else reverseDir = false;

                    //scatteredLayer->Fill(g+1);

                    glayerCounter++;

                    if (leaveLayer)
                    {
                        startMaterial = actualMaterial;
                        scatteredThisLayer = false;
                        hasbeenEvaporized = false;
                        leaveLayer = false;
                        //recordedThisLayer = false;
                        newLayer = true;
                        glayerCounter = 0;
                        if (reverseDir) g--;
                        else g++;
                        stillFirstLayer = false;
                    }

                    gAlt = g;

                    if (g < 0) escapes++;

                    if ((absorbBreak) || (g < 0) || (g == geometries.size()))
                    {
                        if (allTrackFileOutput)
                        {
                            if ((true) && (neutronTrackCoordinatesFullSet.size() > 0))
                            {
                                if (doBatchRun2D) allTrackOutputFile.open(outputFolder + "NeutronTrackData_" + castIntToString(paramInt) + ".dat", ios::out | ios::app);
                                else allTrackOutputFile.open(outputFolder + "NeutronTrackData.dat", ios::out | ios::app);

                                for (int nt = 0; nt <= (neutronTrackCoordinatesFullSet.size() - 9); nt += 9)
                                {
                                    if ((nt > 0) && ((neutronTrackCoordinatesFullSet[nt + 5] != neutronTrackCoordinatesFullSet[nt - 4]) || (neutronTrackCoordinatesFullSet[nt + 7] == 0) || (neutronTrackCoordinatesFullSet[nt + 7] == 1)))
                                    {
                                        allTrackOutputFile << n << "\t" << neutronTrackCoordinatesFullSet.at(nt + 0) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << endl;
                                    }
                                    else
                                    {
                                        if (nt == 0) allTrackOutputFile << n << "\t" << neutronTrackCoordinatesFullSet.at(nt + 0) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 1) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 2) * 0.001 << "\t" << neutronTrackCoordinatesFullSet.at(nt + 3) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 4) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 5) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 6) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 7) << "\t" << neutronTrackCoordinatesFullSet.at(nt + 8) << endl;
                                    }
                                    //previousTrackEnergy = neutronTrackCoordinatesFullSet.at(nt+5);
                                }
                                allTrackOutputFile.close();
                            }
                        }
                    }

                    if (detHitted)
                    {
                        if (otherAllTrackFileOutput)
                        {
                            if ((neutronAbsorbedbyDetector || detHitted) && (neutronTrackCoordinatesFullSet.size() > 0))
                            {
                                double xtJ, ytJ, ztJ, etJ, ztJ2;
                                double metricFac = 0.4; //0.001*squareDim;

                                for (int nt = 0; nt <= (neutronTrackCoordinatesFullSet.size() - 18); nt += 9)
                                {
                                    xtJ = (neutronTrackCoordinatesFullSet.at(nt + 0) * 0.5 + squareDim * 0.15) / squareDim * 700.;
                                    //xtJ2 = neutronTrackCoordinatesFullSet.at(nt+7+0)/squareDim+squareDim/2.;
                                    //ztJ2 = neutronTrackCoordinatesFullSet.at(nt+7+2)/squareDim+squareDim/2.;
                                    ytJ = (neutronTrackCoordinatesFullSet.at(nt + 1) * 0.5 + squareDim * 0.15) / squareDim * 700.;
                                    //ytJ2 = neutronTrackCoordinatesFullSet.at(nt+7+1)/squareDim+squareDim/2.;
                                    //ztJ = (neutronTrackCoordinatesFullSet.at(nt+2)*1.)/squareDim;
                                    ztJ = (neutronTrackCoordinatesFullSet.at(nt + 2) * 0.5 + squareDim * 0.5) / squareDim * 700.;
                                    ztJ2 = (neutronTrackCoordinatesFullSet.at(nt + 9 + 2) * 0.5 + squareDim * 0.5) / squareDim * 700.;

                                    etJ = TMath::Log10(neutronTrackCoordinatesFullSet.at(nt + 5) * 1e6) * 1000.;

                                    thetaTr = neutronTrackCoordinatesFullSet[3 + nt];
                                    phiTr = neutronTrackCoordinatesFullSet[4 + nt];

                                    cosPhiTr = cos(phiTr);
                                    sinPhiTr = sin(phiTr);
                                    tanThetaTr = tan(thetaTr);

                                    trackMetricFactor = fabs(metricFac * cos(thetaTr));

                                    continueTracking = true;
                                    ztrCounter = 0;

                                    for (double ztr = ztJ; continueTracking;)
                                    {
                                        if (ztrCounter > 4000) continueTracking = false;

                                        if (ztrCounter != 0)
                                        {
                                            if (thetaTr < piHalf)
                                            {
                                                ztr += trackMetricFactor;
                                                if (ztr > ztJ2) continueTracking = false;
                                            }
                                            else
                                            {
                                                ztr -= trackMetricFactor;
                                                if (ztr < ztJ2) continueTracking = false;
                                            }

                                            if (!continueTracking) break;
                                        }
                                        ztrCounter++;

                                        xTrack = cosPhiTr * fabs(tanThetaTr * (ztr - ztJ)) + xtJ;
                                        yTrack = sinPhiTr * fabs(tanThetaTr * (ztr - ztJ)) + ytJ;

                                        xCoord.push_back(xTrack);
                                        yCoord.push_back(yTrack);
                                        zCoord.push_back(ztr);
                                        eCoord.push_back(etJ);
                                    }
                                }

                                neutronTrackNeutrons++;

                                if (neutronTrackNeutrons % neutronTrackNeutronsToExport == 0)
                                {
                                    int ntJ = neutronTrackNeutrons / neutronTrackNeutronsToExport;

                                    string folderString = (string)outputFolder+"data/"+(string)castIntToString(ntJ)+"/";
                                    std::replace(folderString.begin(), folderString.end(), '/', '\\');
                                    string foldercommandstring = "mkdir " + folderString;
                                    const char* foldercommandstring2 = foldercommandstring.c_str();
                                    system(foldercommandstring2);

                                    allTrackOutputFile.open(folderString+(string)castIntToString(ntJ)+"-NeutronTrackData"+(string)castIntToString(ntJ)+".json", ios::out | ios::app);

                                    allTrackOutputFile << "{\n\"runNo\":\"0\"," << endl << "\"eventNo\":\"" << ntJ << "\"," << endl;

                                    allTrackOutputFile << "\"x\": [";
                                    for (int cc = 0; cc < xCoord.size(); cc++)
                                    {
                                        allTrackOutputFile << xCoord.at(cc); if (cc != xCoord.size() - 1) allTrackOutputFile << ",";
                                    }
                                    allTrackOutputFile << "]," << endl;

                                    allTrackOutputFile << "\"y\": [";
                                    for (int cc = 0; cc < yCoord.size(); cc++)
                                    {
                                        allTrackOutputFile << yCoord.at(cc); if (cc != yCoord.size() - 1) allTrackOutputFile << ",";
                                    }
                                    allTrackOutputFile << "]," << endl;

                                    allTrackOutputFile << "\"z\": [";
                                    for (int cc = 0; cc < zCoord.size(); cc++)
                                    {
                                        allTrackOutputFile << zCoord.at(cc); if (cc != zCoord.size() - 1) allTrackOutputFile << ",";
                                    }
                                    allTrackOutputFile << "]," << endl;

                                    allTrackOutputFile << "\"q\": [";
                                    for (int cc = 0; cc < eCoord.size(); cc++)
                                    {
                                        allTrackOutputFile << eCoord.at(cc); if (cc != eCoord.size() - 1) allTrackOutputFile << ",";
                                    }
                                    allTrackOutputFile << "]," << endl;

                                    allTrackOutputFile << "\"geom\" : \"protodune\"," << endl;
                                    allTrackOutputFile << "\"type\" : \"wire-cell\"" << endl;

                                    allTrackOutputFile << "}" << endl;

                                    allTrackOutputFile.close();

                                    xCoord.clear(); yCoord.clear(); zCoord.clear(); eCoord.clear();
                                }
                            }
                            neutronTrackCoordinatesFullSet.clear();
                        }
                    }

                    if (absorbBreak)  break;
                }
                //end of For loop of stack

                if (drawSingleNeutronGraphs)
                {
                    TCanvas* cGraph = new TCanvas(string("cGraph "+castLongToString(n)).c_str(),"Neutron Path",1280,1280);
                    CanvasFashion(cGraph); TGraphFashion(neutronPath, "lateral Position x [m]", "Depth [mm]", true); //neutronPath->GetYaxis()->SetRangeUser(-100,5);
                    //neutronPath->GetXaxis()->SetRangeUser(-50,50);
                    //neutronPath->GetYaxis()->SetRangeUser(-1500,50000);
                    //neutronPath->SetMarkerSize(0);
                    neutronPath->Draw("APL");
                    cGraph->SaveAs((string)outputFolder+neutronPicSubPath+"/NeutronPath_"+castLongToString(n)+".png");

                    TCanvas* cMGGraph = new TCanvas(string("cMGGraph "+castLongToString(n)).c_str(),"Neutron Path",1280,1280);
                    CanvasFashion(cMGGraph);
                    for (Int_t k = 0; k < graphs.size(); k++)
                    {
                        multigraph->Add(graphs.at(k));
                    }
                    multigraph->Draw("APL");
                    TMultiGraphFashion(multigraph, "lateral Position x [m]", "Depth [mm]", true);
                    multigraph->Draw("APL");
                    //graphs.back()->Draw("PL");
                    cMGGraph->SaveAs((string)outputFolder+neutronPicSubPath+"/NeutronPathMG_"+castLongToString(n)+".png");
                }

                if (drawSingleNeutronPropagation) { allNeutronCoordinates.push_back(neutronCoordinates); }

                //end of neutron loop
            }

            //end of batch loop
            if ((doBatchRun) || (doDetectorBatchRun) || (doDetectorBatchRun2) || (doDetectorAngleBatchRun) || (doBatchRunDensity) || (doBatchRun2D))
            {
                //TString batchStr = castDoubleToString(xr)+"\t"+castDoubleToString(yr)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27))+absElement->GetBinContent(absElement->FindBin(10)));
                //TString batchStr = castDoubleToString(rHe3)+"\t"+castIntToString(neutrons)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27))+absElement->GetBinContent(absElement->FindBin(10)));
                //TString batchStr = castDoubleToString(paramEnergy)+"\t"+castIntToString(neutrons)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27))+absElement->GetBinContent(absElement->FindBin(10)));

                //TString batchStr = castDoubleToString(paramEnergy)+"\t"+castIntToString(neutrons)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27)))+"\t"+castDoubleToString(absElement->GetBinContent(absElement->FindBin(10)))+"\t"+castDoubleToString(absElementL1->GetBinContent(absElement->FindBin(10)))+"\t"+castDoubleToString(absElementL2->GetBinContent(absElement->FindBin(10)))+"\t"+castDoubleToString(absElementL3->GetBinContent(absElement->FindBin(10)));
                //TString batchStr = castDoubleToString(paramEnergy)+"\t"+castIntToString(neutrons)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27)))+"\t"+
                //        castDoubleToString(absElement->GetBinContent(absElement->FindBin(10)))+"\t"+
                //        castDoubleToString(absElementL1->GetBinContent(absElementL1->FindBin(4)))+"\t"+castDoubleToString(absElementL1->GetBinContent(absElementL1->FindBin(6)))+"\t"+castDoubleToString(absElementL1->GetBinContent(absElementL1->FindBin(8)));
                TString batchStr = castDoubleToString(batchVar) + "\t" + castLongToString(neutrons) + "\t" + castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(28))) + "\t" + castDoubleToString(absElement->GetBinContent(absElement->FindBin(10))) + "\t" + castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27))) + "\t" + castDoubleToString(absElement->GetBinContent(absElement->FindBin(3)));

                //TString batchStr2 = castDoubleToString(batchVar)+"\t"+castLongToString(neutrons)+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(27)))+"\t"+castDoubleToString(absMaterial->GetBinContent(absMaterial->FindBin(32)));
                //for(int bb = 0; bb < absElement->GetNbinsX(); bb++) {absElement->SetBinContent(bb, 0); }
                //for(int bb = 0; bb < absMaterial->GetNbinsX(); bb++) {absMaterial->SetBinContent(bb, 0); }
                //storeToFile(outputFolder, "batchRunEff2.txt", batchStr2, true);

                storeToFile((string)outputFolder, "batchRunEff.txt", batchStr, true);
                if ((doBatchRunDensity) || (doBatchRun2D)) uiM->exportToSave();
            }

        } //End of param loop

        simulationRunning = false;

        time(&end);

        Double_t dif = difftime(end, start);
        cout << "\r" << castFloatToString(100., 6) << " % completed ";
        int secondsEnd = dif;
        int minutesEnd = dif / 60;
        int hoursEnd = minutesEnd / 60;
        cout << endl << "Runtime: " << castDoubleToString(dif) << " s (" << int(hoursEnd) <<" h, "<< int(minutesEnd%60) <<" min, " << int(secondsEnd%60) <<" s)" << endl << endl;

        if (doTheCoastalTransect)
        {
            TString detOutput = castIntToString(detd) + "\t" + castIntToString(detc) + "\t" + castIntToString(detb) + "\t" + castIntToString(deta) + "\t" + castIntToString(det0)
                + "\t" + castIntToString(det1) + "\t" + castIntToString(det2) + "\t" + castIntToString(det3) + "\t" + castIntToString(det4) + "\t" + castIntToString(det5)
                + "\t" + castIntToString(det6) + "\t" + castIntToString(det7) + "\t" + castIntToString(det8) + "\t" + castIntToString(det9) + "\t" + castIntToString(det10) + "\t" + castIntToString(det11)
                + "\t" + castIntToString(det12) + "\t" + castIntToString(det13) + "\t" + castIntToString(det14) + "\t" + castIntToString(det15) + "\t" + castIntToString(det16);

            storeToFile((string)outputFolder, "detHits_"+castFloatToString(soilWaterFrac,4)+".txt", detOutput, true);

            cout << det1 << endl;		cout << det2 << endl;		cout << det3 << endl;		cout << det4 << endl;		cout << det5 << endl;		cout << det6 << endl;
            cout << det7 << endl;		cout << det8 << endl;		cout << det9 << endl;		cout << det10 << endl;		cout << det11 << endl;		cout << det12 << endl;
        }

        uiM->redrawNeutronMap(difftime(end, start) - pauseTime);

        int maxRecNeutronCoord = 0;
        int pointcter = 0;
        int inpScale;

        if (drawSingleNeutronPropagation)
        {
            for (int n = 0; n < allNeutronCoordinates.size(); n++)
            {
                if (allNeutronCoordinates.at(n).size() > maxRecNeutronCoord) maxRecNeutronCoord = allNeutronCoordinates.at(n).size();
            }

            TF1* fWater = new TF1("fWater", "+0*x", -squareDim / 1.85 / 1000., 0);  fWater->SetFillColorAlpha(kBlue, 0.5); fWater->SetLineColorAlpha(kBlue, 0.8); //fWater->SetLineWidth(150);
            TF1* fSoil = new TF1("fSoil", "0+0*x", 0, squareDim / 1.85 / 1000.);  fSoil->SetFillColorAlpha(kOrange, 0.5); fSoil->SetLineColorAlpha(kOrange, 0.8);

            for (int k = 0; ((k < maxRecNeutronCoord) && (k < numberOfFrames)); k++)
            {
                TGraph* neutronDots = new TGraph();
                pointcter = 0;

                TMultiGraph* multigraphNeutrons = new TMultiGraph();

                //put ( 0*n/10<k) for inpScale = 0;
                for (int n = 0; (n < allNeutronCoordinates.size()) && (0 * n / 10 < k); n++)
                {
                    //inpScale = n/10;
                    inpScale = 0;
                    if (allNeutronCoordinates.at(n).size() > k - inpScale)
                    {
                        if ((fabs(allNeutronCoordinates.at(n).at(k - inpScale)[0]) < squareDim / 1.9) && (fabs(allNeutronCoordinates.at(n).at(k - inpScale)[1]) < squareDim / 1.9))
                        {

                            TGraph* neutronMGDots = new TGraph();
                            //scale = (log10(allNeutronCoordinates.at(n).at(k-inpScale)[3])+7.7)/12.;	if (scale<0.001) scale = 0.001;
                            scale = (log10(allNeutronCoordinates.at(n).at(k - inpScale)[3]) + 5.9) / 7.; if (scale < 0.001) scale = 0.001; //for evapEnergies
                            rgbValues = getRGBfromHCL(getLinearC(scale, 0, 240), getScaledColorValue(scale, 0.4, 0.91, 0.9), getScaledColorValue(scale, 2.2, 0.8, 0.9));
                            neutronMGDots->SetMarkerColor(TColor::GetColor(rgbValues.at(0), rgbValues.at(1), rgbValues.at(2)));
                            neutronMGDots->SetPoint(0, (allNeutronCoordinates.at(n).at(k - inpScale)[0]) / 1000., (-1.) * (allNeutronCoordinates.at(n).at(k - inpScale)[2]) / 1000.);

                            neutronMGDots->SetMarkerSize(.9); neutronMGDots->SetMarkerStyle(20);
                            multigraphNeutrons->Add(neutronMGDots);
                            //pointcter++;
                        }
                    }
                }

                //for side view
                neutronDots->SetPoint(pointcter, -squareDim / 1.85 / 1000., 0); neutronDots->SetPoint(pointcter + 1, squareDim / 1.85 / 1000., 0);
                //for square view

                TCanvas* cGraphSingleN = new TCanvas(string("cGraphSingleN "+castIntToString(k)).c_str(),"Neutron Dots",1550,700); //700
                CanvasFashion(cGraphSingleN);
                TGraphFashion(neutronDots, "x [m]", "Depth [m]", true);
                neutronDots->SetMarkerSize(0.7); neutronDots->SetMarkerStyle(20); neutronDots->SetMarkerColor(8);
                //neutronDots->GetXaxis()->SetRangeUser(-squareDim/1.5,squareDim/1.5);
                neutronDots->GetYaxis()->SetRangeUser(-9, 145);
                //neutronDots->GetYaxis()->SetRangeUser(-9,45);
                //neutronDots->Draw("AP");

                TGraph* block = new TGraph(); block->SetFillColor(42);// block->SetFillColor(38);
                block->SetPoint(0, -squareDim / 2. / 1000., 0);  block->SetPoint(1, -squareDim / 7.4 / 1000., 0);
                block->SetPoint(2, -squareDim / 4. / 1000., 0); block->SetPoint(3, -squareDim / 8 * 3 / 1000., 0);

                TGraph* block2 = new TGraph(); block2->SetFillColor(42);
                block2->SetPoint(0, 0., 0); block2->SetPoint(1, squareDim / 7 / 1000., 0);

                neutronDots->SetLineWidth(0); neutronDots->SetMarkerSize(0.001);
                multigraphNeutrons->Add(neutronDots);
                multigraphNeutrons->Draw("APL");
                block->Draw("B1SAME");
                block2->Draw("B1SAME");
                TMultiGraphFashion(multigraphNeutrons, "x [m]", "Depth [m]", true);
                //multigraphNeutrons->GetYaxis()->SetRangeUser(-9,45);
                multigraphNeutrons->GetYaxis()->SetRangeUser(-9, 280);
                multigraphNeutrons->GetXaxis()->SetRangeUser(-600, 120);
                multigraphNeutrons->GetYaxis()->SetLabelOffset(999); multigraphNeutrons->GetYaxis()->SetTitle("");
                multigraphNeutrons->GetYaxis()->SetLabelSize(0);
                multigraphNeutrons->Draw("PL");

                //fSoil->Draw("SAME");
                //fWater->Draw("SAME");

                string foutName;
                if (k<1000) foutName = "NeutronPathSingle_0"+(string)castIntToString(k);
                if (k<100) foutName = "NeutronPathSingle_00"+(string)castIntToString(k);
                if (k<10) foutName = "NeutronPathSingle_000"+(string)castIntToString(k);

                cGraphSingleN->SaveAs((string)outputFolder+neutronPicSubPath+"/"+foutName+".png");

                delete neutronDots;
                if (multigraphNeutrons != 0x0) delete multigraphNeutrons;
            }
        }

        delay(10);
        uiM->exportToSave();
        delay(300);

        if (detOutputFile) detOutputFile.close();
        if (detLayerOutputFile) detLayerOutputFile.close();
        delay(200);

        for (int bn = 0; bn < scatteredSurfaceSpectrumBack->GetNbinsX(); bn++) { scatteredSurfaceSpectrumBack->SetBinContent(bn, 0); }
        for (int bn = 0; bn < scatteredSurfaceSpectrum->GetNbinsX(); bn++) { scatteredSurfaceSpectrum->SetBinContent(bn, 0); }
        for (int bn = 0; bn < cosmicSpectrum->GetNbinsX(); bn++) { cosmicSpectrum->SetBinContent(bn, 0); }

        for (int bn = 0; bn < scatDepth->GetNbinsX(); bn++) { scatDepth->SetBinContent(bn, 0); }
        for (int bn = 0; bn < scatteredSurfaceDepth->GetNbinsX(); bn++) { scatteredSurfaceDepth->SetBinContent(bn, 0); }
        for (int bn = 0; bn < scatteredSurfaceMaxDepth->GetNbinsX(); bn++) { scatteredSurfaceMaxDepth->SetBinContent(bn, 0); }

        if (!noGUIMode)
        {
            uiM->setStatus(1, "Clearing Data");
            cout << "Internal: ";
        }
        histoClearUp(&allTHs2);
        if (!noGUIMode) cout << "External: ";
        histoClearUp(&allTHs);
        if (!noGUIMode) cout << "done" << endl;
        if (!noGUIMode) uiM->setStatus(1, "");
    }

    catch (std::out_of_range& e)
    {
        cout << "Out of Range Exception Termination";
    }

    return true;
}

string tmpStringSimulate = "";

/**
 * Function to be called on clicking the simulate button.
 * Prepares the simulation settings and checks the geometry configuration for consistency
 *
 */
void MainWindow::on_pushButton_Simulate_clicked()
{
    Int_t tp1 = 1, tp2 = 2, mat;

    tmpStringSimulate.reserve(200);

    tmpStringSimulate = ui->lineEdit_InputSpectrumFolder->text().toStdString();
    //    std::replace(tmpStringSimulate.begin(), tmpStringSimulate.end(), '\\', '/');
    inputSpectrumFile = tmpStringSimulate; //ui->lineEdit_InputSpectrumFolder->text().toStdString();

    tmpStringSimulate = ui->lineEdit_OutputFolder->text().toStdString();
    //    std::replace(tmpStringSimulate.begin(), tmpStringSimulate.end(), '\\', '/');
    outputFolder = tmpStringSimulate; //ui->lineEdit_OutputFolder->text().toStdString();
    tmpStringSimulate = ui->lineEdit_CrosssectionFolder->text().toStdString();
    //    std::replace(tmpStringSimulate.begin(), tmpStringSimulate.end(), '\\', '/');
    endfFolder = tmpStringSimulate; // ui->lineEdit_CrosssectionFolder->text().toStdString();

    setupGeometry();

    // checks if geometry is set correctly
    if (model->rowCount() < 2)
    {
        setStatus(1, "Error: Layer Config");
        setStatus(2, "Not enough Layers");
        return;
    }

    if (TMath::Abs(((geometries.at(startingLayer)[4]) / 1000.) - ((geometries.at(groundLayer)[4]) / 1000.)) > 799)
    {
        setStatus(1, "Error: Layer Config");
        setStatus(2, "Source Layer Elevation too high");
        return;
    }

    if (((geometries.at(detectorLayer)[5]) / 1000.) > 40.)
    {
        setStatus(1, "Error: Detector Layer Config");
        setStatus(2, "Detector Layer too thick");
        return;
    }

    if (((geometries.at(startingLayer)[5]) / 1000.) > 250.)
    {
        setStatus(1, "Error: Source Layer Config");
        setStatus(2, "Source Layer too thick");
        return;
    }

    if (startingLayer < 0)
    {
        setStatus(1, "Error: Source Layer Config");
        setStatus(2, "Layer declaration < 0");
        return;
    }


    if (detectorLayer < 0)
    {
        setStatus(1, "Error: Detector Layer Config");
        setStatus(2, "Layer declaration < 0");
        return;
    }

    if (groundLayer < 0)
    {
        setStatus(1, "Error: Ground Layer Config");
        setStatus(2, "Layer declaration < 0");
        return;
    }

    for (int i = 0; i < geometries.size() - 1; i++)
    {
        if ((geometries.at(i)[4] + geometries.at(i)[5] - 1) > geometries.at(i + 1)[4])
        {
            if ((i + 1) == detectorLayer) {}
            else
            {
                tp1 = i + 1;
                tp2 = i + 2;
                string error2 = "Layers "+(string)castIntToString(tp1)+"+"+(string)castIntToString(tp2)+" colliding";
                setStatus(1, "Error: Layer Config");
                setStatus(2, error2);
                return;
            }
        }
    }

    for (int i = 0; i < geometries.size() - 1; i++)
    {
        if (geometries.at(i)[5] < 0)
        {
            tp1 = i + 1;
            string error2 = "Negative Height of Layer "+(string)castIntToString(tp1);
            setStatus(1, "Error: Layer Config");
            setStatus(2, error2);
            return;
        }
    }

    for (int i = 0; i < geometries.size() - 1; i++)
    {
        if ((geometries.at(i)[4] + geometries.at(i)[5] + 0.5) < geometries.at(i + 1)[4])
        {
            if (i + 1 == detectorLayer);
            else
            {
                tp1 = i + 1;
                tp2 = i + 2;
                string error2 = "Gap between Layers "+(string)castIntToString(tp1)+"+"+(string)castIntToString(tp2);
                setStatus(1, "Error: Layer Config");
                setStatus(2, error2);
                return;
            }
        }
    }
    /*
    for (int i = 0; i < maxLayersAllowed; i++)
    {
       additionalDetectorLayers[i] = -1;
    }

    totalAdditionalDetectorLayers = 0;

    // checking whether additional detector layers were defined
    for (int i = 0; i < geometries.size() - 1; i++)
    {
        mat = (int)geometries.at(i)[6];

        if ((mat > 500) && (mat < 600))
        {
            additionalDetectorLayers[i] = totalAdditionalDetectorLayers;

            totalAdditionalDetectorLayers++;

            mat = mat - 500;setupGeometry()
            geometries.at(i)[6] = mat;
        }
    }
    */

    for (int i = 0; i < geometries.size() - 1; i++)
    {
        mat = (int)geometries.at(i)[6];

        tp1 = i + 1;
        if ((mat < 5) || (mat > 50)) // here the materials are hardcoded... could be improved
        {
            string error2 = "At Layer "+(string)castIntToString(tp1);
            setStatus(1, "Error: Wrong Material");
            return;
        }
    }

    if (!noGUIMode) {setStatus(1, "Reading Spectrum Parameter Files");    delay(5);}

    if (loadParamData((string)endfFolder)) return;

    if (!simulationRunning)
    {
        simulationRunning = true;
        cosmicNSimulator(this);
        delay(50);
    }

    ui->neutronCountView->setText("");
    ui->label_timeRemain->setText("");
    ui->label_npers->setText("()");
}


/**
 * Function to be called on clicking the stop button.
 *
 */
void MainWindow::on_pushButton_stop_clicked()
{
    if (simulationRunning)
    {
        stopRunning = true;
        simulationRunning = false;
    }
}


/**
 * Function to be called on changing checkBoxR status, setting the log scale setting.
 * @param arg1
 *
 */
void MainWindow::on_checkBoxR_stateChanged(int arg1)
{
    if (logScaleR) logScaleR = false;
    else logScaleR = true;

    if (logScaleR)
    {
        QSharedPointer<QCPAxisTickerLog> logTicker4(new QCPAxisTickerLog);
        ui->customPlot4->yAxis->setTicker(logTicker4);
        ui->customPlot4->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else
    {
        ui->customPlot4->yAxis->setScaleType(QCPAxis::stLinear);

        QSharedPointer<QCPAxisTickerFixed> linTicker4(new QCPAxisTickerFixed);
        linTicker4->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
        linTicker4->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

        ui->customPlot4->yAxis->setTicker(linTicker4);
    }
    ui->customPlot4->update();
    ui->customPlot4->rescaleAxes();
    ui->customPlot4->replot();
}

/**
 * Function to be called on checking the log checkbox in the spatial view tab for log scale y.
 *
 *@param arg1
 */
void MainWindow::on_checkBoxD_stateChanged(int arg1)
{
    if (logScaleD) logScaleD = false;
    else logScaleD = true;

    if (logScaleD)
    {
        QSharedPointer<QCPAxisTickerLog> logTicker5(new QCPAxisTickerLog);
        ui->customPlot5->yAxis->setTicker(logTicker5);
        ui->customPlot5->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else
    {
        ui->customPlot5->yAxis->setScaleType(QCPAxis::stLinear);

        QSharedPointer<QCPAxisTickerFixed> linTicker5(new QCPAxisTickerFixed);
        linTicker5->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
        linTicker5->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

        ui->customPlot5->yAxis->setTicker(linTicker5);
    }
    ui->customPlot5->update();
    ui->customPlot5->rescaleAxes();
    ui->customPlot5->replot();
}

/**
 * Function to be called on changing the slider for
 * soil water content in the Physical parameters tab.
 *
 * @param position
 */
void MainWindow::on_sliderSoilMoisture_sliderMoved(int position)
{
    soilWaterFrac = (position * 1.) / 200.;
    soilWaterFracVar = (position * 1.) / 200.;
    string posText = castFloatToString((position * 1.) / 2., 4) + " %";
    ui->labelSM->setText(QString::fromStdString(posText));

    fpSoilMoist = (position * 1.) / 200.;

    if (activateFP) rangeIntegral = -1;

    if (activateFP) replotFootprint();
}

/**
 * Function to be called on changing the Air Humidity
 * slider in the Physical parameters tab.
 * @param position
 */
void MainWindow::on_sliderAirHum_sliderMoved(int position)
{
    absHumidityAir = (position * 1.) / 3.;
    absHumidityAirVar = absHumidityAir;

    string posText = castFloatToString(absHumidityAir, 4) + " g/m<sup>3</sup>";
    ui->labelHum->setText(QString::fromStdString(posText));

    rLuftWater = absHumidityAir / 1e6;
    fpHum = 10. / 0.6 * (position * 1.) / 50.;

    if (activateFP) rangeIntegral = -1;

    if (activateFP) replotFootprint();
}

/**
 * Function to be called on changing the
 * maximum number of neutrons for simulation.
 * @param arg1
 */
void MainWindow::on_lineEditNeutrinsTotal_textChanged(const QString& arg1)
{
    long neuronString = arg1.toLong();
    ui->lineEditNeutrinsTotal->setPalette(*paletteB);
    if ((neuronString > 0) && (neuronString < 1e12))   neutrons = neuronString;
    else ui->lineEditNeutrinsTotal->setPalette(*paletteR);
}

/**
 * Function to be called on changing the Atmospheric depth
 * slider in the Physical parameters tab.
 * @param position
 */
void MainWindow::on_sliderAtm1_sliderMoved(int position)
{
    atmDensity = position;
    string posText = (string)castIntToString(position) + " g/cm<sup>2</sup>";
    ui->labelAtm1->setText(QString::fromStdString(posText));
    pressureFac = (position * 1.) / 1020.;
    atmPressure = atmDensity * 100.;
    rLuft = atmPressure / 287.058 / sysTemperature / 1e3;
}


/**
 * Function to be called on changing the cutoff rigidity
 * slider in the Physical parameters tab.
 * @param position
 */
void MainWindow::on_sliderRigidity_sliderMoved(int position)
{
    rigidity = (position * 1.) / 10.;
    string posText = castFloatToString((position * 1.) / 10., 4);
    ui->labelrigidity->setText(QString::fromStdString(posText));
}

/**
 * Function to be called on changing the dimension of
 * the domain in the computational parameters tab.
 * @param arg1
 */
void MainWindow::on_lineEditSquareDim_textChanged(const QString& arg1)
{
    double squareDimSring = arg1.toDouble();
    ui->lineEditSquareDim->setPalette(*paletteB);
    if ((squareDimSring > 0) && (squareDimSring < 1e9))
    {
        squareDim = squareDimSring * 1000.;
        matrixStartX = squareDimSring * (-0.5);
        matrixStartY = squareDimSring * (-0.5);
        matrixMetricFactor = squareDimSring / (inputMatrixPixels * 1.); //meter per pixel
    }
    else ui->lineEditSquareDim->setPalette(*paletteR);
}

/**
 * Function to be called on changing the
 * number of neutrons for refreshing the GUI.
 * @param arg1
 */
void MainWindow::on_lineEditRefresh_textChanged(const QString& arg1)
{
    int refreshCycleString = arg1.toInt();
    ui->lineEditRefresh->setPalette(*paletteB);
    if ((refreshCycleString > 10) && (refreshCycleString < 1e8))   refreshCycle = refreshCycleString;
    else ui->lineEditRefresh->setPalette(*paletteR);
}

/**
 * Function to clear a vector.
 * @param vec
 */
void dataDelete(QVector<double>* vec)
{
    vec->clear();
}

/**
 * Function to be called on clicking the clear button
 * to clear the simulation.
 *
 */
void MainWindow::on_pushButton_Clear_clicked()
{
    if (alreadyStarted)
    {
        if (!noGUIMode) setStatus(1, "Clearing Data");
        delay(15);
        if (liveTHs.size() > 1) histoLiveClearUp(&liveTHs);

        ui->neutronCountView->setText("");
        ui->label_timeRemain->setText("");
        ui->progressBar->setValue(1);
        ui->label_npers->setText("()");
        ui->neutronCountView->setText("");
        ui->label_detectorNs->setText("()");
        if (!noGUIMode) setStatus(1, "");
        if (!noGUIMode) setStatus(2, "");
        redrawNeutronMap(-1);

        ui->detectorLayerRatioBar->setValue(0);
        ui->detectorRatioBar->setValue(0);
    }
}

/**
 * Function to be called on clicking the pause button
 * to pause the simulation.
 *
 */
void MainWindow::on_pushButton_Pause_clicked()
{
    if (simulationRunning)
    {
        if (pausehere)
        {
            ui->pushButton_Pause->setChecked(false); ui->pushButton_Pause->setDown(false); setStatus(1, "");
            pausehere = false;
            stopMode = false;
        }
        else
        {
            ui->pushButton_Pause->setChecked(true);  ui->pushButton_Pause->setDown(true); setStatus(1, "Paused");
            pausehere = true;
            stopMode = true;
        }
    }
}

/**
 * Function to be called on change of the
 * radius of the detector in the detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditDetRad_textChanged(const QString& arg1)
{
    float detRadString = arg1.toFloat();
    ui->lineEditDetRad->setPalette(*paletteB);
    if ((detRadString > -0.01) && (detRadString < 1e4))   detRad = 1000. * detRadString;
    else ui->lineEditDetRad->setPalette(*paletteR);
}

/**
 * Function to be called on change of the x coordinate
 * of the position of the detector in the detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditDetX_textChanged(const QString& arg1)
{
    float detPosString = arg1.toFloat();
    ui->lineEditDetX->setPalette(*paletteB);
    if ((detPosString > -1e5) && (detPosString < 1e5)) { detPosX[0] = 1000. * detPosString; }
    else ui->lineEditDetX->setPalette(*paletteR);
}

/**
 * Function to be called on change of the y coordinate
 * of the position of the detector in the detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditDety_textChanged(const QString& arg1)
{
    float detPosString = arg1.toFloat();
    ui->lineEditDety->setPalette(*paletteB);
    if ((detPosString > -1e5) && (detPosString < 1e5))   detPosY[0] = 1000. * detPosString;
    else ui->lineEditDety->setPalette(*paletteR);
}

/**
 * Function to be called on clicking the Round radio button for the
 * source geometry in the computational parameter tab.
 *
 */
void MainWindow::on_beamRound_clicked()
{
    ui->beamRound->setChecked(true);
    ui->beamSquare->setChecked(false);
    useRadialBeam = true;
    useRectShape = false;
    ui->label_21->setText("Radius");
}

/**
 * Function to be called on clicking the Quadratic radio button for the
 * source geometry in the computational parameter tab.
 *
 */
void MainWindow::on_beamSquare_clicked()
{
    ui->beamRound->setChecked(false);
    ui->beamSquare->setChecked(false);
    useRadialBeam = false;
    useRectShape = true;
    ui->label_21->setText("1/2 Side Length");
}

/**
 * Function to be called on changing the value of the radius in
 * source geometry in the computational parameter tab.
 * @param arg1
 */
void MainWindow::on_lineEditBeamRad_textChanged(const QString& arg1)
{
    double beamRadiusString = arg1.toDouble();
    ui->lineEditBeamRad->setPalette(*paletteB);
    if ((beamRadiusString > 0) && (beamRadiusString < squareDim / 2.))   beamRadius = beamRadiusString * 1000.;
    else ui->lineEditBeamRad->setPalette(*paletteR);
}

/**
 * Function to be called on clicking the river, width radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 *
 */
void MainWindow::on_radioRiver_clicked()
{
    if (!islandSetup)
    {
        islandSetup = true;
        riverSetup = false;
        coastSetup = false;
        lakeSetup = false;
    }
    else
    {
    }
}

/**
 * Function to be called on clicking the coast at x radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 *
 */
void MainWindow::on_radioCoast_clicked()
{
    if (!coastSetup)
    {
        islandSetup = false;
        riverSetup = false;
        coastSetup = true;
        lakeSetup = false;
    }
    else
    {
    }
}

/**
 * Function to be called on clicking the Island, diameter [m] radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 *
 *
 */
void MainWindow::on_radioIsland_clicked()
{
    if (!islandSetup)
    {
        islandSetup = true;
        riverSetup = false;
        coastSetup = false;
        lakeSetup = false;
    }
    else
    {
    }
}

/**
 * Function to be called on clicking the Lake, diameter [m] radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 *
 *
 */
void MainWindow::on_radioLake_clicked()
{
    if (!lakeSetup)
    {
        islandSetup = false;
        riverSetup = false;
        coastSetup = false;
        lakeSetup = true;
    }
    else
    {
    }

}

/**
 * Function to be called on clicking the none radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 */
void MainWindow::on_radioButton_clicked()
{
    islandSetup = false;
    riverSetup = false;
    coastSetup = false;
    lakeSetup = false;
}


/**
 * Function to be called on toggling the river,width radio button
 * for the Topological presets (water, land) in the computational parameter tab.
 * @param checked
 */
void MainWindow::on_radioRiver_toggled(bool checked)
{
}

/**
 * Function to be called on changing the value of the Neutron-sensitive energy higher band limit
 * in the detector and detector layer in the detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditTHLhigh_textChanged(const QString& arg1)
{
    float energyString = arg1.toFloat();
    ui->lineEditTHLhigh->setPalette(*paletteB);
    if ((energyString > 0) && (energyString < 20) && (energyString > energylowTHL))   energyhighTHL = energyString;
    else ui->lineEditTHLhigh->setPalette(*paletteR);
}

/**
 * Function to be called on changing the value of the Neutron-sensitive energy lower band limit
 * in the detector and detector layer in the detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditTHLlow_textChanged(const QString& arg1)
{
    float energyString = arg1.toFloat();
    ui->lineEditTHLlow->setPalette(*paletteB);
    if ((energyString > 0) && (energyString < 20) && (energyString < energyhighTHL))   energylowTHL = energyString;
    else ui->lineEditTHLlow->setPalette(*paletteR);
}

/**
 * Function to be called on clicking the log checkbox in the Range view Tab.
 *
 */
void MainWindow::on_checkBoxRS_clicked()
{
    if (logScaleRS) logScaleRS = false;
    else logScaleRS = true;

    if (logScaleRS)
    {
        QSharedPointer<QCPAxisTickerLog> logTicker6(new QCPAxisTickerLog);
        ui->customPlot6->yAxis->setTicker(logTicker6);
        ui->customPlot6->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else
    {
        ui->customPlot6->yAxis->setScaleType(QCPAxis::stLinear);

        QSharedPointer<QCPAxisTickerFixed> linTicker6(new QCPAxisTickerFixed);
        linTicker6->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssReadability);
        linTicker6->setScaleStrategy(QCPAxisTickerFixed::ssMultiples);

        ui->customPlot6->yAxis->setTicker(linTicker6);
    }
    ui->customPlot6->update();
    ui->customPlot6->rescaleAxes();
    ui->customPlot6->replot();
}

/**
 * Function to be called on changing the value of the river width,
 * for the Topological presets (water, land) in the computational parameter tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_River_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_River->setPalette(*paletteB);
    if ((valueString > 0) && (valueString < 10000) && (valueString < squareDim / 1000.))   riverDiameter = valueString * 1000.;
    else ui->lineEdit_River->setPalette(*paletteR);
}

/**
 * Function to be called on changing the value of the Coast at x [m],
 * for the Topological presets (water, land) in the computational parameter tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_River_2_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_River_2->setPalette(*paletteB);
    if ((valueString > -squareDim / 2000.) && (valueString < squareDim / 2000.))  coastPosition = valueString * 1000.;
    else ui->lineEdit_River_2->setPalette(*paletteR);
}

/**
 * Function to be called on changing the value of the Island, diameter [m],
 * for the Topological presets (water, land) in the computational parameter tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_Island_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_Island->setPalette(*paletteB);
    if ((valueString > 0) && (valueString < 10000) && (valueString < squareDim / 1000.))   islandDiameter = valueString * 1000.;
    else ui->lineEdit_Island->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_Lake_selectionChanged()
{
}

/**
 * Function to be called on changing the value of the Lake, diameter [m],
 * for the Topological presets (water, land) in the computational parameter tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_Lake_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_Lake->setPalette(*paletteB);
    if ((valueString > 0) && (valueString < 10000) && (valueString < squareDim / 1000.))  lakeDiameter = valueString * 1000.;
    else ui->lineEdit_Lake->setPalette(*paletteR);
}

/**
 * Function to be called on clicking the Write Detector Neutron Hits
 * to File checkbox 7for the individual neutron data in the Export tab.
 *
 */
void MainWindow::on_checkBoxFileOutput_clicked()
{
}

/**
 * Function to be called to set the mouse pointer coordinates when clicked on the map in the Bird's Eyes view.
 * @param target
 * @param event
 */
bool MainWindow::eventFilter(QObject* target, QEvent* event)
{
    if (target == ui->customPlot2 && event->type() == QEvent::MouseButtonPress)
    {
        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);

        if (_mouseEvent->button() == Qt::MiddleButton)
        {
            godzillaMode = true;
        }
        if ((pausehere) && (stopMode)) pausehere = false;
    }

    if (target == ui->customPlot2 && event->type() == QEvent::MouseMove)
    {

        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);

        float x = ui->customPlot2->xAxis->pixelToCoord(_mouseEvent->pos().x());
        float y = ui->customPlot2->yAxis->pixelToCoord(_mouseEvent->pos().y());

        xCustomPos = x * 1000.;
        yCustomPos = y * 1000.;
    }

    if (target == ui->customPlot2 && event->type() == QEvent::MouseButtonRelease)
    {
        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);
        if (_mouseEvent->button() == Qt::MiddleButton)
        {
            godzillaMode = false;
        }
        if (stopMode)  pausehere = true;
    }

    return false;
}

/**
 * Function to redraw the Side View for the map, deprecated.
 *
 */
void MainWindow::redrawSideView()
{
    domainLowDrawEdge = domainLowEdge;
    domainUpperDrawEdge = domainUpperEdge;
    domainZLowDrawEdge = domainZLowEdge;
    domainZUpperDrawEdge = domainZUpperEdge;
}

/**
 * Function to redraw the Top View for the neutron density map.
 *
 */
void MainWindow::redrawTopView()
{
    if ((ui->tabWidget_live->currentIndex() == 3) || (!simulationRunning))
    {
        int elementCount = ui->customPlot10->plotLayout()->elementCount();

        for (int i = 0; i < elementCount; i++)
        {
            if (qobject_cast<QCPLayoutGrid*>(ui->customPlot10->plotLayout()->elementAt(i)))  ui->customPlot10->plotLayout()->removeAt(i);
        }
        ui->customPlot10->plotLayout()->simplify();

        QCPLayoutGrid* subLayout = new QCPLayoutGrid;
        ui->customPlot10->plotLayout()->addElement(1, 0, subLayout);
        subLayout->setMargins(QMargins(53, 0, 16, 5));
        QCPColorScale* colorScale = new QCPColorScale(ui->customPlot10);
        subLayout->addElement(0, 0, colorScale);
        colorScale->setType(QCPAxis::atBottom);

        QColor lightBlue = QColor(215, 237, 255);
        colorScale->axis()->setSubTickPen(QPen(lightBlue));
        colorScale->axis()->setTickPen(QPen(lightBlue));
        colorScale->axis()->setTickLabelColor(lightBlue);
        colorScale->axis()->setLabelColor(lightBlue);

        ui->customPlot10->clearPlottables();

        float sliderDetectorValue = ui->horizontalSliderDetector->value();
        float rangeValue = 2. + sliderDetectorValue / 200. * 500.; //in m

        const int unitCells = 30;

        if (ui->horizontalSliderDetector->value() <= 1)
        {
            QCPColorMap* colorMapDetector = new QCPColorMap(ui->customPlot10->xAxis, ui->customPlot10->yAxis);

            colorMapDetector->data()->setSize(250, 250);
            colorMapDetector->data()->setRange(QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001), QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001));

            for (int x = 0; x < 500; x += 2)
            {
                for (int y = 0; y < 500; y += 2)
                {
                    colorMapDetector->data()->setCell((int)x / 2, (int)y / 2, ::detectorOriginMap->GetBinContent(x, y) + ::detectorOriginMap->GetBinContent(x + 1, y) + ::detectorOriginMap->GetBinContent(x, y + 1) + ::detectorOriginMap->GetBinContent(x + 1, y + 1));
                }
            }

            colorMapDetector->setDataRange(QCPRange(0, 2.2));
            colorMapDetector->setGradient(QCPColorGradient::gpCold);

            string detectorLabelMaxText = "-";
            string detectorLabelText = "-";

            ui->label_DetectorRangeMaximum->setText(QString::fromStdString(detectorLabelMaxText));
            ui->label_DetectorRange->setText(QString::fromStdString(detectorLabelText));

            colorScale->setDataRange(QCPRange(0, 2.2));
            colorScale->setGradient(QCPColorGradient::gpCold);
        }
        else
        {
            QCPColorMap* colorMapDetectorRelativeSquare = new QCPColorMap(ui->customPlot10->xAxis, ui->customPlot10->yAxis);

            colorMapDetectorRelativeSquare->data()->setSize(250, 250);
            colorMapDetectorRelativeSquare->data()->setRange(QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001), QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001));

            float allDetectorNeutrons = detectorOriginMap->GetEntries();
            float detectororiginMapValue;
            float detectororiginMapDrawValue, detectororiginMapDrawValueMax = 0;
            bool mapValueFound = false;

            float pixX, pixY, sumMatrixX, sumMatrixY, sumMatrixXp1, sumMatrixYp1;

            TMatrixF sumMatrix(unitCells, unitCells);

            for (int x = 0; x < 500; x += 1)
            {
                pixX = (-squareDim * 0.5 + squareDim * 0.002 * x - detPosX[0]) * 0.001; //relative to Det in m
                for (int y = 0; y < 500; y += 1)
                {
                    detectororiginMapValue = detectorOriginMap->GetBinContent(x, y);
                    if (detectororiginMapValue < 0.1) continue;

                    pixY = (-squareDim * 0.5 + squareDim * 0.002 * y - detPosY[0]) * 0.001;

                    for (int m = 0; m < unitCells; m += 1)
                    {
                        mapValueFound = false;
                        sumMatrixX = -rangeValue + 2. * rangeValue * m / (unitCells * 1.);
                        sumMatrixXp1 = sumMatrixX + 2. * rangeValue / (unitCells * 1.);
                        if ((pixX >= sumMatrixX) && (pixX < sumMatrixXp1))
                        {
                            for (int n = 0; n < unitCells; n += 1)
                            {
                                sumMatrixY = -rangeValue + 2. * rangeValue * n / (unitCells * 1.);
                                sumMatrixYp1 = sumMatrixY + 2. * rangeValue / (unitCells * 1.);
                                if ((pixY >= sumMatrixY) && (pixY < sumMatrixYp1))
                                {
                                    sumMatrix(m, n) = sumMatrix(m, n) + detectororiginMapValue;
                                    mapValueFound = true;
                                    break;
                                }
                            }
                        }
                        if (mapValueFound) break;
                    }
                }
            }

            for (int x = 0; x < 250; x += 1)
            {
                pixX = (-squareDim * 0.5 + squareDim * 0.004 * x - detPosX[0]) * 0.001; //relative to Det in m
                for (int y = 0; y < 250; y += 1)
                {
                    detectororiginMapDrawValue = 0;

                    pixY = (-squareDim * 0.5 + squareDim * 0.004 * y - detPosY[0]) * 0.001;

                    for (int m = 0; m < unitCells; m += 1)
                    {
                        mapValueFound = false;
                        sumMatrixX = -rangeValue + 2. * rangeValue * m / (unitCells * 1.);
                        sumMatrixXp1 = sumMatrixX + 2. * rangeValue / (unitCells * 1.);
                        if ((pixX >= sumMatrixX) && (pixX < sumMatrixXp1))
                        {
                            for (int n = 0; n < unitCells; n += 1)
                            {
                                sumMatrixY = -rangeValue + 2. * rangeValue * n / (unitCells * 1.);
                                sumMatrixYp1 = sumMatrixY + 2. * rangeValue / (unitCells * 1.);
                                if ((pixY >= sumMatrixY) && (pixY < sumMatrixYp1))
                                {
                                    detectororiginMapDrawValue = sumMatrix(m, n);
                                    if (detectororiginMapDrawValue > detectororiginMapDrawValueMax) detectororiginMapDrawValueMax = detectororiginMapDrawValue;
                                    mapValueFound = true;
                                    break;
                                }
                            }
                        }
                        if (mapValueFound) break;
                    }
                    colorMapDetectorRelativeSquare->data()->setCell((int)x, (int)y, detectororiginMapDrawValue / allDetectorNeutrons);
                }
            }

            float detectorDataScaleRange = 1.05 * (0.01 + ui->horizontalSliderDetectorColor->value() / 100. * detectororiginMapDrawValueMax / allDetectorNeutrons);
            colorMapDetectorRelativeSquare->setDataRange(QCPRange(0, detectorDataScaleRange));
            colorMapDetectorRelativeSquare->setGradient(QCPColorGradient::gpJet);

            int lltext1 = detectorDataScaleRange * 100;
            int lltext2 = rangeValue;
            string detectorLabelMaxText = castIntToString(lltext1) + " %";
            string detectorLabelText = castIntToString(lltext2) + " m";

            ui->label_DetectorRangeMaximum->setText(QString::fromStdString(detectorLabelMaxText));
            ui->label_DetectorRange->setText(QString::fromStdString(detectorLabelText));

            colorScale->setDataRange(QCPRange(0, 100 * detectorDataScaleRange));
            colorScale->setGradient(QCPColorGradient::gpJet);
        }

        ui->customPlot10->rescaleAxes();
        ui->customPlot10->update();
        ui->customPlot10->replot();
    }

    if ((ui->tabWidget_live->currentIndex() == 0) || (!simulationRunning))
    {
        ui->customPlot2->clearPlottables();

        QCPColorMap* colorMap = new QCPColorMap(ui->customPlot2->xAxis, ui->customPlot2->yAxis);

        if (plotTopViewLog) { colorMap->setDataScaleType(QCPAxis::stLogarithmic); }

        colorMap->data()->setSize(250, 250);
        colorMap->data()->setRange(QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001), QCPRange(-::squareDim * 0.5 * 0.001, ::squareDim * 0.5 * 0.001));

        double allEntries = 0;
        long allEntriesInt = 0;
        double allEntriesMaximum = 0;
        float actualEntry = 0;
        float minZ = 0;

        if (plotTopViewLog) minZ = 1;

        if (ui->checkBoxGradient2->isChecked())
        {
            TH2F* densityMapCopy = new TH2F();
            TH2F* densityMapRed = new TH2F();

            if (showDensityMap) { densityMapRed = (TH2F*)densityMap->Clone(""); densityMapCopy = (TH2F*)densityMap->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }
            if (showDensityMapIntermediate) { densityMapRed = (TH2F*)densityMapIntermediate->Clone(""); densityMapCopy = (TH2F*)densityMapIntermediate->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }
            if (showDensityMapFast) { densityMapRed = (TH2F*)densityMapFast->Clone(""); densityMapCopy = (TH2F*)densityMapFast->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }
            if (showDensityMapAlbedo) { densityMapRed = (TH2F*)densityMapAlbedo->Clone(""); densityMapCopy = (TH2F*)densityMapAlbedo->Clone(""); densityMapCopy->RebinX(2); densityMapCopy->RebinY(2); densityMapRed->RebinX(2); densityMapRed->RebinY(2); getGradientMatrixFromTH2(densityMapRed, densityMapCopy); }

            for (int x = 0; x < 250; x += 1)
            {
                for (int y = 0; y < 250; y += 1)
                {
                    actualEntry = densityMapCopy->GetBinContent(x, y);
                    colorMap->data()->setCell((int)x, (int)y, actualEntry);
                    allEntries += actualEntry;
                    if (actualEntry > allEntriesMaximum) allEntriesMaximum = actualEntry;
                }
            }

            entryWeight = (allEntries * 1.) / 250. / 250. * 2.; maximumWeight = 2. * allEntriesMaximum;

            delete densityMapCopy;
            delete densityMapRed;
        }
        else
        {
            for (int x = 0; x < 500; x += 2)
            {
                for (int y = 0; y < 500; y += 2)
                {
                    if (showDensityTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityTrackMap->GetBinContent(x, y) + ::densityTrackMap->GetBinContent(x + 1, y) + ::densityTrackMap->GetBinContent(x, y + 1) + ::densityTrackMap->GetBinContent(x + 1, y + 1));
                    if (showDensityIntermediateTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityIntermediateTrackMap->GetBinContent(x, y) + ::densityIntermediateTrackMap->GetBinContent(x + 1, y) + ::densityIntermediateTrackMap->GetBinContent(x, y + 1) + ::densityIntermediateTrackMap->GetBinContent(x + 1, y + 1));
                    if (showDensityFastTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityFastTrackMap->GetBinContent(x, y) + ::densityFastTrackMap->GetBinContent(x + 1, y) + ::densityFastTrackMap->GetBinContent(x, y + 1) + ::densityFastTrackMap->GetBinContent(x + 1, y + 1));
                    if (showDensityAlbedoTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityAlbedoTrackMap->GetBinContent(x, y) + ::densityAlbedoTrackMap->GetBinContent(x + 1, y) + ::densityAlbedoTrackMap->GetBinContent(x, y + 1) + ::densityAlbedoTrackMap->GetBinContent(x + 1, y + 1));

                    if (showDensityEnergyTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityEnergyTrackMap->GetBinContent(x, y) + ::densityEnergyTrackMap->GetBinContent(x + 1, y) + ::densityEnergyTrackMap->GetBinContent(x, y + 1) + ::densityEnergyTrackMap->GetBinContent(x + 1, y + 1));

                    if (showDensityThermalTrackMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityThermalTrackMap->GetBinContent(x, y) + ::densityThermalTrackMap->GetBinContent(x + 1, y) + ::densityThermalTrackMap->GetBinContent(x, y + 1) + ::densityThermalTrackMap->GetBinContent(x + 1, y + 1));
                    if (showDensityMapThermal)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapThermal->GetBinContent(x, y) + ::densityMapThermal->GetBinContent(x + 1, y) + ::densityMapThermal->GetBinContent(x, y + 1) + ::densityMapThermal->GetBinContent(x + 1, y + 1));

                    if (showDensityMap)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMap->GetBinContent(x, y) + ::densityMap->GetBinContent(x + 1, y) + ::densityMap->GetBinContent(x, y + 1) + ::densityMap->GetBinContent(x + 1, y + 1));
                    if (showDensityMapIntermediate)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapIntermediate->GetBinContent(x, y) + ::densityMapIntermediate->GetBinContent(x + 1, y) + ::densityMapIntermediate->GetBinContent(x, y + 1) + ::densityMapIntermediate->GetBinContent(x + 1, y + 1));
                    if (showDensityMapFast)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapFast->GetBinContent(x, y) + ::densityMapFast->GetBinContent(x + 1, y) + ::densityMapFast->GetBinContent(x, y + 1) + ::densityMapFast->GetBinContent(x + 1, y + 1));
                    if (showDensityMapAlbedo)    colorMap->data()->setCell((int)x / 2, (int)y / 2, ::densityMapAlbedo->GetBinContent(x, y) + ::densityMapAlbedo->GetBinContent(x + 1, y) + ::densityMapAlbedo->GetBinContent(x, y + 1) + ::densityMapAlbedo->GetBinContent(x + 1, y + 1));

                }
            }

            if (showDensityTrackMap) { entryWeight = (densityTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityTrackMap->GetBinContent(densityTrackMap->GetMaximumBin()); }
            if (showDensityIntermediateTrackMap) { entryWeight = (densityIntermediateTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityIntermediateTrackMap->GetBinContent(densityIntermediateTrackMap->GetMaximumBin()); }
            if (showDensityFastTrackMap) { entryWeight = (densityFastTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityFastTrackMap->GetBinContent(densityFastTrackMap->GetMaximumBin()); }
            if (showDensityAlbedoTrackMap) { entryWeight = (densityAlbedoTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityAlbedoTrackMap->GetBinContent(densityAlbedoTrackMap->GetMaximumBin()); }

            if (showDensityEnergyTrackMap) { entryWeight = (densityEnergyTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityEnergyTrackMap->GetBinContent(densityEnergyTrackMap->GetMaximumBin()); }

            if (showDensityThermalTrackMap) { entryWeight = (densityThermalTrackMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityThermalTrackMap->GetBinContent(densityThermalTrackMap->GetMaximumBin()); }
            if (showDensityMapThermal) { entryWeight = (densityMapThermal->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityMapThermal->GetBinContent(densityMapThermal->GetMaximumBin()); }

            if (showDensityMap) { entryWeight = (densityMap->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityMap->GetBinContent(densityMap->GetMaximumBin()); }
            if (showDensityMapIntermediate) { entryWeight = (densityMapIntermediate->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityMapIntermediate->GetBinContent(densityMapIntermediate->GetMaximumBin()); }
            if (showDensityMapFast) { entryWeight = (densityMapFast->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityMapFast->GetBinContent(densityMapFast->GetMaximumBin()); }
            if (showDensityMapAlbedo) { entryWeight = (densityMapAlbedo->GetEntries()) / 250. / 250. * 4.; maximumWeight = 4. * densityMapAlbedo->GetBinContent(densityMapAlbedo->GetMaximumBin()); }

            allEntriesInt = entryWeight * 250. * 250. * 4.;
            //string numberofDetectorLayerNs = castLongToString(allEntriesInt);

            totalRealTime = (allEntriesInt*1.)/neutronRealScalingFactor;
            string numberofDetectorLayerNs;
            if (totalRealTime < 1)  numberofDetectorLayerNs = castDoubleToString(totalRealTime*1000.,6)+" ms";
            else  numberofDetectorLayerNs = castDoubleToString(totalRealTime,6)+" s";
            ui->label_detectorLayerNs->setText(QString::fromStdString(numberofDetectorLayerNs));

            allEntriesInt = densityMapHighEnergy->GetEntries();
            if (useExtraCounter) ui->label_detectorLayerNs2->setText(QString::fromStdString((string)castLongToString(allEntriesInt)));
        }

        float colorscaleMax = (ui->horizontalSliderColor->value() * 1.) / 200. * maximumWeight;
        if (colorscaleMax < 1) colorscaleMax = 1;

        colorMap->setGradient(QCPColorGradient::gpJet);

        if (ui->radioButton_NeutronNight->isChecked()) colorMap->setGradient(QCPColorGradient::gpNight);
        if (ui->radioButton_NeutronCold->isChecked()) colorMap->setGradient(QCPColorGradient::gpCold);
        if (ui->radioButton_NeutronPolar->isChecked()) colorMap->setGradient(QCPColorGradient::gpThermal);
        if (ui->radioButton_NeutronHot->isChecked()) colorMap->setGradient(QCPColorGradient::gpHot);
        if (ui->radioButton_NeutronThermal->isChecked()) colorMap->setGradient(QCPColorGradient::gpPolar);
        if (ui->radioButton_NeutronGrayScale->isChecked()) colorMap->setGradient(QCPColorGradient::gpGrayscale);

        if ((!silderColorMoved) && ((!ui->checkBoxClearEveryDisplayRefresh->isChecked()) || (!ui->checkBoxSaveEvery_2->isChecked())))
        {
            if (entryWeight * 1.1 < 5) colorMap->setDataRange(QCPRange(minZ, 5));
            else colorMap->setDataRange(QCPRange(minZ, entryWeight * 1.1));
        }
        else
        {
            if (ui->horizontalSliderColorZero->value() == 0) colorMap->setDataRange(QCPRange(minZ, colorscaleMax));
            else colorMap->setDataRange(QCPRange((ui->horizontalSliderColorZero->value() * 1.) / 200. * colorscaleMax, colorscaleMax));
        }

        if (manualColorZero < minZ) manualColorZero = minZ;

        int lowColorValue = (ui->horizontalSliderColorZero->value() * 1.) / 200. * colorscaleMax;
        int highColorValue = colorscaleMax;

        if (ui->checkBoxManual->isChecked())
        {
            useManualColors = true;
            colorMap->setDataRange(QCPRange(manualColorZero, manualColor));
        }
        else
        {
            useManualColors = false;
            ui->lineEditManualColorZero->setText(QString::fromStdString((string)castIntToString(lowColorValue)));
            ui->lineEditManualColor->setText(QString::fromStdString((string)castIntToString(highColorValue)));
        }

        string colorLabelText = "["+(string)castIntToString(lowColorValue)+"-"+(string)castIntToString(highColorValue)+"]";
        ui->label_ColorRange->setText(QString::fromStdString(colorLabelText));

        ui->customPlot2->update();
        ui->customPlot2->replot();

        if (visualization != 0x0)
        {
            if (visualization->isVisible()) updateEnlargedView = true;
            else  updateEnlargedView = false;
        }
        if (visualization2 != 0x0)
        {
            if (visualization2->isVisible()) updateEnlargedView2 = true;
            else  updateEnlargedView2 = false;
        }
        if (updateEnlargedView)
        {
            redrawEnlargedView();
        }
        if (updateEnlargedView2)
        {
            redrawEnlargedView2();
        }
    }

    if ((ui->tabWidget_live->currentIndex() == 2) || (!simulationRunning))
    {
        if (showDensityTrackMap)            cutView = densityTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityIntermediateTrackMap)cutView = densityIntermediateTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityFastTrackMap)        cutView = densityFastTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityAlbedoTrackMap)      cutView = densityAlbedoTrackMap->ProjectionX("proj", 240, 260);

        if (showDensityEnergyTrackMap)      cutView = densityEnergyTrackMap->ProjectionX("proj", 240, 260);

        if (showDensityThermalTrackMap)     cutView = densityThermalTrackMap->ProjectionX("proj", 240, 260);
        if (showDensityMapThermal)          cutView = densityMapThermal->ProjectionX("proj", 240, 260);

        if (showDensityMap)                 cutView = densityMap->ProjectionX("proj", 240, 260);
        if (showDensityMapIntermediate)     cutView = densityMapIntermediate->ProjectionX("proj", 240, 260);
        if (showDensityMapFast)             cutView = densityMapFast->ProjectionX("proj", 240, 260);
        if (showDensityMapAlbedo)           cutView = densityMapAlbedo->ProjectionX("proj", 240, 260);

        for (int i = 0; (i < cutView->GetNbinsX()) && (i < 500); ++i)
        {
            ::plotGUIxBinsCutView[i] = cutView->GetBinCenter(i + 1);
            ::plotGUIyBinsCutView[i] = cutView->GetBinContent(i + 1);
        }

        ui->customPlot3->graph(0)->setData(::plotGUIxBinsCutView, ::plotGUIyBinsCutView);
        ui->customPlot3->update();
        ui->customPlot3->rescaleAxes();
        ui->customPlot3->replot();
    }
}

/**
 * Function to be called on click of radio Button for the
 * 1 eV-1 keV option for Hit Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_map_clicked()
{
    if (!showDensityMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMapFast = false;
        showDensityMapIntermediate = false;
        showDensityMapAlbedo = false;
        showDensityMap = true;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 1;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the
 * 1 keV-0.5 MeV option for Hit Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapInter_clicked()
{
    if (!showDensityMapIntermediate)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMapFast = false;
        showDensityMap = false;
        showDensityMapAlbedo = false;
        showDensityMapIntermediate = true;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 2;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the
 * 0.5 MeV-10 MeV option for Hit Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapFast_clicked()
{
    if (!showDensityMapFast)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapAlbedo = false;
        showDensityMapFast = true;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 3;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the
 * Detector Selection option for Hit Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapAlbedo_clicked()
{
    if (!showDensityMapAlbedo)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = true;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 4;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the 1 eV-1 keV option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrack_clicked()
{
    if (!showDensityTrackMap)
    {
        showDensityTrackMap = true;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 11;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the 1 keV-0.5 MeV option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrackInter_clicked()
{
    if (!showDensityIntermediateTrackMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = true;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 12;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the 0.5 MeV-10 MeV option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrackFast_clicked()
{
    if (!showDensityFastTrackMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = true;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 13;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the Detector Selection option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrackAlbedo_clicked()
{
    if (!showDensityAlbedoTrackMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = true;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 14;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the Thermal option
 * for Hit Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapThermal_clicked()
{
    if (!showDensityMapThermal)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = true;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 5;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the Thermal option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrackThermal_clicked()
{
    if (!showDensityThermalTrackMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = true;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = false;

        densityMapButtonID = 15;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on click of radio Button for the Energy Dependent option
 * for Track Density in the Display energy window in the Display tab.
 *
 */
void MainWindow::on_radioButton_mapTrackEnergy_clicked()
{
    if (!showDensityEnergyTrackMap)
    {
        showDensityTrackMap = false;
        showDensityIntermediateTrackMap = false;
        showDensityFastTrackMap = false;
        showDensityAlbedoTrackMap = false;

        showDensityThermalTrackMap = false;
        showDensityMapThermal = false;

        showDensityMap = false;
        showDensityMapIntermediate = false;
        showDensityMapFast = false;
        showDensityMapAlbedo = false;

        showDensityEnergyTrackMap = true;

        densityMapButtonID = 10;
    }

    if (alreadyStarted)  redrawTopView();
}

/**
 * Function to be called on change the path of the
 * Input Spectrum Calculation File textbox in the Folder tab.
 *
 */
void MainWindow::on_lineEdit_InputSpectrumFolder_editingFinished()
{
    QString valueString = ui->lineEdit_InputSpectrumFolder->text();
    string textHere = valueString.toStdString();

    std::replace(textHere.begin(), textHere.end(), '\\', '/');
    TString textHere2 = textHere;

    TFile f_here(textHere2);

    ui->lineEdit_InputSpectrumFolder->setPalette(*paletteB);

    if (f_here.IsZombie())
    {
        ifstream f(textHere.c_str());

        ui->lineEdit_InputSpectrumFolder->setPalette(*paletteB);

        if (!(f.good()))
        {
            if ((textHere == "default") || (textHere == "") || (textHere == "N/A") || (textHere == "n/a"))
            {
                ui->lineEdit_InputSpectrumFolder->setPalette(*paletteGray);
                inputSpectrumFile = "";
            }
            else
            {
                ui->lineEdit_InputSpectrumFolder->setPalette(*paletteR);
            }
        }
        else
        {
            if (textHere.length() < 7)  ui->lineEdit_InputSpectrumFolder->setPalette(*paletteR);
            else
            {
                inputSpectrumFile = textHere;
                ui->label_16->setText("Input Spectrum Calculation File (ASCII)");
            }
        }
    }
    else
    {
        inputSpectrumFile = textHere;
        ui->label_16->setText("Input Spectrum Calculation File (ROOT)");
    }
}

/**
 * Function to be called on change the path of the
 * Cross Section Folder (containing the ENDF databases)textbox in the Folder tab.
 *
 */
void MainWindow::on_lineEdit_CrosssectionFolder_editingFinished()
{
    QString valueString = ui->lineEdit_CrosssectionFolder->text();
    string textHere = valueString.toStdString();
    std::replace(textHere.begin(), textHere.end(), '\\', '/');

    ifstream input_streamCheck(textHere + "absorbH1.txt", ios::in);

    ui->lineEdit_CrosssectionFolder->setPalette(*paletteB);

    if (input_streamCheck.fail())  ui->lineEdit_CrosssectionFolder->setPalette(*paletteR);
    else
    {
        endfFolder = textHere;
    }
}

/**
 * Function to be called on change the path of the
 * Output Folder textbox in the Folder tab.
 *
 */
void MainWindow::on_lineEdit_OutputFolder_editingFinished()
{
    QString valueString = ui->lineEdit_OutputFolder->text();
    string textHere = valueString.toStdString();
    std::replace(textHere.begin(), textHere.end(), '\\', '/');

    if (access(textHere.c_str(), 0) == 0)
    {
        ui->lineEdit_OutputFolder->setPalette(*paletteB);
        outputFolder = textHere;
    }
    else
    {
        ui->lineEdit_OutputFolder->setPalette(*paletteR);
    }
}

/**
 * Function to be called on clicking the URANOS icon
 * displays details about the application.
 *
 */
void MainWindow::on_pushButton_about_clicked()
{
    QString messageString = "<p align='center'> This is <br> <br>";
    messageString += "<big>URANOS</big> <br>";
    messageString += "- <br> ";
    messageString += "The Ultra Rapid Neutron-Only Simulation <br> ";
    messageString += "For Environmental Research<br><br> ";

    messageString += "Physikalisches Institut, Heidelberg University, Germany <br>";
    messageString += "UFZ — Helmholtz Centre for Environmental Research, Leipzig, Germany <br> <br> <br> ";
    messageString += "Developed by <br> ";
    messageString += "Markus Köhli (PI Heidelberg)<br> ";
    messageString += "koehli@physi.uni-heidelberg.de <br>";
    messageString += "<br>";
    messageString += "In Collaboration with <br>";
    messageString += "Martin Schrön (UFZ Leipzig)<br> ";
    messageString += "martin.schroen@ufz.de <br>";
    messageString += "<br>";
    messageString += "This software is research in progress, any results should be treated with care, in case contact the author for support.<br>";
    messageString += "<small>If you use URANOS in a scientific publication, please mail a copy to the author.<br>";
    messageString += "This will help to continue the support of URANOS in the future.</small><br> ";
    messageString += "<br> <br>";
    messageString += "For technical support or questions contact<br>";
    messageString += "uranos@physi.uni-heidelberg.de <br> <br>";
    messageString += "Citation: M. Köhli et al., Geosci. Model. Dev., 16, 449-477, 2023 <br><br>";
    messageString+=        "v1.10 (27.02.2023)<br> ";
    messageString+=        "<small>Based on QT 5.14.2, ROOT 6.22.08 and QCustomPlot 2.1.1 (MSVC 2017 32bit)</small> <br>";
    messageString += "<small>(see also attached information)</small> <br><br>";

    messageString += "<font color='#e4e4e4'> All your neutrons <br> are belong to us!<\font>";

    QMessageBox about_box(this);
    about_box.setWindowTitle("About URANOS");
    about_box.setTextInteractionFlags(Qt::TextSelectableByMouse);
    about_box.setText(messageString);
    about_box.setIconPixmap(QPixmap("about.png"));
    about_box.setParent(this);

    about_box.exec();
}

/**
 * Function to be called on checking the transparent checkbox
 * - the neutrons are removed upon detection when non-transparent.
 *
 */
void MainWindow::on_checkBoxTransparent_clicked()
{
    if (ui->checkBoxTransparent->isChecked()) detectorAbsorbing = false;
    else detectorAbsorbing = true;
}

/**
 * Function to be called on checking the basic spectrum checkbox
 * - check to use a mean basic spectrum.
 *
 */
void MainWindow::on_checkBoxBasicSpectrum_clicked()
{
    if (useBasicSpectrum)
    {
        useBasicSpectrum = false;
        ui->sliderRigidity->setEnabled(true);
        int positionSlider = ui->sliderRigidity->value();
        string text = castFloatToString((positionSlider * 1.) / 10., 4);
        ui->labelrigidity->setText(QString::fromStdString(text));
    }
    else
    {
        useBasicSpectrum = true; ui->sliderRigidity->setEnabled(false);  ui->labelrigidity->setText("-");
    }
}

/**
 * Function to activate sky evaporation.
 *
 */
void MainWindow::activateSkyEvaporation()
{
    doSkyEvaporation = true;
    ui->groupBox_Evaporation->setHidden(false);
    ui->checkBoxThermal->setHidden(true);
}

/**
 * Function to activate detector batch run.
 *
 */
void MainWindow::activateDetectorBatchRun()
{
    doDetectorBatchRun = true;
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxFileOutputPDF_2->setChecked(true);
}

void MainWindow::activateDetectorBatchRun2()
{
    doDetectorBatchRun2 = true;
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxFileOutputPDF_2->setChecked(true);
}

void MainWindow::activateDetectorAngleBatchRun()
{
    doDetectorAngleBatchRun = true;
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxFileOutputPDF_2->setChecked(true);
}

/**
 * Function to activate batch run.
 *
 */
void MainWindow::activateBatchRun()
{
    doBatchRun = true;
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxFileOutputPDF_2->setChecked(true);
}

/**
 * Function to activate parameter batch run.
 *@param - parMin
 *@param - parMax
 */
void MainWindow::activateParameterBatchRun(int parMin, int parMax)
{
    doBatchRun2D = true;
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxFileOutputPDF_2->setChecked(true);

    paramMin = parMin;
    paramMax = parMax;
}

/**
 * Function to activate thermal neutron transport.
 *
 */
void MainWindow::activateThermal()
{
    noThermalRegime = false;
    maxScatterings = 1500; // maximum number of scatterings of a neutron until killed

    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxThermal->setChecked(true);

    ui->radioButton_mapThermal->setHidden(false);
    ui->radioButton_mapTrackThermal->setHidden(false);
    ui->checkBoxThermalMap->setHidden(false);
    ui->checkBoxThermalData->setHidden(false);
}

void MainWindow::setConfigFilePath(string pathtoConfigFile)
{
    if (pathtoConfigFile.length() > 1)
    {
        configFilePath = pathtoConfigFile;
        configFilePathConfigured = true;
    }
}

/**
 * Function to disable the GUI and run in command line.
 * @param - path to ConfigFile (uranos.cfg)
 *
 */
void MainWindow::disabledGUIRun(string pathtoConfigFile)
{
    cout << "GUI disabled" << endl;

    noGUIMode = true;
    configFilePath = pathtoConfigFile;

    if (pathtoConfigFile.length() > 1)
    {
        configFilePathConfigured = true;
        cout << "Using Config File " << pathtoConfigFile << endl;
        importSettings();
        //setupImport();
        on_pushButton_LoadGeometry_clicked();
        if (layerMapsImport)
        {
            useImage = true;
            on_checkBox_useImage_clicked();
        }
    }
    else
    {
        importSettings();
        //on_pushButton_ReadGeometry_clicked();
        on_pushButton_LoadGeometry_clicked();
    }

    setupGeometry();

    //cout << "Cosmic Neutron Spectrum Definition..." << endl;

    if (!(loadParamData((string)endfFolder)))
    {
        cout << "done" << endl;
        cosmicNSimulator(this);
    }
    else
    {
        cout << endl << "failed" << endl;
    }
}

/**
 * Function to activate thermal sky evaporation.
 *
 */
void MainWindow::activateThermalSkyEvaporation()
{
    noThermalRegime = false;
    doSkyEvaporation = true;
    highResCalc = true;
    highhighResCalc = true;

    maxScatterings = 1500; // maximum number of scatterings of a neutron until killed

    ui->groupBox_Evaporation->setHidden(false);
    ui->checkBoxThermal->setHidden(false);
    ui->checkBoxThermal->setChecked(true);

    ui->radioButton_mapThermal->setHidden(false);
    ui->radioButton_mapTrackThermal->setHidden(false);
    ui->checkBoxThermalMap->setHidden(false);
    ui->checkBoxThermalData->setHidden(false);

    exportHighResTrackData = true;
    ui->checkBoxHighResTrackingData->setChecked(true);
}

/**
 * Function to activate silent mode with no command line output.
 *
 */
void MainWindow::beSilent()
{
    silentMode = true;
}

/**
 * Function to be called on clicking add layer + button
 * in layer control in the Parameters Tab.
 *
 */
void MainWindow::on_pushButton_AddLayer_clicked()
{
    if (row > 0)   model->insertRow(row, QModelIndex());
    else  model->insertRow(0, QModelIndex());

    ui->spinBox_StartingLayer->setMaximum(ui->spinBox_StartingLayer->maximum() + 1);
    ui->spinBox_DetectorLayer->setMaximum(ui->spinBox_DetectorLayer->maximum() + 1);
    ui->spinBox_GroundLayer->setMaximum(ui->spinBox_GroundLayer->maximum() + 1);

    if (ui->checkBox_useImage->isChecked()) on_checkBox_useImage_clicked();
}

/**
 * Function to be called on clicking remove layer - button
 * in layer control in the Parameters Tab.
 *
 */
void MainWindow::on_pushButton_RemoveLayer_clicked()
{
    if (row > 0)   model->removeRow(row, QModelIndex());
    else   model->removeRow(0, QModelIndex());

    if (ui->checkBox_useImage->isChecked()) on_checkBox_useImage_clicked();
}

/**
 * Function to be called on clicking generate button
 * in layer control in the Parameters Tab.
 *
 */
void MainWindow::on_pushButton_ReadGeometry_clicked()
{
    setGeometry();

    QStandardItem* item;

    QSortFilterProxyModel* proxy1 = new QSortFilterProxyModel();
    proxy1->setSourceModel(model);

    int modelRows = model->rowCount();

    for (int z = 0; z < modelRows - 1; z++)
    {
        model->removeRow(0, QModelIndex());
    }

    for (int i = 0; i < geometries.size() - 1; i++)
    {
        model->insertRow(0, QModelIndex());
    }

    for (int i = 0; i < geometries.size(); i++)
    {

        QString position = QString::number((geometries.at(i)[4]) / 1000.);
        item = model->itemFromIndex(model->index(i, 0, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
        item->setText(position);

        QString height = QString::number((geometries.at(i)[5]) / 1000.);
        item = model->itemFromIndex(model->index(i, 1, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
        item->setText(height);

        QString material = QString::number(geometries.at(i)[6]);
        item = model->itemFromIndex(model->index(i, 2, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
        item->setText(material);

    }
    modelRows = model->rowCount();

    for (int k = geometries.size(); k < modelRows; k++)
    {
        model->removeRow(geometries.size(), QModelIndex());
    }

    ui->spinBox_StartingLayer->setValue(2);
    ui->spinBox_DetectorLayer->setValue(4);
    ui->spinBox_GroundLayer->setValue(6);

    ui->checkBox_useImage->setChecked(false);
    layerMapsImport = false;
    useImage = false;
}

/**
 * Function to be called on clicking save button in
 * Layer Control in the Parameters tab.
 *
 */
void MainWindow::on_pushButton_SaveGeometry_clicked()
{
    setupGeometry();

    ofstream* stream_out;
    stream_out = new ofstream(workFolder + "UranosGeometryConfig.dat", ofstream::out);

    *stream_out << (startingLayer + 1) << endl;
    *stream_out << (detectorLayer + 1) << endl;
    *stream_out << (groundLayer + 1) << endl;

    for (int i = 0; i < geometries.size(); i++)
    {
        if (geometries.at(i)[4] > 10000) {*stream_out <<  setprecision(6) << ((geometries.at(i)[4]) * 0.001);  *stream_out << setprecision(10);}
        else {*stream_out << ((geometries.at(i)[4]) * 0.001);}
        *stream_out << "\t";
        if (geometries.at(i)[5] > 10000) {*stream_out << setprecision(6) << ((geometries.at(i)[5]) * 0.001);  *stream_out << setprecision(10);}
        else { *stream_out << (geometries.at(i)[5]) * 0.001;}
        if (additionalDetectorLayers[i] >= 0) {*stream_out << "\t" << (geometries.at(i)[6] + 500) << endl;}
        else {*stream_out << "\t" << (geometries.at(i)[6]) << endl;}
    }

    stream_out->close();
}

void MainWindow::on_spinBox_StartingLayer_valueChanged(const QString& arg1)
{

}

/**
 * Function to be called on changing the source layer in
 * Layer Control in the Parameters tab.
 * @param arg1
 */
void MainWindow::on_spinBox_StartingLayer_valueChanged(int arg1)
{
    setupGeometry();
    ui->spinBox_StartingLayer->setMaximum(geometries.size());
    startingLayer = arg1 - 1;
}

/**
 * Function to be called on changing the Detector layer in
 * Layer Control in the Parameters tab.
 * @param arg1
 */
void MainWindow::on_spinBox_DetectorLayer_valueChanged(int arg1)
{
    setupGeometry();
    ui->spinBox_DetectorLayer->setMaximum(geometries.size());
    detectorLayer = arg1 - 1;
    detectorHeight = geometries.at(detectorLayer)[4] + 0.5 * geometries.at(detectorLayer)[5];
}

/**
 * Function to be called on changing the Ground layer in
 * Layer Control in the Parameters tab.
 * @param arg1
 */
void MainWindow::on_spinBox_GroundLayer_valueChanged(int arg1)
{
    setupGeometry();
    ui->spinBox_GroundLayer->setMaximum(geometries.size());
    groundLayer = arg1 - 1;
}

/**
 * Function to be called on clicking Load button in
 * Layer Control in the Parameters tab.
 *
 */
void MainWindow::on_pushButton_LoadGeometry_clicked()
{
    QStandardItem* item;

    TString line;
    int linecounter = 0;

    string geometryPath = (string)workFolder+"UranosGeometryConfig.dat";

    ifstream input_stream(geometryPath, ios::in);

    if (!access(geometryPath.c_str(), 0) == 0)
    {
        if (!silentMode)
        {
            if (noGUIMode) { cout<<"No Config File"<<endl; }
        }
        setStatus(1,"No Geometry Config File"); delay(1);
        return;
    }
    int  modelRows = model->rowCount();

    for (int z = 0; z <  modelRows; z++)
    {
        model->removeRow(0, QModelIndex());
    }

    if (!noGUIMode) {setStatus(1,"Loading Geometry File"); delay(1);}
    if (noGUIMode) { cout<<"Loading Geometry Files"<<endl; }

    int a, b, c, i;
    float posInput, heightInput, materialInput;

    while (line.ReadLine(input_stream))
    {
        istrstream stream(line.Data());

        if (linecounter == 0) stream >> a;
        if (linecounter == 1) stream >> b;
        if (linecounter == 2) stream >> c;

        i = linecounter - 3;
        if (linecounter > 2)
        {
            model->insertRow(i, QModelIndex());

            stream >> posInput;
            QString position = QString::number(posInput);
            item = model->itemFromIndex(model->index(i, 0, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
            item->setText(position);

            stream >> heightInput;
            QString height = QString::number(heightInput);
            item = model->itemFromIndex(model->index(i, 1, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
            item->setText(height);

            stream >> materialInput;
            QString material = QString::number(materialInput);
            item = model->itemFromIndex(model->index(i, 2, QModelIndex()));  item->setTextAlignment(Qt::AlignCenter);
            item->setText(material);
        }
        linecounter++;
    }

    modelRows = model->rowCount();

    for (int k = linecounter - 3; k < modelRows; k++)
    {
        model->removeRow(geometries.size(), QModelIndex());
    }

    setupGeometry();

    ui->spinBox_StartingLayer->setMaximum(geometries.size());
    ui->spinBox_DetectorLayer->setMaximum(geometries.size());
    ui->spinBox_GroundLayer->setMaximum(geometries.size());

    ui->spinBox_StartingLayer->setValue(a);
    ui->spinBox_DetectorLayer->setValue(b);
    ui->spinBox_GroundLayer->setValue(c);

    startingLayer = a - 1;
    detectorLayer = b - 1;
    groundLayer = c - 1;

    if (layerMapsImport)
    {
        ui->checkBox_useImage->setChecked(true);

        if (!noGUIMode) {setStatus(2,"Loading Input Definitions"); delay(1);}

        on_checkBox_useImage_clicked();
        checkInputPics();
    }
    else
    {
        ui->checkBox_useImage->setChecked(false);
    }

    if (!noGUIMode) setStatus(1,"");
}

/**
 * Function to be called on changing the
 * Work Config Directory textbox in the Folders tab.
 *
 */
void MainWindow::on_lineEdit_WorkFolder_editingFinished()
{
    QString valueString = ui->lineEdit_WorkFolder->text();
    string textHere = valueString.toStdString();
    std::replace(textHere.begin(), textHere.end(), '\\', '/');

    if ((textHere == "default") || (textHere == ""))
    {
        ui->lineEdit_WorkFolder->setPalette(*paletteGray);
        workFolder = "";
        visualization->setWorkFolder(textHere);
    }
    else
    {
        if (access(textHere.c_str(), 0) == 0)
        {
            ui->lineEdit_WorkFolder->setPalette(*paletteB);
            workFolder = textHere;
            visualization->setWorkFolder(textHere);
        }
        else
        {
            ui->lineEdit_WorkFolder->setPalette(*paletteR);
        }
    }
}

/**
 * Function to be called on clicking the
 * material codes button in the parameters tab.
 *
 */
void MainWindow::on_pushButton_6_clicked()
{
    QString messageString = "Available Materials like <br>  7 = Salt Water <br> 8 = Snow (0.03 g/cm<sup>3</sup>) <br> 9 = Water <br> 10 = Dry Air <br> 11 = Air with Humidity <br> 12 = Quarz <br> 20 = Soil with Soil Moisture <br> 21 = Plants <br> 23 = Cat Litter <br> 24 = Asphalt <br> 25 = HDPE  <br> 26 = Aluminum <br> 27 = Helium-3 <br> 28 = BF<sub>3</sub> <br> 29 = Gd<sub>2</sub>O<sub>3</sub><br>  31 = HDPE   <br> 32 = Stainless Steel     <br>  33 = Methane  <br> <br> See full list in the Manual or in the file MaterialCodes.txt";
    QMessageBox::about(this, tr("Material List"), messageString);
}

/**
 * Function to be called on clicking the
 * radio button fission source.
 *
 */
void MainWindow::on_radioButton_fission_clicked()
{
    doFusion = false;
    doFission = true;
    doAmBe = false;
    doMonoEnergetic = false;
    doThermalSource = false;
    doNoSource = false;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button thermal source.
 *
 */
void MainWindow::on_radioButton_ThermalSource_clicked()
{
    doFusion = false;
    doFission = false;
    doAmBe = false;
    doMonoEnergetic = false;
    doThermalSource = true;
    doNoSource = false;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button monoenergetic source.
 *
 */
void MainWindow::on_radioButton_MonoenergeticSource_clicked()
{
    doFusion = false;
    doFission = false;
    doAmBe = false;
    doMonoEnergetic = true;
    doThermalSource = false;
    doNoSource = false;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button fusion source.
 *
 */
void MainWindow::on_radioButton_fusion_clicked()
{
    doFusion = true;
    doFission = false;
    doAmBe = false;
    doMonoEnergetic = false;
    doThermalSource = false;
    doNoSource = false;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button no source.
 *
 */
void MainWindow::on_radioButton_NoSource_clicked()
{
    doFusion = false;
    doFission = false;
    doAmBe = false;
    doMonoEnergetic = false;
    doThermalSource = false;
    doNoSource = true;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button AmBe source.
 *
 */
void MainWindow::on_radioButton_AmBe_clicked()
{
    doFusion = false;
    doFission = false;
    doAmBe = true;
    doMonoEnergetic = false;
    doThermalSource = false;
    doNoSource = false;
    doModeratedCf = false;
}

/**
 * Function to be called on clicking the
 * radio button Moderated Cf source.
 *
 */
void MainWindow::on_radioButton_ModeratedCf_clicked()
{
    doFusion = false;
    doFission = false;
    doAmBe = false;
    doMonoEnergetic = false;
    doThermalSource = false;
    doNoSource = false;
    doModeratedCf = true;
}


/**
 * Function to be called to recreate the input matrices.
 *
 */
void recreateInputMatrices()
{
    inputPicVector.clear();

    /*
    inputMatrix1 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix1);
    inputMatrix2 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix2);
    inputMatrix3 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix3);
    inputMatrix4 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix4);
    inputMatrix5 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix5);
    inputMatrix6 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix6);
    inputMatrix7 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix7);
    inputMatrix8 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix8);
    inputMatrix9 =  TMatrixF(inputMatrixPixels,inputMatrixPixels); inputPicVector.push_back(inputMatrix9);
    inputMatrix10 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix10);
    inputMatrix11 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix11);
    inputMatrix12 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix12);
    inputMatrix13 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix13);
    inputMatrix14 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix14);
    inputMatrix15 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix15);
    inputMatrix16 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix16);
    inputMatrix17 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix17);
    inputMatrix18 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix18);
    inputMatrix19 =  TMatrixF(inputMatrixPixels,inputMatrixPixels);inputPicVector.push_back(inputMatrix19);
*/
}

/**
 * Function to be called to get vectors for new materials.
 * @param fileString
 */
vector<float> getMaterialVector(string fileString)
{
    vector<float> newMaterial;
    newMaterial.push_back(0);

    TString line;
    string temp;
    float number;

    ifstream input_stream(fileString, ios::in);
    while (line.ReadLine(input_stream))
    {
        istrstream stream(line.Data());
        stream >> temp;
        stream >> number;

        if ((number > 0) && (number < 100))
        {
            newMaterial.push_back(number);
        }
        else
        {
            newMaterial.push_back(0);
        }
    }
    return newMaterial;
}

/**
 * Function to be called to Check the input material definitions.
 *
 */
void checkInputMaterialDefinitions()
{
    ifstream* stream_in;
    bool foundData = false;

    string inputMaterialFileName;
    for (int i = 1; i < maxLayersAllowed; i++)
    {
        inputMaterialFileName = (string)workFolder+"Material"+(string)castIntToString(i)+".dat";
        stream_in = new ifstream(inputMaterialFileName, ifstream::in);
        if (stream_in->good())
        {
            inputMaterials[i - 1] = 1;
        }
        else inputMaterials[i - 1] = 0;
        stream_in->close();
    }

    for (int i = 1; i < maxLayersAllowed; i++)
    {
        if (inputMaterials[i - 1] == 1)
        {
            materialVector.push_back(getMaterialVector(inputMaterialFileName));
            materialVector.end()->at(0) = i;
        }
    }
}

/**
 * Function to be called to check the files possibly containing input pics.
 * inputpics: material definitions
 * inputpics2: density definitions
 * inputpics3: porosity definitions
 *
 */
void MainWindow::checkInputPics()
{
    if (!noGUIMode) {setStatus(1, "");   setStatus(2, "");}
    ifstream* stream_in;
    bool foundData = false;
    TString add = "";

    for (int i = 0; i < maxLayersAllowed; i++) { inputPicSizes[i] = 0; }

    for (int i = 1; i < maxLayersAllowed; i++)
    {
        stream_in = new ifstream((string)workFolder+castIntToString(i)+".dat",ifstream::in);
        if (stream_in->good())
        {
            inputPics[i - 1] = 1; // this is an ASCII material definition map
        }
        else inputPics[i - 1] = 0; // no map
        stream_in->close();

        stream_in = new ifstream((string)workFolder+castIntToString(i)+"d.dat",ifstream::in);
        if (stream_in->good())
        {
            inputPics2[i - 1] = 1; // this is an ASCII density definition map
        }
        else inputPics2[i - 1] = 0;
        stream_in->close();

        stream_in = new ifstream((string)workFolder+castIntToString(i)+"p.dat",ifstream::in);
        if (stream_in->good())
        {
            inputPics3[i - 1] = 1; // this is an ASCII porosity definition map
        }
        else inputPics3[i - 1] = 0;
        stream_in->close();

        QImage matrixImageHere;

        string imageStrHere = (string)workFolder+(string)castIntToString(i)+".png";

        if (matrixImageHere.load(QString::fromStdString(imageStrHere)))
        {
            inputPics[i - 1] = 2; // this is a png material definition map
        }

        imageStrHere = (string)workFolder+(string)castIntToString(i)+"d.png";

        if (matrixImageHere.load(QString::fromStdString(imageStrHere)))
        {
            inputPics2[i - 1] = 2; // this is a png density definition map
        }

        imageStrHere = (string)workFolder+(string)castIntToString(i)+"p.png";

        if (matrixImageHere.load(QString::fromStdString(imageStrHere)))
        {
            inputPics3[i - 1] = 2; // this is a png porosity definition map
        }

        //matrixImageHere;
    }
    int temp;

    TString line;
    int lineC = 0;
    QImage matrixImage;
    string imageStr;
    int matrixWidth;

    for (int i = 0; i < maxLayersAllowed; i++)
    {
        if ((inputPics[i] == 0) && (inputPics2[i] == 0) && (inputPics3[i] == 0))  inputPicSizes[i] = 0;
        else
        {
            temp = i + 1;

            if (inputPics[i] == 1) add = "";
            if (inputPics[i] == 2) add = "";

            if (inputPics[i] == 1)
            {
                lineC = 0;
                ifstream input_stream((string)workFolder+(string)castIntToString(temp)+add+".dat",ios::in);
                while (line.ReadLine(input_stream))
                {
                    istrstream stream(line.Data());
                    lineC++;
                }
                inputPicSizes[i] = lineC;
                //break;
                foundData = true;
            }

            if (inputPics[i] == 2)
            {
                imageStr = (string)workFolder+(string)castIntToString(temp)+add+".png";

                if (matrixImage.load(QString::fromStdString(imageStr)))
                {
                    matrixWidth = matrixImage.width();
                    inputPicSizes[i] = matrixWidth;
                    if (matrixImage.width() != matrixImage.height())
                    {
                        Int_t matrixWidthI = matrixWidth;  Int_t matrixHeightI = matrixImage.height();
                        string imageErrorString = (string)castIntToString(matrixWidthI)+"px  x "+(string)castIntToString(matrixHeightI)+"px";
                        setStatus(1, "Image Error: Asymmetry");  setStatus(2, imageErrorString);
                    }
                    else  foundData = true;
                }
            }

            if (inputPics2[i] == 1) add = "d";
            if (inputPics2[i] == 2) add = "d";

            if (inputPics2[i] == 1)
            {
                lineC = 0;
                ifstream input_stream((string)workFolder+(string)castIntToString(temp)+add+".dat",ios::in);
                while (line.ReadLine(input_stream))
                {
                    istrstream stream(line.Data());
                    lineC++;
                }

                if (inputPicSizes[i] < 1) { inputPicSizes[i] = lineC; foundData = true; }
                else if (inputPicSizes[i] != lineC)
                {
                    setStatus(1, "Image Error: Incongruency");   setStatus(2, "Density and material map size");
                    foundData = false;
                }
                else  foundData = true;
            }

            if (inputPics2[i] == 2)
            {
                imageStr = (string)workFolder+(string)castIntToString(temp)+add+".png";

                if (matrixImage.load(QString::fromStdString(imageStr)))
                {
                    matrixWidth = matrixImage.width();
                    if (matrixImage.width() != matrixImage.height())
                    {
                        Int_t matrixWidthI = matrixWidth;  Int_t matrixHeightI = matrixImage.height();
                        string imageErrorString = (string)castIntToString(matrixWidthI)+"px  x "+(string)castIntToString(matrixHeightI)+"px";
                        setStatus(1, "Image Error: Asymmetry");  setStatus(2, imageErrorString);
                    }
                    else
                        if (inputPicSizes[i] < 1) { inputPicSizes[i] = matrixWidth; foundData = true; }
                        else if (inputPicSizes[i] != matrixWidth)
                        {
                            setStatus(1, "Image Error: Incongruency");   setStatus(2, "Density and material map size");
                            foundData = false;
                        }
                        else  foundData = true;
                }
            }

            if (inputPics3[i] == 1) add = "p";
            if (inputPics3[i] == 2) add = "p";

            if (inputPics3[i] == 1)
            {
                lineC = 0;
                ifstream input_stream((string)workFolder+(string)castIntToString(temp)+add+".dat",ios::in);
                while (line.ReadLine(input_stream))
                {
                    istrstream stream(line.Data());
                    lineC++;
                }

                if (inputPicSizes[i] < 1) { inputPicSizes[i] = lineC; foundData = true; }
                else if (inputPicSizes[i] != lineC)
                {
                    setStatus(1, "Image Error: Incongruency");   setStatus(2, "Porosity and material map size");
                    foundData = false;
                }
                else  foundData = true;
            }

            if (inputPics3[i] == 2)
            {
                imageStr = (string)workFolder+(string)castIntToString(temp)+add+".png";

                if (matrixImage.load(QString::fromStdString(imageStr)))
                {
                    matrixWidth = matrixImage.width();
                    if (matrixImage.width() != matrixImage.height())
                    {
                        Int_t matrixWidthI = matrixWidth;  Int_t matrixHeightI = matrixImage.height();
                        string imageErrorString = (string)castIntToString(matrixWidthI)+"px  x "+(string)castIntToString(matrixHeightI)+"px";
                        setStatus(1, "Image Error: Asymmetry");  setStatus(2, imageErrorString);
                    }
                    else
                        if (inputPicSizes[i] < 1) { inputPicSizes[i] = matrixWidth; foundData = true; }
                        else if (inputPicSizes[i] != matrixWidth)
                        {
                            setStatus(1, "Image Error: Incongruency");   setStatus(2, "Porosity and material map size");
                            foundData = false;
                        }
                        else  foundData = true;
                }
            }
        }
    }

    string picLetter = "";
    inputMatrixDefsShown = false;
    if (inputMatrixDefsShown) cout<<"Input Matrix definitions: ";

    for (int z = 0; z < model->rowCount(); ++z)
    {
        QStandardItem* item = model->itemFromIndex(model->index(z, 3, QModelIndex()));
        item->setTextAlignment(Qt::AlignCenter);
        picLetter = "";

        temp = z + 1;
        if ((inputPics[z] == 0) && (inputPics2[z] == 0) && (inputPics3[z] == 0))
        {
            item->setText(QString::fromStdString((string)""));
        }
        else
        {
            if (inputMatrixDefsShown) cout<<castIntToString(temp);

            if (inputPics[z] == 1)
            {
                picLetter += "m";
                if (inputMatrixDefsShown) cout<<"m";
            }
            if (inputPics2[z] == 1)
            {
                picLetter += "d";
                if (inputMatrixDefsShown) cout<<"d";
            }
            if (inputPics3[z] == 1)
            {
                picLetter += "p";
                if (inputMatrixDefsShown) cout<<"p";
            }
            if (inputPics[z] == 2)
            {
                picLetter += "M";
                if (inputMatrixDefsShown) cout<<"M";
            }
            if (inputPics2[z] == 2)
            {
                picLetter += "D";
                if (inputMatrixDefsShown) cout<<"D";
            }
            if (inputPics3[z] == 2)
            {
                picLetter += "P";
                if (inputMatrixDefsShown) cout<<"P";
            }
            if (inputMatrixDefsShown) cout<<" ";

            item->setText(QString::fromStdString((string)castIntToString(temp) + picLetter + " [" + (string)castIntToString(inputPicSizes[z]) + "]"));
        }
        inputMatrixPixels = inputPicSizes[z]; //doesn't make so much sense
    }
    if (inputMatrixDefsShown) cout<<endl;

    if (foundData)
    {
        haveDifferentSoilMoistures = true;

        matrixStartX = squareDim / 1000. * (-0.5);
        matrixStartY = squareDim / 1000. * (-0.5);
        matrixMetricFactor = squareDim / 1000. / (inputMatrixPixels * 1.); //meter per pixel
    }
    else haveDifferentSoilMoistures = false;
}

/**
 * Function to be called on clicking Use layer maps checkbox in the Parameters tab.
 *
 */
void MainWindow::on_checkBox_useImage_clicked()
{
    bool checkBoxChecked = ui->checkBox_useImage->isChecked();

    if (!(checkBoxChecked))
    {
        useImage = false;
        haveDifferentSoilMoistures = false;
        layerMapsImport = false;
        for (int z = 0; z < model->rowCount(); ++z)
        {
            QStandardItem* item = model->itemFromIndex(model->index(z, 3, QModelIndex()));
            item->setText(QString::fromStdString((string)""));
        }
    }
    else
    {
        layerMapsImport = true;
        useImage = true;
        checkInputPics();
    }
}

/**
 * Function to be called on clicking view layer maps button in the Parameters tab.
 *
 */
void MainWindow::on_pushButton_Show_clicked()
{
    if ((useImage) && (haveDifferentSoilMoistures))
    {
        DialogShowPic* dialogshowpic = new DialogShowPic(this);
        dialogshowpic->show();
    }
}

/**
 * Function to be called on get the T matrix for the input pic definitions.
 *
 */
TMatrixF MainWindow::getTMatrix(int i)
{
    TMatrixF dummyTMatrix(0, 0);

    TString add = "";
    if (inputPics2[i] == 1) add = "d";
    if (inputPics2[i] == 2) add = "d";
    if (inputPics3[i] == 1) add = "p";
    if (inputPics3[i] == 2) add = "p";

    int i2 = i + 1;
    TString str = castIntToString(i2) + add;

    if ((inputPics[i] == 1) || ((inputPics2[i] == 1)) || (inputPics3[i] == 1))
    {
        TMatrixF matr = readmatrix(workFolder, str, "dat", -1, inputPicSizes[i]);
        turnInputMatrix(matr);
        return matr;
    }

    if ((inputPics[i] == 2) || ((inputPics2[i] == 2)) || (inputPics3[i] == 2))
    {
        QRgb pixelValue;
        string imageStr = (string)workFolder+(string)str+".png";
        QImage matrixImage;
        matrixImage.load(QString::fromStdString(imageStr));
        int matrixWidth = matrixImage.width();
        int matrixHeight = matrixImage.height();

        TMatrixF matr(matrixWidth, matrixHeight);
        for (int i = 0; i < matrixWidth; ++i)
        {
            for (int j = 0; j < matrixHeight; ++j)
            {
                pixelValue = matrixImage.pixel(i, matrixHeight - 1 - j);
                matr(i, j) = qGray(pixelValue);
            }
        }
        return matr;
    }

    return  dummyTMatrix;
}

/**
 * Function to be called on changing the upper bound slider in the Parameters tab.
 *
 */
void MainWindow::on_horizontalSliderColor_sliderMoved(int position)
{
    if (alreadyStarted)
    {
        silderColorMoved = true;

        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on changing the lower bound slider in the Parameters tab.
 *  @param position
 */
void MainWindow::on_horizontalSliderColorZero_sliderMoved(int position)
{
    if (alreadyStarted)
    {
        silderColorMoved = true;

        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the radio button Dark Gray Scale
 * for the Neutron Color Scheme in the Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronNight_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the radio button Cold
 * for the Neutron Color Scheme in the Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronCold_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the radio button Polar
 * for the Neutron Color Scheme in the Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronPolar_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the radio button URANOS
 * for the Neutron Color Scheme in the Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronRainbow_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on moving the slider for the footprint moisture.
 * deprecated
 *
 */
void MainWindow::on_horizontalSliderFPMoist_sliderMoved(int position)
{
    fpSoilMoist = (position * 1.) / 200.;
    float fpSoilMoist2 = fpSoilMoist * 100.;
    string text = castFloatToString(fpSoilMoist2, 4) + " %";
    ui->labelFPMoist->setText(QString::fromStdString(text));

    rangeIntegral = -1;

    replotFootprint();
}

/**
 * Function to be called on moving the slider for the footprint humidity.
 * deprecated
 *
 */
void MainWindow::on_horizontalSliderFPHum_sliderMoved(int position)
{
    fpHum = (position * 1.);
    int fpHumInt = fpHum;
    string text = castIntToString(fpHumInt) + " g/cm3";
    ui->labelFPHum->setText(QString::fromStdString(text));

    rangeIntegral = -1;

    replotFootprint();
}

/**
 * Function to be called on checking the checkbox for display the footprint function
 * deprecated
 *
 */
void MainWindow::on_checkBoxFPLog_clicked()
{
    if (activateFP) replotFootprint();
}

/**
 * Function to be called on moving the slider for the footprint moisture.
 * deprecated
 *
 */
void MainWindow::on_horizontalSliderFPHum_2_sliderMoved(int position)
{
    if (!activateFP) return;

    plotGUIxBinsFootprFuncLine[0] = 0;
    plotGUIxBinsFootprFuncLine[1] = position;
    plotGUIxBinsFootprFuncLine[2] = position + 1;
    plotGUIxBinsFootprFuncLine[3] = position + 2;

    integralSliderMoved = true;

    string text = (string)castIntToString(position)+" m";
    ui->label_FPIntegral->setText(QString::fromStdString(text));

    replotFootprint();
}

void MainWindow::on_pushButton_ActivateFP_clicked()
{
    if (activateFP)
    {
        ui->pushButton_ActivateFP->setChecked(false); ui->pushButton_ActivateFP->setDown(false);
        activateFP = false;
    }
    else
    {
        ui->pushButton_ActivateFP->setChecked(true); ui->pushButton_ActivateFP->setDown(true);
        activateFP = true;
        replotFootprint();
    }
}


/**
 * Function to be called on changing the value of Only Record in Material No
 * textbox in the detector layer of detector tab.
 *
 */
void MainWindow::on_lineEditScotoma_editingFinished()
{
    QString arg = ui->lineEditScotoma->text();
    float angleString = arg.toFloat();
    ui->lineEditScotoma->setPalette(*paletteB);
    if ((angleString >= 0) && (angleString <= 360) && (angleString >= 0))   downwardScotomaAngle = (angleString / 2.) / 360. * 2. * TMath::Pi(); // (180-angleString/2.)/360.*2.*TMath::Pi();
    else ui->lineEditScotoma->setPalette(*paletteR);
}

/**
 * Function to be called on  clicking the Gradient view
 * checkbox in the displayed energy window of the display tab.
 *
 */
void MainWindow::on_checkBoxGradient2_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the check box
 * Selected Energy Data in spatial distributions in the export tab.
 *
 */
void MainWindow::on_checkBoxSelectedData_clicked()
{
    if (exportSelectedData) exportSelectedData = false;
    else exportSelectedData = true;
}

/**
 * Function to be called on clicking the check box
 * Fast Neutron Data in spatial distributions in the export tab.
 *
 */
void MainWindow::on_checkBoxFastData_clicked()
{
    if (exportFastData) exportFastData = false;
    else exportFastData = true;
}

/**
 * Function to be called on clicking the check box
 * Intermediate Energy Data in spatial distributions in the export tab.
 *
 */
void MainWindow::on_checkBoxIntermediateData_clicked()
{
    if (exportIntermediateData) exportIntermediateData = false;
    else exportIntermediateData = true;
}

/**
 * Function to be called on clicking the check box
 * Epithermal Data in spatial distributions in the export tab.
 *
 */
void MainWindow::on_checkBoxEpithermalData_clicked()
{
    if (exportEpithermalData) exportEpithermalData = false;
    else exportEpithermalData = true;
}


/**
 * Function to be called on clicking the check box
 * Thermal Energy Data in spatial distributions in the export tab.
 *
 */
void MainWindow::on_checkBoxThermalData_clicked()
{
    if (exportThermalData) exportThermalData = false;
    else exportThermalData = true;
}

/**
 * Function to be called on changing the slider Soil Porosity [Vol%]
 * in the parameters tab.
 *  @param position
 */
void MainWindow::on_sliderSoilPorosity_sliderMoved(int position)
{
    soilSolidFracVar = 1 - (position * 1.) / 100.;

    string text = castFloatToString((position * 1.), 4) + " %";
    ui->labelSP->setText(QString::fromStdString(text));
}

/**
 * Function to be called on clicking the Advanced Analysis Raw Output (ROOT)
 * checkbox in General Options in the export tab.
 * @param checked
 */
void MainWindow::on_checkBoxFileOutputPDF_2_clicked(bool checked)
{
    if (!checked) uranosRootOutput = false;
    else uranosRootOutput = true;
}

/**
 * Function to be called on clicking the Create new folder for every export
 * checkbox in General Options in the export tab.
 * @param checked
 */
void MainWindow::on_checkBoxCreateFolder_clicked(bool checked)
{
    if (!checked) createSeparateFolderEachExport = false;
    else createSeparateFolderEachExport = true;
}

/**
 * Function to be called on clicking the Create new folder for every export
 * checkbox in General Options in the export tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_xPos_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_xPos->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  xPosSource = valueString * 1000.;
    else ui->lineEdit_xPos->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_yPos_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_yPos->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  yPosSource = valueString * 1000.;
    else ui->lineEdit_yPos->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_zPos_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_zPos->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  zPosSource = valueString * 1000.;
    else ui->lineEdit_zPos->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_xSize_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_xSize->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  xSizeSource = valueString * 1000.;
    else ui->lineEdit_xSize->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_ySize_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_ySize->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  ySizeSource = valueString * 1000.;
    else ui->lineEdit_ySize->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_zPos_2_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_zPos_2->setPalette(*paletteB);
    if ((valueString > -0.01) && (valueString < 1000))  radiusSource = valueString * 1000.;
    else ui->lineEdit_zPos_2->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_SourceEnergy_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_SourceEnergy->setPalette(*paletteB);
    if ((valueString > 0) && (valueString < 100) && (true))  sourceEnergy = valueString;
    else ui->lineEdit_SourceEnergy->setPalette(*paletteR);
}

/**
 * Function to be called on clicking the activate thermal transport checkbox
 * checkbox in the Computational parameters tab.
 * @param checked
 */
void MainWindow::on_checkBoxThermal_toggled(bool checked)
{
    if (!checked) noThermalRegime = false;
    else noThermalRegime = true;
}

/**
 * Function to be called on changing the Lower bound
 * textbox value in Neutron Color Scheme in Display tab.
 * @param arg1
 */
void MainWindow::on_lineEditManualColorZero_textChanged(const QString& arg1)
{
    int valueString = arg1.toInt();
    ui->lineEditManualColorZero->setPalette(*paletteB);
    if ((valueString > -100000) && (valueString < 100000000))  manualColorZero = valueString;
    else ui->lineEditManualColorZero->setPalette(*paletteR);

    redrawTopView();
}

/**
 * Function to be called on changing the Upper bound
 * textbox value in Neutron Color Scheme in Display tab.
 * @param arg1
 */
void MainWindow::on_lineEditManualColor_textChanged(const QString& arg1)
{
    int valueString = arg1.toInt();
    ui->lineEditManualColor->setPalette(*paletteB);
    if ((valueString > -1000000) && (valueString < 100000000))  {manualColor = valueString; if (manualColor <= manualColorZero) manualColor = manualColorZero + 1;}
    else ui->lineEditManualColor->setPalette(*paletteR);

    redrawTopView();
}

/**
 * Function to be called on clicking the No Track Recording checkbox in
 * Disable Tracks and only count Layer hits in Computational Speed in the Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxNoTrack_toggled(bool checked)
{
    if (!checked) noTrackRecording = false;
    else noTrackRecording = true;
}

/**
 * Function to be called on clicking the check every checkbox in
 * Special Purpose in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxSaveEvery_toggled(bool checked)
{
    if (!checked) saveEveryXNeutrons = false;
    else saveEveryXNeutrons = true;

    saveEveryXNeutronsPower = ui->spinBoxSaveEvery->value();
}

/**
 * Function to be called on changing the value of Only Record in Material No in
 * Detector Layer in Detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditScotoma_2_textChanged(const QString& arg1)
{
    int valueString = arg1.toInt();
    ui->lineEditScotoma_2->setPalette(*paletteB);
    if ((valueString > -1) && (valueString <= 256))  detectorSensitiveMaterial = valueString;
    else ui->lineEditScotoma_2->setPalette(*paletteR);
}

/**
 * Function to be called on changing the value of Only Record in Material No in
 * Detector Layer in Detector tab.
 *
 */
void MainWindow::on_checkBox_clicked()
{
    if (useDetectorSensitiveMaterial) useDetectorSensitiveMaterial = false;
    else useDetectorSensitiveMaterial = true;
}

/**
 * Function to be called on clicking Record Neutrons only once per Layer checkboxin
 * Exclude Multiple Scatteringin Detector Layer in Detector tab.
 * @param checked
 */
void MainWindow::on_checkBox_NoMultipleScattering_toggled(bool checked)
{
    if (!checked) noMultipleScatteringRecording = false;
    else noMultipleScatteringRecording = true;

    noMultipleScatteringRecording = false;
}

/**
 * Function to be called on clicking Manual checkbox
 * in Neutron Color Scheme in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxManual_toggled(bool checked)
{
    redrawTopView();
    //redrawNeutronMap(0);
}

/**
 * Function to be called on clicking Track all Layers checkbox
 * in Detector Layer in Detector tab.
 * @param checked
 */
void MainWindow::on_checkBox_TrackAllLayers_toggled(bool checked)
{
    if (checked) trackAllLayers = true;
    else trackAllLayers = false;
}

/**
 * Function to be called on clicking Logarthmic checkbox
 * in Neutron Color Scheme in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxLogarithmic_toggled(bool checked)
{
    if (checked) {plotTopViewLog = true; visualization->setplotTopViewLog(true); }
    else { plotTopViewLog = false; visualization->setplotTopViewLog(false); }

    if (checked) {visualization2->setplotTopViewLog2(true); }
    else {visualization2->setplotTopViewLog2(false); }

    redrawTopView();
}

/**
 * Function to be called on clicking Gray Scale radio button
 * in Neutron Color Scheme in Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronGrayScale_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking Hot radio button
 * in Neutron Color Scheme in Display tab.
 *
 */
void MainWindow::on_radioButton_NeutronHot_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking Thermal radio button
 * in Neutron Color Scheme in Display tab .
 *
 */
void MainWindow::on_radioButton_NeutronThermal_clicked()
{
    if (alreadyStarted)
    {
        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking Vertical Cylinder radio button
 * in Detector Type in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_Cylinder_toggled(bool checked)
{
    if (checked)
    {
        useCylindricalDetector = true;
        useSphericalDetector = false;
        useySheetDetector = false;
        usexSheetDetector = false;
        ui->lineEditDetLength->setEnabled(false);
        ui->lineEditDetRad->setEnabled(true);
    }
    else
    {
        //useSphericalDetector = true;
        //useCylindricalDetector = false;
    }
}

/**
 * Function to be called on clicking Sphere radio button
 * in Detector Type in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_Sphere_toggled(bool checked)
{
    if (checked)
    {
        useCylindricalDetector = false;
        useSphericalDetector = true;
        useySheetDetector = false;
        usexSheetDetector = false;
        ui->lineEditDetLength->setEnabled(false);
        ui->lineEditDetRad->setEnabled(true);
    }
    else
    {
        //useSphericalDetector = false;
        //useCylindricalDetector = true;
    }
}

/**
 * Function to be called on clicking Energy Band Model radio button
 * in Precision in Detector Layer in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_detectorLayerEnergyBand_toggled(bool checked)
{
    if (checked)
    {
        useRealisticModelLayer = false;
    }
}

/**
 * Function to be called on clicking Physics Model radio button
 * in Precision in Detector Layer in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_detectorLayerRealistic_toggled(bool checked)
{
    if (checked)
    {
        useRealisticModelLayer = true;
    }
}

/**
 * Function to be called on clicking Energy Band Model radio button
 * in Precision in Detector in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_detectorEnergyBand_toggled(bool checked)
{
    if (checked)
    {
        useRealisticModelDetector = false;
    }
}

/**
 * Function to be called on clicking Physics Model radio button
 * in Precision in Detector in Detector tab.
 * @param checked
 */
void MainWindow::on_radioButton_detectorRealistic_toggled(bool checked)
{
    if (checked)
    {
        useRealisticModelDetector = true;
    }
}


/**
 * Function to be called on clicking Volume Source extended down to the ground radio button
 * in source geometry in computational parameters tab.
 * @param checked
 */
void MainWindow::on_checkBox_VolumeSource_toggled(bool checked)
{
    if (checked)
    {
        useVolumeSource = true;
    }
    else
    {
        useVolumeSource = false;
    }
}

/**
 * Function to be called on clicking Experimental High Energy Cascade Model check box
 * in Computational Model in computational parameters tab.
 * @param checked
 */
void MainWindow::on_checkBox_HEModel_toggled(bool checked)
{
    if (checked)
    {
        useHECascadeModel = true;
    }
    else
    {
        useHECascadeModel = false;
    }
}

/**
 * Function to be called on clicking Thermal Physics Extension Model check box
 * in Computational Model in computational parameters tab.
 * @param checked
 */
void MainWindow::on_checkBox_activateThermal_toggled(bool checked)
{
    if (checked)
    {
        noThermalRegime = false;
    }
    else
    {
        noThermalRegime = true;
    }
}

/**
 * Function to be called on changing Range of Interest slider
 * in Detector Color Scheme in Display tab.
 * @param position
 */
void MainWindow::on_horizontalSliderDetector_sliderMoved(int position)
{
    if (alreadyStarted)
    {
        silderDetectorMoved = true;

        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on changing Maximum Scale slider
 * in Detector Color Scheme in Display tab.
 * @param position
 */
void MainWindow::on_horizontalSliderDetectorColor_sliderMoved(int position)
{
    if (alreadyStarted)
    {
        silderDetectorColorMoved = true;

        scaleFactor = ui->horizontalSliderDetectorColor->value();

        if (simulationRunning)
        {
            if (pausehere)redrawTopView();
        }
        else redrawTopView();
    }
}

/**
 * Function to be called on clicking the Enlarge button (+)
 * in the Birds Eye view and Spectra tab.
 *
 */
void MainWindow::on_pushButton_Enlarge_clicked()
{
    delete visualization;
    visualization = new VisualizationEnlarge(this);
    setConnectionsEnlarge();
    visualization->show();

    updateEnlargedView = true;
    if (alreadyStarted)
    {
        redrawEnlargedView();
    }
}

/**
 * Function to be called on changing Downward Accpetance Angle [0-360°]
 * in the Field of View Limitation in Detector tab.
 *
 */
void MainWindow::on_lineEditAntiScotoma_editingFinished()
{
    QString arg = ui->lineEditAntiScotoma->text();
    float angleString = arg.toFloat();
    ui->lineEditAntiScotoma->setPalette(*paletteB);
    if ((angleString >= 0) && (angleString <= 360) && (angleString >= 0))   downwardAcceptanceAngle = (angleString / 2.) / 360. * 2. * TMath::Pi();
    else ui->lineEditAntiScotoma->setPalette(*paletteR);
}

/**
 * Function to be called on clicking Clear every checkbox
 * in the Special Purpose in Display tab.
 * @param - checked
 */
void MainWindow::on_checkBoxSaveEvery_2_clicked(bool checked)
{
    clearEveryXNeutrons = checked;
}

void MainWindow::on_lineEditClearEveryXNeutrons_editingFinished()
{

}

/**
 * Function to be called on changing Neutrons textbox
 * in the Special Purpose in Display tab.
 * @param arg1
 */
void MainWindow::on_lineEditClearEveryXNeutrons_textChanged(const QString& arg1)
{
    int valueString = arg1.toInt();
    ui->lineEditClearEveryXNeutrons->setPalette(*paletteB);
    if ((valueString > -1) && (valueString < 2e8))  clearEveryXNeutronsNumber = valueString;
    else ui->lineEditClearEveryXNeutrons->setPalette(*paletteR);
}

/**
 * Function to be called on changing Auto Update textbox
 * in the Computational Speed in Display tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_AutoUpdate_textChanged(const QString& arg1)
{
    float numberString = arg1.toFloat();
    ui->lineEdit_AutoUpdate->setPalette(*paletteB);
    if ((numberString > 0) && (numberString <= 360))   refreshTime = numberString;
    else ui->lineEdit_AutoUpdate->setPalette(*paletteR);
}

/**
 * Function to be called on changing Auto Refresh Ratecheckbox
 * in the Computational Speed in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxAutoRefreshRate_toggled(bool checked)
{
    setAutoRefreshRate = checked;
}

/**
 * Function to be called on toggling Auto Refresh checkbox
 * in the Computational Speed in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxClearEveryDisplayRefresh_toggled(bool checked)
{
    setAutoRefreshRateClearing = checked;
}

/**
 * Function to be called on clicking Auto Refresh checkbox
 * in the Computational Speed in Display tab.
 * @param checked
 */
void MainWindow::on_checkBoxClearEveryDisplayRefresh_clicked(bool checked)
{
    //setAutoRefreshRateClearing = checked;
}

/**
 * Function to be called on clicking Undefined Material Warnings
 * checkbox in General Options in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxWarnMAterial_toggled(bool checked)
{
    warnUndefinedMaterial = checked;
}

/**
 * Function to be called on clicking Track Data
 * in Spatial Distributions in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxTrackingData_clicked(bool checked)
{
    exportTrackData = checked;
}


/**
 * Function to be called on clicking High Res Track Data
 * in Spatial Distributions in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxHighResTrackingData_clicked(bool checked)
{
    exportHighResTrackData = checked;
    highResCalc = checked;
}

/**
 * Function to be called on clicking Sheet along y-Axis radio button
 * in Detector Type in Detector tab.
 *
 */
void MainWindow::on_radioButton_ySheet_clicked()
{
    useCylindricalDetector = false;
    useSphericalDetector = false;
    useySheetDetector = true;
    usexSheetDetector = false;
    ui->lineEditDetLength->setEnabled(true);
    ui->lineEditDetRad->setEnabled(false);
}

/**
 * Function to be called on clicking Sheet along x-Axis radio button
 * in Detector Type in Detector tab.
 *
 */
void MainWindow::on_radioButton_xSheet_clicked()
{
    useCylindricalDetector = false;
    useSphericalDetector = false;
    useySheetDetector = false;
    usexSheetDetector = true;
    ui->lineEditDetLength->setEnabled(true);
    ui->lineEditDetRad->setEnabled(false);
}

/**
 * Function to be called on clicking Sheet Length [m] textbox
 * in Detector in Detector tab.
 * @param arg1
 */
void MainWindow::on_lineEditDetLength_textChanged(const QString& arg1)
{
    float detLengthString = arg1.toFloat();
    ui->lineEditDetLength->setPalette(*paletteB);
    if ((detLengthString > -0.01) && (detLengthString < 1e4))   detLength = 1000. * detLengthString;
    else ui->lineEditDetLength->setPalette(*paletteR);
}

void MainWindow::on_lineEdit_zSize_textChanged(const QString& arg1)
{
    float valueString = arg1.toFloat();
    ui->lineEdit_zSize->setPalette(*paletteB);
    if ((valueString > -10000) && (valueString < 10000) && (valueString < 10. * squareDim / 1000.))  zSizeSource = valueString * 1000.;
    else ui->lineEdit_zSize->setPalette(*paletteR);
}

void MainWindow::on_radioButton_TopToBottom_toggled(bool checked)
{
    if (checked) sourceDirection = 1;
}

void MainWindow::on_radioButton_BottomToTop_toggled(bool checked)
{
    if (checked) sourceDirection = 2;
}

void MainWindow::on_radioButton_LeftToRight_toggled(bool checked)
{
    if (checked) sourceDirection = 3;
}

void MainWindow::on_radioButton_RightToLeft_toggled(bool checked)
{
    if (checked) sourceDirection = 4;
}

void MainWindow::on_radioButton_Omni_toggled(bool checked)
{
    if (checked) sourceDirection = 0;
}

void MainWindow::on_checkBox_DomainCutoff_toggled(bool checked)
{
    domainCutoff = checked;
}


/**
 * Function to be called on changing x domain size textbox in
 * Computaional Boundary Conditions in Computational parameters.
 * @param arg1
 */
void MainWindow::on_lineEditDomainFactor_textChanged(const QString& arg1)
{
    //int squareDimSring = arg1.toInt();
    double dimString = arg1.toDouble();
    ui->lineEditDomainFactor->setPalette(*paletteB);
    if ((dimString > 1) && (dimString < 1e9))
    {
        domainCutoffFactor = dimString;
    }
    else ui->lineEditDomainFactor->setPalette(*paletteR);
}

/**
 * Function to be called on clicking Remove if textbox for x domain size in
 * Computaional Boundary Conditions in Computational parameters.
 * @param checked
 */
void MainWindow::on_checkBox_DomainCutoffMeters_toggled(bool checked)
{
    domainCutoff2 = checked;
}

/**
 * Function to be called on changing m beyond domain textbox in
 * Computaional Boundary Conditions in Computational parameters.
 * @param arg1
 */
void MainWindow::on_lineEditDomainMeters_textChanged(const QString& arg1)
{
    int dimString = arg1.toInt();
    ui->lineEditDomainMeters->setPalette(*paletteB);
    if ((dimString > -1) && (dimString < 1e9))
    {
        domainCutoffMeters = dimString;
    }
    else ui->lineEditDomainMeters->setPalette(*paletteR);
}

/**
 * Function to be called on clicking Write Detector Neutron Tracks to File checkbox in
 * Individual Neutron Data in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxFileOutput2_toggled(bool checked)
{
    detTrackFileOutput = checked;
    trackAllLayers = checked;
    ui->checkBox_TrackAllLayers->setChecked(checked);
}

/**
 * Function to be called on clicking Export all Tracks checkbox in
 * Individual Neutron Data in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxExportAllTracks_toggled(bool checked)
{
    allTrackFileOutput = checked;
    trackAllLayers = checked;
    ui->checkBox_TrackAllLayers->setChecked(checked);
}

/**
 * Function to be called on clicking Export all Tracks checkbox in
 * Individual Neutron Data in Export tab.
 * @param checked
 */
void MainWindow::on_checkBox_ReflectiveBoundaries_toggled(bool checked)
{
    ui->checkBox_ReflectiveBoundaries->setChecked(checked);
    reflectiveBoundaries = checked;
}

/**
 * Function to be called on clicking Reflective Boundaries checkbox in
 * Computational Boundary Conditions in Computational Parameter tab.
 * @param checked
 */
void MainWindow::on_checkBox_PeriodicBoundaries_toggled(bool checked)
{
    ui->checkBox_PeriodicBoundaries->setChecked(checked);
    periodicBoundaries = checked;
}

/**
 * Function to be called on changing Detector Energy Calibration File textbox in
 * Folders tab.
 *
 */
void MainWindow::on_lineEdit_DetectorFile_editingFinished()
{
    QString valueString = ui->lineEdit_DetectorFile->text();
    string textHere = valueString.toStdString();

    std::replace(textHere.begin(), textHere.end(), '\\', '/');

    ifstream f(textHere.c_str());

    ui->lineEdit_DetectorFile->setPalette(*paletteB);

    if (!(f.good()))
    {
        if ((textHere == "default") || (textHere == "") || (textHere == "N/A") || (textHere == "n/a"))
        {
            ui->lineEdit_DetectorFile->setPalette(*paletteGray);
            detectorResponseFunctionFile = "";
        }
        else
        {
            ui->lineEdit_DetectorFile->setPalette(*paletteR);
        }
    }
    else
    {
        if (textHere.length() < 7)  ui->lineEdit_DetectorFile->setPalette(*paletteR);
        else
        {
            detectorResponseFunctionFile = textHere;
        }
    }
}

/**
 * Function to be called on changing Detector Energy Calibration File textbox in
 * Folders tab.
 * @param arg1
 */
void MainWindow::on_lineEdit_DetectorFile_textChanged(const QString& arg1)
{
    QString valueString = ui->lineEdit_DetectorFile->text();
    string textHere = valueString.toStdString() + ".2";

    std::replace(textHere.begin(), textHere.end(), '\\', '/');

    ifstream f(textHere.c_str());

    ui->lineEdit_DetectorFile->setPalette(*paletteB);

    if (!(f.good()))
    {
    }
    else
    {
        if (textHere.length() < 7) {}
        else
        {
            ui->label_25->setText("Detector Energy Calibration File (+1)");
            detectorResponseFunctionFile2 = textHere;
            useAdditionalDetectorModel = true;
        }
    }
}

/**
 * Function to be called on clicking Write Detector Neutron Hits to File checkbox in
 * Individual Neutron Data in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxFileOutput_toggled(bool checked)
{
    detFileOutput = checked;
}

/**
 * Function to be called on clicking Write Detector Layer Neutron Hits to File checkbox in
 * Individual Neutron Data in Export tab.
 * @param checked
 */
void MainWindow::on_checkBoxFileOutput3_toggled(bool checked)
{
    detLayerFileOutput = checked;
}

/**
 * Exports the settings of the GUI to the cfg file
 *
 */
void MainWindow::on_pushButton_SaveConfig_clicked()
{
    if (configFilePathConfigured)
    {
        exportSettings(configFilePath);
    }
    else
    {
        exportSettings("");
    }
    if (!noGUIMode) setStatus(1, "Uranos.cfg written");
}

/**
 * Opens the second enlarged visualization window
 *
 */
void MainWindow::on_pushButton_Enlarge2_clicked()
{
    delete visualization2;
    visualization2 = new VisualizationEnlarge2(this);
    visualization2->show();

    updateEnlargedView2 = true;
    if (alreadyStarted)
    {
        redrawEnlargedView2();
    }
}

/**
 * Activates the individual source placement on mouse event
 * @param checked
 */
void MainWindow::setGodzillaMode(bool var)
{
   godzillaMode = var;
}

/**
 * Activates the individual source placement on mouse event
 * @param xc, yc
 */
void MainWindow::setSourcePos(float xc, float yc)
{
   godzillaMode = true;
   xCustomPos = xc;
   yCustomPos = yc;
}

/**
 * Activates the individual source placement on mouse event
 * @param xc
 */
void MainWindow::setSourcePosX(float xc)
{
   //godzillaMode = true;
   if (godzillaMode) xCustomPos = xc;
}

/**
 * Activates the individual source placement on mouse event
 * @param yc
 */
void MainWindow::setSourcePosY(float yc)
{
   //godzillaMode = true;
   if (godzillaMode) yCustomPos = yc;
}

void MainWindow::setPause(bool var)
{
    pausehere = var;
}

void MainWindow::setContinue(bool var)
{
    pausehere = !var;
}

void MainWindow::setStopMode(bool var)
{
    stopMode = var;
}

