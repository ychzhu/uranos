/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>

#include "qcustomplot.h"

#include "Toolkit.h"

#include "TRandom3.h"
#include <string>
#include <vector>
#include "time.h"
#include "TLegend.h"

#include "dialogshowpic.h"
#include "visualizationenlarge.h"
#include "visualizationenlarge2.h"


namespace Ui {
class MainWindow;

}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    QMessageBox* msgBox;
public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


    void setupImport();
    void setupGraph(int index);
    void getGraphConfig();
    void setupRunSpectraGraph(QCustomPlot *customPlot);
    void setupRunBirdsEyeViewGraph(QCustomPlot *customPlot);
    void setupRunHorizSliceXGraph(QCustomPlot *customPlot);
    void setupRunNpassDetGraph(QCustomPlot *customPlot);
    void setupRunDepOfIntGraph(QCustomPlot *customPlot);
    void setupRunNpassDetLyrGraph(QCustomPlot *customPlot);
    void setupRunRelIntVsEGraph(QCustomPlot *customPlot);
    void setupRunOrginOfNGraph(QCustomPlot *customPlot);

    void setupRangeFunction(QCustomPlot *customPlot);

    void setupLiveTHs();

    void setupTable(QTableView* table);

    void redrawEnlargedView();

    void redrawEnlargedView2();

    void redrawNeutronMap(double difftime);

    void buttonClickFunction();

    void activateThermal();

    void setConfigFilePath(string pathtoConfigFile);

    void disabledGUIRun(string pathtoConfigFile);

    void activateSkyEvaporation();

    void activateThermalSkyEvaporation();

    void activateDetectorBatchRun();

    void activateDetectorBatchRun2();

    void activateDetectorAngleBatchRun();

    void activateBatchRun();

    void activateParameterBatchRun(int parMin, int parMax);

    void beSilent();

    void setGeometry();

    void setStatus(int numberLabel, string msg);

    void setupGeometry();

    void exportToSave();

   //static TMatrixF readMatrixPNG(TString folder, TString filename);

    static TMatrixF getTMatrix(int i);

    //bool cosmicNSimulator(MainWindow* uiM);

    void redrawTopView();

    void redrawSideView();

    void formatPlotToColor(QCustomPlot *customPlot);

    void replotFootprint();

    void formatForVectorGraphics();

    void checkInputPics();

    TSpline3* getSplinedEnergyModelFromMatrix(TMatrixF* matrix, bool logXValues, bool normalize);

    TSpline3* getSplinedDetectorEnergyModelFromFile(TString dname, int linesToSkip, bool logXValues, bool normalize);

    bool loadParamData(string fileFolder);

    bool isFinished() const { return isFinishedVar; }

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_spinBox_2_valueChanged(int arg1);

    void on_pushButton_Simulate_clicked();

    void on_pushButton_stop_clicked();

    void on_checkBoxR_stateChanged(int arg1);

    void on_checkBoxD_stateChanged(int arg1);

    void on_sliderSoilMoisture_sliderMoved(int position);

    void on_sliderAirHum_sliderMoved(int position);

    void on_lineEditNeutrinsTotal_textChanged(const QString &arg1);

    void on_sliderAtm1_sliderMoved(int position);

    void on_sliderRigidity_sliderMoved(int position);

    void on_lineEditSquareDim_textChanged(const QString &arg1);

    void on_lineEditRefresh_textChanged(const QString &arg1);

    //void on_pushButtonClear_clicked();

    void on_pushButton_Clear_clicked();

    void on_pushButton_Pause_clicked();

    void on_lineEditDetRad_textChanged(const QString &arg1);

    void on_lineEditDetX_textChanged(const QString &arg1);

    void on_lineEditDety_textChanged(const QString &arg1);

    void on_beamRound_clicked();

    void on_beamSquare_clicked();

    void on_radioRiver_clicked();

    void on_radioCoast_clicked();

    void on_radioIsland_clicked();

    void on_radioLake_clicked();

    void on_lineEditBeamRad_textChanged(const QString &arg1);

    void on_radioRiver_toggled(bool checked);

    void on_radioButton_clicked();

    void on_lineEditTHLhigh_textChanged(const QString &arg1);

    void on_lineEditTHLlow_textChanged(const QString &arg1);

    void on_checkBoxRS_clicked();

    void on_lineEdit_River_textChanged(const QString &arg1);

    void on_lineEdit_River_2_textChanged(const QString &arg1);

    void on_lineEdit_Island_textChanged(const QString &arg1);

    void on_lineEdit_Lake_selectionChanged();

    void on_lineEdit_Lake_textChanged(const QString &arg1);

    void on_checkBoxFileOutput_clicked();

    //void on_radioButton_2_clicked();

    void on_radioButton_map_clicked();

    void on_radioButton_mapInter_clicked();

    void on_radioButton_mapFast_clicked();

    void on_radioButton_mapAlbedo_clicked();


    void on_lineEdit_InputSpectrumFolder_editingFinished();

    void on_lineEdit_CrosssectionFolder_editingFinished();

    void on_lineEdit_OutputFolder_editingFinished();

    void on_pushButton_about_clicked();

    void on_checkBoxTransparent_clicked();

    void on_checkBoxBasicSpectrum_clicked();

    void on_pushButton_AddLayer_clicked();

    void setIntType1();
    void setIntType2();
    void setIntType3();

    bool eventFilter(QObject *target, QEvent *event);

    void setFocus(const QModelIndex &idx);

    QPushButton* typebutton(QModelIndex idx);

    void on_pushButton_RemoveLayer_clicked();

    void on_pushButton_ReadGeometry_clicked();

    void on_pushButton_SaveGeometry_clicked();

    void on_spinBox_StartingLayer_valueChanged(const QString &arg1);

    void on_spinBox_StartingLayer_valueChanged(int arg1);

    void on_spinBox_DetectorLayer_valueChanged(int arg1);

    void on_spinBox_GroundLayer_valueChanged(int arg1);

    void on_pushButton_LoadGeometry_clicked();

    void on_lineEdit_WorkFolder_editingFinished();

    void on_pushButton_6_clicked();

    void on_radioButton_fission_clicked();

    void on_radioButton_fusion_clicked();

    void on_checkBox_useImage_clicked();

    void on_pushButton_Show_clicked();

    void on_horizontalSliderColor_sliderMoved(int position);

    void on_horizontalSliderColorZero_sliderMoved(int position);

    void on_radioButton_NeutronNight_clicked();

    void on_radioButton_NeutronCold_clicked();

    void on_radioButton_NeutronPolar_clicked();

    void on_radioButton_NeutronRainbow_clicked();

    void on_horizontalSliderFPMoist_sliderMoved(int position);

    void on_horizontalSliderFPHum_sliderMoved(int position);

    void on_checkBoxFPLog_clicked();

    void on_horizontalSliderFPHum_2_sliderMoved(int position);

    void on_pushButton_ActivateFP_clicked();

    void on_lineEditScotoma_editingFinished();

    void on_checkBoxGradient2_clicked();

    void on_checkBoxSelectedData_clicked();

    void on_checkBoxFastData_clicked();

    void on_checkBoxIntermediateData_clicked();

    void on_checkBoxEpithermalData_clicked();

    //void on_lineEdit_Porosity_textChanged(const QString &arg1);

    void on_sliderSoilPorosity_sliderMoved(int position);

    void on_checkBoxFileOutputPDF_2_clicked(bool checked);

    void on_checkBoxCreateFolder_clicked(bool checked);

    void on_radioButton_AmBe_clicked();

    void on_lineEdit_xPos_textChanged(const QString &arg1);

    void on_lineEdit_yPos_textChanged(const QString &arg1);

    void on_lineEdit_zPos_textChanged(const QString &arg1);

    void on_lineEdit_zPos_2_textChanged(const QString &arg1);

    void on_checkBoxThermal_toggled(bool checked);

    void on_radioButton_mapTrack_clicked();

    void on_lineEditManualColorZero_textChanged(const QString &arg1);

    void on_lineEditManualColor_textChanged(const QString &arg1);

    void on_radioButton_mapTrackInter_clicked();

    void on_radioButton_mapTrackFast_clicked();

    void on_radioButton_mapTrackAlbedo_clicked();

    void on_checkBoxNoTrack_toggled(bool checked);

    void on_radioButton_mapThermal_clicked();

    void on_radioButton_mapTrackThermal_clicked();

    void on_checkBoxThermalData_clicked();

    void on_checkBoxSaveEvery_toggled(bool checked);

    void on_radioButton_NoSource_clicked();

    void on_lineEdit_xSize_textChanged(const QString &arg1);

    void on_lineEdit_ySize_textChanged(const QString &arg1);

    void on_radioButton_ThermalSource_clicked();

    void on_radioButton_MonoenergeticSource_clicked();

    void on_lineEdit_SourceEnergy_textChanged(const QString &arg1);

    void on_lineEditScotoma_2_textChanged(const QString &arg1);

    void on_checkBox_clicked();

    void on_checkBox_NoMultipleScattering_toggled(bool checked);

    void on_checkBoxManual_toggled(bool checked);

    void on_checkBox_TrackAllLayers_toggled(bool checked);

    void on_checkBoxLogarithmic_toggled(bool checked);

    void on_radioButton_NeutronGrayScale_clicked();

    void on_radioButton_NeutronHot_clicked();

    void on_radioButton_NeutronThermal_clicked();

    void on_radioButton_ModeratedCf_clicked();

    void on_radioButton_Cylinder_toggled(bool checked);

    void on_radioButton_Sphere_toggled(bool checked);

    void on_radioButton_detectorLayerEnergyBand_toggled(bool checked);

    void on_radioButton_detectorLayerRealistic_toggled(bool checked);

    void on_radioButton_detectorEnergyBand_toggled(bool checked);

    void on_radioButton_detectorRealistic_toggled(bool checked);

    void on_checkBox_VolumeSource_toggled(bool checked);

    void on_checkBox_HEModel_toggled(bool checked);

    void on_checkBox_activateThermal_toggled(bool checked);

    void on_horizontalSliderDetector_sliderMoved(int position);

    void on_horizontalSliderDetectorColor_sliderMoved(int position);

    void exportSettings(string str);

    bool importSettings();

    void on_pushButton_Enlarge_clicked();

    void on_lineEditAntiScotoma_editingFinished();

    void on_checkBoxSaveEvery_2_clicked(bool checked);

    void on_lineEditClearEveryXNeutrons_editingFinished();

    void on_lineEditClearEveryXNeutrons_textChanged(const QString &arg1);

    void on_lineEdit_AutoUpdate_textChanged(const QString &arg1);

    void on_checkBoxAutoRefreshRate_toggled(bool checked);

    void on_checkBoxClearEveryDisplayRefresh_toggled(bool checked);

    void on_checkBoxClearEveryDisplayRefresh_clicked(bool checked);

    void on_checkBoxWarnMAterial_toggled(bool checked);

    void on_radioButton_mapTrackEnergy_clicked();



    void on_checkBoxTrackingData_clicked(bool checked);

    void on_checkBoxHighResTrackingData_clicked(bool checked);

    void on_radioButton_ySheet_clicked();

    void on_radioButton_xSheet_clicked();

    void on_lineEditDetLength_textChanged(const QString &arg1);

    void on_lineEdit_zSize_textChanged(const QString &arg1);

    void on_radioButton_TopToBottom_toggled(bool checked);

    void on_radioButton_BottomToTop_toggled(bool checked);

    void on_radioButton_LeftToRight_toggled(bool checked);

    void on_radioButton_RightToLeft_toggled(bool checked);

    void on_radioButton_Omni_toggled(bool checked);

    void on_checkBox_DomainCutoff_toggled(bool checked);

    void on_lineEditDomainFactor_textChanged(const QString &arg1);

    void on_checkBox_DomainCutoffMeters_toggled(bool checked);

    void on_lineEditDomainMeters_textChanged(const QString &arg1);

    void on_checkBoxFileOutput2_toggled(bool checked);

    void on_checkBoxExportAllTracks_toggled(bool checked);

    void on_checkBox_ReflectiveBoundaries_toggled(bool checked);

    void on_checkBox_PeriodicBoundaries_toggled(bool checked);

    void on_lineEdit_DetectorFile_editingFinished();

    void on_lineEdit_DetectorFile_textChanged(const QString &arg1);

    void on_checkBoxFileOutput_toggled(bool checked);

    void on_checkBoxFileOutput3_toggled(bool checked);

    void on_pushButton_SaveConfig_clicked();

    void on_pushButton_Enlarge2_clicked();

private:
    Ui::MainWindow *ui;
    Ui::DialogShowPic *dialogshowpic;
    bool isFinishedVar;
    //Ui::VisualizationEnlarge *visualizationEnlarge;
};




bool cosmicNSimulator(MainWindow *ui);

double rangeFunctionCombined(double* x, double *par);

void delay( int millisecondsToWait );

void deleteStatsTH(TH1 *allTHs);
void logYaxis(TH2* h);
void logXaxis(TH2* h);
void logaxis(TH1* h);
void rebinX(TH1* h) ;
double getWWProb(double density, double atWeight, double cs, double lambda);
double getLfromE(double energy);
double getEfromL(double lambda);
void printTF1(TF1* function, TString outputFolder, TString filename);
bool checkDetectorHit(double x, double y, double z, double deltaZ, double phi, double theta);
double getEvaporationEnergy(double theta, TRandom * r);
double getFissionEnergy(TRandom * r);
vector<float> getThermalPDF(const double nEnergy, const float massElm, const float temperature, TRandom * r);
double getThermalEnergy(TF1 *spectrumFunc, TRandom * r);
double getDistanceToPoint(double stVx, double stVy, double stVz, double theta, double phi, double px, double py, double pz);
Double_t legendrian(double* x, Double_t* par);
Double_t legendrian10fold(double* x, Double_t* par);
string legendrian10folded(Float_t* par);
static int ArraySize(double array[]);
//double legendrianNfold(double* x, Double_t* par);
TMatrixF readSigmaEnergy(TString folder, TString filename);
void modifyCSmatrix(TMatrixF* sigmaMatrix, float factorLowE, float factorHighE);
float endfNumberConv(string str);
vector<TMatrixF> readAngularTabulatedCoefficients(TString folder, TString filename, const int coefficients);
TMatrixF readAngularCoefficients(TString folder, TString filename);
float getIndexHorizontalPosition(const TMatrixF& matrix, int line, double value, bool doLogSearch);
float getIndexPosition(const TMatrixF& matrix, double value, bool doLogSearch);
double calcMeanCS(const TMatrixF& matrix, double energy);
double getHighEnergyCosTheta(const TMatrixF& angleMatrix, const TMatrixF& cumulatedProbMatrix, double energy, double prob);
double getAngleFromCumulativeFunction(const TF1* spectrumFunc, float min, float max, TRandom* r);
TF1* calcMeanAngularDistribution(const TMatrixF& matrix, double energy);
void turnInputMatrix(TMatrixF& inputMatrix);
int changeInputMatrixValue(int value);
void generateNormalizedSpectrum(TH1F* precalculatedSpectrum, int entries, TString file1);
void dataDelete(QVector<double>* vec);



#endif // MAINWINDOW_H
