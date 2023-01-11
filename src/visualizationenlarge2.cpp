/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2022 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#include "visualizationenlarge2.h"
#include "ui_visualizationEnlarge2.h"


// variables for drawing the enlarged visulization view
float squareDimV2 = 10;

int horizontalSliderValueV2 = 50;
int horizontalSliderColorZeroValueV2 = 0;
int manualColorZeroV2 = 0;
int manualColorV2 = 0;
int colorSchemeV2 = 0;
int screenshotcounterV2 = 1;

bool plotTopViewLogV2 = false;
bool silderColorMovedV2 = false;
bool useManualColorV2 = false;

string workFolderV2;

VisualizationEnlarge2::VisualizationEnlarge2(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VisualizationEnlarge2)
{
    ui->setupUi(this);
    setupRunGraph2(ui->customPlotEnlarge2);
}

VisualizationEnlarge2::~VisualizationEnlarge2()
{
    delete ui;
}

void VisualizationEnlarge2::setSquareDimSize2(float squareDimSize)
{
    squareDimV2 = squareDimSize;
}

void VisualizationEnlarge2::sethorizontalSliderValue2(int value)
{
    horizontalSliderValueV2 = value;
}

void VisualizationEnlarge2::sethorizontalSliderColorZeroValue2(int value)
{
    horizontalSliderColorZeroValueV2 = value;
}

void VisualizationEnlarge2::setmanualColorZero2(int value)
{
    manualColorZeroV2 = value;
}

void VisualizationEnlarge2::setmanualColor2(int value)
{
    manualColorV2 = value;
}

void VisualizationEnlarge2::setsilderColorMoved2(bool value)
{
    silderColorMovedV2 = value;
}

void VisualizationEnlarge2::setplotTopViewLog2(bool value)
{
    plotTopViewLogV2 = value;
}

void VisualizationEnlarge2::setuseManualColor2(bool value)
{
    useManualColorV2 = value;
}

void VisualizationEnlarge2::setColorScheme2(int value)
{
    colorSchemeV2 = value;
}

void VisualizationEnlarge2::setWorkFolder2(string text)
{
    workFolderV2 = text;
}

void VisualizationEnlarge2::plotGraph2(TH2F* data, int size, float squareDimSize)
{
    setSquareDimSize2(squareDimSize);

    float entryWeight, maximumWeight;
    float pixelsum;

    QColor lightBlue = QColor (215, 237, 255);

    ui->customPlotEnlarge2->clearPlottables();
    ui->customPlotEnlarge2->plotLayout()->remove(ui->customPlotEnlarge2->plotLayout()->element(0,1));

    int elementCount = ui->customPlotEnlarge2->plotLayout()->elementCount();

    for(int i = 0; i < elementCount; i++)
    {
        if(qobject_cast<QCPLayoutGrid*>(ui->customPlotEnlarge2->plotLayout()->elementAt(i)))  ui->customPlotEnlarge2->plotLayout()->removeAt(i);
    }
    ui->customPlotEnlarge2->plotLayout()->simplify();

    QCPColorMap *colorMap = new QCPColorMap( ui->customPlotEnlarge2->xAxis,  ui->customPlotEnlarge2->yAxis);

    if (size == 500)  colorMap->data()->setSize(500, 500);
    if (size == 1000)  colorMap->data()->setSize(1000, 1000);
    if (size == 1500)  colorMap->data()->setSize(1500, 1500);
    if (size == 2000)  colorMap->data()->setSize(2000, 2000);
    if (size == 4000)  colorMap->data()->setSize(1000, 1000);
    colorMap->data()->setRange(QCPRange(-::squareDimV2/2./1000., ::squareDimV2/2./1000.), QCPRange(-::squareDimV2/2./1000., ::squareDimV2/2./1000.));

    if (plotTopViewLogV2) colorMap->setDataScaleType(QCPAxis::stLogarithmic);
    else colorMap->setDataScaleType(QCPAxis::stLinear);

    QCPColorScale *colorScale = new QCPColorScale(ui->customPlotEnlarge2);

    colorScale->setType(QCPAxis::atRight);
    colorMap->setColorScale(colorScale);
    colorScale->axis()->setLabel("Neutron Hit Density");

    if (plotTopViewLogV2) colorScale->setDataScaleType(QCPAxis::stLogarithmic);
    else colorScale->setDataScaleType(QCPAxis::stLinear);

    colorScale->axis()->setSubTickPen(QPen(lightBlue));
    colorScale->axis()->setTickPen(QPen(lightBlue));
    colorScale->axis()->setTickLabelColor(lightBlue);
    colorScale->axis()->setLabelColor(lightBlue);

    QCPLayoutGrid *subLayout = new QCPLayoutGrid;
    ui->customPlotEnlarge2->plotLayout()->addElement(0, 1, subLayout);

    subLayout->setMargins(QMargins(10, 10, 10, 5));
    subLayout->addElement(0, 0, colorScale);

    ui->customPlotEnlarge2->plotLayout()->element(0,0)->setMaximumSize(1500,1500);

    float minZ = 0;

    if (plotTopViewLogV2) minZ = 1;

    if (true)
    {
        if (size == 500)
        {
            for (int x=0; x<500; x+=1)
            {
                for (int y=0; y<500; y+=1)
                {
                    colorMap->data()->setCell((int)x, (int)y, data->GetBinContent(x,y));
                }
            }
            entryWeight = (data->GetEntries())/(500.*500.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }

        if (size == 1000)
        {
            for (int x=0; x<1000; x+=1)
            {
                for (int y=0; y<1000; y+=1)
                {
                    colorMap->data()->setCell((int)x, (int)y, data->GetBinContent(x,y));
                }
            }
            entryWeight = (data->GetEntries())/(1000.*1000.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }

        if (size == 1500)
        {
            for (int x=0; x<1500; x+=1)
            {
                for (int y=0; y<1500; y+=1)
                {
                    colorMap->data()->setCell((int)x, (int)y, data->GetBinContent(x,y));
                }
            }
            entryWeight = (data->GetEntries())/(1500.*1500.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }

        if (size == 2000)
        {
            for (int x=0; x<2000; x+=1)
            {
                for (int y=0; y<2000; y+=1)
                {
                    colorMap->data()->setCell((int)x, (int)y, data->GetBinContent(x,y));
                }
            }
            entryWeight = (data->GetEntries())/(2000.*2000.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }

        if (size == 2020)
        {
            for (int x=0; x<1000; x+=1)
            {
                for (int y=0; y<1000; y+=1)
                {
                    pixelsum = 0;
                    for (int xadd=0; xadd<2; xadd+=1)
                    {
                        for (int yadd=0; yadd<2; yadd+=1)
                        {
                            pixelsum+=data->GetBinContent(2*x+xadd,2*y+yadd);
                        }
                    }
                    pixelsum = pixelsum/4.;

                    colorMap->data()->setCell((int)x, (int)y, pixelsum);
                }
            }
            entryWeight = (data->GetEntries())/(1000.*1000.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }

        if (size == 4000)
        {
            for (int x=0; x<1000; x+=1)
            {
                for (int y=0; y<1000; y+=1)
                {
                    pixelsum = 0;
                    for (int xadd=0; xadd<4; xadd+=1)
                    {
                        for (int yadd=0; yadd<4; yadd+=1)
                        {
                            pixelsum+=data->GetBinContent(4*x+xadd,4*y+yadd);
                        }
                    }
                    pixelsum = pixelsum/16.;

                    colorMap->data()->setCell((int)x, (int)y, pixelsum);
                }
            }
            entryWeight = (data->GetEntries())/(1000.*1000.)*4.; maximumWeight = 1.* data->GetBinContent(data->GetMaximumBin());
        }
    }

    float colorscaleMax = (horizontalSliderValueV2*1.)/99.*maximumWeight;
    if (colorscaleMax < 1) colorscaleMax = 1;

    if (colorSchemeV2==0) {colorMap->setGradient(QCPColorGradient::gpJet); colorScale->setGradient(QCPColorGradient::gpJet);}
    if (colorSchemeV2==1) {colorMap->setGradient(QCPColorGradient::gpNight); colorScale->setGradient(QCPColorGradient::gpNight);}
    if (colorSchemeV2==2) {colorMap->setGradient(QCPColorGradient::gpCold); colorScale->setGradient(QCPColorGradient::gpCold);}
    if (colorSchemeV2==3) {colorMap->setGradient(QCPColorGradient::gpThermal); colorScale->setGradient(QCPColorGradient::gpThermal);}
    if (colorSchemeV2==4) {colorMap->setGradient(QCPColorGradient::gpHot); colorScale->setGradient(QCPColorGradient::gpHot);}
    if (colorSchemeV2==5) {colorMap->setGradient(QCPColorGradient::gpPolar); colorScale->setGradient(QCPColorGradient::gpPolar);}
    if (colorSchemeV2==6) {colorMap->setGradient(QCPColorGradient::gpGrayscale); colorScale->setGradient(QCPColorGradient::gpGrayscale);}

    colorMap->setInterpolate(true);

    if (manualColorZeroV2 < minZ) {manualColorZeroV2 = minZ;}

    if (!silderColorMovedV2)
    {
        if (entryWeight*1.1 < 5) {colorMap->setDataRange(QCPRange(minZ,5));  colorScale->setDataRange(QCPRange(minZ,5));}
        else {colorMap->setDataRange(QCPRange(minZ,entryWeight*1.1)); colorScale->setDataRange(QCPRange(minZ,entryWeight*1.1));}
    }
    else
    {
        if (horizontalSliderColorZeroValueV2 == 0) { colorMap->setDataRange(QCPRange(minZ , colorscaleMax ));  colorScale->setDataRange(QCPRange(minZ,colorscaleMax));}
        else { colorMap->setDataRange(QCPRange( (horizontalSliderColorZeroValueV2*1.)/99.*colorscaleMax , colorscaleMax ));  colorScale->setDataRange(QCPRange((horizontalSliderColorZeroValueV2*1.)/99.*colorscaleMax,colorscaleMax));}
    }

    if (useManualColorV2)
    {
        colorMap->setDataRange(QCPRange(manualColorZeroV2,manualColorV2));
        colorScale->setDataRange(QCPRange(manualColorZeroV2,manualColorV2));
    }

    ui->customPlotEnlarge2->rescaleAxes();
    ui->customPlotEnlarge2->replot();
}

void VisualizationEnlarge2::setupRunGraph2(QCustomPlot *customPlot)
{
    //QColor lightBlue = QColor ( 234 , 246 , 255);
    QColor lightBlue = QColor (215, 237, 255);
    QColor someBlue = QColor (0, 88, 156);
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
    customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom);

    QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis,customPlot->yAxis);
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->xAxis->setLabel("x [m]");
    customPlot->yAxis->setLabel("y [m]");

    colorMap->data()->setSize(1000, 1000);
    colorMap->data()->setRange(QCPRange(-::squareDimV2/2./1000., ::squareDimV2/2./1000.), QCPRange(-::squareDimV2/2./1000., ::squareDimV2/2./1000.));

    colorMap->setGradient(QCPColorGradient::gpJet);
    colorMap->rescaleDataRange(true);
    QCPRange dataR = colorMap->dataRange(); dataR*=2;
    colorMap->setDataRange(dataR);

    QCPColorScale *colorScale = new QCPColorScale(ui->customPlotEnlarge2);
    ui->customPlotEnlarge2->plotLayout()->addElement(0, 1, colorScale);
    colorScale->setType(QCPAxis::atRight);
    colorMap->setColorScale(colorScale);
    colorScale->axis()->setLabel("Neutron Flux");

    ui->customPlotEnlarge2->plotLayout()->element(0,0)->setMaximumSize(1500,1500);

    colorScale->axis()->setSubTickPen(QPen(lightBlue));
    colorScale->axis()->setTickPen(QPen(lightBlue));
    colorScale->axis()->setTickLabelColor(lightBlue);
    colorScale->axis()->setLabelColor(lightBlue);

    colorScale->setDataRange(dataR);

    customPlot->rescaleAxes();

    customPlot->replot();
}

/*
void VisualizationEnlarge2::on_pushButtonPNG_toggled(bool checked)
{
    //originalPixmap.save(QString::fromStdString(workFolderV)+fileName, format.toAscii());}

*/

void VisualizationEnlarge2::on_pushButtonPNG2_clicked()
{
    QScreen *screen = QGuiApplication::primaryScreen();

    QPixmap originalPixmap;

    originalPixmap = screen->grabWindow(ui->customPlotEnlarge2->winId());

    QString fileName = QString::fromStdString(workFolderV2)+"/neutronHighResExport_"+QString::number(screenshotcounterV2)+".png";

    if (!originalPixmap.save(fileName))
    {
            QMessageBox::warning(this, tr("Save Error"), tr("The image could not be saved to \"%1\".").arg(QDir::toNativeSeparators(fileName)));
    }
    screenshotcounterV2++;
}

