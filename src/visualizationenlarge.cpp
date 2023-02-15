/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#include "visualizationenlarge.h"
#include "ui_visualizationenlarge.h"


// variables for drawing the enlarged visulization view
float squareDimV = 10;

int horizontalSliderValue = 50;
int horizontalSliderColorZeroValue = 0;
int manualColorZeroV = 0;
int manualColorV = 0;
int colorScheme = 0;
int screenshotcounter = 1;

bool plotTopViewLogV = false;
bool silderColorMovedV = false;
bool useManualColor = false;

float xCustomPosV, yCustomPosV;

string workFolderV;

VisualizationEnlarge::VisualizationEnlarge(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VisualizationEnlarge)
{
    ui->setupUi(this);
    setupRunGraph(ui->customPlot);

    ui->customPlot->installEventFilter(this);
}

VisualizationEnlarge::~VisualizationEnlarge()
{
    delete ui;
}

bool VisualizationEnlarge::eventFilter(QObject* target, QEvent* event)
{
    float x;
    float y;

    if (target == ui->customPlot && event->type() == QEvent::MouseButtonPress)
    {
        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);

        if (_mouseEvent->button() == Qt::RightButton)
        {
            x = ui->customPlot->xAxis->pixelToCoord(_mouseEvent->pos().x());
            y = ui->customPlot->yAxis->pixelToCoord(_mouseEvent->pos().y());

            //setCursorPos(x * 1000., y * 1000.);

            emit on_VisualizationEnlargeCursorX(x * 1000.);
            emit on_VisualizationEnlargeCursorY(y * 1000.);
            emit on_VisualizationEnlargePressed(true);
        }
    }

    if (target == ui->customPlot && event->type() == QEvent::MouseMove)
    {
        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);

        x = ui->customPlot->xAxis->pixelToCoord(_mouseEvent->pos().x());
        y = ui->customPlot->yAxis->pixelToCoord(_mouseEvent->pos().y());

        //xCustomPosV = x * 1000.;
        //yCustomPosV = y * 1000.;

        //setCursorPos(x * 1000., y * 1000.);

        emit on_VisualizationEnlargeCursorX(x * 1000.);
        emit on_VisualizationEnlargeCursorY(y * 1000.);
    }

    if (target == ui->customPlot && event->type() == QEvent::MouseButtonRelease)
    {
        QMouseEvent* _mouseEvent = static_cast<QMouseEvent*>(event);

        if (_mouseEvent->button() == Qt::RightButton)
        {
            emit on_VisualizationEnlargeReleased(true);
            emit on_VisualizationEnlargeState(false);
        }
    }

    return false;
}


void VisualizationEnlarge::setSquareDimSize(float squareDimSize)
{
    squareDimV = squareDimSize;
}

void VisualizationEnlarge::sethorizontalSliderValue(int value)
{
    horizontalSliderValue = value;
}

void VisualizationEnlarge::sethorizontalSliderColorZeroValue(int value)
{
    horizontalSliderColorZeroValue = value;
}

void VisualizationEnlarge::setmanualColorZero(int value)
{
    manualColorZeroV = value;
}

void VisualizationEnlarge::setmanualColor(int value)
{
    manualColorV = value;
}

void VisualizationEnlarge::setsilderColorMoved(bool value)
{
    silderColorMovedV = value;
}

void VisualizationEnlarge::setplotTopViewLog(bool value)
{
    plotTopViewLogV = value;
}

void VisualizationEnlarge::setuseManualColor(bool value)
{
    useManualColor = value;
}

void VisualizationEnlarge::setColorScheme(int value)
{
    colorScheme = value;
}

void VisualizationEnlarge::setWorkFolder(string text)
{
    workFolderV = text;
}

void VisualizationEnlarge::plotGraph(TH2F* data, int size, float squareDimSize)
{
    setSquareDimSize(squareDimSize);

    float entryWeight, maximumWeight;
    float pixelsum;

    QColor lightBlue = QColor (215, 237, 255);

    ui->customPlot->clearPlottables();
    ui->customPlot->plotLayout()->remove(ui->customPlot->plotLayout()->element(0,1));

    int elementCount = ui->customPlot->plotLayout()->elementCount();

    for(int i = 0; i < elementCount; i++)
    {
        if(qobject_cast<QCPLayoutGrid*>(ui->customPlot->plotLayout()->elementAt(i)))  ui->customPlot->plotLayout()->removeAt(i);
    }
    ui->customPlot->plotLayout()->simplify();

    QCPColorMap *colorMap = new QCPColorMap( ui->customPlot->xAxis,  ui->customPlot->yAxis);

    if (size == 500)  colorMap->data()->setSize(500, 500);
    if (size == 1000)  colorMap->data()->setSize(1000, 1000);
    if (size == 2000)  colorMap->data()->setSize(2000, 2000);
    if (size == 4000)  colorMap->data()->setSize(1000, 1000);
    colorMap->data()->setRange(QCPRange(-::squareDimV/2./1000., ::squareDimV/2./1000.), QCPRange(-::squareDimV/2./1000., ::squareDimV/2./1000.));

    if (plotTopViewLogV) colorMap->setDataScaleType(QCPAxis::stLogarithmic);
    else colorMap->setDataScaleType(QCPAxis::stLinear);

    QCPColorScale *colorScale = new QCPColorScale(ui->customPlot);

    colorScale->setType(QCPAxis::atRight);
    colorMap->setColorScale(colorScale);
    colorScale->axis()->setLabel("Neutron Hit Density");

    if (plotTopViewLogV) colorScale->setDataScaleType(QCPAxis::stLogarithmic);
    else colorScale->setDataScaleType(QCPAxis::stLinear);

    colorScale->axis()->setSubTickPen(QPen(lightBlue));
    colorScale->axis()->setTickPen(QPen(lightBlue));
    colorScale->axis()->setTickLabelColor(lightBlue);
    colorScale->axis()->setLabelColor(lightBlue);

    QCPLayoutGrid *subLayout = new QCPLayoutGrid;
    ui->customPlot->plotLayout()->addElement(0, 1, subLayout);

    subLayout->setMargins(QMargins(10, 10, 10, 5));
    subLayout->addElement(0, 0, colorScale);

    ui->customPlot->plotLayout()->element(0,0)->setMaximumSize(1000,1000);

    float minZ = 0;

    if (plotTopViewLogV) minZ = 1;

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

    float colorscaleMax = (horizontalSliderValue*1.)/99.*maximumWeight;
    if (colorscaleMax < 1) colorscaleMax = 1;

    if (colorScheme==0) {colorMap->setGradient(QCPColorGradient::gpJet); colorScale->setGradient(QCPColorGradient::gpJet); }
    if (colorScheme==1) {colorMap->setGradient(QCPColorGradient::gpNight); colorScale->setGradient(QCPColorGradient::gpNight); }
    if (colorScheme==2) {colorMap->setGradient(QCPColorGradient::gpCold); colorScale->setGradient(QCPColorGradient::gpCold);}
    if (colorScheme==3) {colorMap->setGradient(QCPColorGradient::gpThermal); colorScale->setGradient(QCPColorGradient::gpThermal);}
    if (colorScheme==4) {colorMap->setGradient(QCPColorGradient::gpHot); colorScale->setGradient(QCPColorGradient::gpHot);}
    if (colorScheme==5) {colorMap->setGradient(QCPColorGradient::gpPolar); colorScale->setGradient(QCPColorGradient::gpPolar);}
    if (colorScheme==6) {colorMap->setGradient(QCPColorGradient::gpGrayscale); colorScale->setGradient(QCPColorGradient::gpGrayscale);}

    colorMap->setInterpolate(true);

    if (manualColorZeroV < minZ) { manualColorZeroV = minZ;}

    if (!silderColorMovedV)
    {
        if (entryWeight*1.1 < 5) {colorMap->setDataRange(QCPRange(minZ,5));  colorScale->setDataRange(QCPRange(minZ,5));}
        else {colorMap->setDataRange(QCPRange(minZ,entryWeight*1.1)); colorScale->setDataRange(QCPRange(minZ,entryWeight*1.1));}
    }
    else
    {
        if (horizontalSliderColorZeroValue == 0) { colorMap->setDataRange(QCPRange( minZ , colorscaleMax ));  colorScale->setDataRange(QCPRange(minZ,colorscaleMax));}
        else { colorMap->setDataRange(QCPRange( (horizontalSliderColorZeroValue*1.)/99.*colorscaleMax , colorscaleMax ));  colorScale->setDataRange(QCPRange((horizontalSliderColorZeroValue*1.)/99.*colorscaleMax,colorscaleMax));}
    }

    if (useManualColor)
    {
        colorMap->setDataRange(QCPRange( manualColorZeroV,manualColorV));
        colorScale->setDataRange(QCPRange(manualColorZeroV,manualColorV));
    }

    ui->customPlot->rescaleAxes();
    ui->customPlot->replot();
}

void VisualizationEnlarge::setupRunGraph(QCustomPlot *customPlot)
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

    colorMap->data()->setSize(500, 500);
    colorMap->data()->setRange(QCPRange(-::squareDimV/2./1000., ::squareDimV/2./1000.), QCPRange(-::squareDimV/2./1000., ::squareDimV/2./1000.));

    colorMap->setGradient(QCPColorGradient::gpJet);
    colorMap->rescaleDataRange(true);
    QCPRange dataR = colorMap->dataRange(); dataR*=2;
    colorMap->setDataRange(dataR);

    QCPColorScale *colorScale = new QCPColorScale(ui->customPlot);
    ui->customPlot->plotLayout()->addElement(0, 1, colorScale);
    colorScale->setType(QCPAxis::atRight);
    colorMap->setColorScale(colorScale);
    colorScale->axis()->setLabel("Neutron Flux");

    ui->customPlot->plotLayout()->element(0,0)->setMaximumSize(1000,1000);

    colorScale->axis()->setSubTickPen(QPen(lightBlue));
    colorScale->axis()->setTickPen(QPen(lightBlue));
    colorScale->axis()->setTickLabelColor(lightBlue);
    colorScale->axis()->setLabelColor(lightBlue);

    colorScale->setDataRange(dataR);

    customPlot->rescaleAxes();

    customPlot->replot();
}

void VisualizationEnlarge::on_pushButtonPNG_toggled(bool checked)
{
    //originalPixmap.save(QString::fromStdString(workFolderV)+fileName, format.toAscii());
}

void VisualizationEnlarge::on_pushButtonPNG_clicked()
{
    QScreen *screen = QGuiApplication::primaryScreen();

    QPixmap originalPixmap;

    originalPixmap = screen->grabWindow(ui->customPlot->winId());

    QString fileName = QString::fromStdString(workFolderV)+"/neutronHitDensityExport_"+QString::number(screenshotcounter)+".png";

    if (!originalPixmap.save(fileName))
    {
            QMessageBox::warning(this, tr("Save Error"), tr("The image could not be saved to \"%1\".").arg(QDir::toNativeSeparators(fileName)));
    }
    screenshotcounter++;
}

/*
/**
 * cursor placement on mouse event
 * @param xc, yc

void VisualizationEnlarge::setCursorPos(float x, float y)
{
   xCustomPosV = x;
   yCustomPosV = y;
}
*/

/*
void VisualizationEnlarge::on_VisualizationEnlarge_cursorX(float x)
{
    xCustomPosV = x;
}
*/
