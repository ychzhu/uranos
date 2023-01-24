/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#include "dialogshowpic.h"
#include "ui_dialogshowpic.h"

#include "qcustomplot.h"
#include "mainwindow.h"
#include "Toolkit.h"

DialogShowPic::DialogShowPic(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogShowPic)
{
    ui->setupUi(this);

    setupShowGraph(ui->widget);
}

int number=1;

void DialogShowPic::setupShowGraph(QCustomPlot *customPlot)
{
   TMatrixF matrix = MainWindow::getTMatrix(number);

   QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis,customPlot->yAxis);
   customPlot->axisRect()->setupFullAxesBox(true);
   customPlot->xAxis->setLabel("x [pixel]");
   customPlot->yAxis->setLabel("y [pixel]");

   //customPlot->addPlottable(colorMap);
   colorMap->data()->setSize(matrix.GetColUpb(), matrix.GetColUpb());
   colorMap->data()->setRange(QCPRange(0, matrix.GetColUpb()), QCPRange(0, matrix.GetColUpb()));

   for (int x=0; x<matrix.GetColUpb(); ++x)
     for (int y=0; y<matrix.GetColUpb(); ++y)
       colorMap->data()->setCell(x, y, matrix(x,y));
   colorMap->setGradient(QCPColorGradient::gpJet);
   colorMap->rescaleDataRange(true);
   customPlot->rescaleAxes();
   colorMap->setDataRange(QCPRange(0,256));
   customPlot->replot();
}


DialogShowPic::~DialogShowPic()
{
    delete ui;
}

void DialogShowPic::replotGraph()
{
     TMatrixF matrix = MainWindow::getTMatrix(number);

     ui->widget->clearPlottables();

     ui->widget->xAxis->setRange(QCPRange(0, matrix.GetColUpb()));
     ui->widget->yAxis->setRange(QCPRange(0, matrix.GetColUpb()));
     QCPColorMap *colorMap = new QCPColorMap( ui->widget->xAxis,  ui->widget->yAxis);

    //Int_t a = matrix.GetColUpb();
    //Int_t b = matrix.GetNcols();

    colorMap->data()->setSize(matrix.GetNcols(), matrix.GetNcols());
    colorMap->data()->setRange(QCPRange(0, matrix.GetColUpb()), QCPRange(0, matrix.GetColUpb()));

     for (int x=0; x<matrix.GetColUpb(); ++x)
       for (int y=0; y<matrix.GetColUpb(); ++y)
         colorMap->data()->setCell(x, y, matrix(x,y));
     colorMap->setGradient(QCPColorGradient::gpCold);
     colorMap->rescaleDataRange(true);
     ui->widget->rescaleAxes();
     colorMap->setDataRange(QCPRange(0,256));
     ui->widget->replot();
}

void DialogShowPic::on_pushButtonPlus_clicked()
{
    number++;

    int newInt = number+1;
    string numberString = castIntToString(newInt);
    ui->labelNo->setText(QString::fromStdString(numberString));
    replotGraph();
}

void DialogShowPic::on_pushButtonMinus_clicked()
{
    if (number>0)
    {
        number--;
        int newInt = number+1;
        string numberString = castIntToString(newInt);
        ui->labelNo->setText(QString::fromStdString(numberString));
        replotGraph();
    }
}


