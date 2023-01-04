/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2022 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#ifndef VISUALIZATIONENLARGE2_H
#define VISUALIZATIONENLARGE2_H

#include <QDialog>

#include "qcustomplot.h"

#include "Toolkit.h"

namespace Ui {
class VisualizationEnlarge2;
}

class VisualizationEnlarge2 : public QDialog
{
    Q_OBJECT
    QMessageBox* msgBox2;

public:
    explicit VisualizationEnlarge2(QWidget *parent = 0);
    ~VisualizationEnlarge2();

    void setSquareDimSize2(float squareDimSize);

    void sethorizontalSliderValue2(int value);

    void sethorizontalSliderColorZeroValue2(int value);

    void setmanualColorZero2(int value);

    void setmanualColor2(int value);

    void setsilderColorMoved2(bool value);

    void setplotTopViewLog2(bool value);

    void setuseManualColor2(bool value);

    void setColorScheme2(int value);

    void setWorkFolder2(string text);

    void plotGraph2(TH2F* data, int size, float squareDimSize);

    void setupRunGraph2(QCustomPlot *customPlot);

private slots:
    /*
    void on_pushButtonPNG_toggled(bool checked);
*/
    void on_pushButtonPNG2_clicked();

private:
    Ui::VisualizationEnlarge2 *ui;

};

#endif // VISUALIZATIONENLARGE2_H
