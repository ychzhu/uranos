/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#ifndef VISUALIZATIONENLARGE_H
#define VISUALIZATIONENLARGE_H

#include <QDialog>

#include <qcustomplot.h>

#include "Toolkit.h"

namespace Ui {
class VisualizationEnlarge;
}

class VisualizationEnlarge : public QDialog
{
    Q_OBJECT
    QMessageBox* msgBox;

public:
    explicit VisualizationEnlarge(QWidget *parent = 0);
    ~VisualizationEnlarge();

    void setSquareDimSize(float squareDimSize);

    void sethorizontalSliderValue(int value);

    void sethorizontalSliderColorZeroValue(int value);

    void setmanualColorZero(int value);

    void setmanualColor(int value);

    void setsilderColorMoved(bool value);

    void setplotTopViewLog(bool value);

    void setuseManualColor(bool value);

    void setColorScheme(int value);

    void setWorkFolder(string text);

    void plotGraph(TH2F* data, int size, float squareDimSize);

    void setupRunGraph(QCustomPlot *customPlot);

    //void setFocus(MainWindow* wF);

signals:

  //void setCursorPos(float xc, float yc);

  void on_VisualizationEnlargeCursorX(float x);

  void on_VisualizationEnlargeCursorY(float y);

  void on_VisualizationEnlargePressed(bool var);

  void on_VisualizationEnlargeReleased(bool var);

  void on_VisualizationEnlargeState(bool var);

private slots:    

    bool eventFilter(QObject *target, QEvent *event);

    void on_pushButtonPNG_toggled(bool checked);

    void on_pushButtonPNG_clicked();

private:
    Ui::VisualizationEnlarge *ui;
};

#endif // VISUALIZATIONENLARGE_H
