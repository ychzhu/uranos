/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2022 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/


#ifndef DIALOGSHOWPIC_H
#define DIALOGSHOWPIC_H

#include <QDialog>
#include "qcustomplot.h"

namespace Ui {
class DialogShowPic;
}

class DialogShowPic : public QDialog
{
    Q_OBJECT

public:
    explicit DialogShowPic(QWidget *parent = 0);
    ~DialogShowPic();

    void setupShowGraph(QCustomPlot *customPlot);
    void replotGraph();

private slots:
    void on_pushButtonPlus_clicked();

    void on_pushButtonMinus_clicked();

private:
    Ui::DialogShowPic *ui;
};

#endif // DIALOGSHOWPIC_H
