/********************************************************************************
** Form generated from reading UI file 'visualizationEnlarge2.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VISUALIZATIONENLARGE2_H
#define UI_VISUALIZATIONENLARGE2_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QPushButton>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_VisualizationEnlarge2
{
public:
    QCustomPlot *customPlotEnlarge2;
    QPushButton *pushButtonPNG2;

    void setupUi(QDialog *VisualizationEnlarge2)
    {
        if (VisualizationEnlarge2->objectName().isEmpty())
            VisualizationEnlarge2->setObjectName(QString::fromUtf8("VisualizationEnlarge2"));
        VisualizationEnlarge2->resize(1712, 1520);
        customPlotEnlarge2 = new QCustomPlot(VisualizationEnlarge2);
        customPlotEnlarge2->setObjectName(QString::fromUtf8("customPlotEnlarge2"));
        customPlotEnlarge2->setGeometry(QRect(-1, 0, 1710, 1530));
        pushButtonPNG2 = new QPushButton(customPlotEnlarge2);
        pushButtonPNG2->setObjectName(QString::fromUtf8("pushButtonPNG2"));
        pushButtonPNG2->setGeometry(QRect(1634, 20, 61, 41));

        retranslateUi(VisualizationEnlarge2);

        QMetaObject::connectSlotsByName(VisualizationEnlarge2);
    } // setupUi

    void retranslateUi(QDialog *VisualizationEnlarge2)
    {
        VisualizationEnlarge2->setWindowTitle(QCoreApplication::translate("VisualizationEnlarge2", "Dialog", nullptr));
        pushButtonPNG2->setText(QCoreApplication::translate("VisualizationEnlarge2", "URANOS\n"
"Export", nullptr));
    } // retranslateUi

};

namespace Ui {
    class VisualizationEnlarge2: public Ui_VisualizationEnlarge2 {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VISUALIZATIONENLARGE2_H
