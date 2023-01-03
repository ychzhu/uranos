/********************************************************************************
** Form generated from reading UI file 'visualizationenlarge.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VISUALIZATIONENLARGE_H
#define UI_VISUALIZATIONENLARGE_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QPushButton>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_VisualizationEnlarge
{
public:
    QCustomPlot *customPlot;
    QPushButton *pushButtonPNG;

    void setupUi(QDialog *VisualizationEnlarge)
    {
        if (VisualizationEnlarge->objectName().isEmpty())
            VisualizationEnlarge->setObjectName(QString::fromUtf8("VisualizationEnlarge"));
        VisualizationEnlarge->resize(1212, 1052);
        VisualizationEnlarge->setMaximumSize(QSize(1250, 1052));
        customPlot = new QCustomPlot(VisualizationEnlarge);
        customPlot->setObjectName(QString::fromUtf8("customPlot"));
        customPlot->setGeometry(QRect(0, 0, 1211, 1052));
        customPlot->setMaximumSize(QSize(1250, 1052));
        pushButtonPNG = new QPushButton(customPlot);
        pushButtonPNG->setObjectName(QString::fromUtf8("pushButtonPNG"));
        pushButtonPNG->setGeometry(QRect(1140, 10, 61, 41));

        retranslateUi(VisualizationEnlarge);

        QMetaObject::connectSlotsByName(VisualizationEnlarge);
    } // setupUi

    void retranslateUi(QDialog *VisualizationEnlarge)
    {
        VisualizationEnlarge->setWindowTitle(QCoreApplication::translate("VisualizationEnlarge", "Enlarged Visualization", nullptr));
        pushButtonPNG->setText(QCoreApplication::translate("VisualizationEnlarge", "URANOS\n"
"Export", nullptr));
    } // retranslateUi

};

namespace Ui {
    class VisualizationEnlarge: public Ui_VisualizationEnlarge {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VISUALIZATIONENLARGE_H
