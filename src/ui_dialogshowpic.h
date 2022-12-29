/********************************************************************************
** Form generated from reading UI file 'dialogshowpic.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOGSHOWPIC_H
#define UI_DIALOGSHOWPIC_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_DialogShowPic
{
public:
    QDialogButtonBox *buttonBox;
    QCustomPlot *widget;
    QPushButton *pushButtonPlus;
    QPushButton *pushButtonMinus;
    QLabel *labelNo;
    QLabel *label;

    void setupUi(QDialog *DialogShowPic)
    {
        if (DialogShowPic->objectName().isEmpty())
            DialogShowPic->setObjectName(QString::fromUtf8("DialogShowPic"));
        DialogShowPic->resize(627, 590);
        buttonBox = new QDialogButtonBox(DialogShowPic);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(430, 550, 171, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        widget = new QCustomPlot(DialogShowPic);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(20, 20, 591, 521));
        pushButtonPlus = new QPushButton(DialogShowPic);
        pushButtonPlus->setObjectName(QString::fromUtf8("pushButtonPlus"));
        pushButtonPlus->setGeometry(QRect(60, 555, 31, 23));
        pushButtonMinus = new QPushButton(DialogShowPic);
        pushButtonMinus->setObjectName(QString::fromUtf8("pushButtonMinus"));
        pushButtonMinus->setGeometry(QRect(20, 555, 31, 23));
        labelNo = new QLabel(DialogShowPic);
        labelNo->setObjectName(QString::fromUtf8("labelNo"));
        labelNo->setGeometry(QRect(163, 559, 46, 13));
        label = new QLabel(DialogShowPic);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(107, 557, 51, 16));

        retranslateUi(DialogShowPic);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogShowPic, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogShowPic, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogShowPic);
    } // setupUi

    void retranslateUi(QDialog *DialogShowPic)
    {
        DialogShowPic->setWindowTitle(QCoreApplication::translate("DialogShowPic", "Input Matrix View", nullptr));
        pushButtonPlus->setText(QCoreApplication::translate("DialogShowPic", "+", nullptr));
        pushButtonMinus->setText(QCoreApplication::translate("DialogShowPic", "-", nullptr));
        labelNo->setText(QCoreApplication::translate("DialogShowPic", "-", nullptr));
        label->setText(QCoreApplication::translate("DialogShowPic", "Layer #:", nullptr));
    } // retranslateUi

};

namespace Ui {
    class DialogShowPic: public Ui_DialogShowPic {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOGSHOWPIC_H
