#-------------------------------------------------
#
# Project created by QtCreator 2016-06-13T21:45:52
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = UranosGUI
TEMPLATE = app

CONFIG += console

#CONFIG   -= shared
#CONFIG   += static
#CONFIG   += staticlib

#QMAKE_LFLAGS		= -static -enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc

RC_FILE = urc.rc

SOURCES += main.cpp\
        mainwindow.cpp \#
    qcustomplot.cpp \
    Toolkit.cpp \
    dialogshowpic.cpp \
    customSplashScreen.cpp \
    visualizationenlarge.cpp


HEADERS  += mainwindow.h \
    qcustomplot.h \
    Toolkit.h \
    dialogshowpic.h \
    customSplashScreen.h \
    visualizationenlarge.h

FORMS    += mainwindow.ui \
    dialogshowpic.ui \
    visualizationenlarge.ui

DEFINES += _CRT_SECURE_NO_WARNINGS
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0

QMAKE_LFLAGS +=  /FORCE
QMAKE_LFLAGS += /INCREMENTAL:NO

win32:CXXFLAGS += /O2

QMAKE_CXXFLAGS += -openmp
QMAKE_CXXFLAGS += -MP


#LIBS += -fopenmp

#LIBS += -LC:/Qt2/5.12.11/msvc2017/ -lQt5Core
win32: LIBS += -L$$PWD/root/lib/ -llibMatrix
win32: LIBS += -L$$PWD/root/lib/ -llibCore
win32: LIBS += -L$$PWD/root/lib/ -llibMathCore
win32: LIBS += -L$$PWD/root/lib/ -llibGui
win32: LIBS += -L$$PWD/root/lib/ -llibHist
#win32: LIBS += -L$$PWD/root/lib/ -llibCint
win32: LIBS += -L$$PWD/root/lib/ -llibGraf
win32: LIBS += -L$$PWD/root/lib/ -llibRIO
win32: LIBS += -L$$PWD/root/lib/ -llibGPad
#win32: LIBS += -L$$PWD/GnuWin32/lib -llibgsl
#win32: LIBS += -L$$PWD/gsl/lib -lgsl

INCLUDEPATH += $$PWD/root/include
#INCLUDEPATH += $$PWD/GnuWin32/include
#INCLUDEPATH += $$PWD/gsl/include

#DEPENDPATH += $$PWD/GnuWin32/lib
#DEPENDPATH += $$PWD/GnuWin32
#DEPENDPATH += $$PWD/GnuWin32/include

#DEPENDPATH += $$PWD/gsl/lib
#DEPENDPATH += $$PWD/gsl
#DEPENDPATH += $$PWD/gsl/include

DEPENDPATH += $$PWD/root/include


#win32:!win32-g++: PRE_TARGETDEPS += $$PWD/gsl/lib/gsl.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/gsl/lib/gsl.lib

#win32: PRE_TARGETDEPS += $$PWD/GnuWin32/lib/libgsl.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/GnuWin32/lib/libgsl.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libMatrix.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libMatrix.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libCore.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libCore.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libMathCore.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libMathCore.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGui.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGui.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libHist.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libHist.a

#win32: PRE_TARGETDEPS += $$PWD/root/lib/libCint.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libCint.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGraf.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGraf.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libRIO.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libRIO.a

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGPad.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGPad.a

#win32: LIBS += -L$$PWD/GnuWin32/lib/ -llibgsl
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/GnuWin32/lib/ -lgsld
#else:unix: LIBS += -L$$PWD/GnuWin32/lib/ -llibgsl

#INCLUDEPATH += $$PWD/GnuWin32/include
#DEPENDPATH += $$PWD/GnuWin32/include

#win32: PRE_TARGETDEPS += $$PWD/GnuWin32/lib/libgsl.lib
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/GnuWin32/lib/libgsld.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/GnuWin32/lib/gsl.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/GnuWin32/lib/gsld.lib
#else:unix: PRE_TARGETDEPS += $$PWD/GnuWin32/lib/libgsl.a

 #   PRECOMPILED_HEADER = pch/precompiled_header.h
 #   CONFIG += precompile_header
