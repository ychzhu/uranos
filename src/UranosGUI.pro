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
    visualizationenlarge.cpp \
    visualizationenlarge2.cpp


HEADERS  += mainwindow.h \
    qcustomplot.h \
    Toolkit.h \
    dialogshowpic.h \
    customSplashScreen.h \
    visualizationenlarge.h \
    visualizationenlarge2.h

FORMS    += mainwindow.ui \
    dialogshowpic.ui \
    visualizationEnlarge2.ui \
    visualizationenlarge.ui

DEFINES += _CRT_SECURE_NO_WARNINGS
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0

CONFIG += static

win32:QMAKE_LFLAGS +=  /FORCE
win32:QMAKE_LFLAGS += /INCREMENTAL:NO

win32:CXXFLAGS += /O2
unix: QMAKE_CXXFLAGS += -Wno-sign-compare
unix: QMAKE_CXXFLAGS += -Wno-unused-parameter
unix: QMAKE_CXXFLAGS += -Wno-unused-variable
unix: QMAKE_CXXFLAGS += -Wno-unused-but-set-variable
unix: QMAKE_CXXFLAGS += -Wno-maybe-uninitialized
unix: QMAKE_CXXFLAGS += -Wno-unused-result
unix: QMAKE_CXXFLAGS += -Wno-unused-function
unix: QMAKE_CXXFLAGS += -Wno-misleading-indentation


#QMAKE_CXXFLAGS += -openmp
#unix: QMAKE_CXXFLAGS += -MM
win32: QMAKE_CXXFLAGS += -MP

#LIBS += -fopenmp

#LIBS += -LC:/Qt2/5.12.11/msvc2017/ -lQt5Core
win32: LIBS += -L$$PWD/root/lib/ -llibMatrix
win32: LIBS += -L$$PWD/root/lib/ -llibCore
win32: LIBS += -L$$PWD/root/lib/ -llibMathCore
win32: LIBS += -L$$PWD/root/lib/ -llibGui
win32: LIBS += -L$$PWD/root/lib/ -llibHist
win32: LIBS += -L$$PWD/root/lib/ -llibGraf
win32: LIBS += -L$$PWD/root/lib/ -llibRIO
win32: LIBS += -L$$PWD/root/lib/ -llibGPad

unix: LIBS += -L$$PWD/root/lib/ -lMatrix
unix: LIBS += -L$$PWD/root/lib/ -lCore
unix: LIBS += -L$$PWD/root/lib/ -lMathCore
unix: LIBS += -L$$PWD/root/lib/ -lGui
unix: LIBS += -L$$PWD/root/lib/ -lHist
unix: LIBS += -L$$PWD/root/lib/ -lGraf
unix: LIBS += -L$$PWD/root/lib/ -lRIO
unix: LIBS += -L$$PWD/root/lib/ -lGpad
unix: LIBS += -L$$PWD/root/lib/ -lThread
unix: LIBS += -L$$PWD/root/lib/ -lImt
unix: LIBS += -L$$PWD/root/lib/ -ltbb

unix: LIBS += -W


INCLUDEPATH += $$PWD/root/include

DEPENDPATH += $$PWD/root/include


win32: PRE_TARGETDEPS += $$PWD/root/lib/libMatrix.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libMatrix.a
#unix: PRE_TARGETDEPS += $$PWD/root/lib/libMatrix.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libCore.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libCore.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libCore.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libMathCore.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libMathCore.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libCore.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGui.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGui.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libGui.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libHist.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libHist.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libHist.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGraf.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGraf.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libGraf.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libRIO.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libRIO.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libRIO.so

win32: PRE_TARGETDEPS += $$PWD/root/lib/libGPad.lib
#else:win32-g++: PRE_TARGETDEPS += $$PWD/root/lib/libGPad.a
#else:unix: PRE_TARGETDEPS += $$PWD/root/lib/libGpad.so
