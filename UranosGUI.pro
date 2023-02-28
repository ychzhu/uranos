#-------------------------------------------------
#
# Project created by QtCreator 2016-06-13T21:45:52
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

win32:TARGET = UranosGUI
unix:TARGET = uranos
TEMPLATE = app

DESTDIR = build/
OBJECTS_DIR = build/.obj
MOC_DIR = build/.moc
RCC_DIR = build/.rcc
UI_DIR = build/.ui

CONFIG += console

#CONFIG   -= shared
#CONFIG   += static
#CONFIG   += staticlib

#QMAKE_LFLAGS		= -static -enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc

win32:RC_FILE = resources/urc.rc

INCLUDEPATH += src/

SOURCES += src/main.cpp\
    src/mainwindow.cpp \
    src/qcustomplot.cpp \
    src/Toolkit.cpp \
    src/dialogshowpic.cpp \
    src/customSplashScreen.cpp \
    src/visualizationenlarge.cpp \
    src/visualizationenlarge2.cpp


HEADERS += src/mainwindow.h \
    src/qcustomplot.h \
    src/Toolkit.h \
    src/dialogshowpic.h \
    src/customSplashScreen.h \
    src/visualizationenlarge.h \
    src/visualizationenlarge2.h

FORMS += ui/dialogshowpic.ui \
    ui/mainwindow.ui \
    ui/visualizationenlarge.ui \
    ui/visualizationEnlarge2.ui

resources.files = resources/about.png \
	resources/splashScreen.png \
	resources/uranos.ico \
	resources/uranos-logo-smallx2.png
resources.prefix = /
RESOURCES = resources

isEmpty(PREFIX) {
	PREFIX = /usr
}

target.path = $$PREFIX/bin
assets.path = $$PREFIX/share/uranos
assets.files = data/ENDFdata.zip data/IncomingSpectrum
INSTALLS += target assets


DEFINES += _CRT_SECURE_NO_WARNINGS
DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0

# CONFIG += static

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


win32: QMAKE_CXXFLAGS += -openmp
win32: QMAKE_LFLAGS += -fopenmp
win32: LIBS += -fopenmp
#unix: QMAKE_CXXFLAGS += -MM
win32: QMAKE_CXXFLAGS += -MP

#LIBS += -fopenmp

#LIBS += -LC:/Qt2/5.12.11/msvc2017/ -lQt5Core
LIBS += -L$$PWD/root/lib/
win32: LIBS += -llibMatrix -llibCore -llibMathCore -llibGui -llibHist -llibGraf -llibRIO -llibGPad
unix: LIBS += -lMatrix -lCore -lMathCore -lGui -lHist -lGraf -lRIO -lGpad -lThread -lImt -ltbb
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
