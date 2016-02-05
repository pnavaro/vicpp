#-------------------------------------------------
#
# Project created by QtCreator 2011-08-19T11:27:29
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = vicpp
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    xdmf.cpp \
    grid.cpp \
    poisson.cpp \
    tsc.cpp \
    vortexring.cpp

HEADERS += \
    xdmf.hpp \
    grid.hpp \
    poisson.hpp \
    tsc.hpp \
    vortexring.hpp


macx: LIBS += -L$$PWD/../../../../../../usr/local/lib/ -lhdf5

INCLUDEPATH += $$PWD/../../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../../usr/local/include

macx: PRE_TARGETDEPS += $$PWD/../../../../../../usr/local/lib/libhdf5.a

macx: LIBS += -L$$PWD/../../../../../../usr/local/lib/ -lfftw3

INCLUDEPATH += $$PWD/../../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../../usr/local/include

macx: PRE_TARGETDEPS += $$PWD/../../../../../../usr/local/lib/libfftw3.a
