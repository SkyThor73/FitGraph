QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++20
QMAKE_CXXFLAGS += -Wa,-mbig-obj

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    fadjusting.cpp \
    fcalculation.cpp \
    fdatapcs.cpp \
    felement.cpp \
    fexperimental.cpp \
    foutput.cpp \
    main.cpp \
    fitgraphmainwindow.cpp \
    qcustomplot.cpp

HEADERS += \
    fadjusting.h \
    fcalculation.h \
    fdatapcs.h \
    felement.h \
    fexperimental.h \
    fitgraphmainwindow.h \
    foutput.h \
    qcustomplot.h

FORMS += \
    fitgraphmainwindow.ui

TRANSLATIONS += \
    FitGraph_ru_RU.ts
CONFIG += lrelease
CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
