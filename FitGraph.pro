QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++20

mingw: QMAKE_CXXFLAGS += -Wa,-mbig-obj
msvc:  QMAKE_CXXFLAGS += /bigobj

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += \
    $$PWD/src \
    $$PWD/3rdparty/qcustomplot \
    $$PWD/3rdparty/Eigen

SOURCES += \
    src/fadjusting.cpp \
    src/fcalculation.cpp \
    src/fdatapcs.cpp \
    src/felement.cpp \
    src/fexperimental.cpp \
    src/foutput.cpp \
    src/main.cpp \
    src/fitgraphmainwindow.cpp \
    3rdparty/qcustomplot/qcustomplot.cpp

HEADERS += \
    src/fadjusting.h \
    src/fcalculation.h \
    src/fdatapcs.h \
    src/felement.h \
    src/fexperimental.h \
    src/fitgraphmainwindow.h \
    src/foutput.h \
    3rdparty/qcustomplot/qcustomplot.h

FORMS += \
    ui/fitgraphmainwindow.ui

TRANSLATIONS += \
    translations/FitGraph_ru_RU.ts
CONFIG += lrelease
CONFIG += embed_translations

RESOURCES += \
    resources.qrc

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target