TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    globals.cpp \
    Main.cpp

HEADERS += \
    globals.h
#QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11

#LIBS+=  /usr/lib/paraview/libvtkCommonCore.so  /usr/lib/paraview/libvtkFiltersCore.so
#LIBS+=  /usr/lib/paraview/libvtkIOImage.so /usr/lib/paraview/libvtkIOCore.so /usr/lib/paraview/libvtkIOXMLParser.so
#LIBS+=  /usr/lib/paraview/libvtkCommonDataModel.so /usr/lib/paraview/libvtkIOXML.so /usr/lib/paraview/libvtkIOCore.so
#LIBS+=  /usr/lib/paraview/libvtkIOExport.so /usr/lib/paraview/libvtkIOLegacy.so
INCLUDEPATH += /usr/include/paraview
#INCLUDEPATH += /usr/lib/paraview

LIBS+=  -lGL -lGLU -lglut -lm
LIBS += -L/usr/lib/paraview -lvtkCommonCore -lvtkFiltersCore -lvtkIOImage -lvtkIOCore -lvtkIOXMLParser -lvtkCommonDataModel -lvtkIOXML -lvtkIOCore -lvtkIOExport -lvtkIOLegacy

#LIBS += -lopenGL32 -lGLU32 -lm
#LIBS += -L$$PWD/my_lib -lglut32

