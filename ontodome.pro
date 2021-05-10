TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
		src/main.cpp \
		src/base/relation.cpp \
		src/base/thing.cpp

HEADERS += \
		src/base/relation.h \
		src/base/thing.h \
		src/models/gasmodels/gasmodel.h \
		src/models/gasmodels/gasmodelcv.h \
		src/models/nanomodels/nucleation/cnt.h \
		src/ontodome.h
