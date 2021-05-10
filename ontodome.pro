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
		src/models/gasmodel.h \
		src/models/gasmodelcv.h \
		src/ontodome.h
