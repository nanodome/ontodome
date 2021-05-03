TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
		src/main.cpp \
		src/base/relation.cpp \
		src/base/thing.cpp

HEADERS += \
		src/base/relation.h \
		src/base/thing.h
