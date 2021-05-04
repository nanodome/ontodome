TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
#		src/gasphase/gasphase.cpp \
#		src/gasphase/gasphasecv.cpp \
		src/main.cpp \
		src/base/relation.cpp \
		src/base/thing.cpp

HEADERS += \
		src/base/relation.h \
		src/base/thing.h \
#		src/gasphase/gasphase.h \
#		src/gasphase/gasphasecv.h \
		src/nanodome.h \
#		src/species/species.h
