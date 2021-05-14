TEMPLATE = app
CONFIG += console c++20
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
		src/main.cpp \
		src/base/relation.cpp \
		src/base/thing.cpp \
		src/tools/clock.cpp \
#		src/tools/rkf45.cpp \
		src/tools/utilities.cpp

HEADERS += \
		src/base/relation.h \
		src/base/thing.h \
#		src/models/gasmodels/gasmodel.h \
        #src/models/gasmodels/gasmodelcv.h \
        #src/models/nanomodels/moments/momentmodelpratsinis.h \
        #src/models/nanomodels/nucleation/cnt.h \
        #src/models/statemodels/stateinterpolator.h \
		src/ontodome.h \
		src/tools/clock.h \
#		src/tools/rkf45.h \
		src/tools/utilities.h
