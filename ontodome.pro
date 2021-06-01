TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
	src/main.cpp \
	src/base/relation.cpp \
	src/base/thing.cpp \
	src/tools/clock.cpp \
	src/tools/utilities.cpp

HEADERS += \
	src/knowledgegenerators/saturationpressurematerialrelation.h \
	src/knowledgegenerators/saturationpressuremodels/saturationpressurepolynomialmodel.h \
	src/knowledgegenerators/surfacetensionmodels/surfacetensionpolynomialmodel.h \
	src/ontodome.h \
	src/tools/clock.h \
	src/tools/utilities.h \
	src/base/datatypes.h \
	src/base/baseclass.h \
	src/base/relation.h \
	src/base/species.h \
	src/base/thing.h \
	src/knowledgegenerators/knowledgegenerators.h \
	src/knowledgegenerators/surfacetensionmaterialrelation.h
#	src/knowledgegenerators/saturationpressurematerialrelation.h
	#src/knowledgegenerators/gasmodels/gasmodel.h \
	#src/knowledgegenerators/gasmodels/gasmodelcv.h \
	#src/knowledgegenerators/nanomodels/moments/momentmodelpratsinis.h \
	#src/knowledgegenerators/nanomodels/nucleation/cnt.h \
	#src/knowledgegenerators/statemodels/stateinterpolator.h \
