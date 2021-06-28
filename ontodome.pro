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
	src/knowledgegenerators/gasmodels/gasmodels.h \
	src/knowledgegenerators/nanomodels/moments/moments.h \
	src/knowledgegenerators/nanomodels/nucleation/classicalnucleationtheory.h \
	src/knowledgegenerators/nanomodels/nucleation/nucleation.h \
	src/knowledgegenerators/saturationpressuremodels/saturationpressurematerialrelation.h \
	src/knowledgegenerators/saturationpressuremodels/saturationpressurepolynomialsoftwaremodel.h \
	src/knowledgegenerators/surfacetensionmodels/surfacetensionmaterialrelation.h \
	src/knowledgegenerators/surfacetensionmodels/surfacetensionpolynomialsoftwaremodel.h \
	src/ontodome.h \
	src/tools/clock.h \
	src/tools/utilities.h \
	src/base/datatypes.h \
	src/base/baseclass.h \
	src/base/relation.h \
	src/base/species.h \
	src/base/thing.h \
	src/knowledgegenerators/knowledgegenerators.h \
	src/knowledgegenerators/surfacetensionmaterialrelation.h \
	src/knowledgegenerators/saturationpressurematerialrelation.h \
	src/knowledgegenerators/nanomodels/moments/momentmodelpratsinis.h \
	src/knowledgegenerators/nanomodels/nucleation/cnt.h \
	#src/knowledgegenerators/statemodels/stateinterpolator.h \
	src/knowledgegenerators/gasmodels/gasmodel.h
