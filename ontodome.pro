TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    #PBM sources
#    src/knowledgegenerators/nanomodels/pbm/aggregate.cpp \
#    src/knowledgegenerators/nanomodels/pbm/particle.cpp \
#    src/knowledgegenerators/nanomodels/pbm/particlebond.cpp \
#    src/knowledgegenerators/nanomodels/pbm/particlephase.cpp \
#    src/knowledgegenerators/nanomodels/pbm/pbmaggregate.cpp \
#    src/knowledgegenerators/nanomodels/pbm/pbmfractalparticlephase.cpp \
#    src/knowledgegenerators/nanomodels/pbm/pbmparticlephase.cpp \
    #Tools and basic sources
    src/base/relation.cpp \
    src/base/thing.cpp \
    src/tools/clock.cpp \
    src/tools/utilities.cpp

HEADERS += \
    #PBM headers
    src/knowledgegenerators/nanomodels/pbm/aggregate.h \
    src/knowledgegenerators/nanomodels/pbm/bond.h \
    src/knowledgegenerators/nanomodels/pbm/collisionalobject.h \
    src/knowledgegenerators/nanomodels/pbm/edge.h \
    src/knowledgegenerators/nanomodels/pbm/fractalaggregate.h \
    src/knowledgegenerators/nanomodels/pbm/mesoobject.h \
    src/knowledgegenerators/nanomodels/pbm/objectcounter.h \
    src/knowledgegenerators/nanomodels/pbm/particle.h \
    src/knowledgegenerators/nanomodels/pbm/particlebond.h \
    src/knowledgegenerators/nanomodels/pbm/particlephase.h \
    src/knowledgegenerators/nanomodels/pbm/pbmaggregate.h \
    src/knowledgegenerators/nanomodels/pbm/pbmfractalparticlephase.h \
    src/knowledgegenerators/nanomodels/pbm/pbmparticlephase.h \
    src/knowledgegenerators/nanomodels/pbm/physicalobject.h \
    #Gas models headers
    src/knowledgegenerators/gasmodels/gasmodels.h \
    src/knowledgegenerators/gasmodels/gasmodel.h \
    #Moment models headers
    src/knowledgegenerators/nanomodels/moments/moments.h \
    src/knowledgegenerators/nanomodels/moments/momentmodelpratsinis.h \
    #Nucleation models headers
    src/knowledgegenerators/nanomodels/nucleation/classicalnucleationtheory.h \
    src/knowledgegenerators/nanomodels/nucleation/nucleation.h \
    src/knowledgegenerators/nanomodels/nucleation/cnt.h \
    #Material relations headers
    src/knowledgegenerators/surfacetensionmaterialrelation.h \
    src/knowledgegenerators/saturationpressurematerialrelation.h \
    src/knowledgegenerators/saturationpressuremodels/saturationpressurematerialrelation.h \
    src/knowledgegenerators/saturationpressuremodels/saturationpressurepolynomialsoftwaremodel.h \
    src/knowledgegenerators/surfacetensionmodels/surfacetensionmaterialrelation.h \
    src/knowledgegenerators/surfacetensionmodels/surfacetensionpolynomialsoftwaremodel.h \
    #Tools and basic headers
    src/ontodome.h \
    src/base/ontodome_base.h \
    src/tools/clock.h \
    src/tools/utilities.h \
    src/base/datatypes.h \
    src/base/baseclass.h \
    src/base/relation.h \
    src/base/species.h \
    src/base/thing.h \
    src/knowledgegenerators/knowledgegenerators.h
    #src/knowledgegenerators/statemodels/stateinterpolator.h \

