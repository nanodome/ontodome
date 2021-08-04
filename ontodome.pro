TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    #Tools and basic sources
    src/base/relation.cpp \
    src/base/thing.cpp \
    src/tools/clock.cpp \
    src/tools/utilities.cpp

HEADERS += \
    #Particle Phase general headers
    src/knowledgegenerators/nanomodels/particles/aggregate/aggregate.h \
    src/knowledgegenerators/nanomodels/particles/bond/bond.h \
    src/knowledgegenerators/nanomodels/particles/base/collisionalobject.h \
    src/knowledgegenerators/nanomodels/particles/bond/edge.h \
    src/knowledgegenerators/nanomodels/particles/aggregate/fractalaggregate.h \
    src/knowledgegenerators/nanomodels/particles/aggregate/rattleaggregate.h \
    src/knowledgegenerators/nanomodels/particles/aggregate/spatialaggregate.h \
    src/knowledgegenerators/nanomodels/particles/base/mesoobject.h \
    src/knowledgegenerators/nanomodels/particles/base/objectcounter.h \
    src/knowledgegenerators/nanomodels/particles/particle/particle.h \
    src/knowledgegenerators/nanomodels/particles/bond/particlebond.h \
    src/knowledgegenerators/nanomodels/particles/particlephase/particlephase.h \
    src/knowledgegenerators/nanomodels/particles/aggregate/pbmaggregate.h \
    src/knowledgegenerators/nanomodels/particles/particlephase/pbmparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/base/physicalobject.h \
    src/knowledgegenerators/nanomodels/particles/particlephase/cgmdparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/particlephase/dynamicparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/particlephase/randomwalkparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/spatial/cell.h \
    src/knowledgegenerators/nanomodels/particles/spatial/grid.h \
    src/knowledgegenerators/nanomodels/particles/tetra_constrainer.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_vertex.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_constrainer.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_constraint.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_link.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_bond.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_constants.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_tetrahedron.h \
    src/knowledgegenerators/nanomodels/particles/spatial/tetra_face.h \
    src/knowledgegenerators/nanomodels/particles/dynamic/dynamicparticle.h \
    src/knowledgegenerators/nanomodels/particles/dynamic/dynamicpoint.h \
    src/knowledgegenerators/nanomodels/particles/spatial/point.h \
    src/knowledgegenerators/nanomodels/particles/utilities/ndm_random.h \
    src/knowledgegenerators/nanomodels/particles/utilities/i_o.h \
    #PBM Headers
    src/knowledgegenerators/nanomodels/particles/PBM/pbmfractalparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/PBM/pbmbccaparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/PBM/pbmdlcaparticlephase.h \
    #CGMD Headers
    src/knowledgegenerators/nanomodels/particles/CGMD/rigidbodylangevinparticlephase.h \
    src/knowledgegenerators/nanomodels/particles/CGMD/constrainedlangevinparticlephase.h \
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

