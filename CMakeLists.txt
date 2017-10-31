CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

PROJECT(GoldEndotherial_IAEA_1)

FIND_PACKAGE(Geant4 REQUIRED ui_all vis_all)
IF(NOT Geant4_FOUND)
    MESSAGE(FATAL_ERROR "Geant4 not found!")
ENDIF(NOT Geant4_FOUND)

INCLUDE(${Geant4_USE_FILE})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
#  ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
FILE(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

ADD_EXECUTABLE(BYGNP BYGNP.cc ${sources} ${headers})
TARGET_LINK_LIBRARIES(BYGNP ${Geant4_LIBRARIES})

set(BYGNP_SCRIPTS
    init_vis.mac
    PSF_-0.05.IAEAphsp
    PSF_-0.05.IAEAheader
)

foreach(_script ${BYGNP_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY)
endforeach()




