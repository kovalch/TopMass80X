SET(Unfold_DIR "${${PROJECT_NAME}_SOURCE_DIR}/../unfolding")
SET( ${PROJECT_NAME}_ROOTDICT_DIR ${CMAKE_BINARY_DIR}/ROOTdicts )
SET( DICT_NAME ${${PROJECT_NAME}_ROOTDICT_DIR}/${PROJECT_NAME}_dict.cc )

SET( ${PROJECT_NAME}_HEADERS_THAT_NEED_ROOTDICT
  TopAnalysis
  diLeptonic_LinkDef )

SET( ${PROJECT_NAME}_LIBSOURCES
  HistoListReader
  TopAnalysis
  ttbarUtils
  analysisStructs
  AnalyzerBaseClass
  AnalyzerControlPlots
  plotterclass
  DilepSVDFunctions
  )

SET( Unfold_LIBSOURCES
  TopSVDFunctions
  TopSVDUnfold
  BaseSVDUnfold
  )

SET( ${PROJECT_NAME}_EXECUTABLES
  load_Analysis
  Histo
  TypicalErrors
  kinRecoEfficienciesAndSF
  createSystCovMatrix
  )


ADD_DEFINITIONS(-Wall -Wextra -pedantic -Wno-vla -O2 -g -std=c++0x -DSTANDALONE_COMPILE)


#####################
# include directories
#####################

INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${ZTopUtils_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ntupleCommon_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( "${Unfold_DIR}" )
INCLUDE_DIRECTORIES( BEFORE "${${PROJECT_NAME}_SOURCE_DIR}/include" ${${PROJECT_NAME}_ROOTDICT_DIR} )


####################
# targets
####################

MESSAGE (STATUS "starting ROOT dict header search:")
FOREACH( filename ${${PROJECT_NAME}_HEADERS_THAT_NEED_ROOTDICT} )
  FIND_FILE( ${filename}_fullPath "${filename}.h"
    PATHS "${${PROJECT_NAME}_SOURCE_DIR}/include"  NO_DEFAULT_PATH )
  MESSAGE (STATUS "  ${filename} found at: ${${filename}_fullPath}")
  SET(ROOT_DICT_HEADERS ${ROOT_DICT_HEADERS} ${${filename}_fullPath} )
  UNSET(${filename}_fullPath CACHE)
ENDFOREACH( filename ${${PROJECT_NAME}_HEADERS_THAT_NEED_ROOTDICT} )


FILE( MAKE_DIRECTORY ${${PROJECT_NAME}_ROOTDICT_DIR} )
ADD_CUSTOM_COMMAND( OUTPUT ${DICT_NAME}
                    COMMAND ${ROOT_CINT}
                    ARGS "-f" "${DICT_NAME}" "-c" "-p" "Math/LorentzVector.h" "Math/PtEtaPhiM4D.h" ${ROOT_DICT_HEADERS}
                    WORKING_DIRECTORY ${${PROJECT_NAME}_ROOTDICT_DIR}
                    DEPENDS ${ROOT_DICT_HEADERS}
                    COMMENT "generating dictionary ${DICT_NAME} from" ${ROOT_DICT_HEADERS})


## Unfortunately this explicit search for .C files is needed due to missing directory structure for code in unfolding/
MESSAGE (STATUS "starting search for Unfolding source files")
FOREACH( filename ${Unfold_LIBSOURCES} )
  FIND_FILE( ${filename}_cfilePath NAMES "${filename}.C"
    PATHS  ${Unfold_DIR}  NO_DEFAULT_PATH )
  MESSAGE (STATUS "  ${filename} found at: ${${filename}_cfilePath}")
  SET(Unfold_LIBSOURCES_full ${Unfold_LIBSOURCES_full} ${${filename}_cfilePath} )
  UNSET(${filename}_cfilePath CACHE)
ENDFOREACH( filename ${${PROJECT_NAME}_LIBSOURCES} )


ADD_LIBRARY( ${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIBSOURCES} ${Unfold_LIBSOURCES_full} ${DICT_NAME} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ZTopUtils_LIBRARIES} ${ntupleCommon_LIBRARIES} )


INSTALL(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION lib
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE  )


FOREACH( exec ${${PROJECT_NAME}_EXECUTABLES} )
  ADD_EXECUTABLE( ${exec} ${exec} )
  TARGET_LINK_LIBRARIES( ${exec} ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ZTopUtils_LIBRARIES} ${ntupleCommon_LIBRARIES} )
  INSTALL(TARGETS ${exec}
    RUNTIME DESTINATION bin )
ENDFOREACH( exec ${${PROJECT_NAME}_EXECUTABLES} )
