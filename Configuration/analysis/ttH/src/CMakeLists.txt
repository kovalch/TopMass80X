SET(CommonAnalysis_DIR "${${PROJECT_NAME}_SOURCE_DIR}/../common")
SET( ${PROJECT_NAME}_ROOTDICT_DIR ${CMAKE_BINARY_DIR}/ROOTdicts )
SET( DICT_NAME ${${PROJECT_NAME}_ROOTDICT_DIR}/${PROJECT_NAME}_dict.cc )

SET( ${PROJECT_NAME}_HEADERS_THAT_NEED_ROOTDICT
  AnalysisBase
  HiggsAnalysis
  classes
  classesFwd
  LinkDef )

SET( ${PROJECT_NAME}_LIBSOURCES
  utils
  CommandLineParameters
  HistoListReader
  higgsUtils
  sampleHelpers
  analysisHelpers
  analysisStructs
  AnalysisBase
  HiggsAnalysis
  analysisUtils
  analysisObjectStructs
  KinReco
  KinematicReconstruction
  KinematicReconstruction_LSroutines
  KinematicReconstruction_MeanSol
  MvaTreeHandler
  MvaTreeAnalyzer
  MvaWeights2d
  MvaReader
  MvaFactory
  MvaValidation
  mvaSetup
  mvaStructs
  AnalysisHistograms
  BasicAnalyzer
  JetChargeAnalyzer
  JetMatchAnalyzer
  JetCategories
  ScaleFactors
  BTagUtils
  Playground
  DijetAnalyzer
  plotterHelpers
  plotterUtils
  Samples
  DyScaleFactors
  EventYields
  Plotter
  RootFileReader
  )

SET( ${PROJECT_NAME}_EXECUTABLES
  load_Analysis
  train_bdt_topSystemJetAssignment
  histo_bdt_topSystemJetAssignment
  histo_mvaVariables_topSystemJetAssignment
  Histo
  )


ADD_DEFINITIONS(-Wall -Wextra -pedantic -Wno-vla -O2 -g -std=c++0x -DSTANDALONE_COMPILE)


#####################
# include directories
#####################

INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${ZTopUtils_INCLUDE_DIRS} )
## TO-DO: use installed header files as soon as available
INCLUDE_DIRECTORIES( "${CommonAnalysis_DIR}/include")
## TO-DO: do proper separation in include and src directories
INCLUDE_DIRECTORIES( BEFORE "${${PROJECT_NAME}_SOURCE_DIR}/include" ${${PROJECT_NAME}_ROOTDICT_DIR} )


####################
# targets
####################

SET(FILE_SEARCH_PATHS "${${PROJECT_NAME}_SOURCE_DIR}/src" "${${PROJECT_NAME}_SOURCE_DIR}/include" "${CommonAnalysis_DIR}/src" "${CommonAnalysis_DIR}/include")

MESSAGE (STATUS "starting ROOT dict header search:")
FOREACH( filename ${${PROJECT_NAME}_HEADERS_THAT_NEED_ROOTDICT} )
  FIND_FILE( ${filename}_fullPath "${filename}.h"
    PATHS  ${FILE_SEARCH_PATHS}  NO_DEFAULT_PATH )
  MESSAGE (STATUS "  ${filename} --> ${${filename}_fullPath}")
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


MESSAGE (STATUS "starting search for source files")
FOREACH( filename ${${PROJECT_NAME}_LIBSOURCES} )
  FIND_FILE( ${filename}_cfilePath NAMES "${filename}.cc"
    PATHS  ${FILE_SEARCH_PATHS}  NO_DEFAULT_PATH )
  MESSAGE (STATUS "  ${filename} --> ${${filename}_cfilePath}")
  SET(${PROJECT_NAME}_LIBSOURCES_full ${${PROJECT_NAME}_LIBSOURCES_full} ${${filename}_cfilePath} )
  UNSET(${filename}_cfilePath CACHE)
ENDFOREACH( filename ${${PROJECT_NAME}_LIBSOURCES} )


ADD_LIBRARY( ${PROJECT_NAME} SHARED ${${PROJECT_NAME}_LIBSOURCES_full} ${DICT_NAME} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ZTopUtils_LIBRARIES} )

INSTALL(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION lib
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE  )


FOREACH( exec ${${PROJECT_NAME}_EXECUTABLES} )
  ADD_EXECUTABLE( ${exec} ${exec} )
  TARGET_LINK_LIBRARIES( ${exec} ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ZTopUtils_LIBRARIES} )
  INSTALL(TARGETS ${exec}
    RUNTIME DESTINATION bin )
ENDFOREACH( exec ${${PROJECT_NAME}_EXECUTABLES} )
