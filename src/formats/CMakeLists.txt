if(WITH_INCHI)
  add_definitions(-DINCHI_LINK_AS_DLL)
  add_subdirectory(inchi)
endif(WITH_INCHI)

set(formats
    APIInterface
    CSRformat
    PQSformat
    MCDLformat
    alchemyformat
    adfformat
    acrformat
    amberformat
    balstformat
    bgfformat
    boxformat
    cacaoformat
    cacheformat
    carformat
    cccformat
    chem3dformat
    chemdrawct
    chemtoolformat
    cifformat
    copyformat
    crkformat
    cssrformat
    dmolformat
    fastsearchformat
    fastaformat
    fchkformat
    featformat
    fhformat
    fingerprintformat
    freefracformat
    gamessformat
    gaussformat
    gausscubeformat
    gausszmatformat
    ghemicalformat
    gromos96format
    hinformat
    jaguarformat
    mdlformat
    mmcifformat
    mmodformat
    mol2format
    moldenformat
    molreport
    mopacformat
    mpdformat
    mpqcformat
    msiformat
    msmsformat
    nwchemformat
    opendxformat
    outformat
    pcmodelformat
    pdbformat
    pngformat
    povrayformat
    pqrformat
    qchemformat
    reportformat
    shelxformat
    smilesformat
    thermoformat
    tinkerformat
    titleformat
    turbomoleformat
    unichemformat
    viewmolformat
    xedformat
    xyzformat
    yasaraformat
    zindoformat
)

# Additional files for PDB and PQR formats
set(pdbformat_additional_sources atomrecord.cpp)
set(pqrformat_additional_sources atomrecord.cpp)

if(Boost_FOUND)
  set(formats
    ${formats}
    rxnformat
    chemdrawcdx
    chemkinformat
    rsmiformat
  )
endif(Boost_FOUND)

if(WITH_INCHI)
  set(formats
    ${formats}
    inchiformat
  )
  set(libs ${libs} inchi)
  set(inchiformat_additional_sources getinchi.cpp)
endif(WITH_INCHI)

if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(xml_formats
    cdxmlformat
    cmlformat
    pubchem
    xmlformat
  )
  if(Boost_FOUND)
    set(xml_formats
        ${xml_formats}
        cmlreactformat
    )
  endif(Boost_FOUND)
endif(LIBXML2_FOUND)

if(NOT WIN32)
  set(formats
    ${formats}
    xtcformat
  )
endif(NOT WIN32)


if(NOT WIN32)
  set(libs ${libs} m)
endif(NOT WIN32)

if(BUILD_SHARED)
  if(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_win32
        )
  else(WIN32)
    set(openbabel_srcs ${openbabel_srcs}
        dlhandler_unix
        )
  endif(WIN32)
endif(BUILD_SHARED)

foreach(format ${formats})
  add_library(${format} ${PLUGIN_TYPE} ${format}.cpp
    ${${format}_additional_sources}
    ${CMAKE_BINARY_DIR}/include/openbabel/babelconfig.h)
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
  install(TARGETS ${format}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION ${LIB}/openbabel/${BABEL_VERSION}
                  ARCHIVE DESTINATION ${LIB}/openbabel/${BABEL_VERSION})
  set_target_properties(${format} PROPERTIES
                        OUTPUT_NAME ${format}
                        PREFIX ""
                        SUFFIX ${MODULE_EXTENSION})
endforeach(format)

foreach(format ${xml_formats})
  add_library(${format} ${PLUGIN_TYPE} xml/${format}.cpp xml/xml.cpp)
  target_link_libraries(${format} ${libs} ${BABEL_LIBRARY})
  install(TARGETS ${format}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION ${LIB}/openbabel/${BABEL_VERSION}
                  ARCHIVE DESTINATION ${LIB}/openbabel/${BABEL_VERSION})
  set_target_properties(${format} PROPERTIES
                        OUTPUT_NAME ${format}
                        PREFIX ""
                        SUFFIX ${MODULE_EXTENSION})
endforeach(format)
