
add_executable(roundtrip roundtrip.cpp)
target_link_libraries(roundtrip ${BABEL_LIBRARY})
install(TARGETS roundtrip
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)

set(tests
    uniqueid
    smiles
    stereo
    tetrahedral
    tetranonplanar
    tetraplanar
#    squareplanar
    cistrans
   )

foreach(test ${tests})
  add_executable(${test}test ${test}test.cpp)
  target_link_libraries(${test}test ${BABEL_LIBRARY})
  add_test(${test}_Test ${CMAKE_BINARY_DIR}/test/${test}test)
  set_tests_properties(${test}_Test PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
endforeach(test ${tests})
