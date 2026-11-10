cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0122 NEW)
project(mmkit VERSION 0.2.1 LANGUAGES CXX)
include(GNUInstallDirs)
include(CTest)
include(CMakePackageConfigHelpers)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
option(MMKIT_BUILD_PYTHON "Build Python bindings" OFF)
option(MMKIT_BUILD_JAVA "Build JNI bindings and JAR" OFF)
option(MMKIT_BUILD_CSHARP "Build C# P/Invoke bindings" OFF)
find_package(RDKit REQUIRED)
add_library(mmkit STATIC mmkit/Formula.cpp mmkit/Compound.cpp mmkit/Adduct.cpp)
add_library(mmkit::mmkit ALIAS mmkit)
target_compile_features(mmkit PUBLIC cxx_std_17)
set_target_properties(mmkit PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mmkit PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(mmkit PRIVATE RDKit::RDKit)
install(TARGETS mmkit EXPORT mmkitTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT cpp)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT cpp)
configure_package_config_file(cmake/mmkitConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/mmkitConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mmkit)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/mmkitConfigVersion.cmake COMPATIBILITY SameMajorVersion)
install(EXPORT mmkitTargets NAMESPACE mmkit:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mmkit COMPONENT cpp)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mmkitConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/mmkitConfigVersion.cmake
  cmake/FindRDKit.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mmkit COMPONENT cpp)

if(MMKIT_BUILD_PYTHON)
  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_native bindings/python/module.cpp)
  target_link_libraries(_native PRIVATE mmkit)
  set_target_properties(_native PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python/mmkit")
  foreach(module __init__ Formula Compound Adduct)
    configure_file(bindings/python/mmkit/${module}.py python/mmkit/${module}.py COPYONLY)
  endforeach()
  if(BUILD_TESTING)
    add_test(NAME python COMMAND ${Python_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_build_test.py
      ${CMAKE_CURRENT_BINARY_DIR}/python)
  endif()
  install(TARGETS _native LIBRARY DESTINATION mmkit COMPONENT python)
endif()

if(MMKIT_BUILD_JAVA OR MMKIT_BUILD_CSHARP)
  find_package(SWIG 4 REQUIRED)
  include(UseSWIG)
endif()
if(MMKIT_BUILD_JAVA)
  find_package(Java 11 REQUIRED COMPONENTS Development)
  find_package(JNI REQUIRED COMPONENTS JVM)
  include(UseJava)
  set(_java_dir ${CMAKE_CURRENT_BINARY_DIR}/java)
  file(MAKE_DIRECTORY ${_java_dir})
  configure_file(bindings/swig/mmkit.i ${CMAKE_CURRENT_BINARY_DIR}/mmkit_java.i COPYONLY)
  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/mmkit_java.i PROPERTIES CPLUSPLUS ON
    USE_TARGET_INCLUDE_DIRECTORIES TRUE SWIG_FLAGS "-package;org.mmkit")
  swig_add_library(mmkit_java TYPE SHARED LANGUAGE java OUTPUT_DIR ${_java_dir}
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/mmkit_java.i)
  target_include_directories(mmkit_java PRIVATE ${JNI_INCLUDE_DIRS})
  target_link_libraries(mmkit_java PRIVATE mmkit)
  set_target_properties(mmkit_java PROPERTIES OUTPUT_NAME mmkit_java)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mmkit-${PROJECT_VERSION}.jar
    COMMAND ${CMAKE_COMMAND} -DJAVA_DIR=${_java_dir} -DJAVAC=${Java_JAVAC_EXECUTABLE}
      -DJAR=${Java_JAR_EXECUTABLE} -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/mmkit-${PROJECT_VERSION}.jar
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildJar.cmake
    DEPENDS mmkit_java ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildJar.cmake)
  add_custom_target(mmkit_jar ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mmkit-${PROJECT_VERSION}.jar)
  install(TARGETS mmkit_java LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT java)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mmkit-${PROJECT_VERSION}.jar DESTINATION ${CMAKE_INSTALL_DATADIR}/java COMPONENT java)
  if(BUILD_TESTING)
    add_test(NAME java COMMAND ${CMAKE_COMMAND} -DJAVAC=${Java_JAVAC_EXECUTABLE} -DJAVA=${Java_JAVA_EXECUTABLE}
      -DJAR_FILE=${CMAKE_CURRENT_BINARY_DIR}/mmkit-${PROJECT_VERSION}.jar
      -DLIB_DIR=$<TARGET_FILE_DIR:mmkit_java> -DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}/tests/java/Smoke.java
      -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/java-test -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestJava.cmake)
  endif()
endif()
if(MMKIT_BUILD_CSHARP)
  set(_cs_dir ${CMAKE_CURRENT_BINARY_DIR}/csharp)
  file(MAKE_DIRECTORY ${_cs_dir})
  configure_file(bindings/swig/mmkit.i ${CMAKE_CURRENT_BINARY_DIR}/mmkit_csharp.i COPYONLY)
  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/mmkit_csharp.i PROPERTIES CPLUSPLUS ON
    USE_TARGET_INCLUDE_DIRECTORIES TRUE SWIG_FLAGS "-namespace;Mmkit;-dllimport;mmkit_csharp")
  swig_add_library(mmkit_csharp TYPE SHARED LANGUAGE csharp OUTPUT_DIR ${_cs_dir}
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/mmkit_csharp.i)
  target_link_libraries(mmkit_csharp PRIVATE mmkit)
  set_target_properties(mmkit_csharp PROPERTIES OUTPUT_NAME mmkit_csharp)
  find_program(DOTNET_EXECUTABLE dotnet REQUIRED)
  configure_file(bindings/csharp/Mmkit.csproj.in ${_cs_dir}/Mmkit.csproj @ONLY)
  add_custom_target(mmkit_nuget ALL
    COMMAND ${DOTNET_EXECUTABLE} pack ${_cs_dir}/Mmkit.csproj -c Release -o ${CMAKE_CURRENT_BINARY_DIR}/packages
      -p:NativeLibrary=$<TARGET_FILE:mmkit_csharp>
    DEPENDS mmkit_csharp)
  if(BUILD_TESTING)
    add_test(NAME csharp COMMAND ${CMAKE_COMMAND} -E env
      "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:mmkit_csharp>:$ENV{LD_LIBRARY_PATH}"
      ${DOTNET_EXECUTABLE} run --project ${CMAKE_CURRENT_SOURCE_DIR}/tests/csharp/Smoke.csproj
      -c Release -p:MmkitProject=${_cs_dir}/Mmkit.csproj -p:NativeLibrary=$<TARGET_FILE:mmkit_csharp>)
  endif()
endif()
if(BUILD_TESTING)
  add_executable(mmkit_test tests/cpp/test_formula.cpp)
  target_link_libraries(mmkit_test PRIVATE mmkit)
  add_test(NAME formula COMMAND mmkit_test)
  add_executable(mmkit_compound_test tests/cpp/test_compound.cpp)
  target_link_libraries(mmkit_compound_test PRIVATE mmkit)
  add_test(NAME compound COMMAND mmkit_compound_test)
  add_executable(mmkit_adduct_test tests/cpp/test_adduct.cpp)
  target_link_libraries(mmkit_adduct_test PRIVATE mmkit)
  add_test(NAME adduct COMMAND mmkit_adduct_test)
endif()
