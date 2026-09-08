cmake_minimum_required(VERSION 3.16)
project(nnls-chroma VERSION 1.2.0 LANGUAGES C CXX)

# NNLS Chroma / Chordino / Tuning Vamp plugins, plus a command-line chord
# extractor and (under Emscripten) a WebAssembly module.
#
#   cmake -B build && cmake --build build && ctest --test-dir build
#
# Options:
#   NNLS_BUILD_PLUGIN   build the nnls-chroma Vamp plugin library (default ON)
#   NNLS_BUILD_CLI      build the chordino command-line tool (default ON)
#   NNLS_BUILD_TESTS    build and register the tests (default ON)
#   NNLS_VAMP_SDK_DIR   path to an unpacked Vamp plugin SDK source tree;
#                       when unset the SDK is fetched with FetchContent

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
if(MSVC)
  # M_PI for the Vamp host SDK's Window.h and the plugin sources.
  add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(NNLS_BUILD_PLUGIN "Build the nnls-chroma Vamp plugin library" ON)
option(NNLS_BUILD_CLI "Build the chordino command-line tool" ON)
option(NNLS_BUILD_TESTS "Build the tests" ON)
set(NNLS_VAMP_SDK_DIR "" CACHE PATH "Vamp plugin SDK source tree (fetched when empty)")

if(EMSCRIPTEN)
  set(NNLS_BUILD_PLUGIN OFF)
  set(NNLS_BUILD_TESTS OFF)
endif()

# ---------------------------------------------------------------- Vamp SDK

if(NNLS_VAMP_SDK_DIR)
  set(VAMP_SDK_SOURCE_DIR "${NNLS_VAMP_SDK_DIR}")
else()
  include(FetchContent)
  FetchContent_Declare(vamp_plugin_sdk
    GIT_REPOSITORY https://github.com/vamp-plugins/vamp-plugin-sdk.git
    GIT_TAG        vamp-plugin-sdk-v2.10
    GIT_SHALLOW    TRUE)
  FetchContent_GetProperties(vamp_plugin_sdk)
  if(NOT vamp_plugin_sdk_POPULATED)
    FetchContent_Populate(vamp_plugin_sdk)
  endif()
  set(VAMP_SDK_SOURCE_DIR "${vamp_plugin_sdk_SOURCE_DIR}")
endif()
message(STATUS "Vamp SDK: ${VAMP_SDK_SOURCE_DIR}")

set(VAMP_SDK_SRC "${VAMP_SDK_SOURCE_DIR}/src/vamp-sdk")
set(VAMP_HOSTSDK_SRC "${VAMP_SDK_SOURCE_DIR}/src/vamp-hostsdk")

# Plugin-side SDK (what the plugin library links).
add_library(vamp-sdk STATIC
  "${VAMP_SDK_SRC}/PluginAdapter.cpp"
  "${VAMP_SDK_SRC}/RealTime.cpp"
  "${VAMP_SDK_SRC}/FFT.cpp")
target_include_directories(vamp-sdk PUBLIC "${VAMP_SDK_SOURCE_DIR}")

# Host-side SDK, split so the WebAssembly build can skip dynamic loading.
add_library(vamp-hostsdk-adapters STATIC
  "${VAMP_HOSTSDK_SRC}/RealTime.cpp"
  "${VAMP_HOSTSDK_SRC}/PluginWrapper.cpp"
  "${VAMP_HOSTSDK_SRC}/PluginBufferingAdapter.cpp"
  "${VAMP_HOSTSDK_SRC}/PluginChannelAdapter.cpp"
  "${VAMP_HOSTSDK_SRC}/PluginInputDomainAdapter.cpp"
  "${VAMP_HOSTSDK_SRC}/PluginSummarisingAdapter.cpp")
target_include_directories(vamp-hostsdk-adapters PUBLIC "${VAMP_SDK_SOURCE_DIR}")
target_compile_definitions(vamp-hostsdk-adapters PUBLIC _VAMP_PLUGIN_IN_HOST_NAMESPACE=1)

if(NOT EMSCRIPTEN)
  add_library(vamp-hostsdk STATIC
    "${VAMP_HOSTSDK_SRC}/Files.cpp"
    "${VAMP_HOSTSDK_SRC}/PluginHostAdapter.cpp"
    "${VAMP_HOSTSDK_SRC}/PluginLoader.cpp"
    "${VAMP_HOSTSDK_SRC}/host-c.cpp")
  target_link_libraries(vamp-hostsdk PUBLIC vamp-hostsdk-adapters ${CMAKE_DL_LIBS})
endif()

# ---------------------------------------------------------------- plugin core

set(NNLS_CORE_SOURCES
  NNLSBase.cpp NNLSChroma.cpp Chordino.cpp Tuning.cpp
  chromamethods.cpp viterbi.cpp nnls.c)

function(nnls_core_flags target)
  target_compile_definitions(${target} PRIVATE
    $<$<C_COMPILER_ID:MSVC>:_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS>
    $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS>)
  target_compile_options(${target} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>)
endfunction()

if(NNLS_BUILD_PLUGIN)
  add_library(nnls-chroma MODULE ${NNLS_CORE_SOURCES} plugins.cpp)
  nnls_core_flags(nnls-chroma)
  target_link_libraries(nnls-chroma PRIVATE vamp-sdk)
  set_target_properties(nnls-chroma PROPERTIES
    PREFIX ""
    OUTPUT_NAME "nnls-chroma"
    CXX_VISIBILITY_PRESET default
    C_VISIBILITY_PRESET default
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/plugin")
  foreach(cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(nnls-chroma PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_${cfg} "${CMAKE_BINARY_DIR}/plugin"
      RUNTIME_OUTPUT_DIRECTORY_${cfg} "${CMAKE_BINARY_DIR}/plugin")
  endforeach()
  if(APPLE)
    set_target_properties(nnls-chroma PROPERTIES SUFFIX ".dylib")
    target_link_options(nnls-chroma PRIVATE
      "-Wl,-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.list")
  elseif(WIN32)
    target_sources(nnls-chroma PRIVATE nnls-chroma.def)
    if(MINGW)
      target_link_options(nnls-chroma PRIVATE -static -static-libgcc -static-libstdc++)
    endif()
  else()
    target_link_options(nnls-chroma PRIVATE
      "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.map")
  endif()
  # Files hosts expect next to the library.
  foreach(f chord.dict nnls-chroma.n3 nnls-chroma.cat)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
      configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${f}" "${CMAKE_BINARY_DIR}/plugin/${f}" COPYONLY)
    endif()
  endforeach()
endif()

# The same plugin classes compiled for in-process use by a host (the CLI,
# the tests and the WebAssembly module).
add_library(chordino-core STATIC ${NNLS_CORE_SOURCES} cli/analyze.cpp)
nnls_core_flags(chordino-core)
target_include_directories(chordino-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(chordino-core PUBLIC vamp-hostsdk-adapters)

# ---------------------------------------------------------------- CLI

if(NNLS_BUILD_CLI AND NOT EMSCRIPTEN)
  add_library(chordino-audiofile STATIC cli/audiofile.cpp)
  target_include_directories(chordino-audiofile PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  target_compile_definitions(chordino-audiofile PRIVATE
    $<$<C_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)

  add_executable(chordino cli/chordino.cpp)
  target_compile_definitions(chordino PRIVATE CHORDINO_VERSION="${PROJECT_VERSION}")
  target_link_libraries(chordino PRIVATE chordino-core chordino-audiofile)
  if(MINGW)
    target_link_options(chordino PRIVATE -static -static-libgcc -static-libstdc++)
  endif()
  install(TARGETS chordino RUNTIME DESTINATION bin)
endif()

# ---------------------------------------------------------------- WebAssembly

if(EMSCRIPTEN)
  add_executable(chordino-wasm wasm/chordino_wasm.cpp)
  target_link_libraries(chordino-wasm PRIVATE chordino-core)
  set_target_properties(chordino-wasm PROPERTIES OUTPUT_NAME "chordino" SUFFIX ".js"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/wasm")
  target_link_options(chordino-wasm PRIVATE
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createChordino -sENVIRONMENT=web,worker,node
    -sALLOW_MEMORY_GROWTH=1 -sFILESYSTEM=0
    "-sEXPORTED_FUNCTIONS=['_chordino_analyze_json','_chordino_free','_malloc','_free']"
    "-sEXPORTED_RUNTIME_METHODS=['cwrap','UTF8ToString','HEAPF32']")
  configure_file(wasm/demo.html "${CMAKE_BINARY_DIR}/wasm/index.html" COPYONLY)
  configure_file(wasm/chordino-wasm.js "${CMAKE_BINARY_DIR}/wasm/chordino-wasm.js" COPYONLY)
  configure_file(wasm/test.mjs "${CMAKE_BINARY_DIR}/wasm/test.mjs" COPYONLY)
endif()

# ---------------------------------------------------------------- tests

if(NNLS_BUILD_TESTS)
  enable_testing()
  add_executable(test_chords tests/test_chords.cpp tests/synth.cpp)
  target_link_libraries(test_chords PRIVATE chordino-core chordino-audiofile)
  add_test(NAME chords COMMAND test_chords "${CMAKE_BINARY_DIR}/progression.wav")

  if(NNLS_BUILD_PLUGIN)
    add_executable(test_plugin tests/test_plugin.cpp tests/synth.cpp)
    target_link_libraries(test_plugin PRIVATE chordino-core vamp-hostsdk)
    add_test(NAME plugin COMMAND test_plugin)
    set_tests_properties(plugin PROPERTIES
      ENVIRONMENT "VAMP_PATH=${CMAKE_BINARY_DIR}/plugin")
  endif()

  if(NNLS_BUILD_CLI)
    add_test(NAME cli-text COMMAND chordino "${CMAKE_BINARY_DIR}/progression.wav")
    set_tests_properties(cli-text PROPERTIES DEPENDS chords
      PASS_REGULAR_EXPRESSION "C\n.*F\n.*G\n.*Am\n")
    add_test(NAME cli-json COMMAND chordino --format json "${CMAKE_BINARY_DIR}/progression.wav")
    set_tests_properties(cli-json PROPERTIES DEPENDS chords
      PASS_REGULAR_EXPRESSION "\"chords\": \\[.*\"label\": \"Am\"")
    add_test(NAME cli-lab COMMAND chordino --format lab "${CMAKE_BINARY_DIR}/progression.wav")
    set_tests_properties(cli-lab PROPERTIES DEPENDS chords
      PASS_REGULAR_EXPRESSION "\tC\n.*\tF\n.*\tG\n.*\tAm\n")
    add_test(NAME cli-chroma COMMAND chordino --chroma "${CMAKE_BINARY_DIR}/progression.wav")
    set_tests_properties(cli-chroma PROPERTIES DEPENDS chords
      PASS_REGULAR_EXPRESSION "time,A,Bb,B,C")
    add_test(NAME cli-tuning COMMAND chordino --tuning "${CMAKE_BINARY_DIR}/progression.wav")
    set_tests_properties(cli-tuning PROPERTIES DEPENDS chords
      PASS_REGULAR_EXPRESSION "(43[89]|44[0-2])\\.[0-9] Hz")
    add_test(NAME cli-version COMMAND chordino --version)
    set_tests_properties(cli-version PROPERTIES PASS_REGULAR_EXPRESSION "chordino ${PROJECT_VERSION}")
  endif()
endif()

# ---------------------------------------------------------------- install

if(NNLS_BUILD_PLUGIN)
  if(APPLE)
    set(NNLS_PLUGIN_INSTALL_DIR "$ENV{HOME}/Library/Audio/Plug-Ins/Vamp" CACHE PATH "Vamp plugin directory")
  elseif(WIN32)
    set(NNLS_PLUGIN_INSTALL_DIR "$ENV{ProgramFiles}/Vamp Plugins" CACHE PATH "Vamp plugin directory")
  else()
    set(NNLS_PLUGIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/vamp" CACHE PATH "Vamp plugin directory")
  endif()
  install(TARGETS nnls-chroma LIBRARY DESTINATION "${NNLS_PLUGIN_INSTALL_DIR}")
  install(FILES chord.dict nnls-chroma.n3 nnls-chroma.cat DESTINATION "${NNLS_PLUGIN_INSTALL_DIR}" OPTIONAL)
endif()
