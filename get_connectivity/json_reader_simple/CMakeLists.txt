cmake_minimum_required(VERSION 3.21..3.26)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0083 NEW)
cmake_policy(SET CMP0104 NEW)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(tallarn LANGUAGES CXX)

# Set default install prefix to the project source directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "..." FORCE)
endif()

# nlohmann_json
if(USE_EXTERNAL_JSON)
  find_package(nlohmann_json 3.11.3 REQUIRED)
else()
  include(FetchContent)
  FetchContent_Declare(
    json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
        SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/json")
  FetchContent_MakeAvailable(json)
endif()


# Add subdirectories for each language component

add_executable(tallarn ${CMAKE_SOURCE_DIR}/source/main.cpp) 
target_include_directories(tallarn PUBLIC ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(tallarn PUBLIC nlohmann_json::nlohmann_json)


include(FeatureSummary)
feature_summary(WHAT ALL)

install(TARGETS tallarn RUNTIME DESTINATION bin)
