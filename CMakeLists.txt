cmake_minimum_required(VERSION 3.1)
project(Triangulation)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_PREFIX_PATH "/usr/lib/x86_64-linux-gnu/qt5/")
set(CMAKE_UIC_EXECUTABLE "/usr/bin/uic")
set(Qt5_DIR "/usr/lib/x86_64-linux-gnu/cmake/Qt5")

# Εύρεση της βιβλιοθήκης CGAL & Qt5
find_package(CGAL REQUIRED COMPONENTS Core Qt5)


# Ενεργοποίηση των warnings για debug 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -Wl,-rpath,/usr/lib/x86_64-linux-gnu")

add_definitions(-DCGAL_USE_BASIC_VIEWER)
find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)

# Εύρεση της βιβλιοθήκης Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS})

# Εύρεση των απαιτούμενων βιβλιοθηκών
include_directories(${CGAL_INCLUDE_DIRS})
link_directories( Qt5::Widgets ${CGAL_LIBRARY_DIRS})

# Προσθήκη των πηγών
add_executable(triangulation src/main.cpp src/triangulation.cpp src/strategies.cpp src/old_triangulation.cpp)

# Προσθήκη του include φακέλου για το json.hpp
target_include_directories(triangulation PRIVATE include)
# Σύνδεση της CGAL και Boost βιβλιοθήκης
target_link_libraries(triangulation CGAL::CGAL ${Boost_LIBRARIES} CGAL::CGAL_Qt5 Qt5::Widgets Qt5::OpenGL)