cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11")

project(vtkdemo)


# Подключаем gmsh

set(GMSH_ROOT "/usr/local/opt/gmsh") 
set(GMSH_INCLUDE_DIR "${GMSH_ROOT}/include")
set(GMSH_LIBRARY "${GMSH_ROOT}/lib/libgmsh.dylib")


# Подключаем vtk
find_package(VTK COMPONENTS CommonCore CommonDataModel IOLegacy REQUIRED)


include_directories(/usr/local/Cellar/nlohmann-json/3.11.3/include)



add_executable(cubic2d cubic2d.cpp)
add_executable(cubic3d cubic3d.cpp)
add_executable(tetr3d tetr3d.cpp)
add_executable(Lab2 Lab2.cpp)


# Аналогично с VTK
set(VTK_DIR "/usr/local/Cellar/vtk/9.4.1_1")
find_package(VTK REQUIRED)
set(VTK_ROOT "/usr/local/Cellar/vtk/9.4.1_1") 
set(VTK_INCLUDE_DIRS "${VTK_ROOT}/include/vtk-9.4")
# ... подключим её со всеми её компонентами её нативными средствами
include(${VTK_USE_FILE})

# если этот вариант не сработал, собирайтесь из исходников
# можно пользоваться инструкциями от gmsh, все библиотеки собираются по единому принципу

# Укажем библиотеки для линковки
# ... тут только VTK
target_link_libraries(cubic2d PRIVATE ${VTK_LIBRARIES})
target_link_libraries(cubic3d PRIVATE ${VTK_LIBRARIES})
target_link_libraries(Lab2 PRIVATE ${VTK_LIBRARIES})

# ... а тут ещё и gmsh
#target_link_libraries(tetr3d ${VTK_LIBRARIES})
target_link_libraries(tetr3d PRIVATE ${VTK_LIBRARIES})
target_link_libraries(Lab2 PRIVATE ${VTK_LIBRARIES})


#target_link_libraries(tetr3d ${GMSH_LIBRARIES})



#gmsh директории

target_include_directories(cubic2d PRIVATE ${GMSH_INCLUDE_DIR}) #Вроде не надо но мало ли
target_include_directories(cubic3d PRIVATE ${GMSH_INCLUDE_DIR}) #Вроде не надо но мало ли
target_include_directories(tetr3d PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(Lab2 PRIVATE ${GMSH_INCLUDE_DIR})





target_link_libraries(cubic2d PRIVATE ${GMSH_LIBRARY}) #Вроде не надо но мало ли
target_link_libraries(cubic3d PRIVATE ${GMSH_LIBRARY}) #Вроде не надо но мало ли
target_link_libraries(tetr3d PRIVATE ${GMSH_LIBRARY})
target_link_libraries(Lab2 PRIVATE ${GMSH_LIBRARY})

#vtk директории

target_include_directories(cubic2d PRIVATE ${VTK_INCLUDE_DIRS})
target_include_directories(cubic3d PRIVATE ${VTK_INCLUDE_DIRS})
target_include_directories(tetr3d PRIVATE ${VTK_INCLUDE_DIRS})

target_include_directories(Lab2 PRIVATE ${VTK_INCLUDE_DIRS})




