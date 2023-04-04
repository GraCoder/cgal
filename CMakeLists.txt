cmake_minimum_required(VERSION 3.8.0)
set(prj_name cgal)
project(cgal)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(Qt5 COMPONENTS Widgets Gui Network REQUIRED)
#set(moc_h
#)
#qt5_wrap_cpp(mocfiles ${moc_h})
#source_group("moc" FILES ${mocfiles})
#qt5_add_resources(mocfiles)

find_package(CGAL)

set(hdr
  Polygon_2_Clipper.h
)

set(src
	main.cpp
)

add_executable(${prj_name} WIN32 ${hdr} ${src} ${mocfiles})

set_target_properties(${prj_name} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${prj_name}>"
                                               VS_DEBUGGER_COMMAND           "$<TARGET_FILE:${prj_name}>"
                                               VS_DEBUGGER_ENVIRONMENT       "PATH=%PATH%;${CMAKE_PREFIX_PATH}/bin")

target_link_libraries(${prj_name} PRIVATE CGAL::CGAL)

#target_link_libraries(${prj_name} Qt5::Widgets Qt5::Gui Qt5::Network)
#target_compile_definitions(${prj_name} PUBLIC QAPPLICATION_CLASS=QApplication)

