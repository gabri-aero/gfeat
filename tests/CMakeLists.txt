message(STATUS "Adding tests")

file(GLOB SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

foreach(SUBDIR ${SUBDIRS})
    if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR} AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/CMakeLists.txt)
        message(STATUS "Adding subdirectory: ${SUBDIR}")
        add_subdirectory(${SUBDIR})
    endif()
endforeach()