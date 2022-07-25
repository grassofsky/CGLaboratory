function(_find_module_libraries module_name use_default_path namelist_var postfix)
    set(lib_path_var ${module_name}_LIBRARIES${postfix})
    set(search_path_var CGLAB_LIB_DIRS${postfix})
    if (${use_default_path})
        find_library(${lib_path_var}
            NAMES ${${namelist_var}}
            PATHS ${${search_path_var}})
    else()
        find_library(${lib_path_var}
            NAMES ${${namelist_var}}
            PATHS ${${search_path_var}}
            NO_DEFAULT_PATH
            NO_CMAKE_FIND_ROOT_PATH)
    endif()
    mark_as_advanced(${lib_path_var})

    if ("${${lib_path_var}}" STREQUAL "${lib_path_var}-NOTFOUND")
        message(FATAL_ERROR "CGLaboratory Library ${_lib_MODULE}${postfix} not found in ${${search_path_var}}" )
    endif ()
endfunction()

#
# example:
# _cglab_find_module_libraries(
#   MODULE
#     GLEW 
#   NAMELIST
#     GLEW
#     glew
#     glew64
#     glew64d)
#
function(_cglab_find_module_libraries)
    cmake_parse_arguments(_lib
    "USE_DEFAULT_PATH"
    "NAMESPACE;MODULE"
    "NAMELIST"
    ${ARGN})

    if((NOT _lib_NAMESPACE) OR "${_lib_NAMESPACE}" STREQUAL "")
        set(_lib_NAMESPACE "VLibs")
    endif()
    set(_target_name "${_lib_NAMESPACE}::${_lib_MODULE}")
    if(WIN32)
        set(_target_property_name IMPORTED_IMPLIB)
    else()
        set(_target_property_name IMPORTED_LOCATION)
    endif()

    add_library(${_target_name} SHARED IMPORTED)

    get_property(_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(NOT _IS_MULTI_CONFIG)
        _find_module_libraries(
            ${_lib_MODULE}
            ${_lib_USE_DEFAULT_PATH}
            _lib_NAMELIST
            "")
        if(CGLAB_VERBOS)
            message(STATUS "FIND ${_lib_MODULE}: " ${${_lib_MODULE}_LIBRARIES})
        endif()
        set(LIBRARY_PATH_${_lib_MODULE}
            "${${_lib_MODULE}_LIBRARIES}" CACHE STRING "Path of ${_lib_MODULE}" FORCE)
        set_target_properties(${_target_name} PROPERTIES
            ${_target_property_name} "${${_lib_MODULE}_LIBRARIES}")
    else()
        _find_module_libraries(
            ${_lib_MODULE}
            ${_lib_USE_DEFAULT_PATH}
            _lib_NAMELIST
            "_RELEASE")
        _find_module_libraries(
            ${_lib_MODULE}
            ${_lib_USE_DEFAULT_PATH}
            _lib_NAMELIST
            "_DEBUG")

        if (CGLAB_VERBOS)
            message(STATUS "FIND ${_lib_MODULE} Debug: " ${${_lib_MODULE}_LIBRARIES_DEBUG})
            message(STATUS "FIND ${_lib_MODULE} Release: " ${${_lib_MODULE}_LIBRARIES_RELEASE})
        endif()
        set(LIBRARY_PATH_${_lib_MODULE}
            "$<$<CONFIG:Debug>:${${_lib_MODULE}_LIBRARIES_DEBUG}>"
            "$<$<CONFIG:Release>:${${_lib_MODULE}_LIBRARIES_RELEASE}>" CACHE STRING "Path of ${_lib_MODULE}" FORCE)

        set_target_properties(${_target_name} PROPERTIES
            ${_target_property_name}_DEBUG "${${_lib_MODULE}_LIBRARIES_DEBUG}"
            ${_target_property_name}_RELEASE "${${_lib_MODULE}_LIBRARIES_RELEASE}")
    endif()
    mark_as_advanced(LIBRARY_PATH_${_lib_MODULE})

    if (CGLAB_VERBOS)
        message(STATUS "LIBRARY_PATH_${_lib_MODULE}: " ${LIBRARY_PATH_${_lib_MODULE}})
    endif()
endfunction()
