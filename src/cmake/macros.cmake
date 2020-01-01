################################################################################
##
## BCL to FASTQ file converter
## Copyright (c) 2007-2015 Illumina, Inc.
##
## This software is covered by the accompanying EULA
## and certain third party copyright/licenses, and any user of this
## source file is bound by the terms therein.
##
################################################################################
##
## file macros.cmake
##
## CMake configuration file for common installation macros.
##
## author Roman Petrovski
## author Mauricio Varea
##
################################################################################

macro(configure_files srcDir destDir pattern)
    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files)

macro(install_files srcDir destDir pattern perm)
    file(GLOB templateFiles ${srcDir}/${pattern})
    file(INSTALL DESTINATION ${destDir} TYPE FILE
         FILES ${templateFiles} PERMISSIONS ${perm})
endmacro(install_files)

macro(configure_files_recursively srcDir destDir pattern)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files_recursively)

macro(install_files_recursively srcDir destDir pattern perm)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        get_filename_component(DIRNAME "${templateFile}" PATH)
        file(INSTALL DESTINATION ${destDir}/${DIRNAME} TYPE FILE
             FILES ${srcDir}/${templateFile} PERMISSIONS ${perm})
    endforeach(templateFile)
endmacro(install_files_recursively)

#   
# Macro to find libraries, with support for static-only search
#
macro(bcl2fastq_find_library name header library)
if    (NOT ${name}_INCLUDE_DIR)
    find_path(${name}_INCLUDE_DIR ${header} 
              HINTS ENV C_INCLUDE_PATH ENV CPATH ENV CPLUS_INCLUDE_PATH)
endif (NOT ${name}_INCLUDE_DIR)
if    (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
    find_library(${name}_LIBRARY 
                 NAMES "${BCL2FASTQ_LIBRARY_PREFIX}${library}${BCL2FASTQ_LIBRARY_SUFFIX}"
                 HINTS ENV LIBRARY_PATH)
endif (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
if(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
    set (HAVE_${name} ${${name}_LIBRARY})
    message (STATUS "Found ${name}  header: ${${name}_INCLUDE_DIR}/${header}")
    message (STATUS "Found ${name} library: ${${name}_LIBRARY}")
endif(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
endmacro(bcl2fastq_find_library)

#   
# Macro to find libraries, with support for static-only search
#
macro(bcl2fastq_find_header_or_die variable file)
find_file(${variable} ${file} HINTS ENV C_INCLUDE_PATH ENV CPATH ENV CPLUS_INCLUDE_PATH)
if    (${variable})
    message(STATUS "${file} found as ${${variable}}")
else  (${variable})
    message(FATAL_ERROR "Required header ${file} not found.")
endif (${variable})
endmacro(bcl2fastq_find_header_or_die)


