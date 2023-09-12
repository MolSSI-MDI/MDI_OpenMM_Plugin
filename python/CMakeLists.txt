set(PACKAGE_NAME openmmmdi)
set(WRAP_FILE ${PACKAGE_NAME}/ExamplePluginWrapper.cpp)
set(MODULE_NAME openmmmdi)

# OPENMM_PYTHON_STAGING_DIR is a staging area for python, swig, and C files in the python package we are making.
set(OPENMM_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python"
    CACHE PATH "Temporary staging area for Python API wrappers")
mark_as_advanced(OPENMM_PYTHON_STAGING_DIR)

# Create package directory structure
file(MAKE_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/${PACKAGE_NAME})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}/mdisimulation.py ${OPENMM_PYTHON_STAGING_DIR}/${PACKAGE_NAME}/mdisimulation.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${PACKAGE_NAME}/__init__.py ${OPENMM_PYTHON_STAGING_DIR}/${PACKAGE_NAME}/__init__.py)

# Execute SWIG to generate source code for the Python module.

add_custom_command(
    OUTPUT "${WRAP_FILE}"
    COMMAND "${SWIG_EXECUTABLE}"
        -python -c++
        -o "${WRAP_FILE}"
        "-I${OPENMM_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/exampleplugin.i"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/exampleplugin.i"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Compile the Python module.

add_custom_target(PythonInstall DEPENDS "${WRAP_FILE}")
set(EXAMPLEPLUGIN_HEADER_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(EXAMPLEPLUGIN_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

set(OPENMM_PYTHON_USER_INSTALL OFF CACHE BOOL
    "Whether to install OpenMM Python binary module into the user site-packages directory")
mark_as_advanced(OPENMM_PYTHON_USER_INSTALL)

if(OPENMM_PYTHON_USER_INSTALL)
    # Install binary module to user location
add_custom_command(TARGET PythonInstall
    COMMAND "${PYTHON_EXECUTABLE}" setup.py build
    COMMAND "${PYTHON_EXECUTABLE}" setup.py install --user
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
else()
add_custom_command(TARGET PythonInstall
    COMMAND "${PYTHON_EXECUTABLE}" setup.py build
    COMMAND "${PYTHON_EXECUTABLE}" setup.py install
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
endif()
