file(GLOB_RECURSE PYXMOLPP_SOURCES
**.cpp
        **.h
)
pybind11_add_module(vcc ${PYXMOLPP_SOURCES})
target_include_directories(vcc PRIVATE ${CMAKE_CURRENT_LIST_DIR}/library)
target_include_directories(vcc PRIVATE ${MY_LIBRARIES})
target_link_libraries(vcc PRIVATE cpp_library)