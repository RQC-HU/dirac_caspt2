message(STATUS "Getting git hash")
message(STATUS "commit_hash.cmake PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message(STATUS "commit_hash.cmake EXECUTABLE_OUTPUT_PATH: ${EXECUTABLE_OUTPUT_PATH}")
message(STATUS "commit_hash.cmake CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
set(git_hash "unknown") # Default value

if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    file(STRINGS "${PROJECT_SOURCE_DIR}/.git/HEAD" git_head)

    if(git_head MATCHES "ref: ")
        string(REGEX REPLACE "ref: " "" git_head "${git_head}")
        file(STRINGS "${PROJECT_SOURCE_DIR}/.git/${git_head}" git_hash LIMIT_COUNT 1)
    else()
        set(git_hash "${git_head}")
    endif()
else()
    set(git_hash "unknown")
endif()

file(WRITE ${EXECUTABLE_OUTPUT_PATH}/.commit_hash "${git_hash}")
