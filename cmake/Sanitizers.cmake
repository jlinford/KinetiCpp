function(target_sanitizer_options target)

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # see https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html#silencing-unsigned-integer-overflow
        # Compile with -g and -fno-omit-frame-pointer to get proper debug information in your binary
        target_compile_options(${target} PRIVATE -g)
        target_compile_options(${target} PRIVATE -O2)
        target_compile_options(${target} PRIVATE -fno-omit-frame-pointer)
        target_compile_options(${target} PRIVATE -Werror -Wall -Wextra)

        target_compile_options(${target} PRIVATE -fsanitize=address)
        target_link_libraries(${target} PRIVATE -fsanitize=address)

        target_compile_options(${target} PRIVATE -fsanitize=undefined,float-divide-by-zero)
        target_link_libraries(${target} PRIVATE -fsanitize=undefined,float-divide-by-zero)

        target_compile_options(${target} PRIVATE -fsanitize=integer)
        target_link_libraries(${target} PRIVATE -fsanitize=integer)

        target_compile_options(${target} PRIVATE -fsanitize=nullability)
        target_link_libraries(${target} PRIVATE -fsanitize=nullability)

    endif()
    
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${target} PRIVATE -g)
        target_compile_options(${target} PRIVATE -O2)
        target_compile_options(${target} PRIVATE -fno-omit-frame-pointer)
        target_compile_options(${target} PRIVATE -Werror -Wall -Wextra)

        # need to use gold linker, otherwise travis gets '/usr/bin/ld: --push-state: unknown option' error
        target_link_libraries(${target} PRIVATE -fuse-ld=gold)

        if (NOT CMAKE_HOST_APPLE) 
            target_compile_options(${target} PRIVATE -fsanitize=undefined,float-divide-by-zero,float-cast-overflow)
            target_link_libraries(${target} PRIVATE -fsanitize=undefined,float-divide-by-zero,float-cast-overflow)

            target_compile_options(${target} PRIVATE -fsanitize=pointer-compare,pointer-subtract)
            target_link_libraries(${target} PRIVATE -fsanitize=pointer-compare,pointer-subtract)

            target_compile_options(${target} PRIVATE -fsanitize=address)
            target_link_libraries(${target} PRIVATE -fsanitize=address)
        endif()
    endif()

endfunction()
