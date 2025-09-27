cmake_minimum_required(VERSION 3.10)

# Define the project name and version
project(EFEM VERSION 1.0 LANGUAGES CXX)


# Optionally, you can specify C++ standard



# Add source files
 add_executable(nfem src/main.cpp 
                     src/efem.cpp
                     src/MeshEngine.cpp )


 target_compile_features(nfem PUBLIC cxx_std_17)
 target_include_directories(nfem PRIVATE  "include" "C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2/include")

target_link_libraries(nfem
    PRIVATE
        "C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2/lib/mkl_rt.lib"
        
)