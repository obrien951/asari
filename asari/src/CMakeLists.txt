cmake_minimum_required(VERSION 2.8.12)
project(asariSpecToChrom)

find_package(pybind11 REQUIRED)
pybind11_add_module(asariSpecToChrom asariSpecToChrom.cpp specToChrom.cc spectrum.cc base64.cc )
target_link_libraries(asariSpecToChrom PRIVATE z)
