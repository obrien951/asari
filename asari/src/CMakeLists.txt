cmake_minimum_required(VERSION 2.8.12)
project(asariSpecToChrom)

find_package(pybind11 REQUIRED)
target_link_libraries(asariSpecToChrom z)
pybind11_add_module(asariSpecToChrom asariSpecToChrom.cpp specToChrom.cc)
