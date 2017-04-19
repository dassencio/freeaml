#!/bin/bash -e

function run_tests
{
    local use_openmp="$1"
    local build_dir="build-openmp-${use_openmp}"

    mkdir "${build_dir}"
    cd "${build_dir}"
    cmake -DUSE_OPENMP=${use_openmp} ..
    make -j4
    make test
    cd ..
    rm -rf "${build_dir}"
}

function build_documentation
{
    local build_dir="build-doc"

    mkdir "${build_dir}"
    cd "${build_dir}"
    cmake ..
    make doc
    cd ..
    rm -rf "${build_dir}"
}

if [ "$TRAVIS_OS_NAME" == "linux" ]
then
    run_tests ON
    run_tests OFF
    build_documentation
else
    run_tests OFF
fi
