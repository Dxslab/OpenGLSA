#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
set -euo pipefail

if [[ $# -ne 1 || ( "${1:-}" != intel && "${1:-}" != gnu ) ]]; then
    echo 'Usage: bash install_deps.sh {intel|gnu}' >&2
    exit 2
fi

case "$1" in
    intel)
        # Honor explicit choices and fill missing compilers from the same suite.
        intel_suite=""
        for compiler in "${FC:-}" "${CC:-}" "${CXX:-}"; do
            case "${compiler##*/}" in
                ifort|icc|icpc) suite=classic ;;
                ifx|icx|icpx) suite=llvm ;;
                *) continue ;;
            esac
            if [[ -n "$intel_suite" && "$intel_suite" != "$suite" ]]; then
                echo 'Choose FC, CC and CXX from the same Intel compiler suite.' >&2
                exit 1
            fi
            intel_suite="$suite"
        done
        if [[ -z "$intel_suite" ]]; then
            if command -v ifx >/dev/null && command -v icx >/dev/null && command -v icpx >/dev/null; then
                intel_suite=llvm
            elif command -v ifort >/dev/null && command -v icc >/dev/null && command -v icpc >/dev/null; then
                intel_suite=classic
            elif [[ -z "${FC:-}" || -z "${CC:-}" || -z "${CXX:-}" ]]; then
                echo 'Load a complete Intel compiler suite or set FC, CC and CXX explicitly.' >&2
                exit 1
            fi
        fi
        if [[ "$intel_suite" == classic ]]; then
            export FC="${FC:-ifort}" CC="${CC:-icc}" CXX="${CXX:-icpc}"
        else
            export FC="${FC:-ifx}" CC="${CC:-icx}" CXX="${CXX:-icpx}"
        fi
        for compiler in "$FC" "$CC" "$CXX"; do
            if ! command -v "$compiler" >/dev/null; then
                echo "Compiler not found: $compiler" >&2
                exit 1
            fi
        done
        ;;
    gnu)
        export FC=gfortran CC=gcc CXX=g++
        ;;
esac

project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
prefix="$project_root/.deps"
work="$project_root/.deps-build"
version=3.9.1
archive="$work/arpack-ng-$version.tar.gz"
mkdir -p "$work"

if [[ ! -f "$archive" ]]; then
    curl --fail --location --retry 3 \
        "https://codeload.github.com/opencollab/arpack-ng/tar.gz/refs/tags/$version" \
        --output "$archive.part"
    mv "$archive.part" "$archive"
fi
tar -xzf "$archive" -C "$work"

cmake -S "$work/arpack-ng-$version" -B "$work/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DMPI=ON
cmake --build "$work/build"
cmake --install "$work/build"

mkdir -p "$prefix/share/licenses/arpack-ng"
cp "$work/arpack-ng-$version/COPYING" "$prefix/share/licenses/arpack-ng/COPYING"
echo "ARPACK/PARPACK installed in $prefix"

# TecIO uses Boost headers only; download a fixed version into the local prefix.
boost_version=1.69.0
boost_dir=boost_1_69_0
boost_archive="$work/$boost_dir.tar.bz2"
if [[ ! -f "$boost_archive" ]]; then
    curl --fail --location --retry 3 \
        "https://archives.boost.io/release/$boost_version/source/$boost_dir.tar.bz2" \
        --output "$boost_archive.part"
    mv "$boost_archive.part" "$boost_archive"
fi
tar -xjf "$boost_archive" -C "$work" "$boost_dir/boost" "$boost_dir/LICENSE_1_0.txt"
mkdir -p "$prefix/include" "$prefix/share/licenses/boost"
cp -pR "$work/$boost_dir/boost" "$prefix/include/"
cp "$work/$boost_dir/LICENSE_1_0.txt" "$prefix/share/licenses/boost/"
boost_include="$prefix/include"
echo "Boost $boost_version headers installed in $boost_include"

tecio_revision=732bcc7f298a2c38b70b65865a2b2db488528649
tecio_archive="$work/tecio-$tecio_revision.tar.gz"
tecio_source="$work/tecio-$tecio_revision/teciosrc"
tecio_build="$work/tecio-build"
if [[ ! -f "$tecio_archive" ]]; then
    curl --fail --location --retry 3 \
        "https://codeload.github.com/Izumiko/tecio/tar.gz/$tecio_revision" \
        --output "$tecio_archive.part"
    mv "$tecio_archive.part" "$tecio_archive"
fi
if [[ ! -f "$tecio_source/CMakeLists.txt" ]]; then
    tar -xzf "$tecio_archive" -C "$work"
fi
cmake -S "$tecio_source" -B "$tecio_build" \
    -DCMAKE_BUILD_TYPE=Release -DBoost_INCLUDE_DIR="$boost_include" \
    -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX"
cmake --build "$tecio_build" --target tecio --parallel "${JOBS:-4}"

# Upstream has no install target. Keep its public interfaces and license together.
mkdir -p "$prefix/lib" "$prefix/include/tecio" "$prefix/share/licenses/tecio"
cp "$tecio_build/libtecio.a" "$prefix/lib/"
for interface in StandardIntegralTypes.h TECIO.h tecio_Exports.h tecio.inc tecio.for tecio.f90; do
    cp "$tecio_source/$interface" "$prefix/include/tecio/"
done
cp "$tecio_source/tecio_license_agreement.txt" "$prefix/share/licenses/tecio/"
echo "TecIO ($tecio_revision, serial teciosrc) installed in $prefix"
