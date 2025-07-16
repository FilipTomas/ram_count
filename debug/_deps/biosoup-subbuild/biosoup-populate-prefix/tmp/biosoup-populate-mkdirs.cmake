# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/ftomas/github/ram/debug/_deps/biosoup-src")
  file(MAKE_DIRECTORY "/home/ftomas/github/ram/debug/_deps/biosoup-src")
endif()
file(MAKE_DIRECTORY
  "/home/ftomas/github/ram/debug/_deps/biosoup-build"
  "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix"
  "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/tmp"
  "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/src/biosoup-populate-stamp"
  "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/src"
  "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/src/biosoup-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/src/biosoup-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/ftomas/github/ram/debug/_deps/biosoup-subbuild/biosoup-populate-prefix/src/biosoup-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
