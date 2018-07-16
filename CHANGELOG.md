# Changelog
All notable changes to this project will be documented in this file.

## [0.9.2]
### Added
- support for Windows build
- continuous integration support for Windows (Appveyor)
- serial file open (for VisIt reader)
- partitioning, compression and particles examples
- autotests for partitioned and compression io mode
- generate midx when local partitioned IDX is used
- particles data read and write (file per process)
- particles two phase I/O (experimental)
- compression and decompression (zfp, experimental)
- comments on most of the code

### Changed
- metadata format 6.1 and support for previous (6)
- use enums for PIDX types
- use literal io mode, endianess
- refactoring examples 
- refactoring aggregation
- license BSD-3
- improved autotesting scripts (configuration) for Travis CI

### Fixed
- raw read performance
- zfp library build and deploy in CMake 
- zfp types redefinition

## [0.9.1]
### Added
- Procedural testing capability using scripts
that generate and validate the data dumped on disk
for different configurations

### Fixed
- Fixed minor bugs in read and write

### Changed
- Folder structure: examples, develop, tools
- Profiling executables are now part of tools

## [0.9.0] - 2015-05-12
### Added
- First release: write and read in IDX format
