# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]
### Added
- compression and decompression (zfp)
- particles data read and write (file per process)
- particles restructuring (experimental)
- support for Windows build
- continuous integration support for Windows (Appveyor)
- serial file open (for VisIt reader)
- partitioning examples
- generate midx when local partitioned IDX is used

### Changed
- metadata format 6.1 and support for previous (6)
- use enums for PIDX types
- use literal io mode, endianess
- refactoring examples 
- improved testing scripts and autotesting for Travis CI

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
