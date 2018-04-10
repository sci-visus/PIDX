# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]
### Added
- compression and decompression (zfp)
- particle data interface
- support for Windows build
- continuous integration support for Windows (Appveyor)
- serial file open (for VisIt reader)

### Changed
- improved testing scripts and autotesting for Travis CI

### Fixed
- raw read performance
- zfp library build and deploy in CMake 

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
