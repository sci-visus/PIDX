Getting Started
===============================================

--------------------------------------
Prerequisites 
--------------------------------------

Download and install cmake from:: 

  http://www.cmake.org/cmake/resources/software.html

Get source code from::

	https://github.com/sci-visus/PIDX


--------------------------------------
Compilation
--------------------------------------

From shell::

	cd <path/to/pidx>
	mkdir -p build 
	cd build
	cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/install/path/PIDX
	make && make install


--------------------------------------
Optional Build Settings
--------------------------------------

Explicitly control whether to build static or shared library::

  -DBUILD_SHARED_LIBS=[ON|OFF]

Enable Examples::

  -DPIDX_BUILD_EXAMPLES=[ON|OFF]

Enable Profiling::

  -DPIDX_BUILD_PROFILE=[ON|OFF]


--------------------------------------
Your Project
--------------------------------------

CMake:
  ``FindPIDX.cmake`` is provided in root installation directory after running ``make install``. It configures PIDX_INCLUDE_DIR and PIDX_LIBRARIES for use by your application.

Makefile:
  Compile with ``-I/install/path/PIDX`` (add to CFLAGS).

  Link with ``-L/install/path/PIDX/lib -lpidx`` (add to LDFLAGS).

--------------------------------------
Documentation
--------------------------------------

A few simple function calls are all that's necessary to instrument your application to utilize PIDX I/O:

Please see the examples in the "examples" folder for a demonstration.

Documentation from the PIDX Tutorial can be found here: <https://sites.google.com/site/bigdatahpc/>
