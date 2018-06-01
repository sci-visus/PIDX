The master branch is currently [![Build Status](https://travis-ci.org/sci-visus/PIDX.svg?branch=master)](https://travis-ci.org/sci-visus/PIDX) on Linux and OSX and <img src="https://ci.appveyor.com/api/projects/status/github/sci-visus/PIDX?svg=true" alt="Project Badge"/> on Windows.

--------------------------------------
PIDX
--------------------------------------

PIDX is an efficient parallel I/O library that reads and writes multiresolution IDX data files.

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

From shell:

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

  Link with ``-L/install/path/PIDX/lib -lpidx -lzfp`` (add to LDFLAGS).

--------------------------------------
Documentation
--------------------------------------

A few simple function calls are all that's necessary to instrument your application to utilize PIDX I/O:

Please see the examples in the "examples" folder and our wiki: <https://github.com/sci-visus/PIDX/wiki>

Documentation from the PIDX Tutorial can be found here: <https://sites.google.com/site/bigdatahpc/>

--------------------------------------
Support
--------------------------------------

**Users forum:**

Please join our forum at <http://forum.visus.org/>

**Modifications:**

We welcome machine-specific modifications and bug fixes: please report issues and send pull requests to <http://github.com/sci-visus/PIDX>.

**Extending:**

If you would like to extend PIDX, please contact <pascucci@sci.utah.edu>.

**Commercial Use:**

If you would like to use PIDX for commercial applications, please contact <pascucci@sci.utah.edu>.


--------------------------------------
Publications
--------------------------------------

Kumar, S., Humphrey, A., Usher, W., Petruzza, S., Peterson, B., Schmidt, J.A., Harris, D., Isaac, B., Thornock, J., Harman, T. and Pascucci, V., 2018, March. Scalable Data Management of the Uintah Simulation Framework for Next-Generation Engineering Problems with Radiation. In Asian Conference on Supercomputing Frontiers (pp. 219-240). Springer, Cham.

Kumar, S., Hoang, D., Petruzza, S., Edwards, J. and Pascucci, V., 2017, December. Reducing network congestion and synchronization overhead during aggregation of hierarchical data. In High Performance Computing (HiPC), 2017 IEEE 24th International Conference on (pp. 223-232). IEEE.

Kumar, S., Edwards, J., Bremer, P.T., Knoll, A., Christensen, C., Vishwanath, V., Carns, P., Schmidt, J.A. and Pascucci, V., 2014, November. Efficient I/O and storage of adaptive-resolution data. In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (pp. 413-423). IEEE Press.

Kumar, S., Saha, A., Vishwanath, V., Carns, P., Schmidt, J.A., Scorzelli, G., Kolla, H., Grout, R., Latham, R., Ross, R. and Papkafa, M.E., 2013, November. Characterization and modeling of pidx parallel I/O for performance optimization. In Proceedings of the International Conference on High Performance Computing, Networking, Storage and Analysis (p. 67). ACM.

Kumar, S., Vishwanath, V., Carns, P., Levine, J.A., Latham, R., Scorzelli, G., Kolla, H., Grout, R., Ross, R., Papka, M.E. and Chen, J., 2012, November. Efficient data restructuring and aggregation for I/O acceleration in PIDX. In High Performance Computing, Networking, Storage and Analysis (SC), 2012 International Conference for (pp. 1-11). IEEE.

Kumar, S., Vishwanath, V., Carns, P., Summa, B., Scorzelli, G., Pascucci, V., Ross, R., Chen, J., Kolla, H. and Grout, R., 2011, September. PIDX: Efficient parallel I/O for multi-resolution multi-dimensional scientific datasets. In Cluster Computing (CLUSTER), 2011 IEEE International Conference on (pp. 103-111). IEEE.

Kumar, S., Pascucci, V., Vishwanath, V., Carns, P., Hereld, M., Latham, R., Peterka, T., Papka, M.E. and Ross, R., 2010, November. Towards parallel access of multi-dimensional, multi-resolution scientific data. In Petascale Data Storage Workshop (PDSW), 2010 5th (pp. 1-5). IEEE.
