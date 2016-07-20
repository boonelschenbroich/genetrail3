GeneTrail2
==========
This is the GeneTrail2 C++ library, a collection of algorithms and statistics
for processing gene expression data. It provides the compute kernels for the
[GeneTrail2 web service](https://genetrail2.bioinf.uni-sb.de).

Citing
------
If you use this code in your publication, please consider citing the following
publication:

> Stöckel D., Kehl T., Trampert P., Schneider L., Backes C., Ludwig N., Gerasch A.,
> Kaufmann M., Gessler M., Graf N., Meese E., Keller A., Lenhof H.-P.,
> *Multi-omics Enrichment Analysis using the GeneTrail2 Web Service*, Bioinformatics 2016,
> [doi: 10.1093/bioinformatics/btv770](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv770)

Installation
------------
Currently we only support building GeneTrail2 under Linux. Windows and MacOS X
are not supported. For compiling GeneTrail2 you need the following software:

- [GCC](https://gcc.gnu.org/) >= 4.9 or [clang](http://clang.llvm.org/) >= 3.5
- [CMake](https://cmake.org/) >= 2.8.12.2
- [Boost](http://www.boost.org/) >= 1.55
- [Eigen](http://eigen.tuxfamily.org/) >= 3.2.3
- [RapidJSON](https://github.com/miloyip/rapidjson) >= 1.0.2
- [GMP](https://gmplib.org/) >= 5.0.0 (optional)
- [Google Test Framework](https://github.com/google/googletest/) >= 1.7.0 (optional)

Create a directory named `build` in the source directory. Inside the build
directory type:

	cmake -DEIGEN3_INCLUDE_DIR=/path/to/eigen3/ -DGTEST_SRC_DIR=/path/to/gtest/ -DCMAKE_BUILD_TYPE=Release ..

If all goes well you can now compile a version of GeneTrail2
by typing:

	make

For installing GeneTrail2 specify the installation path using `-DCMAKE_INSTALL_PATH=`
and type:

	make install

If you built with support for the Google Test Framework, you can run the unit
tests by typing

	make test

A suite of more expensive integration tests can be started by typing

	make integration_test

License
-------
The code of the GeneTrail2 C++ library is licensed under the *GNU Lesser General
Public License 3 (LGPL v3)*. This means you are allowed to use the library in
your code with little restrictions. Especially you do not need to make your
source code available as long as you are just using the library. However, you
are required to open source all your changes to the library itself under the
conditions of the LGPLv3.

The GNU project hosts an extensive [FAQ](https://www.gnu.org/licenses/gpl-faq.html)
covering possible licensing questions.
