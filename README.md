# histo-header
Histogram in C++11, 1D, simple header-only, inspired by R, calculate and optimize breaks automatically. Accepts different precissions.

# Compile the library.
git clone https://github.com/phcerdan/histo-header.git
mkdir histo-header/build
cd histo-header/build
ccmake ../src

The library will be generated in build/lib

You need a compiler with support for -std=c++11, and probably a modern cmake.

# Using it.
You should read the test file as a sample for usage.

# Test
If you are thinking about contributing, you should get gmock.
I have an environment variable called GMOCK_HOME pointing to the build directory, and the cmake reflects that.
