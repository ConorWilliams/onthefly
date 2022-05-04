[![Actions Status](https://github.com/ConorWilliams/onthefly/workflows/MacOS/badge.svg)](https://github.com/ConorWilliams/onthefly/actions)
[![Actions Status](https://github.com/ConorWilliams/onthefly/workflows/Windows/badge.svg)](https://github.com/ConorWilliams/onthefly/actions)
[![Actions Status](https://github.com/ConorWilliams/onthefly/workflows/Ubuntu/badge.svg)](https://github.com/ConorWilliams/onthefly/actions)
[![Actions Status](https://github.com/ConorWilliams/onthefly/workflows/Style/badge.svg)](https://github.com/ConorWilliams/onthefly/actions)
[![Actions Status](https://github.com/ConorWilliams/onthefly/workflows/Intel/badge.svg)](https://github.com/ConorWilliams/onthefly/actions)

# OnTheFly: Off-lattice kinetic Monte-Carlo

## Compilation

### Build and run the standalone target

Use the following command to build and run the executable target.

```bash
cmake -S standalone -B build/standalone
cmake --build build/standalone
./build/standalone/OnTheFly --help
```

### Build and run test suite

Use the following commands from the project's root directory to run the test suite.

```bash
cmake -S test -B build/test
cmake --build build/test
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test

# or simply call the executable: 
./build/test/OnTheFlyTests
```

### Run clang-format

Use the following commands from the project's root directory to check and fix C++ and CMake source style.
This requires _clang-format_, _cmake-format_ and _pyyaml_ to be installed on the current system.

```bash
cmake -S test -B build/test

# view changes
cmake --build build/test --target format

# apply changes
cmake --build build/test --target fix-format
```


### Build the documentation

The documentation is automatically built and [published](https://thelartians.github.io/ModernCppStarter) whenever a [GitHub Release](https://help.github.com/en/github/administering-a-repository/managing-releases-in-a-repository) is created.
To manually build documentation, call the following command.

```bash
cmake -S docs -B build/docs
cmake --build build/docs --target GenerateDocs
# view the docs
open build/docs/doxygen/html/index.html
```

To build the documentation locally, you will need Doxygen, jinja2 and Pygments on installed your system.

### Build everything at once

The project also includes an `all` directory that allows building all targets at the same time.
This is useful during development, as it exposes all subprojects to your IDE and avoids redundant builds of the library.

```bash
cmake -S all -B build
cmake --build build

# run tests
./build/test/OnTheFlyTests
# format code
cmake --build build --target fix-format
# run standalone
./build/standalone/OnTheFly --help
# build docs
cmake --build build --target GenerateDocs
```

### Ccache

Ccache can be enabled by configuring with `-DUSE_CCACHE=<ON | OFF>`.
