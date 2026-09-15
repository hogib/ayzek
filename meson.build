project('ayzek', 'cpp',
  version: '0.1.0',
  default_options: [
    'cpp_std=c++23',
    'warning_level=3',
    'werror=true',
    'b_ndebug=if-release',
  ])

# C++23 modules were tested on this toolchain and do not work through Meson:
# GCC 16.2 and Clang 22.1 both build modules correctly by hand, but Meson 1.12
# only wires dyndep scanning for MSVC. Forcing it on via -fmodules in cpp_args
# reaches ninja with a malformed dyndep ("inputs may not also have inputs").
# Headers until that is fixed upstream; none of the code below depends on it.

cpp = meson.get_compiler('cpp')
threads = dependency('threads')

inc = include_directories('src')
subdir('src')
subdir('tests')
