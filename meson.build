project('ayzek', 'cpp',
  version: '0.1.0',
  default_options: [
    'cpp_std=c++23',
    'warning_level=3',
    'werror=true',
    'b_ndebug=if-release',
  ])

# C++23 modules are not used: Meson 1.12 supports module dependency scanning
# only for MSVC, and enabling -fmodules with GCC 16.2 or Clang 22.1 produces an
# invalid ninja dyndep file.

cpp = meson.get_compiler('cpp')
threads = dependency('threads')

inc = include_directories('src')
subdir('src')
subdir('app')
subdir('tests')
subdir('tools')
