ayzek_lib = static_library('ayzek',
  sources: files('weights.cpp', 'nn.cpp', 'models.cpp', 'dsp.cpp'),
  include_directories: inc,
  dependencies: [threads],
)

ayzek_dep = declare_dependency(
  link_with: ayzek_lib,
  include_directories: inc,
  dependencies: [threads],
)
