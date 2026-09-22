# Benchmarks: bench/ayzek_bench.cpp, docs/impl/10-benchmarks.md.
#
#   meson test -C build-release --benchmark          all groups, results in
#                                                     build-release/bench/*.json
#   build-release/bench/ayzek_bench --help
bench_version = vcs_tag(
  command: ['git', '-C', meson.project_source_root(), 'describe', '--always', '--dirty'],
  input: 'version.hpp.in',
  output: 'bench_version.hpp',
)
ayzek_bench = executable('ayzek_bench', 'ayzek_bench.cpp', bench_version,
  dependencies: [ayzek_dep],
  install: false)

models_dir = meson.project_source_root() / 'models'
foreach g : ['stages', 'layers', 'capacity']
  benchmark(g, ayzek_bench,
    args: ['--models', models_dir, '--group', g,
           '--json', meson.current_build_dir() / g + '.json'],
    timeout: 900)
endforeach

# Keeps the benchmark building and running in the ordinary test suite: a few
# calls of everything, on any build type.
test('bench_quick', ayzek_bench,
  args: ['--models', models_dir, '--group', 'stages,layers,capacity', '--quick',
         '--threads', '2', '--allow-debug'],
  timeout: 300, is_parallel: false)
