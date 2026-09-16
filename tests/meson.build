# Each stage adds its test here. A ring buffer that is subtly wrong drops
# samples silently, so the tests are the point rather than a formality.
foreach t : ['test_ring', 'test_steim', 'test_reorder']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe)
endforeach

# Compared against PyTorch and scipy outputs from tools/export_models.py. They
# skip (exit 77) until those files exist.
foreach t : ['test_dsp', 'test_models']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe, args: [meson.project_source_root()], timeout: 300)
endforeach

# End to end on real data: replay three stations, require the M4.9 to be
# declared and located within 10 km. Skips without the demo data.
test('demo_m49', find_program('demo_check.sh'),
  args: [ayzek_exe, meson.project_source_root()], timeout: 600, is_parallel: false)
