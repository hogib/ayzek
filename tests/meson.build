# Unit tests.
foreach t : ['test_ring', 'test_steim', 'test_reorder']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe)
endforeach

# Comparisons with reference outputs from tools/export_models.py and
# tools/export_magnitude.py. Skipped (exit 77) if those files are missing.
foreach t : ['test_dsp', 'test_models', 'test_magnitude']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe, args: [meson.project_source_root()], timeout: 300)
endforeach

# End-to-end test on the demo data; see demo_check.sh.
test('demo_m49', find_program('demo_check.sh'),
  args: [ayzek_exe, meson.project_source_root()], timeout: 600, is_parallel: false)
