# Unit tests.
foreach t : ['test_ring', 'test_steim', 'test_reorder', 'test_writer', 'test_locate', 'test_trigger', 'test_assess', 'test_seedlink']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe)
endforeach

# Decoding of every supported encoding against ObsPy's reading of the same
# files (tests/reference/mseed, tools/make_mseed_fixtures.py).
test('test_mseed', executable('test_mseed', 'test_mseed.cpp', dependencies: [ayzek_dep]),
  args: [meson.project_source_root()])

# Comparisons with reference outputs from tools/export_models.py,
# tools/export_stalta.py and tools/export_pd.py.
# Skipped (exit 77) if those files are missing.
foreach t : ['test_dsp', 'test_models', 'test_stalta', 'test_transformer', 'test_pd']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe, args: [meson.project_source_root()], timeout: 300)
endforeach

# FilterPicker against Lomax's C function (reference/filterpicker/README.md),
# which is built here only, as he wrote it.
if add_languages('c', required: false, native: false)
  fp_reference = static_library('fp_reference', 'reference/filterpicker/fp_reference.c',
    override_options: ['werror=false', 'warning_level=0'], c_args: ['-w'])
  exe = executable('test_filterpicker', 'test_filterpicker.cpp',
    dependencies: [ayzek_dep], link_with: fp_reference)
  test('test_filterpicker', exe, args: [meson.project_source_root()], timeout: 300)
endif

# End-to-end tests on the demo data; see demo_check.sh.
test('demo_m49', find_program('demo_check.sh'),
  args: [ayzek_exe, meson.project_source_root(), 'demo', '0'], timeout: 600, is_parallel: false)
test('demo_ko_m49', find_program('demo_check.sh'),
  args: [ayzek_exe, meson.project_source_root(), 'demo_ko', '1'], timeout: 600, is_parallel: false)

# Live SeedLink input against tools/seedlink_mock.py on demo_ko; see
# seedlink_check.sh. Native builds only: under qemu ayzek would not keep up
# with the mock server.
if not meson.is_cross_build()
  foreach mode : ['live', 'reconnect']
    test('seedlink_' + mode, find_program('seedlink_check.sh'),
      args: [ayzek_exe, meson.project_source_root(), mode], timeout: 1800, is_parallel: false)
  endforeach
endif
