# Unit tests. When b_sanitize includes undefined, AYZEK_UBSAN tells
# test_mseed_fuzz to keep clear of a known overflow (see there).
unit_args = '@0@'.format(get_option('b_sanitize')).contains('undefined') ? ['-DAYZEK_UBSAN'] : []
foreach t : ['test_ring', 'test_steim', 'test_mseed_fuzz', 'test_reorder', 'test_writer', 'test_locate', 'test_trigger', 'test_assess', 'test_network']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep], cpp_args: unit_args)
  test(t, exe)
endforeach

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
