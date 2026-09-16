# Each stage adds its test here. A ring buffer that is subtly wrong drops
# samples silently, so the tests are the point rather than a formality.
foreach t : ['test_ring', 'test_steim']
  exe = executable(t, t + '.cpp', dependencies: [ayzek_dep])
  test(t, exe)
endforeach
