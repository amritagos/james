project('james', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

# Library dependencies
graphlib_dep = dependency('graph_lib', fallback : ['graph_lib', 'graphlib_dep'])

inc = include_directories('james/include')
james_dep = declare_dependency(include_directories : inc, dependencies: graphlib_dep)

tests = [
  ['Test_Dummy', 'test/test_dummy.cpp'],
]

test_inc = []
test_inc += inc

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

fmt_dep = dependency('fmt')

foreach t : tests
  exe = executable(t.get(0), t.get(1),
    dependencies : [Catch2, fmt_dep, graphlib_dep],
    include_directories : test_inc
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach