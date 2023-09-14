# TODO:
# - move to fortran_std=f2018; this requires changes in libkordersim, and also replacing isnan() by ieee_is_nan() in various MEX (or use -fall-intrinsics)
# - use buildtype=debugoptimized by default? (and in the preprocessor too)
# - with -Dprefer_static=true, under Octave/Windows, we are now linking the compiler libs (libgcc, libstdc++, libgfortran, libquadmath, libssp, libgomp) statically (contrary to what we were doing with autotools). In theory this is better, because the compiler used for creating the Octave binary may be different from the one used for creating our MEX. Check that this is ok
# - configuration option to disable documentation
# - configuration option to disable preprocessor build
# - add -Wold-style-cast C++ flag except when building flex-generated files
# - determine minimal meson version required, and declare it in the project() function (here and the preprocessor); if possible at low cost, diminish the requirement. NB: Ubuntu “jammy” 22.04 has 0.61, Debian “bullseye” 11 has 0.56 (but bullseye-backports has 1.0), Ubuntu ”focal” 20.04 has 0.53

project('dynare',
        'cpp', 'fortran', 'c',
        version : '6-unstable',
        default_options : [ 'cpp_std=gnu++20', 'fortran_std=none', 'c_std=gnu17', 'warning_level=2' ])

add_global_arguments('-Wimplicit-interface', '-Wno-compare-reals', language : 'fortran')
add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

cpp_compiler = meson.get_compiler('cpp')
fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')

### Preprocessor

subdir('preprocessor/src')
subdir('preprocessor/doc')

### Generated M-file

sed_exe = find_program('sed')
custom_target(output : 'dynare_version.m', input : 'matlab/dynare_version.m.in',
              command : [ sed_exe, 's/@PACKAGE_VERSION@/' + meson.project_version() + '/', '@INPUT@' ],
              capture : true,
              build_by_default : true) # FIXME: This option can be removed when “install” is set to true

### MEX files

mex_incdir = include_directories('mex/sources')

## Various dependencies

if host_machine.system() != 'windows'
  dl_dep = dependency('dl')
else
  # Under Windows, we don’t use dlopen but rely on LoadLibrary
  dl_dep = []
endif

openmp_dep = dependency('openmp')
gsl_dep = dependency('gsl')

if get_option('build_for') == 'octave'
  matio_dep = dependency('matio')
else
  # We don’t use MatIO under MATLAB
  matio_dep = []
endif

pthread_t_sizeof = c_compiler.sizeof('pthread_t', prefix : '#include <pthread.h>')
# TODO: when requiring meson ⩾ 1.2, incorporate pthread.F08 into the dependency object with “extra_files” option
# (and thus remove the pthread_fortran_iface variable)
pthread_fortran_dep = declare_dependency(compile_args : '-DSIZEOF_PTHREAD_T=' + pthread_t_sizeof.to_string())

## Determine MEX compilation options

if get_option('build_for') == 'matlab'
  matlab_path = get_option('matlab_path')
  if matlab_path == ''
    error('The “matlab_path” option must be provided when doing a MATLAB build')
  endif
  matlab_version = run_command('scripts/get-matlab-version', matlab_path, check : true).stdout().strip()
  matlab_minimal_version = [ '8.3', 'R2014a' ]

  if matlab_version.version_compare('<' + matlab_minimal_version[0])
    error('MATLAB is too old (version ' + matlab_version + '), please upgrade to version ' + matlab_minimal_version[0] + ' (' + matlab_minimal_version[1] + ') at least.')
  endif

  matlab_version_hex = run_command('scripts/get-matlab-version', '--hex', matlab_path, check : true).stdout().strip()
  matlab_exe = find_program(matlab_path / 'bin' / 'matlab', required : not meson.is_cross_build())

  if host_machine.system() == 'linux' and host_machine.cpu_family() == 'x86_64'
    mexext = 'mexa64'
    matlab_arch = 'glnxa64'
    export_file = matlab_path / 'extern/lib/glnxa64/mexFunction.map'
    export_link_arg = '-Wl,--version-script,' + export_file
  elif host_machine.system() == 'windows' and host_machine.cpu_family() == 'x86_64'
    mexext = 'mexw64'
    matlab_arch = 'win64'
    export_file = meson.current_source_dir() / 'mex/build/matlab/mex.def'
    export_link_arg = export_file
    if get_option('build_for') == 'matlab'
      arch_fortran_args = [ '-fno-underscoring' ]
    endif
  elif host_machine.system() == 'darwin'
    if host_machine.cpu_family() == 'x86_64'
      mexext = 'mexmaci64'
      matlab_arch = 'maci64'
    elif host_machine.cpu_family() == 'aarch64'
      mexext = 'mexmaca64'
      matlab_arch = 'maca64'
    else
      error('Unsupported platform')
    endif
    export_file = meson.current_source_dir() / 'mex/build/matlab/mexFunction-MacOSX.map'
    export_link_arg = '-Wl,-exported_symbols_list,' + export_file
  else
    error('Unsupported platform')
  endif

  matlab_defs = [ '-DMATLAB_MEX_FILE', '-DMATLAB_VERSION=' + matlab_version_hex, '-DMEXEXT=".' + mexext + '"' ]
  matlab_incdir = include_directories(matlab_path / 'extern/include')

  mex_kwargs = { 'name_prefix' : '',
                 'name_suffix' : mexext,
                 'include_directories' : [ mex_incdir, matlab_incdir ],
                 'cpp_args' : matlab_defs,
                 'fortran_args' : matlab_defs + [ '-fexceptions' ] + get_variable('arch_fortran_args', []),
                 'c_args' : matlab_defs + [ '-fexceptions' ],
                 'link_args' : [ export_link_arg, '-L' + (matlab_path / 'bin' / matlab_arch), '-lmx', '-lmex', '-lmat' ],
                 'link_depends' : export_file }

  # For unit tests
  exe_rpath = matlab_path / 'bin' / matlab_arch
  exe_link_args = [ '-L' + exe_rpath, '-lmx', '-lmex', '-lmat' ]

  # No need to use find_library() for the following libraries, since they are always shipped with MATLAB
  blas_dep = declare_dependency(link_args : '-lmwblas')
  lapack_dep = declare_dependency(link_args : '-lmwlapack', dependencies : blas_dep)
  umfpack_dep = declare_dependency(link_args : '-lmwumfpack', dependencies : blas_dep)
  ut_dep = declare_dependency(link_args : '-lut')

  slicot_dep = declare_dependency(dependencies : [ fortran_compiler.find_library('slicot64_pic'), blas_dep, lapack_dep ])
else # Octave build
  octave_exe = find_program('octave', required : not meson.is_cross_build())
  mkoctfile_exe = find_program('mkoctfile')
  octave_minimal_version = '6.2.0'
  octave_version = run_command(mkoctfile_exe, '-v', check : true).stdout().replace('mkoctfile, version ', '')

  if octave_version.version_compare('<' + octave_minimal_version)
    error('Octave is too old (version ' + octave_version + '), please upgrade to version ' + octave_minimal_version + ' at least.')
  endif

  octave_incflags = run_command(mkoctfile_exe, '-p', 'INCFLAGS', check : true).stdout().split()
  octlibdir = run_command(mkoctfile_exe, '-p', 'OCTLIBDIR', check : true).stdout().strip()

  # Determine whether to link MEX files against the Octave libraries. mkoctfile
  # no longer does this by default but in practice it is needed for Windows and
  # macOS.
  octave_libs = run_command(mkoctfile_exe, '-p', 'OCTAVE_LIBS', check : true).stdout().split()
  if host_machine.system() == 'windows' or host_machine.system() == 'darwin'
    # Under Windows, --enable-link-all-dependencies is hardcoded in src/mkoctfile.cc.in.
    # Under macOS, the Homebrew formula passes --enable-link-all-dependencies
    # to the configure script.
    octave_link_args = [ '-L' + octlibdir ] + octave_libs
  else
    octave_link_args = []
  endif

  # For unit tests
  exe_rpath = octlibdir
  exe_link_args = [ '-L' + octlibdir ] + octave_libs

  octave_defs = [ '-DOCTAVE_MEX_FILE', '-DMEXEXT="mex"' ]

  mex_kwargs = { 'name_prefix' : '',
                 'name_suffix' : 'mex',
                 'include_directories' : [ mex_incdir ],
                 'cpp_args' : octave_incflags + octave_defs,
                 'fortran_args' : octave_incflags + octave_defs,
                 'c_args' : octave_incflags + octave_defs,
                 'link_args' : octave_link_args}

  # The -L argument is useful when cross-compiling.
  blas_dep = declare_dependency(link_args : [ '-L' + (octlibdir / '../..') ] + run_command(mkoctfile_exe, '-p', 'BLAS_LIBS', check : true).stdout().split())
  lapack_dep = declare_dependency(link_args : run_command(mkoctfile_exe, '-p', 'LAPACK_LIBS', check : true).stdout().split(),
                                  dependencies : blas_dep)

  # Create a dependency object for UMFPACK.
  # The dependency returned by find_library('umfpack') is not enough, because we also want the define
  # that indicates the location of umfpack.h, so we construct a new dependency object.
  if cpp_compiler.has_header('suitesparse/umfpack.h', args : octave_incflags)
    umfpack_def = '-DHAVE_SUITESPARSE_UMFPACK_H'
  elif cpp_compiler.has_header('umfpack.h', args : octave_incflags)
    umfpack_def = '-DHAVE_UMFPACK_H'
  else
    error('Can’t find umfpack.h')
  endif
  # Do not enforce static linking even if prefer_static is true, since that library is shipped
  # with Octave.
  # The “dirs” argument is useful when cross-compiling.
  umfpack_dep_tmp = cpp_compiler.find_library('umfpack', dirs : octlibdir / '../..', static : false)
  umfpack_dep = declare_dependency(compile_args : umfpack_def, dependencies : [ umfpack_dep_tmp, blas_dep ])

  # This library does not exist under Octave
  ut_dep = []

  # First look for libslicot, then if needed fallback on libslicot_pic
  slicot_dep_tmp = fortran_compiler.find_library('slicot', required : false)
  if not slicot_dep_tmp.found()
    slicot_dep_tmp = fortran_compiler.find_library('slicot_pic')
  endif
  slicot_dep = declare_dependency(dependencies : [ slicot_dep_tmp, blas_dep, lapack_dep ])
endif

# When static linking is preferred, try to statically link the compiler libraries
if get_option('prefer_static')
  static_flags_pre = [ '-static-libgcc', '-static-libstdc++' ]
  static_flags_post = []
  if host_machine.system() == 'windows'
    # Under Debian 12, libgfortran.a and libquadmath.a are not compiled with -fPIC, so can’t be linked in a MEX.
    # Under macOS, -static-libgcc implies -static-libgfortran (see gfortran -dumpspecs), and for libquadmath
    # we use a hack with a local copy of libquadmath.a.
    # -static-libquadmath was introduced in GCC 13. Until we require the latter, use a hack.
    static_flags_pre += [ '-static-libgfortran', '-Wl,-Bstatic,--whole-archive', '-lquadmath', '-Wl,-Bdynamic,--no-whole-archive' ]

    # Hack to avoid dynamically linking against libwinpthread DLL (which is
    # pulled in by libstdc++, even without using threads, since we are using
    # the POSIX threads version of MinGW).
    static_flags_pre += [ '-Wl,-Bstatic,--whole-archive', '-lwinpthread', '-Wl,-Bdynamic,--no-whole-archive' ]

    # Hack for libssp, which is pulled in by -fstack-protector (curiously only
    # on some MEX files), see windows/build.sh. Note that the link against
    # libssp should not happen with compilers from MSYS2, see:
    # https://www.msys2.org/news/#2022-10-10-libssp-is-no-longer-required
    # But it happens with Debian’s cross compilers (as of Debian “bookworm”
    # 12). Also note that the -lssp must come by the end of the link command
    # (otherwise it will have to be enclosed within --whole-archive).
    static_flags_post = [ '-Wl,-Bstatic', '-lssp', '-Wl,-Bdynamic' ]
  endif

  mex_kwargs += { 'link_args' : static_flags_pre + mex_kwargs.get('link_args', []) + static_flags_post }

  # NB: constructing a dependency object with link_args : ['-Wl,-Bstatic', '-lgomp', '-Wl,-Bdynamic'] does not work,
  # because it reorders the three arguments and puts -lgomp at the end
  openmp_dep_tmp = cpp_compiler.find_library('gomp', static : true)
  openmp_dep = declare_dependency(dependencies : [ openmp_dep, openmp_dep_tmp ])
endif

# For use when creating intermediate static libraries to be incorporated in MEX files
static_library_kwargs = mex_kwargs + { 'name_prefix' : [], 'name_suffix' : [], 'link_args' : [], 'pic' : true }

mex_blas_fortran_iface = [ 'mex/sources/matlab_mex.F08', 'mex/sources/blas_lapack.F08' ]
pthread_fortran_iface = [ 'mex/sources/pthread.F08' ]


## Various core MEX

shared_module('mjdgges', [ 'mex/sources/mjdgges/mjdgges.F08' ] + mex_blas_fortran_iface, kwargs : mex_kwargs, dependencies : lapack_dep)

shared_module('num_procs', 'mex/sources/num_procs/num_procs.cc', kwargs : mex_kwargs)

perfect_foresight_problem_src = [ 'mex/sources/perfect_foresight_problem/perfect_foresight_problem.cc',
                                  'mex/sources/perfect_foresight_problem/DynamicModelCaller.cc' ]
shared_module('perfect_foresight_problem', perfect_foresight_problem_src, kwargs : mex_kwargs, dependencies : [ openmp_dep, dl_dep ])

block_trust_region_src = [ 'mex/sources/block_trust_region/dulmage_mendelsohn.f08',
                           'mex/sources/block_trust_region/matlab_fcn_closure.F08',
                           'mex/sources/block_trust_region/trust_region.f08',
                           'mex/sources/block_trust_region/mexFunction.f08' ] + mex_blas_fortran_iface
shared_module('block_trust_region', block_trust_region_src, kwargs : mex_kwargs, dependencies : lapack_dep)

bytecode_src = [ 'mex/sources/bytecode/bytecode.cc',
                 'mex/sources/bytecode/Interpreter.cc',
                 'mex/sources/bytecode/Mem_Mngr.cc',
                 'mex/sources/bytecode/SparseMatrix.cc',
                 'mex/sources/bytecode/Evaluate.cc',
                 'mex/sources/bytecode/BasicSymbolTable.cc' ]
preprocessor_headers_dep = declare_dependency(include_directories : include_directories('preprocessor/src'))
shared_module('bytecode', bytecode_src, kwargs : mex_kwargs, dependencies : [ umfpack_dep, ut_dep, preprocessor_headers_dep ])

shared_module('sparse_hessian_times_B_kronecker_C', 'mex/sources/kronecker/sparse_hessian_times_B_kronecker_C.cc',
              kwargs : mex_kwargs, dependencies : openmp_dep)

# TODO: A_times_B_kronecker_C does not depend on LAPACK, but since the
# interfaces to both BLAS and LAPACK are in the same source file, the
# dependency must be added. Think about splitting into two files.
# NB: The problem does not appear with libkordersim because since it is a library,
# the LAPACK stuff is never pulled in.
shared_module('A_times_B_kronecker_C', [ 'mex/sources/kronecker/A_times_B_kronecker_C.f08' ] + mex_blas_fortran_iface,
              kwargs : mex_kwargs, dependencies : [ blas_dep, lapack_dep ])

shared_module('cycle_reduction', [ 'mex/sources/cycle_reduction/mexFunction.f08' ] + mex_blas_fortran_iface,
              kwargs : mex_kwargs, dependencies : [ blas_dep, lapack_dep ])

shared_module('logarithmic_reduction', [ 'mex/sources/logarithmic_reduction/mexFunction.f08' ] + mex_blas_fortran_iface,
              kwargs : mex_kwargs, dependencies : [ blas_dep, lapack_dep ])

shared_module('disclyap_fast', [ 'mex/sources/disclyap_fast/disclyap_fast.f08' ] + mex_blas_fortran_iface,
              kwargs : mex_kwargs, dependencies : [ blas_dep, lapack_dep ])

# TODO: Same remark as A_times_B_kronecker_C
shared_module('riccati_update', [ 'mex/sources/riccati_update/mexFunction.f08' ] + mex_blas_fortran_iface,
              kwargs : mex_kwargs, dependencies : [ blas_dep, lapack_dep ])

qmc_sequence_src = [ 'mex/sources/sobol/qmc_sequence.cc',
                     'mex/sources/sobol/sobol.f08' ]
# Hack for statically linking libgfortran
# Since qmc_sequence is a mix of C++ and Fortran, the linker invoked is the C++ one.
# Meson then rightly adds -lgfortran, but the -static-libgfortran flag does not work.
qmc_sequence_mex_kwargs = mex_kwargs
if get_option('prefer_static') and host_machine.system() == 'windows'
  qmc_sequence_mex_kwargs += { 'link_args' : qmc_sequence_mex_kwargs.get('link_args') + [ '-Wl,-Bstatic', '-lgfortran', '-Wl,-Bdynamic' ] }
endif
shared_module('qmc_sequence', qmc_sequence_src, kwargs : qmc_sequence_mex_kwargs, dependencies : [ blas_dep, openmp_dep ])

shared_module('kalman_steady_state', 'mex/sources/kalman_steady_state/kalman_steady_state.cc', kwargs : mex_kwargs, dependencies : slicot_dep)


## k-order simulation stuff

kordersim_src = [ 'mex/sources/libkordersim/pascal.f08',
                  'mex/sources/libkordersim/sort.f08',
                  'mex/sources/libkordersim/partitions.f08',
                  'mex/sources/libkordersim/tensors.f08',
                  'mex/sources/libkordersim/simulation.f08',
                  'mex/sources/libkordersim/struct.f08' ] + mex_blas_fortran_iface + pthread_fortran_iface
kordersim_lib = static_library('kordersim', kordersim_src, kwargs : static_library_kwargs, dependencies : [ blas_dep, pthread_fortran_dep ])

shared_module('folded_to_unfolded_dr', 'mex/sources/folded_to_unfolded_dr/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('k_order_mean', 'mex/sources/k_order_mean/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('k_order_simul', 'mex/sources/k_order_simul/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)

shared_module('local_state_space_iteration_2', 'mex/sources/local_state_space_iterations/local_state_space_iteration_2.cc', kwargs : mex_kwargs, dependencies : openmp_dep)
shared_module('local_state_space_iteration_3', 'mex/sources/local_state_space_iterations/local_state_space_iteration_3.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('local_state_space_iteration_k', 'mex/sources/local_state_space_iterations/local_state_space_iteration_k.f08', kwargs : mex_kwargs, link_with : kordersim_lib)


## k-order resolution stuff

korder_src = [ 'mex/sources/libkorder/kord/approximation.cc',
               'mex/sources/libkorder/kord/decision_rule.cc',
               'mex/sources/libkorder/kord/dynamic_model.cc',
               'mex/sources/libkorder/kord/faa_di_bruno.cc',
               'mex/sources/libkorder/kord/first_order.cc',
               'mex/sources/libkorder/kord/korder.cc',
               'mex/sources/libkorder/kord/korder_stoch.cc',
               'mex/sources/libkorder/kord/journal.cc',
               'mex/sources/libkorder/sylv/BlockDiagonal.cc',
               'mex/sources/libkorder/sylv/GeneralMatrix.cc',
               'mex/sources/libkorder/sylv/GeneralSylvester.cc',
               'mex/sources/libkorder/sylv/IterativeSylvester.cc',
               'mex/sources/libkorder/sylv/KronUtils.cc',
               'mex/sources/libkorder/sylv/KronVector.cc',
               'mex/sources/libkorder/sylv/QuasiTriangular.cc',
               'mex/sources/libkorder/sylv/QuasiTriangularZero.cc',
               'mex/sources/libkorder/sylv/SchurDecomp.cc',
               'mex/sources/libkorder/sylv/SchurDecompEig.cc',
               'mex/sources/libkorder/sylv/SimilarityDecomp.cc',
               'mex/sources/libkorder/sylv/SylvException.cc',
               'mex/sources/libkorder/sylv/SylvMatrix.cc',
               'mex/sources/libkorder/sylv/SylvParams.cc',
               'mex/sources/libkorder/sylv/SymSchurDecomp.cc',
               'mex/sources/libkorder/sylv/TriangularSylvester.cc',
               'mex/sources/libkorder/sylv/Vector.cc',
               'mex/sources/libkorder/tl/equivalence.cc',
               'mex/sources/libkorder/tl/fine_container.cc',
               'mex/sources/libkorder/tl/fs_tensor.cc',
               'mex/sources/libkorder/tl/gs_tensor.cc',
               'mex/sources/libkorder/tl/int_sequence.cc',
               'mex/sources/libkorder/tl/kron_prod.cc',
               'mex/sources/libkorder/tl/normal_moments.cc',
               'mex/sources/libkorder/tl/permutation.cc',
               'mex/sources/libkorder/tl/ps_tensor.cc',
               'mex/sources/libkorder/tl/pyramid_prod.cc',
               'mex/sources/libkorder/tl/pyramid_prod2.cc',
               'mex/sources/libkorder/tl/rfs_tensor.cc',
               'mex/sources/libkorder/tl/sparse_tensor.cc',
               'mex/sources/libkorder/tl/stack_container.cc',
               'mex/sources/libkorder/tl/symmetry.cc',
               'mex/sources/libkorder/tl/t_container.cc',
               'mex/sources/libkorder/tl/t_polynomial.cc',
               'mex/sources/libkorder/tl/tensor.cc',
               'mex/sources/libkorder/tl/tl_static.cc',
               'mex/sources/libkorder/tl/twod_matrix.cc',
               'mex/sources/libkorder/utils/pascal_triangle.cc',
               'mex/sources/libkorder/utils/int_power.cc',
               'mex/sources/libkorder/utils/sthread.cc',
               'mex/sources/libkorder/k_ord_dynare.cc',
               'mex/sources/libkorder/dynamic_dll.cc',
               'mex/sources/libkorder/dynamic_m.cc' ]
korder_incdir = include_directories('mex/sources/libkorder', 'mex/sources/libkorder/tl', 'mex/sources/libkorder/sylv',
                                    'mex/sources/libkorder/kord', 'mex/sources/libkorder/utils')
korder_lib = static_library('korder', korder_src,
                            kwargs : static_library_kwargs + { 'include_directories' : static_library_kwargs.get('include_directories') + korder_incdir},
                            dependencies : [ blas_dep, lapack_dep, dl_dep ])

korder_mex_kwargs = mex_kwargs + { 'include_directories' : mex_kwargs.get('include_directories') + korder_incdir}
shared_module('gensylv', 'mex/sources/gensylv/gensylv.cc', kwargs : korder_mex_kwargs, link_with : korder_lib)
shared_module('k_order_perturbation', 'mex/sources/k_order_perturbation/k_order_perturbation.cc', kwargs : korder_mex_kwargs, link_with : korder_lib)
k_order_welfare_src = [ 'mex/sources/k_order_welfare/k_order_welfare.cc',
                        'mex/sources/k_order_welfare/approximation_welfare.cc',
                        'mex/sources/k_order_welfare/k_ord_objective.cc',
                        'mex/sources/k_order_welfare/objective_m.cc' ]
shared_module('k_order_welfare', k_order_welfare_src, kwargs : korder_mex_kwargs, link_with : korder_lib)

# Unit tests

korder_test_kwargs = { 'include_directories' : [ mex_incdir, korder_incdir ],
                       'link_args' : exe_link_args,
                       'build_rpath' : exe_rpath,
                       'link_with' : korder_lib }

korder_sylv_test_exe = executable('korder_sylv_test', [ 'mex/sources/libkorder/sylv/tests/MMMatrix.cc',
                                                        'mex/sources/libkorder/sylv/tests/tests.cc'],
                                  kwargs : korder_test_kwargs)
test('korder_sylv', korder_sylv_test_exe, workdir : meson.current_source_dir() / 'mex/sources/libkorder/sylv/tests', timeout : 300)

korder_tl_test_exe = executable('korder_tl_test', [ 'mex/sources/libkorder/tl/tests/factory.cc',
                                                    'mex/sources/libkorder/tl/tests/monoms.cc',
                                                    'mex/sources/libkorder/tl/tests/tests.cc'],
                                kwargs : korder_test_kwargs)
test('korder_tl', korder_tl_test_exe, timeout : 300)

korder_kord_test_exe = executable('korder_kord_test', 'mex/sources/libkorder/kord/tests/tests.cc', kwargs : korder_test_kwargs)
test('korder_kord', korder_kord_test_exe, timeout : 1500)


## MS-SBVAR stuff

ms_sbvar_src = [ 'contrib/ms-sbvar/utilities_dw/arrays/dw_array.c',
                 'contrib/ms-sbvar/utilities_dw/arrays/dw_matrix_array.c',
                 'contrib/ms-sbvar/utilities_dw/ascii/dw_ascii.c',
                 'contrib/ms-sbvar/utilities_dw/ascii/dw_parse_cmd.c',
                 'contrib/ms-sbvar/utilities_dw/elliptical/dw_elliptical.c',
                 'contrib/ms-sbvar/utilities_dw/error/dw_error.c',
                 'contrib/ms-sbvar/utilities_dw/histogram/dw_histogram.c',
                 'contrib/ms-sbvar/utilities_dw/math/dw_math.c',
                 'contrib/ms-sbvar/utilities_dw/matrix/dw_matrix.c',
                 'contrib/ms-sbvar/utilities_dw/matrix/bmatrix.c',
                 'contrib/ms-sbvar/utilities_dw/sort/dw_matrix_sort.c',
                 'contrib/ms-sbvar/utilities_dw/stat/dw_rand_gsl.c',
                 'contrib/ms-sbvar/utilities_dw/stat/dw_matrix_rand.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_switch.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_switchio.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_dirichlet_restrictions.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_metropolis_theta.c',
                 'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARbase.c',
                 'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARio.c',
                 'mex/sources/ms-sbvar/modify_for_mex.cc' ]
ms_sbvar_defs = [ '-DSTRUCTURED_COLUMN_MAJOR' ]
ms_sbvar_incdir = include_directories('contrib/ms-sbvar/utilities_dw/include', 'contrib/ms-sbvar/switch_dw/switching', 'mex/sources/ms-sbvar')
ms_sbvar_lib = static_library('ms_sbvar', ms_sbvar_src,
                              kwargs : static_library_kwargs + { 'c_args' : static_library_kwargs.get('c_args') + ms_sbvar_defs,
                                                                 'cpp_args' : static_library_kwargs.get('cpp_args') + ms_sbvar_defs,
                                                                 'include_directories' : static_library_kwargs.get('include_directories') + ms_sbvar_incdir },
                              dependencies : [ blas_dep, lapack_dep, gsl_dep, matio_dep, ut_dep ])

# NB: in the following, mex_top_level.cc is put as extra source and not in ms_sbvar_lib,
# to ensure that mexFunction() is present in the MEX
mex_ms_sbvar_kwargs = mex_kwargs + { 'c_args' : mex_kwargs.get('c_args') + ms_sbvar_defs,
                                     'cpp_args' : mex_kwargs.get('cpp_args') + ms_sbvar_defs,
                                     'include_directories' : mex_kwargs.get('include_directories') + ms_sbvar_incdir,
                                     'link_with' : ms_sbvar_lib,
                                     'sources' : 'mex/sources/ms-sbvar/mex_top_level.cc' }
ms_sbvar_create_init_file_src = [ 'contrib/ms-sbvar/switch_dw/state_space/sbvar/create_init_file.c',
                                  'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARio_matlab.c' ]
shared_module('ms_sbvar_create_init_file', ms_sbvar_create_init_file_src, kwargs : mex_ms_sbvar_kwargs, dependencies : matio_dep)
ms_sbvar_command_line_src = [ 'contrib/ms-sbvar/switch_dw/switching/dw_switch_opt.c',
                              'contrib/ms-sbvar/switch_dw/switching/dw_mdd_switch.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/dw_sbvar_command_line.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_estimate.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_simulate.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_probabilities.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_mdd.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_forecast.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_variance_decomposition.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_impulse_responses.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/dw_csminwel.c' ]
shared_module('ms_sbvar_command_line', ms_sbvar_command_line_src, kwargs : mex_ms_sbvar_kwargs, dependencies : [ gsl_dep, matio_dep ])


### Integration tests

# Create a test driver (in bash) that takes as arguments:
# - unique name of the test
# - get_option('build_for')
# - the path to the MATLAB/Octave executable: get_option('build_for') == 'matlab' ? matlab_exe.full_path() : octave_exe.full_path()
# - the version of MATLAB/Octave (needed under MATLAB for determining the batch options)
# - a list of .mod files that must be executed sequentially (i.e. a flattened dependency tree)
# - a separator (e.g. "--")
# - a list of extra files on which the .mod files depend (.inc, .mat, …)
#
# The test driver would do the following:
# - create a temporary directory (with mktemp -d or its equivalent on Windows/macOS), of the form dynare-$(testname).XXXXXX
# - copy the .mod files and the extra files into that temporary directory, keeping the original directory structure (a test dependency may be in another directory)
# - set the environment variable DYNARE_BUILD_DIR
# - run MATLAB/Octave on the .mod files sequentially (using a thin .m wrapper for printing the stack trace in case of error)
# - return the correct exitcode
#
# In the present file, the tests list could be organized in a array of arrays. Each inner array would correspond to a test, and contain:
# - the unique name of the test
# - an array containing the list of .mod files
# - an array containing the list of extra files
# - an array indicating the suite(s) to which this test belongs
# - an optional timeout value?
#
# This tests list would be used in a foreach loop around the test() command.
#
# Other items:
# - Decide whether to automatically delete the temporary directories created by the test driver. Cons: removes test data useful for debugging. Pro: avoid excessive disk filling.
# - See what to do with xvfb-run (see #1892). Maybe try to detect it from meson.build, and pass it optionally to the test driver script
# - Include the k-order unit tests in some suite(s)
# - See if the tests can be automatically disabled if the MATLAB/Octave executable is not found (using the disabler option of find_program?)
# - Deal with other types of tests (integration tests implemented as pure .m scripts, unit tests for .m files)
