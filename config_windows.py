BUILDDIR        = '#build/release'
DISTDIR         = '#dist'
CXX             = 'cl'
CC              = 'cl'
# /O2=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
# EDIT: MSVC2017 bugs Mitsuba when using /O2 because of /Ob2 - replaced it with the full set and /Ob2 for /Ob1
# EDIT: /Og is deprecated too
CXXFLAGS        = ['/nologo', #'/Og',
'/Oi', '/Ot', '/Oy', '/Ob1', #'/Ob2',
'/Gs', '/GF', '/Gy', #'/O2',
'/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'SINGLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D', 'MTS_SSE', '/D', 'MTS_HAS_COHERENT_RT', '/D', '_CONSOLE', '/D', 'NDEBUG', '/D', 'OPENEXR_DLL', '/openmp']
# '/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'DOUBLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D', '_CONSOLE', '/D', 'NDEBUG', '/D', 'OPENEXR_DLL', '/openmp']
SHCXXFLAGS      = CXXFLAGS
TARGET_ARCH     = 'x86_64'
MSVC_VERSION    = '15.0'
LINKFLAGS       = ['/nologo', '/SUBSYSTEM:CONSOLE', '/MACHINE:X64', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT', '/MANIFEST']
BASEINCLUDE     = ['#include', '#dependencies/include']
BASELIB         = ['msvcrt', 'ws2_32', 'Half', 'zlib']
BASELIBDIR      = ['#dependencies/lib']
OEXRINCLUDE     = ['#dependencies/include/openexr']
OEXRLIB         = ['IlmImf', 'IlmThread', 'Iex', 'zlib']
BOOSTLIB        = ['boost_system-vc141-mt-1_64', 'boost_filesystem-vc141-mt-1_64', 'boost_thread-vc141-mt-1_64']
XERCESLIB       = ['xerces-c_3']
PNGLIB          = ['libpng16']
JPEGLIB         = ['jpeg']
GLLIB           = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS         = ['/D', 'GLEW_MX']
SHLIBPREFIX     = 'lib'
SHLIBSUFFIX     = '.dll'
LIBSUFFIX       = '.lib'
PROGSUFFIX      = '.exe'
QTINCLUDE       = ['#dependencies/include']
QTDIR           = '#dependencies'
FFTWLIB         = ['libfftw3-3']
