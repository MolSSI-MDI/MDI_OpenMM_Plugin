from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
openmmmdi_header_dir = '@EXAMPLEPLUGIN_HEADER_DIR@'
openmmmdi_library_dir = '@EXAMPLEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

class CMakeExtension(Extension):
    def __init__(self, name, source_dir=''):
        Extension.__init__(self, name, sources=[])
        self.source_dir = os.path.abspath(source_dir)

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='openmmmdi._openmmmdi',
                      sources=['openmmmdi/ExamplePluginWrapper.cpp'],
                      libraries=['OpenMM', 'ExamplePlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), os.path.join(os.path.join(openmm_dir, 'include'),'MDI_Library'), openmmmdi_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), os.path.join(os.path.join(openmm_dir, 'lib'),'MDI_Library'), openmmmdi_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='openmmmdi',
      version='1.0',
      packages=['openmmmdi'],
      ext_modules=[CMakeExtension('mdi'), extension],
     )
