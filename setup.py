import distutils.core as dc
import distutils.extension as dex
from Cython.Distutils import build_ext

dc.setup(
    name = 'ljlib',
    ext_modules=[
        dex.Extension('ljlib',
                      sources=['ljlib.pyx', 'ljforces.cpp'],
                      extra_compile_args=['-std=c++11'],
                      language='c++')
        ],
    cmdclass = {'build_ext': build_ext}
)
