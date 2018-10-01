from distutils.core import setup, Extension
 
tt_irt1 = Extension('tt_irt1', sources = ['tt_irt1_int32.c'])
 
setup (name = 'tt_irt',
        version = '1.0',
        description = 'Inverse Rosenblatt transform in TT format',
        ext_modules = [tt_irt1],
        py_modules = ['tt_irt'])

