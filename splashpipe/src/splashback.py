# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_splashback')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_splashback')
    _splashback = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_splashback', [dirname(__file__)])
        except ImportError:
            import _splashback
            return _splashback
        try:
            _mod = imp.load_module('_splashback', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _splashback = swig_import_helper()
    del swig_import_helper
else:
    import _splashback
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class dp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, dp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, dp, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _splashback.new_dp()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _splashback.delete_dp
    __del__ = lambda self: None

    def assign(self, value):
        return _splashback.dp_assign(self, value)

    def value(self):
        return _splashback.dp_value(self)

    def cast(self):
        return _splashback.dp_cast(self)
    if _newclass:
        frompointer = staticmethod(_splashback.dp_frompointer)
    else:
        frompointer = _splashback.dp_frompointer
dp_swigregister = _splashback.dp_swigregister
dp_swigregister(dp)

def dp_frompointer(t):
    return _splashback.dp_frompointer(t)
dp_frompointer = _splashback.dp_frompointer

class splashback(_object):
    """Proxy of C++ splashback class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, splashback, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, splashback, name)
    __repr__ = _swig_repr
    __swig_setmethods__["tree"] = _splashback.splashback_tree_set
    __swig_getmethods__["tree"] = _splashback.splashback_tree_get
    if _newclass:
        tree = _swig_property(_splashback.splashback_tree_get, _splashback.splashback_tree_set)

    def __init__(self, *args, **kwargs):
        """__init__(splashback self, float xrmin=0.5, float xrmax=15.0, int xrbins=15, char * xoutfile, float xmag_limit=25.0, float zmax=0.5, int Njack=25, bool deproject=False, int colored=0, bool verbose=False) -> splashback"""
        this = _splashback.new_splashback(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _splashback.delete_splashback
    __del__ = lambda self: None

    def allocate_lens_memory(self, xNcen):
        """allocate_lens_memory(splashback self, int xNcen) -> int"""
        return _splashback.splashback_allocate_lens_memory(self, xNcen)


    def process_lens(self, xra, xdec, xzred, xjackreg, wt=1.0):
        """process_lens(splashback self, float xra, float xdec, float xzred, int xjackreg, float wt=1.0) -> int"""
        return _splashback.splashback_process_lens(self, xra, xdec, xzred, xjackreg, wt)


    def finalize_lenses(self):
        """finalize_lenses(splashback self) -> int"""
        return _splashback.splashback_finalize_lenses(self)


    def process_source(self, sra, sdec, smag, disable_magcheck, color=1.E30):
        """process_source(splashback self, float sra, float sdec, float smag, bool disable_magcheck, float color=1.E30) -> int"""
        return _splashback.splashback_process_source(self, sra, sdec, smag, disable_magcheck, color)


    def finalize_results(self, writeok=False):
        """finalize_results(splashback self, bool writeok=False) -> int"""
        return _splashback.splashback_finalize_results(self, writeok)


    def test_searchrecord(self):
        """test_searchrecord(splashback self) -> int"""
        return _splashback.splashback_test_searchrecord(self)


    def deprojection_kernel(self, x):
        """deprojection_kernel(splashback self, double x) -> double"""
        return _splashback.splashback_deprojection_kernel(self, x)

splashback_swigregister = _splashback.splashback_swigregister
splashback_swigregister(splashback)

# This file is compatible with both classic and new-style classes.


