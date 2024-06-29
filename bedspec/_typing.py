from types import UnionType
from typing import Callable
from typing import Type
from typing import Union
from typing import _BaseGenericAlias  # type: ignore[attr-defined]
from typing import get_args
from typing import get_origin


class MethodType:
    def __init__(self, func: Callable, obj: object) -> None:
        self.__func__ = func
        self.__self__ = obj

    def __call__(self, *args: object, **kwargs: object) -> object:
        func = self.__func__
        obj = self.__self__
        return func(obj, *args, **kwargs)


class classmethod_generic:
    def __init__(self, f: Callable) -> None:
        self.f = f

    def __get__(self, obj: object, cls: object | None = None) -> Callable:
        if cls is None:
            cls = type(obj)
        method = MethodType(self.f, cls)
        method._generic_classmethod = True  # type: ignore[attr-defined]
        return method


def __getattr__(self: object, name: str | None = None) -> object:
    if hasattr(obj := orig_getattr(self, name), "_generic_classmethod"):
        obj.__self__ = self
    return obj


orig_getattr = _BaseGenericAlias.__getattr__
_BaseGenericAlias.__getattr__ = __getattr__


def is_union(annotation: Type) -> bool:
    """Test if we have a union type annotation or not."""
    return get_origin(annotation) in {Union, UnionType}


def is_optional(annotation: Type) -> bool:
    """Return if this type annotation is optional (a union type with None) or not."""
    return is_union(annotation) and type(None) in get_args(annotation)


def singular_non_optional_type(annotation: Type) -> Type:
    """Return the non-optional version of a singular type annotation."""
    if not is_optional(annotation):
        return annotation

    not_none: list[Type] = [arg for arg in get_args(annotation) if arg is not type(None)]
    if len(not_none) == 1:
        return not_none[0]
    else:
        raise TypeError(f"Complex non-optional types are not supported! Found: {not_none}")
