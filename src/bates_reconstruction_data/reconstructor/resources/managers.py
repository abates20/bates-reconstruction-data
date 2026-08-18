from typing import Generic, TypeVar, Union, Callable, Optional
import os
from pathlib import Path
from functools import wraps

import cobra

from ..resources.files import get_resource_path
from ..utils import load_json


T = TypeVar("T")
U = TypeVar("U")
PathLike = Union[str, bytes, os.PathLike]
FileLoader = Callable[[PathLike], T]


class SizedDict(dict):

    def __init__(self, maxsize: Optional[int] = None):
        if maxsize is not None and maxsize < 1:
            raise ValueError("maxsize parameter in SizedDict must be >= 1.")
        super().__init__()
        self._size = maxsize

    def __setitem__(self, key, value):
        if self._size:
            while len(self) >= self._size:
                oldest_key = next(iter(self))
                super().__delitem__(oldest_key)
        super().__setitem__(key, value)


def _cache(maxsize: int = 1, try_copy: bool = True):

    def decorator(func):
        cache = SizedDict(maxsize)

        @wraps(func)
        def wrapper(*args, **kwargs):
            key = (args, tuple(sorted(kwargs.items())))
            if key in cache:
                result = cache[key]
            else:
                result = func(*args, **kwargs)
            cache[key] = result

            if try_copy and hasattr(result, "copy"):
                result = result.copy()
            
            return result
        
        return wrapper
    return decorator


class Resource(Generic[T, U]):

    def __init__(
        self,
        value: Union[T, U],
        *,
        load_hook: Optional[Callable[[U], T]] = None,
    ):
        self._value = value
        self._hook = load_hook

    def load(self) -> T:
        return self._hook(self._value) if self._hook else self._value


class ResourceManager(Generic[T]):

    def __init__(
        self,
        resource_name: str = "resource",
        *,
        load_hook: Optional[Callable[[PathLike], T]] = None
    ):
        self._dict: dict[str, Resource[T, PathLike]] = {}
        self._hook = load_hook
        self._resource_name = resource_name

    @property
    def registered_names(self) -> list[str]:
        return list(self._dict)

    def register(self, key: str, resource: T):
        self._dict[key] = Resource(resource)

    def register_file(self, key: str, resource_path: PathLike):
        self._dict[key] = Resource(resource_path, load_hook=self._hook)

    @_cache(maxsize=1, try_copy=True)
    def load(self, key: str) -> T:
        if key not in self._dict:
            options = ", ".join(repr(x) for x in self.registered_names)
            raise KeyError(
                f"{repr(key)} is not a registered {self._resource_name}. Available options are {options}."
            )
        return self._dict[key].load()


# Create a manager for universal models and register the models
universal_manager: ResourceManager[cobra.Model] = ResourceManager(
    "universal model",
    load_hook=cobra.io.read_sbml_model
)
universal_manager.register_file(
    "universal_prokaryote",
    get_resource_path("universal_prokaryote.sbml.gz")
)


# Create a manager for template models and register the models
template_manager: ResourceManager[cobra.Model] = ResourceManager(
    "template model",
    load_hook=cobra.io.read_sbml_model
)
template_manager.register_file(
    "gram_negative",
    get_resource_path("gram_neg_template.sbml")
)
template_manager.register_file(
    "gram_positive",
    get_resource_path("gram_pos_template.sbml")
)


# Create a manager for the DIAMOND databases and register any databases
diamond_db_manager: ResourceManager[Path] = ResourceManager("DIAMOND database")
diamond_db_manager.register(
    "kegg_prokaryote",
    get_resource_path("screened_kegg_prokaryotes_pep_db.dmnd")
)


# Create a manager for gene -> rxn mappings
gene_rxn_map_manager: ResourceManager[dict[str, list[str]]] = ResourceManager(
    "gene-to-reaction map",
    load_hook=load_json
)
gene_rxn_map_manager.register_file(
    "kegg_seed",
    get_resource_path("gene_modelseed.json.gz")
)


# Create a manager for gene -> name mappings
gene_name_manager: ResourceManager[dict[str, str]] = ResourceManager(
    "gene-name map",
    load_hook=load_json
)
gene_name_manager.register_file(
    "kegg",
    get_resource_path("gene_names.json.gz")
)
