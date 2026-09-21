from typing import Optional, TypeAlias, Callable, Any
import logging
import itertools
import re
import functools

import requests

from bates_reconstruction_data.utils import limit, parallel


_Parser: TypeAlias = Callable[[dict[str, Any], str, str], dict[str, Any]]
_SKIP: _Parser = lambda dct, field, value: dct

# API rate limit values
# (KEGG limit is technically 3 per second but doing 2 per second to be safe)
KEGG_N = 3
KEGG_T = 1.5


@limit.rate_limit(KEGG_N, KEGG_T)
def request(url: str, method: str = "GET", session: Optional[requests.Session] = None) -> requests.Response:
    """
    Sends a request to the provided url, rate-limited to comply with KEGG's
    policy (no more than 3 requests per second).
    """
    logging.debug(f"KEGG API call - {url}")
    if session is not None:
        return session.request(method, url)
    else:
        return requests.request(method, url)


def batched_request(
    url_template: str,
    iterable,
    batch_size: int = 10,
    method: str = "GET",
    threads: int = 1,
    session: Optional[requests.Session] = None,
    show_progress: bool = True
) -> list[requests.Response]:
    if session is None:
        session = requests.Session()
        cleanup = True
    else:
        cleanup = False

    results: list[requests.Response] = parallel.apply(
        request,
        [url_template.format("+".join(b)) for b in batched(iterable, batch_size)],
        itertools.repeat(method),
        itertools.repeat(session),
        threads=threads,
        show_progress=show_progress
    )

    if cleanup:
        session.close()

    return results


def batched(iterable, n):
    """
    Yield batches of size `n` from an iterable.
    """
    it = iter(iterable)
    while True:
        batch = tuple(itertools.islice(it, n))
        if not batch:
            return
        yield batch


class FlatFileReader:

    EOF_SYMBOL = "///"
    _parsers: dict[str, _Parser] = {}

    @classmethod
    def parse(cls, flatfile: str, default_parser: _Parser = _SKIP) -> list[dict]:
        entries = []
        current_entry: dict[str, Any] = {}
        current_field = None
        parser = default_parser

        for line in flatfile.splitlines():
            if line == cls.EOF_SYMBOL:
                entries.append(current_entry)

                # Reset everything for new entry
                current_entry = {}
                current_field = None
                parser = default_parser
                continue

            if line == "":
                continue # skip empty lines

            field_name, rest = re.split(r"\s{2,12}", line, maxsplit=1)
            if field_name != "":
                current_field = field_name
                parser = cls.get_parser(field_name, default_parser)

            current_entry = parser(current_entry, current_field, rest) # ty:ignore[invalid-argument-type]
        return entries

    @classmethod
    def register(cls, name: str):
        def decorator(func: _Parser):

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)

            cls._parsers[name] = wrapper

            return wrapper
        return decorator

    @classmethod
    def get_parser(cls, name: str, default: _Parser) -> _Parser:
        return cls._parsers.get(name, default)


@FlatFileReader.register("ENTRY")
def parse_entry(fields_dct: dict[str, Any], field: str, value: str):
    entry_id, *_, database = value.split()
    fields_dct[field] = {"entry_id": entry_id, "database": database}
    return fields_dct


@FlatFileReader.register("SYMBOL")
def parse_symbol(fields_dct: dict[str, Any], field: str, value: str):
    fields_dct[field] = value
    return fields_dct


@FlatFileReader.register("NAME")
def parse_name(fields_dct: dict[str, Any], field: str, value: str):
    fields_dct[field] = value
    return fields_dct


@FlatFileReader.register("PATHWAY")
def parse_pathway(fields_dct: dict[str, Any], field: str, value: str):
    pathway_list: list[dict[str, str]] = fields_dct.get(field, [])

    pattern = re.compile(r"(map\d{5})\s+(.*)")
    m = pattern.match(value)
    if m is None:
        raise Exception(f"Error parsing KEGG pathway: {value}")

    pid, pname = m.groups()
    pathway_list.append({"pathway_id": pid, "name": pname})
    fields_dct[field] = pathway_list
    return fields_dct


@FlatFileReader.register("MODULE")
def parse_reaction(fields_dct: dict[str, Any], field: str, value: str):
    mod_list: list[dict[str, str]] = fields_dct.get(field, [])

    pattern = re.compile(r"(M\d{5})\s+(.*)")
    m = pattern.match(value)
    if m is None:
        raise Exception(f"Error parsing KEGG module: {value}")

    mid, mname = m.groups()
    mod_list.append({"module_id": mid, "name": mname})
    fields_dct[field] = mod_list
    return fields_dct


@FlatFileReader.register("REACTION")
def parse_reaction(fields_dct: dict[str, Any], field: str, value: str):
    rxn_list: list[dict[str, str]] = fields_dct.get(field, [])

    pattern = re.compile(r"(R\d{5})\s+(.*)")
    m = pattern.match(value)
    if m is None:
        raise Exception(f"Error parsing KEGG reaction: {value}")
    rid, rname = m.groups()
    rxn_list.append({"reaction_id": rid, "name": rname})
    fields_dct[field] = rxn_list
    return fields_dct


@FlatFileReader.register("BRITE")
def parse_brite(fields_dct: dict[str, Any], field: str, value: str):
    brite_list: list[str] = fields_dct.get(field, [])
    brite_list.append(value)
    fields_dct[field] = brite_list
    return fields_dct


@FlatFileReader.register("GENES")
def parse_gene(fields_dct: dict[str, Any], field: str, value: str):
    gene_list: list[dict] = fields_dct.get(field, [])

    organism, gene_ids = re.split(r":\s+", value)
    gene_dct = {"org_id": organism, "org_genes": []}
    for gene in gene_ids.split():
        symbol_pattern = re.compile(r"(\(.+\))")
        gene_id = f"{organism.lower()}:{symbol_pattern.sub('', gene)}"
        gene_name = ""
        m = symbol_pattern.search(gene)
        if m:
            gene_name = m.group().strip("()")
        gene_dct["org_genes"].append({"gene_id": gene_id, "symbol": gene_name})

    gene_list.append(gene_dct)
    fields_dct[field] = gene_list
    return fields_dct


@FlatFileReader.register("REFERENCE")
def parse_reference(fields_dct: dict[str, Any], field: str, value: str):
    if value.startswith("PMID:"):
        ...
    return fields_dct
