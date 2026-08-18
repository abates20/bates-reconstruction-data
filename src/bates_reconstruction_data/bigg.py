"""
bigg
----

Module for using the BiGG API.
"""


from typing import Optional

import requests

from bates_reconstruction_data.utils import limit


BASE_URL = "http://bigg.ucsd.edu/api/v2"


@limit.rate_limit(10, 1) # BiGG asks that API calls be limited to 10 requests per second
def _call_api(url, session: Optional[requests.Session] = None) -> requests.Response:
    if session is None:
        return requests.get(url)
    return session.get(url)


def get_model_metadata(model_id: Optional[str] = None, session: Optional[requests.Session] = None):
    url = f"{BASE_URL}/models" if model_id is None else f"{BASE_URL}/models/{model_id}"
    response = _call_api(url, session)
    return response.json()

def get_model(model_id: str, session: Optional[requests.Session] = None):
    url = f"{BASE_URL}/models/{model_id}/download"
    response = _call_api(url, session)
    return response.json()

def get_gene_details(model_id: str, gene_id: str, session: Optional[requests.Session] = None):
    url = f"{BASE_URL}/models/{model_id}/genes/{gene_id}"
    response = _call_api(url, session)
    return response.json()

def get_genes(model_id: str, session: Optional[requests.Session] = None):
    url = f"{BASE_URL}/models/{model_id}/genes"
    response = _call_api(url, session)
    return response.json()

def get_universal_metabolite(met_id: str, session: Optional[requests.Session] = None):
    url = f"{BASE_URL}/universal/metabolites/{met_id}"
    response = _call_api(url, session)
    return response.json()
