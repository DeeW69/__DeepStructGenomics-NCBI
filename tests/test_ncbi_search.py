"""Bounded NCBI search and preview integrity, without live network access."""
import hashlib
import json

import pytest
import requests

from deepstructgenomics.data_sources.ncbi_client import NCBIClient, normalize_user_sequence
from deepstructgenomics.gui.services import preview_ncbi, run_rna


@pytest.fixture
def responses(monkeypatch):
    calls, queue = [], []

    def request(self, method, url, **kwargs):
        calls.append((url, kwargs["params"]))
        payload = queue.pop(0)
        if isinstance(payload, Exception):
            raise payload
        result = requests.Response()
        result.status_code = 200
        result._content = json.dumps(payload).encode()
        return result

    monkeypatch.setattr(requests.Session, "request", request)
    monkeypatch.setattr("deepstructgenomics.data_sources.ncbi_client.time.sleep", lambda _: None)
    return queue, calls


def test_search_batches_summaries_and_preserves_pagination(responses):
    queue, calls = responses
    queue.extend([{"esearchresult": {"count": "40", "idlist": ["10", "11"]}},
                  {"result": {uid: {"accessionversion": f"NR_{uid}.2", "slen": 123,
                                    "organism": "Homo sapiens", "biomol": "ncRNA", "title": "Example"}
                              for uid in ("10", "11")}}])
    result = NCBIClient().search_sequences("HOTAIR", organism="Homo sapiens", start=20, page_size=2)
    assert result["total"] == 40 and result["start"] == 20
    assert [r["accession"] for r in result["rows"]] == ["NR_10.2", "NR_11.2"]
    assert len(calls) == 2 and calls[1][1]["id"] == "10,11"
    assert calls[0][1]["retstart"] == 20 and calls[0][1]["retmax"] == 2
    assert '"Homo sapiens"[Organism]' in result["query"]
    assert "biomol_ncrna[PROP]" in result["query"]


def test_empty_search_does_not_fetch_summaries(responses):
    queue, calls = responses
    queue.append({"esearchresult": {"count": "0", "idlist": [],
                                   "errorlist": {"phrasesnotfound": ["unknown"], "fieldsnotfound": []}}})
    result = NCBIClient().search_sequences("unknown", rna_only=False)
    assert result["rows"] == [] and len(calls) == 1
    assert result["unmatched_terms"] == ["unknown"] and "biomol" not in result["query"]


@pytest.mark.parametrize("payload", [[], {"error": "failed"}, {"esearchresult": {}},
    {"esearchresult": {"count": "bad", "idlist": []}},
    {"esearchresult": {"count": "1", "idlist": [None]}},
    {"esearchresult": {"errorlist": {"fieldsnotfound": ["bad"]}}}])
def test_malformed_search_is_recoverable(responses, payload):
    responses[0].append(payload)
    with pytest.raises(ValueError):
        NCBIClient().search_sequences("BRCA1")


def test_timeout_is_not_an_empty_result(responses):
    responses[0].append(requests.Timeout())
    with pytest.raises(requests.Timeout):
        NCBIClient().search_sequences("BRCA1")


def test_missing_summary_is_not_silently_dropped(responses):
    responses[0].extend([{"esearchresult": {"count": "1", "idlist": ["1"]}}, {"result": {}}])
    with pytest.raises(ValueError, match="incompl"):
        NCBIClient().search_sequences("BRCA1")


@pytest.mark.parametrize("parameters", [{"term": " "}, {"term": "x", "page_size": 51},
    {"term": "x", "start": -1}, {"term": "x", "organism": 'bad"[all]'}])
def test_invalid_input_does_not_access_network(responses, parameters):
    with pytest.raises(ValueError):
        NCBIClient().search_sequences(**parameters)
    assert responses[1] == []


def test_preview_cache_and_analysis_keep_ncbi_provenance(tmp_path, responses):
    queue, calls = responses
    # FASTA is returned by a separate endpoint; use real client caching below.
    queue.extend([{"result": {"1": {"title": "Reference"}}}])
    client_method = NCBIClient._call_endpoint
    from unittest.mock import patch

    def endpoint(self, endpoint, parameters):
        if endpoint == "efetch.fcgi":
            return ">NR_TEST.1 Reference\nGGGGAAAACCCC\n"
        return client_method(self, endpoint, parameters)

    with patch.object(NCBIClient, "_call_endpoint", endpoint):
        preview = preview_ncbi(tmp_path, {"accession": "NR_TEST.1"})
        assert preview["metadata"]["provenance"]["source"] == "NCBI"
        assert preview["sha256"] == hashlib.sha256(b"GGGGAAAACCCC").hexdigest()
        entry = run_rna(tmp_path, {"accession": preview["accession"], "expected_ncbi_sha256": preview["sha256"]})
    assert entry["label"] == "NR_TEST.1" and len(calls) == 1


def test_changed_preview_stops_before_prediction(tmp_path, monkeypatch):
    monkeypatch.setattr(NCBIClient, "fetch_sequence", lambda *args, **kwargs: normalize_user_sequence("ACGU", "test"))
    monkeypatch.setattr("deepstructgenomics.gui.services.DeepStructPipeline.run_and_export",
                        lambda *args: pytest.fail("Changed reference must not be analyzed"))
    with pytest.raises(ValueError, match="depuis"):
        run_rna(tmp_path, {"accession": "NR_TEST.1", "expected_ncbi_sha256": "different"})
