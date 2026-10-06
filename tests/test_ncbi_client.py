"""Exercise NCBI response handling without opening a network connection."""

from __future__ import annotations

import json
from unittest.mock import Mock, call

import pytest
import requests

from deepstructgenomics.config import NCBIConfig
from deepstructgenomics.data_sources.ncbi_client import NCBIClient


ACCESSION = "NC_TEST.1"
FASTA = f">{ACCESSION} DNA reference\r\nat gc\r\n"


def response(body: str, status: int = 200) -> requests.Response:
    result = requests.Response()
    result.status_code = status
    result.encoding = "utf-8"
    result._content = body.encode("utf-8")
    return result


@pytest.fixture(autouse=True)
def block_network(monkeypatch):
    def reject_network(*args, **kwargs):
        pytest.fail("NCBI client tests must not make real HTTP requests")

    monkeypatch.setattr(requests.Session, "request", reject_network)
    monkeypatch.setattr(requests.adapters.HTTPAdapter, "send", reject_network)


@pytest.fixture
def mock_ncbi(monkeypatch):
    def install(*replies):
        pending = iter(replies)

        def request(method, url, **kwargs):
            reply = next(pending, None)
            if reply is None:
                pytest.fail(f"Unexpected additional request: {method} {url}")
            if isinstance(reply, Exception):
                raise reply
            reply.url = url
            reply.request = requests.Request(method, url, params=kwargs.get("params")).prepare()
            return reply

        mocked = Mock(side_effect=request)
        monkeypatch.setattr(requests.Session, "request", mocked)
        return mocked

    return install


def test_fetch_preserves_ncbi_metadata_and_sends_config_to_both_endpoints(mock_ncbi):
    summary = {
        "uid": "42",
        "accessionversion": ACCESSION,
        "title": "Metadata title",
        "organism": "Test organism",
    }
    mocked = mock_ncbi(
        response(FASTA),
        response(json.dumps({"result": {"uids": ["42"], "42": summary}})),
    )
    config = NCBIConfig(
        base_url="https://ncbi.example.test/eutils",
        database="nucleotide",
        tool_name="offline-test",
        email="maintainer@example.test",
        api_key="test-key",
        timeout=7,
    )

    record = NCBIClient(config).fetch_sequence(ACCESSION)

    assert record.identifier == ACCESSION
    assert record.description == "DNA reference"
    assert record.sequence == "AUGC"
    assert record.metadata == summary
    common = {
        "db": "nucleotide",
        "id": ACCESSION,
        "tool": "offline-test",
        "email": "maintainer@example.test",
        "api_key": "test-key",
    }
    assert mocked.call_args_list == [
        call(
            "GET", f"{config.base_url}/efetch.fcgi",
            params={**common, "rettype": "fasta", "retmode": "text"},
            timeout=7, allow_redirects=True,
        ),
        call(
            "GET", f"{config.base_url}/esummary.fcgi",
            params={**common, "retmode": "json"},
            timeout=7, allow_redirects=True,
        ),
    ]


def test_missing_description_uses_summary_title_and_optional_credentials_are_omitted(mock_ncbi):
    mocked = mock_ncbi(
        response(f">{ACCESSION}\nATGC\n"),
        response(json.dumps({"result": {"uids": ["42"], "42": {"title": "Summary title"}}})),
    )

    record = NCBIClient().fetch_sequence(ACCESSION)

    assert record.description == "Summary title"
    assert mocked.call_count == 2
    for recorded in mocked.call_args_list:
        assert "email" not in recorded.kwargs["params"]
        assert "api_key" not in recorded.kwargs["params"]
        assert recorded.kwargs["timeout"] == NCBIConfig().timeout


@pytest.mark.parametrize(
    "body",
    ["", " \r\n\t", f">{ACCESSION}\n", "NCBI error: record unavailable", f">{ACCESSION}\nATGN\n"],
    ids=["empty", "whitespace", "no-bases", "not-fasta", "ambiguous-base"],
)
def test_invalid_fasta_stops_before_summary_request(mock_ncbi, body):
    mocked = mock_ncbi(response(body))

    with pytest.raises(ValueError):
        NCBIClient().fetch_sequence(ACCESSION)

    assert mocked.call_count == 1
    assert mocked.call_args.args[1].endswith("/efetch.fcgi")


@pytest.mark.parametrize(
    "body",
    [
        "", " \r\n", "<html>temporarily unavailable</html>",
        "null", "[]", "{}", '{"result": null}', '{"result": []}',
        '{"result": {}}', '{"result": {"uids": []}}',
        '{"result": {"uids": ["42"], "42": null}}',
        '{"result": {"uids": ["42"], "42": []}}',
    ],
    ids=[
        "empty", "whitespace", "invalid-json", "null-payload", "list-payload", "missing-result",
        "null-result", "list-result", "empty-result", "empty-uids", "null-summary", "list-summary",
    ],
)
def test_unusable_summary_keeps_sequence_with_empty_metadata(mock_ncbi, body):
    mocked = mock_ncbi(response(FASTA), response(body))

    record = NCBIClient().fetch_sequence(ACCESSION)

    assert record.identifier == ACCESSION
    assert record.description == "DNA reference"
    assert record.sequence == "AUGC"
    assert record.metadata == {}
    assert mocked.call_count == 2


@pytest.mark.parametrize("endpoint", ["efetch", "esummary"])
@pytest.mark.parametrize("status", [400, 404, 429, 500])
def test_http_failure_is_propagated_without_retry_or_later_request(mock_ncbi, endpoint, status):
    # Bodies are deliberately invalid: status errors must precede parsing.
    failed = response("upstream error", status=status)
    replies = [response(FASTA), failed] if endpoint == "esummary" else [failed]
    mocked = mock_ncbi(*replies)

    with pytest.raises(requests.HTTPError) as error:
        NCBIClient().fetch_sequence(ACCESSION)

    assert error.value.response is failed
    assert error.value.response.status_code == status
    assert error.value.response.url.endswith(f"/{endpoint}.fcgi")
    assert mocked.call_count == len(replies)


@pytest.mark.parametrize("endpoint", ["efetch", "esummary"])
def test_timeout_is_propagated_without_retry_or_later_request(mock_ncbi, endpoint):
    timeout = requests.Timeout("NCBI did not respond in time")
    replies = [response(FASTA), timeout] if endpoint == "esummary" else [timeout]
    mocked = mock_ncbi(*replies)

    with pytest.raises(requests.Timeout) as error:
        NCBIClient().fetch_sequence(ACCESSION)

    assert error.value is timeout
    assert mocked.call_count == len(replies)
    assert mocked.call_args.args[1].endswith(f"/{endpoint}.fcgi")
