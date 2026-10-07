import json
from unittest.mock import Mock

import pytest

from deepstructgenomics.data_sources.ncbi_client import NCBIClient


def client(path):
    result = NCBIClient(cache_dir=path)
    result._call_endpoint = Mock(return_value=">test description\nACGT")
    result._fetch_summary = Mock(return_value={"title": "test"})
    return result


def test_cache_hit_refresh_and_provenance(tmp_path):
    first = client(tmp_path)
    record = first.fetch_sequence("test")
    assert record.sequence == "ACGU"
    assert not record.metadata["provenance"]["cache_hit"]
    second = client(tmp_path)
    cached = second.fetch_sequence("test")
    assert cached.metadata["provenance"]["cache_hit"]
    assert cached.metadata["provenance"]["retrieved_at"] == record.metadata["provenance"]["retrieved_at"]
    second._call_endpoint.assert_not_called()
    second.fetch_sequence("test", refresh=True)
    second._call_endpoint.assert_called_once()


def test_corruption_and_failed_refresh(tmp_path):
    ncbi = client(tmp_path)
    ncbi.fetch_sequence("../../test")
    path = next(tmp_path.glob("*.json"))
    payload = json.loads(path.read_text())
    payload["record"]["sequence"] = "AAAA"
    path.write_text(json.dumps(payload))
    assert ncbi.fetch_sequence("../../test").sequence == "ACGU"
    before = path.read_bytes()
    ncbi._call_endpoint.side_effect = ValueError("failed")
    with pytest.raises(ValueError):
        ncbi.fetch_sequence("../../test", refresh=True)
    assert path.read_bytes() == before
