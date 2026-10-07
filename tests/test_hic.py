import pytest

from deepstructgenomics.analysis.hic import summarize_contacts


def write_contacts(tmp_path, rows):
    path = tmp_path / "contacts.tsv"
    path.write_text("chrom1\tstart1\tchrom2\tstart2\tcount\n" + rows)
    return path


def test_summary_counts_and_diagonal(tmp_path):
    path = write_contacts(tmp_path, "chr1\t0\tchr1\t0\t3\nchr1\t0\tchr1\t10\t5\nchr1\t10\tchr2\t0\t2\n")
    result = summarize_contacts(path, assembly="synthetic", bin_size=10, top_k=1)
    assert result["total_weight"] == 10
    assert result["cis_fraction"] == .8
    assert result["diagonal_weight"] == 3
    assert result["bin_coverage"][0]["weight"] == 8
    assert result["top_contacts"][0]["count"] == 5
    assert len(result["top_contacts"]) == 1


@pytest.mark.parametrize("rows", [
    "", "chr1\t0\tchr1\t10\tnan\n", "chr1\t0\tchr1\t10\t-1\n",
    "chr1\t1\tchr1\t10\t1\n", "chr1\t-10\tchr1\t10\t1\n",
    "chr1\t0\tchr1\t10\t1\nchr1\t10\tchr1\t0\t1\n",
    "chr1\t0\tchr1\t10\n",
])
def test_invalid_contacts(tmp_path, rows):
    with pytest.raises(ValueError):
        summarize_contacts(write_contacts(tmp_path, rows), assembly="synthetic", bin_size=10)


def test_zero_counts(tmp_path):
    result = summarize_contacts(write_contacts(tmp_path, "chr1\t0\tchr1\t10\t0\n"),
                                assembly="synthetic", bin_size=10)
    assert result["cis_fraction"] is None
