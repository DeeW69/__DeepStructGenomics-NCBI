"""Isolated compute process. Requests travel over stdin, never shell arguments."""
import json
import sys


def main():
    try:
        from .services import run_hic, run_rna, search_ncbi, preview_ncbi
        request = json.load(sys.stdin)
        operation = {"rna": run_rna, "hic": run_hic, "ncbi_search": search_ncbi,
                     "ncbi_preview": preview_ncbi}[request["kind"]]
        result = operation(request["workspace"], request["parameters"])
        print(json.dumps(result, ensure_ascii=True))
        return 0
    except Exception as exc:
        import requests
        message = ("NCBI inaccessible : vérifier l'accession et la connexion."
                   if isinstance(exc, requests.RequestException) else str(exc))
        print(json.dumps({"error": message}, ensure_ascii=True))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
