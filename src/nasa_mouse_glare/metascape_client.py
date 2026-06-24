"""Small Metascape web-client wrapper for NASA mouse GLARE gene lists.

Metascape does not document these web-app endpoints as a stable public API.
This module mirrors the browser workflow so repeated GLARE analyses can be run
without manual uploads, while keeping the integration isolated and easy to
replace if Metascape changes its frontend contract.
"""

from __future__ import annotations

import argparse
import http.cookiejar
import json
import mimetypes
import time
import uuid
from pathlib import Path
from typing import Any
from urllib import parse, request
from urllib.error import HTTPError, URLError


DEFAULT_BASE_URL = "https://metascape.org/gp_server"
DEFAULT_ANALYSIS_SPECIES = 10090
DEFAULT_INPUT_SPECIES = 10090
DEFAULT_REPORT_FILES = [
    "metascape_result.xlsx",
    "Enrichment_GO/_FINAL_GO.csv",
    "Enrichment_GO/GO_AllLists.csv",
    "Enrichment_heatmap/HeatmapSelectedGO.csv",
    "Enrichment_PPI/GO_MCODE.csv",
    "Enrichment_PPI/_FINAL_MCODE.csv",
    "Enrichment_PPI/MCODE.csv",
    "Enrichment_QC/GO_PaGenBase.csv",
    "Enrichment_QC/GO_TRRUST.csv",
]


class MetascapeError(RuntimeError):
    """Raised when the Metascape web client receives an unusable response."""


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def read_ids(path: str | Path) -> list[str]:
    text = Path(path).read_text(encoding="utf-8")
    return [item for item in _split_ids(text) if item]


def _split_ids(text: str) -> list[str]:
    normalized = text.replace(",", "\n").replace(";", "\n").replace("\t", "\n")
    return [line.strip() for line in normalized.splitlines() if line.strip()]


def _json_loads(value: Any) -> Any:
    if isinstance(value, str):
        return json.loads(value)
    return value


def _dedupe(values: list[Any]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value).strip()
        if text and text not in seen:
            result.append(text)
            seen.add(text)
    return result


def selected_term_ids(tree: list[dict[str, Any]]) -> list[Any]:
    ids: list[Any] = []

    def visit(node: dict[str, Any]) -> None:
        children = node.get("items") or node.get("children") or []
        if node.get("checked") and not children and "id" in node:
            ids.append(node["id"])
        for child in children:
            if isinstance(child, dict):
                visit(child)

    for root in tree:
        visit(root)
    return ids


class MetascapeClient:
    def __init__(self, base_url: str = DEFAULT_BASE_URL, timeout: int = 120):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        cookie_jar = http.cookiejar.CookieJar()
        self.opener = request.build_opener(request.HTTPCookieProcessor(cookie_jar))

    def _url(self, path: str, query: dict[str, Any] | None = None) -> str:
        url = f"{self.base_url}/{path.lstrip('/')}"
        if query:
            url += "?" + parse.urlencode(query)
        return url

    def _open(self, req: request.Request) -> bytes:
        try:
            with self.opener.open(req, timeout=self.timeout) as response:
                return response.read()
        except HTTPError as exc:
            body = exc.read().decode("utf-8", "replace")
            raise MetascapeError(
                f"Metascape HTTP {exc.code} for {req.full_url}: {body[:1000]}"
            ) from exc
        except URLError as exc:
            raise MetascapeError(f"Metascape request failed: {exc}") from exc

    def get_text(self, path: str, query: dict[str, Any] | None = None) -> str:
        req = request.Request(self._url(path, query), method="GET")
        return self._open(req).decode("utf-8", "replace")

    def get_json(self, path: str, query: dict[str, Any] | None = None) -> Any:
        return json.loads(self.get_text(path, query))

    def post_json(self, path: str, payload: dict[str, Any]) -> Any:
        data = json.dumps(payload).encode("utf-8")
        req = request.Request(
            self._url(path),
            data=data,
            method="POST",
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
        )
        text = self._open(req).decode("utf-8", "replace")
        return json.loads(text) if text.strip() else {}

    def post_multipart(
        self,
        path: str,
        fields: dict[str, Any],
        file_field: str,
        file_path: str | Path,
    ) -> Any:
        file_path = Path(file_path)
        boundary = f"----nasaMouseMetascape{uuid.uuid4().hex}"
        chunks: list[bytes] = []
        for key, value in fields.items():
            chunks.extend(
                [
                    f"--{boundary}\r\n".encode("utf-8"),
                    (
                        f'Content-Disposition: form-data; name="{key}"\r\n\r\n'
                    ).encode("utf-8"),
                    str(value).encode("utf-8"),
                    b"\r\n",
                ]
            )
        content_type = mimetypes.guess_type(file_path.name)[0] or "text/plain"
        chunks.extend(
            [
                f"--{boundary}\r\n".encode("utf-8"),
                (
                    f'Content-Disposition: form-data; name="{file_field}"; '
                    f'filename="{file_path.name}"\r\n'
                ).encode("utf-8"),
                f"Content-Type: {content_type}\r\n\r\n".encode("utf-8"),
                file_path.read_bytes(),
                b"\r\n",
                f"--{boundary}--\r\n".encode("utf-8"),
            ]
        )
        body = b"".join(chunks)
        req = request.Request(
            self._url(path),
            data=body,
            method="POST",
            headers={
                "Content-Type": f"multipart/form-data; boundary={boundary}",
                "Accept": "application/json",
            },
        )
        text = self._open(req).decode("utf-8", "replace")
        return json.loads(text) if text.strip() else {}

    def create_session(self) -> str:
        session_id = self.get_text(
            "get_session_id",
            {"old_session_id": "null", "idsurl": "undefined"},
        ).strip().strip('"')
        if not session_id:
            raise MetascapeError("Metascape did not return a session ID")
        return session_id

    def upload_gene_lists(
        self,
        session_id: str,
        gene_lists_csv: str | Path,
        multiple_list: bool = True,
    ) -> dict[str, Any]:
        return self.post_multipart(
            "upload_excel",
            {
                "session_id": session_id,
                "isMultipleList": str(bool(multiple_list)).lower(),
            },
            "files",
            gene_lists_csv,
        )

    def apply_species(
        self,
        session_id: str,
        input_species: int,
        analysis_species: int,
        multiple_list: bool = True,
        input_column: str = "Gene",
    ) -> dict[str, Any]:
        return self.post_json(
            "apply_species",
            {
                "specifiedSpeciesOption": {
                    "session_id": session_id,
                    "multipleList": {"value": bool(multiple_list)},
                    "input_column": input_column,
                    "specifiedSpecies": {
                        "analysisSpecies": int(analysis_species),
                        "inputSpecies": int(input_species),
                    },
                }
            },
        )

    def get_color_list(self, session_id: str) -> list[str]:
        result = self.get_json("get_enrichment_color_list", {"session_id": session_id})
        return list(result.get("color_list", []))

    def convert_background(
        self,
        session_id: str,
        background_ids: list[str],
    ) -> dict[str, Any]:
        return self.post_json(
            "convert_background_list_2_gene_id",
            {"session_id": session_id, "ids": background_ids},
        )

    def term_membership_count(
        self,
        analysis_species: int,
        data: list[str],
    ) -> list[dict[str, Any]]:
        result = self.post_json(
            "termMembershipCount",
            {"analysisSpecies": int(analysis_species), "data": data},
        )
        return _json_loads(result["result"])

    def run_enrichment(
        self,
        session_id: str,
        data: list[str],
        checked_nodes: list[Any],
        backgroundlist: list[str] | None,
        *,
        multiple_list: bool = True,
        p_cutoff: float = 0.01,
        min_overlap: str = "3",
        min_enrichment: float = 1.5,
        cluster_number_cutoff: int = 20,
        disable_ppi: bool = False,
        first_call: int = 1,
    ) -> Any:
        options: dict[str, Any] = {
            "wherePutResult": "CURRENT_SHEET",
            "hasHeader": False,
            "oldHeader": "Gene_ID",
            "one2many": "KEEP_FIRST_ONLY",
            "minOverlap": str(min_overlap),
            "pCutoff": float(p_cutoff),
            "minEnrichment": float(min_enrichment),
            "gpec": False,
            "go_selective": False,
            "qc_ignore": False,
            "ppi_database_type": "PHYSICAL_CORE",
            "clusterNumberCutoff": int(cluster_number_cutoff),
            "checkedNodes": checked_nodes,
            "session_id": session_id,
            "isMultipleList": bool(multiple_list),
            "isExpressAnalysis": False,
            "isL1000": False,
            "ppiOption": {
                "disablePPI": bool(disable_ppi),
                "minSize": "3",
                "maxSize": "500",
            },
        }
        if backgroundlist is not None:
            options["backgroundlist"] = backgroundlist
        return self.post_json(
            "enrichmentanalysismultiplelist",
            {
                "options": options,
                "data": data,
                "async": True,
                "first_call": int(first_call),
            },
        )

    def job_status(self, session_id: str, custom: str) -> dict[str, Any]:
        return self.get_json(
            "get_job_status",
            {"custom": custom, "session_id": session_id},
        )

    def wait_for_job(
        self,
        session_id: str,
        custom: str,
        timeout_seconds: int,
        poll_interval: int,
    ) -> dict[str, Any]:
        deadline = time.time() + timeout_seconds
        last_status: dict[str, Any] = {}
        while time.time() < deadline:
            last_status = self.job_status(session_id, custom)
            complete = last_status.get("complete")
            total = last_status.get("total")
            status = last_status.get("last_act") or last_status.get("next_action")
            if total and complete == total:
                return last_status
            log(
                f"Metascape {custom} status: "
                f"{complete}/{total} {status or ''}".strip()
            )
            time.sleep(poll_interval)
        raise MetascapeError(
            f"Timed out waiting for Metascape {custom} job. "
            f"Last status: {last_status}"
        )

    def make_report(self, session_id: str) -> str:
        return self.get_text("make_analysis_report", {"session_id": session_id})

    def report_information(self, session_id: str) -> dict[str, Any]:
        return self.get_json("get_report_information", {"session_id": session_id})

    def wait_for_report(
        self,
        session_id: str,
        timeout_seconds: int,
        poll_interval: int,
    ) -> dict[str, Any]:
        deadline = time.time() + timeout_seconds
        last: dict[str, Any] = {}
        while time.time() < deadline:
            last = self.report_information(session_id)
            files = set(last.get("all_files", []))
            if "all.zip" in files or "metascape_result.xlsx" in files:
                return last
            log("Metascape report is still being generated")
            time.sleep(poll_interval)
        raise MetascapeError(
            f"Timed out waiting for Metascape report. Last response: {last}"
        )

    def download_file(self, session_id: str, file_name: str, output_path: str | Path) -> None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        url = self._url(
            "get_file",
            {"session_id": session_id, "file_name": file_name, "rename": "true"},
        )
        req = request.Request(url, method="GET")
        output_path.write_bytes(self._open(req))


def converted_gene_ids(applied_species_response: dict[str, Any]) -> list[str]:
    rows = _json_loads(applied_species_response["content"])
    return _dedupe([row.get("Gene", "") for row in rows])


def final_background(
    converted_background: dict[str, Any],
    input_genes: list[str],
) -> list[str]:
    server_background = _dedupe(converted_background.get("gene_ids", []))
    return _dedupe(server_background + [g for g in input_genes if g not in server_background])


def compact_report_files(report_info: dict[str, Any], include_zip: bool) -> list[str]:
    available = set(report_info.get("all_files", []))
    files = [
        name
        for name in DEFAULT_REPORT_FILES
        if name.split("/")[0] in available or name in available
    ]
    if include_zip and "all.zip" in available:
        files.append("all.zip")
    return files


def download_report_files(
    client: MetascapeClient,
    session_id: str,
    files: list[str],
    output_dir: Path,
) -> list[dict[str, str]]:
    warnings: list[dict[str, str]] = []
    for file_name in files:
        destination = output_dir / file_name
        log(f"Downloading {file_name}")
        try:
            client.download_file(session_id, file_name, destination)
        except MetascapeError as exc:
            message = str(exc)
            if "HTTP 404" not in message:
                raise
            log(f"Skipping missing optional Metascape file: {file_name}")
            warnings.append({"file": file_name, "reason": "HTTP 404"})
    if warnings:
        pd_rows = ["file\treason"] + [
            f"{row['file']}\t{row['reason']}" for row in warnings
        ]
        (output_dir / "download_warnings.tsv").write_text(
            "\n".join(pd_rows) + "\n",
            encoding="utf-8",
        )
    return warnings


def submit(args: argparse.Namespace) -> None:
    output_dir_template = str(args.output_dir) if args.output_dir else None
    client = MetascapeClient(args.base_url, timeout=args.request_timeout)

    session_id = client.create_session()
    report_url = f"https://metascape.org/gp/index.html#/reportfinal/{session_id}"
    if output_dir_template:
        output_dir = Path(output_dir_template.format(session_id=session_id))
    else:
        output_dir = Path("outputs/metascape_runs") / session_id
    output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Created Metascape session {session_id}")

    upload_response = client.upload_gene_lists(session_id, args.gene_lists)
    (output_dir / "upload_response_summary.json").write_text(
        json.dumps(
            {
                "session_id": session_id,
                "found_header": upload_response.get("found_header"),
                "isMultipleList": upload_response.get("isMultipleList"),
                "guess_type": upload_response.get("guess_type"),
                "filename": upload_response.get("filename"),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    log("Uploaded gene-list CSV")

    applied = client.apply_species(
        session_id,
        input_species=args.input_species,
        analysis_species=args.analysis_species,
    )
    data = converted_gene_ids(applied)
    if not data:
        raise MetascapeError("Species application returned no converted gene IDs")
    log(f"Applied species and retained {len(data):,} converted input genes")

    background = None
    background_summary: dict[str, Any] | None = None
    if args.background:
        background_ids = read_ids(args.background)
        converted = client.convert_background(session_id, background_ids)
        background = final_background(converted, data)
        background_summary = {
            "submitted_background_ids": len(background_ids),
            "recognized_background_ids": len(converted.get("gene_ids", [])),
            "unrecognized_background_ids": len(converted.get("wrong_input_ids", [])),
            "final_background_ids": len(background),
            "first_unrecognized": converted.get("wrong_input_ids", [])[:20],
        }
        (output_dir / "background_summary.json").write_text(
            json.dumps(background_summary, indent=2) + "\n",
            encoding="utf-8",
        )
        log(
            "Prepared custom background: "
            f"{background_summary['final_background_ids']:,} final IDs"
        )

    tree = client.term_membership_count(args.analysis_species, data)
    checked_nodes = selected_term_ids(tree)
    if not checked_nodes:
        raise MetascapeError("Metascape returned no checked ontology terms")
    log(f"Selected {len(checked_nodes):,} ontology/source terms")

    run_summary = {
        "session_id": session_id,
        "report_url": report_url,
        "gene_lists": str(args.gene_lists),
        "background": str(args.background) if args.background else None,
        "analysis_species": args.analysis_species,
        "input_species": args.input_species,
        "converted_input_genes": len(data),
        "checked_nodes": len(checked_nodes),
        "background_summary": background_summary,
        "p_cutoff": args.p_cutoff,
        "min_overlap": args.min_overlap,
        "min_enrichment": args.min_enrichment,
        "disable_ppi": args.disable_ppi,
    }

    if args.prepare_only:
        (output_dir / "metascape_run_summary.json").write_text(
            json.dumps(run_summary, indent=2) + "\n",
            encoding="utf-8",
        )
        log(f"Prepared session only. Report URL will be {report_url}")
        return

    client.run_enrichment(
        session_id,
        data,
        checked_nodes,
        background,
        p_cutoff=args.p_cutoff,
        min_overlap=args.min_overlap,
        min_enrichment=args.min_enrichment,
        cluster_number_cutoff=args.cluster_number_cutoff,
        disable_ppi=args.disable_ppi,
    )
    timeout_seconds = int(args.timeout_minutes * 60)
    client.wait_for_job(session_id, "GO", timeout_seconds, args.poll_interval)
    if not args.disable_ppi:
        client.wait_for_job(session_id, "PPI", timeout_seconds, args.poll_interval)

    log("Requesting Metascape report generation")
    client.make_report(session_id)
    report_info = client.wait_for_report(
        session_id, timeout_seconds, args.poll_interval
    )
    (output_dir / "AnalysisReport.html").write_text(
        report_info.get("analysisReportHtml", ""),
        encoding="utf-8",
    )
    (output_dir / "report_information.json").write_text(
        json.dumps(
            {
                key: value
                for key, value in report_info.items()
                if key != "analysisReportHtml"
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    if not args.no_download:
        download_report_files(
            client,
            session_id,
            compact_report_files(report_info, args.include_zip),
            output_dir,
        )

    (output_dir / "metascape_run_summary.json").write_text(
        json.dumps(run_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    log(f"Metascape report: {report_url}")
    log(f"Outputs: {output_dir}")


def status(args: argparse.Namespace) -> None:
    client = MetascapeClient(args.base_url, timeout=args.request_timeout)
    report_info = client.report_information(args.session_id)
    if not args.full:
        report_info = {
            key: value
            for key, value in report_info.items()
            if key != "analysisReportHtml"
        }
    print(json.dumps(report_info, indent=2))


def download(args: argparse.Namespace) -> None:
    client = MetascapeClient(args.base_url, timeout=args.request_timeout)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_info = client.report_information(args.session_id)
    (output_dir / "AnalysisReport.html").write_text(
        report_info.get("analysisReportHtml", ""),
        encoding="utf-8",
    )
    files = args.files or compact_report_files(report_info, args.include_zip)
    download_report_files(client, args.session_id, files, output_dir)
    (output_dir / "report_information.json").write_text(
        json.dumps(
            {
                key: value
                for key, value in report_info.items()
                if key != "analysisReportHtml"
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Automate the Metascape web workflow for GLARE gene lists."
    )
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--request-timeout", type=int, default=180)
    subparsers = parser.add_subparsers(dest="command", required=True)

    submit_parser = subparsers.add_parser(
        "submit", help="Upload gene lists, run enrichment, and download results."
    )
    submit_parser.add_argument("--gene-lists", required=True, type=Path)
    submit_parser.add_argument("--background", type=Path)
    submit_parser.add_argument("--output-dir", type=Path)
    submit_parser.add_argument("--input-species", type=int, default=DEFAULT_INPUT_SPECIES)
    submit_parser.add_argument(
        "--analysis-species", type=int, default=DEFAULT_ANALYSIS_SPECIES
    )
    submit_parser.add_argument("--p-cutoff", type=float, default=0.01)
    submit_parser.add_argument("--min-overlap", default="3")
    submit_parser.add_argument("--min-enrichment", type=float, default=1.5)
    submit_parser.add_argument("--cluster-number-cutoff", type=int, default=20)
    submit_parser.add_argument("--disable-ppi", action="store_true")
    submit_parser.add_argument("--prepare-only", action="store_true")
    submit_parser.add_argument("--no-download", action="store_true")
    submit_parser.add_argument("--include-zip", action="store_true")
    submit_parser.add_argument("--timeout-minutes", type=float, default=90)
    submit_parser.add_argument("--poll-interval", type=int, default=30)
    submit_parser.set_defaults(func=submit)

    status_parser = subparsers.add_parser("status", help="Print report JSON.")
    status_parser.add_argument("--session-id", required=True)
    status_parser.add_argument("--full", action="store_true")
    status_parser.set_defaults(func=status)

    download_parser = subparsers.add_parser(
        "download", help="Download compact files for an existing session."
    )
    download_parser.add_argument("--session-id", required=True)
    download_parser.add_argument("--output-dir", required=True, type=Path)
    download_parser.add_argument("--files", nargs="*")
    download_parser.add_argument("--include-zip", action="store_true")
    download_parser.set_defaults(func=download)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
