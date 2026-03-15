#!/usr/bin/env python3
"""Publish a dataset to Zenodo via the REST API.

Creates a new deposit, uploads files, sets metadata, and publishes.
Metadata is read from a YAML config file specific to each dataset.

Prerequisites:
    pip install requests pyyaml

    Generate an API token at:
    https://zenodo.org/account/settings/applications/
    (Personal access tokens -> New token -> deposit:actions + deposit:write)

    Store the token in ~/.zenodo_token or pass via --token.

Usage:
    python pipeline/publish_to_zenodo.py DATASET_DIR [--token TOKEN] [--dry-run]

Example:
    python pipeline/publish_to_zenodo.py top2018_full/
    python pipeline/publish_to_zenodo.py top2018_mc/ --dry-run

The DATASET_DIR must contain:
    zenodo.yaml         - metadata config
    zenodo_upload/      - directory with files to upload
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import requests
    import yaml
except ImportError:
    print(
        "ERROR: requires 'requests' and 'pyyaml'. "
        "Install with: pip install requests pyyaml",
        file=sys.stderr,
    )
    sys.exit(1)

ZENODO_API = "https://zenodo.org/api"


def load_token(token_arg: str | None) -> str:
    """Load Zenodo API token from argument or ~/.zenodo_token."""
    if token_arg:
        return token_arg
    token_file = Path.home() / ".zenodo_token"
    if token_file.exists():
        return token_file.read_text().strip()
    print(
        "ERROR: No API token. Pass --token or create ~/.zenodo_token",
        file=sys.stderr,
    )
    sys.exit(1)


def create_deposit(token: str) -> tuple[int, str]:
    """Create a new empty deposit. Returns (deposit_id, bucket_url)."""
    r = requests.post(
        f"{ZENODO_API}/deposit/depositions",
        params={"access_token": token},
        json={},
        headers={"Content-Type": "application/json"},
    )
    r.raise_for_status()
    data = r.json()
    return data["id"], data["links"]["bucket"]


def upload_file(bucket_url: str, filepath: Path, token: str) -> None:
    """Upload a single file to the deposit bucket."""
    print(f"  Uploading {filepath.name} ({filepath.stat().st_size / 1e6:.1f} MB)...")
    with open(filepath, "rb") as fh:
        r = requests.put(
            f"{bucket_url}/{filepath.name}",
            data=fh,
            params={"access_token": token},
        )
    r.raise_for_status()
    print(f"    done (md5:{r.json()['checksum'].replace('md5:', '')})")


def set_metadata(deposit_id: int, metadata: dict, token: str) -> None:
    """Set deposit metadata."""
    r = requests.put(
        f"{ZENODO_API}/deposit/depositions/{deposit_id}",
        params={"access_token": token},
        data=json.dumps({"metadata": metadata}),
        headers={"Content-Type": "application/json"},
    )
    r.raise_for_status()


def publish(deposit_id: int, token: str) -> str:
    """Publish the deposit. Returns the DOI."""
    r = requests.post(
        f"{ZENODO_API}/deposit/depositions/{deposit_id}/actions/publish",
        params={"access_token": token},
    )
    r.raise_for_status()
    return r.json()["doi"]


def build_metadata(config: dict) -> dict:
    """Build Zenodo metadata dict from YAML config."""
    # Required fields
    meta: dict = {
        "title": config["title"],
        "upload_type": "dataset",
        "description": config["description"],
        "creators": config["creators"],
        "access_right": "open",
        "license": config.get("license", "cc-by-4.0"),
        "version": config.get("version", "1.0.0"),
    }

    # Optional fields
    if "keywords" in config:
        meta["keywords"] = config["keywords"]
    if "notes" in config:
        meta["notes"] = config["notes"]
    if "publication_date" in config:
        meta["publication_date"] = config["publication_date"]
    if "related_identifiers" in config:
        meta["related_identifiers"] = config["related_identifiers"]
    if "references" in config:
        meta["references"] = config["references"]
    if "contributors" in config:
        meta["contributors"] = config["contributors"]
    if "subjects" in config:
        meta["subjects"] = config["subjects"]

    return meta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Publish a dataset to Zenodo via REST API"
    )
    parser.add_argument(
        "dataset_dir", type=Path,
        help="Dataset directory containing zenodo.yaml and zenodo_upload/",
    )
    parser.add_argument(
        "--token", type=str, default=None,
        help="Zenodo API token (or store in ~/.zenodo_token)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate config and list files without uploading",
    )
    args = parser.parse_args()

    # Load config
    config_path = args.dataset_dir / "zenodo.yaml"
    if not config_path.exists():
        print(f"ERROR: {config_path} not found", file=sys.stderr)
        sys.exit(1)

    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    # Find upload files
    upload_dir = args.dataset_dir / "zenodo_upload"
    if not upload_dir.is_dir():
        print(f"ERROR: {upload_dir} not found", file=sys.stderr)
        sys.exit(1)

    files = sorted(upload_dir.iterdir())
    if not files:
        print(f"ERROR: no files in {upload_dir}", file=sys.stderr)
        sys.exit(1)

    # Build metadata
    metadata = build_metadata(config)

    print(f"Dataset: {config['title']}")
    print(f"Version: {metadata.get('version', '?')}")
    print(f"Files ({len(files)}):")
    for f in files:
        print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")

    if args.dry_run:
        print("\nMetadata:")
        print(json.dumps(metadata, indent=2))
        print("\n[DRY RUN] Would create deposit and upload above files.")
        return

    # Load token and publish
    token = load_token(args.token)

    print("\nCreating deposit...")
    deposit_id, bucket_url = create_deposit(token)
    print(f"  Deposit ID: {deposit_id}")

    print("\nUploading files...")
    for f in files:
        upload_file(bucket_url, f, token)

    print("\nSetting metadata...")
    set_metadata(deposit_id, metadata, token)

    print("\nPublishing...")
    doi = publish(deposit_id, token)
    print(f"\nPublished! DOI: {doi}")
    print(f"URL: https://doi.org/{doi}")


if __name__ == "__main__":
    main()
