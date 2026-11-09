#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.13"
# dependencies = []
# ///

import tomllib
from pathlib import Path


def main() -> None:
    file_path = Path(__file__).parent.parent / "pyproject.toml"
    with open(file_path, "rb") as f:
        data = tomllib.load(f)
    authors = []
    for author in data["project"]["authors"]:
        authors.append(author["email"])
    print(", ".join(authors))


if __name__ == "__main__":
    main()
