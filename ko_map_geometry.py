"""Read KO KGML graphics without fetching XML resources or inventing coordinates.

Each ortholog entry is a feature group, potentially containing several KOs.  Its
graphics remain a single biological feature group in the renderer and reports.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


MAX_KGML_BYTES = 32 * 1024 * 1024
KO_PATTERN = re.compile(r"ko:(K\d{5})", re.IGNORECASE)
MAP_PATTERN = re.compile(r"(?:ko|map)(\d{5})", re.IGNORECASE)
# Official KGML files carry this declaration.  Remove it before XML parsing;
# neither this DTD nor any other remote resource is opened.
OFFICIAL_DOCTYPE = re.compile(
    r'<!DOCTYPE\s+pathway\s+SYSTEM\s+([\'"])'
    r'https?://www\.(?:kegg|genome)\.jp/kegg/xml/KGML_v[0-9]+(?:\.[0-9]+)+_\.dtd\1\s*>',
    re.IGNORECASE,
)


def normalize_map_id(value: str) -> str:
    """Canonicalize a display mapNNNNN or KO reference koNNNNN to koNNNNN."""
    if not isinstance(value, str) or not (match := MAP_PATTERN.fullmatch(value.strip())):
        raise ValueError("KEGG map ID must be koNNNNN or mapNNNNN")
    return "ko" + match.group(1)


def _number(value: str | None, field: str, *, positive: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid KGML graphic {field}") from None
    if not math.isfinite(number) or abs(number) > 1_000_000 or positive and number <= 0:
        raise ValueError(f"Invalid KGML graphic {field}")
    return number


def parse_kgml_text(text: str, map_id: str | None = None) -> dict:
    """Parse a local UTF-8 KGML document into JSON-safe original geometry.

    Only the canonical KEGG external DOCTYPE is stripped, without resolving it.
    All entity declarations, internal subsets, alternate DTDs, namespaces and
    XML include mechanisms are rejected.  Links in KGML remain inert metadata.
    """
    if not isinstance(text, str) or len(text.encode("utf-8")) > MAX_KGML_BYTES:
        raise ValueError("KGML exceeds the 32 MiB local map size limit")
    safe_text = OFFICIAL_DOCTYPE.sub("", text)
    if re.search(r"<!\s*(?:DOCTYPE|ENTITY)\b|\bxmlns\s*=|\bxmlns:", safe_text, re.IGNORECASE):
        raise ValueError("KGML DTDs, entities and external XML references are not supported")
    try:
        root = ET.fromstring(safe_text)
    except ET.ParseError as error:
        raise ValueError(f"Invalid KGML XML: {error}") from None
    if root.tag != "pathway" or any("}" in element.tag or ":" in element.tag for element in root.iter()):
        raise ValueError("Expected a plain KGML pathway document without XML includes")
    source_name = root.get("name", "").removeprefix("path:")
    if not re.fullmatch(r"ko\d{5}", source_name):
        raise ValueError("Use the KO reference KGML (path:koNNNNN), not organism-specific or map KGML")
    canonical = normalize_map_id(source_name)
    if root.get("org", "ko") != "ko" or root.get("number", canonical[2:]) != canonical[2:]:
        raise ValueError("KGML pathway name, organism and number disagree")
    if map_id is not None and normalize_map_id(map_id) != canonical:
        raise ValueError("KGML pathway does not match the selected map ID")

    entries, all_bounds, seen_ids, warnings = [], [], set(), []
    for element in root.findall("entry"):
        entry_id = element.get("id", "")
        if not re.fullmatch(r"\d+", entry_id) or entry_id in seen_ids:
            raise ValueError("KGML entry IDs must be unique nonnegative integers")
        seen_ids.add(entry_id)
        names = element.get("name", "").split()
        ko_ids = []
        if element.get("type") == "ortholog":
            for name in names:
                match = KO_PATTERN.fullmatch(name)
                if match and (ko := match.group(1).upper()) not in ko_ids:
                    ko_ids.append(ko)
        graphics = []
        for index, graphic in enumerate(element.findall("graphics"), 1):
            shape = graphic.get("type", "rectangle")
            if shape not in ("line", "rectangle", "roundrectangle", "circle"):
                warnings.append(f"Entry {entry_id}: unsupported graphics type {shape!r} omitted")
                continue
            geometry = {"id": f"{entry_id}:{index}", "type": shape,
                        "label": graphic.get("name", ""), "coords": []}
            if shape == "line":
                pieces = graphic.get("coords", "").split(",")
                if len(pieces) < 4 or len(pieces) % 2:
                    raise ValueError(f"KGML line {geometry['id']} needs at least two coordinate pairs")
                numbers = [_number(piece, "coords") for piece in pieces]
                geometry["coords"] = [numbers[pos:pos + 2] for pos in range(0, len(numbers), 2)]
                xs, ys = numbers[::2], numbers[1::2]
                geometry.update(x=(min(xs) + max(xs)) / 2, y=(min(ys) + max(ys)) / 2,
                                width=max(xs) - min(xs), height=max(ys) - min(ys))
                all_bounds.append((min(xs), max(xs), min(ys), max(ys)))
            else:
                geometry.update(x=_number(graphic.get("x"), "x"), y=_number(graphic.get("y"), "y"),
                                width=_number(graphic.get("width"), "width", positive=True),
                                height=_number(graphic.get("height"), "height", positive=True))
                x, y, width, height = (geometry[key] for key in ("x", "y", "width", "height"))
                all_bounds.append((x - width / 2, x + width / 2, y - height / 2, y + height / 2))
            graphics.append(geometry)
        entries.append({"id": entry_id, "type": element.get("type", ""), "names": names,
                        "ko_ids": ko_ids, "graphics": graphics,
                        "components": [child.get("id", "") for child in element.findall("component")],
                        "reaction": element.get("reaction", "")})
    ko_entries = [entry for entry in entries if entry["ko_ids"] and entry["graphics"]]
    if not ko_entries:
        raise ValueError("KGML contains no valid KO ortholog graphics to position differential results")
    # Relationships are preserved as metadata.  Only explicit graphics supply
    # drawable paths: no centre-to-centre or force-layout edges are generated.
    relations = [{**item.attrib, "subtypes": [dict(child.attrib) for child in item.findall("subtype")]}
                 for item in root.findall("relation")]
    reactions = [{**item.attrib,
                  "substrates": [dict(child.attrib) for child in item.findall("substrate")],
                  "products": [dict(child.attrib) for child in item.findall("product")]}
                 for item in root.findall("reaction")]
    return {"schema_version": 1, "source": "ko_network", "map_id": canonical,
            "title": root.get("title", canonical), "coordinate_system": "KGML_original_y_down",
            "source_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            "entries": entries, "relations": relations, "reactions": reactions,
            "bounds": {"xmin": min(row[0] for row in all_bounds), "xmax": max(row[1] for row in all_bounds),
                       "ymin": min(row[2] for row in all_bounds), "ymax": max(row[3] for row in all_bounds)},
            "counts": {"entries": len(entries), "ko_feature_groups": len(ko_entries),
                       "ko_graphics": sum(len(entry["graphics"]) for entry in ko_entries),
                       "unique_kos": len({ko for entry in ko_entries for ko in entry["ko_ids"]})},
            "warnings": warnings}


def parse_kgml(path: str | Path, map_id: str | None = None) -> dict:
    path = Path(path)
    if path.stat().st_size > MAX_KGML_BYTES:
        raise ValueError("KGML exceeds the 32 MiB local map size limit")
    try:
        return parse_kgml_text(path.read_text(encoding="utf-8-sig"), map_id)
    except UnicodeError:
        raise ValueError("KGML must be UTF-8 XML") from None


def write_geometry(kgml_path: str | Path, output_path: str | Path, map_id: str | None = None) -> dict:
    geometry = parse_kgml(kgml_path, map_id)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(geometry, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return geometry


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kgml_path")
    parser.add_argument("output_path")
    parser.add_argument("--map-id")
    args = parser.parse_args()
    write_geometry(args.kgml_path, args.output_path, args.map_id)
