#!/usr/bin/env python3
"""
Extract CHARMM-GUI membrane-system metadata into a text report.

Expected files are auto-detected recursively:
    input.config.dat
    step3_nlipids_upper.prm
    step3_nlipids_lower.prm
    step4_components.str
    step4.3_ion.str
    step2_orient.out
    gromacs/step5_input.gro

Python requirement: Python 3.9+ (Standard Library only)
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

NUM_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
NUM_RE = re.compile(NUM_PATTERN)

DEFAULT_WATER_NAMES = {
    "SOL", "WAT", "HOH", "TIP3", "TIP3P", "TIP4", "TIP4P",
    "TIP5", "TIP5P", "SPC", "SPCE", "SWM4", "OH2"
}

DEFAULT_ION_NAMES = {
    "SOD", "NA", "NA+", "POT", "K", "K+",
    "CLA", "CL", "CL-", "CAL", "CA", "CA2", "CA2+",
    "MG", "MG2", "MG2+", "ZN", "ZN2", "ZN2+",
    "CES", "CS", "LIT", "LI", "RB",
    "BAR", "BA", "BA2", "CD", "CD2", "CU", "CU1", "CU2",
    "FE", "FE2", "FE3", "MN", "MN2", "AL", "AL3",
    "F", "F-", "BR", "BR-", "IOD", "I", "I-"
}

COMMON_LIPID_NAMES = {
    "POPC", "POPE", "POPG", "POPS", "POPA",
    "DOPC", "DOPE", "DOPG", "DOPS", "DOPA",
    "DPPC", "DPPE", "DPPG", "DPPS", "DPPA",
    "DMPC", "DMPE", "DMPG", "DMPS", "DMPA",
    "DSPC", "DSPE", "DSPG", "DSPS",
    "DLPC", "DLPE", "DLPG", "DLPS",
    "SOPC", "SOPE", "SOPG", "SOPS",
    "PAPC", "PAPE", "PAPG", "PAPS",
    "SAPC", "SDPC", "SLPC", "DAPC", "DUPC",
    "DEPC", "DHPC", "DIPC",
    "PSM", "SSM", "NSM", "SM",
    "CHOL", "CHL1", "CHL",
    "TOCL1", "TOCL2", "TLCL1", "TLCL2", "CARD",
    "CER", "CER1", "CER2", "CER3", "CER160", "CER180",
    "GLCER", "BGLC", "GM1", "GM3",
    "LPA", "LPC", "LPE", "LPG", "LPS",
    "PIP", "PIP2", "PIP3", "PI", "POPI",
    "DAG", "TAG"
}

TARGET_FILES = {
    "config": "input.config.dat",
    "lipid_upper": "step3_nlipids_upper.prm",
    "lipid_lower": "step3_nlipids_lower.prm",
    "components": "step4_components.str",
    "ions": "step4.3_ion.str",
    "orientation": "step2_orient.out",
    "gro": "step5_input.gro",
}


@dataclass(frozen=True)
class Evidence:
    path: Path
    line_number: int
    text: str


@dataclass
class Assignment:
    key: str
    value: str
    line_number: int
    raw: str


@dataclass
class LipidInfo:
    counts: dict[str, int] = field(default_factory=dict)
    total: Optional[int] = None
    evidence: list[Evidence] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class IonInfo:
    types: dict[str, str] = field(default_factory=dict)
    count_evidence: list[Evidence] = field(default_factory=list)
    concentration_evidence: list[Evidence] = field(default_factory=list)
    type_evidence: list[Evidence] = field(default_factory=list)


@dataclass
class GroInfo:
    title: str
    declared_atom_count: int
    parsed_atom_count: int
    residue_counts: Counter
    residue_atom_counts: Counter
    raw_box_line: str
    raw_box_values: list[str]
    vectors: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    vector_tokens: tuple[
        tuple[str, str, str],
        tuple[str, str, str],
        tuple[str, str, str],
    ]
    inferred_box_type: str
    vector_lengths: tuple[float, float, float]
    gamma_degrees: float
    xy_area_nm2: float
    volume_nm3: float
    warnings: list[str] = field(default_factory=list)


def natural_key(value: str):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    ]


def normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def strip_inline_comment(line: str) -> str:
    result = line
    for marker in ("!", "#"):
        pos = result.find(marker)
        if pos >= 0:
            result = result[:pos]
    return result.strip()


def parse_assignments(path: Path) -> list[Assignment]:
    assignments = []
    for line_number, raw in enumerate(read_text(path).splitlines(), 1):
        stripped = strip_inline_comment(raw)
        if not stripped or stripped.startswith("*"):
            continue

        match = re.match(
            r"^\s*set\s+([A-Za-z_][A-Za-z0-9_.()\[\]-]*)\s*(?:=\s*)?(.+?)\s*$",
            stripped,
            re.IGNORECASE,
        )
        if not match:
            match = re.match(
                r"^\s*([A-Za-z_][A-Za-z0-9_.()\[\]-]*)\s*(?:=|:)\s*(.+?)\s*$",
                stripped,
            )
        if not match:
            match = re.match(
                r"^\s*([A-Za-z_][A-Za-z0-9_.()\[\]+\--]*)\s+([^\s].*?)\s*$",
                stripped,
            )

        if match:
            assignments.append(
                Assignment(
                    key=match.group(1).strip(),
                    value=match.group(2).strip().strip("\"'"),
                    line_number=line_number,
                    raw=raw.rstrip(),
                )
            )
    return assignments


def numeric_values(value: str) -> list[float]:
    values = []
    for item in NUM_RE.findall(value):
        try:
            values.append(float(item.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return values


def integer_values(value: str) -> list[int]:
    result = []
    for number in numeric_values(value):
        if math.isfinite(number) and abs(number - round(number)) < 1.0e-8:
            result.append(int(round(number)))
    return result


def first_integer(value: str) -> Optional[int]:
    values = integer_values(value)
    return values[0] if values else None


def unique_evidence(records: Iterable[Evidence]) -> list[Evidence]:
    result = []
    seen = set()
    for record in records:
        key = (record.path.resolve(), record.line_number, record.text)
        if key not in seen:
            seen.add(key)
            result.append(record)
    return result


def detect_files(root: Path) -> tuple[dict[str, Optional[Path]], list[str]]:
    warnings = []
    all_files = [path for path in root.rglob("*") if path.is_file()]
    by_name = defaultdict(list)

    for path in all_files:
        by_name[path.name.lower()].append(path)

    selected: dict[str, Optional[Path]] = {}

    for label, expected_name in TARGET_FILES.items():
        candidates = by_name.get(expected_name.lower(), [])
        if not candidates:
            selected[label] = None
            continue

        if label == "gro":
            candidates.sort(
                key=lambda p: (
                    0 if "gromacs" in [part.lower() for part in p.parts] else 1,
                    len(p.parts),
                    natural_key(str(p)),
                )
            )
        else:
            candidates.sort(key=lambda p: (len(p.parts), natural_key(str(p))))

        selected[label] = candidates[0]
        if len(candidates) > 1:
            warnings.append(
                f"Multiple {expected_name} files found; selected {relative_path(candidates[0], root)}."
            )

    return selected, warnings


def is_species_token(value: str) -> bool:
    token = value.strip().strip("\"'").upper()
    if not re.fullmatch(r"[A-Z][A-Z0-9+\-]{1,15}", token):
        return False

    excluded = {
        "YES", "NO", "TRUE", "FALSE", "NONE", "NULL",
        "UPPER", "LOWER", "TOTAL", "TYPE", "NAME",
        "LIPID", "LIPIDS", "NLIPID", "NLIPIDS",
        "WATER", "ION", "IONS", "PROTEIN", "PEPTIDE",
        "SIZE", "XSIZE", "YSIZE", "BOX", "XBOX", "YBOX",
        "AREA", "RADIUS", "SEED", "COUNT", "NUMBER",
        "RECTANGULAR", "HEXAGONAL", "TRICLINIC",
    }
    if token in excluded or token.startswith(("NLIPID", "NUM", "COUNT")):
        return False
    return True


def words_from_value(value: str) -> list[str]:
    return re.findall(r"[A-Za-z][A-Za-z0-9+\-]{1,15}", value)


def parse_lipid_file(path: Optional[Path]) -> LipidInfo:
    info = LipidInfo()
    if path is None:
        return info

    assignments = parse_assignments(path)
    indexed_species: dict[str, tuple[str, Assignment]] = {}
    indexed_counts: dict[str, tuple[int, Assignment]] = {}
    species_lists: list[tuple[list[str], Assignment]] = []
    count_lists: list[tuple[list[int], Assignment]] = []
    direct_counts: list[tuple[str, int, Assignment]] = []
    total_declarations: list[tuple[int, Assignment]] = []

    for assignment in assignments:
        key = normalized_key(assignment.key)
        words = [
            word.upper()
            for word in words_from_value(assignment.value)
            if is_species_token(word)
        ]
        integers = integer_values(assignment.value)

        count_match = (
            re.fullmatch(r"(?:n|num|count)lipid(\d+)", key)
            or re.fullmatch(r"lipid(\d+)(?:n|num|count)", key)
        )
        species_match = (
            re.fullmatch(r"lipid(\d+)", key)
            or re.fullmatch(r"lipid(?:type|name|resname)(\d+)", key)
            or re.fullmatch(r"(?:type|name|resname)lipid(\d+)", key)
            or re.fullmatch(r"lipid(\d+)(?:type|name|resname)", key)
        )

        if count_match and integers:
            indexed_counts[count_match.group(1)] = (integers[0], assignment)
            continue

        if species_match and words:
            indexed_species[species_match.group(1)] = (words[0], assignment)
            continue

        if key in {
            "lipid", "lipids", "lipidtype", "lipidtypes",
            "lipidname", "lipidnames", "resname", "resnames"
        } and words:
            species_lists.append((words, assignment))

        if key in {
            "nlipid", "nlipids", "numlipid", "numlipids",
            "lipidcount", "lipidcounts"
        } and integers:
            count_lists.append((integers, assignment))
            total_declarations.append((sum(integers), assignment))

        direct_match = re.fullmatch(
            r"(?:n|num|count)[_-]?([A-Za-z][A-Za-z0-9+\-]{1,15})",
            assignment.key,
            re.IGNORECASE,
        )
        if direct_match and integers:
            species = direct_match.group(1).upper()
            if is_species_token(species):
                direct_counts.append((species, integers[0], assignment))
                continue

        if (
            integers
            and assignment.key.upper() == assignment.key
            and is_species_token(assignment.key)
        ):
            direct_counts.append(
                (assignment.key.upper(), integers[0], assignment)
            )

    paired_records = []
    for index in sorted(set(indexed_species).intersection(indexed_counts), key=natural_key):
        species, s_assign = indexed_species[index]
        count, c_assign = indexed_counts[index]
        paired_records.append((species, count, s_assign, c_assign))

    for species_entry, count_entry in zip(species_lists, count_lists):
        s_vals, s_assign = species_entry
        c_vals, c_assign = count_entry
        if len(s_vals) == len(c_vals):
            for sp, cnt in zip(s_vals, c_vals):
                paired_records.append((sp, cnt, s_assign, c_assign))
        elif len(s_vals) == 1 and len(c_vals) == 1:
            paired_records.append((s_vals[0], c_vals[0], s_assign, c_assign))

    for species, count, s_assign, c_assign in paired_records:
        info.counts[species] = info.counts.get(species, 0) + count
        info.evidence.append(Evidence(path, s_assign.line_number, s_assign.raw))
        info.evidence.append(Evidence(path, c_assign.line_number, c_assign.raw))

    for species, count, assignment in direct_counts:
        if species not in info.counts:
            info.counts[species] = count
            info.evidence.append(Evidence(path, assignment.line_number, assignment.raw))
        elif info.counts[species] != count:
            info.warnings.append(
                f"{path.name}: conflicting counts detected for {species}: {info.counts[species]} and {count}."
            )

    if info.counts:
        info.total = sum(info.counts.values())
    elif total_declarations:
        info.total = total_declarations[0][0]
        assignment = total_declarations[0][1]
        info.evidence.append(Evidence(path, assignment.line_number, assignment.raw))

    if info.total is not None and total_declarations:
        declared = total_declarations[0][0]
        if info.counts and declared != info.total:
            info.warnings.append(
                f"{path.name}: declared total {declared} differs from parsed per-species total {info.total}."
            )

    info.evidence = unique_evidence(info.evidence)
    return info


def collect_config_evidence(path: Optional[Path]) -> dict[str, list[Evidence]]:
    categories = {
        "box_type": [],
        "box_xy": [],
        "water_thickness": [],
        "final_box": [],
    }
    if path is None:
        return categories

    assignments = {item.line_number: item for item in parse_assignments(path)}
    x_keys = {"boxx", "xbox", "xsize", "sizex", "xdim", "xdimension", "xlength", "lengthx", "boxlengthx"}
    y_keys = {"boxy", "ybox", "ysize", "sizey", "ydim", "ydimension", "ylength", "lengthy", "boxlengthy"}
    xy_keys = {"boxxy", "xybox", "xysize", "boxxysize", "lateralsize", "lateralbox", "boxlateral", "xy"}

    for line_number, raw in enumerate(read_text(path).splitlines(), 1):
        low = raw.lower()
        assignment = assignments.get(line_number)
        key = normalized_key(assignment.key) if assignment else ""
        evidence = Evidence(path, line_number, raw.rstrip())

        if (
            key in {"boxtype", "boxshape", "pbcshape", "pbctype", "geometry"}
            or re.search(r"\bbox\s*(?:type|shape)\b", low)
            or re.search(r"\b(?:type|shape)\s*(?:of\s+)?(?:the\s+)?box\b", low)
        ):
            categories["box_type"].append(evidence)

        if (
            key in x_keys or key in y_keys or key in xy_keys
            or re.search(r"\bbox.*\b(?:xy|x/y|lateral)\b", low)
            or re.search(r"\b(?:xy|x/y|lateral).*\bbox\b", low)
        ):
            if not re.search(r"\b(?:final|calculated|computed)\b", low):
                categories["box_xy"].append(evidence)

        if (
            ("water" in key and any(t in key for t in ("thick", "height", "layer", "padding", "distance")))
            or key in {"waterthickness", "watthick", "hydrationthickness"}
            or re.search(r"\bwater.*\b(?:thickness|height|padding|layer|distance)\b", low)
            or re.search(r"\b(?:thickness|height|padding|distance).*\bwater\b", low)
        ):
            categories["water_thickness"].append(evidence)

        if (
            re.search(r"\b(?:final|calculated|computed).*\bbox\b", low)
            or re.search(r"\bbox.*\b(?:final|calculated|computed)\b", low)
            or (any(t in key for t in ("final", "calc", "computed")) and "box" in key)
        ):
            categories["final_box"].append(evidence)

    for k in categories:
        categories[k] = unique_evidence(categories[k])
    return categories


def collect_calculated_box_evidence(paths: Iterable[Optional[Path]]) -> list[Evidence]:
    records = []
    for path in paths:
        if path is None:
            continue
        for line_number, raw in enumerate(read_text(path).splitlines(), 1):
            low = raw.lower()
            if (
                re.search(r"\b(?:final|calculated|computed).*\bbox\b", low)
                or re.search(r"\bbox.*\b(?:final|calculated|computed)\b", low)
                or (
                    "box" in low
                    and re.search(r"\b(?:x|y|xy|x/y|lateral)\b", low)
                    and re.search(r"\b(?:size|length|dimension)\b", low)
                )
            ):
                records.append(Evidence(path, line_number, raw.rstrip()))
    return unique_evidence(records)


def collect_component_declarations(path: Optional[Path]) -> list[Evidence]:
    records = []
    if path is None:
        return records

    component_tokens = {
        "wat", "water", "sol", "tip3", "lipid", "prot",
        "protein", "peptide", "ion", "sod", "cla", "pot", "atom", "residue", "molecule"
    }

    for assignment in parse_assignments(path):
        key = normalized_key(assignment.key)
        count = first_integer(assignment.value)
        if count is None:
            continue

        relevant = any(token in key for token in component_tokens)
        count_style = key.startswith("n") or "count" in key or "number" in key or "num" in key

        if relevant and count_style:
            records.append(Evidence(path, assignment.line_number, assignment.raw))

    return unique_evidence(records)


def parse_ion_information(paths: Iterable[Optional[Path]]) -> IonInfo:
    info = IonInfo()
    assignments_by_path = []

    for path in paths:
        if path is None:
            continue
        assignments = parse_assignments(path)
        assignments_by_path.append((path, assignments))

        for line_number, raw in enumerate(read_text(path).splitlines(), 1):
            low = raw.lower()
            if re.search(r"\b(?:concentration|conc\.?|molarity|molar|salinity)\b", low):
                info.concentration_evidence.append(Evidence(path, line_number, raw.rstrip()))

    for path, assignments in assignments_by_path:
        for assignment in assignments:
            key = normalized_key(assignment.key)
            words = [word.upper() for word in words_from_value(assignment.value) if is_species_token(word)]
            if not words:
                continue

            label = None
            if key in {"pos", "positive", "cation", "posion", "positiveion"}:
                label = "positive"
            elif key in {"neg", "negative", "anion", "negion", "negativeion"}:
                label = "negative"
            elif "ion" in key and any(t in key for t in ("type", "name", "resname")):
                label = assignment.key
            elif key in {"ion", "ions", "iontype"}:
                label = assignment.key

            if label:
                info.types[label] = words[0]
                info.type_evidence.append(Evidence(path, assignment.line_number, assignment.raw))

    discovered_types = set(info.types.values())

    for path, assignments in assignments_by_path:
        for assignment in assignments:
            key = normalized_key(assignment.key)
            count = first_integer(assignment.value)
            if count is None:
                continue

            is_count = False
            if key in {
                "npos", "npositive", "ncation", "numpositive",
                "nneg", "nnegative", "nanion", "numnegative",
                "nion", "nions", "numion", "numions", "ioncount"
            }:
                is_count = True

            direct_match = re.fullmatch(r"(?:n|num|count)([a-z][a-z0-9+\-]{1,15})", key)
            if direct_match:
                possible_species = direct_match.group(1).upper()
                if possible_species in DEFAULT_ION_NAMES or possible_species in discovered_types:
                    is_count = True

            if is_count:
                info.count_evidence.append(Evidence(path, assignment.line_number, assignment.raw))

    info.type_evidence = unique_evidence(info.type_evidence)
    info.count_evidence = unique_evidence(info.count_evidence)
    info.concentration_evidence = unique_evidence(info.concentration_evidence)
    return info


def collect_orientation(path: Optional[Path]) -> dict[str, list[Evidence]]:
    result = {
        "tilt": [],
        "insertion_depth": [],
        "other_ppm": [],
    }
    if path is None:
        return result

    for line_number, raw in enumerate(read_text(path).splitlines(), 1):
        low = raw.lower()
        record = Evidence(path, line_number, raw.rstrip())

        if "tilt" in low:
            result["tilt"].append(record)

        if (
            "insertion depth" in low
            or re.search(r"\bdepth\b", low)
            or re.search(r"\bz[- ]?(?:position|coordinate|center)\b", low)
            or re.search(r"\bdistance.*\bmembrane\b", low)
        ):
            result["insertion_depth"].append(record)

        if (
            "ppm" in low
            or "orientation" in low
            or "hydrophobic thickness" in low
            or "membrane position" in low
        ):
            result["other_ppm"].append(record)

    for k in result:
        result[k] = unique_evidence(result[k])
    return result


def collect_area_evidence(paths: Iterable[Optional[Path]]) -> tuple[list[Evidence], list[Evidence]]:
    strong, weak = [], []
    for path in paths:
        if path is None:
            continue
        for line_number, raw in enumerate(read_text(path).splitlines(), 1):
            low = raw.lower()
            normalized = normalized_key(raw)
            record = Evidence(path, line_number, raw.rstrip())

            strong_match = (
                re.search(r"cross[- ]?section(?:al)?", low)
                or ("area" in low and re.search(r"\b(?:protein|peptide|solute|excluded?|exclusion)\b", low))
                or ("area" in normalized and any(t in normalized for t in ("prot", "pept", "cross", "excl")))
                or (re.search(r"\b(?:exclude|exclusion|excluded)\b", low) and re.search(r"\blipid", low) and re.search(r"\b(?:area|radius|distance)\b", low))
            )
            if strong_match:
                strong.append(record)
            elif "area" in low:
                weak.append(record)

    return unique_evidence(strong), unique_evidence(weak)


def vector_norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(sum(c * c for c in vector))


def dot(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return sum(a * b for a, b in zip(first, second))


def cross(first: tuple[float, float, float], second: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def parse_gro(path: Path) -> GroInfo:
    lines = read_text(path).splitlines()
    if len(lines) < 3:
        raise ValueError("GRO file has fewer than three lines.")

    title = lines[0].rstrip()
    try:
        declared_atom_count = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError("Invalid GRO atom-count line.") from exc

    atom_start = 2
    atom_end = atom_start + declared_atom_count
    if len(lines) <= atom_end:
        raise ValueError(
            f"GRO declares {declared_atom_count} atoms, but file is truncated before box vectors."
        )

    atom_lines = lines[atom_start:atom_end]
    residue_counts = Counter()
    residue_atom_counts = Counter()
    last_residue = None
    warnings = []

    for index, line in enumerate(atom_lines, 1):
        if len(line) < 15:
            warnings.append(f"Atom line {index} is shorter than standard GRO format.")

        resid_text = line[0:5].strip()
        resname = line[5:10].strip().upper()

        try:
            resid = int(resid_text)
        except ValueError:
            fields = line.split()
            if not fields:
                continue
            match = re.match(r"(\d+)([A-Za-z0-9+\-]+)", fields[0])
            if not match:
                continue
            resid = int(match.group(1))
            resname = match.group(2).upper()

        residue_key = (resid, resname)
        if residue_key != last_residue:
            residue_counts[resname] += 1
            last_residue = residue_key

        residue_atom_counts[resname] += 1

    box_line = lines[atom_end].strip()
    raw_values = box_line.split()

    if len(raw_values) not in {3, 9}:
        raise ValueError(f"GRO box line has {len(raw_values)} values; expected 3 or 9.")

    try:
        values = [float(v.replace("D", "E").replace("d", "e")) for v in raw_values]
    except ValueError as exc:
        raise ValueError("Could not parse GRO box-vector values.") from exc

    zero = "0"
    if len(values) == 3:
        a = (values[0], 0.0, 0.0)
        b = (0.0, values[1], 0.0)
        c = (0.0, 0.0, values[2])
        token_vectors = (
            (raw_values[0], zero, zero),
            (zero, raw_values[1], zero),
            (zero, zero, raw_values[2]),
        )
    else:
        a = (values[0], values[3], values[4])
        b = (values[5], values[1], values[6])
        c = (values[7], values[8], values[2])
        token_vectors = (
            (raw_values[0], raw_values[3], raw_values[4]),
            (raw_values[5], raw_values[1], raw_values[6]),
            (raw_values[7], raw_values[8], raw_values[2]),
        )

    off_diagonal = (abs(a[1]), abs(a[2]), abs(b[0]), abs(b[2]), abs(c[0]), abs(c[1]))
    lengths = (vector_norm(a), vector_norm(b), vector_norm(c))
    tolerance = 1.0e-8

    if any(val > tolerance for val in off_diagonal):
        inferred_type = "triclinic"
    elif abs(lengths[0] - lengths[1]) < tolerance and abs(lengths[1] - lengths[2]) < tolerance:
        inferred_type = "cubic/orthorhombic"
    else:
        inferred_type = "orthorhombic"

    denom = lengths[0] * lengths[1]
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, dot(a, b) / denom)))) if denom > 0.0 else float("nan")
    xy_area = vector_norm(cross(a, b))
    volume = abs(dot(a, cross(b, c)))

    extra_lines = [l for l in lines[atom_end + 1:] if l.strip()]
    if extra_lines:
        warnings.append(f"{len(extra_lines)} non-empty line(s) occur after the GRO box line.")

    return GroInfo(
        title=title,
        declared_atom_count=declared_atom_count,
        parsed_atom_count=len(atom_lines),
        residue_counts=residue_counts,
        residue_atom_counts=residue_atom_counts,
        raw_box_line=box_line,
        raw_box_values=raw_values,
        vectors=(a, b, c),
        vector_tokens=token_vectors,
        inferred_box_type=inferred_type,
        vector_lengths=lengths,
        gamma_degrees=gamma,
        xy_area_nm2=xy_area,
        volume_nm3=volume,
        warnings=warnings,
    )


def read_mdp_settings(path: Path) -> dict[str, str]:
    settings = {}
    for raw in read_text(path).splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        settings[k.strip().lower()] = v.strip()
    return settings


def collect_existing_commands(root: Path, gro_path: Optional[Path]) -> list[tuple[Path, str]]:
    candidates = []
    search_dirs = {root}
    if gro_path is not None:
        search_dirs.add(gro_path.parent)

    for directory in search_dirs:
        if not directory.exists():
            continue
        candidates.extend(directory.glob("README*"))
        candidates.extend(directory.glob("readme*"))
        candidates.extend(directory.glob("*.sh"))

    commands = []
    seen_paths = set()

    for path in sorted(set(candidates), key=lambda p: natural_key(str(p))):
        if not path.is_file() or path.resolve() in seen_paths:
            continue
        seen_paths.add(path.resolve())

        logical_lines = []
        buf = ""
        for raw in read_text(path).splitlines():
            stripped = raw.strip()
            if not stripped:
                continue
            buf = buf + " " + stripped if buf else stripped
            if buf.endswith("\\"):
                buf = buf[:-1].rstrip()
                continue
            logical_lines.append(buf)
            buf = ""
        if buf:
            logical_lines.append(buf)

        for line in logical_lines:
            cleaned = re.sub(r"^\s*\$\s*", "", line).strip()
            if re.search(r"(^|\s)(?:gmx(?:_mpi)?|grompp|mdrun|mpirun.*gmx|charmm|namd\d*|./\S+)", cleaned, re.IGNORECASE):
                commands.append((path, cleaned))

    return commands[:250]


def generate_gromacs_commands(root: Path, gro_path: Optional[Path]) -> tuple[Optional[Path], list[Path], list[str]]:
    if gro_path is None:
        return None, [], []

    run_dir = gro_path.parent
    mdp_files = sorted(run_dir.glob("*.mdp"), key=lambda p: natural_key(p.name))
    if not mdp_files:
        return run_dir, [], []

    topol_path = run_dir / "topol.top"
    topology = topol_path.name if topol_path.exists() else "topol.top"

    commands = []
    prev_gro = gro_path.name
    prev_cpt = None
    ref_gro = gro_path.name

    for mdp in mdp_files:
        settings = read_mdp_settings(mdp)
        stem = mdp.stem
        define = settings.get("define", "")
        continuation = settings.get("continuation", "no").lower()
        refcoord_scaling = settings.get("refcoord_scaling", "").lower()

        req_ref = "POSRES" in define.upper() or refcoord_scaling not in {"", "no", "none"}

        cmd = ["gmx", "grompp", "-f", mdp.name, "-o", f"{stem}.tpr", "-c", prev_gro]
        if req_ref:
            cmd.extend(["-r", ref_gro])
        if continuation in {"yes", "true"} and prev_cpt:
            cmd.extend(["-t", prev_cpt])
        cmd.extend(["-p", topology])

        commands.append(" ".join(shlex.quote(p) for p in cmd))
        commands.append(f"gmx mdrun -deffnm {shlex.quote(stem)}")

        prev_gro = f"{stem}.gro"
        prev_cpt = f"{stem}.cpt"

    return run_dir, mdp_files, commands


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "N/A"
    return f"{value:.10g}"


def append_evidence(output: list[str], records: list[Evidence], root: Path, empty_msg: str = "NOT FOUND", max_recs: int = 80):
    if not records:
        output.append(f"  {empty_msg}")
        return
    for r in records[:max_recs]:
        output.append(f"  - {relative_path(r.path, root)}:{r.line_number}: {r.text.strip()}")
    if len(records) > max_recs:
        output.append(f"  ... {len(records) - max_recs} additional matching lines omitted.")


def create_report(
    root: Path,
    selected: dict[str, Optional[Path]],
    detection_warnings: list[str],
    extra_water_names: set[str],
    extra_ion_names: set[str],
    extra_lipid_names: set[str],
    include_run_template: bool,
) -> str:
    warnings = list(detection_warnings)
    config_evidence = collect_config_evidence(selected["config"])

    upper = parse_lipid_file(selected["lipid_upper"])
    lower = parse_lipid_file(selected["lipid_lower"])
    warnings.extend(upper.warnings)
    warnings.extend(lower.warnings)

    component_declarations = collect_component_declarations(selected["components"])
    ion_info = parse_ion_information([selected["ions"], selected["config"], selected["components"]])
    orientation = collect_orientation(selected["orientation"])

    area_strong, area_weak = collect_area_evidence([
        selected["config"], selected["lipid_upper"], selected["lipid_lower"],
        selected["components"], selected["ions"], selected["orientation"]
    ])
    calc_box_evidence = collect_calculated_box_evidence([
        selected["config"], selected["lipid_upper"], selected["lipid_lower"], selected["components"]
    ])

    gro = None
    if selected["gro"] is not None:
        try:
            gro = parse_gro(selected["gro"])
            warnings.extend(gro.warnings)
        except Exception as exc:
            warnings.append(f"Could not parse {selected['gro'].name}: {exc}")

    water_names = {n.upper() for n in DEFAULT_WATER_NAMES | extra_water_names}
    ion_names = {n.upper() for n in DEFAULT_ION_NAMES | extra_ion_names}
    ion_names.update(v.upper() for v in ion_info.types.values())

    parsed_lipid_names = set(upper.counts) | set(lower.counts)
    lipid_names = {n.upper() for n in (COMMON_LIPID_NAMES | parsed_lipid_names | extra_lipid_names)}

    output = []
    output.append("=" * 79)
    output.append("CHARMM-GUI SYSTEM INFORMATION REPORT")
    output.append("=" * 79)
    output.append(f"Search root: {root.resolve()}")
    output.append("")

    output.append("1. DETECTED INPUT FILES")
    output.append("-" * 79)
    for label, exp_name in TARGET_FILES.items():
        p = selected[label]
        output.append(f"  {exp_name:<32} {relative_path(p, root) if p else 'NOT FOUND'}")
    output.append("")

    output.append("2. BOX AND SOLVATION SETTINGS")
    output.append("-" * 79)
    output.append("Specified box type:")
    append_evidence(output, config_evidence["box_type"], root)
    output.append("\nSpecified box XY dimensions/settings:")
    append_evidence(output, config_evidence["box_xy"], root)
    output.append("\nSpecified water thickness/padding:")
    append_evidence(output, config_evidence["water_thickness"], root)
    output.append("\nCalculated/final box declarations found in builder files:")
    append_evidence(output, unique_evidence(config_evidence["final_box"] + calc_box_evidence), root)

    output.append("\nDefinitive final box from step5_input.gro:")
    if gro:
        la, lb, lc = gro.vector_lengths
        output.append(f"  Inferred box type: {gro.inferred_box_type}")
        output.append(f"  Lattice lengths |a| x |b| x |c|: {format_float(la)} x {format_float(lb)} x {format_float(lc)} nm")
        output.append(f"  Lattice lengths in Angstrom: {format_float(la * 10.0)} x {format_float(lb * 10.0)} x {format_float(lc * 10.0)} A")
        output.append(f"  Final XY/lateral size: |a| x |b| = {format_float(la)} x {format_float(lb)} nm")
        output.append(f"  Final XY/lateral size: {format_float(la * 10.0)} x {format_float(lb * 10.0)} A")
        output.append(f"  Angle gamma(a,b): {format_float(gro.gamma_degrees)} degrees")
        output.append(f"  XY lattice area |a x b|: {format_float(gro.xy_area_nm2)} nm^2 ({format_float(gro.xy_area_nm2 * 100.0)} A^2)")
        output.append(f"  Box volume: {format_float(gro.volume_nm3)} nm^3")
    else:
        output.append("  NOT AVAILABLE: step5_input.gro was not parsed.")

    output.append("\n3. LIPID COUNTS BY LEAFLET AND SPECIES")
    output.append("-" * 79)
    species = sorted(set(upper.counts) | set(lower.counts), key=natural_key)
    if species:
        output.append(f"  {'Species':<14}{'Upper':>10}{'Lower':>10}{'Planned':>12}{'Final GRO':>12}{'Delta':>10}")
        output.append(f"  {'-' * 14}{'-' * 10}{'-' * 10}{'-' * 12}{'-' * 12}{'-' * 10}")
        for lip in species:
            u_cnt = upper.counts.get(lip, 0)
            l_cnt = lower.counts.get(lip, 0)
            pln = u_cnt + l_cnt
            if gro:
                f_cnt = gro.residue_counts.get(lip, 0)
                d_text = f"{f_cnt - pln:+d}"
                f_text = str(f_cnt)
            else:
                f_text, d_text = "-", "-"
            output.append(f"  {lip:<14}{u_cnt:>10}{l_cnt:>10}{pln:>12}{f_text:>12}{d_text:>10}")
    else:
        output.append("  Per-species lipid counts: NOT FOUND")

    u_tot = str(upper.total) if upper.total is not None else "NOT FOUND"
    l_tot = str(lower.total) if lower.total is not None else "NOT FOUND"
    output.append(f"\n  Upper-leaflet total: {u_tot}\n  Lower-leaflet total: {l_tot}")
    if upper.total is not None and lower.total is not None:
        output.append(f"  Planned lipid total: {upper.total + lower.total}")

    output.append("\nUpper-leaflet source records:")
    append_evidence(output, upper.evidence, root)
    output.append("\nLower-leaflet source records:")
    append_evidence(output, lower.evidence, root)

    output.append("\n4. CONFIRMED FINAL COMPONENT COUNTS")
    output.append("-" * 79)
    if gro:
        f_lipids = {k: v for k, v in gro.residue_counts.items() if k in lipid_names}
        f_waters = {k: v for k, v in gro.residue_counts.items() if k in water_names}
        f_ions = {k: v for k, v in gro.residue_counts.items() if k in ion_names}

        output.append(f"  Lipid molecules: {sum(f_lipids.values())}")
        for k in sorted(f_lipids, key=natural_key):
            output.append(f"    {k:<12} {f_lipids[k]}")

        output.append(f"  Water molecules: {sum(f_waters.values())}")
        for k in sorted(f_waters, key=natural_key):
            output.append(f"    {k:<12} {f_waters[k]}")

        output.append(f"  Ions: {sum(f_ions.values())}")
        for k in sorted(f_ions, key=natural_key):
            output.append(f"    {k:<12} {f_ions[k]}")

        classified = sum(f_lipids.values()) + sum(f_waters.values()) + sum(f_ions.values())
        output.append(f"  Other residue records (protein/peptide/other): {sum(gro.residue_counts.values()) - classified}")
    else:
        output.append("  NOT AVAILABLE: final GRO file was not parsed.")

    output.append("\nComponent-count declarations from step4_components.str:")
    append_evidence(output, component_declarations, root)

    output.append("\n5. WATER INFORMATION")
    output.append("-" * 79)
    if gro:
        w_counts = {k: v for k, v in gro.residue_counts.items() if k in water_names}
        output.append(f"  Confirmed number of water molecules: {sum(w_counts.values())}")
        output.append(f"  Water residue names: {', '.join(sorted(w_counts)) if w_counts else 'None found'}")
    else:
        output.append("  Confirmed number of water molecules: NOT AVAILABLE")
    output.append("\nSpecified water thickness source:")
    append_evidence(output, config_evidence["water_thickness"], root)

    output.append("\n6. ION TYPE, COUNT, AND CONCENTRATION")
    output.append("-" * 79)
    output.append("Declared ion types:")
    if ion_info.types:
        for lbl, ityp in ion_info.types.items():
            output.append(f"  {lbl}: {ityp}")
    else:
        output.append("  NOT FOUND")

    output.append("\nIon type source records:")
    append_evidence(output, ion_info.type_evidence, root)
    output.append("\nDeclared ion count records:")
    append_evidence(output, ion_info.count_evidence, root)
    output.append("\nDeclared ion concentration records:")
    append_evidence(output, ion_info.concentration_evidence, root)

    output.append("\nConfirmed ion counts from final GRO:")
    if gro:
        gro_ions = {k: v for k, v in gro.residue_counts.items() if k in ion_names}
        if gro_ions:
            for k in sorted(gro_ions, key=natural_key):
                output.append(f"  {k}: {gro_ions[k]}")
            output.append(f"  Total ions: {sum(gro_ions.values())}")
        else:
            output.append("  No recognized ion residues found.")
    else:
        output.append("  NOT AVAILABLE")

    output.append("\n7. PEPTIDE ORIENTATION / PPM INFORMATION")
    output.append("-" * 79)
    output.append("Tilt angle records:")
    append_evidence(output, orientation["tilt"], root)
    output.append("\nInsertion depth / membrane-position records:")
    append_evidence(output, orientation["insertion_depth"], root)
    output.append("\nOther PPM/orientation records:")
    append_evidence(output, orientation["other_ppm"], root)

    output.append("\n8. PEPTIDE CROSS-SECTIONAL / LIPID-EXCLUSION AREA")
    output.append("-" * 79)
    if area_strong:
        append_evidence(output, area_strong, root)
    elif area_weak:
        output.append("  No unambiguous cross-sectional area found. Other candidate lines:")
        append_evidence(output, area_weak, root)
    else:
        output.append("  NOT FOUND")

    output.append("\n9. TOTAL ATOM COUNT AND ALL RESIDUE COUNTS")
    output.append("-" * 79)
    if gro:
        output.append(f"  GRO title: {gro.title}")
        output.append(f"  Declared total atom count: {gro.declared_atom_count}")
        output.append(f"  Parsed atom-line count: {gro.parsed_atom_count}")
        output.append(f"  Total contiguous residue records: {sum(gro.residue_counts.values())}")
        output.append(f"\n  {'Residue':<14}{'Residues/Molecules':>20}{'Atoms':>14}")
        output.append(f"  {'-' * 14}{'-' * 20}{'-' * 14}")
        for r in sorted(gro.residue_counts, key=natural_key):
            output.append(f"  {r:<14}{gro.residue_counts[r]:>20}{gro.residue_atom_counts[r]:>14}")
    else:
        output.append("  NOT AVAILABLE")

    output.append("\n10. EXACT FINAL BOX VECTORS")
    output.append("-" * 79)
    if gro:
        ta, tb, tc = gro.vector_tokens
        a, b, c = gro.vectors
        output.append(f"  Raw GRO box line: {gro.raw_box_line}")
        output.append(f"\n  GROMACS lattice vectors in nm:\n    a = ({ta[0]}, {ta[1]}, {ta[2]})\n    b = ({tb[0]}, {tb[1]}, {tb[2]})\n    c = ({tc[0]}, {tc[1]}, {tc[2]})")
        output.append(f"\n  Lattice vectors converted to Angstrom:\n    a = ({', '.join(format_float(v * 10.0) for v in a)})\n    b = ({', '.join(format_float(v * 10.0) for v in b)})\n    c = ({', '.join(format_float(v * 10.0) for v in c)})")
    else:
        output.append("  NOT AVAILABLE")

    output.append("\n11. RUN INSTRUCTIONS / STEP SEQUENCE")
    output.append("-" * 79)
    conceptual = [
        ("config", "Input/build configuration"),
        ("orientation", "PPM/orientation result"),
        ("lipid_upper", "Upper-leaflet lipid placement"),
        ("lipid_lower", "Lower-leaflet lipid placement"),
        ("components", "Component assembly definition"),
        ("ions", "Ion replacement"),
        ("gro", "Final assembled coordinates"),
    ]
    seq_idx = 1
    for k, desc in conceptual:
        if selected[k]:
            output.append(f"  {seq_idx}. {desc}: {relative_path(selected[k], root)}")
            seq_idx += 1

    cmds = collect_existing_commands(root, selected["gro"])
    output.append("\nCommands extracted from README/shell files:")
    if cmds:
        for p, c in cmds:
            output.append(f"  [{relative_path(p, root)}] {c}")
    else:
        output.append("  No explicit command lines detected.")

    r_dir, mdps, gen_cmds = generate_gromacs_commands(root, selected["gro"])
    output.append("\nDetected GROMACS MDP stage order:")
    if mdps:
        for i, m in enumerate(mdps, 1):
            output.append(f"  {i}. {m.name}")
    else:
        output.append("  No MDP files detected beside step5_input.gro.")

    if include_run_template and r_dir and gen_cmds:
        output.append("\nGenerated GROMACS command template:")
        output.append(f"  cd {shlex.quote(relative_path(r_dir, root) or '.')}")
        for c in gen_cmds:
            output.append(f"  {c}")

    output.append("\n12. WARNINGS / EXTRACTION NOTES")
    output.append("-" * 79)
    if warnings:
        for w in dict.fromkeys(warnings):
            output.append(f"  - {w}")
    else:
        output.append("  No parser warnings.")

    output.append("")
    return "\n".join(output)


def parse_name_set(value: str) -> set[str]:
    return {item.strip().upper() for item in value.split(",") if item.strip()} if value else set()


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract CHARMM-GUI membrane metadata into a report.")
    parser.add_argument("root", nargs="?", default=".", help="Project directory to search recursively.")
    parser.add_argument("-o", "--output", default="system_information.txt", help="Output report file path.")
    parser.add_argument("--water-resnames", default="", help="Additional water residue names (comma-separated).")
    parser.add_argument("--ion-resnames", default="", help="Additional ion residue names (comma-separated).")
    parser.add_argument("--lipid-resnames", default="", help="Additional lipid residue names (comma-separated).")
    parser.add_argument("--no-run-template", action="store_true", help="Do not generate GROMACS execution template.")

    args = parser.parse_args()
    root = Path(args.root).expanduser().resolve()

    if not root.is_dir():
        print(f"ERROR: Search root is not a valid directory: {root}", file=sys.stderr)
        return 2

    selected, detection_warnings = detect_files(root)

    report_content = create_report(
        root=root,
        selected=selected,
        detection_warnings=detection_warnings,
        extra_water_names=parse_name_set(args.water_resnames),
        extra_ion_names=parse_name_set(args.ion_resnames),
        extra_lipid_names=parse_name_set(args.lipid_resnames),
        include_run_template=not args.no_run_template,
    )

    output_path = Path(args.output).expanduser()
    if not output_path.is_absolute():
        output_path = root / output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report_content, encoding="utf-8")

    print(f"Report successfully written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())