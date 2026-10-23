import argparse
import csv
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import date, timedelta
from pathlib import Path


XLSX_NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
SKIP_DIRS = {"omezarr", "work", "logs", "extracted_tifs", "scene_inputs", "pyramid", "metadata", "additional data"}
NON_SAMPLE_HEADERS = {"date", "time", "supposed_tim", "supposed_time", "step", "chemical"}
END_SAMPLE_HEADERS = {"amount", "temp", "where", "comments"}
TREATMENT_TOKENS = {"t", "treatment", "treatmen"}
CONTROL_TOKENS = {"c", "control", "scrambled"}

FIELDNAMES = [
    "run",
    "dataset_name",
    "input_file",
    "source_file",
    "format",
    "run_date",
    "sample",
    "condition",
    "sheet_found",
    "sample_match",
]


def sanitize_name(value):
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_")
    return cleaned or "osfish_lif"


def compact(value):
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def tokens(value):
    return [token for token in re.split(r"[^a-z0-9]+", str(value or "").lower()) if token]


def column_index(cell_ref):
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    index = 0
    for letter in letters:
        index = index * 26 + (ord(letter) - ord("A") + 1)
    return index - 1


def read_xlsx(path):
    """Return {sheet_name: [[cell, ...], ...]} using only the standard library."""
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", XLSX_NS):
                shared.append("".join(t.text or "" for t in item.iter(f"{{{XLSX_NS['m']}}}t")))

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {rel.get("Id"): rel.get("Target") for rel in rels}

        sheets = {}
        for sheet in workbook.find("m:sheets", XLSX_NS):
            target = targets[sheet.get(f"{{{XLSX_NS['r']}}}id")].lstrip("/")
            if not target.startswith("xl/"):
                target = f"xl/{target}"
            root = ET.fromstring(archive.read(target))
            rows = []
            for row in root.iter(f"{{{XLSX_NS['m']}}}row"):
                values = {}
                for cell in row.findall("m:c", XLSX_NS):
                    value = cell.find("m:v", XLSX_NS)
                    cell_type = cell.get("t")
                    if cell_type == "s" and value is not None:
                        text = shared[int(value.text)]
                    elif cell_type == "inlineStr":
                        text = "".join(t.text or "" for t in cell.iter(f"{{{XLSX_NS['m']}}}t"))
                    else:
                        text = value.text if value is not None else ""
                    values[column_index(cell.get("r"))] = (text or "").strip()
                width = max(values) + 1 if values else 0
                rows.append([values.get(index, "") for index in range(width)])
            sheets[sheet.get("name")] = rows
    return sheets


def normalize_date(value):
    value = str(value or "").strip()
    if not value:
        return ""
    match = re.fullmatch(r"(\d{1,2})\.(\d{1,2})\.(\d{4})", value)
    if match:
        day, month, year = match.groups()
        return f"{year}-{int(month):02d}-{int(day):02d}"
    try:
        serial = float(value)
    except ValueError:
        return ""
    if serial > 20000:
        return (date(1899, 12, 30) + timedelta(days=int(serial))).isoformat()
    return ""


def sheet_info(rows):
    header = rows[0] if rows else []
    lowered = [cell.lower() for cell in header]
    start = lowered.index("step") + 1 if "step" in lowered else 0
    samples = []
    for cell in header[start:]:
        if cell.lower() in END_SAMPLE_HEADERS:
            break
        if cell and cell.lower() not in NON_SAMPLE_HEADERS:
            samples.append(cell)
    run_date = ""
    for row in rows[1:]:
        run_date = normalize_date(row[0] if row else "")
        if run_date:
            break
    return {"samples": samples, "run_date": run_date}


def numbered_sample(stem, run, samples):
    """Match Zeiss-style names like SF02_03_01 (sample 3, image 1) to a column s3_..."""
    rest = stem[len(run):] if stem.lower().startswith(run.lower()) else stem
    match = re.match(r"[^0-9]*0*(\d+)[^0-9]", rest + "_")
    if not match:
        return ("", "")
    for sample in samples:
        if tokens(sample)[:1] == [f"s{int(match.group(1))}"]:
            return (sample, "number")
    return ("", "")


def match_sample(lif_stem, samples):
    """Pick the sheet sample column that best matches a file name."""
    lif_compact = compact(lif_stem)
    lif_tokens = set(tokens(lif_stem))
    best = ("", "")
    best_score = 0
    for sample in samples:
        sample_compact = compact(sample)
        sample_tokens = set(tokens(sample))
        if not sample_compact:
            continue
        if sample_compact in lif_compact:
            score, how = 1000 + len(sample_compact), "name"
        elif sample_tokens and sample_tokens <= lif_tokens:
            score, how = 500 + len(sample_tokens), "tokens"
        else:
            continue
        if score > best_score:
            best, best_score = (sample, how), score
    return best


def condition_from(*values):
    for value in values:
        value_tokens = set(tokens(value))
        is_treatment = bool(value_tokens & TREATMENT_TOKENS)
        is_control = bool(value_tokens & CONTROL_TOKENS)
        if is_treatment and not is_control:
            return "treatment"
        if is_control and not is_treatment:
            return "control"
    return ""


def file_format(path):
    name = path.name.lower()
    if name.endswith(".lif"):
        return "lif"
    if name.endswith(".czi"):
        return "czi"
    if name.endswith(".ome.tif") or name.endswith(".ome.tiff"):
        return "ome-tiff"
    return ""


def file_stem(path):
    return re.sub(r"(?i)\.ome$", "", path.stem)


def run_of(relative):
    return relative.parts[0] if len(relative.parts) > 1 else file_stem(relative)


def find_inputs(runs_dir):
    """LIF and CZI files everywhere; OME-TIFF only in runs without LIFs.

    Plain TIFFs and OME-TIFFs next to LIFs are exports of LIF scenes, so they are
    skipped. A raw OME-TIFF tile is skipped when its *_Merged version exists.
    """
    candidates = []
    for path in sorted(runs_dir.rglob("*")):
        relative = path.relative_to(runs_dir)
        if not path.is_file() or any(part.lower() in SKIP_DIRS for part in relative.parts[:-1]):
            continue
        kind = file_format(path)
        if kind:
            candidates.append((path, relative, kind))

    runs_with_lif = {run_of(relative) for _, relative, kind in candidates if kind == "lif"}
    selected = []
    for path, relative, kind in candidates:
        if kind == "ome-tiff":
            if run_of(relative) in runs_with_lif:
                continue
            merged = [path.with_name(f"{file_stem(path)}{sep}Merged{ext}") for sep in ("_", " ") for ext in (".ome.tif", ".ome.tiff")]
            if any(candidate.exists() for candidate in merged):
                continue
        selected.append((path, relative, kind))
    return drop_cross_run_copies(selected)


def drop_cross_run_copies(selected):
    """Keep one copy of files duplicated into another run (same folder, name and size)."""
    def key(item):
        path, relative, _ = item
        return (path.parent.name, path.name, path.stat().st_size)

    runs_by_key = {}
    folder_counts = {}
    for item in selected:
        runs_by_key.setdefault(key(item), set()).add(run_of(item[1]))
        folder = (run_of(item[1]), item[0].parent.name)
        folder_counts[folder] = folder_counts.get(folder, 0) + 1

    kept = []
    for item in selected:
        runs = runs_by_key[key(item)]
        run = run_of(item[1])
        if len(runs) > 1:
            folder = item[0].parent.name
            best = max(sorted(runs), key=lambda name: folder_counts.get((name, folder), 0))
            if run != best:
                print(f"Skipping {item[1]}: same file as in run {best}", file=sys.stderr)
                continue
        kept.append(item)
    return kept


def run_alias(folder_run, stems, sheets_by_name):
    """A folder without its own sheet (e.g. 13052025) belongs to the sheet all its files are named after."""
    if folder_run.lower() in sheets_by_name:
        return folder_run
    owners = {
        name for name in sheets_by_name
        if stems and all(stem.lower().startswith(f"{name}_") for stem in stems)
    }
    return next(iter(owners)) if len(owners) == 1 else folder_run


def dataset_names(entries):
    """Name each LIF <run>_<stem>, adding parent folders only when stems collide."""
    by_key = {}
    for entry in entries:
        by_key.setdefault((entry["run"], entry["stem"]), []).append(entry)
    names = {}
    for (run, stem), group in by_key.items():
        for entry in group:
            parts = entry["parents"] + [stem] if len(group) > 1 else [stem]
            name = "_".join(parts)
            if not compact(name).startswith(compact(run)):
                name = f"{run}_{name}"
            names[entry["input_file"]] = sanitize_name(name)
    if len(set(names.values())) != len(names):
        raise RuntimeError("Generated dataset names are not unique; rename LIF files or edit the inventory")
    return names


def main():
    parser = argparse.ArgumentParser(
        description="Inventory all LIF/CZI/OME-TIFF files below a runs directory and link them to samples.xlsx."
    )
    parser.add_argument("--runs-dir", required=True)
    parser.add_argument("--samples-xlsx", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--runs", default="", help="Optional comma-separated run names to include.")
    args = parser.parse_args()

    runs_dir = Path(args.runs_dir)
    sheets = read_xlsx(args.samples_xlsx) if args.samples_xlsx else {}
    sheets_by_name = {name.lower(): sheet_info(rows) for name, rows in sheets.items()}
    selected_runs = {run.strip().lower() for run in args.runs.split(",") if run.strip()}

    inputs = find_inputs(runs_dir)
    stems_by_folder = {}
    for path, relative, _ in inputs:
        stems_by_folder.setdefault(run_of(relative), []).append(file_stem(path))
    sheet_names = {name.lower(): name for name in sheets}
    aliases = {}
    for folder_run, stems in stems_by_folder.items():
        alias = run_alias(folder_run, stems, sheets_by_name)
        aliases[folder_run] = sheet_names.get(alias, alias)
        if aliases[folder_run] != folder_run:
            print(f"Folder {folder_run} has no sheet; using sheet {aliases[folder_run]} (all files are named {aliases[folder_run]}_*)", file=sys.stderr)

    entries = []
    for path, relative, kind in inputs:
        run = aliases[run_of(relative)]
        if selected_runs and run.lower() not in selected_runs:
            continue
        entries.append(
            {
                "run": run,
                "stem": file_stem(path),
                "parents": list(relative.parts[1:-1]),
                "input_file": str(path),
                "source_file": str(relative),
                "format": kind,
            }
        )

    names = dataset_names(entries)
    rows = []
    for entry in entries:
        info = sheets_by_name.get(entry["run"].lower())
        sample, how = match_sample(entry["stem"], info["samples"]) if info else ("", "")
        if info and not sample:
            sample, how = numbered_sample(entry["stem"], entry["run"], info["samples"])
        rows.append(
            {
                "run": entry["run"],
                "dataset_name": names[entry["input_file"]],
                "input_file": entry["input_file"],
                "source_file": entry["source_file"],
                "format": entry["format"],
                "run_date": info["run_date"] if info else "",
                "sample": sample,
                "condition": condition_from(entry["stem"], sample, *entry["parents"]),
                "sheet_found": "TRUE" if info else "FALSE",
                "sample_match": how,
            }
        )

    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    unmatched = sum(1 for row in rows if not row["sample"])
    missing_sheets = sorted({row["run"] for row in rows if row["sheet_found"] == "FALSE"})
    print(f"Wrote {len(rows)} files to {args.output} ({unmatched} without a sample match)", file=sys.stderr)
    if missing_sheets:
        print(f"Runs without a samples.xlsx sheet: {', '.join(missing_sheets)}", file=sys.stderr)


if __name__ == "__main__":
    main()
