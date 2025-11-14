#!/usr/bin/env python3
"""Utility helpers for quickly inspecting and transforming photometry tables."""

from __future__ import annotations

import argparse
import configparser
import csv
import json
import math
import statistics
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


COLOR_KEY_MAP = {
    "bv": "Tv_bv",
    "vr": "Tv_vr",
    "vi": "Tv_vi",
}


def read_table(table_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with table_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
    except FileNotFoundError:
        sys.stderr.write(f"error: file not found -> {table_path}\n")
        sys.exit(1)

    if not rows:
        sys.stderr.write(f"error: file is empty -> {table_path}\n")
        sys.exit(1)

    header = reader.fieldnames or []
    return header, rows


def collect_star_labels(header: Iterable[str]) -> list[str]:
    labels = []
    for column in header:
        if column.startswith("Source_AMag_"):
            labels.append(column.removeprefix("Source_AMag_"))
    # Preserve natural ordering by preferring target/check before comparisons.
    def sort_key(label: str) -> tuple[int, str]:
        if label == "T1":
            return (0, label)
        if label == "T2":
            return (1, label)
        if label.startswith("C"):
            return (2, label)
        return (3, label)

    return sorted(labels, key=sort_key)


def report_star_presence(header: Iterable[str]) -> None:
    has_t1 = "Source_AMag_T1" in header
    has_t2 = "Source_AMag_T2" in header
    comparison_columns = sorted(col for col in header if col.startswith("Source_AMag_C"))

    print(f"Target T1 present: {'yes' if has_t1 else 'no'}")
    print(f"Check star T2 present: {'yes' if has_t2 else 'no'}")
    print(f"Comparison stars found: {len(comparison_columns)}")
    if comparison_columns:
        ids = ", ".join(comparison_columns)
        print(f"Comparison columns: {ids}")


def extract_target_position(rows: Iterable[Mapping[str, str]]) -> tuple[float, float]:
    for row in rows:
        try:
            ra_val = float(row["RA_T1"])
            dec_val = float(row["DEC_T1"])
        except (KeyError, TypeError, ValueError):
            continue
        else:
            # RA appears to be stored in hours; convert to degrees.
            ra_deg = ra_val * 15.0 if ra_val <= 24.0 else ra_val
            return ra_deg, dec_val
    raise RuntimeError("unable to locate numeric RA_T1/DEC_T1 entries in table")


def format_ra_deg(ra_deg: float) -> str:
    total_seconds = ra_deg / 15.0 * 3600.0
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds % 60
    return f"{hours:02d}h {minutes:02d}m {seconds:05.2f}s"


def format_dec_deg(dec_deg: float) -> str:
    sign = "-" if dec_deg < 0 else "+"
    abs_deg = abs(dec_deg)
    degrees = int(abs_deg)
    arcminutes = int((abs_deg - degrees) * 60)
    arcseconds = (abs_deg - degrees - arcminutes / 60) * 3600
    return f"{sign}{degrees:02d}° {arcminutes:02d}' {arcseconds:04.1f}\""


def query_aavso_comparisons(
    ra_deg: float,
    dec_deg: float,
    radius_arcmin: float,
    mag_limit: float,
) -> dict:
    # AAVSO VSP API expects specific format - using 'fov' parameter for field of view
    # Note: API may return HTML error page if parameters are incorrect
    params = {
        "format": "json",
        "fov": str(radius_arcmin),
        "maglimit": str(mag_limit),
        "ra": str(ra_deg),
        "dec": str(dec_deg),
    }
    url = "https://apps.aavso.org/vsp/api/chart/?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "stdphot-tool/1.0",
            "Accept": "application/json",
        }
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            content = response.read().decode('utf-8')
            # Check if we got HTML instead of JSON (API error)
            if content.strip().startswith('<!DOCTYPE') or content.strip().startswith('<html'):
                sys.stderr.write(f"error: AAVSO API returned HTML instead of JSON\n")
                sys.stderr.write(f"This usually means invalid parameters or API endpoint\n")
                sys.stderr.write(f"URL requested: {url}\n")
                sys.exit(1)
            try:
                return json.loads(content)
            except json.JSONDecodeError as json_err:
                sys.stderr.write(f"error: failed to parse AAVSO response as JSON\n")
                sys.stderr.write(f"Response content (first 500 chars): {content[:500]}\n")
                sys.stderr.write(f"JSON error: {json_err}\n")
                sys.exit(1)
    except urllib.error.HTTPError as exc:
        sys.stderr.write(f"error: AAVSO VSP API returned HTTP {exc.code}: {exc.reason}\n")
        try:
            error_content = exc.read().decode('utf-8')
            sys.stderr.write(f"Error response: {error_content[:500]}\n")
        except Exception:
            pass
        sys.exit(1)
    except urllib.error.URLError as exc:
        sys.stderr.write(f"error: failed to reach AAVSO VSP API: {exc}\n")
        sys.exit(1)


def emit_aavso_summary(payload: Mapping[str, object]) -> None:
    comp_stars = payload.get("photometry")
    if not isinstance(comp_stars, list):
        print("No comparison star data returned from AAVSO.")
        return

    if not comp_stars:
        print("No comparison stars found in search area.")
        return

    print(f"AAVSO comparison stars retrieved: {len(comp_stars)}")
    for entry in comp_stars:
        if not isinstance(entry, Mapping):
            continue
        try:
            # Extract label and AUID
            label = entry.get("label", "?")
            auid = entry.get("auid", "")
            
            # Extract RA/Dec (in sexagesimal format from AAVSO)
            ra_str = entry.get("ra", "")
            dec_str = entry.get("dec", "")
            
            # Extract bands list
            bands = entry.get("bands", [])
            if not isinstance(bands, list):
                continue
            
            # Build a dict of band -> magnitude
            band_dict = {}
            for band_entry in bands:
                if isinstance(band_entry, Mapping):
                    band_name = band_entry.get("band")
                    mag_value = band_entry.get("mag")
                    error_value = band_entry.get("error")
                    if band_name and mag_value is not None:
                        band_dict[band_name] = {
                            "mag": float(mag_value),
                            "error": float(error_value) if error_value is not None else None
                        }
            
            # Get V magnitude
            v_mag = band_dict.get("V", {}).get("mag", float("nan"))
            v_err = band_dict.get("V", {}).get("error")
            
            # Calculate B-V color if both available
            b_v_color = "?"
            if "B" in band_dict and "V" in band_dict:
                b_mag = band_dict["B"]["mag"]
                v_mag_val = band_dict["V"]["mag"]
                b_v_color = f"{b_mag - v_mag_val:.3f}"
            
            # Calculate V-I color if both available
            v_i_color = "?"
            if "V" in band_dict and "Ic" in band_dict:
                v_mag_val = band_dict["V"]["mag"]
                i_mag = band_dict["Ic"]["mag"]
                v_i_color = f"{v_mag_val - i_mag:.3f}"
            
        except (TypeError, ValueError, KeyError) as e:
            continue
        
        # Format output line
        v_str = f"{v_mag:.3f}" if not (v_mag != v_mag) else "  -  "
        err_str = f"±{v_err:.3f}" if v_err is not None else ""
        print(f"  {label:>4s} ({auid:15s})  V={v_str}{err_str:8s}  B-V={b_v_color:>6s}  V-I={v_i_color:>6s}  {ra_str} {dec_str}")


def parse_sexagesimal_ra(ra_str: str) -> float:
    """Convert RA from HH:MM:SS.SS to decimal degrees."""
    parts = ra_str.split(':')
    if len(parts) != 3:
        return float('nan')
    hours = float(parts[0])
    minutes = float(parts[1])
    seconds = float(parts[2])
    return (hours + minutes / 60.0 + seconds / 3600.0) * 15.0


def parse_sexagesimal_dec(dec_str: str) -> float:
    """Convert Dec from ±DD:MM:SS.S to decimal degrees."""
    sign = 1.0
    if dec_str.startswith('-'):
        sign = -1.0
        dec_str = dec_str[1:]
    elif dec_str.startswith('+'):
        dec_str = dec_str[1:]
    
    parts = dec_str.split(':')
    if len(parts) != 3:
        return float('nan')
    degrees = float(parts[0])
    minutes = float(parts[1])
    seconds = float(parts[2])
    return sign * (degrees + minutes / 60.0 + seconds / 3600.0)


def angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Calculate angular separation in arcminutes using the haversine formula."""
    ra1_rad = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad = math.radians(ra2)
    dec2_rad = math.radians(dec2)
    
    delta_ra = ra2_rad - ra1_rad
    delta_dec = dec2_rad - dec1_rad
    
    a = math.sin(delta_dec / 2)**2 + math.cos(dec1_rad) * math.cos(dec2_rad) * math.sin(delta_ra / 2)**2
    c = 2 * math.asin(math.sqrt(a))
    
    # Return in arcminutes
    return math.degrees(c) * 60.0


def match_comparison_stars(
    rows: list[dict[str, str]],
    aavso_stars: list[dict],
    target_ra: float,
    target_dec: float,
    match_tolerance_arcmin: float = 2.0
) -> dict[str, dict]:
    """
    Match comparison stars from the photometry table to AAVSO catalog.
    
    Returns a dict mapping star labels (C3, C4, etc.) to AAVSO star data.
    """
    # Extract comparison star positions from first row of photometry
    if not rows:
        return {}
    
    first_row = rows[0]
    comp_positions = {}
    
    # Find all comparison stars
    for key in first_row.keys():
        if key.startswith("RA_C") and key[4:].isdigit():
            label = key[3:]  # Extract "C3", "C4", etc.
            try:
                ra_val = float(first_row[key])
                dec_val = float(first_row[f"DEC_{label}"])
                # Convert RA hours to degrees if needed
                ra_deg = ra_val * 15.0 if ra_val <= 24.0 else ra_val
                comp_positions[label] = (ra_deg, dec_val)
            except (KeyError, ValueError, TypeError):
                continue
    
    # Also get T2 position for matching
    try:
        t2_ra = float(first_row["RA_T2"])
        t2_dec = float(first_row["DEC_T2"])
        t2_ra_deg = t2_ra * 15.0 if t2_ra <= 24.0 else t2_ra
        comp_positions["T2"] = (t2_ra_deg, t2_dec)
    except (KeyError, ValueError, TypeError):
        pass
    
    # Match each comparison star to AAVSO catalog
    matches = {}
    
    for label, (comp_ra, comp_dec) in comp_positions.items():
        best_match = None
        best_distance = float('inf')
        
        for aavso_star in aavso_stars:
            aavso_ra_str = aavso_star.get("ra", "")
            aavso_dec_str = aavso_star.get("dec", "")
            
            aavso_ra = parse_sexagesimal_ra(aavso_ra_str)
            aavso_dec = parse_sexagesimal_dec(aavso_dec_str)
            
            if math.isnan(aavso_ra) or math.isnan(aavso_dec):
                continue
            
            distance = angular_separation(comp_ra, comp_dec, aavso_ra, aavso_dec)
            
            if distance < match_tolerance_arcmin and distance < best_distance:
                best_distance = distance
                best_match = {
                    **aavso_star,
                    "separation_arcmin": distance
                }
        
        if best_match:
            matches[label] = best_match
    
    return matches


def compute_transformation(
    rows: list[dict[str, str]],
    matches: dict[str, dict],
    color_basis: str = "vi"
) -> dict:
    """
    Compute zero-point transformation from instrumental to standard V magnitudes.
    
    Uses matched comparison stars to determine the transformation:
    V_standard = V_instrumental + ZP + k * color
    
    Returns dict with zero_point, color_coefficient, rms, num_stars used.
    """
    # Collect data points from all observations
    data_points = []
    
    for row in rows:
        for label, aavso_data in matches.items():
            # Get instrumental magnitude
            inst_mag_col = f"Source_AMag_{label}"
            if inst_mag_col not in row:
                continue
            
            try:
                inst_mag = float(row[inst_mag_col])
            except (ValueError, TypeError):
                continue
            
            # Get standard V magnitude from AAVSO
            bands = aavso_data.get("bands", [])
            std_v_mag = None
            std_b_mag = None
            std_i_mag = None
            
            for band in bands:
                if isinstance(band, dict):
                    band_name = band.get("band")
                    mag_val = band.get("mag")
                    if band_name == "V" and mag_val is not None:
                        std_v_mag = float(mag_val)
                    elif band_name == "B" and mag_val is not None:
                        std_b_mag = float(mag_val)
                    elif band_name == "Ic" and mag_val is not None:
                        std_i_mag = float(mag_val)
            
            if std_v_mag is None:
                continue
            
            # Calculate color index
            color_index = None
            if color_basis == "bv" and std_b_mag is not None:
                color_index = std_b_mag - std_v_mag
            elif color_basis == "vi" and std_i_mag is not None:
                color_index = std_v_mag - std_i_mag
            
            if color_index is not None:
                # Delta mag = standard - instrumental
                delta_mag = std_v_mag - inst_mag
                data_points.append({
                    "label": label,
                    "inst_mag": inst_mag,
                    "std_mag": std_v_mag,
                    "color": color_index,
                    "delta": delta_mag
                })
    
    if len(data_points) < 2:
        return {
            "zero_point": float('nan'),
            "color_coefficient": float('nan'),
            "rms": float('nan'),
            "num_stars": 0,
            "num_obs": 0
        }
    
    # Fit a line: delta_mag = ZP + k * color
    # Using simple least squares
    n = len(data_points)
    sum_color = sum(p["color"] for p in data_points)
    sum_delta = sum(p["delta"] for p in data_points)
    sum_color_sq = sum(p["color"]**2 for p in data_points)
    sum_color_delta = sum(p["color"] * p["delta"] for p in data_points)
    
    # Solve for slope (k) and intercept (ZP)
    denom = n * sum_color_sq - sum_color**2
    if abs(denom) < 1e-10:
        # All colors are the same, just compute mean zero point
        zero_point = sum_delta / n
        color_coeff = 0.0
    else:
        color_coeff = (n * sum_color_delta - sum_color * sum_delta) / denom
        zero_point = (sum_delta - color_coeff * sum_color) / n
    
    # Calculate RMS residual
    residuals = []
    for p in data_points:
        predicted = zero_point + color_coeff * p["color"]
        residuals.append(p["delta"] - predicted)
    
    rms = math.sqrt(sum(r**2 for r in residuals) / len(residuals)) if residuals else float('nan')
    
    # Count unique stars
    unique_stars = len(set(p["label"] for p in data_points))
    
    return {
        "zero_point": zero_point,
        "color_coefficient": color_coeff,
        "rms": rms,
        "num_stars": unique_stars,
        "num_obs": len(data_points),
        "data_points": data_points
    }


def apply_transformation_to_targets(
    rows: list[dict[str, str]],
    transformation: dict,
    aavso_stars: list[dict],
    matches: dict[str, dict],
    color_basis: str = "vi"
) -> dict[str, list[tuple[int, float]]]:
    """
    Apply the computed transformation to T1 and T2 to get standard V magnitudes.
    
    Returns dict with "T1" and "T2" keys containing lists of (row_index, V_standard) tuples.
    
    For T1 color: use median color of matched comparison stars.
    For T2 color: attempt to find T2 in AAVSO catalog or use median comparison color.
    """
    transformed = {"T1": [], "T2": []}
    
    zp = transformation["zero_point"]
    k = transformation["color_coefficient"]
    
    if math.isnan(zp):
        return transformed
    
    # Compute median color from matched comparison stars
    comp_colors = []
    for label, match_data in matches.items():
        bands = match_data.get("bands", [])
        band_dict = {}
        for band in bands:
            if isinstance(band, dict):
                band_name = band.get("band")
                mag_val = band.get("mag")
                if band_name and mag_val is not None:
                    band_dict[band_name] = float(mag_val)
        
        if color_basis == "bv" and "B" in band_dict and "V" in band_dict:
            comp_colors.append(band_dict["B"] - band_dict["V"])
        elif color_basis == "vr" and "V" in band_dict and "R" in band_dict:
            comp_colors.append(band_dict["V"] - band_dict["R"])
        elif color_basis == "vi" and "V" in band_dict and "Ic" in band_dict:
            comp_colors.append(band_dict["V"] - band_dict["Ic"])
    
    if not comp_colors:
        sys.stderr.write("warning: no comparison star colors available, using color=0.0 for targets\n")
        median_color = 0.0
    else:
        median_color = statistics.median(comp_colors)
    
    # Use median comparison color for T1
    t1_color = median_color
    
    # Try to find T2 in AAVSO catalog for its color
    # For now, also use median color as we'd need T2 position to match
    t2_color = median_color
    
    # Process each observation
    for row_index, row in enumerate(rows):
        # Transform T1
        try:
            t1_inst = float(row["Source_AMag_T1"])
            t1_std = t1_inst + zp + k * t1_color
            transformed["T1"].append((row_index, t1_std))
        except (KeyError, ValueError, TypeError):
            pass
        
        # Transform T2
        try:
            t2_inst = float(row["Source_AMag_T2"])
            t2_std = t2_inst + zp + k * t2_color
            transformed["T2"].append((row_index, t2_std))
        except (KeyError, ValueError, TypeError):
            pass
    
    return transformed


def load_coefficients(path: Path) -> Mapping[str, float]:
    parser = configparser.ConfigParser()
    with path.open(encoding="utf-8") as handle:
        parser.read_file(handle)

    if "Coefficients" not in parser:
        sys.stderr.write(f"error: section [Coefficients] not found in {path}\n")
        sys.exit(1)

    coeffs = {}
    for key, value in parser["Coefficients"].items():
        try:
            coeffs[key] = float(value)
        except ValueError as exc:  # pragma: no cover - defensive
            sys.stderr.write(f"error: could not parse coefficient {key}={value!r}: {exc}\n")
            sys.exit(1)
    return coeffs


def load_colors(path: Path, color_basis: str) -> Mapping[str, float]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    colors: dict[str, float] = {}
    for label, value in payload.items():
        if isinstance(value, Mapping):
            if color_basis not in value:
                sys.stderr.write(
                    f"error: no '{color_basis}' entry for star {label!r} in {path}\n"
                )
                sys.exit(1)
            colors[label.upper()] = float(value[color_basis])
        else:
            colors[label.upper()] = float(value)
    return colors


def apply_v_transformation(
    rows: list[dict[str, str]],
    star_labels: Iterable[str],
    coeffs: Mapping[str, float],
    colors: Mapping[str, float],
    color_basis: str,
) -> dict[str, list[float]]:
    coeff_key = COLOR_KEY_MAP[color_basis]
    if coeff_key not in coeffs:
        sys.stderr.write(
            f"error: coefficient '{coeff_key}' not present in coefficients file\n"
        )
        sys.exit(1)

    slope = coeffs[coeff_key]
    transformed: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        for label in star_labels:
            color_value = colors.get(label.upper())
            column = f"Source_AMag_{label}"
            if color_value is None or column not in row:
                continue
            try:
                mag = float(row[column])
            except (TypeError, ValueError):  # skip non-numeric entries
                continue
            transformed[label].append(mag + slope * color_value)
    return transformed


def emit_transformation_summary(
    transformed: Mapping[str, list[float]],
    color_basis: str,
) -> None:
    if not transformed:
        print("No transformed magnitudes were produced (check color mapping).")
        return

    print(f"Transformed V magnitudes (basis={color_basis.upper()}):")
    for label in sorted(transformed):
        values = transformed[label]
        if not values:
            continue
        mean_mag = statistics.fmean(values)
        print(f"  {label}: n={len(values)} mean={mean_mag:.4f} min={min(values):.4f} max={max(values):.4f}")


def write_augmented_table(
    output_path: Path,
    header: list[str],
    rows: list[dict[str, str]],
    transformed: Mapping[str, list[float]],
) -> None:
    fieldnames = list(header)
    augment_columns = []
    for label in sorted(transformed):
        column = f"V_Transformed_{label}"
        augment_columns.append(column)
        if column not in fieldnames:
            fieldnames.append(column)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row_index, row in enumerate(rows):
            row_out = dict(row)
            for label in sorted(transformed):
                values = transformed[label]
                column = f"V_Transformed_{label}"
                try:
                    value = values[row_index]
                except IndexError:
                    value = ""
                else:
                    value = f"{value:.6f}"
                row_out[column] = value
            writer.writerow(row_out)

    added = ", ".join(augment_columns)
    print(f"Wrote transformed data with columns [{added}] to {output_path}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("photometry", help="path to the photometry .tbl file")
    parser.add_argument(
        "--coeffs",
        help="path to transformation coefficients (.txt/.ini)" ,
        type=Path,
    )
    parser.add_argument(
        "--color-basis",
        choices=sorted(COLOR_KEY_MAP.keys()),
        help="color index to use with the V-band transformation (default: vi if available)",
    )
    parser.add_argument(
        "--colors",
        type=Path,
        help=(
            "JSON file mapping star IDs (e.g. 'T1', 'C3') to color indices. "
            "Values can either be numeric or nested dictionaries keyed by the selected color basis."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional path to write an augmented photometry table with transformed V magnitudes",
    )
    parser.add_argument(
        "--fetch-aavso",
        action="store_true",
        help="query the AAVSO VSP service using the target position",
    )
    parser.add_argument(
        "--aavso-radius",
        type=float,
        default=60.0,
        help="search radius in arcminutes for AAVSO comparison stars (default: 60)",
    )
    parser.add_argument(
        "--aavso-mag-limit",
        type=float,
        default=15.0,
        help="magnitude limit for AAVSO comparison stars (default: 15.0)",
    )
    parser.add_argument(
        "--match-tolerance",
        type=float,
        default=2.0/60.0,
        help="maximum separation in arcminutes for star matching (default: 2 arcsec = 0.033 arcmin)",
    )
    parser.add_argument(
        "--transform",
        action="store_true",
        help="compute transformation from instrumental to standard V magnitudes",
    )
    return parser


def main(argv: list[str]) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv[1:])

    table_path = Path(args.photometry).expanduser().resolve()
    header, rows = read_table(table_path)
    report_star_presence(header)

    if args.fetch_aavso:
        try:
            ra_deg, dec_deg = extract_target_position(rows)
        except RuntimeError as exc:
            sys.stderr.write(f"error: {exc}\n")
            sys.exit(1)
        print(
            f"Target coordinates (T1): RA={ra_deg:.6f}° ({format_ra_deg(ra_deg)}), "
            f"Dec={dec_deg:.6f}° ({format_dec_deg(dec_deg)})"
        )
        payload = query_aavso_comparisons(
            ra_deg,
            dec_deg,
            radius_arcmin=args.aavso_radius,
            mag_limit=args.aavso_mag_limit,
        )
        emit_aavso_summary(payload)
        
        # Perform star matching and transformation if requested
        if args.transform:
            aavso_stars = payload.get("photometry", [])
            if not aavso_stars:
                sys.stderr.write("error: no AAVSO stars available for transformation\n")
                sys.exit(1)
            
            print("\n=== Star Matching ===")
            matches = match_comparison_stars(
                rows, 
                aavso_stars, 
                ra_deg, 
                dec_deg,
                match_tolerance_arcmin=args.match_tolerance
            )
            
            if not matches:
                sys.stderr.write("error: no comparison stars matched between table and AAVSO\n")
                sys.exit(1)
            
            print(f"Matched {len(matches)} comparison stars:")
            for label, match_data in sorted(matches.items()):
                aavso_label = match_data.get("label", "?")
                aavso_auid = match_data.get("auid", "?")
                separation = match_data.get("separation_arcmin", 0)
                
                # Get V mag
                bands = match_data.get("bands", [])
                v_mag = None
                for band in bands:
                    if isinstance(band, dict) and band.get("band") == "V":
                        v_mag = band.get("mag")
                        break
                v_str = f"{v_mag:.3f}" if v_mag is not None else "  -  "
                
                print(f"  {label:>3s} -> AAVSO {aavso_label:>4s} ({aavso_auid:15s})  V={v_str}  sep={separation:.2f}\"")
            
            # Determine color basis
            color_basis = args.color_basis or "vi"
            
            print(f"\n=== Transformation (color basis: {color_basis.upper()}) ===")
            transformation = compute_transformation(rows, matches, color_basis=color_basis)
            
            if transformation["num_stars"] == 0:
                sys.stderr.write("error: insufficient data for transformation\n")
                sys.exit(1)
            
            zp = transformation["zero_point"]
            k = transformation["color_coefficient"]
            rms = transformation["rms"]
            n_stars = transformation["num_stars"]
            n_obs = transformation["num_obs"]
            
            print(f"Zero point (ZP):        {zp:+.4f} mag")
            print(f"Color coefficient (k):  {k:+.4f} mag/{color_basis.upper()}")
            print(f"RMS residual:           {rms:.4f} mag")
            print(f"Number of comp stars:   {n_stars}")
            print(f"Number of observations: {n_obs}")
            
            # Apply to targets
            print("\n=== Transformed Target Magnitudes ===")
            target_mags = apply_transformation_to_targets(
                rows, 
                transformation, 
                aavso_stars,
                matches,
                color_basis=color_basis
            )
            
            # Compute median color from comparison stars for display
            comp_colors = []
            for label, match_data in matches.items():
                bands = match_data.get("bands", [])
                band_dict = {}
                for band in bands:
                    if isinstance(band, dict):
                        band_name = band.get("band")
                        mag_val = band.get("mag")
                        if band_name and mag_val is not None:
                            band_dict[band_name] = float(mag_val)
                
                if color_basis == "bv" and "B" in band_dict and "V" in band_dict:
                    comp_colors.append(band_dict["B"] - band_dict["V"])
                elif color_basis == "vr" and "V" in band_dict and "R" in band_dict:
                    comp_colors.append(band_dict["V"] - band_dict["R"])
                elif color_basis == "vi" and "V" in band_dict and "Ic" in band_dict:
                    comp_colors.append(band_dict["V"] - band_dict["Ic"])
            
            median_color = statistics.median(comp_colors) if comp_colors else 0.0
            
            # Display summary statistics
            for target in ["T1", "T2"]:
                if target_mags[target]:
                    mags = [mag for _, mag in target_mags[target]]
                    mean_mag = statistics.fmean(mags)
                    std_dev = statistics.stdev(mags) if len(mags) > 1 else 0.0
                    print(f"\n{target} Summary:")
                    print(f"  Mean V:     {mean_mag:.4f} mag")
                    print(f"  Std Dev:    {std_dev:.4f} mag")
                    print(f"  N obs:      {len(mags)}")
                    print(f"  Color used: {median_color:+.3f} {color_basis.upper()} (median of comp stars)")
                    
                    # Show individual measurements
                    print(f"\n{target} Individual Measurements:")
                    print(f"  {'Obs#':>5s}  {'V_inst':>8s}  {'V_std':>8s}")
                    for row_idx, v_std in target_mags[target]:
                        try:
                            v_inst = float(rows[row_idx][f"Source_AMag_{target}"])
                            print(f"  {row_idx+1:5d}  {v_inst:8.4f}  {v_std:8.4f}")
                        except (KeyError, ValueError, IndexError):
                            pass

    # If no coefficients were provided, we're done after reporting presence.
    if not args.coeffs:
        return

    color_basis = args.color_basis or ("vi" if "vi" in COLOR_KEY_MAP else None)
    if not color_basis:
        parser.error("unable to determine which color basis to use; please supply --color-basis")

    if not args.colors:
        parser.error("--colors is required when applying transformations")

    coeffs = load_coefficients(args.coeffs.expanduser().resolve())
    colors = load_colors(args.colors.expanduser().resolve(), color_basis)

    star_labels = collect_star_labels(header)
    transformed = apply_v_transformation(rows, star_labels, coeffs, colors, color_basis)
    emit_transformation_summary(transformed, color_basis)

    if args.output:
        write_augmented_table(args.output.expanduser().resolve(), header, rows, transformed)


if __name__ == "__main__":
    main(sys.argv)
