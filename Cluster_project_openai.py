#!/usr/bin/env python3
"""Identify open-cluster members and estimate their astrometric properties.

The script queries Gaia DR3, identifies members with DBSCAN in proper motion
and parallax, applies a minimum membership score, and writes catalogues and
diagnostic diagrams.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


LOG = logging.getLogger("gaia_clusters")
GAIA_TABLE = "gaiadr3.gaia_source_lite"
QUERY_VERSION = "2"
MAX_G_MAG = 18.0
PARALLAX_OFFSET_MAS = 0.017
REQUIRED_COLUMNS = {
    "source_id",
    "ra",
    "dec",
    "parallax",
    "parallax_error",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
}
OPTIONAL_COLUMNS = {
    "ruwe",
    "parallax_over_error",
    "visibility_periods_used",
    "phot_bp_rp_excess_factor",
    "phot_g_mean_flux",
    "phot_g_mean_flux_error",
}


@dataclass(frozen=True)
class ClusterConfig:
    """Fixed starting values; all scientific cuts remain CLI-configurable."""

    name: str
    ra_deg: float
    dec_deg: float
    radius_deg: float
    expected_pmra: float
    expected_pmdec: float
    expected_parallax: float
    av_prior: float
    feh_prior: float


CLUSTERS: dict[str, ClusterConfig] = {
    "pleiades": ClusterConfig(
        "Pleiades", 56.75, 24.1167, 2.0, 19.15, -45.72, 7.35, 0.12, 0.00
    ),
    "praesepe": ClusterConfig(
        "Praesepe", 130.05, 19.67, 2.0, -35.70, -12.80, 5.40, 0.08, 0.12
    ),
    "m67": ClusterConfig(
        "M67", 132.825, 11.800, 2.0, -10.90, -2.90, 1.14, 0.13, 0.03
    ),
}


def require_runtime() -> dict[str, Any]:
    """Load optional scientific dependencies only when executing the pipeline."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from astroquery.gaia import Gaia
        from astropy.table import Table
        from sklearn.cluster import DBSCAN
        from sklearn.preprocessing import RobustScaler
    except ImportError as error:
        raise RuntimeError(
            "Missing a required package. Install: "
            "pip install astroquery astropy numpy pandas scipy scikit-learn matplotlib"
        ) from error
    return {
        "pd": pd,
        "plt": plt,
        "Gaia": Gaia,
        "Table": Table,
        "DBSCAN": DBSCAN,
        "RobustScaler": RobustScaler,
    }


def cache_paths(cache_dir: Path, cluster: ClusterConfig) -> tuple[Path, Path]:
    stem = cluster.name.lower().replace(" ", "_")
    return cache_dir / f"{stem}_raw.ecsv", cache_dir / f"{stem}_query.json"


def fetch_catalogue(
    cluster: ClusterConfig,
    cache_dir: Path,
    refresh: bool,
    offline: bool,
    runtime: dict[str, Any],
) -> tuple[Any, dict[str, Any]]:
    """Return cached data or query the Gaia DR3 source-lite table."""
    table_path, metadata_path = cache_paths(cache_dir, cluster)
    table_class = runtime["Table"]
    if table_path.exists() and metadata_path.exists() and not refresh:
        metadata = json.loads(metadata_path.read_text())
        query = metadata.get("query", {})
        cache_matches = (
            metadata.get("catalog") == GAIA_TABLE
            and metadata.get("query_version") == QUERY_VERSION
            and query.get("radius_deg") == cluster.radius_deg
            and query.get("max_g_mag") == MAX_G_MAG
        )
        if cache_matches:
            LOG.info("Using cached catalogue for %s", cluster.name)
            return table_class.read(table_path, format="ascii.ecsv"), metadata
        LOG.info("Ignoring stale catalogue cache for %s", cluster.name)
    if offline:
        raise FileNotFoundError(
            f"No current Gaia cache for {cluster.name}. Re-run without --offline to query Gaia."
        )

    columns = [
        "source_id", "ra", "dec", "parallax", "parallax_error", "pmra",
        "pmra_error", "pmdec", "pmdec_error", "phot_g_mean_mag",
        "phot_bp_mean_mag", "phot_rp_mean_mag", "ruwe", "parallax_over_error",
        "phot_bp_rp_excess_factor",
    ]
    adql = f"""
        SELECT {", ".join(columns)}
        FROM {GAIA_TABLE}
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {cluster.ra_deg}, {cluster.dec_deg}, {cluster.radius_deg})
        ) = 1
        AND phot_g_mean_mag <= {MAX_G_MAG}
    """
    LOG.info(
        "Querying %s (%s degree cone, G <= %s) via Gaia %s",
        cluster.name,
        cluster.radius_deg,
        MAX_G_MAG,
        GAIA_TABLE,
    )
    try:
        job = runtime["Gaia"].launch_job_async(adql, dump_to_file=False)
        table = job.get_results()
    except Exception as error:
        raise RuntimeError(f"Gaia query failed for {cluster.name}: {error}") from error

    metadata = {
        "cluster": asdict(cluster),
        "source": "Gaia Archive",
        "catalog": GAIA_TABLE,
        "query": {
            "ra_deg": cluster.ra_deg,
            "dec_deg": cluster.dec_deg,
            "radius_deg": cluster.radius_deg,
            "max_g_mag": MAX_G_MAG,
        },
        "query_version": QUERY_VERSION,
        "selected_columns": columns,
        "adql": " ".join(adql.split()),
        "row_count": len(table),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
    }

    cache_dir.mkdir(parents=True, exist_ok=True)
    table.write(table_path, format="ascii.ecsv", overwrite=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return table, metadata


def to_dataframe(table: Any, runtime: dict[str, Any]) -> Any:
    dataframe = table.to_pandas()
    for column in dataframe.columns:
        dataframe[column] = dataframe[column].replace([np.inf, -np.inf], np.nan)
    return dataframe


def prepare_catalogue(dataframe: Any) -> Any:
    """Add quality flags and derived astrometry without discarding the raw sample."""
    required_finite = list(REQUIRED_COLUMNS - {"source_id"})
    valid = dataframe[required_finite].notna().all(axis=1)
    valid &= dataframe["parallax_error"] > 0
    valid &= dataframe["pmra_error"] > 0
    valid &= dataframe["pmdec_error"] > 0
    if "parallax_over_error" in dataframe:
        valid &= dataframe["parallax_over_error"] >= 5.0
    else:
        valid &= dataframe["parallax"] / dataframe["parallax_error"] >= 5.0
    if "ruwe" in dataframe:
        valid &= dataframe["ruwe"].fillna(np.inf) < 1.4
    if "visibility_periods_used" in dataframe:
        valid &= dataframe["visibility_periods_used"].fillna(0) >= 8
    if "phot_bp_rp_excess_factor" in dataframe:
        colour = dataframe["phot_bp_mean_mag"] - dataframe["phot_rp_mean_mag"]
        lower = 1.0 + 0.015 * colour**2
        upper = 1.3 + 0.06 * colour**2
        excess = dataframe["phot_bp_rp_excess_factor"]
        valid &= excess.isna() | ((excess > lower) & (excess < upper))

    dataframe = dataframe.copy()
    dataframe["quality_flag"] = valid
    dataframe["bp_rp"] = dataframe["phot_bp_mean_mag"] - dataframe["phot_rp_mean_mag"]
    dataframe["parallax_corrected"] = dataframe["parallax"] + PARALLAX_OFFSET_MAS
    return dataframe


def select_members(
    dataframe: Any,
    cluster: ClusterConfig,
    eps: float,
    min_samples: int,
    runtime: dict[str, Any],
) -> tuple[Any, dict[str, float]]:
    """Find the astrometric component closest to the configured cluster locus."""
    result = dataframe.copy()
    result["cluster_label"] = -1
    result["member_score"] = 0.0
    result["membership_probability"] = 0.0
    result["member_flag"] = False
    usable = result["quality_flag"].to_numpy()
    if usable.sum() < min_samples:
        raise RuntimeError(f"Only {usable.sum()} quality sources remain; cannot cluster {cluster.name}.")

    features = result.loc[usable, ["pmra", "pmdec", "parallax_corrected"]].to_numpy()
    scaler = runtime["RobustScaler"]()
    scaled = scaler.fit_transform(features)
    labels = runtime["DBSCAN"](eps=eps, min_samples=min_samples).fit_predict(scaled)
    result.loc[usable, "cluster_label"] = labels
    expected_scaled = scaler.transform(
        [[cluster.expected_pmra, cluster.expected_pmdec, cluster.expected_parallax]]
    )[0]
    candidates = [label for label in np.unique(labels) if label >= 0]
    if not candidates:
        raise RuntimeError(
            f"No DBSCAN component found for {cluster.name}. Try increasing --dbscan-eps."
        )

    def label_cost(label: int) -> float:
        component = scaled[labels == label]
        centre = np.median(component, axis=0)
        expected_distance = float(np.linalg.norm(centre - expected_scaled))
        sky = result.loc[usable].iloc[np.flatnonzero(labels == label)][["ra", "dec"]].to_numpy()
        sky_spread = float(np.nanmedian(np.hypot(sky[:, 0] - cluster.ra_deg, sky[:, 1] - cluster.dec_deg)))
        return expected_distance + 0.10 * sky_spread - 0.001 * len(component)

    chosen_label = min(candidates, key=label_cost)
    member_positions = np.flatnonzero(labels == chosen_label)
    component = scaled[member_positions]
    centre = np.median(component, axis=0)
    scatter = np.maximum(1.4826 * np.median(np.abs(component - centre), axis=0), 0.10)
    distance = np.sqrt(np.sum(((scaled - centre) / scatter) ** 2, axis=1))
    scores = np.exp(-0.5 * distance**2)
    usable_indices = result.index[usable]
    result.loc[usable_indices, "member_score"] = scores
    result.loc[usable_indices, "membership_probability"] = scores
    result.loc[usable_indices, "member_flag"] = (labels == chosen_label) & (scores >= 0.5)
    parameters = {
        "label": float(chosen_label),
        "pmra": float(np.median(result.loc[result.member_flag, "pmra"])),
        "pmdec": float(np.median(result.loc[result.member_flag, "pmdec"])),
        "parallax": float(np.median(result.loc[result.member_flag, "parallax_corrected"])),
        "n_members": float(result.member_flag.sum()),
    }
    return result, parameters


def estimate_distance(members: Any) -> tuple[float, float, float]:
    parallaxes = members["parallax_corrected"].to_numpy()
    errors = members["parallax_error"].to_numpy()
    weights = 1.0 / np.maximum(errors, 0.02) ** 2
    parallax = float(np.average(parallaxes, weights=weights))
    parallax_error = float(max(1.0 / np.sqrt(weights.sum()), 0.015))
    if parallax <= 0:
        raise RuntimeError("Member parallax is non-positive; no reliable distance modulus is available.")
    distance_pc = 1000.0 / parallax
    distance_modulus = 5.0 * np.log10(distance_pc) - 5.0
    return parallax, parallax_error, float(distance_modulus)


def add_hr_quantities(dataframe: Any, distance_modulus: float, av: float) -> Any:
    result = dataframe.copy()
    # Broad-band Gaia extinction coefficients used as a transparent first-pass model.
    result["abs_g"] = result["phot_g_mean_mag"] - distance_modulus - 0.86 * av
    result["bp_rp_0"] = result["bp_rp"] - 0.42 * av
    return result


def plot_diagnostics(dataframe: Any, cluster: ClusterConfig, output_dir: Path, runtime: dict[str, Any]) -> None:
    plt = runtime["plt"]
    output_dir.mkdir(parents=True, exist_ok=True)
    quality = dataframe[dataframe.quality_flag]
    members = dataframe[dataframe.member_flag]
    figure, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    axes[0, 0].scatter(dataframe.ra, dataframe.dec, s=1, c="0.75", label="all")
    axes[0, 0].scatter(members.ra, members.dec, s=4, c="tab:red", label="members")
    axes[0, 0].set(xlabel="RA (deg)", ylabel="Dec (deg)", title=f"{cluster.name}: sky position")
    axes[0, 0].legend(markerscale=3)
    axes[0, 1].scatter(quality.pmra, quality.pmdec, s=1, c="0.65")
    axes[0, 1].scatter(members.pmra, members.pmdec, s=4, c="tab:red")
    axes[0, 1].set(xlabel="pmRA (mas/yr)", ylabel="pmDec (mas/yr)", title="Proper motion")
    axes[1, 0].hist(quality.parallax_corrected, bins=80, color="0.65")
    axes[1, 0].hist(members.parallax_corrected, bins=40, color="tab:red", alpha=0.8)
    axes[1, 0].set(xlabel="Corrected parallax (mas)", ylabel="Sources", title="Parallax")
    axes[1, 1].scatter(quality.pmra, quality.parallax_corrected, s=1, c="0.65")
    axes[1, 1].scatter(members.pmra, members.parallax_corrected, s=4, c="tab:red")
    axes[1, 1].set(xlabel="pmRA (mas/yr)", ylabel="Corrected parallax (mas)", title="Astrometric locus")
    figure.savefig(output_dir / f"{cluster.name.lower()}_diagnostics.png", dpi=180)
    plt.close(figure)


def plot_hr_diagram(
    dataframe: Any,
    cluster: ClusterConfig,
    output_dir: Path,
    runtime: dict[str, Any],
) -> None:
    plt = runtime["plt"]
    quality = dataframe[dataframe.quality_flag]
    members = dataframe[dataframe.member_flag]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for axis in axes:
        axis.scatter(quality.bp_rp_0, quality.abs_g, s=1, c="0.80", rasterized=True)
        axis.invert_yaxis()
        axis.set(xlabel=r"$(G_{BP} - G_{RP})_0$", ylabel=r"$M_G$", title=cluster.name)
    points = axes[0].scatter(
        members.bp_rp_0, members.abs_g, c=members.membership_probability, s=5, cmap="viridis"
    )
    figure.colorbar(points, ax=axes[0], label="membership probability")
    axes[1].scatter(members.bp_rp_0, members.abs_g, s=5, c="tab:red", rasterized=True)
    figure.savefig(output_dir / f"{cluster.name.lower()}_hr_diagram.png", dpi=180)
    plt.close(figure)


def plot_color_diagram(
    dataframe: Any,
    cluster: ClusterConfig,
    output_dir: Path,
    runtime: dict[str, Any],
) -> None:
    """Plot apparent Gaia colour versus G magnitude."""
    plt = runtime["plt"]
    quality = dataframe[dataframe.quality_flag]
    members = dataframe[dataframe.member_flag]
    figure, axis = plt.subplots(figsize=(7, 6), constrained_layout=True)
    axis.scatter(quality.bp_rp, quality.phot_g_mean_mag, s=1, c="0.80", rasterized=True, label="quality sample")
    axis.scatter(
        members.bp_rp,
        members.phot_g_mean_mag,
        c=members.membership_probability,
        s=6,
        cmap="viridis",
        rasterized=True,
        label="members",
    )
    axis.invert_yaxis()
    axis.set(
        xlabel=r"$G_{BP} - G_{RP}$",
        ylabel=r"$G$ (mag)",
        title=f"{cluster.name}: Gaia colour–magnitude diagram",
    )
    axis.legend(markerscale=2)
    figure.savefig(output_dir / f"{cluster.name.lower()}_color_diagram.png", dpi=180)
    plt.close(figure)


def estimate_turnoff(members: Any) -> dict[str, float]:
    """Report a reproducible turn-off proxy for the selected members."""
    candidates = members[(members.bp_rp_0 > -0.2) & (members.bp_rp_0 < 1.3)]
    if len(candidates) < 5:
        candidates = members
    brightest = candidates.nsmallest(max(1, min(5, len(candidates))), "abs_g")
    return {
        "turnoff_abs_g_proxy": float(brightest.abs_g.median()),
        "turnoff_bp_rp_proxy": float(brightest.bp_rp_0.median()),
    }


def write_summary(rows: list[dict[str, Any]], output_dir: Path, runtime: dict[str, Any]) -> None:
    dataframe = runtime["pd"].DataFrame(rows)
    dataframe.to_csv(output_dir / "cluster_summary.csv", index=False)
    lines = ["# Gaia Cluster Analysis Summary", "", dataframe.to_markdown(index=False), ""]
    (output_dir / "report.md").write_text("\n".join(lines))


def process_cluster(cluster: ClusterConfig, arguments: argparse.Namespace, runtime: dict[str, Any]) -> dict[str, Any]:
    table, metadata = fetch_catalogue(
        cluster,
        arguments.cache_dir,
        arguments.refresh,
        arguments.offline,
        runtime,
    )
    dataframe = prepare_catalogue(to_dataframe(table, runtime))
    dataframe, membership = select_members(
        dataframe, cluster, arguments.dbscan_eps, arguments.min_samples, runtime
    )
    members = dataframe[dataframe.member_flag].copy()
    if len(members) < 5:
        raise RuntimeError(f"Only {len(members)} members selected for {cluster.name}.")
    parallax, parallax_error, distance_modulus = estimate_distance(members)
    dataframe = add_hr_quantities(dataframe, distance_modulus, cluster.av_prior)
    members = dataframe[dataframe.member_flag].copy()
    plot_diagnostics(dataframe, cluster, arguments.output_dir, runtime)
    plot_hr_diagram(dataframe, cluster, arguments.output_dir, runtime)
    plot_color_diagram(dataframe, cluster, arguments.output_dir, runtime)
    turnoff = estimate_turnoff(members)
    member_path = arguments.output_dir / f"{cluster.name.lower()}_members.csv"
    members.to_csv(member_path, index=False)
    (arguments.output_dir / f"{cluster.name.lower()}_metadata.json").write_text(
        json.dumps({"query": metadata, "membership": membership, "turnoff": turnoff}, indent=2)
        + "\n"
    )
    return {
        "cluster": cluster.name,
        "n_queried": len(dataframe),
        "n_quality": int(dataframe.quality_flag.sum()),
        "n_members": int(len(members)),
        "parallax_mas": parallax,
        "parallax_error_mas": parallax_error,
        "distance_modulus": distance_modulus,
        "av": cluster.av_prior,
        "feh": cluster.feh_prior,
        **turnoff,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", nargs="+", choices=sorted(CLUSTERS), default=sorted(CLUSTERS))
    parser.add_argument("--cache-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--refresh", action="store_true", help="Ignore cached Gaia results")
    parser.add_argument("--offline", action="store_true", help="Require cached Gaia results")
    parser.add_argument("--dbscan-eps", type=float, default=0.55)
    parser.add_argument("--min-samples", type=int, default=12)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    logging.basicConfig(level=logging.DEBUG if arguments.verbose else logging.INFO, format="%(levelname)s: %(message)s")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    runtime = require_runtime()
    summaries: list[dict[str, Any]] = []
    for name in arguments.clusters:
        summaries.append(process_cluster(CLUSTERS[name], arguments, runtime))
    write_summary(summaries, arguments.output_dir, runtime)
    LOG.info("Completed %d cluster(s). Results: %s", len(summaries), arguments.output_dir)


if __name__ == "__main__":
    main()