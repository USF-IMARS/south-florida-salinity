#!/usr/bin/env python3
"""Sample NOAA BlueTopo bathymetry onto a regular lon/lat grid (no full tile download)."""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import boto3
import geopandas as gpd
import numpy as np
import rasterio
from botocore import UNSIGNED
from botocore.config import Config
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject
from shapely.geometry import box

BUCKET = "noaa-ocs-nationalbathymetry-pds"
TILE_SCHEME_KEY = "BlueTopo/_BlueTopo_Tile_Scheme/BlueTopo_Tile_Scheme_20260626_132625.gpkg"
ELEVATION_BAND = 1
LAND_GEOJSON_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "master/geojson/ne_10m_land.geojson"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sample NOAA BlueTopo elevation onto a lon/lat grid."
    )
    parser.add_argument("lon_min", type=float)
    parser.add_argument("lat_min", type=float)
    parser.add_argument("lon_max", type=float)
    parser.add_argument("lat_max", type=float)
    parser.add_argument("--n-lon", type=int, required=True)
    parser.add_argument("--n-lat", type=int, required=True)
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output GeoTIFF (EPSG:4326, band 1 = BlueTopo elevation).",
    )
    parser.add_argument(
        "--cache-dir",
        default="data/bathymetry/bluetopo_cache",
        help="Directory for cached tile scheme.",
    )
    return parser.parse_args()


def s3_client():
    return boto3.client("s3", config=Config(signature_version=UNSIGNED))


def tile_scheme_path(cache_dir: Path) -> Path:
    return cache_dir / "BlueTopo_Tile_Scheme.gpkg"


def load_tile_scheme(cache_dir: Path) -> gpd.GeoDataFrame:
    scheme_path = tile_scheme_path(cache_dir)
    scheme_path.parent.mkdir(parents=True, exist_ok=True)
    if not scheme_path.exists():
        print(f"Downloading tile scheme from s3://{BUCKET}/{TILE_SCHEME_KEY}")
        s3_client().download_file(BUCKET, TILE_SCHEME_KEY, str(scheme_path))
    gdf = gpd.read_file(scheme_path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf


def tiles_for_bbox(
    gdf: gpd.GeoDataFrame,
    lon_min: float,
    lat_min: float,
    lon_max: float,
    lat_max: float,
    *,
    resolution: str | None = None,
) -> list[str]:
    aoi = box(lon_min, lat_min, lon_max, lat_max)
    hits = gdf[gdf.intersects(aoi)].copy()
    if hits.empty:
        raise RuntimeError("No BlueTopo tiles intersect the requested bounding box.")

    if resolution is not None:
        urls = (
            hits[hits["Resolution"] == resolution]["GeoTIFF_Link"]
            .dropna()
            .unique()
            .tolist()
        )
        return urls

    preferred = hits[hits["Resolution"] == "4m"]["GeoTIFF_Link"].dropna().unique().tolist()
    if preferred:
        return preferred

    fallback = hits["GeoTIFF_Link"].dropna().unique().tolist()
    if not fallback:
        raise RuntimeError("No BlueTopo GeoTIFF URLs found for this area.")
    return fallback


def mosaic_bluetopo(
    tile_urls: list[str],
    lon_min: float,
    lat_min: float,
    lon_max: float,
    lat_max: float,
    n_lon: int,
    n_lat: int,
    *,
    resampling: Resampling = Resampling.nearest,
) -> tuple[np.ndarray, rasterio.transform.Affine]:
    dst_crs = "EPSG:4326"
    dst_transform = from_bounds(lon_min, lat_min, lon_max, lat_max, n_lon, n_lat)
    dst_shape = (n_lat, n_lon)
    dst = np.full(dst_shape, np.nan, dtype="float32")

    for url in tile_urls:
        vsicurl = f"/vsicurl/{url}"
        temp = np.full(dst_shape, np.nan, dtype="float32")
        with rasterio.open(vsicurl) as src:
            reproject(
                source=rasterio.band(src, ELEVATION_BAND),
                destination=temp,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=resampling,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
        fill = np.isnan(dst) & ~np.isnan(temp)
        dst[fill] = temp[fill]
        overlap = ~np.isnan(dst) & ~np.isnan(temp)
        if overlap.any():
            dst[overlap] = temp[overlap]

    return dst, dst_transform


def fill_missing_elevation(
    elevation: np.ndarray,
    fill_urls: list[str],
    dst_transform: rasterio.transform.Affine,
    dst_crs: str,
    *,
    resampling: Resampling = Resampling.nearest,
) -> np.ndarray:
    if not fill_urls or not np.isnan(elevation).any():
        return elevation

    dst = elevation.copy()
    dst_shape = dst.shape
    missing = np.isnan(dst)
    for url in fill_urls:
        vsicurl = f"/vsicurl/{url}"
        temp = np.full(dst_shape, np.nan, dtype="float32")
        with rasterio.open(vsicurl) as src:
            reproject(
                source=rasterio.band(src, ELEVATION_BAND),
                destination=temp,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                resampling=resampling,
                src_nodata=np.nan,
                dst_nodata=np.nan,
            )
        fill = missing & ~np.isnan(temp)
        if not fill.any():
            continue
        dst[fill] = temp[fill]
        missing = np.isnan(dst)

    return dst


def land_geojson_path(cache_dir: Path) -> Path:
    return cache_dir / "ne_10m_land.geojson"


def load_land_polygons(cache_dir: Path) -> gpd.GeoDataFrame:
    path = land_geojson_path(cache_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        print(f"Downloading Natural Earth land polygons from {LAND_GEOJSON_URL}")
        urllib.request.urlretrieve(LAND_GEOJSON_URL, path)
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf


def apply_land_mask(
    elevation: np.ndarray,
    transform: rasterio.transform.Affine,
    cache_dir: Path,
    lon_min: float,
    lat_min: float,
    lon_max: float,
    lat_max: float,
) -> np.ndarray:
    land_gdf = load_land_polygons(cache_dir).to_crs("EPSG:4326")
    aoi = box(lon_min, lat_min, lon_max, lat_max)
    hits = land_gdf[land_gdf.intersects(aoi)]
    if hits.empty:
        return elevation, np.zeros(elevation.shape, dtype=np.uint8)

    land_raster = rasterize(
        [(geom, 1) for geom in hits.geometry],
        out_shape=elevation.shape,
        transform=transform,
        fill=0,
        all_touched=True,
        dtype=np.uint8,
    )
    masked = elevation.copy()
    on_land = land_raster == 1
    masked[on_land] = np.nan
    print(f"Natural Earth land mask cleared {int(on_land.sum())} cell(s)")
    return masked, land_raster


def main() -> None:
    args = parse_args()
    lon_min, lon_max = sorted((args.lon_min, args.lon_max))
    lat_min, lat_max = sorted((args.lat_min, args.lat_max))
    cache_dir = Path(args.cache_dir)
    output_path = Path(args.output)

    gdf = load_tile_scheme(cache_dir)
    urls_4m = tiles_for_bbox(gdf, lon_min, lat_min, lon_max, lat_max, resolution="4m")
    urls_8m = tiles_for_bbox(gdf, lon_min, lat_min, lon_max, lat_max, resolution="8m")
    if not urls_4m and not urls_8m:
        raise RuntimeError("No BlueTopo 4 m or 8 m tiles found for this area.")
    print(
        f"Sampling BlueTopo at {args.n_lon} x {args.n_lat}: "
        f"{len(urls_4m)} tile(s) at 4 m"
        + (f", {len(urls_8m)} tile(s) at 8 m for gap fill" if urls_8m else "")
    )

    elevation, transform = mosaic_bluetopo(
        urls_4m or urls_8m,
        lon_min,
        lat_min,
        lon_max,
        lat_max,
        args.n_lon,
        args.n_lat,
        resampling=Resampling.nearest,
    )
    if urls_4m and urls_8m:
        elevation = fill_missing_elevation(
            elevation,
            urls_8m,
            transform,
            "EPSG:4326",
            resampling=Resampling.nearest,
        )

    elevation, land_raster = apply_land_mask(
        elevation,
        transform,
        cache_dir,
        lon_min,
        lat_min,
        lon_max,
        lat_max,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "height": args.n_lat,
        "width": args.n_lon,
        "count": 2,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": np.nan,
    }
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(elevation, 1)
        dst.write(land_raster.astype("float32"), 2)
        dst.set_band_description(1, "BlueTopo elevation (NAVD88, m)")
        dst.set_band_description(2, "Natural Earth land mask (1=land)")

    if not np.isfinite(elevation).any():
        raise RuntimeError("BlueTopo sampling returned no valid elevation values.")

    valid = int(np.isfinite(elevation).sum())
    total = elevation.size
    land = int(np.sum(elevation >= 0))
    unknown = int(np.sum(~np.isfinite(elevation)))
    print(
        "Elevation range:",
        float(np.nanmin(elevation)),
        "to",
        float(np.nanmax(elevation)),
    )
    print(
        f"Coverage: {valid}/{total} cells ({100 * valid / total:.1f}%), "
        f"land (elev>=0): {land}, unknown/land-masked: {unknown}"
    )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
