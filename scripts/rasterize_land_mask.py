#!/usr/bin/env python3
"""Rasterize Natural Earth land polygons onto a lon/lat grid."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from shapely.geometry import box

LAND_GEOJSON_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "master/geojson/ne_10m_land.geojson"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rasterize land mask onto a lon/lat grid.")
    parser.add_argument("lon_min", type=float)
    parser.add_argument("lat_min", type=float)
    parser.add_argument("lon_max", type=float)
    parser.add_argument("lat_max", type=float)
    parser.add_argument("--n-lon", type=int, required=True)
    parser.add_argument("--n-lat", type=int, required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--cache-dir",
        default="data/bathymetry/bluetopo_cache",
        help="Directory for cached Natural Earth land polygons.",
    )
    return parser.parse_args()


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


def main() -> None:
    args = parse_args()
    lon_min, lon_max = sorted((args.lon_min, args.lon_max))
    lat_min, lat_max = sorted((args.lat_min, args.lat_max))
    cache_dir = Path(args.cache_dir)
    output_path = Path(args.output)

    land_gdf = load_land_polygons(cache_dir).to_crs("EPSG:4326")
    aoi = box(lon_min, lat_min, lon_max, lat_max)
    hits = land_gdf[land_gdf.intersects(aoi)]
    transform = from_bounds(lon_min, lat_min, lon_max, lat_max, args.n_lon, args.n_lat)
    land_raster = rasterize(
        [(geom, 1) for geom in hits.geometry],
        out_shape=(args.n_lat, args.n_lon),
        transform=transform,
        fill=0,
        all_touched=True,
        dtype=np.uint8,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "height": args.n_lat,
        "width": args.n_lon,
        "count": 1,
        "dtype": "uint8",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": 0,
    }
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(land_raster, 1)

    print(f"Wrote land mask with {int(land_raster.sum())} land cell(s) to {output_path}")


if __name__ == "__main__":
    main()
