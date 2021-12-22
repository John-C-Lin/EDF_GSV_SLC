#!/usr/bin/env python

import os
from typing import List, NamedTuple, Tuple

import fastkml as kml
import numpy as np
import numpy.typing as npt
import pandas as pd
from shapely.geometry import Point, Polygon
from shapely import vectorized


def median(x: pd.Series) -> float:
    return x.median(skipna=True)


def count(x: pd.Series) -> float:
    return (~x.isna()).sum()


# def floor(x: float, digits: int = 0) -> float:
#     """Return floor to specified number of digits."""
#     factor = 10 ** digits
#     return np.floor(x * factor) / factor


# def ceil(x: float, digits: int = 0) -> float:
#     """Return ceiling to specified number of digits."""
#     factor = 10 ** digits
#     return np.ceil(x * factor) / factor


# def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
#     return (df - df.min()) / (df.max() - df.min())


# def get_point_grid(
#     latitude_min: float, latitude_max: float, longitude_min: float, longitude_max: float
# ) -> List[Point]:
#     longitudes = np.arange(longitude_min, longitude_max, 1e-3)
#     latitudes = np.arange(latitude_min, latitude_max, 1e-3)

#     points = []
#     for longitude in longitudes:
#         for latitude in latitudes:
#             points.append(Point(longitude, latitude))

#     return points


def get_polygons(filename: str) -> List[kml.Placemark]:
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} missing")

    with open(filename, "rb") as f:
        xml = f.read()

    k = kml.KML()
    k.from_string(xml)

    feature = next(k.features())
    if type(feature) is not kml.Document:
        raise ValueError("Top level feature is not kml.Document")

    document = feature
    polygons_with_background_drives = list(document.features())

    names = set([polygon.name for polygon in polygons_with_background_drives])
    names = names - {"Emigration_Canyon", "Magna_Loop", "Historic"}

    polygons = []
    for polygon in polygons_with_background_drives:
        if polygon.name in names and type(polygon) is kml.Placemark:
            polygons.append(polygon)

    return polygons


# def get_polygon_by_location(
#     latitude: float, longitude: float, polygons: List[kml.Placemark]
# ) -> str:
#     for polygon in polygons:
#         geometry = polygon.geometry
#         if geometry.contains()


# def get_points_in_polygons(points: List[Point], polygons: List[Polygon]) -> List[dict]:
#     polygon_points = {polygon.name: [] for polygon in polygons}
#     for polygon in polygons:
#         geometry = polygon.geometry
#         for point in points:
#             if geometry.contains(point):
#                 longitude = np.round(point.x)
#                 latitude = np.round(point.y)
#                 polygon_points[polygon.name].append(
#                     {"longitude": longitude, "latitude": latitude}
#                 )
#     return polygon_points


# def rle(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
#     """Run length encoding to identify repeated values

#     Parameters
#     ----------
#     x : np.ndarray
#         vector to encode

#     Returns
#     -------
#     Tuple[start, length, value]
#         run length encoded values and positions
#     """
#     x = np.asarray(x)
#     n = len(x)
#     if n == 0:
#         return (
#             np.array([], dtype=int),
#             np.array([], dtype=int),
#             np.array([], dtype=x.dtype),
#         )

#     start = np.r_[0, np.flatnonzero(~np.isclose(x[1:], x[:-1])) + 1]
#     length = np.diff(np.r_[start, n])
#     return (start, length, x[start])


def haversine_distance(
    points: npt.ArrayLike,
    prev_points: npt.ArrayLike,
    radius_meters: float = 6371008.8,
) -> np.ndarray:
    """Return distance between each set of lonlat points."""
    points_radians = np.array(points) * np.pi / 180
    prev_points_radians = np.array(prev_points) * np.pi / 180

    dlongitude_radians = points_radians[:, 0] - prev_points_radians[:, 0]
    dlatitude_radians = points_radians[:, 1] - prev_points_radians[:, 1]

    a = (
        np.sin(dlatitude_radians / 2) ** 2
        + np.cos(points_radians[:, 0])
        * np.cos(prev_points_radians[:, 0])
        * np.sin(dlongitude_radians / 2) ** 2
    )
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return radius_meters * c


def get_speed(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike, time: npt.ArrayLike
) -> np.ndarray:
    """Calculate approximate speed in meters per second from location timeseries."""
    order = np.argsort(time)
    if any(np.diff(order) < 0):
        raise ValueError("records must be sorted in time")

    longitude = np.array(longitude)
    latitude = np.array(latitude)
    time = np.array(time)

    points = np.column_stack([longitude[1:], latitude[1:]])
    prev_points = np.column_stack([longitude[0:-1], latitude[0:-1]])

    distance_meters = haversine_distance(points, prev_points)
    interval_seconds = np.diff(time) / np.timedelta64(1, "s")

    speed_meters_per_second = distance_meters / interval_seconds
    invalid_mask = (
        (interval_seconds > 2)
        | (distance_meters > 80)
        | (speed_meters_per_second > 40)  # 90 mph
    )
    speed_meters_per_second[invalid_mask] = np.nan
    return np.insert(speed_meters_per_second, 0, np.nan)


observations = pd.read_feather("data/tables/4_primary_keys_uuids.feather")

observations = observations.groupby("system_id").apply(
    lambda group: group.assign(
        speed_meters_per_second=get_speed(
            group.longitude.values, group.latitude.values, group.time.values
        )
    )
)

observations.to_feather("data/receptors/observations.feather")


receptor_list = []
for index, group in observations.groupby("uuid"):
    uuid = index
    system_id = group.system_id.iloc[0]
    n = len(group)

    polygon_id = group.polygon.value_counts().index[0]

    receptor = {
        "time": group.time.median(),
        "system_id": system_id,
        "uuid": uuid,
        "n": n,
        "longitude": median(group.longitude),
        "latitude": median(group.latitude),
        "o3_ppb": median(group.o3_ppb),
        "o3_ppb_n": count(group.o3_ppb),
        "ch4d_ppm": median(group.ch4d_ppm),
        "ch4d_ppm_ex": median(group.ch4d_ppm_ex),
        "ch4d_ppm_base": median(group.ch4d_ppm_base),
        "ch4d_ppm_n": count(group.ch4d_ppm),
        "co2d_ppm": median(group.co2d_ppm),
        "co2d_ppm_ex": median(group.co2d_ppm_ex),
        "co2d_ppm_base": median(group.co2d_ppm_base),
        "co2d_ppm_n": count(group.co2d_ppm),
        "bc_ngm3": median(group.bc_ngm3),
        "bc_ngm3_ex": median(group.bc_ngm3_ex),
        "bc_ngm3_base": median(group.bc_ngm3_base),
        "bc_ngm3_n": count(group.bc_ngm3),
        "co_ppb": median(group.co_ppb),
        "co_ppb_ex": median(group.co_ppb_ex),
        "co_ppb_base": median(group.co_ppb_base),
        "co_ppb_n": count(group.co_ppb),
        "no_ppb": median(group.no_ppb),
        "no_ppb_ex": median(group.no_ppb_ex),
        "no_ppb_base": median(group.no_ppb_base),
        "no_ppb_n": count(group.no_ppb),
        "no2_ppb": median(group.no2_ppb),
        "no2_ppb_ex": median(group.no2_ppb_ex),
        "no2_ppb_base": median(group.no2_ppb_base),
        "no2_ppb_n": count(group.no2_ppb),
        "pm25_ugm3": median(group.pm25_ugm3),
        "pm25_ugm3_ex": median(group.pm25_ugm3_ex),
        "pm25_ugm3_base": median(group.pm25_ugm3_base),
        "pm25_ugm3_n": count(group.pm25_ugm3),
    }
    receptor_list.append(receptor)

receptors = pd.DataFrame(receptor_list)
# receptors = pd.read_feather("data/receptors/receptors.feather")

polygons = get_polygons("data/polygons/Priority_1_v2019-11-24.kml")

for polygon in polygons:
    mask = vectorized.contains(
        polygon.geometry, receptors.longitude, receptors.latitude
    )
    receptors.loc[mask, "polygon_id"] = polygon.name  # type: ignore


receptors.to_feather("data/receptors/receptors.feather")
