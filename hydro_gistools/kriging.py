#!/usr/bin/env python
# coding: utf-8

"""
该模块用于进行克里金插值

- `ok` - 即Ordinary Kriging，普通克里金插值方法
"""


import os
import numpy as np

# import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import matplotlib.pyplot as plt
import geopandas as gpd
from osgeo import osr
from osgeo import gdal

import math
import warnings
import io


def ok(
    src,
    des,
    field,
    rows=250,
    columns=250,
    auto_cell_size=False,
    variogram_model="spherical",
    coordinates_type="euclidean",
    n_closest_points=12,
    save_asc=False,
    asc_filename="ok.asc",
):
    """
    即Ordinary Kriging，普通克里金插值方法

    References: [OrdinaryKriging](https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html)

    Args:
        src (str): 必选，源文件名，即参与运算的Point图层，文件格式为.shp，注意：图层坐标最好不要为空或未知
        des (str): 必选，目标文件名，即运算得到的栅格图像，文件格式为.tif
        field (str): 必选，源文件中参与计算的字段名称
        rows (int): 可选，生成图像的行数
        columns (int): 可选，生成图像的列数
        auto_cell_size (bool): 可选，是否自动调整行列数使网格尺寸接近正方形
        variogram_model (str): 可选，变异函数，可以为linear, power, gaussian, spherical, exponential
        coordinates_type (str): 可选，坐标系类型，如果源数据坐标为地理坐标，则为geographic；如果源数据坐标为投影坐标，则为euclidean
        n_closest_points (int): 可选，计算采用临近点的数量
        save_asc (bool) 可选，是否保存中间数据，即.asc格式图像数据
        asc_filename (str) 可选，中间数据.asc格式文件名

    """

    in_point = gpd.GeoDataFrame.from_file(src)

    lons = in_point.geometry.x.values
    lats = in_point.geometry.y.values
    values = in_point[field].values

    box = in_point.geometry.total_bounds
    if auto_cell_size:
        rows, columns = cell_size(box, rows, columns)

    grid_lon = np.linspace(box[0], box[2], columns)
    grid_lat = np.linspace(box[1], box[3], rows)

    OK = OrdinaryKriging(
        lons,
        lats,
        values,
        nlags=6,
        variogram_model=variogram_model,
        coordinates_type=coordinates_type,
    )

    z, ss = OK.execute(
        "grid", grid_lon, grid_lat, backend="C", n_closest_points=n_closest_points
    )

    asc_path = os.path.dirname(os.path.abspath(asc_filename))
    if not os.path.exists(asc_path):
        os.mkdir(asc_path)

    write_asc_grid(grid_lon, grid_lat, z, filename=asc_filename)

    write_asc_tiff(asc_filename, des, in_point.crs.srs, save_asc)


def write_asc_tiff(asc_filename, des, srs, save_asc):
    drv = gdal.GetDriverByName("GTiff")
    ds_in = gdal.Open(asc_filename)
    ds_out = drv.CreateCopy(des, ds_in)
    ds_out.SetProjection(srs)
    ds_in = None
    ds_out = None

    if save_asc and os.path.exists(asc_filename):
        os.remove(asc_filename)


def cell_size(box, rows, columns):
    delta_row = (box[2] - box[0]) / rows
    delta_column = (box[3] - box[1]) / columns
    delta = min(delta_row, delta_column)

    r = (box[2] - box[0]) / delta
    c = (box[3] - box[1]) / delta
    rows = int(r)
    columns = int(c)

    while math.fabs(r - rows) / rows > 0.01 or math.fabs(c - columns) / columns >= 0.01:
        if math.fabs(r - rows) / rows > 0.01:
            rows = rows + 1

        if math.fabs(c - columns) / columns > 0.01:
            columns = columns + 1

        delta_row = (box[2] - box[0]) / rows
        delta_column = (box[3] - box[1]) / columns
        delta = min(delta_row, delta_column)

        r = (box[2] - box[0]) / delta
        c = (box[3] - box[1]) / delta
        rows = int(r)
        columns = int(c)

    return rows, columns


def write_asc_grid(x, y, z, filename="output.asc", no_data=-999.0, style=1):
    if np.ma.is_masked(z):
        z = np.array(z.tolist(no_data))

    x = np.squeeze(np.array(x))
    y = np.squeeze(np.array(y))
    z = np.squeeze(np.array(z))
    nrows = z.shape[0]
    ncols = z.shape[1]

    if z.ndim != 2:
        raise ValueError("Two-dimensional grid is required to write *.asc grid.")
    if x.ndim > 1 or y.ndim > 1:
        raise ValueError(
            "Dimensions of X and/or Y coordinate arrays are not "
            "as expected. Could not write *.asc grid."
        )
    if z.shape != (y.size, x.size):
        warnings.warn(
            "Grid dimensions are not as expected. "
            "Incorrect *.asc file generation may result.",
            RuntimeWarning,
        )
    if np.amin(x) != x[0] or np.amin(y) != y[0]:
        warnings.warn(
            "Order of X or Y coordinates is not as expected. "
            "Incorrect *.asc file generation may result.",
            RuntimeWarning,
        )

    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])
    if not np.isclose(abs((x[-1] - x[0]) / (x.shape[0] - 1)), dx) or not np.isclose(
        abs((y[-1] - y[0]) / (y.shape[0] - 1)), dy
    ):
        raise ValueError(
            "X or Y spacing is not constant; *.asc grid cannot be written."
        )
    cellsize = -1
    if style == 2:
        if dx != dy:
            raise ValueError(
                "X and Y spacing is not the same. "
                "Cannot write *.asc file in the specified format."
            )
        cellsize = dx

    xllcenter = x[0]
    yllcenter = y[0]

    # Note that these values are flagged as -1. If there is a problem in trying
    # to write out style 2, the -1 value will appear in the output file.
    xllcorner = -1
    yllcorner = -1
    if style == 2:
        xllcorner = xllcenter - dx / 2.0
        yllcorner = yllcenter - dy / 2.0

    with io.open(filename, "w") as f:
        if style == 1:
            f.write("NCOLS          " + "{:<n}".format(ncols) + "\n")
            f.write("NROWS          " + "{:<n}".format(nrows) + "\n")
            f.write("XLLCENTER      " + "{:<f}".format(xllcenter) + "\n")
            f.write("YLLCENTER      " + "{:<f}".format(yllcenter) + "\n")
            f.write("DX             " + "{:<f}".format(dx) + "\n")
            f.write("DY             " + "{:<f}".format(dy) + "\n")
            f.write("NODATA_VALUE   " + "{:<f}".format(no_data) + "\n")
        elif style == 2:
            f.write("NCOLS          " + "{:<n}".format(ncols) + "\n")
            f.write("NROWS          " + "{:<n}".format(nrows) + "\n")
            f.write("XLLCORNER      " + "{:<f}".format(xllcorner) + "\n")
            f.write("YLLCORNER      " + "{:<f}".format(yllcorner) + "\n")
            f.write("CELLSIZE       " + "{:<f}".format(cellsize) + "\n")
            f.write("NODATA_VALUE   " + "{:<f}".format(no_data) + "\n")
        else:
            raise ValueError("style kwarg must be either 1 or 2.")

        for m in range(z.shape[0] - 1, -1, -1):
            for n in range(z.shape[1]):
                f.write("{:<16.2f}".format(z[m, n]))
            if m != 0:
                f.write("\n")
