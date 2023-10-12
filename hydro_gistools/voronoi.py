#!/usr/bin/env python
# coding: utf-8

"""
该模块用于生成泰森多边形

- `voronoi_from_shp` - 只通过shp格式点图层生成泰森多边形
- `voronoi_by_geometry` - 根据已有图形范围生成泰森多边形

"""


import os
import whitebox
from functools import partial
from .common import my


def voronoi_from_shp(src, des, data_dir="."):
    """
    通过点直接生成泰森多边形

    References: [VoronoiDiagram](https://whiteboxgeo.com/manual/wbt_book/available_tools/gis_analysis.html?highlight=voro#voronoidiagram)

    Args:
        src (str): 必选，源文件名，即参与运算的Point图层，文件格式为.shp
        des (str): 必选，目标文件名，即运算得到的Polygon图层，文件格式为.shp

    """

    my.data_dir = os.path.abspath(data_dir)
    src = os.path.abspath(src)
    des = os.path.abspath(des)

    wbt = whitebox.WhiteboxTools()

    wbt.voronoi_diagram(src, des, callback=my.my_callback)


import numpy as np
import geopandas as gpd
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union


def voronoi_by_geometry(points, geometry=None, box=None):
    """
    根据点和已知的范围生成泰森多边形

    参考: [Voronoi](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Voronoi.html)

    Args:
        points (GeoDataFrame): 必选，源文件名，即参与运算的Point图层
        geometry (GeoSeries或Shapely.Polygon|MultiPolygon): 可选，目标区域
        box (list): 可选，矩形目标区域，格式如[115,38,136,54],box与geometry二选一即可

    Returns:
        vor_polys (GeoDataFrame): 生成的泰森多边形
    """

    minx, miny, maxx, maxy = points.total_bounds
    if isinstance(geometry, (Polygon, MultiPolygon)):
        minx = geometry.bounds[0] if geometry.bounds[0] < minx else minx
        miny = geometry.bounds[1] if geometry.bounds[1] < miny else miny
        maxx = geometry.bounds[2] if geometry.bounds[2] > maxx else maxx
        maxy = geometry.bounds[3] if geometry.bounds[3] > maxy else maxy

    if isinstance(geometry, gpd.GeoSeries):
        geometry = geometry.to_list()
        geometry = unary_union(geometry)
        minx = geometry.bounds[0] if geometry.bounds[0] < minx else minx
        miny = geometry.bounds[1] if geometry.bounds[1] < miny else miny
        maxx = geometry.bounds[2] if geometry.bounds[2] > maxx else maxx
        maxy = geometry.bounds[3] if geometry.bounds[3] > maxy else maxy

    if box:
        minx = box[0] if box[0] < minx else minx
        miny = box[1] if box[1] < miny else miny
        maxx = box[2] if box[2] > maxx else maxx
        maxy = box[3] if box[3] > maxy else maxy

    deltax = maxx - minx
    deltay = maxy - miny

    pnts = np.dstack((points.geometry.x.to_numpy(), points.geometry.y.to_numpy()))[0]
    gn_pnts = np.concatenate(
        [
            pnts,
            np.array(
                [
                    [100 * deltax + maxx, 100 * deltay + maxy],
                    [100 * deltax + maxx, -100 * deltay + miny],
                    [-100 * deltax + minx, (maxy + miny) / 2],
                ]
            ),
        ]
    )
    vor = Voronoi(gn_pnts)

    bnd_poly = Polygon([[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]])
    if isinstance(geometry, (Polygon, MultiPolygon)):
        bnd_poly = geometry
    if box:
        bnd_poly = Polygon(
            [[box[0], box[1]], [box[2], box[1]], [box[2], box[4]], [box[0], box[4]]]
        )

    vor_polys = []
    for i in range(len(gn_pnts) - 3):
        vor_poly = [vor.vertices[v] for v in vor.regions[vor.point_region[i]]]
        i_cell = bnd_poly.intersection(Polygon(vor_poly))
        if i_cell.geom_type == "Polygon":
            vor_polys.append(Polygon(list(i_cell.exterior.coords[:-1])))
        if i_cell.geom_type == "MultiPolygon":
            for polygon in i_cell:
                vor_polys.append(Polygon(list(polygon.exterior.coords[:-1])))

    vor_polys = gpd.GeoSeries(vor_polys, crs=points.crs)
    vor_polys = vor_polys.to_frame(name="geometry")

    vor_polys = gpd.sjoin(left_df=vor_polys, right_df=points, predicate="intersects")

    return vor_polys
