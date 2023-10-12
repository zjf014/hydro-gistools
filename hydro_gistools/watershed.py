#!/usr/bin/env python
# coding: utf-8

"""
该模块用于划分子流域，目前大尺度流域自动划分还有问题

- `watershed` - 自动划分子流域

"""


import os
import whitebox
from .common import my


def watershed(
    dem, src, des, data_dir=".", save_intermediate_data=False, intermediate_data_dir="."
):
    """
    划分子流域

    Todo:
        - 填缺: [fill_missing_data](https://whiteboxgeo.com/manual/wbt_book/available_tools/geomorphometric_analysis.html?#fillmissingdata)
        - 填洼: [fill_depressions](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#FillDepressions)
        - 流向: [d8_pointer](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#D8Pointer)
        - 流量累积: [d8_flow_accumulation](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#D8FlowAccumulation)
        - 河网分支分级: [strahler_order_basins](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#StrahlerOrderBasins)
        - 出水口校正: [snap_pour_points](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#SnapPourPoints)
        - 生成子流域: [watershed](https://whiteboxgeo.com/manual/wbt_book/available_tools/hydrological_analysis.html#Watershed)
        - 转化为矢量图层: [raster_to_vector_polygons](https://whiteboxgeo.com/manual/wbt_book/available_tools/data_tools.html?#rastertovectorpolygons)

    Args:
        dem (str): 必选，DEM文件名，格式为.tif
        src (str): 必选，源文件名，即参与运算的Point图层，文件格式为.shp
        des (str): 必选，目标文件名，即运算得到的Polygon图层，文件格式为.shp
        save_intermediate_data (bool): 可选，是否保存中间过程数据
        intermediate_data_dir (str): 可选，中间过程数据保存目录

    """

    my.data_dir = os.path.abspath(data_dir)
    dem = os.path.abspath(dem)
    src = os.path.abspath(src)
    des = os.path.abspath(des)

    intermediate_data_dir = os.path.abspath(intermediate_data_dir)
    if os.path.isdir(intermediate_data_dir) and not os.path.exists(
        intermediate_data_dir
    ):
        os.mkdir(intermediate_data_dir)

    wbt = whitebox.WhiteboxTools()

    try:
        wbt.fill_missing_data(
            dem,
            os.path.join(intermediate_data_dir, "fillmissingdata.tif"),
            filter=11,
            weight=2.0,
            no_edges=True,
            callback=my.my_callback,
        )

        wbt.fill_depressions_wang_and_liu(
            os.path.join(intermediate_data_dir, "fillmissingdata.tif"),
            os.path.join(intermediate_data_dir, "filldepression.tif"),
            fix_flats=True,
            flat_increment=None,
            # max_depth=None,
            callback=my.my_callback,
        )

        wbt.d8_pointer(
            os.path.join(intermediate_data_dir, "filldepression.tif"),
            os.path.join(intermediate_data_dir, "d8pointer.tif"),
            esri_pntr=False,
            callback=my.my_callback,
        )

        wbt.d8_flow_accumulation(
            os.path.join(intermediate_data_dir, "d8pointer.tif"),
            os.path.join(intermediate_data_dir, "d8flowaccumulation.tif"),
            out_type="cells",
            log=False,
            clip=False,
            pntr=True,
            esri_pntr=False,
            callback=my.my_callback,
        )

        wbt.strahler_order_basins(
            os.path.join(intermediate_data_dir, "d8pointer.tif"),
            os.path.join(intermediate_data_dir, "d8flowaccumulation.tif"),
            os.path.join(intermediate_data_dir, "strahler.tif"),
            esri_pntr=False,
            callback=my.my_callback,
        )

        wbt.snap_pour_points(
            src,
            os.path.join(intermediate_data_dir, "d8flowaccumulation.tif"),
            os.path.join(intermediate_data_dir, "snapped.shp"),
            500,
            callback=my.my_callback,
        )

        wbt.watershed(
            os.path.join(intermediate_data_dir, "d8pointer.tif"),
            os.path.join(intermediate_data_dir, "snapped.shp"),
            os.path.join(intermediate_data_dir, "watershed.tif"),
            esri_pntr=False,
            callback=my.my_callback,
        )

        wbt.raster_to_vector_polygons(
            os.path.join(intermediate_data_dir, "watershed.tif"),
            des,
            callback=my.my_callback,
        )
    finally:
        os.chdir(my.data_dir)

        if not save_intermediate_data:
            delete_intermediate_data(intermediate_data_dir)


def delete_intermediate_data(data_dir):
    if os.path.exists(os.path.join(data_dir, "fillmissingdata.tif")):
        os.remove(os.path.join(data_dir, "fillmissingdata.tif"))
    if os.path.exists(os.path.join(data_dir, "filldepression.tif")):
        os.remove(os.path.join(data_dir, "filldepression.tif"))
    if os.path.exists(os.path.join(data_dir, "d8pointer.tif")):
        os.remove(os.path.join(data_dir, "d8pointer.tif"))
    if os.path.exists(os.path.join(data_dir, "d8flowaccumulation.tif")):
        os.remove(os.path.join(data_dir, "d8flowaccumulation.tif"))
    if os.path.exists(os.path.join(data_dir, "watershed.tif")):
        os.remove(os.path.join(data_dir, "watershed.tif"))
    if os.path.exists(os.path.join(data_dir, "strahler.tif")):
        os.remove(os.path.join(data_dir, "strahler.tif"))
    if os.path.exists(os.path.join(data_dir, "snapped.shp")):
        files = os.listdir(data_dir)
        for file in files:
            if file.split(".")[0] == "snapped":
                os.remove(os.path.join(data_dir, file))
