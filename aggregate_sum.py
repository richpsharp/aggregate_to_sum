"""Aggregate sum."""
import argparse
import glob
import logging
import os
import sys

from osgeo import gdal
import numpy
import pygeoprocessing

LOGGER = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.DEBUG,
    format=('%(message)s'),
    stream=sys.stdout)
logging.getLogger('ecoshard').setLevel(logging.INFO)
LOGGER = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Aggregate by sum.')
    parser.add_argument(
        'target_size', help=(
            "target pixel size in degrees"), type=float)
    parser.add_argument(
        'filepath', nargs='+', help='Files/patterns to ecoshard.')

    args = parser.parse_args()

    gtiff_driver = gdal.GetDriverByName('GTiff')

    LOGGER.debug('filepath: %s', args.filepath)
    for glob_pattern in args.filepath:
        LOGGER.debug(glob_pattern)
        for file_path in glob.glob(glob_pattern):
            LOGGER.debug(args.target_size)
            target_path = (
                'sum_aggregate_to_%f_%s' % (
                    float(args.target_size), os.path.basename(file_path)))
            LOGGER.debug('making %s', target_path)
            LOGGER.debug(file_path)
            base_info = pygeoprocessing.get_raster_info(file_path)
            base_cols, base_rows = base_info['raster_size']
            base_gt = base_info['geotransform']
            target_gt = [
                base_gt[0], float(args.target_size), 0,
                base_gt[1], 0, -float(args.target_size)]

            base_inv_gt = gdal.InvGeoTransform(base_gt)
            target_inv_gt = gdal.InvGeoTransform(target_gt)

            nodata = base_info['nodata'][0]
            n_cols = int(base_gt[1] / target_gt[1] * base_cols)
            n_rows = int(base_gt[5] / target_gt[5] * base_rows)

            base_raster = gdal.OpenEx(file_path, gdal.OF_RASTER)
            base_band = base_raster.GetRasterBand(1)

            target_raster = gtiff_driver.Create(
                target_path, n_cols, n_rows, 1,
                base_info['datatype'],
                options=(
                    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
                    'BLOCKXSIZE=256', 'BLOCKYSIZE=256'))
            target_band = target_raster.GetRasterBand(1)
            target_band.SetNoDataValue(nodata)
            target_band.Fill(nodata)
            target_raster.SetProjection(base_info['projection'])
            target_raster.SetGeoTransform(target_gt)

            for target_i in range(n_cols):
                LOGGER.debug(
                    '%.2f%% complete for %s', target_i/n_cols*100, target_path)
                for target_j in range(n_rows):
                    target_x, target_y = gdal.ApplyGeoTransform(
                        target_gt, target_i, target_j)
                    base_i, base_j = [
                        int(x) for x in gdal.ApplyGeoTransform(
                            base_inv_gt, target_x, target_y)]

                    target_x_p1, target_y_p1 = gdal.ApplyGeoTransform(
                        target_gt, target_i+1, target_j+1)
                    base_i_p1, base_j_p1 = [
                        int(x) for x in gdal.ApplyGeoTransform(
                            base_inv_gt, target_x_p1, target_y_p1)]

                    win_xsize = int(base_i_p1-base_i)
                    win_ysize = int(base_j_p1-base_j)

                    LOGGER.debug(
                        'B: %d %d, %d %d %d %d', base_cols, base_rows, base_i, base_j, win_xsize, win_ysize)
                    if win_xsize + base_i >= base_rows:
                        win_xsize = base_rows-base_i
                    if win_ysize + base_j >= base_cols:
                        win_ysize = base_cols-base_j
                    LOGGER.debug(
                        'A: %d %d, %d %d %d %d', base_cols, base_rows, base_i, base_j, win_xsize, win_ysize)

                    base_array = base_band.ReadAsArray(
                        xoff=int(base_i), yoff=int(base_j),
                        win_xsize=int(base_i_p1-base_i),
                        win_ysize=int(base_j_p1-base_j))

                    masked_array = base_array[
                        ~numpy.isclose(base_array, nodata)]

                    if masked_array.size > 1:
                        array_sum = numpy.sum(masked_array)
                        LOGGER.debug(array_sum)
                        target_band.WriteArray(
                            [[array_sum]], target_i, target_j)
