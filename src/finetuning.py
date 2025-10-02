import pickle
import os
import fnmatch
import time
import sys
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import warnings

import SimpleITK as sitk # type: ignore
from tifffile import TiffFile # type: ignore
from tifffile.tiffcomment import tiffcomment
import typing as tp
from pathlib import Path
import zlib
from xml.etree import ElementTree
import napari # type: ignore
import vispy.color
from shapely.geometry import mapping, Polygon
from shapely.ops import unary_union
from scipy.ndimage import rotate

from src.alignment import remove_params, extract_downscaled_images
from src.alignment import get_cells_coordinates_SPIAT_CellType, get_annotations
from src.alignment import alignment_report

from src.save_images import panels_name_alignment, generate_channels_list, opal_to_rgb, getLabels

from src.alignment import merge_annotations, get_gdf, transform_annotation, rotate_coordinates_angle, scale_multipolygon_coordinates, plot_multipolygon
from src.alignment import transform_filter_coordinates, filter_coordinates
from src.alignment import save_tables, plot_rasters




def finetuning(id, meshsize, downscaled_images_path, coordinate_tables, visualization, scans_path, resolution_micron, metric, pixel_size_raster_micron, alpha_red):
    ## Load and read the pickle file
    with open(downscaled_images_path, "rb") as file:
        downscaled_images = pickle.load(file)
    ## Get the parameters
    annotations_paths = downscaled_images["params"].get("annotations_paths")
    annotations_names_AnalysisArea = downscaled_images["params"].get("annotations_names_AnalysisArea")
    annotations_names_empty = downscaled_images["params"].get("annotations_names_empty")
    annotations_names_artefacts = downscaled_images["params"].get("annotations_names_artefacts")
    HALO_rotation_path = downscaled_images["params"].get("HALO_rotation_path")
    panels_all = downscaled_images["params"].get("panels")
    output_path = downscaled_images["params"].get("output_path")
    ## Remove "params" from downscaled_images
    downscaled_images = remove_params(downscaled_images)
    ## Get reference panel
    reference_panel = panels_all[1]
    ## Get the downscaled images of the given patient id
    downscaled_images_id = downscaled_images[id]
    
    print("---------------")
    print(f'Patient: {id}')
    metadata_images = {}
    cell_coordinates = {}
    data_frame_cells = {}
    metric_ms_alignment = {}
    img_resize_alignment = {}
    outTx_Rigid_alignment = {}
    outTx_Bspline_alignment = {}
    artefacts_empty_alignment = {}
    analysis_area_alignment = {}
    manual_empty_alignment = {}
    scans_alignment_dict = panels_name_alignment(panels_all, scans_path)

    # Iterate over the alignment with the downscaled images and the files to the full resolution .qptiff files for the corresponding alignment
    for i_ms, ((name_alignment, scans_alignment_paths), (_, downscaled_images_id_name_alignment)) in enumerate(zip(scans_alignment_dict.items(), downscaled_images_id.items())):
        print("-----")
        ## Load downscaled images
        print("Loading downscaled images...")
        panels_alignment = name_alignment.split("_")
        metadata_images, img1, img2, img1_resize, img2_resize = extract_downscaled_images(downscaled_images_id_name_alignment, panels_alignment, 
                                                                                            name_alignment, metadata_images, id)
        ## Alignment
        print(f"Alignment {name_alignment}...")
        spline_order = 3
        (outTx_Bspline_dict, simg2_dict, execution_time_dict, 
            metric_values_dict, outTx_Rigid, simg1_Rigid, 
            simg2_Rigid, nda_Rigid) = alignment_2panels_specific_mesh_size(img1, img2, meshsize[i_ms], metric, spline_order)
        ## Load coordinate table and annotations of each panel
        print(f'Loading tables and annotations...')
        for panel in panels_alignment:
            # Get the panel path
            index = panels_all.index(panel)
            coordinate_table = coordinate_tables[index]
            # Get the file path corresponding to id
            csv_file_path = next((os.path.join(coordinate_table, f) for f in os.listdir(coordinate_table) if id in f and f.endswith(".csv") and not f.startswith(".")), None)
            ## Get the coordinates table and annotations
            scale_percent = metadata_images[name_alignment][f'scale_percent_{panel}']
            image_shape = metadata_images[name_alignment][f'image_shape_{panel}']
            cell_coordinates, data_frame_cells = get_cells_coordinates_SPIAT_CellType(csv_file_path, id, panel, cell_coordinates, data_frame_cells, resolution_micron, HALO_rotation_path, scale_percent, image_shape)
            if annotations_paths:
                annotations = annotations_paths[index]
                annotation_file_path = next(
                    (os.path.join(annotations, f) for f in os.listdir(annotations) 
                    if id in f and (f.endswith(".annotations") or f.endswith(".geojson")) and not f.startswith(".")), 
                    None
                )
                artefacts_empty_alignment, analysis_area_alignment = get_annotations(annotation_file_path, panel, artefacts_empty_alignment, analysis_area_alignment, annotations_names_empty, annotations_names_artefacts, annotations_names_AnalysisArea)
        ## Create the alignment report
        print("Creating the alignment report...")
        outTx_Rigid_alignment, outTx_Bspline_alignment, img_resize_alignment, metric_ms_alignment = alignment_report(id, name_alignment, panels_alignment, cell_coordinates, metadata_images, output_path, 
                                                                                                alpha_red, img1_resize, img2_resize, simg2_Rigid, outTx_Rigid, outTx_Bspline_dict, simg2_dict, 
                                                                                                execution_time_dict, metric_values_dict, spline_order, resolution_micron, pixel_size_raster_micron, 
                                                                                                metric_ms_alignment, img_resize_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment)
    

        ## double cursor visualization
        scan_path1, scan_path2 = [path for path in scans_alignment_paths]
        warnings.filterwarnings("ignore")
        if visualization == "all":
            manual_empty_alignment = mirrored_cursor_visu(scan_path1, scan_path2, id, metadata_images, img1_resize, img2_resize, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, manual_empty_alignment, analysis_area_alignment, artefacts_empty_alignment)            
        elif visualization == name_alignment:
            manual_empty_alignment = mirrored_cursor_visu(scan_path1, scan_path2, id, metadata_images, img1_resize, img2_resize, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, manual_empty_alignment, analysis_area_alignment, artefacts_empty_alignment)

    response = input("Are you sure you want to rewrite the file of the merged annotations with the specified mesh size for this patient's panel alignments? (Y/n): ").strip().lower()
    if response in ['no', 'n']:
        print("Exiting the script without rewriting the annotations file.")
        sys.exit()  
    inner_boundary = None
    if annotations_paths:
        ## Transform and merge the annotations
        print("-----")
        print("Merging the annotations...") 
        mask = merge_annotations(id, artefacts_empty_alignment, analysis_area_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, 
                            img_resize_alignment, metric_ms_alignment, metadata_images, panels_all, reference_panel, resolution_micron, output_path)
        
        if visualization != "0":
            mask = merge_manual_empty_with_mask(manual_empty_alignment, mask, output_path, id, resolution_micron)
    else:
        if manual_empty_alignment:
            panels_tr = [key for key in manual_empty_alignment.keys() if "in" in key]
            manual_empty_alignment_list = []
            for panel in panels_tr:
                manual_empty_alignment_list.append(manual_empty_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))

            ## Filter out empty GeoSeries
            manual_empty_alignment_list = [g for g in manual_empty_alignment_list if not g.empty]
            ## Merge the remaining GeoSeries
            merged_polygons = gpd.GeoSeries(unary_union([g.unary_union for g in manual_empty_alignment_list]))

            # Create a Polygon representing the inner boundary
            inner_boundary = merged_polygons.unary_union
            mask = "inner_boundary"
            print(mask)
            
        else:
            mask = "No_annotations"

    ## Transform and filter coordinates
    print("Transform, filer and merge coordinates...")
    merged_cell_coordinates = pd.DataFrame()
    for name_alignment in downscaled_images_id.keys():
        panels_alignment = name_alignment.split("_")
        merged_cell_coordinates = transform_filter_coordinates(metadata_images, cell_coordinates, data_frame_cells, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, mask, resolution_micron, merged_cell_coordinates, inner_boundary)
    merged_cell_coordinates = filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, resolution_micron, merged_cell_coordinates, inner_boundary)
    ## Save new coordinates table
    save_tables(merged_cell_coordinates, output_path, id)
    ## Plot cell coordinates before and after alignment
    plot_rasters(data_frame_cells, merged_cell_coordinates, cell_coordinates, pixel_size_raster_micron, output_path, panels_all, id, resolution_micron)




def alignment_2panels_specific_mesh_size(img1, img2, meshsize, metric, spline_order):
    ## Rigid alignment
    execution_time_dict = {}
    metric_values_dict = {}
    global metric_values_list
    metric_values_list= []
    start_time = time.time() 
    cimg, simg1_Rigid, simg2_Rigid, outTx_Rigid = imageRegRigid(img1,img2)
    end_time = time.time() 
    execution_time = end_time - start_time  
    # print(f"Execution time for Rigid: {execution_time} seconds")
    nda_Rigid = sitk.GetArrayFromImage(cimg)
    # plt.imshow(nda)
    
    execution_time_dict[f"{metric}_Rigid"] = execution_time
    metric_values_dict[f"{metric}_Rigid"] = metric_values_list

    
    ## B-spline alignment
    # Define variables
    simg2_dict = {}
    outTx_Bspline_dict = {}
    # print("Alignment on mesh size =",meshsize)
    
    transformDomainMeshSize = [meshsize, meshsize]
    metric_values_list = []

    img1 = sitk.GetImageFromArray(simg1_Rigid)
    img2 = sitk.GetImageFromArray(simg2_Rigid)

    img1 = sitk.Cast(img1, sitk.sitkFloat32)
    img2 = sitk.Cast(img2, sitk.sitkFloat32)

    start_time = time.time() 
    cimg, simg1, simg2, outTx_Bspline = imageRegBspline(img1, img2, transformDomainMeshSize, spline_order, metric)
    end_time = time.time() 
    execution_time = end_time - start_time 
    # print(f"Execution time for deformation degree {deformation_degree}: {execution_time} seconds")

    execution_time_dict[f"{metric}_{meshsize}"] = execution_time
    metric_values_dict[f"{metric}_{meshsize}"] = metric_values_list
    simg2_dict[f"{metric}_{meshsize}"] = simg2
    outTx_Bspline_dict[f"{metric}_{meshsize}"] = outTx_Bspline

    simg1 = sitk.GetArrayFromImage(simg1)
    simg2 = sitk.GetArrayFromImage(simg2)
    nda = sitk.GetArrayFromImage(cimg)

    return outTx_Bspline_dict, simg2_dict, execution_time_dict, metric_values_dict, outTx_Rigid, simg1_Rigid, simg2_Rigid, nda_Rigid



def command_iteration(method):
    global metric_values_list
    metric_values_list.append(method.GetMetricValue())
    '''
    print(
        f"{method.GetOptimizerIteration():3} "
        + f"= {method.GetMetricValue():10.5f}"
    )
    '''

def imageRegRigid(fixed,moving):

    R = sitk.ImageRegistrationMethod()

    R.SetMetricAsCorrelation()

    R.SetOptimizerAsRegularStepGradientDescent(
        learningRate=2.0,
        minStep=1e-4,
        numberOfIterations=1000,
        gradientMagnitudeTolerance=1e-8,
    )
    R.SetOptimizerScalesFromIndexShift()

    tx = sitk.CenteredTransformInitializer(
        fixed, moving, sitk.Euler2DTransform()
    )
    
    
    R.SetInitialTransform(tx)

    R.SetInterpolator(sitk.sitkLinear)

    R.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(R))

    outTx = R.Execute(fixed, moving)

    # print("-------")
    # print(outTx)
    # print(f"Optimizer stop condition: {R.GetOptimizerStopConditionDescription()}")
    # print(f" Iteration: {R.GetOptimizerIteration()}")
    # print(f" Metric value: {R.GetMetricValue()}")

    # sitk.WriteTransform(outTx, "output.tfm")

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed)
    resampler.SetInterpolator(sitk.sitkLinear)
    # resampler.SetDefaultPixelValue(100)
    resampler.SetTransform(outTx)

    out = resampler.Execute(moving)

    simg1 = sitk.Cast(sitk.RescaleIntensity(fixed), sitk.sitkUInt8)
    simg2 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
    cimg = sitk.Compose(simg1, simg2, simg1 // 2.0 + simg2 // 2.0)
    
    simg1 = sitk.GetArrayFromImage(simg1)
    simg2 = sitk.GetArrayFromImage(simg2)

    return cimg, simg1, simg2, outTx

def imageRegBspline(fixed,moving,transformDomainMeshSize,spline_order,metric):
    # print(transformDomainMeshSize)
    tx = sitk.BSplineTransformInitializer(fixed, transformDomainMeshSize, spline_order)

    # print("Spline order:")
    # print(tx.GetOrder())

    R = sitk.ImageRegistrationMethod()
    # Correlation between pixel values of 2 images
    if metric == "Mutual information":
        R.SetMetricAsMattesMutualInformation()
    elif metric == "Correlation":
        R.SetMetricAsCorrelation()

    R.SetOptimizerAsLBFGSB(
        gradientConvergenceTolerance=1e-10,
        numberOfIterations=1000,
        maximumNumberOfCorrections=100,
        maximumNumberOfFunctionEvaluations=1000,
        costFunctionConvergenceFactor=1e7,
    )
    R.SetInitialTransform(tx, True)
    R.SetInterpolator(sitk.sitkLinear)
    
    R.AddCommand(sitk.sitkIterationEvent, lambda: command_iteration(R))

    outTx = R.Execute(fixed, moving)

    # print("-------")
    # print(outTx)
    # print(f"Optimizer stop condition: {R.GetOptimizerStopConditionDescription()}")
    # print(f" Iteration: {R.GetOptimizerIteration()}")
    # print(f" Metric value: {R.GetMetricValue()}")

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed)
    resampler.SetInterpolator(sitk.sitkLinear)
    # resampler.SetDefaultPixelValue(100)
    resampler.SetTransform(outTx)

    out = resampler.Execute(moving)
    
    simg1 = sitk.Cast(sitk.RescaleIntensity(fixed), sitk.sitkUInt8)
    simg2 = sitk.Cast(sitk.RescaleIntensity(out), sitk.sitkUInt8)
    cimg = sitk.Compose(simg1, simg2, simg1 // 2.0 + simg2 // 2.0)
    
    return cimg, simg1, simg2, outTx



def load_QPTIFF_high_res(path_scan):
    tif = TiffFile(path_scan)
    tif_tags = {tag.name: tag.value for tag in tif.pages[0].tags.values()}
    xml_content = tif_tags.get("ImageDescription", "").replace("\r\n", "").replace("\t", "")
    try:
        channels = generate_channels_list(xml_content, tif_tags)
    except ValueError as e:
        channel_list, channel_name_dictionary = getLabels(tif_tags)
        # Build channels in same format as JSON-based structure
        channels = []
        for i, label in enumerate(channel_list):
            # Extract fluorophore name from label: e.g., "Opal 570" from "Opal 570 (channel 1)"
            name_part = label.split(" (")[0]
            rgb = opal_to_rgb(name_part)

            channels.append({
                "id": i + 1,
                "name": label,
                "rgb": rgb
            })
    # Load the second compressed image channels in .pages
    img_data = tif.series[0].levels[2].asarray()
    sizeY_fullres2 = tif_tags['ImageLength']
    sizeY_compressed2 = img_data[0].shape[0]

    return channels, img_data, sizeY_compressed2, sizeY_fullres2

def load_TIF_high_res(path_scan):
    tif = TiffFile(path_scan, is_ome=False)
    xml = tiffcomment(path_scan)
    root = ElementTree.fromstring(xml.replace("\n", "").replace("\t", ""))
    sizeY = list(dict.fromkeys([x.get("sizeY") for x in root.iter() if x.tag == "dimension"]))
    # Load the second compressed image channels in .levels
    img_data = tif.series[0].levels[2].asarray()
    sizeY_fullres = int(sizeY[0])
    sizeY_compressed = int(sizeY[2])
    channels = [
        {"id": elem.get("id"), "name": elem.get("name"), "rgb": elem.get("rgb")}
        for elem in root.iter() if elem.tag == "channel"
    ]



    return channels, img_data, sizeY_compressed, sizeY_fullres


def get_channel_names_from_indica_xml(file: Path) -> list[str]:
    """
    Extract the name of channels in IndicaLabs TIFF file.

    Parameters
    ----------
    file: pathlib.Path
        Path to TIFF file.
    """
    xml = tiffcomment(file)
    root = ElementTree.fromstring(xml.replace("\n", "").replace("\t", ""))
    # channels = [x.attrib for x in root.iter() if x.tag == 'channel']
    channels = [x.get("name") for x in root.iter() if x.tag == "channel"]
    return channels

def load_high_res_imgs(folder_path1, folder_path2, id, panels):
    # Define path of the qptiff or tiff
    path1 = "\n".join([os.path.join(folder_path1, f) for f in os.listdir(folder_path1) if not f.startswith('.') and (fnmatch.fnmatch(f.lower(), f"*{id}*.qptiff") or fnmatch.fnmatch(f.lower(), f"*{id}*.tif"))])
    path2 = "\n".join([os.path.join(folder_path2, f) for f in os.listdir(folder_path2) if not f.startswith('.') and (fnmatch.fnmatch(f.lower(), f"*{id}*.qptiff") or fnmatch.fnmatch(f.lower(), f"*{id}*.tif"))])

    # Read files whether it is qptiff or tif
    ext_to_loader = {
        ".qptiff": load_QPTIFF_high_res,
        ".tif": load_TIF_high_res,
    }
    def get_loader(path: str):
        for ext, loader in ext_to_loader.items():
            if path.lower().endswith(ext):
                return loader, ext
        raise ValueError(f"Unsupported file type: {path}")
    # Get loaders + extensions
    loader1, ext1 = get_loader(path1)
    loader2, ext2 = get_loader(path2)
    # Flags
    tif_1 = ext1 == ".tif"
    tif_2 = ext2 == ".tif"
    # Load data
    channels1, img1_data, sizeY_compressed1, sizeY_fullres1 = loader1(path1)
    channels2, img2_data, sizeY_compressed2, sizeY_fullres2 = loader2(path2)

    scale_percent1 = sizeY_compressed1/sizeY_fullres1
    scale_percent2 = sizeY_compressed2/sizeY_fullres2

    '''
    # Saving full resolution images
    np.save(dir_res + f"{id}/{folder_name_alignment}/results/img1_data_{id}_panel_{panels[0]}.npy", img1_data)
    np.save(dir_res + f"{id}/{folder_name_alignment}/results/img2_data_{id}_panel_{panels[0]}.npy", img2_data)
    '''

    return img1_data, channels1, img2_data, channels2, tif_1, tif_2



def mirrored_cursor_visu(scan_path1, scan_path2, id, metadata_images, img1_resize, img2_resize, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, manual_empty_alignment, analysis_area_alignment = None, artefacts_empty_alignment = None, shapes_data_alignment=None, shapes_transformed=None, full_res=None):
    panels = name_alignment.split('_')
    img1_data, channels1, img2_data, channels2, tif_1, tif_2 = load_high_res_imgs(scan_path1, scan_path2, id, panels)
    scale_percent1 = metadata_images[name_alignment][f"scale_percent_{panels_alignment[0]}"]
    scale_percent2 = metadata_images[name_alignment][f"scale_percent_{panels_alignment[1]}"]
    crop_coords1 = metadata_images[name_alignment][f'crop_coords_{panels_alignment[0]}']
    crop_coords2 = metadata_images[name_alignment][f'crop_coords_{panels_alignment[1]}']
    image_shape1 = metadata_images[name_alignment][f'image_shape_{panels_alignment[0]}']
    image_shape2 = metadata_images[name_alignment][f'image_shape_{panels_alignment[1]}']
    scale_percent1_napari = image_shape1[0] / img1_data.shape[1]
    scale_percent2_napari = image_shape2[0] / img2_data.shape[1]
    manual_alignment_displacement = metadata_images[name_alignment]['manual_alignment_displacement']
    manual_alignment_rotation = metadata_images[name_alignment]['manual_alignment_rotation']
    manual_alignment_rotation_shape = metadata_images[name_alignment]['manual_alignment_rotation_shape']
    image_shape_manual_alignment = metadata_images[name_alignment]['image_shape_manual_alignment']
    print(f"Mirrored cursor visualization for alignment {name_alignment}")
    outTx_Rigid = outTx_Rigid_alignment[panels_alignment[0]]
    outTx_Bspline = outTx_Bspline_alignment[panels_alignment[0]]
    outTx_Rigid = outTx_Rigid_alignment.get(panels_alignment[0])
    outTx_Bspline = outTx_Bspline_alignment.get(panels_alignment[0])

    if outTx_Rigid is None or outTx_Bspline is None:
        raise ValueError(f"Missing data for key: {panels_alignment[0]}")



    
    ## Image low resolution but fast
    # By default
    ## Image full resolution but lags
    # full_res=True
    
    def int_to_rgb(value):
        value = int(value)
        r = (value >> 16) & 255
        g = (value >> 8) & 255
        b = value & 255
        return [r / 255.0, g / 255.0, b / 255.0] 


    viewer1 = napari.Viewer()
    viewer1.title = panels[0]

    ## Rotate image in for better visualization of the alignment

    def rotate_image_with_offset(img, angle):
        """
        Rotate a high-dimensional image (C, H, W) with expand=True and return the rotated image
        and the offset needed to keep the original center fixed.

        Parameters:
            img (np.ndarray): Image of shape (C, H, W)
            angle (float): Rotation angle in degrees

        Returns:
            rotated_img (np.ndarray): Rotated image with expand=True
            offset (tuple): (offset_y, offset_x) to re-center image
        """
        channels, height, width = img.shape
        img_shape = img.shape[1:]

        # Rotate each channel separately with expand=True
        rotated_channels = []
        print('rotating image')
        for c in range(channels):
            rotated = rotate(img[c], angle, reshape=True, order=1) 
            rotated_channels.append(rotated)

        # Stack back into shape (C, new_H, new_W)
        rotated_img = np.stack(rotated_channels, axis=0)

        # Original and new center
        orig_center = np.array([height // 2, width // 2])
        new_center = np.array([rotated_img.shape[1] // 2, rotated_img.shape[2] // 2])

        # Calculate offset to re-align rotated image center with original center
        offset = orig_center - new_center 

        return rotated_img, img_shape, offset
    
    img1_data, img1_data_shape, img1_data_offset = rotate_image_with_offset(img1_data, manual_alignment_rotation)

    # Set contrast limits based on whether it is .tif file
    contrast_limit_1 = [0, 50] if tif_1 else [0, 255]
    contrast_limit_2 = [0, 50] if tif_2 else [0, 255]

    for i, channel in enumerate(channels1):
        rgb_values = int_to_rgb(channel['rgb'])
        colorMap = vispy.color.Colormap([[0.0, 0.0, 0.0], rgb_values])
        
        if full_res:
            rlayer = viewer1.add_image(
                [img1_data[i],
                img1_data[i][::4, ::4],
                img1_data[i][::8, ::8]],
                name=channel['name'],
                contrast_limits=contrast_limit_1
            )
        else:
            rlayer = viewer1.add_image(
                img1_data[i],
                name=channel['name'],
                contrast_limits=contrast_limit_1
            )
            
        rlayer.blending = 'additive'
        rlayer.colormap = colorMap
        #rlayer.gamma = gamma_value
        
    viewer2 = napari.Viewer()
    viewer2.title = panels[1]
    for i, channel in enumerate(channels2):
        rgb_values = int_to_rgb(channel['rgb'])
        colorMap = vispy.color.Colormap([[0.0, 0.0, 0.0], rgb_values])
        if full_res:
            rlayer = viewer2.add_image(
                [img2_data[i],
                img2_data[i][::4, ::4],
                img2_data[i][::8, ::8]],
                name=channel['name'], 
                contrast_limits=contrast_limit_2
            )
        else:
            rlayer = viewer2.add_image(
                img2_data[i],
                name=channel['name'], 
                contrast_limits=contrast_limit_2
            )
            
        rlayer.blending = 'additive'
        rlayer.colormap = colorMap
        #rlayer.gamma = gamma_value
    
    # Create a function to synchronize zooming and center between viewer1 and viewer2
    def sync_zoom_center(viewer_src, viewer_dest):
        with viewer_dest.dims.events.ndisplay.blocker(), viewer_dest.camera.events.center.blocker():
            viewer_dest.dims.ndisplay = viewer_src.dims.ndisplay
            viewer_dest.camera.zoom = viewer_src.camera.zoom
            viewer_dest.camera.center = viewer_src.camera.center
    
    # Connect the zoom and center synchronization function to the events of viewer1
    viewer1.camera.events.zoom.connect(lambda event: sync_zoom_center(viewer1, viewer2))
    viewer1.camera.events.center.connect(lambda event: sync_zoom_center(viewer1, viewer2))
    
    @viewer1.mouse_move_callbacks.append
    def mouse_move_callback(viewer, event):
        try:
            # Get the current cursor position
            points1_napari = event.position
            # Transform to the coordinate system of sitk
            points1_sitk = [points1_napari[1], points1_napari[0]]

            ## Transform back the cursor to the coordinate system of the original img1_data (not rotated)
            # Step 1: subtract offset (shift back from new center to old center)
            points1_sitk += img1_data_offset[::-1]  # Reverse (dy, dx) → (dx, dy)
            # Step 2: rotate in the opposite direction to undo rotation
            points1_sitk = rotate_coordinates_angle(points1_sitk, manual_alignment_rotation, img1_data_shape)

            # Transform to the coordinate system of the resized image (the one used for the registration)
            scaled_points1 = tuple(value * scale_percent1_napari for value in points1_sitk)
            # Shift coords if cropping has been done on the first image
            scaled_points1 = shift_tuple_coordinates(scaled_points1, crop_coords1)
            # Shift and rotate coords based on pre-alignment
            scaled_points1 = rotate_coordinates_angle(scaled_points1, -manual_alignment_rotation, manual_alignment_rotation_shape)
            #scaled_points1 = shift_tuple_coordinates(scaled_points1, tuple(np.negative(manual_alignment_rotation_offset)))

            scaled_points1 = shift_tuple_coordinates(scaled_points1, tuple(np.negative(manual_alignment_displacement)))
            
            # Translate the cursor position to the corresponding point in img2_data
            scaled_pointsTR_tmp = outTx_Bspline.TransformPoint(scaled_points1)
            scaled_pointsTR = outTx_Rigid.TransformPoint(scaled_pointsTR_tmp)


            # Shift coords back if cropping has been done on the second image
            scaled_pointsTR = shift_back_tuple_coordinates(scaled_pointsTR, crop_coords2)
            # Transform back to the coordinate system of the full image (the one used for visualization)
            pointsTR_sitk = tuple(value / scale_percent2_napari for value in scaled_pointsTR)
            # Transform back to the coordinate system of napari viewer
            pointsTR_napari = [pointsTR_sitk[1], pointsTR_sitk[0]]
            # Update the translated points in viewer2
            viewer2.layers['Translated Points'].data = pointsTR_napari
            # Calculate the zoom factor between viewer1 and viewer2
            zoom_factor = viewer2.camera.zoom / viewer1.camera.zoom
            # Zoom to the corresponding position in viewer2
            viewer2.camera.center = pointsTR_napari
            viewer2.camera.zoom = zoom_factor * viewer1.camera.zoom
    
        except Exception as e:
            print(f"An error occurred: {e}")
    
    
    # Create an empty layer for translated points in viewer2
    viewer2.add_points([], face_color='red', size=10, name='Translated Points')

    def add_polygons_to_viewer(geometry_list, scale_percent, scale_percent_napari, viewer, viewer_id, edge_color, face_color, opacity, name):
        all_coords = []
        for single_polygon in geometry_list:
            if single_polygon.is_empty:
                continue

            single_polygon = single_polygon.buffer(0)

            if single_polygon.geom_type == 'Polygon':
                polygons = [single_polygon]
            elif single_polygon.geom_type == 'MultiPolygon':
                polygons = list(single_polygon.geoms)
            else:
                continue

            for poly in polygons:
                coords = (np.array(poly.exterior.coords) * scale_percent) / scale_percent_napari
                coords = coords[:, [1, 0]]  # Invert x/y for Napari
                if viewer_id == "viewer1":
                    orig_center = np.array([img1_data_shape[0] / 2, img1_data_shape[1] / 2])  # (y, x)
                    new_center = orig_center - np.array(img1_data_offset) 

                    # Step 1: Translate to origin (center of original image)
                    coords_centered = coords - orig_center

                    # Step 2: Rotate around origin
                    theta = np.deg2rad(manual_alignment_rotation)
                    rot_matrix = np.array([
                        [np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]
                    ])
                    rotated_coords = coords_centered @ rot_matrix.T

                    # Step 3: Translate into rotated image space
                    coords = rotated_coords + new_center



                    




                if coords.shape[0] > 2:
                    all_coords.append(coords)

        if all_coords:
            viewer.add_shapes(
                all_coords,
                shape_type='polygon',
                edge_color=edge_color,
                face_color=face_color,
                opacity=opacity,
                edge_width=10,
                name=name
            )

    if analysis_area_alignment:
        for geometry, scale_percent, scale_percent_napari, viewer, viewer_id in [
            (analysis_area_alignment[panels_alignment[0]], scale_percent1, scale_percent1_napari, viewer1, "viewer1"),
            (analysis_area_alignment[panels_alignment[1]], scale_percent2, scale_percent2_napari, viewer2, "viewer2"),
        ]:
            add_polygons_to_viewer(
                geometry,
                scale_percent,
                scale_percent_napari,
                viewer,
                viewer_id,
                edge_color='green',
                face_color='gray',
                opacity=0.3,
                name="Analysis area"
            )


    if artefacts_empty_alignment:
        for geometry, scale_percent, scale_percent_napari, viewer, viewer_id in [
            (artefacts_empty_alignment[panels_alignment[0]], scale_percent1, scale_percent1_napari, viewer1, "viewer1"),
            (artefacts_empty_alignment[panels_alignment[1]], scale_percent2, scale_percent2_napari, viewer2, "viewer2"),
        ]:
            add_polygons_to_viewer(
                geometry,
                scale_percent,
                scale_percent_napari,
                viewer,
                viewer_id,
                edge_color='transparent',
                face_color='red',
                opacity=0.5,
                name="Artefacts/ Empty area"
            )

    if shapes_data_alignment:
        viewer2.add_shapes(shapes_data_alignment, shape_type='polygon')

    '''
    if shapes_data_alignment and shapes_transformed:
        viewer1.add_shapes(shapes_data_alignment, shape_type='polygon')
        viewer2.add_shapes(shapes_transformed, shape_type='polygon')
    '''
    
    # Display the viewer windows
    viewer1.show()
    viewer2.show()
    napari.run()

    # Compute the maximum indices for each dimension
    max_indices1 = [image_shape1[i] - 1 for i in range(len(image_shape1))]

    ## save manual annotations (clipping area)
    if 'Shapes' in viewer1.layers:
        polygons = []
        for polygon in viewer1.layers['Shapes'].data:
            polygon_sitk = transform_polygon_to_sitk(polygon, img1_data_offset, manual_alignment_rotation, img1_data_shape)
            polygons.append(Polygon(polygon_sitk))
        manual_empty_geometries = gpd.GeoSeries(polygons)

        ## Scale annotations to downscaled image resolution
        manual_empty_array = get_gdf(manual_empty_geometries, scale_percent1_napari, crop_coords1, image_shape1, img1_resize, c="B")
        # First scale with scale percent of the second downscaling (for napari visualization) and second scale with the scale percent of the full resolution image to get back to the resolution of the coordinates in pixel
        manual_empty_alignment[f'{panels_alignment[0]}_in_{panels_alignment[1]}_panel'] = transform_annotation(manual_empty_array, scale_percent1_napari, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, outTx_Rigid, outTx_Bspline, img2_resize)


    return manual_empty_alignment


# Functions that convert coordinates between napari and sitk systems
## reference point (origin) is not the same in both systems
def transform_polygon_to_sitk(polygon,img1_data_offset, manual_alignment_rotation, img1_data_shape):
    try:
        polygon1_sitk = []
        for points1_napari in polygon:
            # Transform to the coordinate system of sitk
            points1_sitk = [points1_napari[1], points1_napari[0]]

            ## Transform back the cursor to the coordinate system of the original img1_data (not rotated)
            # Step 1: subtract offset (shift back from new center to old center)
            points1_sitk += img1_data_offset[::-1]  # Reverse (dy, dx) → (dx, dy)
            # Step 2: rotate in the opposite direction to undo rotation
            points1_sitk = rotate_coordinates_angle(points1_sitk, manual_alignment_rotation, img1_data_shape)

            polygon1_sitk.append(points1_sitk)

        return np.array(polygon1_sitk)

    except Exception as e:
        print(f"An error occurred: {e}")


def merge_manual_empty_with_mask(manual_empty_alignment, mask, output_path, id, resolution_micron):
    panels_tr = [key for key in manual_empty_alignment.keys() if "in" in key]
    manual_empty_alignment_list = []
    for panel in panels_tr:
        manual_empty_alignment_list.append(manual_empty_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))

    ## Filter out empty GeoSeries
    manual_empty_alignment_list = [g for g in manual_empty_alignment_list if not g.empty]
    ## Merge the remaining GeoSeries
    merged_polygons = gpd.GeoSeries(unary_union([g.unary_union for g in manual_empty_alignment_list]))

    # Create a Polygon representing the inner boundary
    inner_boundary = merged_polygons.unary_union

    # Construct the mask by taking the difference between the outer boundary and inner boundary
    if inner_boundary:
        mask = mask.difference(inner_boundary)
    else:
        print("No empty or artefact area.")
        mask = mask

    # Plot the MultiPolygon
    # plot_multipolygon(mask)
    
    ## Convert mask to micro meters
    scaled_mask = scale_multipolygon_coordinates(mask, resolution_micron)
    geojson_dict = mapping(scaled_mask)
    with open(output_path + f"Alignment/merged_annotations/{id}_mask.geojson", "w") as f:
        json.dump(geojson_dict, f, indent=4)

    return mask

def shift_tuple_coordinates(coords, crop_coords):
    """Applies the shift induced by cropping to coordinates as is done to the image."""
    x, y = coords
    shift_x, shift_y = crop_coords  
    shifted_coords = (x - shift_x, y - shift_y)  

    return shifted_coords

def shift_back_tuple_coordinates(coords, crop_coords):
    """Applies the INVERSE of the shift induced by cropping to coordinates to get back to the original image."""
    x, y = coords
    shift_x, shift_y = crop_coords  
    shifted_coords = (x + shift_x, y + shift_y) 

    return shifted_coords
