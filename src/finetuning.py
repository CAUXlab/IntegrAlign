import pickle
import os
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
import napari # type: ignore
import vispy.color
from shapely.geometry import mapping, Polygon
from shapely.ops import unary_union

from src.alignment import remove_params, extract_downscaled_images
from src.alignment import get_cells_coordinates_SPIAT_CellType, get_annotations
from src.alignment import alignment_report

from src.save_images import panels_name_alignment

from src.alignment import merge_annotations, get_gdf, transform_annotation, rotate_coordinates_angle, scale_multipolygon_coordinates, plot_multipolygon
from src.alignment import transform_filter_coordinates, filter_coordinates
from src.alignment import save_tables, plot_rasters




def finetuning(id, meshsize, downscaled_images_path, coordinate_tables, annotations_tables, visualization, scans_path, resolution_micron, metric, pixel_size_raster_micron, alpha_red):
    ## Load and read the pickle file
    with open(downscaled_images_path, "rb") as file:
        downscaled_images = pickle.load(file)
    ## Get the parameters
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

    # iterate over the alignment with the downscaled images and the files to the full resolution .qptiff files for the corresponding alignment
    for i_ms, ((name_alignment, scans_alignment_paths), (_, downscaled_images_id_name_alignment)) in enumerate(zip(scans_alignment_dict.items(), downscaled_images_id.items())):
    # for name_alignment, folder_alignment_paths in panel_alignment_dict.items():
    # for name_alignment, downscaled_images_id_name_alignment in downscaled_images_id.items():
        # downscaled_images_id_name_alignment = downscaled_images_id[name_alignment]
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
            annotations = annotations_tables[index]
            # Get the file path corresponding to id
            csv_file_path = next((os.path.join(coordinate_table, f) for f in os.listdir(coordinate_table) if id in f and f.endswith(".csv") and not f.startswith(".")), None)
            geojson_file_path = next((os.path.join(annotations, f) for f in os.listdir(annotations) if id in f and f.endswith(".geojson") and not f.startswith(".")), None)
            ## Get the coordinates table and annotations
            cell_coordinates, data_frame_cells = get_cells_coordinates_SPIAT_CellType(csv_file_path, panel, cell_coordinates, data_frame_cells, resolution_micron)
            artefacts_empty_alignment, analysis_area_alignment = get_annotations(geojson_file_path, panel, artefacts_empty_alignment, analysis_area_alignment)
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
        sys.exit()  # Exit the script immediately

    ## Transform and merge the annotations
    print("-----")
    print("Merging the annotations...") 
    mask = merge_annotations(id, artefacts_empty_alignment, analysis_area_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, 
                        img_resize_alignment, metric_ms_alignment, metadata_images, panels_all, reference_panel, resolution_micron, output_path)
    
    if visualization != "0":
        mask = merge_manual_empty_with_mask(manual_empty_alignment, mask, output_path, id, resolution_micron)

    ## Transform and filter coordinates
    print("Transform, filer and merge coordinates...")
    merged_cell_coordinates = pd.DataFrame()
    for name_alignment in downscaled_images_id.keys():
        panels_alignment = name_alignment.split("_")
        merged_cell_coordinates = transform_filter_coordinates(metadata_images, cell_coordinates, data_frame_cells, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, mask, resolution_micron, merged_cell_coordinates)
    merged_cell_coordinates = filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, resolution_micron, merged_cell_coordinates)
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
    start_time = time.time()  # Record start time
    cimg, simg1_Rigid, simg2_Rigid, outTx_Rigid = imageRegRigid(img1,img2)
    end_time = time.time()  # Record end time
    execution_time = end_time - start_time  # Calculate execution time
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

    start_time = time.time()  # Record start time
    cimg, simg1, simg2, outTx_Bspline = imageRegBspline(img1, img2, transformDomainMeshSize, spline_order, metric)
    end_time = time.time()  # Record end time
    execution_time = end_time - start_time  # Calculate execution time
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
    # correlation statistique entre les valeurs des pixels de deux images
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





def getLabels(tif_tags, nb_channels):
    substr = "<ScanColorTable-k>"
    start = 0
    if tif_tags['ImageDescription'].find(substr, start) != -1:  
        strings = []
        while True:
            start = tif_tags['ImageDescription'].find(substr, start)
            if start == -1: # when '<ScanColorTable-k>' is not found
                break
            string = tif_tags['ImageDescription'][start+18:start+26]
            strings.append(string)
            start += 1
        marquages = []
        wl = []
        for s in strings:
            if s.startswith('DAPI'):
                marquages.append(s[0:4])
                wl.append('450')
            if s.startswith('Opal'):
                marquages.append(s)
                wl.append(s[5:])
                
        dictionary = {key: value for key, value in enumerate(wl)}
        # change to detailled list
        channel_list = [f'{value} (channel {key})' for key, value in enumerate(marquages)]
    else:
        dictionary = {'DAPI': '450'}
        channel_list = ['450 (channel DAPI)']

    return channel_list, dictionary

def load_QPTIFF_high_res(path):
    # Load qptiff file
    with TiffFile(path) as tif:
        # Get tags from DAPI channel in .pages
        tif_tags = {}
        for tag in tif.pages[0].tags.values():
            name, value = tag.name, tag.value
            tif_tags[name] = value
        # get nb of channels
        last_indexes = 10
        im_sizes = []
        for im in tif.pages[-last_indexes:]:  
            im_sizes.append(len(im.asarray()))
        from collections import Counter
        count_dict = Counter(im_sizes)
        nb_channels = max(count_dict.values())
        # Load array with every channel
        index_first_channel_second_downscaling = nb_channels*2+1 # 8 channels, 1 low res RGB, 8 channels (first downscaling of the pyramid representation) to get to the first channel of the second downscaling
        index_last_channel_second_downscaling = nb_channels*3+1 # 8 channels, 1 low res RGB, 8 channels (first downscaling of the pyramid representation), 8 channels (second downscaling of the pyramid representation) to get to the last channel of the second downscaling
        img_data = np.stack([tif.pages[i].asarray() for i in range(index_first_channel_second_downscaling, index_last_channel_second_downscaling)], axis=0)
        # img_data = tif.series[0].asarray()
        
    channel_list, channel_name_dictionary = getLabels(tif_tags, nb_channels)
    
    return img_data, tif_tags, channel_list, channel_name_dictionary


def load_high_res_imgs(folder_path1, folder_path2, id, panels):
    # Define path of the qptiff
    path1 = folder_path1 + id + f"_panel_{panels[0]}.unmixed.qptiff"
    path2 = folder_path2 + id + f"_panel_{panels[1]}.unmixed.qptiff"

    img1_data, tif_tags1, channel_list1, channel_name_dictionary1 = load_QPTIFF_high_res(path1)
    img2_data, tif_tags2, channel_list2, channel_name_dictionary2 = load_QPTIFF_high_res(path2)
    
    

    '''
    # Saving full resolution images
    np.save(dir_res + f"{id}/{folder_name_alignment}/results/img1_data_{id}_panel_{panels[0]}.npy", img1_data)
    np.save(dir_res + f"{id}/{folder_name_alignment}/results/img2_data_{id}_panel_{panels[0]}.npy", img2_data)
    '''

    return img1_data, tif_tags1, channel_list1, channel_name_dictionary1, img2_data, tif_tags2, channel_list2, channel_name_dictionary2



def mirrored_cursor_visu(scan_path1, scan_path2, id, metadata_images, img1_resize, img2_resize, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, manual_empty_alignment, analysis_area_alignment = None, artefacts_empty_alignment = None, shapes_data_alignment=None, shapes_transformed=None, full_res=None):
    print(f"Loading .qptiff images...")
    panels = name_alignment.split('_')
    img1_data, tif_tags1, channel_list1, channel_name_dictionary1, img2_data, tif_tags2, channel_list2, channel_name_dictionary2 = load_high_res_imgs(scan_path1, scan_path2, id, panels)
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
    

    viewer1 = napari.Viewer()
    viewer1.title = panels[0]
    for i in range(len(channel_list1)):
        wavelength = int(list(channel_name_dictionary1.values())[i])
        rgb_values = wavelength_to_rgb(wavelength)
        colorMap = vispy.color.Colormap([[0.0, 0.0, 0.0], rgb_values])
        if full_res:
            rlayer = viewer1.add_image([img1_data[i],
                                       img1_data[i][::4, ::4],
                                       img1_data[i][::8, ::8]],
                                       name=channel_list1[i], contrast_limits=[0, 255])
        else:
            rlayer = viewer1.add_image(img1_data[i],
                           name=channel_list1[i], contrast_limits=[0, 255])
        rlayer.blending = 'additive'
        rlayer.colormap = colorMap
        
    viewer2 = napari.Viewer()
    viewer2.title = panels[1]
    for i in range(len(channel_list2)):
        wavelength = int(list(channel_name_dictionary2.values())[i])
        rgb_values = wavelength_to_rgb(wavelength)
        colorMap = vispy.color.Colormap([[0.0, 0.0, 0.0], rgb_values])
        if full_res:
            rlayer = viewer2.add_image([img2_data[i],
                                       img2_data[i][::4, ::4],
                                       img2_data[i][::8, ::8]],
                                       name=channel_list2[i], contrast_limits=[0, 255])
        else:
            rlayer = viewer2.add_image(img2_data[i],
                           name=channel_list2[i], contrast_limits=[0, 255])
        
        rlayer.blending = 'additive'
        rlayer.colormap = colorMap
    
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
            # Transform to the coordinate system of the resized image (the one used for the registration)
            scaled_points1 = tuple(value * scale_percent1_napari for value in points1_sitk)
            # Shift coords if cropping has been done on the first image
            scaled_points1 = shift_tuple_coordinates(scaled_points1, crop_coords1)
            # Shift and rotate coords based on pre-alignment
            scaled_points1 = rotate_coordinates_angle(scaled_points1, -manual_alignment_rotation, image_shape_manual_alignment)
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

    if analysis_area_alignment:
        geometry = analysis_area_alignment[panels_alignment[0]]
        # Loop through all geometries in the GeoSeries
        for single_polygon in geometry:
            # Ensure it's a valid Polygon
            if single_polygon.geom_type == 'Polygon':
                # Extract exterior coordinates
                single_polygon = single_polygon.buffer(0)
                if single_polygon.geom_type == 'Polygon':
                    # Still a Polygon after buffering
                    coords = (np.array(single_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer1.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")
                elif single_polygon.geom_type == 'MultiPolygon':
                    # It turned into a MultiPolygon after buffering
                    for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                        coords = (np.array(sub_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                        # Invert x and y coordinates (napari coordinate system)
                        coords = coords[:, [1, 0]]
                        # Add the sub-polygon to the viewer
                        if coords.shape[0] > 2:
                            viewer1.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")
            
            # Handling for MultiPolygon geometry (if not already handled)
            elif single_polygon.geom_type == 'MultiPolygon':
                for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                    coords = (np.array(sub_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the sub-polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer1.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")


        geometry = analysis_area_alignment[panels_alignment[1]]
        # Loop through all geometries in the GeoSeries
        for single_polygon in geometry:
            # Ensure it's a valid Polygon
            if single_polygon.geom_type == 'Polygon':
                # Extract exterior coordinates
                single_polygon = single_polygon.buffer(0)
                if single_polygon.geom_type == 'Polygon':
                    # Still a Polygon after buffering
                    coords = (np.array(single_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer2.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")
                elif single_polygon.geom_type == 'MultiPolygon':
                    # It turned into a MultiPolygon after buffering
                    for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                        coords = (np.array(sub_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                        # Invert x and y coordinates (napari coordinate system)
                        coords = coords[:, [1, 0]]
                        # Add the sub-polygon to the viewer
                        if coords.shape[0] > 2:
                            viewer2.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")
            
            # Handling for MultiPolygon geometry (if not already handled)
            elif single_polygon.geom_type == 'MultiPolygon':
                for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                    coords = (np.array(sub_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the sub-polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer2.add_shapes([coords], shape_type='polygon', edge_color='red', face_color='gray', opacity=0.3, edge_width=10, name="Analysis area")

    '''
    if artefacts_empty_alignment:
        geometry = artefacts_empty_alignment[panels_alignment[0]]
        # Loop through all geometries in the GeoSeries
        for single_polygon in geometry:
            # Ensure it's a valid Polygon
            if single_polygon.geom_type == 'Polygon':
                # Extract exterior coordinates
                single_polygon = single_polygon.buffer(0)
                if single_polygon.geom_type == 'Polygon':
                    # Still a Polygon after buffering
                    coords = (np.array(single_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer1.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")
                elif single_polygon.geom_type == 'MultiPolygon':
                    # It turned into a MultiPolygon after buffering
                    for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                        coords = (np.array(sub_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                        # Invert x and y coordinates (napari coordinate system)
                        coords = coords[:, [1, 0]]
                        # Add the sub-polygon to the viewer
                        if coords.shape[0] > 2:
                            viewer1.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")
            
            # Handling for MultiPolygon geometry (if not already handled)
            elif single_polygon.geom_type == 'MultiPolygon':
                for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                    coords = (np.array(sub_polygon.exterior.coords) * scale_percent1)/scale_percent1_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the sub-polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer1.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")


        geometry = artefacts_empty_alignment[panels_alignment[1]]
        # Loop through all geometries in the GeoSeries
        for single_polygon in geometry:
            # Ensure it's a valid Polygon
            if single_polygon.geom_type == 'Polygon':
                # Extract exterior coordinates
                single_polygon = single_polygon.buffer(0)
                if single_polygon.geom_type == 'Polygon':
                    # Still a Polygon after buffering
                    coords = (np.array(single_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer2.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")
                elif single_polygon.geom_type == 'MultiPolygon':
                    # It turned into a MultiPolygon after buffering
                    for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                        coords = (np.array(sub_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                        # Invert x and y coordinates (napari coordinate system)
                        coords = coords[:, [1, 0]]
                        # Add the sub-polygon to the viewer
                        if coords.shape[0] > 2:
                            viewer2.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")
            
            # Handling for MultiPolygon geometry (if not already handled)
            elif single_polygon.geom_type == 'MultiPolygon':
                for sub_polygon in single_polygon.geoms:  # Use .geoms to iterate over the individual polygons
                    coords = (np.array(sub_polygon.exterior.coords) * scale_percent2)/scale_percent2_napari
                    # Invert x and y coordinates (napari coordinate system)
                    coords = coords[:, [1, 0]]
                    # Add the sub-polygon to the viewer
                    if coords.shape[0] > 2:
                        viewer2.add_shapes([coords], shape_type='polygon', edge_color='white', face_color='white', opacity=0.6, edge_width=1, name="Artefacts/ Empty area")
    '''


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
            polygon_sitk = transform_polygon_to_sitk(polygon)
            polygons.append(Polygon(polygon_sitk))
        manual_empty_geometries = gpd.GeoSeries(polygons)

        ## Scale annotations to downscaled image resolution
        manual_empty_array = get_gdf(manual_empty_geometries, scale_percent1_napari, crop_coords1, image_shape1, img1_resize, c="B")
        # First scale with scale percent of the second downscaling (for napari visualization) and second scale with the scale percent of the full resolution image to get back to the resolution of the coordinates in pixel
        manual_empty_alignment[f'{panels_alignment[0]}_in_{panels_alignment[1]}_panel'] = transform_annotation(manual_empty_array, scale_percent1_napari, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, image_shape_manual_alignment, outTx_Rigid, outTx_Bspline, img2_resize)


    return manual_empty_alignment


# Functions that convert coordinates between napari and sitk systems
## reference point (origin) is not the same in both systems
def transform_polygon_to_sitk(polygon):
    try:
        polygon1_sitk = []
        for points1_napari in polygon:
            # Transform to the coordinate system of sitk
            points1_sitk = [points1_napari[1], points1_napari[0]]
            polygon1_sitk.append(points1_sitk)

        return np.array(polygon1_sitk)

    except Exception as e:
        print(f"An error occurred: {e}")



"""
This function is based on code originally written in Fortran by Dan Bruton,
and now available here: http://www.midnightkite.com/color.html.
There are a number of implementations floating around in various languages (Pascal, R etc.),
often with the logic tweaked slightly.
"""
def wavelength_to_rgb(nm):

    gamma = 0.8
    max_intensity = 255
    factor = 0

    rgb = {"R": 0, "G": 0, "B": 0}

    if 380 <= nm <= 439:
        rgb["R"] = -(nm - 440) / (440 - 380)
        rgb["G"] = 0.0
        rgb["B"] = 1.0
    elif 440 <= nm <= 489:
        rgb["R"] = 0.0
        rgb["G"] = (nm - 440) / (490 - 440)
        rgb["B"] = 1.0
    elif 490 <= nm <= 509:
        rgb["R"] = 0.0
        rgb["G"] = 1.0
        rgb["B"] = -(nm - 510) / (510 - 490)
    elif 510 <= nm <= 579:
        rgb["R"] = (nm - 510) / (580 - 510)
        rgb["G"] = 1.0
        rgb["B"] = 0.0
    elif 580 <= nm <= 644:
        rgb["R"] = 1.0
        rgb["G"] = -(nm - 645) / (645 - 580)
        rgb["B"] = 0.0
    elif 645 <= nm <= 780:
        rgb["R"] = 1.0
        rgb["G"] = 0.0
        rgb["B"] = 0.0

    if 380 <= nm <= 419:
        factor = 0.3 + 0.7 * (nm - 380) / (420 - 380)
    elif 420 <= nm <= 700:
        factor = 1.0
    elif 701 <= nm <= 780:
        factor = 0.3 + 0.7 * (780 - nm) / (780 - 700)

    if rgb["R"] > 0:
        rgb["R"] = int(max_intensity * ((rgb["R"] * factor) ** gamma))
    else:
        rgb["R"] = 0

    if rgb["G"] > 0:
        rgb["G"] = int(max_intensity * ((rgb["G"] * factor) ** gamma))
    else:
        rgb["G"] = 0

    if rgb["B"] > 0:
        rgb["B"] = int(max_intensity * ((rgb["B"] * factor) ** gamma))
    else:
        rgb["B"] = 0

    return (rgb["R"]/255, rgb["G"]/255, rgb["B"]/255)



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
    shifted_coords = (x - shift_x, y - shift_y)  # Return a tuple instead of an array

    return shifted_coords

def shift_back_tuple_coordinates(coords, crop_coords):
    """Applies the INVERSE of the shift induced by cropping to coordinates to get back to the original image."""
    x, y = coords
    shift_x, shift_y = crop_coords  
    shifted_coords = (x + shift_x, y + shift_y)  # Return a tuple instead of an array

    return shifted_coords
