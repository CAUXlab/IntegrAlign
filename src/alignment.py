import pickle
import SimpleITK as sitk # type: ignore
import matplotlib.pyplot as plt
import os
import numpy as np
import time
from tqdm import tqdm # type: ignore
import matplotlib.patches as mpatches
from io import BytesIO
import base64
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.stats import pearsonr
from rasterio.transform import from_origin # type: ignore
from pathlib import Path
from shapely.geometry import mapping, Polygon, MultiPolygon, Point
import geopandas as gpd
from shapely.ops import unary_union
import json
import shapely.vectorized


def alignment(downscaled_images_path, coordinate_tables, annotations_tables, resolution_micron, number_ms, metric, pixel_size_raster_micron, alpha_red):
    ## Load and read the pickle file
    with open(downscaled_images_path, "rb") as file:
        downscaled_images = pickle.load(file)

    ## Get the parameters
    common_ids = downscaled_images["params"].get("common_ids")
    panels_all = downscaled_images["params"].get("panels")
    # print(panels_all)
    output_path = downscaled_images["params"].get("output_path")
    turn_img_id_dict = downscaled_images["params"].get("turn_img_id_dict")
    
    ## Remove "params" from downscaled_images
    downscaled_images = remove_params(downscaled_images)
    ## Create alignment folder to store alignment reports
    Path(output_path+'Alignment/reports').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/merged_annotations').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/merged_tables').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/plots').mkdir(parents=True, exist_ok=True)
    
    
    ## Get reference panel
    reference_panel = panels_all[1]

    for id, downscaled_images_id in downscaled_images.items():
        print("---------------")
        print(f'Patient: {id}')
        metadata_images = {}
        cell_coordinates = {}
        data_frame_cells = {}
        metric_ms_alignment = {}
        img_resize_alignment = {}
        outTx_Rigid_alignment = {}
        outTx_Bspline_alignment = {}
        outTx_Bspline_inverse_alignment = {}
        artefacts_empty_alignment = {}
        analysis_area_alignment = {}
        for name_alignment, downscaled_images_id_name_alignment in downscaled_images_id.items():
            print("-----")
            print(f"Alignment {name_alignment}...")
            ## Load downscaled images
            print("Loading downscaled images...")
            panels_alignment = name_alignment.split("_")
            metadata_images, img1, img2, img1_resize, img2_resize = extract_downscaled_images(downscaled_images_id_name_alignment, panels_alignment, 
                                                                                              name_alignment, turn_img_id_dict, metadata_images, id)
            ## Alignment
            spline_order = 3
            (outTx_Bspline_dict, simg2_dict, execution_time_dict, 
             metric_values_dict, outTx_Rigid, simg1_Rigid, 
             simg2_Rigid, nda_Rigid) = alignment_2panels(name_alignment, id, img1, img2, number_ms, metric, spline_order)
            ## Load coordinate table and annotations of each panel
            print(f'Loading tables and annotations...')
            for panel in panels_alignment:
                # Get the panel path
                index = panels_all.index(panel)
                coordinate_table = coordinate_tables[index]
                annotations = annotations_tables[index]
                # Get the file path corresponding to id
                csv_file_path = next((os.path.join(coordinate_table, f) for f in os.listdir(coordinate_table) if id in f and f.endswith(".csv")), None)
                geojson_file_path = next((os.path.join(annotations, f) for f in os.listdir(annotations) if id in f and f.endswith(".geojson")), None)
                ## Get the coordinates table and annotations
                cell_coordinates, data_frame_cells = get_cells_coordinates_SPIAT_CellType(csv_file_path, panel, cell_coordinates, data_frame_cells, resolution_micron)
                artefacts_empty_alignment, analysis_area_alignment = get_annotations(geojson_file_path, panel, artefacts_empty_alignment, analysis_area_alignment)
            ## Create the alignment report
            print("Creating the alignment report...")
            outTx_Rigid_alignment, outTx_Bspline_alignment, img_resize_alignment, metric_ms_alignment = alignment_report(id, name_alignment, panels_alignment, cell_coordinates, metadata_images, output_path, 
                                                                                                    alpha_red, img1_resize, img2_resize, simg2_Rigid, outTx_Rigid, outTx_Bspline_dict, simg2_dict, 
                                                                                                    execution_time_dict, metric_values_dict, spline_order, resolution_micron, pixel_size_raster_micron, 
                                                                                                    metric_ms_alignment, img_resize_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment)
        ## Transform and merge the annotations
        print("-----")
        print("Merging the annotations...")    
        mask = merge_annotations(id, artefacts_empty_alignment, analysis_area_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, 
                          img_resize_alignment, metric_ms_alignment, metadata_images, panels_all, reference_panel, resolution_micron, output_path)
        
        
        ## Transform and filter coordinates
        print("Transform, filer and merge coordinates...")
        merged_cell_coordinates = pd.DataFrame()
        for name_alignment in downscaled_images_id.keys():
            panels_alignment = name_alignment.split("_")
            merged_cell_coordinates = transform_filter_coordinates(metadata_images, cell_coordinates, data_frame_cells, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, mask, resolution_micron, merged_cell_coordinates)
        merged_cell_coordinates = filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, merged_cell_coordinates)
        ## Save new coordinates table
        save_tables(merged_cell_coordinates, output_path, id)
        ## Plot cell coordinates before and after alignment
        plot_rasters(data_frame_cells, merged_cell_coordinates, cell_coordinates, pixel_size_raster_micron, output_path, panels_all, id, resolution_micron)












def remove_params(d):
    ''' Recursive function to remove "params" from the dictionnary downscaled_images. '''
    if isinstance(d, dict):
        # Remove "params" key and recursively clean nested dictionaries
        return {k: remove_params(v) for k, v in d.items() if k != "params"}
    elif isinstance(d, list):
        # If it's a list, recursively clean each item
        return [remove_params(item) for item in d]
    return d

def extract_downscaled_images(downscaled_images_id_name_alignment, panels, name_alignment, turn_img_id_dict, metadata_images, id):
    ## Load images
    img1 = downscaled_images_id_name_alignment.get("img1")
    img2 = downscaled_images_id_name_alignment.get("img2")
    img1_resize = downscaled_images_id_name_alignment.get("img1_resize")
    img2_resize = downscaled_images_id_name_alignment.get("img2_resize")
    panel1, panel2 = panels
    ## Get metadata
    metadata_images[name_alignment] = {
        f"tif_tags_{panel1}": downscaled_images_id_name_alignment.get("tif_tags1"),
        f"channel_list_{panel1}": downscaled_images_id_name_alignment.get("channel_list1"),
        f"channel_name_dictionary_{panel1}": downscaled_images_id_name_alignment.get("channel_name_dictionary1"),
        f"scale_percent_{panel1}": downscaled_images_id_name_alignment.get("scale_percent1"),
        f"image_shape_{panel1}": img1_resize.shape,
        f"turn_img_{panel1}": turn_img_id_dict[id + "_" + panel1],
        f"tif_tags_{panel2}": downscaled_images_id_name_alignment.get("tif_tags2"),
        f"channel_list_{panel2}": downscaled_images_id_name_alignment.get("channel_list2"),
        f"channel_name_dictionary_{panel2}": downscaled_images_id_name_alignment.get("channel_name_dictionary2"),
        f"scale_percent_{panel2}": downscaled_images_id_name_alignment.get("scale_percent2"),
        f"image_shape_{panel2}": img2_resize.shape,
        f"turn_img_{panel2}": turn_img_id_dict[id + "_" + panel2]
    }

    return metadata_images, img1, img2, img1_resize, img2_resize


def alignment_2panels(name_alignment, id, img1, img2, number_ms, metric, spline_order):
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
    deformation_degree_all = [i for i in range(1, number_ms + 1)]
    # print("Deformation degree:",deformation_degree_all)
    
    # Alignment on multiple mesh sizes
    for deformation_degree in tqdm(deformation_degree_all, desc="Alignment on different mesh size", unit="ms"):
        transformDomainMeshSize = [deformation_degree, deformation_degree]
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
    
        execution_time_dict[f"{metric}_{deformation_degree}"] = execution_time
        metric_values_dict[f"{metric}_{deformation_degree}"] = metric_values_list
        simg2_dict[f"{metric}_{deformation_degree}"] = simg2
        outTx_Bspline_dict[f"{metric}_{deformation_degree}"] = outTx_Bspline
    
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
        numberOfIterations=5000,
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


def get_cells_coordinates_SPIAT_CellType(csv_file_path, panel, cell_coordinates, data_frame_cells, resolution_micron):
    df = pd.read_csv(csv_file_path)
    
    new_df = pd.DataFrame()
    # Compute the center coordinates and store them in the new DataFrame
    ## Convert pixel coordinates to micrometers
    new_df['x (µm)'] = df['x']
    new_df['y (µm)'] = df['y']

    new_df['x'] = df['x'] * resolution_micron
    new_df['y'] = df['y'] * resolution_micron
        
    # new_df['x'] = ((df['XMin_pixel'] + df['XMax_pixel']) / 2)
    # new_df['y'] = ((df['YMin_pixel'] + df['YMax_pixel']) / 2)
    
    new_df['Phenotype'] = df['Phenotype']
    new_df['Cell_type'] = df['Cell_type']
    new_df['Object Id'] = df['ID']
    new_df['Classifier Label'] = df['Tissue_Category']

    '''
    print(f"\nTotal number of cells: {len(new_df)}")
    print("--")
    cell_type_counts = df['Cell_type'].value_counts()
    for cell_type, count in cell_type_counts.items():
        print(f"{cell_type}: {count} cells")
        # cell_coordinates[f'{panel}_panel_{cell_type}'] = new_df[new_df['Cell_type'] == cell_type][['x', 'y']].to_numpy()
    '''
    
    cell_coordinates[f'{panel}_panel_DAPI'] = new_df[['x', 'y']].to_numpy()
    data_frame_cells[f'{panel}_panel_df'] = new_df

    return cell_coordinates, data_frame_cells


def get_annotations(geojson_file_path, panel, artefacts_empty_alignment, analysis_area_alignment):
    ## Read annotations file
    gdf = gpd.read_file(geojson_file_path)
    def safe_json_load(x):
        try:
            return json.loads(x)
        except json.JSONDecodeError:
            return {}  # Return an empty dictionary for invalid rows
    gdf['classification'] = gdf['classification'].apply(safe_json_load)
    # print(gdf['classification'].apply(lambda x: x.get('name')).unique())
    ## Convert LineString to Polygon
    gdf_poly = gdf['geometry'].apply(lambda geom: Polygon(geom))
    '''
    fig, ax = plt.subplots(figsize=(6, 6)) 
    gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    '''

    ## Check for the annotation types (non-case sensitive)
    Annotations_classification = ['Artefacts', 'Manual_Artefacts', 'Empty']
    artefacts_empty_gdf = gdf[gdf['classification'].apply(lambda x: x.get('name').lower() in [item.lower() for item in Annotations_classification])]
    # print(artefacts_empty_gdf['classification'].apply(lambda x: x.get('name')).unique())
    ## Convert LineString to Polygon
    artefacts_empty_gdf_poly = artefacts_empty_gdf['geometry'].apply(lambda geom: Polygon(geom))
    '''
    fig, ax = plt.subplots(figsize=(6, 6)) 
    artefacts_empty_gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    plt.show()
    '''
    artefacts_empty_alignment[panel] = artefacts_empty_gdf_poly

    ## Get the analysis area annotation
    analysis_area_gdf = gdf[gdf['classification'].apply(lambda x: x.get('name') in ['Analysis_area'])]
    analysis_area_gdf['classification'].apply(lambda x: x.get('name')).unique()
    # Convert LineString to Polygon
    analysis_area_gdf_poly = analysis_area_gdf['geometry'].apply(lambda geom: Polygon(geom))
    '''
    fig, ax = plt.subplots(figsize=(6, 6)) 
    analysis_area_gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    plt.show()
    '''
    analysis_area_alignment[panel] = analysis_area_gdf_poly

    return artefacts_empty_alignment, analysis_area_alignment




def alignment_report(id, name_alignment, panels, cell_coordinates, metadata_images, output_path, alpha_red, img1_resize, img2_resize, simg2_Rigid, outTx_Rigid, outTx_Bspline_dict, simg2_dict, execution_time_dict, metric_values_dict, spline_order, resolution_micron, pixel_size_raster_micron, metric_ms_alignment, img_resize_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment):
    ## Plot rigid alignment
    plt.figure(figsize=(5, 2), dpi=300)
    plt.imshow(img1_resize, cmap="Blues")
    plt.imshow(simg2_Rigid, cmap="Reds", alpha=alpha_red)

    blue_patch = mpatches.Patch(color='lightblue', label='panel '+panels[0])
    red_patch = mpatches.Patch(color='tomato', label='panel '+panels[1])
    plt.legend(handles=[blue_patch, red_patch], loc='upper right', fontsize=4)
    plt.title("RIGID", fontsize=6)
    plt.axis('off')

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_visual_alignment_RIGID = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()

    ## Plot mesh size alignment
    num_images = len(simg2_dict)
    num_cols = num_images  
    num_rows = 1

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))

    if num_images > 1:
        axes = axes.flatten()
    else:
        axes = [axes]  

    for ax, (key, simg2_Bspline) in zip(axes, simg2_dict.items()):
        # Display simg1_Rigid with 'Blues' colormap
        ax.imshow(img1_resize, cmap="Blues")
        # Overlay simg2_Rigid with 'Reds' colormap and alpha=0.2
        ax.imshow(sitk.GetArrayFromImage(simg2_Bspline), cmap="Reds", alpha=alpha_red)
        # Add legend color
        blue_patch = mpatches.Patch(color='lightblue', label='panel '+panels[0])
        red_patch = mpatches.Patch(color='tomato', label='panel '+panels[1])
        ax.legend(handles=[blue_patch, red_patch], loc='upper right', fontsize=9)
        # Set title from the key
        ax.set_title(key, fontsize=10)
        ax.axis('off')

    for i in range(num_images, num_rows * num_cols):
        axes[i].axis('off')

    plt.tight_layout()
    # plt.savefig(dir_res + f"{id}/{folder_name_alignment}/simg2_ALLms.png", dpi=300)

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_visual_alignment = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()


    ## Plot execution time and pixel wise correlation at convergence
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    x_C = []
    y_C = []
    for key, value in execution_time_dict.items():
        if key.startswith("Correlation"):
            if key.endswith("Rigid"):
                x_C.append("RIGID")
                y_C.append(value)
            else:
                x_C.append(int(key.split('_')[-1]))
                y_C.append(value)
    ax1.plot(x_C, y_C, marker='o', markersize = 4, linestyle='-', color='blue')

    ax1.set_xlabel('Mesh Size', size=16)
    ax1.set_ylabel('Execution Time (seconds)', size=16)


    x_C = []
    y_C = []
    for key, value in metric_values_dict.items():  
        if key.startswith("Rigid"):
            x_C.append("RIGID")
            y_C.append(-value[-1])
        if key.startswith("Correlation"):
            x_C.append(key.split('_')[-1])
            y_C.append(-value[-1])
    ax2.plot(x_C, y_C, marker='o', markersize = 3, linestyle='-', color='blue')

    ax2.set_ylim(0, 1)
    ax2.set_xlabel('Mesh Size', size=16)
    ax2.set_ylabel('- Pixel-wise correlation \nat convergence', size=16)

    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_execTime_correlAtConvergence = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()

    ## Plot jacobian
    # Add optional treshold for better visualisation of the volumetric changes (determinant of the jacobian)
    # This treshold should be representative of the tissue that should not be altered
    # If this treshold isn't define it will plot the image with gray gradient instead of binary pixel values (only white and black)
    thresh_pixelSIMG2 = 60.0  # Example threshold value
    
    simg2 = sitk.GetArrayFromImage(simg2_dict[list(outTx_Bspline_dict.keys())[0]])
    simg2_tr = np.clip((simg2 >= thresh_pixelSIMG2) * 255, 0, 255)
    plt.imshow(simg2_tr, cmap='gray')
    plt.axis('off')

    fig, axs = plt.subplots(2, len(outTx_Bspline_dict), figsize=(len(outTx_Bspline_dict) * 6, 12))

    for i, (metric_ms, outTx_Bspline) in enumerate(outTx_Bspline_dict.items()):
        simg2 = sitk.GetArrayFromImage(simg2_dict[metric_ms])
        transformDomainMeshSize = [int(metric_ms.split('_')[1]), int(metric_ms.split('_')[1])]

        array_contour, xx, yy, u, v, jacobian_det_np_arr, array = grid_deform_detJacobian(simg2, transformDomainMeshSize, spline_order, outTx_Bspline)

        # Plot moving image and grid
        axs[0, i].imshow(simg2, cmap='Reds')
        axs[0, i].imshow(array_contour, alpha=0.3, cmap='Reds_r')
        axs[0, i].quiver(xx, yy, -u, -v, color='black', units='xy', angles='xy', scale_units='xy', scale=1, width=10)
        axs[0, i].scatter(xx, yy, s=1)
        axs[0, i].axis('off')

        # Plot Jacobian determinant map
        colors = [
            [0.0, '#2390ff'],  # Red
            [0.25, '#09537d'],
            [0.5, '#000000'],  # Black
            [0.75, '#ca4f04'],
            [1.0, '#ff1b37'],  # Blue
        ]
        custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
        jacobian_det_np_arr = np.nan_to_num(jacobian_det_np_arr, nan=1.0, posinf=1.5, neginf=0.5)
        im = axs[1, i].imshow(jacobian_det_np_arr, cmap=custom_cmap, vmin=0.5, vmax=1.5)
        cbar = plt.colorbar(im, ax=axs[1, i], orientation='vertical', label='det(Jac)')
        axs[1, i].imshow(array_contour, alpha=0.2, cmap='gray')
        if thresh_pixelSIMG2:
            simg2_tr = np.clip((simg2 >= thresh_pixelSIMG2) * 255, 0, 255)
            axs[1, i].imshow(simg2_tr, cmap='gray', alpha=0.4)
        else:
            axs[1, i].imshow(simg2, cmap='gray', alpha=0.4)
        axs[1, i].axis('off')

    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_deformation = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()

    
    pixel_size_raster_fullres = pixel_size_raster_micron * resolution_micron
    pixel_size_raster_fullres

    scale_percent1 = metadata_images[name_alignment][f'scale_percent_{panels[0]}']
    scale_percent2 = metadata_images[name_alignment][f'scale_percent_{panels[1]}']
    image_shape1 = metadata_images[name_alignment][f'image_shape_{panels[0]}']
    image_shape2 = metadata_images[name_alignment][f'image_shape_{panels[1]}']
    operation1 = metadata_images[name_alignment][f'turn_img_{panels[0]}']
    operation2 = metadata_images[name_alignment][f'turn_img_{panels[1]}']
    coords1 = cell_coordinates[f'{panels[0]}_panel_DAPI']
    coords2 = cell_coordinates[f'{panels[1]}_panel_DAPI']
    label1 = f'DAPI from panel {panels[0]}'
    label2 = f'DAPI {panels[1]}'
    plot_correlation, correlations = correlation_rasters_report(outTx_Bspline_dict, outTx_Rigid, scale_percent1, scale_percent2, coords1, coords2, label1, label2, pixel_size_raster_fullres, image_shape1, image_shape2, operation1, operation2)


    best_ms = max(correlations, key=correlations.get)

    num_images = len(simg2_dict)
    num_cols = num_images  
    num_rows = 1  

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))

    if num_images > 1:
        axes = axes.flatten()
    else:
        axes = [axes]  

    for ax, (key, simg2_Bspline) in zip(axes, simg2_dict.items()):
        # Display simg1_Rigid with 'Blues' colormap
        ax.imshow(img1_resize, cmap="Blues")
        # Overlay simg2_Rigid with 'Reds' colormap and alpha=0.2
        ax.imshow(sitk.GetArrayFromImage(simg2_Bspline), cmap="Reds", alpha=alpha_red)
        # Add legend color
        blue_patch = mpatches.Patch(color='lightblue', label='panel '+panels[0])
        red_patch = mpatches.Patch(color='tomato', label='panel '+panels[1])
        ax.legend(handles=[blue_patch, red_patch], loc='upper right', fontsize=9)
        # Set title from the key
        if key == best_ms:
            ax.set_title(key, fontsize=14, color='red', fontweight='bold')
        else:
            ax.set_title(key, fontsize=10)
        ax.axis('off')

    for i in range(num_images, num_rows * num_cols):
        axes[i].axis('off')

    plt.tight_layout()
    # plt.savefig(dir_res + f"{id}/{folder_name_alignment}/simg2_ALLms.png", dpi=300)

    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_visual_alignment = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()

    with open(output_path + f"Alignment/reports/{id}_{name_alignment}_report_alignment.html", 'w') as file:
        file.write("""
        <html>
        <head>
            <title>Validation of the alignment</title>

        </head>
        <body>
        
            <h1>Validation of the alignment</h1>
        """)
        
        file.write(f"""
            <div style="text-align: left;">
                <img src="data:image/png;base64,{plot_execTime_correlAtConvergence}" alt="Correlation at convergence" style="max-width:60%">
        """)
        
        file.write("""
            <h2>Visual alignment</h2>
            <p>Check the alignment of the downsampled images and the deformations for different mesh sizes</p>
        """)
            
        file.write(f"""
            <div>
                <img src="data:image/png;base64,{plot_visual_alignment_RIGID}" alt="Visual alignment plot RIGID" style="max-width:80%; float: left;">
                <img src="data:image/png;base64,{plot_visual_alignment}" alt="Visual alignment plot" style="max-width:300%; float: left;">
            </div>
        """)
        file.write(f"""
            <div>
                <img src="data:image/png;base64,{plot_deformation}" alt="Visual deformation" style="max-width:300%; float: left;">
            </div>
        """)
        file.write(f"""
            <div>
                <img src="data:image/png;base64,{plot_correlation}" alt="Visual correlation" style="max-width:300%; float: left;">
            </div>
        """)

        
        file.write("""
        </body>
        </html>
        """)



    # Chose mesh size based on the report
    metric_ms = best_ms
    metric_ms_alignment[panels[0]] = metric_ms

    outTx_Rigid_alignment[panels[0]] = outTx_Rigid
    outTx_Bspline_alignment[panels[0]] = outTx_Bspline_dict[metric_ms]
    img_resize_alignment[panels[0]] = img1_resize
    img_resize_alignment[panels[1]] = img2_resize

    return outTx_Rigid_alignment, outTx_Bspline_alignment, img_resize_alignment, metric_ms_alignment



def grid_deform_detJacobian(simg2, transformDomainMeshSize, spline_order, outTx_Bspline):
    ## Visualization of the b-spline deformation
    # Create grid to visualize the deformation in the image space
    grid_spacing = simg2.shape[1]/20, simg2.shape[0]/20
    grid = sitk.GridSource(
        outputPixelType=sitk.sitkFloat32,
        size=(simg2.shape[1], simg2.shape[0]),
        sigma=(0.5, 0.5),
        gridSpacing=(grid_spacing[0], grid_spacing[1]),
        gridOffset=(0, 0),
        spacing=(1, 1),
    )
    array = sitk.GetArrayViewFromImage(grid)
    
    # get mesh_size and number of control points
    mesh_size = np.array(transformDomainMeshSize)
    ctrl_pts = np.array(mesh_size) + spline_order
    
    # transform grid with outTx_Bspline parameters
    transform = sitk.BSplineTransformInitializer(grid, mesh_size.tolist(), spline_order)
    transform.SetParameters(outTx_Bspline.GetParameters())
    
    # Create array with the contour of the deformation
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(grid)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(100)
    resampler.SetOutputPixelType(sitk.sitkFloat32)
    resampled = resampler.Execute(grid)
    array_contour = sitk.GetArrayViewFromImage(resampled)
    
    # Visualize control points
    x_coeff, y_coeff = outTx_Bspline.GetCoefficientImages()
    grid_origin = x_coeff.GetOrigin()
    grid_spacing = x_coeff.GetSpacing()
    
    x = np.linspace(grid_origin[0], grid_origin[0] + (ctrl_pts[0] - 1) * grid_spacing[0], ctrl_pts[0])
    y = np.linspace(grid_origin[1], grid_origin[1] + (ctrl_pts[1] - 1) * grid_spacing[1], ctrl_pts[1])
    xx, yy = np.meshgrid(x, y)
    
    # Visualize transform on control points
    uv = np.reshape(np.array(transform.GetParameters()), (ctrl_pts[0], ctrl_pts[1], 2), order='F')
    u, v = uv[..., 0].T, uv[..., 1].T
    # Transpose y_values to make sure it has the same shape as x_values
    y_values = yy-v
    x_values = xx-u
    result_array = np.dstack((x_values, y_values))


    # Compute displacment field
    displacement_field = sitk.TransformToDisplacementField(transform,
                                      sitk.sitkVectorFloat64,
                                      grid.GetSize(),
                                      grid.GetOrigin(),
                                      grid.GetSpacing(),
                                      grid.GetDirection())
    displacement_field_np_arr = sitk.GetArrayViewFromImage(displacement_field)
    # Compute the jacobian determinant volume
    jacobian_det_volume = sitk.DisplacementFieldJacobianDeterminant(displacement_field)
    jacobian_det_np_arr = sitk.GetArrayViewFromImage(jacobian_det_volume)
    
    return array_contour, xx, yy, u, v, jacobian_det_np_arr, array


def correlation_rasters_report(outTx_Bspline_dict, outTx_Rigid, scale_percent1, scale_percent2, coords1, coords2, label1, label2, pixel_size_raster_fullres, image_shape1, image_shape2, operation1, operation2):
    correlations = {}
    num_images = len(outTx_Bspline_dict)
    num_cols = num_images
    num_rows = 1
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))
    if num_images > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    
    for ax, (ms, outTx_Bspline) in zip(axes, outTx_Bspline_dict.items()):
        # print(ms)

        coords1_tr = transform_coords(coords1, outTx_Bspline, outTx_Rigid, scale_percent1, scale_percent2, image_shape1, image_shape2, operation1, operation2)
    
        # Créer les rasters avec une taille fixe
        min_x = min(np.min(coords1_tr[:, 0]), np.min(coords2[:, 0])) - 1000
        max_x = max(np.max(coords1_tr[:, 0]), np.max(coords2[:, 0])) + 1000
        min_y = min(np.min(coords1_tr[:, 1]), np.min(coords2[:, 1])) - 1000
        max_y = max(np.max(coords1_tr[:, 1]), np.max(coords2[:, 1])) + 1000
        raster1 = create_raster_from_points(coords1_tr, pixel_size=pixel_size_raster_fullres, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
        raster2 = create_raster_from_points(coords2, pixel_size=pixel_size_raster_fullres, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
        
        correlation, p_value = pearsonr(raster1.flatten(), raster2.flatten())
        correlations[ms] = correlation

        
        new_array = np.zeros_like(raster1, dtype=np.uint8)
        new_array[(raster1 == 1) & (raster2 == 1)] = 3
        new_array[(raster1 == 1) & (raster2 == 0)] = 1
        new_array[(raster1 == 0) & (raster2 == 1)] = 2
        rgb_image = np.zeros((*new_array.shape, 3), dtype=np.uint8)
        colors = [[255, 255, 255], [0, 0, 255], [255, 0, 0], [0, 255, 0]]
        for i in range(4):
            rgb_image[new_array == i] = colors[i]
            
        colors_normalized = [[r / 255, g / 255, b / 255] for [r, g, b] in colors]
        legend1_patches = mpatches.Patch(color=colors_normalized[1], label=label1)
        legend2_patches = mpatches.Patch(color=colors_normalized[2], label=label2)
        
        ax.imshow(rgb_image)
        ax.axis('off')  # Turn off the axis
        ax.invert_yaxis()
        ax.legend(handles=[legend1_patches, legend2_patches], loc='upper left')
        ax.text(0.5, 0, f"Pearson correlation = {np.round(correlation, 3)}", transform=ax.transAxes, wrap=True, horizontalalignment='center', fontsize=12)
        #ax.figtext(0.5, 0, f"Pearson correlation = {np.round(correlation, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    
    
    for i in range(num_images, num_rows * num_cols):
        axes[i].axis('off')
    
    plt.tight_layout()
    # plt.savefig(dir_res + f"{id}/{folder_name_alignment}/simg2_ALLms.png", dpi=300)
    
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300)
    buffer.seek(0)
    plot_correlation = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    plt.close()

    return plot_correlation, correlations
    
def transform_coords(coords1, outTx_Bspline, outTx_Rigid, scale_percent1, scale_percent2, image_shape1, image_shape2, operation1, operation2):
    scaled_points1 = coords1*scale_percent1
    # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
    scaled_points1 = rotate_coordinates(scaled_points1, operation1, image_shape1)
    
    # Transform scaled points 1 in the coordinate system of slide 2
    scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points1])
    scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])

    # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
    scaled_pointsTR = rotate_coordinates(scaled_pointsTR, operation2, image_shape2)
    # Scale back to the full size image
    coords1_tr = scaled_pointsTR/scale_percent1

    return coords1_tr

# Create a raster from a set of points
def create_raster_from_points(points, pixel_size, min_x, max_x, min_y, max_y):
    width = int(np.ceil((max_x - min_x) / pixel_size))
    height = int(np.ceil((max_y - min_y) / pixel_size))
    
    transform = from_origin(min_x, max_y, pixel_size, pixel_size)
    
    raster_array = np.zeros((height, width), dtype=np.uint8)
    
    for point in points:
        x, y = point
        col = int((x - min_x) / pixel_size)
        row = int((max_y - y) / pixel_size)
        raster_array[row, col] = 1
    
    return raster_array

def rotate_coordinates(coords, operation, image_shape):
    """Applies the same transformation to coordinates as is done to the image."""
    height, width = image_shape[:2]
    
    if operation == 3:
        # Rotate counterclockwise by 90 degrees: (x, y) -> (-y, x)
        new_coords = np.column_stack((width - coords[:, 1], coords[:, 0]))
    elif operation == 2:
        # Flip left-right and up-down: (x, y) -> (-x, -y)
        new_coords = np.column_stack((width - coords[:, 0], height - coords[:, 1]))
    elif operation == 1:
        # Rotate clockwise by 90 degrees: (x, y) -> (y, -x)
        new_coords = np.column_stack((coords[:, 1], height - coords[:, 0]))
    else:
        # No operation, return coordinates as is
        new_coords = coords
    
    return new_coords



def merge_annotations(id, artefacts_empty_alignment, analysis_area_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, img_resize_alignment, metric_ms_alignment, metadata_images, panels_all, reference_panel, resolution_micron, output_path):
    for panel in panels_all:
        if panel is not reference_panel:
            ## Get annotations 
            artefacts_empty_gdf_poly = artefacts_empty_alignment[panel]
            analysis_area_gdf_poly = analysis_area_alignment[panel]
            ## Recover name alignment to get metadata
            name_alignments = metadata_images.keys()
            name_alignment = next(filter(lambda name: panel in name, name_alignments))
            scale_percent1 = metadata_images[name_alignment][f'scale_percent_{panel}']
            image_shape1 = metadata_images[name_alignment][f'image_shape_{panel}']
            operation1 = metadata_images[name_alignment][f'turn_img_{panel}']
            img_resize1 = img_resize_alignment[panel]
            ## Scale annotations to downscaled image resolution
            artefacts_empty_array = scale_gdf(artefacts_empty_gdf_poly, scale_percent1, operation1, image_shape1, img_resize1, c="B")
            analysis_area_array = scale_gdf(analysis_area_gdf_poly, scale_percent1, operation1, image_shape1, img_resize1, c="B")
            ## Transform scaled annotations to downscaled reference image
            scale_percent2 = metadata_images[name_alignment][f'scale_percent_{reference_panel}']
            image_shape2 = metadata_images[name_alignment][f'image_shape_{reference_panel}']
            operation2 = metadata_images[name_alignment][f'turn_img_{reference_panel}']
            img_resize2 = img_resize_alignment[reference_panel]
            metric_ms = metric_ms_alignment[panel]
            outTx_Rigid = outTx_Rigid_alignment[panel]
            outTx_Bspline = outTx_Bspline_alignment[panel]
            artefacts_empty_alignment[f'{panel}_in_{reference_panel}_panel'] = transform_annotation(artefacts_empty_array, scale_percent1, scale_percent2, image_shape1, image_shape2, operation1, operation2, metric_ms, outTx_Rigid, outTx_Bspline, img_resize2)
            analysis_area_alignment[f'{panel}_in_{reference_panel}_panel'] = transform_annotation(analysis_area_array, scale_percent1, scale_percent2, image_shape1, image_shape2, operation1, operation2, metric_ms, outTx_Rigid, outTx_Bspline, img_resize2)


    ## Fix invalid geometries in artefacts and empty annotations
    panels_tr = [key for key in analysis_area_alignment.keys() if "in" in key]
    panels_tr.append(reference_panel)
    analysis_area_alignment_list = []
    artefacts_empty_alignment_list = []
    for panel in panels_tr:
        analysis_area_alignment_list.append(analysis_area_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))
        artefacts_empty_alignment_list.append(artefacts_empty_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))

    ## Filter out empty GeoSeries
    analysis_area_alignment_list = [g for g in analysis_area_alignment_list if not g.empty]
    artefacts_empty_alignment_list = [g for g in artefacts_empty_alignment_list if not g.empty]

    ## Get the intersection of all analysis_area annotations
    # Step-by-step intersection
    intersections = intersect_polygon_lists(analysis_area_alignment_list[0], analysis_area_alignment_list[1])
    # Add the analysis_area annotation of the third panel if necessary
    if len(panels_all) == 3:
        intersections = intersect_polygon_lists(intersections, analysis_area_alignment_list[2])
    # Convert final intersections to GeoSeries
    intersections_polygons = gpd.GeoSeries(intersections)
    '''
    fig, ax = plt.subplots(figsize=(6, 6)) 
    intersections_polygons.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    plt.show()
    '''

    ## Merge the remaining GeoSeries
    merged_polygons = gpd.GeoSeries(unary_union([g.unary_union for g in artefacts_empty_alignment_list]))
    '''
    fig, ax = plt.subplots(figsize=(6, 6))
    merged_polygons.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    plt.show()
    '''

    # Create a Polygon representing the outer boundary
    outer_boundary = intersections_polygons.unary_union
    # Create a Polygon representing the inner boundary
    inner_boundary = merged_polygons.unary_union
    # Construct the mask by taking the difference between the outer boundary and inner boundary
    mask = outer_boundary.difference(inner_boundary)
    
    # Plot the MultiPolygon
    # plot_multipolygon(scaled_mask)
    
    ## Convert mask to micro meters
    scaled_mask = scale_multipolygon_coordinates(mask, resolution_micron)
    geojson_dict = mapping(scaled_mask)
    with open(output_path + f"Alignment/merged_annotations/{id}_mask.geojson", "w") as f:
        json.dump(geojson_dict, f, indent=4)

    return mask



def scale_gdf(gdf, scale_percent, operation, image_shape, img_resize, c = "R"):
    annotation = []
    for geometry in gdf.geometry:
        if geometry.geom_type == 'Polygon':
            annotation.append(np.array(geometry.exterior.coords))
            for hole in geometry.interiors:
                annotation.append(np.array(hole.coords))
        elif geometry.geom_type == 'MultiPolygon':
            for part in geometry:
                annotation.append(np.array(part.exterior.coords))
                for hole in part.interiors:
                    annotation.append(np.array(hole.coords))
    
    # Convert to polygon object
    polygons = []
    for array in annotation:
        scaled_points = array*scale_percent
        # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
        scaled_points = rotate_coordinates(scaled_points, operation, image_shape)
        # Convert array to list of tuples
        points = [(x, y) for x, y in scaled_points]
        
        # Create Polygon object
        polygons.append(Polygon(points))
    '''
    if c == "R":
        plt.imshow(img_resize, cmap="Reds")
    elif c == "B":
        plt.imshow(img_resize, cmap="Blues")
    # Plot the polygon
    for polygon in polygons:
        # Plot the exterior of the polygon as white
        plt.fill(*polygon.exterior.xy, color='gray', alpha=0.5)
    plt.show()
    '''

    return annotation



def transform_annotation(annotations, scale_percent1, scale_percent2, image_shape1, image_shape2, operation1, operation2, metric_ms, outTx_Rigid, outTx_Bspline, img_resize):
    # print(metric_ms)
    
    # Transform Shapes from IF2 to IF3 coordinate system
    annotation_TR= []
    for shape in annotations:
        scaled_points = shape*scale_percent1
        # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
        scaled_points = rotate_coordinates(scaled_points, operation1, image_shape1)
        
        # Transform scaled points 1 in the coordinate system of slide 2
        scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points])
        scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])

        # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
        scaled_pointsTR = rotate_coordinates(scaled_pointsTR, operation2, image_shape2)
        # Scale back to the full size image
        annotation_TR.append(scaled_pointsTR/scale_percent2)
    
    # Convert to polygon object
    polygonsTR = []
    for array in annotation_TR:
        scaled_points = array*scale_percent2
        # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
        scaled_points = rotate_coordinates(scaled_points, operation2, image_shape2)
        # Convert array to list of tuples
        points = [(x, y) for x, y in scaled_points]
        # Create Polygon object
        polygonsTR.append(Polygon(points))
    '''
    plt.imshow(img_resize, cmap="Reds")
    # Plot the polygon
    for polygon in polygonsTR:
        # Plot the exterior of the polygon as white
        plt.fill(*polygon.exterior.xy, color='gray', alpha=0.5)
    plt.show()
    '''

    polygons = [Polygon(coords) for coords in annotation_TR]
    geo_series = gpd.GeoSeries(polygons)

    return geo_series


# Function to find intersections between two lists of polygons
def intersect_polygon_lists(polygons1, polygons2):
    intersections = []
    for poly1 in polygons1:
        for poly2 in polygons2:
            intersection = poly1.intersection(poly2)
            if not intersection.is_empty:
                intersections.append(intersection)
    return intersections



def plot_multipolygon(multipolygon, linewidth=0.001):
    fig, ax = plt.subplots()

    if isinstance(multipolygon, MultiPolygon):
        polygons = multipolygon.geoms
    else:
        polygons = [multipolygon]

    for polygon in polygons:
        x, y = polygon.exterior.xy
        ax.plot(x, y, 'b', linewidth=linewidth)
        ax.fill(x, y, alpha=0.8, fc='lightblue', ec='blue')
        
        for interior in polygon.interiors:
            x, y = interior.xy
            ax.plot(x, y, 'r', linewidth=linewidth)
            ax.fill(x, y, alpha=0.9, fc='white', ec='red')

    ax.set_aspect('equal')
    ax.invert_yaxis()  # Invert the y-axis
    ax.axis('off')
    plt.show()


def scale_multipolygon_coordinates(multipolygon, resolution_um):
    # List to store the scaled polygons
    scaled_polygons = []
    
    # Iterate through each polygon in the multipolygon using multipolygon.geoms
    for polygon in multipolygon.geoms:
        # Get the exterior and interior coordinates of the polygon
        exterior_coords = np.array(polygon.exterior.coords)
        
        # Scale the exterior coordinates by dividing by resolution_um
        scaled_exterior = exterior_coords / resolution_um
        
        # For interior rings (holes) in the polygon
        scaled_interiors = []
        for interior in polygon.interiors:
            interior_coords = np.array(interior.coords)
            scaled_interior = interior_coords / resolution_um
            scaled_interiors.append(scaled_interior)
        
        # Create a new polygon with scaled coordinates
        scaled_polygon = Polygon(scaled_exterior, scaled_interiors)
        
        # Add the scaled polygon to the list
        scaled_polygons.append(scaled_polygon)
    
    # Return a new multipolygon with the scaled polygons
    return MultiPolygon(scaled_polygons)


def transform_filter_coordinates(metadata_images, cell_coordinates, data_frame_cells, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, mask, resolution_micron, merged_cell_coordinates):
    print("----------------------")
    print(panels_alignment[0])
    scale_percent1 = metadata_images[name_alignment][f'scale_percent_{panels_alignment[0]}']
    image_shape1 = metadata_images[name_alignment][f'image_shape_{panels_alignment[0]}']
    operation1 = metadata_images[name_alignment][f'turn_img_{panels_alignment[0]}']

    scale_percent2 = metadata_images[name_alignment][f'scale_percent_{panels_alignment[1]}']
    image_shape2 = metadata_images[name_alignment][f'image_shape_{panels_alignment[1]}']
    operation2 = metadata_images[name_alignment][f'turn_img_{panels_alignment[1]}']

    outTx_Rigid = outTx_Rigid_alignment[panels_alignment[0]]
    outTx_Bspline = outTx_Bspline_alignment[panels_alignment[0]]

    coords = cell_coordinates[f'{panels_alignment[0]}_panel_DAPI']
    ids = np.arange(len(coords))


    scaled_points = coords*scale_percent1
    # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
    scaled_points = rotate_coordinates(scaled_points, operation1, image_shape1)

    # Transform scaled points 1 in the coordinate system of slide 2
    scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points])
    scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])
    
    # Transform scaled coordinates regarding the induced rotation that as been done to the corresponding image to facilitate the alignment
    scaled_pointsTR = rotate_coordinates(scaled_pointsTR, operation2, image_shape2)
    # Scale back to the full size image
    coords_tr = scaled_pointsTR/scale_percent2
    cell_coordinates[f'{panels_alignment[0]}_panel_DAPI_transformed'] = coords_tr

    # Create boolean mask indicating which points are inside any polygon in the mask
    mask_contains = shapely.vectorized.contains(mask, coords_tr[:, 0], coords_tr[:, 1])
    # Filter points based on the boolean mask
    filtered_coords_tr = coords_tr[mask_contains]
    filtered_ids = ids[mask_contains]
    
    print('Nb cells before :', len(coords_tr))
    print('Nb cells after filrtering :', len(filtered_coords_tr))

    # Plot coordinates before filtering
    plt.figure(figsize=(5, 5))
    plt.scatter(coords_tr[:, 0], coords_tr[:, 1], c='blue', label='Before', alpha=0.5, s=1)

    # Plot coordinates after filtering
    filtered_coords = np.array(filtered_coords_tr)
    plt.scatter(filtered_coords[:, 0], filtered_coords[:, 1], c='red', label='After Filter', s=0.01)

    # Add legend and labels
    plt.legend()
    #plt.grid(True)
    plt.gca().invert_yaxis()

    # Original DataFrame
    df_cells = data_frame_cells[f'{panels_alignment[0]}_panel_df']
    # Filter the DataFrame by `filtered_ids`
    filtered_df = df_cells.loc[filtered_ids, ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']].copy()
    filtered_df[['x (µm)', 'y (µm)']] = filtered_coords_tr / resolution_micron
    filtered_df = filtered_df[['x (µm)', 'y (µm)', 'Cell_type', 'Phenotype', 'Classifier Label', 'Object Id']]
    filtered_df

    filtered_df['Panel'] = panels_alignment[0]
    merged_cell_coordinates = pd.concat([merged_cell_coordinates, filtered_df], ignore_index=True)
    
    return merged_cell_coordinates


def filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, merged_cell_coordinates):
    print("----------------------")
    print(panels_alignment[1])
    coords = cell_coordinates[f'{panels_alignment[1]}_panel_DAPI']
    # Generate IDs
    ids = np.arange(len(coords))

    # Create boolean mask indicating which points are inside any polygon in the mask
    mask_contains = shapely.vectorized.contains(mask, coords[:, 0], coords[:, 1])
    # Filter points based on the boolean mask
    filtered_coords = coords[mask_contains]
    filtered_ids = ids[mask_contains]

    print('Nb cells before :', len(coords))
    print('Nb cells after filrtering :', len(filtered_coords))

    # Plot coordinates before filtering
    plt.figure(figsize=(5, 5))
    plt.scatter(coords[:, 0], coords[:, 1], c='blue', label='Before', alpha=0.5, s=1)

    # Plot coordinates after filtering
    plt.scatter(filtered_coords[:, 0], filtered_coords[:, 1], c='red', label='After Filter', s=0.01)

    # Add legend and labels
    plt.legend()
    #plt.grid(True)
    plt.gca().invert_yaxis()

    # Original DataFrame
    df_cells = data_frame_cells[f'{panels_alignment[1]}_panel_df']
    # Filter the DataFrame by `filtered_ids`
    filtered_df = df_cells.loc[filtered_ids, ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']].copy()
    filtered_df[['x (µm)', 'y (µm)']] = filtered_coords / 2.01294825
    filtered_df = filtered_df[['x (µm)', 'y (µm)', 'Cell_type', 'Phenotype', 'Classifier Label', 'Object Id']]
    filtered_df

    filtered_df['Panel'] = panels_alignment[1]
    merged_cell_coordinates = pd.concat([merged_cell_coordinates, filtered_df], ignore_index=True)

    return merged_cell_coordinates


def save_tables(merged_cell_coordinates, output_path, id):
    merged_cell_coordinates["label"] = merged_cell_coordinates.index
    merged_cell_coordinates.rename(columns={"Classifier Label": "Classifier_Label"}, inplace=True)
    merged_cell_coordinates.drop(columns=["Object Id"], inplace=True)
    print("Combined phenotypes:")
    print(np.unique(merged_cell_coordinates['Cell_type']))

    merged_cell_coordinates.to_csv(output_path + f"Alignment/merged_tables/{id}_merged_cell_coordinates.csv", index=False)



def plot_rasters(data_frame_cells, merged_cell_coordinates, cell_coordinates, pixel_size_raster_micron, output_path, panels_all, id, resolution_micron):
    ## Plot the raster before alignment
    coords1 = data_frame_cells[f'{panels_all[0]}_panel_df'][["x (µm)", "y (µm)"]].to_numpy()
    coords2 = data_frame_cells[f'{panels_all[1]}_panel_df'][["x (µm)", "y (µm)"]].to_numpy()
    coords3 = data_frame_cells[f'{panels_all[2]}_panel_df'][["x (µm)", "y (µm)"]].to_numpy()
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords1[:, 0]), np.min(coords2[:, 0]), np.min(coords3[:, 0])) - 1000
    max_x = max(np.max(coords1[:, 0]), np.max(coords2[:, 0]), np.max(coords3[:, 0])) + 1000
    min_y = min(np.min(coords1[:, 1]), np.min(coords2[:, 1]), np.min(coords3[:, 1])) - 1000
    max_y = max(np.max(coords1[:, 1]), np.max(coords2[:, 1]), np.max(coords3[:, 1])) + 1000
    # Create rasters with a fixed size
    raster1 = create_raster_from_points(coords1, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster2 = create_raster_from_points(coords2, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster3 = create_raster_from_points(coords3, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    # Compute the pearson correlation between rasters
    correlation_1_2, p_value = pearsonr(raster1.flatten(), raster2.flatten())
    correlation_1_3, p_value = pearsonr(raster1.flatten(), raster3.flatten())
    correlation_2_3, p_value = pearsonr(raster2.flatten(), raster3.flatten())
    # Plot the raster
    plt.figure(figsize=(4, 4))
    colors1 = ['white', 'blue']
    custom_cmap1 = plt.cm.colors.ListedColormap(colors1)
    colors2 = ['white', 'green']
    custom_cmap2 = plt.cm.colors.ListedColormap(colors2)
    colors3 = ['white', 'red']
    custom_cmap3 = plt.cm.colors.ListedColormap(colors3)
    legend1_patches = mpatches.Patch(color=colors1[1], label='DAPI T')
    legend2_patches = mpatches.Patch(color=colors2[1], label='DAPI DC')
    legend3_patches = mpatches.Patch(color=colors3[1], label='DAPI TLS')
    plt.imshow(raster1, cmap=custom_cmap1, extent=[min_x, max_x, min_y, max_y], aspect='auto')
    plt.imshow(raster2, cmap=custom_cmap2, extent=[min_x, max_x, min_y, max_y], alpha = 0.5, aspect='auto')
    plt.imshow(raster3, cmap=custom_cmap3, extent=[min_x, max_x, min_y, max_y], alpha = 0.4, aspect='auto')
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.legend(handles=[legend1_patches, legend2_patches, legend3_patches], loc='upper right')
    plt.figtext(0.5, 0, f"Correlation 1-2 : {np.round(correlation_1_2, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.05, f"Correlation 1-3 : {np.round(correlation_1_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.1, f"Correlation 2-3 : {np.round(correlation_2_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(output_path + f"Alignment/plots/{id}_raster_before_alignment.png", dpi=300, bbox_inches='tight')
    plt.close()

    ## Plot the common raster after alignment
    coords1 = cell_coordinates[f'{panels_all[0]}_panel_DAPI_transformed'] / resolution_micron
    coords2 = cell_coordinates[f'{panels_all[1]}_panel_DAPI'] / resolution_micron
    coords3 = cell_coordinates[f'{panels_all[2]}_panel_DAPI_transformed'] / resolution_micron
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords1[:, 0]), np.min(coords2[:, 0]), np.min(coords3[:, 0])) - 1000
    max_x = max(np.max(coords1[:, 0]), np.max(coords2[:, 0]), np.max(coords3[:, 0])) + 1000
    min_y = min(np.min(coords1[:, 1]), np.min(coords2[:, 1]), np.min(coords3[:, 1])) - 1000
    max_y = max(np.max(coords1[:, 1]), np.max(coords2[:, 1]), np.max(coords3[:, 1])) + 1000
    # Create rasters with a fixed size  
    raster1 = create_raster_from_points(coords1, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster2 = create_raster_from_points(coords2, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster3 = create_raster_from_points(coords3, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    # Compute the pearson correlation between rasters
    correlation_1_2, p_value = pearsonr(raster1.flatten(), raster2.flatten())
    correlation_1_3, p_value = pearsonr(raster1.flatten(), raster3.flatten())
    correlation_2_3, p_value = pearsonr(raster2.flatten(), raster3.flatten())
    # Create a new array with the same shape initialized to zeros
    new_array = np.zeros_like(raster1, dtype=np.uint8)
    # Apply the rules to create the new array
    new_array[(raster1 == 1) & (raster2 == 1) & (raster3 == 1)] = 3
    new_array[(raster1 == 1) & (raster2 == 1) & (raster3 == 0)] = 2
    new_array[(raster1 == 0) & (raster2 == 1) & (raster3 == 1)] = 2
    new_array[(raster1 == 1) & (raster2 == 0) & (raster3 == 0)] = 1
    new_array[(raster1 == 0) & (raster2 == 1) & (raster3 == 0)] = 1
    new_array[(raster1 == 0) & (raster2 == 0) & (raster3 == 1)] = 1
    # Create an RGB image with the same height and width as new_array
    rgb_image = np.zeros((*new_array.shape, 3), dtype=np.uint8)
    # Assign specific colors based on the pixel values
    colors = [[255, 255, 255], [255, 0, 0], [255, 165, 0], [0, 255, 0]]
    for i in range(4):
        rgb_image[new_array == i] = colors[i]
    colors_normalized = [[r / 255, g / 255, b / 255] for [r, g, b] in colors]
    legend1_patches = mpatches.Patch(color=colors_normalized[1], label='Only one')
    legend2_patches = mpatches.Patch(color=colors_normalized[2], label='Part common')
    legend3_patches = mpatches.Patch(color=colors_normalized[3], label='Common')
    # Plot the RGB image
    plt.imshow(rgb_image)
    plt.axis('off')  # Turn off the axis
    plt.gca().invert_yaxis()
    plt.legend(handles=[legend1_patches, legend2_patches, legend3_patches], loc='upper right')
    plt.figtext(0.5, 0, f"Correlation 1-2 : {np.round(correlation_1_2, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.05, f"Correlation 1-3 : {np.round(correlation_1_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.1, f"Correlation 2-3 : {np.round(correlation_2_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(output_path + f"Alignment/plots/{id}_common_raster.png", dpi=300, bbox_inches='tight')
    plt.close()


    ## Plot the common raster after alignment and filtering with common annotations
    coords1 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[0]}'][["x (µm)", "y (µm)"]].to_numpy()
    coords2 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[1]}'][["x (µm)", "y (µm)"]].to_numpy()
    coords3 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[2]}'][["x (µm)", "y (µm)"]].to_numpy()
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords1[:, 0]), np.min(coords2[:, 0]), np.min(coords3[:, 0])) - 1000
    max_x = max(np.max(coords1[:, 0]), np.max(coords2[:, 0]), np.max(coords3[:, 0])) + 1000
    min_y = min(np.min(coords1[:, 1]), np.min(coords2[:, 1]), np.min(coords3[:, 1])) - 1000
    max_y = max(np.max(coords1[:, 1]), np.max(coords2[:, 1]), np.max(coords3[:, 1])) + 1000
    # Create rasters with a fixed size  
    raster1 = create_raster_from_points(coords1, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster2 = create_raster_from_points(coords2, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster3 = create_raster_from_points(coords3, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    # Compute the pearson correlation between rasters
    correlation_1_2, p_value = pearsonr(raster1.flatten(), raster2.flatten())
    correlation_1_3, p_value = pearsonr(raster1.flatten(), raster3.flatten())
    correlation_2_3, p_value = pearsonr(raster2.flatten(), raster3.flatten())
    # Create a new array with the same shape initialized to zeros
    new_array = np.zeros_like(raster1, dtype=np.uint8)
    # Apply the rules to create the new array
    new_array[(raster1 == 1) & (raster2 == 1) & (raster3 == 1)] = 3
    new_array[(raster1 == 1) & (raster2 == 1) & (raster3 == 0)] = 2
    new_array[(raster1 == 0) & (raster2 == 1) & (raster3 == 1)] = 2
    new_array[(raster1 == 1) & (raster2 == 0) & (raster3 == 0)] = 1
    new_array[(raster1 == 0) & (raster2 == 1) & (raster3 == 0)] = 1
    new_array[(raster1 == 0) & (raster2 == 0) & (raster3 == 1)] = 1
    # Create an RGB image with the same height and width as new_array
    rgb_image = np.zeros((*new_array.shape, 3), dtype=np.uint8)
    # Assign specific colors based on the pixel values
    colors = [[255, 255, 255], [255, 0, 0], [255, 165, 0], [0, 255, 0]]
    for i in range(4):
        rgb_image[new_array == i] = colors[i] 
    colors_normalized = [[r / 255, g / 255, b / 255] for [r, g, b] in colors]
    legend1_patches = mpatches.Patch(color=colors_normalized[1], label='Only one')
    legend2_patches = mpatches.Patch(color=colors_normalized[2], label='Part common')
    legend3_patches = mpatches.Patch(color=colors_normalized[3], label='Common')
    # Plot the raster
    plt.imshow(rgb_image)
    plt.axis('off')  # Turn off the axis
    plt.gca().invert_yaxis()
    plt.legend(handles=[legend1_patches, legend2_patches, legend3_patches], loc='upper right')
    plt.figtext(0.5, 0, f"Correlation 1-2 : {np.round(correlation_1_2, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.05, f"Correlation 1-3 : {np.round(correlation_1_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.1, f"Correlation 2-3 : {np.round(correlation_2_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(output_path + f"Alignment/plots/{id}_common_raster_filtered_annotations.png", dpi=300, bbox_inches='tight')
    plt.close()
    

    ## Plot the raster after alignment
    coords1_tr = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[0]}'][["x (µm)", "y (µm)"]].to_numpy()
    coords2 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[1]}'][["x (µm)", "y (µm)"]].to_numpy()
    coords3_tr = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[2]}'][["x (µm)", "y (µm)"]].to_numpy()
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords1_tr[:, 0]), np.min(coords2[:, 0]), np.min(coords3_tr[:, 0])) - 1000
    max_x = max(np.max(coords1_tr[:, 0]), np.max(coords2[:, 0]), np.max(coords3_tr[:, 0])) + 1000
    min_y = min(np.min(coords1_tr[:, 1]), np.min(coords2[:, 1]), np.min(coords3_tr[:, 1])) - 1000
    max_y = max(np.max(coords1_tr[:, 1]), np.max(coords2[:, 1]), np.max(coords3_tr[:, 1])) + 1000
    # Create rasters with a fixed size
    raster1 = create_raster_from_points(coords1_tr, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster2 = create_raster_from_points(coords2, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    raster3 = create_raster_from_points(coords3_tr, pixel_size=pixel_size_raster_micron, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
    # Compute the pearson correlation between rasters
    correlation_1_2, p_value = pearsonr(raster1.flatten(), raster2.flatten())
    correlation_1_3, p_value = pearsonr(raster1.flatten(), raster3.flatten())
    correlation_2_3, p_value = pearsonr(raster2.flatten(), raster3.flatten())
    # Plot the rasters
    plt.figure(figsize=(4, 4))
    legend1_patches = mpatches.Patch(color=colors1[1], label='DAPI T')
    legend2_patches = mpatches.Patch(color=colors2[1], label='DAPI DC')
    legend3_patches = mpatches.Patch(color=colors3[1], label='DAPI TLS')
    plt.imshow(raster1, cmap=custom_cmap1, extent=[min_x, max_x, min_y, max_y], aspect='auto')
    plt.imshow(raster2, cmap=custom_cmap2, extent=[min_x, max_x, min_y, max_y], alpha = 0.5, aspect='auto')
    plt.imshow(raster3, cmap=custom_cmap3, extent=[min_x, max_x, min_y, max_y], alpha = 0.3, aspect='auto')
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.legend(handles=[legend1_patches, legend2_patches, legend3_patches], loc='upper right')
    plt.figtext(0.5, 0, f"Correlation 1-2 : {np.round(correlation_1_2, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.05, f"Correlation 1-3 : {np.round(correlation_1_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.5, -0.1, f"Correlation 2-3 : {np.round(correlation_2_3, 3)}", wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(output_path + f"Alignment/plots/{id}_raster_alignment.png", dpi=300, bbox_inches='tight')
    plt.close()



