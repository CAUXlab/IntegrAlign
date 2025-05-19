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
from shapely.geometry import box
from shapely.plotting import plot_polygon
import xml.etree.ElementTree as ET
import geopandas as gpd
from shapely.ops import unary_union
import json
import shapely.vectorized
import matplotlib.colors as mcolors


def alignment(downscaled_images_path, coordinate_tables, resolution_micron, number_ms, metric, pixel_size_raster_micron, alpha_red):
    ## Load and read the pickle file
    with open(downscaled_images_path, "rb") as file:
        downscaled_images = pickle.load(file)

    ## Get the parameters
    annotations_paths = downscaled_images["params"].get("annotations_paths")
    annotations_names_AnalysisArea = downscaled_images["params"].get("annotations_names_AnalysisArea")
    annotations_names_empty = downscaled_images["params"].get("annotations_names_empty")
    annotations_names_artefacts = downscaled_images["params"].get("annotations_names_artefacts")
    common_ids = downscaled_images["params"].get("common_ids")
    panels_all = downscaled_images["params"].get("panels")
    # print(panels_all)
    output_path = downscaled_images["params"].get("output_path")

    ## Remove "params" from downscaled_images
    downscaled_images = remove_params(downscaled_images)
    ## Create alignment folder to store alignment reports
    Path(output_path+'Alignment/reports').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/merged_annotations').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/merged_tables').mkdir(parents=True, exist_ok=True)
    Path(output_path+'Alignment/plots').mkdir(parents=True, exist_ok=True)
    
    
    ## Get reference panel
    reference_panel = panels_all[1]
    breaking_id = True
    for id, downscaled_images_id in downscaled_images.items():
        '''
        print(id)
        if id == "11005":
            breaking_id = False
        if breaking_id:
            continue
        '''
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
        for name_alignment, downscaled_images_id_name_alignment in downscaled_images_id.items():
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
             simg2_Rigid, nda_Rigid) = alignment_2panels(name_alignment, id, img1, img2, number_ms, metric, spline_order)
            ## Load coordinate table and annotations of each panel
            print(f'Loading tables and annotations...')
            for panel in panels_alignment:
                # Get the panel path
                index = panels_all.index(panel)
                coordinate_table = coordinate_tables[index]
                # Get the file path corresponding to id
                csv_file_path = next((os.path.join(coordinate_table, f) for f in os.listdir(coordinate_table) if id in f and f.endswith(".csv") and not f.startswith(".")), None)
                ## Get the coordinates table
                cell_coordinates, data_frame_cells = get_cells_coordinates_SPIAT_CellType(csv_file_path, panel, cell_coordinates, data_frame_cells, resolution_micron)
                ## Get the annotations if given
                if annotations_paths:
                    annotations = annotations_paths[index]
                    # geojson_file_path = next((os.path.join(annotations, f) for f in os.listdir(annotations) if id in f and f.endswith(".geojson") and not f.startswith(".")), None)
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
        if annotations_paths:
            ## Transform and merge the annotations
            print("-----")
            print("Merging the annotations...")    
            mask = merge_annotations(id, artefacts_empty_alignment, analysis_area_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, 
                            img_resize_alignment, metric_ms_alignment, metadata_images, panels_all, reference_panel, resolution_micron, output_path)
        else:
            mask = "No_annotations"
            
        # If the intersection of the analysis area is empty
        # Doesn't save merge coordinates, tables or plots
        # Continue to the next iteration
        if mask is None:
            continue  # Skip to the next iteration
        
        ## Transform and filter coordinates
        print("Transform, filter and merge coordinates...")
        merged_cell_coordinates = pd.DataFrame()
        for name_alignment in downscaled_images_id.keys():
            panels_alignment = name_alignment.split("_")
            merged_cell_coordinates = transform_filter_coordinates(metadata_images, cell_coordinates, data_frame_cells, name_alignment, panels_alignment, outTx_Rigid_alignment, outTx_Bspline_alignment, mask, resolution_micron, merged_cell_coordinates)
        merged_cell_coordinates = filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, resolution_micron, merged_cell_coordinates)
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

def extract_downscaled_images(downscaled_images_id_name_alignment, panels, name_alignment, metadata_images, id):
    ## Load images
    img1 = downscaled_images_id_name_alignment.get("img1")
    img2 = downscaled_images_id_name_alignment.get("img2")
    img1_resize = downscaled_images_id_name_alignment.get("img1_resize")
    img2_resize = downscaled_images_id_name_alignment.get("img2_resize")
    panel1, panel2 = panels
    ## Get metadata
    metadata_images[name_alignment] = {
        f"channels_{panel1}": downscaled_images_id_name_alignment.get("channels1"),
        f"scale_percent_{panel1}": downscaled_images_id_name_alignment.get("scale_percent1"),
        f"crop_coords_{panel1}": downscaled_images_id_name_alignment.get("crop_coords1"),
        f"img_resize_ori_{panel1}": downscaled_images_id_name_alignment.get("img1_resize_ori"),
        f"image_shape_{panel1}": downscaled_images_id_name_alignment.get("img1_shape_ori"),
        f"manual_alignment_displacement": downscaled_images_id_name_alignment.get("manual_alignment_displacement"),
        f"manual_alignment_rotation": downscaled_images_id_name_alignment.get("manual_alignment_rotation_angle"),
        f"manual_alignment_rotation_shape": downscaled_images_id_name_alignment.get("manual_alignment_rotation_shape"),
        f"image_shape_manual_alignment": img1_resize.shape,

        f"channels_{panel2}": downscaled_images_id_name_alignment.get("channels2"),
        f"scale_percent_{panel2}": downscaled_images_id_name_alignment.get("scale_percent2"),
        f"crop_coords_{panel2}": downscaled_images_id_name_alignment.get("crop_coords2"),
        f"img_resize_ori_{panel2}": downscaled_images_id_name_alignment.get("img2_resize_ori"),
        f"image_shape_{panel2}": downscaled_images_id_name_alignment.get("img2_shape_ori"),
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


def get_cells_coordinates_SPIAT_CellType(csv_file_path, panel, cell_coordinates, data_frame_cells, resolution_micron):
    df = pd.read_csv(csv_file_path)
    
    new_df = pd.DataFrame()

    required_columns = ['Phenotype', 'Classifier.Label', 'Object.Id', 'x', 'y']

    required_columns_SPIAT = [
    'x', 'y', 'Phenotype', 'Cell_type', 'ID', 'Tissue_Category',
    'XMin_pixel', 'XMax_pixel', 'YMin_pixel', 'YMax_pixel' 
    ]

    required_columns_ROI = ['x_', 'y_', 'x', 'y']

    required_columns_HALO = ['XMin', 'XMax', 'YMin', 'YMax']

    # Check if it is a table from SPIAT with cell types or not
    if all(col in df.columns for col in required_columns):
        # Store the coordinates in the new DataFrame
        ## Convert micrometers coordinates to pixel
        
        '''
        new_df['x'] = (df['XMin'] + df['XMax'])/2
        new_df['y'] = (df['YMin'] + df['YMax'])/2
        '''

        new_df['x'] = df['x']
        new_df['y'] = df['y']

        new_df['x (micron)'] = new_df['x'] / resolution_micron
        new_df['y (micron)'] = new_df['y'] / resolution_micron
            
        # new_df['x'] = ((df['XMin_pixel'] + df['XMax_pixel']) / 2)
        # new_df['y'] = ((df['YMin_pixel'] + df['YMax_pixel']) / 2)
        
        new_df['Phenotype'] = df['Phenotype']
        new_df['Object.Id'] = df['Object.Id']
        new_df['Classifier.Label'] = df['Classifier.Label']

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

    # Check if it is a table from SPIAT with cell types or not
    elif all(col in df.columns for col in required_columns_SPIAT):
        # Store the coordinates in the new DataFrame
        ## Convert micrometers coordinates to pixel
        new_df['x (micron)'] = df['x']
        new_df['y (micron)'] = df['y']

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

    # Check if all required ROI columns are present in the DataFrame
    elif all(col in df.columns for col in required_columns_ROI):
        new_df = df.copy()

        new_df['x (micron)'] = df['x_']
        new_df['y (micron)'] = df['y_']


        cell_coordinates[f'{panel}_panel_DAPI'] = new_df[['x', 'y']].to_numpy()
        data_frame_cells[f'{panel}_panel_df'] = new_df

    # Check if all required HALO columns are present in the DataFrame
    elif all(col in df.columns for col in required_columns_HALO):
        new_df = df.copy()

        new_df['x'] = (df['XMin'] + df['XMax'])/2
        new_df['y'] = (df['YMin'] + df['YMax'])/2

        new_df['x (micron)'] = new_df['x'] / resolution_micron
        new_df['y (micron)'] = new_df['y'] / resolution_micron

        cell_coordinates[f'{panel}_panel_DAPI'] = new_df[['x', 'y']].to_numpy()
        data_frame_cells[f'{panel}_panel_df'] = new_df

    else:
        raise ValueError("The tables are not in the correct format.\n Required columns : ['Phenotype', 'Classifier.Label', 'Object.Id', 'x', 'y'] or simply ['XMin', 'XMax', 'YMin', 'YMax']")


    return cell_coordinates, data_frame_cells


def get_annotations(annotation_file_path, panel, artefacts_empty_alignment, analysis_area_alignment, annotations_names_Empty, annotations_names_Artefacts, annotations_names_AnalysisArea):
    if annotation_file_path.endswith(".annotations"):
        # Parse the XML file
        tree = ET.parse(annotation_file_path)
        root = tree.getroot()

        if annotations_names_AnalysisArea:
            # Filter analysis area
            gdf = get_gdf_from_annot(root)
            analysis_area_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.lower() for item in annotations_names_AnalysisArea))]
            if analysis_area_gdf["geometry"].is_empty.all():
                print(f"No annotations named {annotations_names_AnalysisArea}")
            else:
                analysis_area_alignment[panel] = analysis_area_gdf["geometry"]
                if annotations_names_Empty or annotations_names_Artefacts:
                    annotations_names_Empty_Artefacts = annotations_names_Empty + annotations_names_Artefacts
                    gdf = get_gdf_from_annot(root)
                    artefacts_empty_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.lower() for item in annotations_names_Empty_Artefacts))]
                    if artefacts_empty_gdf["geometry"].is_empty.all():
                        print(f"No annotations named {annotations_names_Empty_Artefacts}")
                    else:
                        artefacts_empty_alignment[panel] = artefacts_empty_gdf["geometry"]
        else:
            if annotations_names_Empty:
                ## Get the Empty areas
                gdf = get_gdf_from_annot(root, NegativeROA = 0)
                empty_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.lower() for item in annotations_names_Empty))]
                # Find the index of the largest shape
                largest_geometry_index = empty_gdf["geometry"].area.idxmax()
                '''
                # Get the largest shape
                largest_geometry_gdf = artefacts_empty_gdf.loc[[largest_geometry_index]]
                analysis_area_alignment[panel] = largest_geometry_gdf["geometry"]
                '''
                # Remove the largest geometry (full slide annotation) from the empty shapes
                artefacts_empty_alignment[panel] = empty_gdf.drop(largest_geometry_index)["geometry"]

                ## Get the Analysis area (here = Negative to empty areas)
                gdf = get_gdf_from_annot(root, NegativeROA = 1)
                analysis_area_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.lower() for item in annotations_names_Empty))]
                analysis_area_alignment[panel] = analysis_area_gdf["geometry"]
            else:
                print(f"No annotations named {annotations_names_Empty}")
            if annotations_names_Artefacts:
                ## Get the Artefacts areas
                gdf = get_gdf_from_annot(root)
                artefacts_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.lower() for item in annotations_names_Artefacts))]
                artefacts_empty_alignment[panel] = pd.concat([
                    artefacts_empty_alignment[panel],
                    artefacts_gdf["geometry"]
                ])
            else:
                print(f"No annotations named {annotations_names_Empty}")











    if annotation_file_path.endswith(".geojson"):
        ## Read annotations file
        gdf = gpd.read_file(annotation_file_path)
        def safe_json_load(x):
            if isinstance(x, dict):  # If it's already a dictionary, return it as is
                return x
            elif isinstance(x, str):  # If it's a string, try to parse it as JSON
                try:
                    return json.loads(x)
                except json.JSONDecodeError:
                    return {}  # Return an empty dictionary for invalid rows
            else:  # If it's neither, return an empty dictionary
                return {}
        gdf['classification'] = gdf['classification'].apply(safe_json_load)
        '''
        # print(gdf['classification'].apply(lambda x: x.get('name')).unique())
        ## Convert LineString to Polygon
        gdf_poly = gdf['geometry'].apply(lambda geom: Polygon(list(geom.coords) + [geom.coords[0]]))
        fig, ax = plt.subplots(figsize=(6, 6)) 
        gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
        ax.axis('off')
        ax.invert_yaxis()
        '''

        ## Check for the annotation types (non-case sensitive)
        # annotations_empty = ['Artefacts', 'Manual_Artefacts', 'Empty', 'No_tissue']
        # artefacts_empty_gdf = gdf[gdf['classification'].apply(lambda x: x.get('name').lower() in [item.lower() for item in annotations_empty])]
        artefacts_empty_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.get('name').lower() for item in annotations_names_Empty))]

        # print(artefacts_empty_gdf['classification'].apply(lambda x: x.get('name')).unique())
        ## Convert LineString to Polygon
        artefacts_empty_gdf_poly = artefacts_empty_gdf['geometry'].apply(lambda geom: Polygon(list(geom.coords) + [geom.coords[0]]))
        '''
        fig, ax = plt.subplots(figsize=(6, 6)) 
        artefacts_empty_gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
        ax.axis('off')
        ax.invert_yaxis()
        plt.show()
        '''
        artefacts_empty_alignment[panel] = artefacts_empty_gdf_poly

        ## Get the analysis area annotation
        # analysis_area_gdf = gdf[gdf['classification'].apply(lambda x: x.get('name') in ['Analysis_area'])]
        analysis_area_gdf = gdf[gdf['classification'].apply(lambda x: any(item.lower() in x.get('name').lower() for item in annotations_names_AnalysisArea))]
        # Convert LineString to Polygon
        analysis_area_gdf_poly = analysis_area_gdf['geometry'].apply(lambda geom: Polygon(list(geom.coords) + [geom.coords[0]]))
        '''
        fig, ax = plt.subplots(figsize=(6, 6)) 
        analysis_area_gdf_poly.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
        ax.axis('off')
        ax.invert_yaxis()
        plt.show()
        '''
        analysis_area_alignment[panel] = analysis_area_gdf_poly

    return artefacts_empty_alignment, analysis_area_alignment

def get_gdf_from_annot(root, NegativeROA = "all"):
    # List to store annotation polygons
    geometries = []
    classification_labels = []
    
    for annotation in root.findall(".//Annotation"):
        annotation_name = annotation.get("Name")
        # Extract coordinates of the region
        if NegativeROA == 0:
            regions = annotation.findall(".//Region[@NegativeROA='0']")
        elif NegativeROA == 1:
            regions = annotation.findall(".//Region[@NegativeROA='1']")
        else:
            regions = annotation.findall(".//Region")
        for region in regions:
            vertices = region.findall(".//V")
            points = [(int(vertex.get("X")), int(vertex.get("Y"))) for vertex in vertices]

            # Ensure at least 3 points to form a polygon
            if len(points) >= 3:
                polygon = Polygon(points)
                geometries.append(polygon)
                classification_labels.append(annotation_name)
    
    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame({"classification": classification_labels, "geometry": geometries})
    '''
        # Check if annotation matches classification categories (case-insensitive)
        if annotation_name and annotation_name.lower() in annotations_classification:
            print(f"Matched classification: {annotation_name}")
    
            # Extract coordinates of the region
            regions = annotation.findall(".//Region")
            
            for region in regions:
                vertices = region.findall(".//V")
                points = [(int(vertex.get("X")), int(vertex.get("Y"))) for vertex in vertices]
    
                # Ensure at least 3 points to form a polygon
                if len(points) >= 3:
                    polygon = Polygon(points)
                    geometries.append(polygon)
                    classification_labels.append(annotation_name)
    
    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame({"classification": classification_labels, "geometry": geometries})
    '''

    return gdf


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
    # Ensure axs is always a 2D array
    if len(outTx_Bspline_dict) == 1:
        axs = axs.reshape(2, 1)
    

    for i, (metric_ms, outTx_Bspline) in enumerate(outTx_Bspline_dict.items()):
        simg2 = sitk.GetArrayFromImage(simg2_dict[metric_ms])
        transformDomainMeshSize = [int(metric_ms.split('_')[1]), int(metric_ms.split('_')[1])]

        axs = grid_deform_detJacobian(simg2, transformDomainMeshSize, spline_order, outTx_Bspline, i, axs, thresh_pixelSIMG2, simg2_tr)

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
    crop_coords1 = metadata_images[name_alignment][f'crop_coords_{panels[0]}']
    crop_coords2 = metadata_images[name_alignment][f'crop_coords_{panels[1]}']
    manual_alignment_displacement = metadata_images[name_alignment]['manual_alignment_displacement']
    manual_alignment_rotation = metadata_images[name_alignment]['manual_alignment_rotation']
    manual_alignment_rotation_shape = metadata_images[name_alignment]['manual_alignment_rotation_shape']
    image_shape_manual_alignment = metadata_images[name_alignment]['image_shape_manual_alignment']
    img1_resize_ori = metadata_images[name_alignment][f'img_resize_ori_{panels[0]}']
    img2_resize_ori = metadata_images[name_alignment][f'img_resize_ori_{panels[1]}']
    coords1 = cell_coordinates[f'{panels[0]}_panel_DAPI']
    coords2 = cell_coordinates[f'{panels[1]}_panel_DAPI']
    label1 = f'DAPI from panel {panels[0]}'
    label2 = f'DAPI {panels[1]}'
    plot_correlation, correlations = correlation_rasters_report(outTx_Bspline_dict, outTx_Rigid, scale_percent1, scale_percent2, coords1, coords2, label1, label2, pixel_size_raster_fullres, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, img1_resize, img2_resize, img1_resize_ori, img2_resize_ori)


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

    if len(outTx_Bspline_dict) == 1:
        file_name = output_path + f"Alignment/reports/{id}_{name_alignment}_report_alignment_refined.html"
    else:
        file_name = output_path + f"Alignment/reports/{id}_{name_alignment}_report_alignment.html"

    with open(file_name, 'w') as file:
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


def grid_deform_detJacobian(simg2, transformDomainMeshSize, spline_order, outTx_Bspline, i, axs, thresh_pixelSIMG2, simg2_tr):
    grid_spacing = simg2.shape[1] / 20, simg2.shape[0] / 20
    grid = sitk.GridSource(
        outputPixelType=sitk.sitkFloat32,
        size=(simg2.shape[1], simg2.shape[0]),
        sigma=(0.5, 0.5),
        gridSpacing=(grid_spacing[0], grid_spacing[1]),
        gridOffset=(0, 0),
        spacing=(1, 1),
    )

    array = sitk.GetArrayViewFromImage(grid)
    mesh_size = np.array(transformDomainMeshSize)
    ctrl_pts = np.array(mesh_size) + spline_order

    transform = sitk.BSplineTransformInitializer(grid, mesh_size.tolist(), spline_order)
    transform.SetParameters(outTx_Bspline.GetParameters())

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(grid)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(100)
    resampler.SetOutputPixelType(sitk.sitkFloat32)
    resampled = resampler.Execute(grid)

    array_contour = sitk.GetArrayViewFromImage(resampled)

    #
    x_coeff, y_coeff = outTx_Bspline.GetCoefficientImages()
    grid_origin = x_coeff.GetOrigin()
    grid_spacing = x_coeff.GetSpacing()
    x = np.linspace(grid_origin[0], grid_origin[0] + (ctrl_pts[0] - 1) * grid_spacing[0], ctrl_pts[0])
    y = np.linspace(grid_origin[1], grid_origin[1] + (ctrl_pts[1] - 1) * grid_spacing[1], ctrl_pts[1])
    xx, yy = np.meshgrid(x, y)

    uv = np.reshape(np.array(transform.GetParameters()), (ctrl_pts[0], ctrl_pts[1], 2), order='F')
    u, v = uv[..., 0].T, uv[..., 1].T
    y_values = yy - v
    x_values = xx - u
    result_array = np.dstack((x_values, y_values))

    displacement_field = sitk.TransformToDisplacementField(transform, sitk.sitkVectorFloat64,
                                                           grid.GetSize(), grid.GetOrigin(),
                                                           grid.GetSpacing(), grid.GetDirection())
    displacement_field_np_arr = sitk.GetArrayViewFromImage(displacement_field)
    jacobian_det_volume = sitk.DisplacementFieldJacobianDeterminant(displacement_field)
    jacobian_det_np_arr = sitk.GetArrayViewFromImage(jacobian_det_volume)

    ## Plot
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

    return axs


def correlation_rasters_report(outTx_Bspline_dict, outTx_Rigid, scale_percent1, scale_percent2, coords1, coords2, label1, label2, pixel_size_raster_fullres, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, img1_resize, img2_resize, img1_resize_ori, img2_resize_ori):
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

        coords1_tr = transform_coords(coords1, outTx_Bspline, outTx_Rigid, scale_percent1, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, img1_resize, img2_resize, img1_resize_ori, img2_resize_ori)
    
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
    
def transform_coords(coords1, outTx_Bspline, outTx_Rigid, scale_percent1, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, img1_resize, img2_resize, img1_resize_ori, img2_resize_ori):
    scaled_points1 = coords1*scale_percent1

    '''
    # Load and display the image
    plt.imshow(img1_resize_ori, cmap="Blues")

    # Plot points with increased size (adjust "s" for desired size)
    plt.scatter(scaled_points1[:, 0], scaled_points1[:, 1], 
                color='blue', marker='o', s=0.001, label="Scaled Points")  # Increased size
    plt.show()
    '''

    # Shift coords if cropping has been done on the first image
    scaled_points1 = shift_coordinates(scaled_points1, crop_coords1)
    # Shift and rotate coords based on pre-alignment
    scaled_points1 = rotate_coordinates_angle(scaled_points1, -manual_alignment_rotation, manual_alignment_rotation_shape)
    scaled_points1 = shift_coordinates(scaled_points1, tuple(np.negative(manual_alignment_displacement)))

    '''
    # Load and display the image
    plt.imshow(img1_resize, cmap="Blues")
    # Plot points with increased size (adjust "s" for desired size)
    plt.scatter(scaled_points1[:, 0], scaled_points1[:, 1], 
                color='blue', marker='o', s=0.001, label="Scaled Points")  # Increased size
    plt.show()
    '''

    
    # Transform scaled points 1 in the coordinate system of slide 2
    scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points1])
    scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])


    '''
    # Load and display the image
    plt.imshow(img2_resize, cmap="Reds")
    # Plot points with increased size (adjust "s" for desired size)
    plt.scatter(scaled_pointsTR[:, 0], scaled_pointsTR[:, 1], 
                color='blue', marker='o', s=0.001, label="Scaled Points")  # Increased size
    plt.show()
    '''

    # Shift coords back if cropping has been done on the second image
    scaled_pointsTR = shift_back_coordinates(scaled_pointsTR, crop_coords2)

    '''
    # Load and display the image
    plt.imshow(img2_resize_ori, cmap="Reds")
    # Plot points with increased size (adjust "s" for desired size)
    plt.scatter(scaled_pointsTR[:, 0], scaled_pointsTR[:, 1], 
                color='blue', marker='o', s=0.001, label="Scaled Points")  # Increased size
    plt.show()
    '''
    
    # Scale back to the full size image
    coords1_tr = scaled_pointsTR/scale_percent2

    

    return coords1_tr

def shift_coordinates(coords, crop_coords):
    """Applies the shift induced by cropping to coordinates as is done to the image."""
    shift_x, shift_y = crop_coords  
    shifted_coords = np.array([[point[0] - shift_x, point[1] - shift_y] for point in coords])

    return shifted_coords

def shift_back_coordinates(coords, crop_coords):
    """Applies the INVERSE of the shift induced by cropping to coordinates to get back to the original image."""
    shift_x, shift_y = crop_coords  
    shifted_coords = np.array([[point[0] + shift_x, point[1] + shift_y] for point in coords])

    return shifted_coords

def rotate_coordinates_angle(coordinates, angle, img_shape):
    """
    Rotate coordinates to match the image rotation with reshape=False
    """
    # Convert angle to radians
    theta = np.radians(angle)
    
    # Rotation matrix for 2D points
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                                [np.sin(theta), np.cos(theta)]])
    
    # Calculate the center of the image
    center_x, center_y = img_shape[1] / 2, img_shape[0] / 2
    
    # Translate coordinates to the image center
    coords_centered = coordinates - np.array([center_x, center_y])
    
    # Apply the rotation matrix to the centered coordinates
    rotated_coords_centered = np.dot(coords_centered, rotation_matrix.T)
    
    # Translate the coordinates back by adding the center back
    rotated_coords = rotated_coords_centered + np.array([center_x, center_y])
    
    return rotated_coords

def shift_and_rotate_coordinates(coords, crop_coords, angle):
    """Applies the shift induced by cropping and rotates coordinates around the origin."""
    shift_x, shift_y = crop_coords  
    shifted_coords = np.array([[point[0] - shift_x, point[1] - shift_y] for point in coords])
    
    # Convert angle to radians
    theta = np.radians(angle)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    
    # Rotation matrix
    rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])
    
    # Apply rotation
    rotated_coords = np.dot(shifted_coords, rotation_matrix.T)
    
    return rotated_coords


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
            crop_coords1 = metadata_images[name_alignment][f'crop_coords_{panel}']
            img_resize1 = img_resize_alignment[panel]
            manual_alignment_displacement = metadata_images[name_alignment]['manual_alignment_displacement']
            manual_alignment_rotation = metadata_images[name_alignment]['manual_alignment_rotation']
            manual_alignment_rotation_shape = metadata_images[name_alignment]['manual_alignment_rotation_shape']
            image_shape_manual_alignment = metadata_images[name_alignment]['image_shape_manual_alignment']
            ## Scale annotations to downscaled image resolution
            artefacts_empty_array = get_gdf(artefacts_empty_gdf_poly, scale_percent1, crop_coords1, image_shape1, img_resize1, c="B")
            analysis_area_array = get_gdf(analysis_area_gdf_poly, scale_percent1, crop_coords1, image_shape1, img_resize1, c="B")
            ## Transform scaled annotations to downscaled reference image
            scale_percent2 = metadata_images[name_alignment][f'scale_percent_{reference_panel}']
            image_shape2 = metadata_images[name_alignment][f'image_shape_{reference_panel}']
            crop_coords2 = metadata_images[name_alignment][f'crop_coords_{reference_panel}']
            img_resize2 = img_resize_alignment[reference_panel]
            outTx_Rigid = outTx_Rigid_alignment[panel]
            outTx_Bspline = outTx_Bspline_alignment[panel]
            artefacts_empty_alignment[f'{panel}_in_{reference_panel}_panel'] = transform_annotation(artefacts_empty_array, scale_percent1, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, outTx_Rigid, outTx_Bspline, img_resize2)
            analysis_area_alignment[f'{panel}_in_{reference_panel}_panel'] = transform_annotation(analysis_area_array, scale_percent1, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, outTx_Rigid, outTx_Bspline, img_resize2)


    ## Fix invalid geometries in artefacts and empty annotations
    panels_tr = [key for key in analysis_area_alignment.keys() if "in" in key]
    panels_tr.append(reference_panel)
    analysis_area_alignment_list = []
    artefacts_empty_alignment_list = []
    for panel in panels_tr:
        analysis_area_alignment_list.append(analysis_area_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))
        artefacts_empty_alignment_list.append(artefacts_empty_alignment[panel].apply(lambda geom: geom if geom.is_valid else geom.buffer(0)))
    '''
    fig, ax = plt.subplots(figsize=(6, 6)) 
    analysis_area_alignment_list[0].plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    analysis_area_alignment_list[1].plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    analysis_area_alignment_list[2].plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
    ax.axis('off')
    ax.invert_yaxis()
    plt.show()
    '''
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

    # Check if intersections_polygons is empty
    if intersections_polygons.empty:
        print("No intersections found between the analysis areas of each panel. One of the alignment might be suboptimal.")
        return None  # Indicate that this iteration should be skipped

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
    if inner_boundary:
        mask = outer_boundary.difference(inner_boundary)
    else:
        print("No empty or artefact area.")
        mask = outer_boundary

    '''
    ## Check if there is an analysis area annotation
    if analysis_area_alignment_list:
    else:
        ## Merge the remaining GeoSeries
        merged_polygons = gpd.GeoSeries(unary_union([g.unary_union for g in artefacts_empty_alignment_list]))

        fig, ax = plt.subplots(figsize=(6, 6))
        merged_polygons.plot(ax=ax, color='blue', alpha=0.5, aspect='equal')
        ax.axis('off')
        ax.invert_yaxis()
        plt.show()

        # Define an outer boundary with the size of the reference image
        img_fullres_shape2 = tuple(int(dim / scale_percent2) for dim in image_shape2)

        height, width = img_fullres_shape2
        outer_boundary = gpd.GeoSeries(box(0, 0, width, height)).unary_union

        # Create a Polygon representing the inner boundary
        inner_boundary = merged_polygons.unary_union

        # Construct the mask by taking the difference between the outer boundary and inner boundary
        if inner_boundary:
            mask = outer_boundary.difference(inner_boundary)
        else:
            print("No empty or artefact area.")
            mask = outer_boundary
    '''
    
    # Plot the MultiPolygon
    # plot_multipolygon(mask)
    
    ## Convert mask to micro meters
    scaled_mask = scale_multipolygon_coordinates(mask, resolution_micron)
    geojson_dict = mapping(scaled_mask)
    with open(output_path + f"Alignment/merged_annotations/{id}_mask.geojson", "w") as f:
        json.dump(geojson_dict, f, indent=4)

    return mask



def get_gdf(gdf, scale_percent, crop_coords, image_shape, img_resize = None, c = "R"):
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
        # Shift coords if cropping has been done on the first image
        scaled_points = shift_coordinates(scaled_points, crop_coords)
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



def transform_annotation(annotations, scale_percent1, scale_percent2, image_shape1, image_shape2, crop_coords1, crop_coords2, manual_alignment_displacement, manual_alignment_rotation, manual_alignment_rotation_shape, outTx_Rigid, outTx_Bspline, img_resize):
    
    # Transform Shapes from IF2 to IF3 coordinate system
    annotation_TR= []
    for shape in annotations:

        scaled_points = shape*scale_percent1
        # Shift coords if cropping has been done on the first image
        scaled_points = shift_coordinates(scaled_points, crop_coords1)
        # Shift and rotate coords based on pre-alignment
        scaled_points = rotate_coordinates_angle(scaled_points, -manual_alignment_rotation, manual_alignment_rotation_shape)
        scaled_points = shift_coordinates(scaled_points, tuple(np.negative(manual_alignment_displacement)))

        
        # Transform scaled points 1 in the coordinate system of slide 2
        scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points])
        scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])

        # Shift coords back if cropping has been done on the second image
        scaled_pointsTR = shift_back_coordinates(scaled_pointsTR, crop_coords2)
        # Scale back to the full size image
        annotation_TR.append(scaled_pointsTR/scale_percent2)

    
    '''
    # Convert to polygon object
    polygonsTR = []
    for array in annotation_TR:
        scaled_points = array*scale_percent2
        # Shift coords if cropping has been done on the first image
        scaled_points = shift_coordinates(scaled_points, crop_coords2)
        # Convert array to list of tuples
        points = [(x, y) for x, y in scaled_points]
        # Create Polygon object
        polygonsTR.append(Polygon(points))

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
    crop_coords1 = metadata_images[name_alignment][f'crop_coords_{panels_alignment[0]}']
    manual_alignment_displacement = metadata_images[name_alignment]['manual_alignment_displacement']
    manual_alignment_rotation = metadata_images[name_alignment]['manual_alignment_rotation']
    manual_alignment_rotation_shape = metadata_images[name_alignment]['manual_alignment_rotation_shape']
    image_shape_manual_alignment = metadata_images[name_alignment]['image_shape_manual_alignment']
    

    scale_percent2 = metadata_images[name_alignment][f'scale_percent_{panels_alignment[1]}']
    image_shape2 = metadata_images[name_alignment][f'image_shape_{panels_alignment[1]}']
    crop_coords2 = metadata_images[name_alignment][f'crop_coords_{panels_alignment[1]}']

    outTx_Rigid = outTx_Rigid_alignment[panels_alignment[0]]
    outTx_Bspline = outTx_Bspline_alignment[panels_alignment[0]]

    coords = cell_coordinates[f'{panels_alignment[0]}_panel_DAPI']
    ids = np.arange(len(coords))


    scaled_points = coords*scale_percent1
    # Shift coords if cropping has been done on the first image
    scaled_points = shift_coordinates(scaled_points, crop_coords1)
    # Shift and rotate coords based on pre-alignment
    scaled_points = rotate_coordinates_angle(scaled_points, -manual_alignment_rotation, manual_alignment_rotation_shape)
    scaled_points = shift_coordinates(scaled_points, tuple(np.negative(manual_alignment_displacement)))

    # Transform scaled points 1 in the coordinate system of slide 2
    scaled_pointsTR_tmp = np.array([list(outTx_Bspline.TransformPoint(point)) for point in scaled_points])
    scaled_pointsTR = np.array([list(outTx_Rigid.TransformPoint(point)) for point in scaled_pointsTR_tmp])
    
    # Shift coords back if cropping has been done on the second image
    scaled_pointsTR = shift_back_coordinates(scaled_pointsTR, crop_coords2)
    # Scale back to the full size image
    coords_tr = scaled_pointsTR/scale_percent2
    cell_coordinates[f'{panels_alignment[0]}_panel_DAPI_transformed'] = coords_tr

    if mask == "No_annotations":
        filtered_coords_tr = coords_tr
        filtered_ids = ids
    else:
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
    required_columns = ['Phenotype', 'Object.Id', 'Classifier.Label']
    required_columns_SPIAT = ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']
    required_columns_ROI = ['x_', 'y_', 'x', 'y']
    required_columns_HALO = ['XMin', 'XMax', 'YMin', 'YMax']
    
    # Check if all required columns are present in df
    if all(col in df_cells.columns for col in required_columns):
        filtered_df = df_cells.loc[filtered_ids, ['Phenotype', 'Object.Id', 'Classifier.Label']].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords_tr / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)', 'Phenotype', 'Classifier.Label', 'Object.Id']]

    elif all(col in df_cells.columns for col in required_columns_SPIAT):
        filtered_df = df_cells.loc[filtered_ids, ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords_tr / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)', 'Cell_type', 'Phenotype', 'Classifier Label', 'Object Id']]

    elif all(col in df_cells.columns for col in required_columns_ROI):
        # Get all columns in df_cells except for 'x' and 'y'
        columns_without_xy = [col for col in df_cells.columns if col not in ['x', 'y', 'x (micron)', 'y (micron)', 'XMin_pixel', 'XMax_pixel', 'YMin_pixel', 'YMax_pixel']]
        filtered_df = df_cells.loc[filtered_ids, columns_without_xy].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords_tr / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)'] + columns_without_xy]

    elif all(col in df_cells.columns for col in required_columns_HALO):
        # Get all columns in df_cells except for 'x' and 'y'
        columns_without_xy = [col for col in df_cells.columns if col not in ['x', 'y', 'XMin', 'XMax', 'YMin', 'YMax']]
        filtered_df = df_cells.loc[filtered_ids, columns_without_xy].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords_tr / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)'] + columns_without_xy]
    filtered_df['Panel'] = panels_alignment[0]
    merged_cell_coordinates = pd.concat([merged_cell_coordinates, filtered_df], ignore_index=True)

    
    return merged_cell_coordinates


def filter_coordinates(cell_coordinates, panels_alignment, mask, data_frame_cells, resolution_micron, merged_cell_coordinates):
    print("----------------------")
    print(panels_alignment[1])
    coords = cell_coordinates[f'{panels_alignment[1]}_panel_DAPI']
    # Generate IDs
    ids = np.arange(len(coords))

    
    if mask == "No_annotations":
        filtered_coords = coords
        filtered_ids = ids
    else:
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
    required_columns = ['Phenotype', 'Object.Id', 'Classifier.Label']
    required_columns_SPIAT = ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']
    required_columns_ROI = ['x_', 'y_', 'x', 'y']
    required_columns_HALO = ['XMin', 'XMax', 'YMin', 'YMax']

    # Check if all required columns are present in df
    if all(col in df_cells.columns for col in required_columns):
        filtered_df = df_cells.loc[filtered_ids, ['Phenotype', 'Object.Id', 'Classifier.Label']].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)', 'Phenotype', 'Classifier.Label', 'Object.Id']]

    elif all(col in df_cells.columns for col in required_columns_SPIAT):
        filtered_df = df_cells.loc[filtered_ids, ['Cell_type', 'Phenotype', 'Object Id', 'Classifier Label']].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)', 'Cell_type', 'Phenotype', 'Classifier Label', 'Object Id']]

    elif all(col in df_cells.columns for col in required_columns_ROI):
        # Get all columns in df_cells except for 'x' and 'y'
        columns_without_xy = [col for col in df_cells.columns if col not in ['x', 'y', 'x (micron)', 'y (micron)', 'XMin_pixel', 'XMax_pixel', 'YMin_pixel', 'YMax_pixel']]
        filtered_df = df_cells.loc[filtered_ids, columns_without_xy].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)'] + columns_without_xy]

    elif all(col in df_cells.columns for col in required_columns_HALO):
        # Get all columns in df_cells except for 'x' and 'y'
        columns_without_xy = [col for col in df_cells.columns if col not in ['x', 'y', 'XMin', 'XMax', 'YMin', 'YMax']]
        filtered_df = df_cells.loc[filtered_ids, columns_without_xy].copy()
        filtered_df[['x (micron)', 'y (micron)']] = filtered_coords / resolution_micron
        filtered_df = filtered_df[['x (micron)', 'y (micron)'] + columns_without_xy]

    filtered_df['Panel'] = panels_alignment[1]
    merged_cell_coordinates = pd.concat([merged_cell_coordinates.loc[:,~merged_cell_coordinates.columns.duplicated()].reset_index(drop=True), filtered_df.loc[:,~filtered_df.columns.duplicated()].reset_index(drop=True)], ignore_index=True)
    # merged_cell_coordinates = pd.concat([merged_cell_coordinates, filtered_df], ignore_index=True)

    return merged_cell_coordinates


def save_tables(merged_cell_coordinates, output_path, id):
    if 'Cell_type' in merged_cell_coordinates.columns:
        merged_cell_coordinates["label"] = merged_cell_coordinates.index
        merged_cell_coordinates.rename(columns={"Classifier Label": "Classifier_Label"}, inplace=True)
        # merged_cell_coordinates.drop(columns=["Object Id"], inplace=True)
        print("Combined phenotypes:")
        print(np.unique(merged_cell_coordinates['Cell_type']))
    # put panel after transformed coordinates for readability with HALO files
    panel = merged_cell_coordinates.pop('Panel')
    merged_cell_coordinates.insert(2, 'Panel', panel)

    merged_cell_coordinates.to_csv(output_path + f"Alignment/merged_tables/{id}_merged_cell_coordinates.csv", index=False)



def plot_rasters(data_frame_cells, merged_cell_coordinates, cell_coordinates, pixel_size_raster_micron, output_path, panels_all, id, resolution_micron):
    '''
    ## Plot the raster before alignment
    if len(panels_all)==3:
        coords1 = data_frame_cells[f'{panels_all[0]}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
        coords2 = data_frame_cells[f'{panels_all[1]}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
        coords3 = data_frame_cells[f'{panels_all[2]}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
    if len(panels_all)==2:
        coords1 = data_frame_cells[f'{panels_all[0]}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
        coords2 = data_frame_cells[f'{panels_all[1]}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
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
    '''

    ## Plot the raster before alignment
    # Extract coordinates
    coords = [
        data_frame_cells[f'{panel}_panel_df'][["x (micron)", "y (micron)"]].to_numpy()
        for panel in panels_all
    ]
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords_set[:, 0]) for coords_set in coords) - 1000
    max_x = max(np.max(coords_set[:, 0]) for coords_set in coords) + 1000
    min_y = min(np.min(coords_set[:, 1]) for coords_set in coords) - 1000
    max_y = max(np.max(coords_set[:, 1]) for coords_set in coords) + 1000
    # Create rasters with a fixed size
    rasters = [
        create_raster_from_points(
            coords_set, 
            pixel_size=pixel_size_raster_micron, 
            min_x=min_x, 
            max_x=max_x, 
            min_y=min_y, 
            max_y=max_y
        )
        for coords_set in coords
    ]
    # Compute the pearson correlation between rasters
    correlations = []
    correlation_texts = []
    for i in range(len(rasters)):
        for j in range(i + 1, len(rasters)):
            corr = pearsonr(rasters[i].flatten(), rasters[j].flatten())[0]
            correlations.append(corr)
            correlation_texts.append(f"Correlation {panels_all[i]}-{panels_all[j]}: {corr:.2f}")
    # Plot the rasters
    # Define colors and labels dynamically based on panels_all
    colors = ["blue", "green", "red"][:len(panels_all)]
    labels = list(map("DAPI {}".format, panels_all))
    plt.figure(figsize=(4, 4))
    for i, raster in enumerate(rasters):
        cmap = plt.cm.colors.ListedColormap(["white", colors[i]])
        plt.imshow(raster, cmap=cmap, extent=[min_x, max_x, min_y, max_y], alpha=0.5, aspect='auto')
    legend_patches = [
        mpatches.Patch(color=colors[i], label=labels[i])
        for i in range(len(rasters))
    ]
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.legend(handles=legend_patches, loc="upper right")
    # Add correlations values
    for idx, text in enumerate(correlation_texts):
        plt.figtext(0.5, 0.05 + 0.05 * idx, text, fontsize=12, wrap=True, horizontalalignment="center")
    # Save the plot
    plt.savefig(output_path + f"Alignment/plots/{id}_raster_before_alignment.png", dpi=300, bbox_inches="tight")
    plt.close()


    '''
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
    '''


    ## Plot the common raster after alignment
    coords = [
        cell_coordinates[f'{panel}_panel_DAPI_transformed'] / resolution_micron if i % 2 == 0 
        else cell_coordinates[f'{panel}_panel_DAPI'] / resolution_micron
        for i, panel in enumerate(panels_all)
    ]
    # Define the range of coordinates to cover both sets of points
    all_coords = [coords[i] for i in range(len(coords))]
    min_x = min(np.min(c[:, 0]) for c in all_coords) - 1000
    max_x = max(np.max(c[:, 0]) for c in all_coords) + 1000
    min_y = min(np.min(c[:, 1]) for c in all_coords) - 1000
    max_y = max(np.max(c[:, 1]) for c in all_coords) + 1000
    # Create rasters with a fixed size  
    rasters = [
        create_raster_from_points(c, pixel_size=pixel_size_raster_micron, 
                                min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
        for c in all_coords
    ]
    # Compute the pearson correlation between rasters
    correlations = []
    correlation_texts = []
    for i in range(len(rasters)):
        for j in range(i + 1, len(rasters)):
            corr = pearsonr(rasters[i].flatten(), rasters[j].flatten())[0]
            correlations.append(corr)
            correlation_texts.append(f"Correlation {panels_all[i]}-{panels_all[j]}: {corr:.2f}")
    # Create a new array with the same shape initialized to zeros
    new_array = np.zeros_like(rasters[0], dtype=np.uint8)
    # Apply the rules to create the new array
    new_array[np.logical_and.reduce([(r == 1) for r in rasters])] = 3
    for i, r in enumerate(rasters):
        other_rasters = [r for j, r in enumerate(rasters) if j != i]
        new_array[np.logical_and.reduce([r == 1 for r in other_rasters]) & (r == 0)] = 1
    # Create an RGB image with the same height and width as new_array
    rgb_image = np.zeros((*new_array.shape, 3), dtype=np.uint8)
    colors = [[255, 255, 255], [255, 0, 0], [255, 165, 0], [0, 255, 0]]
    for i in range(len(colors)):
        rgb_image[new_array == i] = colors[i]
    colors_normalized = [[r / 255, g / 255, b / 255] for [r, g, b] in colors]
    legends = [
        mpatches.Patch(color=colors_normalized[1], label='Only one'),
        mpatches.Patch(color=colors_normalized[2], label='Part common'),
        mpatches.Patch(color=colors_normalized[3], label='Common')
    ]
    # Plot the RGB image
    plt.imshow(rgb_image)
    plt.axis('off')  # Turn off the axis
    plt.gca().invert_yaxis()
    plt.legend(handles=legends, loc='upper right')
    # Add correlations values
    for idx, text in enumerate(correlation_texts):
        plt.figtext(0.5, 0.05 + 0.05 * idx, text, fontsize=12, wrap=True, horizontalalignment="center")
    plt.savefig(output_path + f"Alignment/plots/{id}_common_raster.png", dpi=300, bbox_inches='tight')
    plt.close()
        


    '''
    ## Plot the common raster after alignment and filtering with common annotations
    coords1 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[0]}'][["x (micron)", "y (micron)"]].to_numpy()
    coords2 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[1]}'][["x (micron)", "y (micron)"]].to_numpy()
    coords3 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[2]}'][["x (micron)", "y (micron)"]].to_numpy()

    
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
    '''

    ## Plot the common raster after alignment and filtering with common annotations
    coords = [
        merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panel}'][["x (micron)", "y (micron)"]].to_numpy()
        for panel in panels_all
    ]
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(c[:, 0]) for c in coords) - 1000
    max_x = max(np.max(c[:, 0]) for c in coords) + 1000
    min_y = min(np.min(c[:, 1]) for c in coords) - 1000
    max_y = max(np.max(c[:, 1]) for c in coords) + 1000
    # Create rasters with a fixed size
    rasters = [
        create_raster_from_points(c, pixel_size=pixel_size_raster_micron, 
                                min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)
        for c in coords
    ]
    # Compute the Pearson correlation between rasters
    correlations = []
    correlation_texts = []
    for i in range(len(rasters)):
        for j in range(i + 1, len(rasters)):
            corr = pearsonr(rasters[i].flatten(), rasters[j].flatten())[0]
            correlations.append(corr)
            correlation_texts.append(f"Correlation {panels_all[i]}-{panels_all[j]}: {corr:.2f}")
    # Create a new array with the same shape initialized to zeros
    new_array = np.zeros_like(rasters[0], dtype=np.uint8)
    # Apply the rules to create the new array
    new_array[np.logical_and.reduce([r == 1 for r in rasters])] = 3
    for i, r in enumerate(rasters):
        other_rasters = [r for j, r in enumerate(rasters) if j != i]
        new_array[np.logical_and.reduce([r == 1 for r in other_rasters]) & (r == 0)] = 1
    # Create an RGB image with the same height and width as new_array
    rgb_image = np.zeros((*new_array.shape, 3), dtype=np.uint8)
    # Assign specific colors based on the pixel values
    colors = [[255, 255, 255], [255, 0, 0], [255, 165, 0], [0, 255, 0]]
    for i in range(len(colors)):
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
    # Add correlations values
    for idx, text in enumerate(correlation_texts):
        plt.figtext(0.5, 0.05 + 0.05 * idx, text, fontsize=12, wrap=True, horizontalalignment="center")
    plt.savefig(output_path + f"Alignment/plots/{id}_common_raster_filtered_annotations.png", dpi=300, bbox_inches='tight')
    plt.close()

    
    '''
    ## Plot the raster after alignment
    coords1_tr = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[0]}'][["x (micron)", "y (micron)"]].to_numpy()
    coords2 = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[1]}'][["x (micron)", "y (micron)"]].to_numpy()
    coords3_tr = merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panels_all[2]}'][["x (micron)", "y (micron)"]].to_numpy()
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
    '''


    ## Plot the raster before alignment
    # Extract coordinates
    coords = [
        merged_cell_coordinates[merged_cell_coordinates['Panel'] == f'{panel}'][["x (micron)", "y (micron)"]].to_numpy()
        for panel in panels_all
    ]
    # Define the range of coordinates to cover both sets of points
    min_x = min(np.min(coords_set[:, 0]) for coords_set in coords) - 1000
    max_x = max(np.max(coords_set[:, 0]) for coords_set in coords) + 1000
    min_y = min(np.min(coords_set[:, 1]) for coords_set in coords) - 1000
    max_y = max(np.max(coords_set[:, 1]) for coords_set in coords) + 1000
    # Create rasters with a fixed size
    rasters = [
        create_raster_from_points(
            coords_set, 
            pixel_size=pixel_size_raster_micron, 
            min_x=min_x, 
            max_x=max_x, 
            min_y=min_y, 
            max_y=max_y
        )
        for coords_set in coords
    ]
    # Compute the pearson correlation between rasters
    correlations = []
    correlation_texts = []
    for i in range(len(rasters)):
        for j in range(i + 1, len(rasters)):
            corr = pearsonr(rasters[i].flatten(), rasters[j].flatten())[0]
            correlations.append(corr)
            correlation_texts.append(f"Correlation {panels_all[i]}-{panels_all[j]}: {corr:.2f}")
    # Plot the rasters
    # Define colors and labels dynamically based on panels_all
    colors = ["blue", "green", "red"][:len(panels_all)]
    labels = list(map("DAPI {}".format, panels_all))
    plt.figure(figsize=(4, 4))
    for i, raster in enumerate(rasters):
        cmap = plt.cm.colors.ListedColormap(["white", colors[i]])
        plt.imshow(raster, cmap=cmap, extent=[min_x, max_x, min_y, max_y], alpha=0.5, aspect='auto')
    legend_patches = [
        mpatches.Patch(color=colors[i], label=labels[i])
        for i in range(len(rasters))
    ]
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.legend(handles=legend_patches, loc="upper right")
    # Add correlations values
    for idx, text in enumerate(correlation_texts):
        plt.figtext(0.5, 0.05 + 0.05 * idx, text, fontsize=12, wrap=True, horizontalalignment="center")
    # Save the plot
    plt.savefig(output_path + f"Alignment/plots/{id}_raster_alignment.png", dpi=300, bbox_inches="tight")
    plt.close()




