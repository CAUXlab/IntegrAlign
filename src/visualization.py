import re 
import os
from tqdm import tqdm # type: ignore

from tifffile import TiffFile # type: ignore
from tifffile.tiffcomment import tiffcomment
from xml.etree import ElementTree
import cv2 # type: ignore

import matplotlib.pyplot as plt # type: ignore
from io import BytesIO
import base64
import json
import numpy as np
import scipy.ndimage
from collections import Counter
import geopandas as gpd


from src.alignment import get_annotations, get_gdf, get_gdf_from_annot

def visualization(scans_paths, annotations, panels, output_path):
    """Generate a visualization report of all the slides for each patients."""
    ## Get the ids of every patients
    ids_panels = []
    for folder_scan, panel in zip(scans_paths, panels):
        ids = find_ids(folder_scan)
        # print(f"There is {len(ids)} patients in the panel {panel}")
        ids_panels.append(ids)

    ## Sort the common ids between the panels
    common_ids = set.intersection(*map(set, ids_panels))
    # sorted_common_ids = sorted(common_ids, key=lambda x: int(x))
    sorted_common_ids = sorted(common_ids, key=lambda x: (x.isdigit(), int(x) if x.isdigit() else x))
    print(f"There is {len(sorted_common_ids)} patients that are common in every panel")

    ## Get the correct annotations names 
    print("Reading annotation files...")
    if annotations:
        all_annotation_names = set()
        for annotations_panel in annotations:
            files = [f for f in os.listdir(annotations_panel) if (f.endswith(".annotations") or f.endswith(".geojson")) and not f.startswith(".")]
            for id in sorted_common_ids:
                annotation_file_path = next(
                    (os.path.join(annotations_panel, f) for f in files if id in f),
                    None
                )
                annotation_names = get_annotation_classification_names(annotation_file_path)
                all_annotation_names.update(annotation_names)

        all_annotation_names = sorted(all_annotation_names)

        print("\nFound the following unique annotation names:")
        for idx, name in enumerate(all_annotation_names):
            print(f"{idx + 1}: {name}")

        # Prompt user to choose corresponding categories by index
        analysis_area_indices = input("\nEnter numbers (comma-separated) for 'AnalysisArea' annotations: ")
        artefact_area_indices = input("Enter numbers (comma-separated) for 'Artefacts' annotations: ")
        empty_area_indices = input("Enter numbers (comma-separated) for 'Empty' annotations: ")

        # Convert input to selected names
        def get_selected_names(indices_str):
            indices = [int(i.strip()) - 1 for i in indices_str.split(",") if i.strip().isdigit()]
            return [all_annotation_names[i] for i in indices if 0 <= i < len(all_annotation_names)]

        annotations_names_AnalysisArea = get_selected_names(analysis_area_indices)
        annotations_names_Artefacts = get_selected_names(artefact_area_indices)
        annotations_names_Empty = get_selected_names(empty_area_indices)

        print("\n✅ Selected AnalysisArea names:", annotations_names_AnalysisArea)
        print("✅ Selected Artefacts names:", annotations_names_Artefacts)
        print("✅ Selected Empty names:", annotations_names_Empty)
    else:
        annotations_names_Empty = []
        annotations_names_Artefacts = []
        annotations_names_AnalysisArea = []
    



    ## Save the report with the compressed images of every common ids
    figs = []
    for id in tqdm(sorted_common_ids, desc="Processing images", unit="ID"):
        fig = plot_panels(scans_paths, annotations, id, panels, annotations_names_Empty, annotations_names_Artefacts, annotations_names_AnalysisArea)  # Generate figure
        figs.append(fig)
    save_html_report(figs, sorted_common_ids, output_path + "images_report.html")

    ## Save parameters for next steps
    dict_params = {"scans_paths":scans_paths, "annotations_paths":annotations, 
                   "annotations_names_empty":annotations_names_Empty, "annotations_names_artefacts":annotations_names_Artefacts, "annotations_names_AnalysisArea":annotations_names_AnalysisArea, 
                   "common_ids":sorted_common_ids, "panels":panels, "output_path":output_path}
    params_json = json.dumps(dict_params,indent=4)
    with open(output_path + "params.json","w") as outfile:
        outfile.write(params_json)
    print(f'Params file saved to {output_path + "params.json"}')


def find_ids(folder_path):
    # id_regex = re.compile(r'^(\d+)')
    id_regex = re.compile(r'^[^_]+')
    
    ids = []
    for filename in os.listdir(folder_path):
        if filename.endswith(('.tif', '.qptiff')):
            match = id_regex.match(filename)
            if match:
                ids.append(match.group())
    
    return ids

import xml.etree.ElementTree as ET

def get_annotation_classification_names(annotation_file_path):
    """
    Parse an .annotations file and return all unique classification names.

    Parameters:
        annotation_file_path (str): Path to the .annotations XML file.

    Returns:
        List[str]: A list of unique classification names found in the file.
    """
    if annotation_file_path.endswith(".annotations"):
        # Parse the XML file
        tree = ET.parse(annotation_file_path)
        root = tree.getroot()

        # Convert XML to GeoDataFrame using external function
        gdf = get_gdf_from_annot(root)

    elif annotation_file_path.endswith(".geojson"):
        ## Read annotations file
        gdf = gpd.read_file(annotation_file_path)
        def extract_name(val):
            if isinstance(val, dict):
                return val.get("name")
            elif isinstance(val, str):
                try:
                    import json
                    data = json.loads(val)
                    return data.get("name")
                except Exception:
                    return None
            else:
                return None
        gdf['classification'] = gdf['classification'].apply(extract_name)
    else:
        raise ValueError("File must have a .annotations extension")

    # Extract and return unique classification names
    return sorted(set(name for name in gdf['classification'].dropna()))


def plot_panels(scans_paths, annotations_paths, id, panels, annotations_names_empty, annotations_names_artefacts, annotations_names_AnalysisArea):
    # print("------------------------")
    # print("Patient ", id)
    imgs_resized = []
    annotations_resized = []
    if annotations_paths:
        ## Get all compressed images and annotations
        for scan, annotations, panel in zip(scans_paths, annotations_paths, panels):
            # Define scan path
            path = next((os.path.join(scan, f) for f in os.listdir(scan) if id in f and f.endswith(('.tif', '.qptiff')) and not f.startswith('.')), None)
            # print(f"{panel}:", path)
            # Load image
            img_resized, scale_percent = load_compressed_img(path)
            # Scale images to 8bit
            img_resized = cv2.convertScaleAbs(img_resized)
            imgs_resized.append(img_resized)

            if annotations_paths:
                # Get annotations
                artefacts_empty_alignment = {}
                analysis_area_alignment = {}
                # geojson_file_path = next((os.path.join(annotations, f) for f in os.listdir(annotations) if id in f and f.endswith(".geojson") and not f.startswith(".")), None)
                annotation_file_path = next(
                    (os.path.join(annotations, f) for f in os.listdir(annotations) 
                    if id in f and (f.endswith(".annotations") or f.endswith(".geojson")) and not f.startswith(".")), 
                    None
                )

                artefacts_empty_alignment, analysis_area_alignment = get_annotations(annotation_file_path, panel, artefacts_empty_alignment, analysis_area_alignment, annotations_names_empty, annotations_names_artefacts, annotations_names_AnalysisArea)

                analysis_area_gdf_poly = analysis_area_alignment[panel]
                analysis_area_array = get_gdf(analysis_area_gdf_poly, scale_percent, crop_coords = (0, 0), image_shape = img_resized.shape, img_resize = img_resized)
                annotations_resized.append([area * scale_percent for area in analysis_area_array])
        ## Plot the compressed images
        fig, axes = plt.subplots(1, len(panels), figsize=(10, 3))
        for i, (img_resized, panel, annotations) in enumerate(zip(imgs_resized, panels, annotations_resized)):
            axes[i].imshow(img_resized, cmap="gray")  # Display the image
            axes[i].set_title(panel, fontsize=6)  # Set panel title
            axes[i].axis('off')  # Turn off axes
            if annotations_paths:
                # Plot the annotations
                for annotation_array in annotations:
                    axes[i].plot(annotation_array[:, 0], annotation_array[:, 1], color='red', linewidth=0.5)
        plt.tight_layout()
        plt.close(fig)
    else:
        ## Get all compressed images
        for scan, panel in zip(scans_paths, panels):
            # Define scan path
            path = next((os.path.join(scan, f) for f in os.listdir(scan) if id in f and f.endswith(('.tif', '.qptiff')) and not f.startswith('.')), None)
            # print(f"{panel}:", path)
            # Load image
            img_resized, scale_percent = load_compressed_img(path)
            # Scale images to 8bit
            img_resized = cv2.convertScaleAbs(img_resized)

            imgs_resized.append(img_resized)
            ## Plot the compressed images
        fig, axes = plt.subplots(1, len(panels), figsize=(10, 3))
        for i, (img_resized, panel) in enumerate(zip(imgs_resized, panels)):
            axes[i].imshow(img_resized, cmap="gray")
            axes[i].set_title(panel, fontsize=6)
            axes[i].axis('off')
        plt.tight_layout()
        plt.close(fig)
    
    return fig


def save_html_report(figs, ids, output_html="report.html"):
    """
    Save all matplotlib figures to an HTML report.
    """
    html_content = "<html><head><title>Image Report</title></head><body>"
    html_content += "<h1>Image Report</h1>"

    for fig, id in zip(figs, ids):
        # Save the figure to a BytesIO buffer
        buffer = BytesIO()
        fig.savefig(buffer, format="png", bbox_inches="tight")
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.read()).decode("utf-8")
        buffer.close()
        plt.close(fig)

        # Add image to the HTMLBytesIO
        html_content += f"<h2>Patient ID: {id}</h2>"
        html_content += f"<img src='data:image/png;base64,{img_base64}'><br>"

    html_content += "</body></html>"

    # Write to file
    with open(output_html, "w") as f:
        f.write(html_content)

    print(f"Report saved to {output_html}")


def load_compressed_img(path_scan):
    if path_scan.endswith('.qptiff'):
        print(path_scan)
        tif = TiffFile(path_scan)
    
        ## Load the most compressed DAPI image in .pages
        # get nb of channels
        last_indexes = 10
        im_sizes = []
        for im in tif.pages[-last_indexes:]:  
            im_sizes.append(len(im.asarray()))
            
        count_dict = Counter(im_sizes)  # Count occurrences of each size
        nb_channels = max(count_dict.values())  # Most frequent count
        most_frequent_size = max(count_dict, key=count_dict.get)  # Get the most frequent size
        # Find indices of channels
        indices = [i for i, size in enumerate(im_sizes) if size == most_frequent_size]
        indices  = [x - last_indexes for x in indices]
        index_DAPI = min(indices)
        img_resize = tif.pages[index_DAPI].asarray()
        # Get scale_percent
        tif_tags= {}
        for tag in tif.pages[0].tags.values():
            name, value = tag.name, tag.value
            tif_tags[name] = value
        scale_percent = img_resize.shape[0] / tif_tags['ImageLength']

    elif path_scan.endswith('.tif'):
        tif = TiffFile(path_scan, is_ome=False)
        # Load the most compressed DAPI image in .levels
        img_resize = tif.series[0].levels[-1].asarray()[0]
        # img_resized = rotate_array(img_resized, 180)  # Rotate
        # Get scale percent from full res image and compressed image sizes
        xml = tiffcomment(path_scan)
        root = ElementTree.fromstring(xml.replace("\n", "").replace("\t", ""))
        sizeX = list(dict.fromkeys([x.get("sizeX") for x in root.iter() if x.tag == "dimension"]))
        scale_percent = int(sizeX[-1])/int(sizeX[0])
    
    return img_resize, scale_percent



def rotate_array(array, angle):
    """
    Rotates a 2D NumPy array by a specified angle without modifying values due to interpolation.
    
    Parameters:
        array (numpy.ndarray): Input 2D array.
        angle (float): Rotation angle in degrees (from -180 to 180).
    
    Returns:
        numpy.ndarray: Rotated array.
    """
    return scipy.ndimage.rotate(array, angle, reshape=False, mode='nearest', order=0)



