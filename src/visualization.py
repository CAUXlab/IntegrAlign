import re 
import os
from tqdm import tqdm # type: ignore

from tifffile import TiffFile # type: ignore
import cv2 # type: ignore

import matplotlib.pyplot as plt # type: ignore
from io import BytesIO
import base64
import json
import numpy as np

from src.alignment import get_annotations, get_gdf

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

    ## Save the report with the compressed images of every common ids
    figs = []
    for id in tqdm(sorted_common_ids, desc="Processing images", unit="ID"):
        fig = plot_panels(scans_paths, annotations, id, panels)  # Generate figure
        figs.append(fig)
    save_html_report(figs, sorted_common_ids, output_path + "images_report.html")

    ## Save parameters for next steps
    dict_params = {"folder_paths":scans_paths, "common_ids":sorted_common_ids, "panels":panels, "output_path":output_path}
    params_json = json.dumps(dict_params,indent=4)
    with open(output_path + "params.json","w") as outfile:
        outfile.write(params_json)
    print(f'Params file saved to {output_path + "params.json"}')


def find_ids(folder_path):
    # id_regex = re.compile(r'^(\d+)')
    id_regex = re.compile(r'^[^_]+')
    
    ids = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.qptiff'):
            match = id_regex.match(filename)
            if match:
                ids.append(match.group())
    
    return ids


def plot_panels(scans_paths, annotations_paths, id, panels):
    # print("------------------------")
    # print("Patient ", id)
    imgs_resized = []
    annotations_resized = []
    if annotations_paths:
        ## Get all compressed images and annotations
        for scan, annotations, panel in zip(scans_paths, annotations_paths, panels):
            # Define scan path
            path = next((os.path.join(scan, f) for f in os.listdir(scan) if id in f and f.endswith(".qptiff")), None)
            # print(f"{panel}:", path)
            # Load image
            diff_index_pages = 10
            img_resized, scale_percent = load_compressed_img(path, diff_index_pages)
            # Scale images to 8bit
            img_resized = cv2.convertScaleAbs(img_resized)
            imgs_resized.append(img_resized)

            # Get annotations
            artefacts_empty_alignment = {}
            analysis_area_alignment = {}
            geojson_file_path = next((os.path.join(annotations, f) for f in os.listdir(annotations) if id in f and f.endswith(".geojson") and not f.startswith(".")), None)
            artefacts_empty_alignment, analysis_area_alignment = get_annotations(geojson_file_path, panel, artefacts_empty_alignment, analysis_area_alignment)

            analysis_area_gdf_poly = analysis_area_alignment[panel]
            analysis_area_array = get_gdf(analysis_area_gdf_poly, scale_percent, operation = 0, image_shape = img_resized.shape, img_resize = img_resized)
            annotations_resized.append([area * scale_percent for area in analysis_area_array])
        ## Plot the compressed images
        fig, axes = plt.subplots(1, len(panels), figsize=(10, 3))
        for i, (img_resized, panel, annotations) in enumerate(zip(imgs_resized, panels, annotations_resized)):
            axes[i].imshow(img_resized, cmap="gray")  # Display the image
            axes[i].set_title(panel, fontsize=6)  # Set panel title
            axes[i].axis('off')  # Turn off axes
            # Plot the annotations
            for annotation_array in annotations:
                axes[i].plot(annotation_array[:, 0], annotation_array[:, 1], color='red', linewidth=0.5)
        plt.tight_layout()
        plt.close(fig)
    else:
        ## Get all compressed images
        for scan, panel in zip(scans_paths, panels):
            # Define scan path
            path = next((os.path.join(scan, f) for f in os.listdir(scan) if id in f and f.endswith(".qptiff")), None)
            # print(f"{panel}:", path)
            # Load image
            diff_index_pages = 10
            img_resized, scale_percent = load_compressed_img(path, diff_index_pages)
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


def load_compressed_img(path_qptiff, diff_index_pages):
    tif = TiffFile(path_qptiff)
    # Load the most compressed DAPI image in .pages
    most_comp_DAPI_index = len(tif.pages) - diff_index_pages
    img_resized = tif.pages[most_comp_DAPI_index].asarray()

    tif_tags= {}
    for tag in tif.pages[0].tags.values():
        name, value = tag.name, tag.value
        tif_tags[name] = value
    scale_percent = img_resized.shape[0] / tif_tags['ImageLength']
    
    return img_resized, scale_percent
