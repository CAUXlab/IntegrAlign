import re 
import os
from tqdm import tqdm # type: ignore

from tifffile import TiffFile # type: ignore
import cv2 # type: ignore

import matplotlib.pyplot as plt # type: ignore
from io import BytesIO
import base64
import json

def visualization(scans_paths, panels, output_path):
    """Generate a visualization report of all the slides for each patients."""
    ## Get the ids of every patients
    ids_panels = []
    for folder_scan, panel in zip(scans_paths, panels):
        ids = find_ids(folder_scan)
        # print(f"There is {len(ids)} patients in the panel {panel}")
        ids_panels.append(ids)

    ## Sort the common ids between the panels
    common_ids = set.intersection(*map(set, ids_panels))
    sorted_common_ids = sorted(common_ids, key=lambda x: int(x))
    print(f"There is {len(sorted_common_ids)} patients that are common in every panel")

    ## Save the report with the compressed images of every common ids
    figs = []
    for id in tqdm(sorted_common_ids, desc="Processing images", unit="ID"):
        fig = plot_panels(scans_paths, id, panels)  # Generate figure
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


def plot_panels(folder_paths, id, panels):
    # print("------------------------")
    # print("Patient ", id)
    imgs_resized = []
    ## Get all compressed images
    for folder, panel in zip(folder_paths, panels):
        # Define scan path
        path = next((os.path.join(folder, f) for f in os.listdir(folder) if id in f and f.endswith(".qptiff")), None)
        # print(f"{panel}:", path)
        # Load image
        diff_index_pages = 10
        img_resized = load_compressed_img(path, diff_index_pages)
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
    
    return img_resized
