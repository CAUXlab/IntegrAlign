import json
import numpy as np
from tqdm import tqdm # type: ignore
import pickle

from tifffile import TiffFile # type: ignore
import cv2 # type: ignore
import SimpleITK as sitk # type: ignore

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill
from openpyxl.formatting.rule import FormulaRule


def save_downscaled_images(params_file_path, excluded_ids, rotation_params):
    ## Get the parameters for alignment
    folder_paths, common_ids, panels, output_path, turn_img_id_dict = get_parameters(params_file_path, rotation_params, excluded_ids)
    ## Get the name of the panels for each alignment (with reference panel at the end)
    data_dict = {}
    panel_alignment_dict = panels_name_alignment(panels, folder_paths)

    unique_name_alignments = set()
    for id in tqdm(common_ids, desc="Loading downscaled images", unit="Patient"):
        data_dict[id] = {}  # Initialize sub-dictionary for this patient
        for name_alignment, folder_alignment_paths in panel_alignment_dict.items():
            # Add the name_alignment to the set (sets automatically handle uniqueness)
            unique_name_alignments.add(name_alignment)
            folder_path1, folder_path2 = [path for path in folder_alignment_paths]
            # Load imgs downscaled and metadata
            (tif_tags1, channel_list1, channel_name_dictionary1, scale_percent1, 
            tif_tags2, channel_list2, channel_name_dictionary2, scale_percent2, 
            img1, img2, img1_resize, img2_resize) = load_downscaled_imgs(folder_path1, folder_path2, id, name_alignment, turn_img_id_dict)
            data_dict[id][name_alignment] = {
                "tif_tags1": tif_tags1,
                "channel_list1": channel_list1,
                "channel_name_dictionary1": channel_name_dictionary1,
                "scale_percent1": scale_percent1,
                "tif_tags2": tif_tags2,
                "channel_list2": channel_list2,
                "channel_name_dictionary2": channel_name_dictionary2,
                "scale_percent2": scale_percent2,
                "img1": img1,
                "img2": img2,
                "img1_resize": img1_resize,
                "img2_resize": img2_resize,
            }
    # Add params
    data_dict["params"] = {"common_ids" : common_ids, "panels" : panels, "output_path" : output_path, "turn_img_id_dict" : turn_img_id_dict}
    # Save the dictionary to a file
    with open(output_path + "downscaled_images.pkl", "wb") as file:
        pickle.dump(data_dict, file)
    print(f'Data saved to {output_path + "downscaled_images.pkl"}')

    ## Create DataFrame for common_ids with empty status columns for manual input
    # Initialize the data dictionary with 'Sample id' and 'Excluded' columns
    data = {
        'Sample id': common_ids,
        'Excluded': ['No'] * len(common_ids)
    }
    # Create a column for each unique name_alignment and add it to the data dictionary
    for name_alignment in unique_name_alignments:
        # Dynamically create a new column name for each unique 'name_alignment'
        column_name = f'Status {name_alignment} alignment'
        data[column_name] = [None] * len(common_ids)  # Initialize with None for manual input
    # Add the 'Comments' column after all the Status columns
    data['Comments'] = [None] * len(common_ids)
    df = pd.DataFrame(data)

    ## Create DataFrame for excluded_ids with empty status columns for manual input
    excluded_data = {
        'Sample id': excluded_ids,
        'Excluded': ['Yes']*len(excluded_ids)
    }
    # Create a column for each unique name_alignment and add it to the data dictionary
    for name_alignment in unique_name_alignments:
        # Dynamically create a new column name for each unique 'name_alignment'
        column_name = f'Status {name_alignment} alignment'
        excluded_data[column_name] = [None] * len(excluded_ids)  # Initialize with None for manual input
    # Add the 'Comments' column after all the Status columns
    excluded_data['Comments'] = [None] * len(excluded_ids)

    # Concatenate the common_ids DataFrame with the excluded_data DataFrame at the bottom
    df = pd.concat([df, pd.DataFrame(excluded_data)], ignore_index=True)

    df['Sample id'] = '"' + df['Sample id'].astype(str) + '"'

    # Save to Excel (instead of CSV) to enable styling
    excel_file_path = output_path + "Alignments_validated.xlsx"
    df.to_excel(excel_file_path, index=False, engine='openpyxl')

    # Load the Excel workbook and sheet
    wb = load_workbook(excel_file_path)
    ws = wb.active

    # Define fill colors for gray, light green, and light red
    gray_fill = PatternFill(start_color="D3D3D3", end_color="D3D3D3", fill_type="solid")
    light_green_fill = PatternFill(start_color="A8E6A1", end_color="A8E6A1", fill_type="solid")  # Light green
    light_red_fill = PatternFill(start_color="FFCCCB", end_color="FFCCCB", fill_type="solid")  # Light red

    # Apply gray fill to all cells in the row if 'Excluded' is 'Yes'
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, min_col=1, max_col=5):  # Columns 1-5 are all the columns (Sample id, Excluded, Status, Comments)
        if row[1].value == 'Yes':  # Check if 'Excluded' (2nd column) is 'Yes'
            for cell in row:  # Apply the gray fill to all cells in the row
                cell.fill = gray_fill

    # Create conditional formatting for 'Status first alignment' and 'Status second alignment' columns
    # Light green for 'Y'
    green_rule = FormulaRule(formula=['$C2="Y"'], fill=light_green_fill)  # Column C = 'Status first alignment'
    ws.conditional_formatting.add('C2:C{}'.format(ws.max_row), green_rule)

    green_rule_2 = FormulaRule(formula=['$D2="Y"'], fill=light_green_fill)  # Column D = 'Status second alignment'
    ws.conditional_formatting.add('D2:D{}'.format(ws.max_row), green_rule_2)

    # Light red for 'N'
    red_rule = FormulaRule(formula=['$C2="N"'], fill=light_red_fill)  # Column C = 'Status first alignment'
    ws.conditional_formatting.add('C2:C{}'.format(ws.max_row), red_rule)

    red_rule_2 = FormulaRule(formula=['$D2="N"'], fill=light_red_fill)  # Column D = 'Status second alignment'
    ws.conditional_formatting.add('D2:D{}'.format(ws.max_row), red_rule_2)

    # Adjust the column widths for better readability of the headers
    ws.column_dimensions['C'].width = 25  # Status first alignment
    ws.column_dimensions['D'].width = 25  # Status second alignment
    ws.column_dimensions['E'].width = 60  # Status second alignment

    # Save the modified workbook with conditional formatting and adjusted column widths
    wb.save(excel_file_path)

    print(f'Data saved to {excel_file_path}')









        
def get_parameters(params_file_path, rotation_params, excluded_ids):
    ''' Get the parameters define in the visualization step. '''
    # Load parameters
    with open(params_file_path, "r") as infile:
        params_dict = json.load(infile)
    folder_paths = params_dict.get("folder_paths")
    common_ids = params_dict.get("common_ids")
    panels = params_dict.get("panels")
    output_path = params_dict.get("output_path")

    ## Exclude specific ids
    common_ids = [id for id in common_ids if id not in excluded_ids]

    ## Define which image should be turn 90 degrees clock wise (1), 180 degrees (2) or 90 degrees counter clock wise (3)
    turn_img_id_dict = {f"{id_}_{panel}": 0 for panel in panels for id_ in common_ids}
    for rota in rotation_params:
        # Split the string to get the key and value
        key, value = rota.rsplit('_', 1)
        # Assign the new value to correponding key
        turn_img_id_dict[key] = int(value)  # Convert string to integer
    
    return folder_paths, common_ids, panels, output_path, turn_img_id_dict

def panels_name_alignment(panels, folder_paths):
    ''' Get the name of the panels and corresponding folder for each alignment. '''
    # If there is 2 panels in total then the reference panel is at the end
    if len(panels) == 2:
        panel_alignment_dict = {"_".join(panels): folder_paths}
    # If there is 3 panels in total then the reference panel is in the middle
    elif len(panels) == 3:
        panel_alignment_dict = {
            "_".join([panels[0], panels[1]]): [folder_paths[0], folder_paths[1]],
            "_".join([panels[2], panels[1]]): [folder_paths[2], folder_paths[1]],
        }
    
    return panel_alignment_dict

def load_downscaled_imgs(folder_path1, folder_path2, id, folder_name_alignment, turn_img_id_dict):
    # Define path of the qptiff
    panels = folder_name_alignment.split('_')
    path1 = folder_path1 + id + f"_panel_{panels[0]}.unmixed.qptiff"
    path2 = folder_path2 + id + f"_panel_{panels[1]}.unmixed.qptiff"
    # print("path1:", path1)
    # print("path2:", path2)
    
    ## Load images
    diff_index_pages = 10
    tif_tags1, channel_list1, channel_name_dictionary1, img1_resize, scale_percent1 = get_data_alignment(path1, diff_index_pages)
    tif_tags2, channel_list2, channel_name_dictionary2, img2_resize, scale_percent2 = get_data_alignment(path2, diff_index_pages)
    # print("Compression")
    comp_lvl1 = int((1/scale_percent1))
    comp_lvl2 = int((1/scale_percent2))
    # print(f"img1: 1/{comp_lvl1}ème")
    # print(f"img2: 1/{comp_lvl2}ème")
    
    # Check if both images are the same size
    if comp_lvl2 < comp_lvl1:
        # print("Not same compression level")
        tif = TiffFile(path1)
        while comp_lvl2 < comp_lvl1:
            diff_index_pages += 8
            most_comp_DAPI_index = len(tif.pages) - diff_index_pages
            img1_resize = tif.pages[most_comp_DAPI_index].asarray()
            scale_percent1 = img1_resize.shape[0] / tif_tags1['ImageLength']
            comp_lvl1 = int((1/scale_percent1))
            # print("Compression")
            # print(f"img1: 1/{comp_lvl1}ème")
            # print(f"img2: 1/{comp_lvl2}ème")
    
    if comp_lvl1 < comp_lvl2:
        # print("Not same compression level")
        tif = TiffFile(path2)
        while comp_lvl1 < comp_lvl2:
            diff_index_pages += 8
            most_comp_DAPI_index = len(tif.pages) - diff_index_pages
            img2_resize = tif.pages[most_comp_DAPI_index].asarray()
            scale_percent2 = img2_resize.shape[0] / tif_tags2['ImageLength']
            comp_lvl2 = int((1/scale_percent2))
            # print("Compression")
            # print(f"img1: 1/{comp_lvl1}ème")
            # print(f"img2: 1/{comp_lvl2}ème")

    # Scale images to 8bit
    img1_8bit = cv2.convertScaleAbs(img1_resize)
    img2_8bit = cv2.convertScaleAbs(img2_resize)

    # Rotate or not the first image
    img1_8bit = transform_image(img1_8bit, turn_img_id_dict[id + "_" + panels[0]])
    img2_8bit = transform_image(img2_8bit, turn_img_id_dict[id + "_" + panels[1]])

    # Convert to sitk image format
    img1 = sitk.GetImageFromArray(img1_8bit)
    img2 = sitk.GetImageFromArray(img2_8bit)
    img1 = sitk.Cast(img1, sitk.sitkFloat32)
    img2 = sitk.Cast(img2, sitk.sitkFloat32)
    
    # Save images
    img1_arr = sitk.GetArrayFromImage(img1)
    img2_arr = sitk.GetArrayFromImage(img2)

    return tif_tags1, channel_list1, channel_name_dictionary1, scale_percent1, tif_tags2, channel_list2, channel_name_dictionary2, scale_percent2, img1, img2, img1_arr, img2_arr


def get_data_alignment(path_qptiff, diff_index_pages):
    tif = TiffFile(path_qptiff)
    # Get tags from DAPI channel in .pages
    tif_tags= {}
    for tag in tif.pages[0].tags.values():
        name, value = tag.name, tag.value
        tif_tags[name] = value
    channel_list, channel_name_dictionary = getLabels(tif_tags)
    # Load the most compressed DAPI image in .pages
    most_comp_DAPI_index = len(tif.pages) - diff_index_pages
    img_resize = tif.pages[most_comp_DAPI_index].asarray()
    scale_percent = img_resize.shape[0] / tif_tags['ImageLength']
    
    return tif_tags, channel_list, channel_name_dictionary, img_resize, scale_percent
    

def transform_image(img, operation):
    """Applies the transformation to the image based on the operation code."""
    if operation == 1:
        # Rotate counterclockwise by 90 degrees (equiv. to clockwise -90 degrees)
        return np.rot90(img, -1)
    elif operation == 2:
        # Flip left-right and then up-down
        return np.fliplr(np.flipud(img))
    elif operation == 3:
        # Rotate clockwise by 90 degrees
        return np.rot90(img)
    else:
        # No operation or undefined operation, return image as is
        return img
    

def getLabels(tif_tags):
    substr = "<ScanColorTable-k>"
    start = 0
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

    return channel_list, dictionary