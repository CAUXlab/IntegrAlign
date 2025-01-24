# IntegrAlign

A comprehensive tool for multi-immunofluorescence panel integration through image alignement

## Prerequisites

Ensure you have Python 3.9.13 installed on your system. You can download it from the [official Python website](https://www.python.org/downloads/).

## Installation

We recommend running IntegrAlign in a virtual environment to prevent conflicts with other packages or system dependencies. Here's how you can do it:

1. Install `virtualenv` package, if it's not already installed:

```bash
pip install virtualenv
```

2. Navigate to the directory where you want to create the virtual environment and run:

```bash
virtualenv IntegrAlign_env
```

3. To activate the virtual environment, use:
   
On Windows:

```bash
IntegrAlign_env\Scripts\activate
```

On Unix or MacOS:

```bash
source IntegrAlign_env/bin/activate
```

4. Navigate to the IntegrAlign directory (with setup.py). Then install the dependencies by running:

```bash
pip install .
```

## Quick start

### 1 Visualization

Once you have activated your virtual environment and installed the necessary dependencies, you can run the first step of the tool by using:

```bash
python main_IntegrAlign.py visualize --scans "panel_T/SCANS/" "panel_DC/SCANS/" "panel_TLS/SCANS/" --annotations "T/Annotations/" "DC/Annotations/" "TLS/Annotations/" --panels T DC TLS --output "output/"
```

This steps will create a report file ```images_report.html``` that plots every panels for each patients. This will help to determine which slides cannot be aligned and if some panels need to be rotated for best alignment.

<ins>Options:</ins>

--scans - 2 or 3 paths to the folders of each panel that contains the .qptiff files of every patients.

--annotations - Paths to 2 or 3 panel folders of the annotation files (.geojson). Optional.

--panels - Panel names.

--output - Path to the main output folder.

Notes : **The order is crucial** for --scans, --annotations and --panels options. It should be consistent at every step and the reference panel has to be in the middle for alignment of 3 panels or in last position for alignment of 2 panels.

### 2 Save downscaled images

This steps will preprocess the data (downscaling and rotation) and save the images without excluded patients. This will generate a file containing downscaled images, allowing for alignment processing without needing access to the ```.qptiff``` files.

```bash
python main_IntegrAlign.py saveimgs --params "output/params.json" --exclude 02006 06001 08006 --rotate 01008_DC_2
```

<ins>Options:</ins>

--params - Path to the params file generated in the Visualization step (located in the defined output folder).

--exclude - Patient IDs to exclude from alignment.

--rotate - Image to rotate with the rotation parameter. Rotation parameters: turn 90 degrees clock wise (1), 180 degrees (2) or 90 degrees counter clock wise (3).

Notes : You provide strings to the --rotate option. For example, "01008_DC_2" will rotate the scan of the DC panel for the patient with ID 01008 by 180 degrees.

### 3 Alignment

This is the main step of this pipeline, here we align all downscaled images, generate an html report for each alignment and merged tables with cells' coordinates of each panel.

Indeed this step run different alignment for different mesh size and then the optimal mesh size is chosen based on the correlation result between the rasters of the cells coordinates (all cells from DAPI staining). This way we avoid the use of a mesh size that would induce deformations within the tissue.

This can be run in a cluster computing environment using the ```params.json```, ```downscaled_images.pkl``` and the coordinate tables files.

```bash
python main_IntegrAlign.py align --dwnscimg "output/downscaled_images.pkl" --tables "T/Cell_positions/" "DC/Cell_positions/" "TLS/Cell_positions/" --annotations "T/Annotations/" "DC/Annotations/" "TLS/Annotations/"
```

<ins>Options:</ins>

--dwnscimg - Paths to the downscaled images file (.pkl).

--tables - Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.

--annotations - Paths to 2 or 3 panel folders of the annotation files (.geojson). Keep the same order used in the visualize step.

--resolution - Resolution in µm (optional). Value by default : 2.012948251135

--maxMS - Max mesh size number (optional). Value by default : 10

--metric - Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"

--raster - Pixel size of the raster in micron (optional). Value by default : 30µm

--alpha - Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4

Notes : If there is no intersections found between the analysis areas of each panel, one of the alignment is considered suboptimal, then it will continue to the next patient without saving merged annotations, tables or plots.
Also, the given parameters have to be in the same panel order as in step 1 Visualization.

### 4 Fine tuning

This step is used to refine the alignment(s) for a specific patient. For example, if you think the chosen mesh size is not optimal, if you want to change any parameter for better visualization (alpha, raster size) or if you need a more precise look over the alignment using mirrored cursor visualization.

This will write, with the specified **mesh size** (resolution of the deformation), a "_refined" alignment report(s) and rewrite the merged coordinates tables and the merged annotations files. 

You can also, if needed, validate the alignment with the mirrored cursor visualization and manually clip areas where the alignment is inadequate. This step of validation **requires the paths to the .qptiff files** of the scans.

```bash
python main_IntegrAlign.py finetuning --id 02005 --meshsize 6 8 --dwnscimg "output/downscaled_images.pkl" --tables "T/Cell_positions/" "DC/Cell_positions/" "TLS/Cell_positions/" --annotations "T/Annotations/" "DC/Annotations/" "TLS/Annotations/" --visualization all --scans "panel_T/SCANS/" "panel_DC/SCANS/" "panel_TLS/SCANS/"
```

<ins>Options:</ins>

--id - Patient ID.

--meshsize - Mesh size of each alignment. \nExemple : 6 8  will align the first panel to the middle panel using a mesh size of 6 and align the last panel to the middle panel using a mesh size of 8

--dwnscimg - Paths to the downscaled images file (.pkl).

--tables - Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.

--annotations - Paths to 2 or 3 panel folders of the annotation files (.geojson). Keep the same order used in the visualize step.

--visualization - Enable mirrored cursor visualization for manual quality control (QC) of alignments, folder paths to the scans is required. Options: '0' (no visualization, default), '1' (visualize the first alignment), '2' (visualize the second alignment), 'all' (visualize all alignments).

--scans - Paths to 2 or 3 panel folders of the scans (.qptiff). Should be ordered with the reference panel in the middle for 3 panels or on the right for 2 panels. THE FILES SHOULD BE THE SAME AS THE ONES USED IN THE BATCH ALIGNMENT BUT CAN BE IN ANOTHER PLACE, if you want to do an alignment with new slides relaunch the full pipeline (with batch alignement) so it uses the new downscaled slides (verify the rotation of the slide). 

--resolution - Resolution in µm (optional). Value by default : 2.012948251135

--metric - Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"

--raster - Pixel size of the raster in micron (optional). Value by default : 30µm

--alpha - Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4

Notes : The given parameters have to be in the same panel order as in step 1 Visualization.

## What is a good alignment?

You have all the information needed to check if the alignment performed well in the corresponding html report located in Alignment/report/ from your output folder. The chosen optimal mesh size is shown in red.

Here there is detailled visualization for each alignment (with different mesh size). The most important images will be the red and blue superposition of both images and the jacobian determinant. 

The first one will show you if the images align well, for example here the alignment between T and DC seems well aligned. 
![image](https://github.com/user-attachments/assets/8bd0e81c-f8d0-4f2e-9c98-caa6986b5526)

The jacobian will show you if deformation is done inside the tissue, here there is no deformation inside the tissue. 
![image](https://github.com/user-attachments/assets/03d1e12f-ca0b-4858-b6f2-ed0203b7278e)

For the alignment between TLS and DC otherwise we can see the blue image doesn't fit well with the red image inside the tissue.
![image](https://github.com/user-attachments/assets/61fa4a95-855e-4114-bd92-783878f355aa)

We also have huge deformation inside the tissue for most of mesh sizes (here mesh size at 8, 9 and 10 ), this alignment did not performed well and need to be excluded in further analysis.
![image](https://github.com/user-attachments/assets/4e47c660-97c0-40a5-b128-5b96f9ab59a5)

Indeed when we visualize the images in higher resolution with the mirrored cursor we can observe that the TLS tissue exhibits irregular tearing throughout, which is inconsistent with the tearing seen in the DC panel. This make the TLS panel impossible to align with the DC panel in our case since the tissue structure is completly different. 

This explains the importance of serial slides being cut in a consistent manner in order to have the same tissue structure.

## Exiting the virtual environment
When you are done, you can deactivate the virtual environment by simply running:

```bash
deactivate
```

## Dependencies
IntegrAlign relies on several Python libraries, listed in the setup.py file:


- ttkthemes==2.1
- customtkinter
- napari[all]
- tk
- opencv-python
- SimpleITK
- scikit-image
- matplotlib
- imagecodecs
- rasterio
- shapely
- geopandas



