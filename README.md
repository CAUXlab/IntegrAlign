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

### 0 Inputs

IntegrAlign requires three main types of files to align multi-IF slides:

1. Scan Images

- Recommended format: .qptiff
Use this format when
.phenotyping has been performed in inForm,
.or when phenotyping has been performed in HALO and unmixing was done directly using the Vectra Polaris scanner.

- Alternative format: .tif (HALO)
Use this when phenotyping was performed in HALO and unmixing was handled by inForm.

2. Annotation Files

- Supported formats: .annotations or .geojson

- Important: To avoid mismatches, ensure annotation names are consistent across patients or use the --names argument to provide all possible names for each categorie.

3. Coordinate Tables

- Format: .csv

- Requirements:

For HALO outputs: must include XMin, XMax, YMin, and YMax columns (pixel coordinates).

For inForm outputs: must include x and y columns (coordinates in microns).

- Important: Avoid raw HALO / Inform outputs with Intensity and Classication columns as **you cannot infer functional markers between aligned panel since cells aren't exactly the same between serial slides**. Instead chose the corresponding cell types to avoid double positive cells and preprocess the tables into final tables with :
HALO : XMin, XMax, YMin, YMax, cell_id, cell_type (from lineage markers), and phenotype (all positive marker in order to have the functional information).
Inform : x (in µm), y (in µm) , cell_id, cell_type, and phenotype.

### 1 Visualization

Once you have activated your virtual environment and installed the necessary dependencies, you can run the first step of the tool by using:

```bash
python main_IntegrAlign.py visualize --scans "panel_T/SCANS/" "panel_DC/SCANS/" "panel_TLS/SCANS/" --annotations "T/Annotations/" "DC/Annotations/" "TLS/Annotations/" --namesEmpty Empty No_tissue --namesArtefacts Artefacts Manual_Artefacts --namesAnalysisArea Analysis_Area --panels T DC TLS --output "output/"
```

This steps will create a report file ```images_report.html``` that plots every panels for each patients. This will help identify which slides cannot be aligned, saving time when processing large patient cohorts.

<ins>Options:</ins>

--scans - 2 or 3 paths to the folders of each panel that contains the .qptiff files of every patients.

--annotations - Paths to 2 or 3 panel folders of the annotation files (.geojson). Optional.

--namesEmpty - Name(s) of the empty areas in the annotations files.

--namesArtefacts - Name(s) of the Artefacts areas in the annotations files. 

--namesAnalysisArea - Name(s) of the Analysis areas in the annotations files.

--panels - Panel names.

--output - Path to the main output folder.

Notes : **The order is crucial** for --scans, --annotations and --panels options. It should be consistent at every step and the reference panel has to be in the middle for alignment of 3 panels or in last position for alignment of 2 panels.

### 2 Save downscaled images

This step will preprocess the data by downscaling the paired images for alignment at the same compression level. Also you can appply cropping (slides with amygdales for control) and manual pre alignment (improve final results) and save the resulting images for the patients that have not been excluded. This will generate a file containing the downscaled images, allowing for alignment processing without needing access to the ```.qptiff``` files. An excel file (```Alignments_validated.xlsx```) will also be generated to track the validation status of each registration after the alignment step completed.

```bash
python main_IntegrAlign.py saveimgs --params "output/params.json" --exclude 02006 06001 08006
```

<ins>Options:</ins>

--params - Path to the params file generated in the Visualization step (located in the defined output folder).

--exclude - Patient IDs to exclude from alignment.

--brightness - Brightness factor for the cropping and manual alignment visualization (optional). Value by default : 1


### 3 Alignment

This is the main step of this pipeline, here we align all downscaled images, generate an html report for each alignment and merged tables with cells' coordinates of each panel.

This step run different alignment for different mesh size and then the optimal mesh size is chosen based on the correlation result between the rasters of the cells coordinates (all cells from DAPI staining). This way we avoid the use of a mesh size that would induce deformations within the tissue.

This can be run in a cluster computing environment using the ```downscaled_images.pkl``` and the coordinates tables files.

```bash
python main_IntegrAlign.py align --dwnscimg "output/downscaled_images.pkl" --tables "T/Cell_positions/" "DC/Cell_positions/" "TLS/Cell_positions/"
```

<ins>Options:</ins>

--dwnscimg - Paths to the downscaled images file (.pkl).

--tables - Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.

--resolution - Resolution in µm (optional). Value by default : 2.012948251135

--maxMS - Max mesh size number (optional). Value by default : 10

--metric - Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"

--raster - Pixel size of the raster in micron (optional). Value by default : 30µm

--alpha - Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4

Notes : The given parameters have to be in the same panel order as in step 1 Visualization.

### 4 Fine tuning

This step allows for detailed visualization and fine-tuning of the alignment(s) for a specific patient. For example, if you think the chosen mesh size is not optimal, if you want to change any parameter for better visualization (alpha, raster size) or if you need a more precise look over the alignment using mirrored cursor visualization.

This will write, with the specified **mesh size** (resolution of the deformation), a "_refined" alignment report(s) and rewrite the merged coordinates tables and the merged annotations files, if specified. 

You can also validate the alignment with the mirrored cursor visualization and manually clip areas where the alignment is inadequate. This step of validation **requires the paths to the .qptiff files** of the scans.

```bash
python main_IntegrAlign.py finetuning --id 02005 --meshsize 6 8 --dwnscimg "output/downscaled_images.pkl" --tables "T/Cell_positions/" "DC/Cell_positions/" "TLS/Cell_positions/" --visualization all --scans "panel_T/SCANS/" "panel_DC/SCANS/" "panel_TLS/SCANS/"
```

<ins>Options:</ins>

--id - Patient ID.

--meshsize - Mesh size of each alignment. \nExemple : 6 8  will align the first panel to the middle panel using a mesh size of 6 and align the last panel to the middle panel using a mesh size of 8

--dwnscimg - Paths to the downscaled images file (.pkl).

--tables - Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.

--visualization - Enable mirrored cursor visualization for manual quality control (QC) of alignments, folder paths to the scans is required. Options: '0' (no visualization, default), '1' (visualize the first alignment), '2' (visualize the second alignment), 'all' (visualize all alignments).

--scans - Paths to 2 or 3 panel folders of the scans (.qptiff). Should be ordered with the reference panel in the middle for 3 panels or on the right for 2 panels. THE FILES SHOULD BE THE SAME AS THE ONES USED IN THE BATCH ALIGNMENT BUT CAN BE IN ANOTHER PLACE, if you want to do an alignment with new slides relaunch the full pipeline (with batch alignement) so it uses the new downscaled slides (verify the rotation of the slide). 

--resolution - Resolution in µm (optional). Value by default : 2.012948251135

--metric - Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"

--raster - Pixel size of the raster in micron (optional). Value by default : 30µm

--alpha - Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4

Notes : The given parameters have to be in the same panel order as in step 1 Visualization.

### 5 Outputs

You can keep track of the alignments using the Alignments_validated.xlsx file located in your output directory. This file is generated during the second step (Save downscaled images).

In the Alignment/merged_tables/ directory, you will find the merged coordinate tables for each patient. These tables contain the concatenated data from all panels, with coordinates transformed accordingly to the reference panel (except for the reference slide, whose coordinates remain unchanged) and in µm. A "Panel" column is also added to indicate which cells belong to which panel. If you need unique cell IDs across panels, you can concatenate the cell ID column with the panel name (using an underscore for example).

If the column names in the panel tables match, they will be merged. Any columns unique to specific panels will appear as NaN or empty cells in the resulting CSV for panels that don't include them.

## What is a good alignment?

You have all the information needed to check if the alignment performed well in the corresponding html report located in Alignment/report/ from your output folder. The chosen optimal mesh size is shown in red.

Here there is detailled visualization for each alignment (with different mesh size). The most important images will be the red and blue superposition of both images and the jacobian determinant. 

The first one will show you if the images align well, for example here the alignment between T and DC seems well aligned. 
![overlay_good_alignment](images/overlay_good_alignment.png)

The jacobian will show you if deformation is done inside the tissue, here there is no deformation inside the tissue. 
![jacobian_good_alignment](images/jacobian_good_alignment.png)

For the alignment between TLS and DC otherwise we can see the blue image doesn't fit well with the red image inside the tissue.
![overlay_bad_alignment](images/overlay_bad_alignment.png)

We also have huge deformation inside the tissue for most of mesh sizes (here mesh size at 8, 9 and 10 ), this alignment did not performed well and need to be excluded in further analysis.
![jacobian_bad_alignment](images/jacobian_bad_alignment.png)

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



