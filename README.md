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
python main_IntegrAlign.py visualize --folders "panel_T/SCANS/" "panel_DC/SCANS/" "panel_TLS/SCANS/" --panels T DC TLS --output "output_path/"
```

This steps will create a report file ```images_report.html``` that plots every panels for each patients. This will help to determine which slides cannot be aligned and if some panels need to be rotated for best alignment.

<ins>Options:</ins>

--folders - 2 or 3 paths to the folders of each panel that contains the .qptiff files of every patients.

--panels - Panel names.

--output - Path to the main output folder.

Notes : **The order is crucial** for --folders and --panels options. It should be consistent at every step and the reference panel has to be in the middle for alignment of 3 panels or in last position for alignment of 2 panels.

### 2 Save downscaled images

This steps will preprocess the data (downscaling and rotation) and save the images without excluded patients. This will generate a file containing downscaled images, allowing for alignment processing without needing access to the ```.qptiff``` files.

```bash
python main_IntegrAlign.py saveimgs --params "main/output/params.json" --exclude 02006 06001 08006 --rotate 01008_DC_2
```

<ins>Options:</ins>

--params - Path to the params file generated in the Visualization step (located in the defined output folder).

--exclude - Patient IDs to exclude from alignment.

--rotate - Image to rotate with the rotation parameter. Rotation parameters: turn 90 degrees clock wise (1), 180 degrees (2) or 90 degrees counter clock wise (3).

Notes : You can provide strings to the --rotate option. For example, "01008_DC_2" will rotate the scan of the DC panel for the patient with ID 01008 by 180 degrees.

### 3 Alignment

This is the main step of this pipeline, here we align all downscaled images, generate an html report for each alignment and merged tables with cells' coordinates of each panel.
This can be run in a cluster computing environment using the ```params.json```, ```downscaled_images.pkl``` and the coordinates tables.

```bash
python main_IntegrAlign.py align --dwnscimg "main/output/downscaled_images.pkl" --tables "T/Cell_positions/" "DC/Cell_positions/" "TLS/Cell_positions/" --annotations "T/Annotations/" "DC/Annotations/" "TLS/Annotations/" --maxMS 3
```

<ins>Options:</ins>

--

Notes : The given parameters have to be in the same panel order as in step 1.

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



