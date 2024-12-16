import os
import argparse

def parse_arguments():
    """Define arguments."""
    parser = argparse.ArgumentParser(description="IntegrAlign: A comprehensive tool for multi-immunofluorescence panel integration through image alignement")
    subparsers = parser.add_subparsers(dest="task", required=True, help="Choose a task to perform")

    # Preprocessing Task
    preprocess_parser = subparsers.add_parser("visualize", help="Generate visualization report of the scans")
    preprocess_parser.add_argument('--folders', nargs='+', required=True,
                                    help="Paths to 2 or 3 panel folders of the scans (.qptiff). Should be ordered with the reference panel in the middle for 3 panels or on the right for 2 panels.")
    preprocess_parser.add_argument('--panels', nargs='+', required=True,
                                    help="List of the panel names. Keep the same order as folders.")
    preprocess_parser.add_argument('--output', required=True, help="Path to save the preprocessing report")

    # saveImages Task
    alignment_parser = subparsers.add_parser("saveimgs", help="Save images and metadata for alignment step")
    alignment_parser.add_argument('--params', required=True,
                                   help="Paths to the params file")
    alignment_parser.add_argument('--exclude', type=str, nargs='+', default=[],
                                   help="List of patient IDs to exclude from alignment")
    alignment_parser.add_argument('--rotate', type=str, nargs='+', default=[],
                                   help="String of the image to rotate with the rotation parameter. \nExemple : \"01008_DC_2\" will turn 180 degrees the scan of the panel DC for the patient with id 01008. \nRotation parameters: turn 90 degrees clock wise (1), 180 degrees (2) or 90 degrees counter clock wise (3)")
    
    # Alignment Task
    alignment_parser = subparsers.add_parser("align", help="Align panels, generate QC report and a table with combined cell types")
    alignment_parser.add_argument('--dwnscimg', required=True,
                                   help="Paths to the downscaled images file (.pkl)")
    alignment_parser.add_argument('--tables', nargs='+', required=True,
                                   help="Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.")
    alignment_parser.add_argument('--annotations', nargs='+', required=True,
                                   help="Paths to 2 or 3 panel folders of the annotation files (.geojson). Keep the same order used in the visualize step.")
    alignment_parser.add_argument('--resolution', type=float, default=2.012948251135, 
                                   help="Resolution in µm (optional). Value by default : 2.012948251135")
    alignment_parser.add_argument('--maxMS', type=int, default=10, 
                                   help="Max mesh size number (optional). Value by default : 10")
    alignment_parser.add_argument('--metric', type=str, default="Correlation", 
                                   help='Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"')
    alignment_parser.add_argument('--raster', type=int, default=30, 
                                   help="Pixel size of the raster in micron (optional). Value by default : 30µm")
    alignment_parser.add_argument('--alpha', type=int, default=0.4, 
                                   help="Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4")

    # Validation Task
    validation_parser = subparsers.add_parser("validation", help="Validate the data for a specific patient")
    validation_parser.add_argument('--patient', required=True, help="Patient ID to validation")
    validation_parser.add_argument('--meshsize', type=int, required=True, help="Mesh size")
    validation_parser.add_argument('--visualization', action='store_true',
                                     help="Enable mirrored cursor visualization for manual QC")
    validation_parser.add_argument('--output', required=True, help="Path to save combined coordinates")

    return parser.parse_args()

'''
def validate_folder(folder_paths):
    """Validate folder structure."""

    for folder in folder_paths:
        scans_path = os.path.join(folder, "SCANS")
        annotations_path = os.path.join(folder, "Annotations")
        if not os.path.isdir(scans_path) or not os.path.isdir(annotations_path):
            raise FileNotFoundError(f"Required subfolders missing in {folder}: SCANS/, Annotations/")
'''


def validate_nb_panels(value):
    if len(value) not in [2, 3]:
        raise argparse.ArgumentTypeError("There should be exactly 2 or 3 panels (folders and panels).")
    return value
