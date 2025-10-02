import os
import argparse

def parse_arguments():
    """Define arguments."""
    parser = argparse.ArgumentParser(description="IntegrAlign: A comprehensive tool for multi-immunofluorescence panel integration through image alignement")
    subparsers = parser.add_subparsers(dest="task", required=True, help="Choose a task to perform")

    # Preprocessing Task
    preprocess_parser = subparsers.add_parser("visualize", help="Generate visualization report of the scans")
    preprocess_parser.add_argument('--scans', nargs='+', required=True,
                                    help="Paths to 2 or 3 panel folders of the scans (.qptiff). Should be ordered with the reference panel in the middle for 3 panels or on the right for 2 panels.")
    preprocess_parser.add_argument('--annotations', nargs='+', default=None,
                                    help="Paths to 2 or 3 panel folders of the annotation files (.annotation or .geojson). Keep the same order used in the visualize step.")
    preprocess_parser.add_argument('--panels', nargs='+', required=True,
                                    help="List of the panel names. Keep the same order as folders.")
    preprocess_parser.add_argument('--output', required=True, help="Path to save the preprocessing report")

    # saveImages Task
    saveimgs_parser = subparsers.add_parser("saveimgs", help="Save images and metadata for alignment step")
    saveimgs_parser.add_argument('--params', required=True,
                                   help="Paths to the params file")
    saveimgs_parser.add_argument('--exclude', type=str, nargs='+', default=[],
                                   help="List of patient IDs to exclude from alignment")
    saveimgs_parser.add_argument('--brightness', type=float, default=1, 
                                   help="Brightness factor for the cropping and manual alignment visu (optional). Value by default : 1")
    saveimgs_parser.add_argument('--HALOrotation', type=str, default=None, 
                                   help='Path to an Excel (.xlsx) file containing rotation information for each slide (optional). Use this when slides were rotated in HALO, causing output coordinates to be rotated but unchanged image file. The file must include ID_patient and Panel columns.')
                                   
    # Alignment Task
    alignment_parser = subparsers.add_parser("align", help="Align panels, generate QC report and a table with combined cell types")
    alignment_parser.add_argument('--dwnscimg', required=True,
                                   help="Paths to the downscaled images file (.pkl)")
    alignment_parser.add_argument('--tables', nargs='+', required=True,
                                   help="Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.")
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
    

    # finetuning Task
    finetuning_parser = subparsers.add_parser("finetuning", help="Validate the data for a specific patient")
    finetuning_parser.add_argument('--id', required=True, help="Patient ID")
    finetuning_parser.add_argument('--meshsize', type=int, nargs='+', required=True, help="Mesh size of each alignment. \nExemple : 6 8  will align the first panel to the middle panel using mesh size of 6 and align the last panel to the middle panel using mesh size of 8")
    finetuning_parser.add_argument('--dwnscimg', required=True,
                                    help="Paths to the downscaled images file (.pkl)")
    finetuning_parser.add_argument('--tables', nargs='+', required=True,
                                    help="Paths to 2 or 3 panel folders of the coordinate tables (.csv). Keep the same order used in the visualize step.")
    finetuning_parser.add_argument('--visualization', type=str, choices=["0", "1", "2", "all"], default="0",
                                    help=("Enable mirrored cursor visualization for manual quality control (QC) of alignments. "
                                        "Options: '0' (no visualization, default), '1' (visualize the first alignment), "
                                        "'2' (visualize the second alignment), 'all' (visualize all alignments)."))
    finetuning_parser.add_argument('--scans', nargs='+',
                                    help="Paths to 2 or 3 panel folders of the scans (.qptiff). THE FILES SHOULD BE THE SAME AS THE ONES USED IN THE BATCH ALIGNMENT BUT CAN BE IN ANOTHER PLACE, if you want to do an alignment with new slides relaunch the full pipeline (with batch alignement) so it uses the new downscaled slides (need to verify the rotation of the slide). Should be ordered with the reference panel in the middle for 3 panels or on the right for 2 panels.")
    finetuning_parser.add_argument('--resolution', type=float, default=2.012948251135, 
                                    help="Resolution in µm (optional). Value by default : 2.012948251135")
    finetuning_parser.add_argument('--metric', type=str, default="Correlation", 
                                    help='Similarity Metric (optional). Reflects the relationship between the intensities of the images, see https://simpleitk.readthedocs.io/en/master/registrationOverview.html . Value by default : "Correlation"')
    finetuning_parser.add_argument('--raster', type=int, default=30, 
                                     help="Pixel size of the raster in micron (optional). Value by default : 30µm")
    finetuning_parser.add_argument('--alpha', type=int, default=0.4, 
                                    help="Transparency of the reference panel in red for visualization of the alignment (optional). Value by default : 0.4")

    return parser.parse_args()


def validate_nb_panels(value):
    if len(value) not in [2, 3]:
        raise argparse.ArgumentTypeError("There should be exactly 2 or 3 panels (folders and panels).")
    return value
