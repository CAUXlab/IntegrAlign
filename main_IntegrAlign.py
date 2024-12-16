from src.utils import parse_arguments
from src.visualization import visualization
from src.save_images import save_downscaled_images
from src.alignment import alignment
from src.validation import validation

### EXEMPLES ###
'''
python main_IntegrAlign.py visualize --folders "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/" --panels T DC TLS --output "/Users/leohermet/Documents/IntegrAlign/main/output/"
python main_IntegrAlign.py saveimgs --params "/Users/leohermet/Documents/IntegrAlign/main/output/params.json" --exclude 02006 06001 08006 --rotate 01008_DC_2
python main_IntegrAlign.py align --dwnscimg "/Users/leohermet/Documents/IntegrAlign/main/output/downscaled_images.pkl" --tables "/Users/leohermet/Documents/SPIAT/BI3/T/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/DC/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/TLS/Output/cell_positions_data_clean_CellType/" --annotations "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/Annotations_geojson_export_2024-08-20/1724135774_Annotations/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/annotations_finales/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/Geojson_Annotations/Annotations_geojson/" --maxMS 3

'''

if __name__ == "__main__":
    args = parse_arguments()

    # Process the task
    if args.task == "visualize":
        # Validate folder structure
        # validate_folder(args.folders)
        visualization(args.folders, args.panels, args.output)

    elif args.task == "saveimgs":
        save_downscaled_images(args.params, args.exclude, args.rotate)

    elif args.task == "align":
        alignment(args.dwnscimg, args.tables, args.annotations, args.resolution, args.maxMS, args.metric, args.raster, args.alpha)

    elif args.task == "validate":
        validation(args.patient, args.meshsize, args.visualization, args.output)