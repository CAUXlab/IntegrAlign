from src.utils.args import parse_arguments
from src.visualization import visualization
from src.save_images import save_downscaled_images
from src.alignment import alignment
from src.finetuning import finetuning

### EXEMPLES ###
'''
python main_IntegrAlign.py visualize --scans "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/SCANS/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/SCANS/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/SCANS/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/Annotations" --panels T DC TLS --output "/Users/leohermet/Downloads/IntegrAlign-main/output/"
python main_IntegrAlign.py saveimgs --params "/Users/leohermet/Downloads/IntegrAlign-main/output/params.json" --rotate 01008_DC_2
python main_IntegrAlign.py align --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/Cell_coordinates/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/Cell_coordinates/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/Cell_coordinates/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/Annotations" --maxMS 3
python main_IntegrAlign.py finetuning --id 02005 --meshsize 3 3 --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/Cell_coordinates/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/Cell_coordinates/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/Cell_coordinates/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/Annotations" --visualization all --scans "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_T/SCANS/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_DC/SCANS/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data/panel_TLS/SCANS/"



python main_IntegrAlign.py visualize --scans "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/" --panels T DC TLS --output "/Users/leohermet/Downloads/IntegrAlign-main/output/"
python main_IntegrAlign.py saveimgs --params "/Users/leohermet/Downloads/IntegrAlign-main/output/params.json" --exclude 02006 06001 08006 --rotate 01008_DC_2
python main_IntegrAlign.py align --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Users/leohermet/Documents/SPIAT/BI3/T/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/DC/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/TLS/Output/cell_positions_data_clean_CellType/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --maxMS 10
python main_IntegrAlign.py finetuning --id 16005 --meshsize 1 8 --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Users/leohermet/Documents/SPIAT/BI3/T/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/DC/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/TLS/Output/cell_positions_data_clean_CellType/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --visualization all --scans "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/"

python main_IntegrAlign.py visualize --scans "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --panels T DC TLS --output "/Users/leohermet/Downloads/IntegrAlign-main/output/"
python main_IntegrAlign.py saveimgs --params "/Users/leohermet/Downloads/IntegrAlign-main/output/params.json" --exclude 02006 06001 08006 --rotate 01008_DC_2
python main_IntegrAlign.py align --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Users/leohermet/Documents/SPIAT/BI3/T/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/DC/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/TLS/Output/cell_positions_data_clean_CellType/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --maxMS 10
python main_IntegrAlign.py finetuning --id 16005 --meshsize 1 8 --dwnscimg "/Users/leohermet/Downloads/IntegrAlign-main/output/downscaled_images.pkl" --tables "/Users/leohermet/Documents/SPIAT/BI3/T/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/DC/Output/cell_positions_data_clean_CellType/" "/Users/leohermet/Documents/SPIAT/BI3/TLS/Output/cell_positions_data_clean_CellType/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --visualization all --scans "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/"




python main_IntegrAlign.py saveimgs --params "/Users/leohermet/Downloads/IntegrAlign-main/output/params.json" --exclude 01001 01003 01008 01011 01013 01018 01022 01024 02005 02006 02009 06001 06002 06004 06006 06007 07003 08004 08006 08009 08011 08012 08014 09001 --rotate 01008_DC_2
python main_IntegrAlign.py visualize --scans "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_T/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_DC/SCANS/" "/Volumes/Bureautique/Equipes/LICL/MIFIS/Projets/BREAST IMMUN 3 (BI3)/panel_TLS/SCANS/" --annotations "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_T/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_DC/Annotations/" "/Volumes/DD_Léo_H/integrAlign/pipeline/data_all/panel_TLS/Annotations/" --panels T DC TLS --output "/Users/leohermet/Downloads/IntegrAlign-main/output/"

'''

if __name__ == "__main__":
    args = parse_arguments()

    # Process the task
    if args.task == "visualize":
        visualization(args.scans, args.annotations, args.panels, args.output)

    elif args.task == "saveimgs":
        save_downscaled_images(args.params, args.exclude)

    elif args.task == "align":
        alignment(args.dwnscimg, args.tables, args.annotations, args.resolution, args.maxMS, args.metric, args.raster, args.alpha)

    elif args.task == "finetuning":
        # Ensure scans is provided if visualization is specified
        if args.visualization != "0" and args.scans is None:
            raise ValueError("--scans is required when --visualization is set to a value other than '0'.")

        finetuning(args.id, args.meshsize, args.dwnscimg, args.tables, args.annotations, args.visualization, args.scans, args.resolution, args.metric, args.raster, args.alpha)