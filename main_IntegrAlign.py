from src.utils.args import parse_arguments
from src.visualization import visualization
from src.save_images import save_downscaled_images
from src.alignment import alignment
from src.finetuning import finetuning

if __name__ == "__main__":
    args = parse_arguments()

    if args.task == "visualize":
        visualization(args.scans, args.annotations, args.panels, args.output)

    elif args.task == "saveimgs":
        save_downscaled_images(args.params, args.exclude, args.brightness, args.HALOrotation)

    elif args.task == "align":
        alignment(args.dwnscimg, args.tables, args.resolution, args.maxMS, args.metric, args.raster, args.alpha)

    elif args.task == "finetuning":
        if args.visualization != "0" and args.scans is None:
            raise ValueError("--scans is required when --visualization is set to a value other than '0'.")

        finetuning(args.id, args.meshsize, args.dwnscimg, args.tables, args.visualization, args.scans, args.resolution, args.metric, args.raster, args.alpha)