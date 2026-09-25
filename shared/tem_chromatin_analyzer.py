import json
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname
from warnings import warn

import numpy as np
import pandas as pd
import skimage
import tifffile as tiff
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  # Need to bypass default limit to get image to load


def process_folder(source_dir, output_dir, **kwargs):

    source_dir = Path(source_dir)

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Process folder to look for qpproj files
    qpproj_path_list = list(source_dir.rglob("*.qpproj"))

    all_data = []

    for file in qpproj_path_list:
        image_uris = get_image_uri_from_qupath_project(file)

        rel_path = (file.parent).relative_to(source_dir)

        for image_name, image_uri in image_uris.items():
            # Look for the mask file
            mask_path = next(
                (source_dir / "export").glob(image_name + "-labels.png"), None
            )

            if not mask_path:
                warn(f"[{image_name}] No label file found.")
                continue

            data = process_image(image_uri, mask_path, output_dir / rel_path)

            if data is not None:
                all_data.append(data)

    df = pd.DataFrame(all_data)
    df.to_csv(output_dir / "summary_data.csv", index=False)


def process_image(filepath, maskpath, output_dir, radius=25, write_csv=True):

    filepath = Path(filepath)

    print(filepath)

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Read the image
    with Image.open(filepath) as image:
        image = np.array(image.convert("L"))

    # Read the nuclear mask from QuPath
    mask_nucl = skimage.io.imread(maskpath)

    if not np.any(mask_nucl):
        warn(f"{filepath.stem}: No nuclear mask was found in label file.")
        return

    # Create a blank array to hold the filtered image
    # image_filt = np.empty_like(image, dtype=np.uint8)

    # ndimage.gaussian_filter(image, sigma=7.0, output=image_filt)
    # ndimage.median_filter(image, size=(10, 10), output=image_filt)

    mask_hchromatin = segment_hchromatin(image)

    # Clean up mask
    mask_hchromatin = mask_hchromatin & mask_nucl

    # Now segment the ring
    mask_nucl_ring = (
        mask_nucl - skimage.morphology.isotropic_erosion(mask_nucl, radius=radius)
    ) > 0

    # from matplotlib import pyplot as plt

    # plt.imshow(mask_nucl_ring)
    # plt.show()
    # plt.close()

    mask_hchromatin_ring = mask_nucl_ring & mask_hchromatin

    # ---Measure data---
    data = {
        "file": str(filepath),
    }

    # -- Chromatin ring properties --

    props_hchromatin = skimage.measure.regionprops(
        mask_hchromatin_ring.astype(np.uint8)
    )
    props_ring = skimage.measure.regionprops(mask_nucl_ring.astype(np.uint8))

    data["hchromatin_area"] = props_hchromatin[0].area
    data["ring_area"] = props_ring[0].area
    data["average_hchromatin_area"] = data["hchromatin_area"] / data["ring_area"]

    # -- Nuclear properties --
    props_nucleus = skimage.measure.regionprops(mask_nucl.astype(np.uint8))

    data["nuclear_area"] = props_nucleus[0].area
    data["nuclear_equiv_diameter_area"] = props_nucleus[0].equivalent_diameter_area
    data["nuclear_eccentricity"] = props_nucleus[0].eccentricity

    data["nuclear_circularity"] = (4 * np.pi * props_nucleus[0].area) / (
        props_nucleus[0].perimeter ** 2
    )

    data["nuclear_solidity"] = props_nucleus[0].solidity

    hull_mask = props_nucleus[0].image_convex.astype(np.uint8)
    hull_perimeter = skimage.measure.perimeter(hull_mask)
    data["nuclear_tortuosity"] = hull_perimeter / props_nucleus[0].perimeter

    # ---Save data---

    # Write raw masks to file
    fn = filepath.stem
    skimage.io.imsave(
        output_dir / (fn + "_hchromatin_full_mask.tif"),
        mask_hchromatin,
        check_contrast=False,
    )
    skimage.io.imsave(
        output_dir / (fn + "_hchromatin_ring_mask.tif"),
        mask_hchromatin_ring,
        check_contrast=False,
    )

    images = (image, mask_nucl, mask_hchromatin_ring)
    channel_names = ["TEM", "Nuclear mask", "Heterochromatin mask"]

    write_combined_image_stack(images, output_dir, channel_names=channel_names)

    if write_csv:
        df = pd.DataFrame(data, index=[0])
        df.to_csv(output_dir / (fn + "_data.csv"), index=False)

    return data


def segment_hchromatin(image):

    # Measure the dark regions (heterochromatin) in the image
    th = skimage.filters.threshold_li(image)

    mask_hchromatin = image < (0.7 * th)

    return mask_hchromatin


def write_combined_image_stack(images, output_dir, channel_names=None):

    try:
        stack = np.stack(images, axis=0)
    except ValueError as e:
        print(f"Stacking failed: {e}")
        for f, img in zip(image_list, images):
            print(f"{f}: shape = {img.shape}, dtype = {img.dtype}")
        return

    if channel_names is None:
        channel_names = [f"Channel {i + 1}" for i in range(stack.shape[0])]

    # print(channel_names)
    # exit()

    tile_size = 256
    subresolutions = 2

    with tiff.TiffWriter(output_dir / "combined.ome.tif", bigtiff=True) as tif:
        options = {
            "photometric": "minisblack",
            "tile": (tile_size, tile_size),
            "compression": "lzw",
            "resolutionunit": "CENTIMETER",
        }

        tif.write(
            stack,
            subifds=subresolutions,
            metadata={"axes": "CYX", "Channel": {"Name": channel_names}},
            **options,
        )

        for level in range(subresolutions):
            mag = 2 ** (level + 1)
            tif.write(
                stack[:, ::mag, ::mag],
                subfiletype=1,  # FILETYPE.REDUCEDIMAGE
                **options,
            )


def get_image_uri_from_qupath_project(qpproj_path):

    qpproj_path = Path(qpproj_path)

    # Get the QuPath project file in the parent directory
    if qpproj_path.is_dir():
        project_file = next(qpproj_path.glob("*.qpproj"), None)
    elif qpproj_path.is_file():
        project_file = qpproj_path

    if project_file is None:
        raise FileNotFoundError(
            f"Could not find a QuPath project (.qpproj) file at path: {qpproj_path}"
        )

    # Read QuPath project and append the image files to a dictionary. The image file
    # name is used as a key so it is easy to match with the label filename later.
    image_uri_dict = {}
    with open(project_file, "r", encoding="utf-8") as f:
        project_data = json.load(f)

        for entry in project_data.get("images", []):
            image_name = entry.get("imageName")
            server_builder = entry.get("serverBuilder", {})
            uri = server_builder.get("builder", {}).get("uri")

            # print(f"{image_name}: Server {uri}")

            # Only map it if both the name and URI exist
            if image_name and uri:
                image_name_clean = Path(image_name).stem

                # Clean the uri to remove the "file:///"
                parsed_url = urlparse(uri)
                local_filepath = url2pathname(parsed_url.path)

                image_uri_dict[image_name_clean] = local_filepath

    return image_uri_dict
