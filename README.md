# LargeVolumeTransformix
Transforming giga- to teravoxel high-resolution volumes with low-resolution registration results. 

# Background
Virtual histology based on hard X-ray microtomography produces volumetric imaging with isotropic pixel sizes down to and below one micron. As field-of-view increases, the resulting datasets grow in size from the gigabyte to terabyte (and even petabyte) scale. Registering and transforming these datasets requires a multi-resolution distributed approach.

As a first step, we developed a pipeline for transforming large virtual histology images using low-resolution registrations. In this case, we collected microtomography of an entire mouse brain with isotropic 0.65 µm pixel size and developed a pipeline to reconstruct the data, see [G. Rodgers et al. " Mosaic microtomography of a full mouse brain with sub-micron pixel size "](https://doi.org/10.1117/12.2633556). This data should be transformed into the Allen Mouse Brain Common Coordinate Framework.

The approach is described in more detail in [C. Tanner et al. "Registration of microtomography images: challenges and approaches"](https://doi.org/10.1117/12.2633922). Briefly, the coordinates region defined in the fixed image is transformed to the moving frame, a bounding box is defined to encompass this skewed region, then the high resolution data is loaded and warped, see image below. For large volumes, the fixed volume is divided into grids, each of which can fit into memory and can be transformed in parallel, then the final volume is formed by combining the transformed sub-regions.

![Coordinate transforms for transformix](https://github.com/grodgers1/LargeVolumeTransformix/blob/main/example/figures/fig_coordtransform.png)

# Usage
Set parameters and filenames in 'DataParameterDefinition.m'. You will need to adjust paths to the data, the resolution and region of interest which should be processed, and the maximum filesize your computer can handle. You will also need to adjust the path to the initial transform for the affine transform given, if relative paths change. More information is found in `DataParameterDefinition.m`.

Run the script `TransformToAtlas.m'. It will first state the selected main parameters.

## Requirements
The pipeline relies on [elastix and transformix](https://elastix.lumc.nl/) are for registration and transformation.

## Data
### Low-resolution volumes for registration
After publication, the datasets will hopefully be made publicly available.

For now, the 32x binned reco and the atlas volumes are on our group's storage:
`/storage/groups/bmc/shared_projects/LargeVolumeTransformix/registration/volumes/`

The atlas volumes were downloaded from [https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/) then re-oriented using the following code:

```
[avgt, ~] = nrrdread([atlas_directory 'average_template_10.nrrd']); # from file exchange 
avgt = permute(avgt,[1,3,2]);
avgt = int16(double(avgt)-2^15);
avgt = flip(avgt,3);
avgt = rot90(avgt,1);

```

For registration, the pixel sizes were divided by 100. The pixel size adjustment allows for non-rigid registration with a bending energy penalty, as the original pixel size results in values below the detection limit.

### High-resolution reconstructions
After publication, the datasets will hopefully be made publicly available.

For now, the 8x binned recos are available on our group's storage:
`/storage/groups/bmc/shared_projects/LargeVolumeTransformix/reconstructions/mouse4_eth/`
