# LargeVolumeTransformix
Transforming giga- to teravoxel high-resolution volumes with low-resolution registration results

## Requirements
[elastix and transformix](https://elastix.lumc.nl/) are needed.

## Data
### Low-resolution volumes for registration
After publication, the datasets will hopefully be made publicly available.

For now, the 32x binned reco and the atlas volumes are on our group's storage:
`/storage/groups/bmc/shared_projects/LargeVolumeTransformix/registration/volumes/`

The atlas volumes were re-oriented and the pixel sizes were divided by 100. The pixel size adjustment allows for non-rigid registration with a bending energy penalty, as the original pixel size results in values below the detection limit.

### High-resolution reconstructions
After publication, the datasets will hopefully be made publicly available.

For now, the 8x binned recos are available on our group's storage:
`/storage/groups/bmc/shared_projects/LargeVolumeTransformix/reconstructions/mouse4_eth/`

## Usage
Follow the script `TransformToAtlas.m`

You will need to adjust paths to the data. You will also need to adjust paths to the initial transform for the affine transform given. More information is found in `TransformToAtlas.m`.

