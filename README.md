# allen CCF tools

Some code to work with the Allen Inst CCF data, specifically the 10µm voxel 2016 version. 

## Usage example:
```
>> tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
>> av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
>> st = loadStructureTree('structure_tree_safe.csv'); % a table of what all the labels mean
>> f = allenAtlasBrowser(tv, av, st);
```

## Requirements
You need the npy-matlab repository to load the data: https://github.com/kwikteam/npy-matlab

You also need the data files. See //zserver/Lab/Atlas/allenCCF or, if you don't have access to that, contact me at nick.steinmetz@gmail.com. 

## Note about annotation volume
The original volume has numbers that correspond to the "id" field in the structure tree, but since I wanted to make a colormap for these, I re-indexed the annotation volume by the row number of the structure tree. So in this version the values correspond to "index"+1. This also allows using uint16 datatype, cutting file size in half. See setup_utils.m.


## (unsorted comments)
## Plot wire mesh of brain (example usage in script_sliceMovie)
bregma = allenCCFbregma();
isBrain = av>1; % >0 for original av, >1 for by_index
gridIn3D(double(isBrain), 0.5, 50, bregma);
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30    25]);
