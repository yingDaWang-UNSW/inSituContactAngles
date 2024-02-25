# inSituContactAngles
In Situ Contact Angle Measurements on Voxelized Images

This codebase is companion to the publication "In situ Characterization of Heterogeneous Surface Wetting in Porous Materials" by YD Wang et al./, published in Advances in Colloids and Interface Sciences.

It computes the local contact angle - mapped to the contact loop of fluid clusters in segmented voxelised images containing 2 fluid phases and 1 solid phase. 

Users can compute the contact angle distributions of their whole domain, and/or interpolate the spatial contact angle map over their solid domain.

This codebase uses 3 methods of measuring the insitu contact angle:

1. Geometric method (Alratrout, 2018) https://doi.org/10.1016/j.advwatres.2017.07.018
2. Topological method (Sun, 2020) https://doi.org/10.1016/j.jcis.2020.05.076
3. Extended Topological method (this study)

To run this code and its example, see the following steps:


0. Users should first open the provided 3D tiff image file "segmentedRockTestTheta30.tif" in imageJ or another visualiser to better understand the data that this codebase works with as input 

1. Open matlab and navigate to the inSituContactAngles folder

2. Run "runThisFile.m" in Matlab. This will plot the contact angle distributions as the computation runs, but will also at the end dump a database "thetaData0.mat". This contains all the information required to reconstruct the results and perform further analysis.

3. Run "generateSpatialMap.m" in Matlab. This will take "thetaData0.mat" and perform spatial analysis. It will create both 3D tiff image files containing the contact angle in the voxels corresponding to the contact loops of the fluid-solid phases as well as the 3D interpolation of the values in these contact loops over the entire solid phase.

