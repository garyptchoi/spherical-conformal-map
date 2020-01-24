spherical_conformal_map: Conformally map a genus-0 closed triangle mesh to the unit sphere

This code computes spherical conformal (i.e. angle-preserving) maps of genus-0 closed triangle meshes using the linear method in [1], which has been applied for human brain mapping, texture mapping, surface registration, cardiac mapping and so on.
Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following paper:
[1] P. T. Choi, K. C. Lam, and L. M. Lui, 
    "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
    SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.

Copyright (c) 2013-2018, Gary Pui-Tung Choi
https://scholar.harvard.edu/choi

===============================================================



Usage:
map = spherical_conformal_map(v,f)


Input:
v: nv x 3 vertex coordinates of a genus-0 triangle mesh
f: nf x 3 triangulations of a genus-0 triangle mesh

Output:
map: nv x 3 vertex coordinates of the spherical conformal map