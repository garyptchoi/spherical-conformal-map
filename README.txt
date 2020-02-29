spherical_conformal_map: Conformally map a genus-0 closed triangle mesh to the unit sphere

This code computes spherical conformal (i.e. angle-preserving) maps of genus-0 closed triangle meshes using the linear method in [1], which has been applied for human brain mapping, texture mapping, surface registration, cardiac mapping and so on.
Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following papers:
[1] P. T. Choi, K. C. Lam, and L. M. Lui, 
    "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
    SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.

(For mobius_area_correction_spherical)
[2] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
    "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
    SIAM Journal on Imaging Sciences, 2020.

Copyright (c) 2013-2020, Gary Pui-Tung Choi
https://scholar.harvard.edu/choi

===============================================================
Usage (see demo.m):
map = spherical_conformal_map(v,f)

Input:
v: nv x 3 vertex coordinates of a genus-0 triangle mesh
f: nf x 3 triangulations of a genus-0 triangle mesh

Output:
map: nv x 3 vertex coordinates of the spherical conformal map

===============================================================
Some possible extensions are also provided (see demo_extension.m):
- spherical_area_preserving_map: 
We replace the cotangent Laplacian in our proposed method by the locally authalic Laplacian. Note that this approach aims to reduce the local area distortion, while the global area distortion is not necessarily minimized.

- iterative_spherical_area_preserving_map:
We use the above modified method to obtain an initial spherical map, and then solve an energy minimization problem to obtain a more area-preserving map. Again, note that the energy only focuses on minimizing the local area distortion, while the global area distortion is not necessarily minimized.

- mobius_area_correction_spherical:
We further reduce the global area distortion of a spherical conformal parameterization while maintaining the conformality, using the Mobius area correction method in [2]. 
[2] G. P. T. Choi, Y. Leung-Liu, X. Gu, and L. M. Lui, 
    "Parallelizable global conformal parameterization of simply-connected surfaces via partial welding."
    SIAM Journal on Imaging Sciences, 2020.
