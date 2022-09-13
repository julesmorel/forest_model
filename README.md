# Forest models reconstruction method

**Reconstruction of 3D forest models from segmented point clouds**

-----------------
## Overview

This method is a complete pipeline of forest mockups reconstruction from wood/leaf segmented point clouds.
It is made of the following steps:
1. Filtering of the leaves points.
2. Outliers removal for an easier clustering.
3. Clustering of the complete scan to isolate the trees.
4. Trees mesh model reconstruction

-----------------
## Setup

The method relies on [Adtree](https://github.com/tudelft3d/AdTree), [TreeQSM](https://github.com/InverseTampere/TreeQSM) or [aRchi](https://github.com/umr-amap/aRchi) to reconstruct 3D trees mesh models from point clouds representing single trees 

### Requirements
* Linux (tested on Ubuntu 20.04.4 LTS)
* PCL 1.11

### Install
Install this library by running the following command:
```shell
cmake .
make
```  

-----------------
## Modeling of forest plots

The method starts by splitting the input point cloud into clusters, each of those clusters corresponding to individual tree. Then the user runs one of the proposed methods (AdTree, TreeQSM or aRchi) in order to build a 3D mockup of each individual tree.

###	Clustering

Because the tree reconstruction method uses point cloud of isolated trees, the input point cloud needs to be at first splitted into clusters.
A statistical removal filter has been added in order to filter the LiDAR noise for better identifying the clusters.
Some of the clusters may correspond to lower vegetation objects; those clusters are filtered out by checking their size along the Z axis. 

To cluster the point cloud, run the following command:

```bash
./segment_terrain.sh INPUT_FILE OUTPUT_DIR SOOMTHNESS_TOLERANCE CURVATURE_THRESHOLD MIN_CLUSTER_SIZE MAX_CLUSTER_SIZE MIN_SIZE MEAN_K STD_DEV_MUL_THRESH
```

Where:

* INPUT_FILE is an ASCII file containing points stored in the following format: X Y Z label (leaf label=0).
* OUTPUT_DIR is the output directory where every tree point clouds will be stored.
* SOOMTHNESS_TOLERANCE: If the deviation between points normals is less than the smoothness tolerance (in degrees) then they are suggested to be in the same cluster 
* CURVATURE_THRESHOLD: if this curvature at the current point is less than the curvature threshold then the algorithm will continue the growth of the cluster using the newly added point.
* MIN_CLUSTER_SIZE and MAX_CLUSTER_SIZE are the size limits of the clusters.
* MIN_SIZE is the minimum size along the Z axis for a cluster to be considered.
* MEAN_K and STD_DEV_MUL_THRESH are the number of neighbors to analyze for each point and the standard deviation multiplier.

We usually use:
```bash
./reconstruction INPUT_FILE OUTPUT_DIR 45. 10. 1000 1000000000 2.0 8 0.5
```

###	Reconstruction

This package contains 3 reconstruction methods designed to produce 3D mesh models from individual scanned trees.  

[**AdTree**]

First, compile AdTree by running the following command:
```bash
cd AdTree
mkdir build
cd build
cmake ..
make
```

Once the input point cloud has been clustered, simply call Adtree in batch mode:
```bash
./Adtree INPUT_DIR OUTPUT_DIR
```
* INPUT_DIR is the directory where every tree point clouds are stored.
* OUTPUT_DIR is the directory containing the resulting tree mockups.

[**TreeQSM**]

To use TreeQSM, call in Matlab the function `run.m` as follow:
```bash
run(LIST_OF_FILES,OUTPUT_DIR)
```
* LIST_OF_FILES is a glob containing the ASCII files resulting from the clustering
* OUTPUT_DIR is the directory containing the resulting tree mockups.

For instance, one can run:

```bash
run("~/GEDI/GEDI008/ascii/*.xyz","~/GEDI/GEDI008/qsm/")
```

[**aRchi**]

To use aRchi, run the following command:
```bash
Rscript reconstruction.R INPUT_DIR OUTPUT_DIR
```
* INPUT_DIR is the directory where every tree point clouds are stored.
* OUTPUT_DIR is the directory containing the resulting tree mockups.

For instance, one can run:
```bash
Rscript reconstruction.R ~/GEDI/GEDI009/trees ~/GEDI/GEDI009/aRchi
```
