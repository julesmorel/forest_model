# Forest models reconstruction method

**Reconstruction of 3D forest models from segmented point clouds**

-----------------
## Overview

This method is a complete pipeline of forest mockups reconstruction from wood/leaf segmented point clouds.
It is made of the following steps:
1. Filtering of the leaves points.
2. Outliers removal for an easier clustering.
3. Clustering of the complete scan to separate the trees.
4. Trees mesh model reconstruction

-----------------
## Setup

The method relies on [Adtree](https://github.com/tudelft3d/AdTree) to reconstruct 3D trees mesh models from point clouds representing single trees

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
## Reconstruction

The method starts by splitting the input point cloud into clusters, each of those clusters corresponding to individual tree. Then it run ADtree in Batch processing mode in order to build a 3D mockup of each individual tree.

###	Clustering

Because the tree reconstruction method uses point cloud of isolated trees, the input point cloud needs to be at first splitted into clusters.
A statistical removal filter has been added in order to filter the LiDAR noise for better identifying the clusters.
Some of the clusters may correspond to lower vegetation objects; those clusters are filtered out by checking their size along the Z axis. 

To cluster the point cloud, run the following command:

```bash
./segment_terrain.sh INPUT_FILE OUTPUT_DIR CLUSTER_TOLERANCE MIN_CLUSTER_SIZE MAX_CLUSTER_SIZE MIN_SIZE MEAN_K STD_DEV_MUL_THRESH
```

Where:

* INPUT_FILE is an ASCII file containing points stored in the following format: X Y Z label (leaf label=0).
* OUTPUT_DIR is the output directory where every tree point clouds will be stored.
* CLUSTER_TOLERANCE is the spatial cluster tolerance.
* MIN_CLUSTER_SIZE and MAX_CLUSTER_SIZE are the size limits of the clusters.
* MIN_SIZE is the minimum size along the Z axis for a cluster to be considered.
* MEAN_K and STD_DEV_MUL_THRESH are the number of neighbors to analyze for each point and the standard deviation multiplier.

We usually use:
```bash
./reconstruction INPUT_FILE OUTPUT_DIR 0.05 1000 1000000000 2.0 8 0.5
```

###	Adtree

Once the input point cloud has been clustered, simply call Adtree in batch mode:
```bash
./Adtree INPUT_DIR OUTPUT_DIR
```
* INPUT_DIR is the directory where every tree point clouds are stored.
* OUTPUT_DIR is the directory containing the resulting tree mockups.