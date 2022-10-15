import logging
import numpy as np
import subprocess

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        # logging.FileHandler("debug.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('slim')

def reconstruction_with_aRchi(segmentation_result_file, output_directory_ascii, output_directory_obj, smooth_tol=45., curv_th=10., cluster_size_min=1000, cluster_size_max=1000000000, minsize_cluster_Z=2., mean_k=8, sd_mul_dev=0.5, archi_D=0.3, archi_cl_dist=0.1, archi_max_d=0.05):
    '''Assign class=2 to the ground points, class=3 to the wood points and class=1 to the remaining point in the .laz files provided
    
    Parameters
    ----------
    segmentation_result_file : string
        the ASCII file resulting from the deep learning model segmentation
    output_directory_ascii : string
        this directory stores the point clouds of the clusters
    output_directory_obj : list of string
        this directory stores the tree models in the OBJ format
    smooth_tol : float
        If the deviation between points normals is less than the smoothness tolerance (in degrees) then they are suggested to be in the same cluster
    curv_th : float
        if this curvature at the current point is less than the curvature threshold then the algorithm will continue the growth of the cluster using the newly added point.
    cluster_size_min : int
        minimum number of points to form a cluster  
    cluster_size_max : int
        maximum number of points to form a cluster      
    minsize_cluster_Z : float
        the minimum size along the Z axis to form a cluster
    mean_k : int
        the number of neighbors to analyze for each point 
    sd_mul_dev : float
        the standard deviation multiplier         
    archi_D : float
        the distance of research for point neighborhood
    archi_cl_dist : float
       the clustering distance 
    archi_max_d : float
        the maximum searching distance for skeleton building                             
    '''

    #Clustering
    cmd = ['./reconstruction', segmentation_result_file, output_directory_ascii, str(smooth_tol), str(curv_th), str(cluster_size_min), str(cluster_size_max), str(minsize_cluster_Z), str(mean_k), str(sd_mul_dev)]
    print(' '.join(cmd))
    subprocess.check_call(' '.join(cmd), shell=True)

    #Reconstruction of tree models from clusters
    cmd = ['Rscript', 'reconstruction.R', output_directory_ascii, output_directory_obj, str(archi_D), str(archi_cl_dist), str(archi_max_d)]
    print(' '.join(cmd))
    subprocess.check_call(' '.join(cmd), shell=True)