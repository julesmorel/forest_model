import logging
import numpy as np
import subprocess
import click

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        # logging.FileHandler("debug.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('slim')

@click.command()
@click.option('--segmentation_result_file', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False), help='the ASCII file resulting from the deep learning model segmentation.')
@click.option('--output_directory_ascii', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), help='this directory stores the point clouds of the clusters.')
@click.option('--output_directory_obj', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), help='this directory stores the tree models in the OBJ format.')
@click.option('--smooth_tol', type=float, default=45., help='If the deviation between points normals is less than the smoothness tolerance (in degrees) then they are suggested to be in the same cluster.')
@click.option('--curv_th', type=float, default=10., help='If this curvature at the current point is less than the curvature threshold then the algorithm will continue the growth of the cluster using the newly added point.')
@click.option('--cluster_size_min', type=int, default=1000, help='minimum number of points to form a cluster  ')
@click.option('--cluster_size_max', type=int, default=1000000000, help='maximum number of points to form a cluster   ')
@click.option('--minsize_cluster_Z', type=float, default=2., help='the minimum size along the Z axis to form a cluster')
@click.option('--mean_k', type=int, default=8, help=' the number of neighbors to analyze for each point for the outlier filter')
@click.option('--sd_mul_dev', type=float, default=0.5, help='the standard deviation multiplier for the outlier filter')
@click.option('--archi_D', type=float, default=0.3, help='the distance of research for point neighborhood')
@click.option('--archi_cl_dist', type=float, default=0.1, help='the clustering distance ')
@click.option('--archi_max_d', type=float, default=0.05, help='the maximum searching distance for skeleton building ')
def reconstruction_with_aRchi(segmentation_result_file, output_directory_ascii, output_directory_obj, smooth_tol, curv_th, cluster_size_min, cluster_size_max, minsize_cluster_Z, mean_k, sd_mul_dev, archi_D, archi_cl_dist, archi_max_d):
    """Assign class=2 to the ground points, class=3 to the wood points and class=1 to the remaining point in the .laz files provided"""

    #Clustering
    cmd = ['./reconstruction', segmentation_result_file, output_directory_ascii, str(smooth_tol), str(curv_th), str(cluster_size_min), str(cluster_size_max), str(minsize_cluster_Z), str(mean_k), str(sd_mul_dev)]
    print(' '.join(cmd))
    subprocess.check_call(' '.join(cmd), shell=True)

    #Reconstruction of tree models from clusters
    cmd = ['Rscript', 'reconstruction.R', output_directory_ascii, output_directory_obj, str(archi_D), str(archi_cl_dist), str(archi_max_d)]
    print(' '.join(cmd))
    subprocess.check_call(' '.join(cmd), shell=True)

if __name__ == '__main__':
    reconstruction_with_aRchi()    