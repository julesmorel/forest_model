#include<iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/filters/passthrough.h>
#include <pcl/common/common.h>
#include <pcl/filters/extract_indices.h>
#include <boost/make_shared.hpp>

#include "util/pointCloudFileReader.h"
#include "util/offsetManager.h"


#include <sys/stat.h>

int main(int argc, char *argv[]){
  std::string directoryOut,filenameIn;
  double clusterTolerance = 1.;
  int minClusterSize=10;
  int maxClusterSize=1000;
  double minSize=2.0;
  int meanK=16;
  double stddevMulThresh=1.0;

  if (argc == 9) {
    filenameIn = argv[1];
    directoryOut = argv[2];
    clusterTolerance = std::stod (argv[3]);
    minClusterSize = std::stoi (argv[4]);
    maxClusterSize = std::stoi (argv[5]);
    minSize = std::stod (argv[6]);
    meanK = std::stoi (argv[7]);
    stddevMulThresh = std::stod (argv[8]);
  }else{
    std::cout<<"Please specify a file to process"<<std::endl;
    exit(0);
  }

  pcl::PointCloud<pcl::PointXYZI> points;  
  //read the offset for the current file if it exists
  offsetManager offsetM(filenameIn);
  double offset_x = offsetM.getOffsetX();
  double offset_y = offsetM.getOffsetY();
  //read the points in the file 
  pointCloudFileReader::read(filenameIn,points,offset_x,offset_y);
  std::cout<<"Input point cloud: "<<points.size()<<" points"<<std::endl;
  //store of the offset if it was not stored before
  offsetM.setOffset(offset_x,offset_y);

  //Filtering the points labelized as "leaf"
  pcl::PointCloud<pcl::PointXYZI> pointsFiltered;  
  pcl::PassThrough<pcl::PointXYZI> pass;
  pass.setInputCloud (points.makeShared());
  pass.setFilterFieldName ("intensity");
  pass.setFilterLimits (0.001, 1.0);
  pass.filter (pointsFiltered);
  std::cout<<"Wood point cloud: "<<pointsFiltered.size()<<" points"<<std::endl;

  //Filtering the outliers to clean up the cloud in order to better separate the clusters
  pcl::PointCloud<pcl::PointXYZI> ptsFilteredOutliers;
  pcl::StatisticalOutlierRemoval<pcl::PointXYZI> sor;
  sor.setInputCloud (pointsFiltered.makeShared());
  sor.setMeanK (meanK);
  sor.setStddevMulThresh (stddevMulThresh);
  sor.filter (ptsFilteredOutliers);
  std::cout<<pointsFiltered.size()-ptsFilteredOutliers.size()<<" points filtered (StatisticalOutlierRemoval)"<<std::endl;

  // Creating the KdTree object for the search method of the extraction
  pcl::search::KdTree<pcl::PointXYZI>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZI>);
  tree->setInputCloud (ptsFilteredOutliers.makeShared());

  //Creating the clusters of points
  std::vector<pcl::PointIndices> cluster_indices;
  pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec;
  ec.setClusterTolerance (clusterTolerance); 
  ec.setMinClusterSize (minClusterSize);
  ec.setMaxClusterSize (maxClusterSize);
  ec.setSearchMethod (tree);
  ec.setInputCloud (ptsFilteredOutliers.makeShared());
  ec.extract (cluster_indices);

  std::cout<<cluster_indices.size()<<" clusters detected"<<std::endl;

  //Saving the clusters as point clouds in the given directory
  int clusterCounter=0;
  for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin (); it != cluster_indices.end (); ++it)
  {
    int idCluster = it-cluster_indices.begin();

    //check if the cluster is tall enough
    pcl::PointCloud<pcl::PointXYZI> pointsCluster;  
    for(int i=0;i<it->indices.size();i++){
      pointsCluster.push_back(ptsFilteredOutliers.at(cluster_indices.at(idCluster).indices[i]));
    }
    pcl::PointXYZI minPt, maxPt;
    pcl::getMinMax3D (pointsCluster, minPt, maxPt);
    double deltaZ = maxPt.z - minPt.z;

    if(deltaZ>minSize){     
      std::string filename = directoryOut + "/tree_" + std::to_string(idCluster) + ".xyz";
      std::ofstream outfile;
      outfile.open(filename, std::ios_base::app);
      for(int i=0;i<it->indices.size();i++){
        pcl::PointXYZI currentPt = ptsFilteredOutliers.at(cluster_indices.at(idCluster).indices[i]);
        outfile <<currentPt.x<<" "<<currentPt.y<<" "<<currentPt.z<<std::endl;
      }
      clusterCounter++;
    }
  }
  std::cout<<clusterCounter<<" tree clusters saved"<<std::endl;
}
