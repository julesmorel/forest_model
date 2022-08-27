#include<iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/filters/passthrough.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/common/common.h>
#include <pcl/filters/extract_indices.h>
#include <boost/make_shared.hpp>

#include "util/pointCloudFileReader.h"
#include "util/offsetManager.h"


#include <sys/stat.h>

int main(int argc, char *argv[]){
  std::string directoryOut,filenameIn;
  double smoothnessThreshold = 45.;
  double curvatureThreshold = 10.;
  int minClusterSize=10;
  int maxClusterSize=1000;
  double minSize=2.0;
  int meanK=16;
  double stddevMulThresh=1.0;

  if (argc == 10) {
    filenameIn = argv[1];
    directoryOut = argv[2];
    smoothnessThreshold = std::stod (argv[3]);
    curvatureThreshold = std::stod (argv[4]);
    minClusterSize = std::stoi (argv[5]);
    maxClusterSize = std::stoi (argv[6]);
    minSize = std::stod (argv[7]);
    meanK = std::stoi (argv[8]);
    stddevMulThresh = std::stod (argv[9]);
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

  //Estimating the normals of the remaining points
  std::vector<pcl::PointIndices> cluster_indices;
  pcl::search::Search<pcl::PointXYZI>::Ptr treeN (new pcl::search::KdTree<pcl::PointXYZI>);
  pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
  pcl::NormalEstimation<pcl::PointXYZI, pcl::Normal> normal_estimator;
  normal_estimator.setSearchMethod (treeN);
  normal_estimator.setInputCloud (ptsFilteredOutliers.makeShared());
  normal_estimator.setKSearch (10);
  normal_estimator.compute (*normals);

  //Creating the clusters of points
  pcl::RegionGrowing<pcl::PointXYZI, pcl::Normal> reg;
  reg.setMinClusterSize (minClusterSize);
  reg.setMaxClusterSize (maxClusterSize);
  reg.setSearchMethod (treeN);
  reg.setNumberOfNeighbours (10);
  reg.setInputCloud (ptsFilteredOutliers.makeShared());
  reg.setInputNormals (normals);
  reg.setSmoothnessThreshold (smoothnessThreshold / 180.0 * M_PI);
  reg.setCurvatureThreshold (curvatureThreshold);
  reg.extract (cluster_indices);

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