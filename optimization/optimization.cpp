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

#include "../util/pointCloudFileReader.h"
#include "../util/offsetManager.h"


#include <sys/stat.h>

int main(int argc, char *argv[]){
  std::string directoryOut,filenameIn;
  double smoothnessThreshold = 45.;
  double curvatureThreshold = 10.;
  int minClusterSize=1000;
  int maxClusterSize=10000000000;
  double minSize=2.0;
  int meanK=8;
  double stddevMulThresh=0.5;

  if (argc == 2) {
    filenameIn = argv[1];
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

  //Optimization parameters
  double smoothnessMin = 0.;
  double smoothnessMax = 0.;
  double smoothnessStep = 2.;

  double curvatureMin = 0.0;
  double curvatureMax = 1.0;
  double curvatureStep = 0.05;

  for(int m=0;m<(int)((smoothnessMax-smoothnessMin)/smoothnessStep);m++){
    for(int n=0;n<(int)((curvatureMax-curvatureMin)/curvatureStep);n++){
        double smoothness=smoothnessMin+m*smoothnessStep;
        double curvature=curvatureMin+n*curvatureStep;

        pcl::RegionGrowing<pcl::PointXYZI, pcl::Normal> reg;
        reg.setMinClusterSize (minClusterSize);
        reg.setMaxClusterSize (maxClusterSize);
        reg.setSearchMethod (treeN);
        reg.setNumberOfNeighbours (10);
        reg.setInputCloud (ptsFilteredOutliers.makeShared());
        reg.setInputNormals (normals);
        reg.setSmoothnessThreshold (smoothness / 180.0 * M_PI);
        reg.setCurvatureThreshold (curvature);
        reg.setSmoothModeFlag(true);
        reg.setCurvatureTestFlag(true);
        reg.extract (cluster_indices);
        
        int clusterCounter=0;
        int clusterNb=0;
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
                clusterCounter+=pointsCluster.size();
                clusterNb++;
            }
        }
        double ratio=((double)clusterCounter)/((double)pointsFiltered.size());
        std::string filename = "optimResult.txt";
        std::ofstream outfile;
        outfile.open(filename, std::ios_base::app);
        outfile <<smoothness<<" "<<curvature<<" "<<ratio<<" "<<clusterNb<<std::endl;
    }
  }
}