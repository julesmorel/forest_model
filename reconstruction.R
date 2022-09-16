remotes::install_github('umr-amap/aRchi')
library("aRchi")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("5 arguments must be supplied (input directory and output directory, the distance of research for point neighborhood, the clustering distance and the maximum searching distance for skeleton building.).n", call.=FALSE)
} else {
   inDir=args[1];
   outDir=args[2];
   D_=args[3];
   cl_dist_=args[4];
   max_d_=args[5];
}

listInputFiles = list.files(inDir)

for (file in listInputFiles) {
  print(file)
  tryCatch(
    expr = {
      input = paste(inDir, fsep = .Platform$file.sep, file, sep = "")
      outputfile=gsub(file, pattern=".xyz", replacement=".obj")
      output = paste(outDir, fsep = .Platform$file.sep, outputfile, sep = "")
      points = read.table(input)
      ar = aRchi::build_aRchi()
      ar = aRchi::add_pointcloud(ar,point_cloud = points)
      ar = skeletonize_pc(ar, D = D_, cl_dist = cl_dist_, max_d = max_d_)
      ar = smooth_skeleton(ar)
      ar =add_radius(ar)
      
      mesh=QSM2mesh(ar@QSM, tmesh = TRUE, sides = 16)
      rgl::shade3d(mesh)
      rgl::writeOBJ(output)
      rgl::close3d()
    },
    error = function(e){ 
      print(e)
    },
    finally = {}
  )
}
