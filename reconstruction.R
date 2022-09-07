remotes::install_github('umr-amap/aRchi')
library("aRchi")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two arguments must be supplied (input directory and output directory).n", call.=FALSE)
} else {
   inDir=args[1];
   outDir=args[2];
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
      ar = skeletonize_pc(ar)
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
    warning = function(w){
      print(w)
    },
    finally = {}
  )
}
