function run(inputFiles,outputDir)
listing = dir(inputFiles);

for i=1:length(listing)
    file=listing(i).name;
    [pathstr, name, ext] = fileparts(file);
    directory=listing(i).folder;
    filename=directory+"/"+file;

    points = load(filename, "-ASCII");
    inputs = define_input(points,1,1,1);
    inputs.GrowthVolCor=1;
    inputs.plot=0;
    inputs.savetxt=0;
    inputs.savemat=0;
    inputs.saveobj=0;
    inputs.disp=0;
    QSM = treeqsm(points,inputs);

    outputFilename=outputDir+name+".obj";
    QSMCyl = QSMBCylindrical(QSM.cylinder.start,QSM.cylinder.axis,QSM.cylinder.length,QSM.cylinder.radius,QSM.cylinder.parent,QSM.cylinder.branch,QSM.cylinder.added);
    QSMCyl.export( ...
        'OBJ', ...
        outputFilename, ...
        'TriangleStem', ...
        'FaceCount',[5 10], ...
        'Closed' ...
    );
end    
end