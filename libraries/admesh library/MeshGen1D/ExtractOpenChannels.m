function [FD,A] = ExtractOpenChannels(DEMFile)

DEM = GRIDobj(DEMFile);
DEMf = fillsinks(DEM);

FD = FLOWobj(DEMf);
A = flowacc(FD);

end