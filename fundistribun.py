import mdtraj as mt
import numpy as np
import os
import csv
#import matplotlib.pyplot as plt
basename="apoCzrA"
filenames=["apoCzrA.2.nc","apoCzrA.3.nc","apoCzrA.4.nc","apoCzrA.5.nc","apoCzrA.6.nc","apoCzrA.7.nc","apoCzrA.8.nc","apoCzrA.9.nc"]

iniframe=mt.load_netcdf('../'+filenames[0], top='../'+basename+'.prmtop', frame=1).topology
amides=list(iniframe.select("name N"))
namides=len(amides)
waterO=iniframe.select("water and name O")
hists=[0]*namides
chunkinfo=[]
chunksize=2000
for filename in filenames:
    for chunk in mt.iterload('../'+filename, chunk=chunksize,top='../'+basename+'.prmtop'):
        chunkinfo.append(chunk.n_frames/chunksize)
        for amide in amides:
            pairs=[[amide,x] for x in waterO]
            indice=amides.index(amide)
            if type(hists[indice]) == int:
               hists[indice]=mt.compute_rdf(chunk, pairs, r_range=None, bin_width=0.005, n_bins=None, periodic=True, opt=True)[1]
               xaxis=mt.compute_rdf(chunk, pairs, r_range=None, bin_width=0.005, n_bins=None, periodic=True, opt=True)[0]
               continue
            hists[indice]=hists[indice]+chunkinfo[-1]*mt.compute_rdf(chunk, pairs, r_range=None, bin_width=0.005, n_bins=None, periodic=True, opt=True)[1]
    print("finish file "+filename)
totchunks=sum(chunkinfo)
print(totchunks)
hist_scale=[]
for hist in hists:
    hist_scale.append(hist/totchunks)
print(hist_scale)

printable=[list(xaxis)]+[list(x) for x in hist_scale]
#with open('histograms.txt', 'a') as f1:
with open("histograms.dat", "w", newline="") as f:
    writer = csv.writer(f,delimiter=' ', lineterminator='\n')
    writer.writerow(["#distance"]+["resid"+str(amides.index(amide)+1) for amide in amides])
    for v in zip(*printable):
       writer.writerow(v)
