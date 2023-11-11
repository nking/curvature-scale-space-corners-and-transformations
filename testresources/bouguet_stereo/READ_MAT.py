# from https://github.com/malcolmreynolds/calib_bouguet
# dependencies: python3, opencv-python
## copy calib_bouguet/build/load.py into this directory 
import load as calib_bouguet

results = calib_bouguet.load_calib_results('full path to Calib_Results_right.mat')

print("intrinsics:")
print(results.intrinsics.matrix)
print("extrinsics:")
#dir(results.extrinsics)
for i in range(14) :
  print("i=", i, " ", results.extrinsics.__getitem__(i))

