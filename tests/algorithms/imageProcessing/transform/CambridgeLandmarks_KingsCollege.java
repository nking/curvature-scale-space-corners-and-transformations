package algorithms.imageProcessing.transform;

/**
 
PoseNet: A Convolutional Network for Real-Time 6-DOF Camera Relocalization
Alex Kendall, Matthew Grimes, Roberto Cipolla
Proceedings of the International Conference on Computer Vision, 2015
http://mi.eng.cam.ac.uk/projects/relocalisation/

Description:
King's College scene from Cambridge Landmarks, a large scale outdoor visual relocalisation dataset taken around Cambridge University. Contains original video, with extracted image frames labelled with their 6-DOF camera pose and a visual reconstruction of the scene. If you use this data, please cite our paper: Alex Kendall, Matthew Grimes and Roberto Cipolla "PoseNet: A Convolutional Network for Real-Time 6-DOF Camera Relocalization." Proceedings of the International Conference on Computer Vision (ICCV), 2015. 

Contents:
Images (.png), Video (.mp4), Text (.txt), Reconstructions (.nvm can be opened with Visual SFM, http://ccwu.me/vsfm/ )

using features from 2 images:
#ImageFile,          Camera Position [X Y Z W P Q R]
seq7/frame00039.png  -13.182436 -18.442966 1.724784 0.678307 0.569025 -0.308146 0.348075
seq8/frame00074.png   13.236864 -24.132282 1.828876 0.769271 0.637893 -0.030718 0.019230

Spatial Extent (m) of entire KC dataset is 140 x 40m
Dist. to Conv. 2-3m, 2.70-3 degrees
* 
* License:
* https://creativecommons.org/licenses/by-nc-sa/2.0/uk/
* Attribution-NonCommercial-ShareAlike 2.0 UK: England & Wales (CC BY-NC-SA 2.0 UK)
* 
* Publication Reference:
http://mi.eng.cam.ac.uk/projects/relocalisation/
* 
@author nichole
 */
public class CambridgeLandmarks_KingsCollege {
    
}
