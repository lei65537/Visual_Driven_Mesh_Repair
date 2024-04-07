# Repairing Man-Made Meshes via Visual Driven Global Optimization with Minimum Intrusion
![](docs/teaser.png)
Lei Chu, Hao Pan, Yang Liu, Wenping Wang. ACM Transactions on Graphics (SIGGRAPH Asia 2019).
## Dependency
      Boost(1.68)  https://sourceforge.net/projects/boost/files/boost-binaries/1.68.0/
      CGAL(4.12.1)  https://github.com/CGAL/cgal/releases
      MOSEK(8.1)  https://www.mosek.com/downloads/8.1.0.81/
      Eigen3(included)
      libigl(included)
      Freeglut(3.0.0)(included)
      GLM(0.9.9.3)(included)
      GLEW(2.1.0)(included)
      
   
## Installation
Our project was developed with VS2015(x64) on Windows10. Following is how to install on Windows.
* Clone to local repository: git clone https://github.com/lei65537/Visual_Driven_Mesh_Repair.git
* Compile the code with CMake(>=3.3)
* Please set(under windows) MOSEK_INCLUDE_DIR as MOSEK_ROOT_PATH\8\tools\platform\win64x86\h and MOSEK_LIBRARIES as MOSEK_ROOT_PATH/8/tools/platform/win64x86/bin/mosek64_8_1.lib
* The license of MOSEK is required.
* Please tick the _Advanced_ and _Grouped_ in CMake
* Please copy the dlls(including _libgmp-10.dll_ and _libmpfr-4.dll_ of CGAL) to the folder of the executable before running

## Citation
If you use our code or model, please cite our paper:

      @Article{Chu:2019:MeshRepair, 
        Title = {Repairing Man-Made Meshes via Visual Driven Global Optimization with Minimum Intrusion}, 
        Author = {Chu, Lei and Pan, Hao and Liu, Yang and Wang, Wenping}, 
        Journal = {ACM Trans. Graph. (SIGGRAPH ASIA)}, 
        issue_date = {November 2019},
        Year = {2019}, 
        Number = {6}, 
        Volume = {38},
        Pages={158:1--158:18},
        numpages = {18},
        DOI={https://doi.org/10.1145/3355089.3356507},
        Publisher = {ACM} 
      }

## Usage
### Input & Output
We support .off and .obj for input and the output is .obj. There will be three output meshes for each input model.
* _vd.obj -- the result mesh after visual driven processing
* _opt.obj -- the result manifold mesh after gloabl topology optimization
* _opt_HR.obj -- the more concise manifold mesh with hidden patch removal

### Call the Program
MeshViewer.exe PATH_TO_INPUT_FILE Eps_Visual Csl_min W1 W2 Num_threads

* PATH_TO_INPUT_FILE - the path to a triangle mesh in .off or .obj
* Eps_Visual and Csl_min - parameters that control the visual driven processing step
* W1 and W2 - the weights of E_ori and E_reg in the manifold mesh reconstruction step
* Num_threads - controls the number of cpu threads used

Suggested default parameter values are given in the following example and discussed in the paper/supplemental document.

### Test Examples

For example: MeshViewer.exe .\\PATH_TO_Test_examples\\airplane_0627.off 0.0056 0.08 5 1 8

## Processed Datasets

We provide 10k processed meshes of ModelNet and ShapeNet by our program. The meshes are named according to their sources in the respective datasets.
* ModelNet processed (1.6GB): https://connecthkuhk-my.sharepoint.com/:u:/g/personal/wenhua00_connect_hku_hk/ETGIi2UtszRBjzd2P65sPg0BMDkBXFz8nqIAXynuC49H3g?e=3kliil
* ShapeNet processed (5GB): https://connecthkuhk-my.sharepoint.com/:u:/g/personal/wenhua00_connect_hku_hk/EYJOeGqcIslJgqzLE4mAoFIB3ky6FPn2uRWK9RMGOEU2Rw?e=YzUDkm

## License
The code and data are released under the MIT license.
## Acknowledgement
