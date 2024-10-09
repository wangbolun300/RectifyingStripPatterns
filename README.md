# Webs via Propagation
<img src="./data/fig/ppl.jpg"  />  

This is an implementation of the paper "Computational design of asymptotic geodesic hybrid gridshells via propagation algorithms, CAD 2024". The code is implemented in C++, and tested on MSVC and GCC. This is developed based on the framework of [Rectifying Strip Patterns](https://github.com/wangbolun300/RectifyingStripPatterns), but we enhanced the abilities for gridshell design. Using our code, you can start from a existing regular web (the distribution of the vertices is in the form of a $m\times n$ grid), a strip (constructed by $2\times n$ vertices), or even a polyline ($n$ vertices), and gradually expand the web into a AAG, AGG, or GGG until reaching the desired size. 

## Citation

If you use our code in your project, please consider citing the [original paper](TODO):

```bibtex
@article{Wang:2024:WebsPropagation,
    title        = {Computational design of asymptotic geodesic hybrid gridshells via propagation
algorithms},
    author       = {Wang, Bolun and Almaskin, Maryam and Pottmann, Helmut},
    year         = 2024,
    month        = xx,
    journal      = {Computer-Aided Design},
    volume       = XX,
    number       = XX,
    articleno    = xx,
    numpages     = XX
}
```

## Compiling Instruction 
To compile the code, first you need to install CMake (https://cmake.org/), 
To build the executable on Linux or macOS:
```sh
cd WebsViaPropagation/
mkdir build
cd build
cmake ../  -DCMAKE_BUILD_TYPE=Release
make
```
To build it on MSVC, you need to download [openmesh](https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh/-/jobs/156362/artifacts/raw/OpenMesh-9.1.0.zip) into the "external" folder, unzip it into "external/openmesh/OpenMesh-9.1.0" and build it into "external/openmesh/OpenMesh-9.1.0/build".

Then double-click the executable `lsc_devbin`, you will see a pop-out window whose buttons link the functions that generate all the results in our paper. 


## Usage
Some useful shortcuts:
* `i` is to invert the normal directions of the viewer.
* `d` is to enable double-sided rendering.
* `x` is to delete the selected mesh in the `Mesh Management` list.

Some tips:
* Please always check if there is any new mesh added into the `Mesh Management` list. Sometimes after pressing the buttons, there is no change in the main view window, but the new mesh is already added into the list. To see the renderings of the new results you need to make other meshes invisible by clicking the checker boxes on the left of the mesh files.
* Please always check the information printed in the terminal, it will tell you the numerical errors, how to read/write files after pressing the corresponding buttons, etc.


### The GUI.
<img src="./data/fig/GUI.png" style="zoom:80%;" /> 

Ignoring the buttons inherited from Rectifying Strip Patterns, the buttons or inputs that are useful for our project are marked in the 4 squares in the figure. They are
1. the parameters. The `weight boundary` is $\lambda_{appro}$ in our paper. `weight laplacian` = $\lambda_{fair}$, `weight pseudo-geodesic` is the weights for $E_{ggg}$, $E_{agg}$ or $E_{aag}$, which is a constant number 1.0 in the paper. `WeightAngle` = $\lambda_{guide}$ for the guide curves.
2. the file list. Every time when the optimization is finished, the result file is shown here. 
3. the buttons for loading, propagating, and optimizing webs. 
4. the `Save Mesh` button. Clicking on the target mesh in the `Mesh Management` list, and save the mesh by clicking `Save Mesh`.

Next is the basic pipeline for the usage.
#### Load files. 
Click on `ReadPlyObj`, you will be able to read a initial curve where you can start your propagation. If you want to read a strip (the size is $2\times n$) or a web (the size is $m\times n$), you need to first specify the parameter `vinrow` = $n$, which is the number of vertices added to the existing web in each propagation. Then click on `easyCheckXXX`, where `XXX` is either `AGG`, `AAG` or `GGG`.  
#### Propagate. 
For AAG, first `AAGPlanIntersect` to compute the rectifying strip, then click on `AAGAdjustInit` 3~5 times to update the propagated vertices, each clicking does one round of coordinate-decent optimization, and click on `AAGPropagate` to finish the propagation. For AGG, first `AGGCurve2Strip`, then 3~5 times of `AGGAdjustInit`, finally `AGGPropagate`. For GGG, First `GGGurve2Strip`, then `GGGPropagate`, finally `OptGGG`. Note that when you are using `easyCheckAAG` to load a strip or a web as input, you need to click on `Calibrate` before `AAGPlanIntersect`.

The parameters `AggPar1`~`AggPar5` are for slightly changing the first strip constructed by loading an initial curve using `ReadPlyObj`. Please read the code to understand how these 5 parameters work.
#### optimize.
* To load a reference curve for AGG, use `ReadRefCurve`. If you have loaded a curve as the initialization, you can also apply rigid transformations to the reference curve by playing with the parameters `CurveAngle` and `CurveRotate`, and use `TransRotateAGG` to apply the transformations. The final curve can be saved using `SaveTransRotate`. If the reference curve is too sparse, you can use `upsampleCurve` and `smoothCurve` to increase the density before saving the curve.
* `ChooseDiag`. The first thing is to choose which diagonal is the third family before global optimization. The different diagonals are due to specifying $P_i$ and $P_{i-1}$ propagate vertex $P_{i+n}$ or vertex $P_{i-1+n}$. The former corresponds to `D0`, and the later corresponds to `D1`. Unfortunatelly, our propagation algorithms don't provide selection from the two diagonals: for GGG, we provide `D0`, and for AGG and AAG, we provide `D1`. Please choose the correct diagonal before global optimization. You can plot the diagonals on the quad mesh using `DrawDiags`. We believe that it won't be hard for the developers to complete the selection of diagonals, since the changes are minor.
* By default, the shape approximate the first curve. If you loaded a strip or a web, and you want to approximate the first strip of it, please check the box `ApproStrip`
* Run optimization by clicking `OptAAG`, `OptAGG` or `OptGGG`, and check the printed information in the terminal to check the quality.
### Interactive Design
* Isometric deformation. 1. Load a web using `vinrow` and `easyCheckIso`. 2. Choose the web, hold number key "2", and use the left button of the mouse to choose a vertex. 3. Hold number key "5", select a position using the left button of the mouse. 4. run optimization by `OptIso`. Note that approximation to the chosen point and the first curve / strip use the same weight `\lambda_{appro}`. If you want to cancel the approximation to the first strip while keeping approximating to the selected point, please  go into the cpp files change the code.
* Steering the guide curve for AGG. First, select a point on the curve using `selectCurvePtID` and `selectCurvePt`. Then hold the number key "5", select a position using the left button of the mouse. After clicking on `EditCurve`, you will see the optimized curve. You can save it by `SaveEditedCurve`.
* Interactive design of AAG and AGG and GGG is similar to the isometric deformation, but with different loading and optimization buttons.




## Some Test Data
We provide one initial curve such that you can test propagation with loading polylines. We provide one AAG, one AGG (along with the guide curve) and one GGG for you to test the propagation and optimizaiton algorithms.
