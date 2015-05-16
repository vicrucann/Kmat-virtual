<2014>, <vicrucann@gmail.com>

To calculate stability of K parameters using the affine compensation  method, run *run_generate_synthetic.m*    
Extract the deviation when noise=0 and distortion=0 (refer to *_std variables at the end).  

Before running the main script, it is necessary to build ellipse center extractor (ellipse-params) and put the executable at /matlab/exe/ folder:  
/Kmat-virtual/matlab/exe/ellipse-params (this file will be required by detect_centers.m function) 

The source can be found:  
git clone https://github.com/vicrucann/ellipseC.git  

Compilation should be done using cmake and make commands
