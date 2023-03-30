# Deconvolution
Wiener &amp; Gauss Filter

The full analysis workflow can be exeuted using these commands:
 - Digitizer:
   - Run the following commands to output: opdetraw_ideal_gen/hist.root & opdetraw_template_gen/hist.root
 ``` shell
 $ lar -c digitizer_ideal_run.fcl dune1x2x6_optical_tutorial_sim_gen_3GeV.root
 ```
 ``` shell
 $ lar -c digitizer_template_run.fcl dune1x2x6_optical_tutorial_sim_gen_3GeV.root
 ``` 
 
 - Deconvolution:
   - Run the following commands to output: deconv_gauss_gen/hist.root & deconv_wiener_gen/hist.root
``` shell
$ lar -c deconvolution_gauss_run.fcl opdetraw_template_gen.root
```
``` shell
$ lar -c deconvolution_wiener_run.fcl opdetraw_template_gen.root
```      
 - OpHitFinder:
   - Run the following commands to output: ophitspe_ideal_gen/hist.root & ophitspe_raw_gen/hist.root & ophitspe_gauss_gen/hist.root & ophitspe_wiener_gen/hist.root
``` shell
$ lar -c ophitfinder_ideal_run.fcl opdetraw_ideal_gen.root
```
``` shell
$ lar -c ophitfinder_raw_run.fcl opdetraw_template_gen.root
```
``` shell
$ lar -c ophitfinder_gauss_run.fcl deconv_gauss_gen.root
```
``` shell 
$ lar -c deconvolution_wiener_run.fcl deconv_wiener_gen.root
```

 Then all the *_hist.root files can be exported to be used in the analysis notebook.

Enjoy!
