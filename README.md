# Deconvolution
Wiener &amp; Gauss Filter

The full analysis workflow can be exeuted using these commands:
 - Digitizer:
 ``` bash
 $ lar -c digitizer_ideal_run.fcl dune1x2x6_optical_tutorial_sim_gen_3GeV.root
 $ lar -c digitizer_template_run.fcl dune1x2x6_optical_tutorial_sim_gen_3GeV.root
 ```    
     - This should output: opdetraw_ideal_gen/hist.root & opdetraw_template_gen/hist.root
 
 - Deconvolution:
   - $ lar -c deconvolution_gauss_run.fcl opdetraw_template_gen.root
     - This should output: deconv_gauss_gen.root & deconv_gauss_hist.root
 
   - $ lar -c deconvolution_wiener_run.fcl opdetraw_template_gen.root
     - This should output: deconv_wiener_gen.root & deconv_wiener_hist.root
 
 - OpHitFinder:
   - $ lar -c ophitfinder_ideal_run.fcl opdetraw_ideal_gen.root
     - This should output: ophitspe_ideal_gen.root & ophitspe_ideal_hist.root

   - $ lar -c ophitfinder_raw_run.fcl opdetraw_template_gen.root
     - This should output: ophitspe_raw_gen.root & ophitspe_raw_hist.root
 
   - $ lar -c ophitfinder_gauss_run.fcl deconv_gauss_gen.root
     - This should output: ophitspe_gauss_gen.root & ophitspe_gauss_hist.root
 
   - $ lar -c deconvolution_wiener_run.fcl deconv_wiener_gen.root
     - This should output: ophitspe_wiener_gen.root & ophitspe_wiener_hist.root
 
 Then all the *_hist.root files can be exported to be used in the analysis notebook.

Enjoy!
