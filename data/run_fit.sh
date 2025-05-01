valor-run_fit \
   --config /sps/juno/ycao/VALOR/softw/JUNO/vmaster/config/analysis/Atmo_Juno_nuebar_2d.xml \
   --config-datasetB /sps/juno/ycao/VALOR/softw/JUNO/vmaster/config/analysis/Atmo_Juno_nuebar_2d_NO.xml \
   --physics-hypothesis void \
   --param-values-option prefit-errors \
   --osc-prior-file priors/void.xml \
   --load-prior-from-file - \
   --fitter minuit+with_simplex_fallback \
   --mesg-thresholds verbosity.xml \
   --without-stat-fluct  \
   --tag numu_test_NO \
   --run 1 \
   --ntoymc 1 \

 