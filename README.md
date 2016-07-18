# FRB-model-selection

Code to implement robust model selection for linear regression models based on MM-estimators and the 
[Fast and Robust Bootstrap] (http://dx.doi.org/10.1214/aos/1021379865)
as described in [Salibian-Barrera, M. and Van Aelst, S. (2008)](http://dx.doi.org/10.1016/j.csda.2008.05.007). 


You will need to create a dynamic library from the code in `FRB-model-selection.c` using,
for example, the following command in your shell:
```R
R CMD SHLIB FRB-model-selection.c
```
Note that if you are running Windows, you will need to have installed the [RTools package](https://cran.r-project.org/bin/windows/Rtools/)
from CRAN. 

The file `FRB-model-selection-example.R` contains a script to apply this method to the 
well-known Boston Housing data set. 

