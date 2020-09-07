# TaskRelated

The code in this repository contains all the MATLAB software used to analyze the data and generate the figures for the following publication: 
Roth, ZN, Ryoo, M, and Merriam, EP (2020). Task-related activity in human visual cortex. (DOI: XXX).

This code depends heavily on two other software packages: [mrTools](https://github.com/justingardner/mrTools), and [mgl](https://github.com/justingardner/mgl).
Additional packages used are: [gru](https://github.com/justingardner/gru), 
and [Circular Statistics Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).

To generate the figures in the paper, download the data from OSF ([DOI: 10.17605/OSF.IO/CBJQ6](https://doi.org/10.17605/osf.io/cbjq6)) into a single folder, and run the following functions. 
If things don’t work out, please contact Zvi Roth ([zvi.roth@mail.huji.ac.il](mailto:zvi.roth@mail.huji.ac.il)) and Eli Merriam ([elisha.merriam@nih.gov](mailto:elisha.merriam@nih.gov)).

Order in which to run scripts:

(retrieve data)
* savePupilData.m
* saveTaskEyeData.m
* saveTaskData.m
* saveTaskData_localizer.m
* saveControlData.m

(generate figures)
* fig2.m
* fig3.m
* fig4.m 
* fig5.m
* fig6.m
* fig7.m
* figS1B.m
* figS2.m
* figS3.m

(statistical significance of results)
* permutationTest.m
* permutationTest_pupil.m
