# OTC

This repository includes the R package to perform optimal transport colocalization analyses of super resolution data. 

## Reproducibility
To reproduce a figure in the paper Colocalization for Super-Resolution Microscopy via Optimal Transport (C. Tameling et al, 2021) just run the script located under Code with the figure name.

## Installation
To use the OTC package on your own data, we recommend to install CPLEX to significantly speed up the computation. Therefore sign up for the acemdic initiative at IBM and download CPLEX. Aferwards download the transport package from CRAN and adapt the makevars file to include CPLEX. With this liniking install the transport package and afterwards the OTC package from this repository. We tested the package under Ubuntu and Windows. 

## Usage

To use the OTC package on your own data, please follow the description and have a look at the demo.R file.

1. Compute transport plans: To do so,  use the function calculate_tplans(*data_path, picsA, picsB, random_sections, n_random_sections, output_path, output_name*) from the OTC package. Here, *data\_path* is the directory where the images that shall be analyzed are located, *picsA* is a list of filenames of TIFF images of protein A and *picsB* for images of protein B,  respectively. Make sure to set *random_sections* to TRUE if your images are bigger than 128x128 pixels and select how many random samples you want to draw (*n_random_sections*). The calculated transport plans will be saved in one file with name *output_name* in the directory specified in *output_path*.
2. Calculate the OTC curves from the transport plans: User the function evaluate_tplans(*data_path, data_list, pxsize, dim, output_path, output_name*). Here, *data_path* is the directory where the file containing all transport plans is located, *data_list* is the list of filenames for these transport plan files, which provides a list of all pairwise OTC curves which will be included into one plot. The size of one pixel in nm and the dimension of the evaluated images needs to be specified in *pxsize* and *dim* respectively, *output_path* and *output_name* are again the directory and the name for the results file. This function returns a data frame with OTC curves.
3. Plot OTC curves: Use plot_otc_curves(*otc_curves, output_path, output_name*) for OTC package. The data frame with OTC curves entering the variable *otc_curves* are the returning value from step 2 go into *otc_curves*. The  figure is automatically saved under *output_name* in the directory *output_path*.

