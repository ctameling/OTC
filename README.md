# OTC

This repository includes all code and data to reproduce all figures in the paper.

To use the OTC package on your own data, please follow the description.

1. Compute transport plans: To do so,  use the function calculate_tplans(*data_path, picsA, picsB, random_sections, n_random_sections, output_path, output_name*) from the OTC package. Here, *data\_path* is the directory where the images that shall be analyzed are located, *picsA* is a list of filenames of TIFF images of protein A and *picsB* for images of protein B,  respectively. Make sure to set *random_sections* to TRUE if your images are bigger than 128x128 pixels and select how many random samples you want to draw (*n_random_sections*). The calculated transport plans will be saved in one file with name *output_name* in the directory specified in *output_path*.
2. Calculate the OTC curves from the transport plans: User the function evaluate_tplans(*data_path, data_list, pxsize, dim, output_path, output_name*). Here, *data_path* is the directory where the file containing all transport plans is located, *data_list* is the list of filenames for these transport plan files, which provides a list of all pairwise OTC curves which will be included into one plot. The size of one pixel in nm and the dimension of the evaluated images needs to be specified in *pxsize* and *dim* respectively, *output_path* and *output_name* are again the directory and the name for the results file. This function returns a data frame with OTC curves.
3. Plot OTC curves: Use plot_otc_curves(*otc_curves, output_path, output_name*) for OTC package. The data frame with OTC curves entering the variable *otc_curves* are the returning value from step 2 go into *otc_curves*. The  figure is automatically saved under *output_name* in the directory *output_path*.

