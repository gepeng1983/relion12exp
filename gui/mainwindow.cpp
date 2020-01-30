/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#include "mainwindow.h"

RelionMainWindow::RelionMainWindow(int w, int h, const char* title):Fl_Window(w,h,title)
{

	// Initialisation
	run_button = NULL;
	print_CL_button = NULL;
	cite_button = NULL;
	ori_w = w;
	ori_h = h;

	// Height of each entry
    step_y = 22;

    color(GUI_BACKGROUND_COLOR);
    int menuheight = 30;
    int tabheight = 25;
    menubar = new Fl_Menu_Bar(0, 0, w, menuheight);
    menubar->add("File/Load settings",  FL_ALT+'l', cb_menubar_load, this);
    menubar->add("File/Save settings",  FL_ALT+'s', cb_menubar_save, this);
    menubar->add("File/Reactivate Run",  FL_ALT+'r', cb_menubar_reactivate_runbutton, this);
    menubar->add("File/About",  FL_ALT+'a', cb_menubar_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_menubar_quit, this);
    current_y += menuheight;

    // Dropdown for runtype choice
    choice_runtype = new Fl_Choice(150, 4, 177, 22);
    choice_runtype->menu(runtype_options);
    choice_runtype->label("Run type:");
    // By default new run
    run_type = PREPROCESS;
    choice_runtype->picked(&runtype_options[0]);
    choice_runtype->callback(cb_menu_runtype, this);
    menu_runtype = choice_runtype;
    menu_runtype->color(GUI_BACKGROUND_COLOR);

    // Dropdown for new/continue choice
    choice_continue = new Fl_Choice(330, 4, 167, 22);
    choice_continue->menu(continue_options);
    // By default new run
    is_continue = false;
    choice_continue->picked(&continue_options[0]);
    choice_continue->callback(cb_menu_continue, this);
    menu_continue = choice_continue;
    menu_continue->color(GUI_BACKGROUND_COLOR);


    restart_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    restart_group->end();

    // Set up tabs
    begin();
    tabs = new Fl_Tabs(0, menuheight, w, h);

    tab1 = new Fl_Group(0, menuheight + tabheight, w, h, "I/O");
    tab1->end();

    tab2 = new Fl_Group(0, menuheight + tabheight, w, h, "CTF");
    tab2->hide();
    // fill tab here
    ctf_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    ctf_group->end();
    ctffind_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    ctffind_group->end();
    tab2->end();

    tab3 = new Fl_Group(0, menuheight + tabheight, w, h, "Optimisation");
    tab3->hide();
    star_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    extract_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    movie_extract_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    movie_extract_group->end();
    movie_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    movie_group->end();
    extract_group->end();
    star_group->end();
    tab3->end();

    tab4 = new Fl_Group(0, menuheight + tabheight, w, h, "Sampling");
    tab4->hide();
    // fill tab here
    localsearch_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    localsearch_group->end();
    rescale_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    rescale_group->end();
    norm_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    norm_group->end();

    tab4->end();


    tab5 = new Fl_Group(0, menuheight + tabheight, w, h, "Running");
    tab5->hide();
    // fill tab here
    queue_group = new Fl_Group(0, menuheight + tabheight, XCOL5, h);
    queue_group->end();
    tab5->end();

    tabs->end();

    menubar->color(GUI_BACKGROUND_COLOR);
    tabs->color(GUI_BACKGROUND_COLOR);
    tab1->color(GUI_BACKGROUND_COLOR);
    tab2->color(GUI_BACKGROUND_COLOR);
    tab3->color(GUI_BACKGROUND_COLOR);
    tab4->color(GUI_BACKGROUND_COLOR);
    tab5->color(GUI_BACKGROUND_COLOR);
    tabs->selection_color(GUI_BACKGROUND_COLOR);
    tab1->selection_color(GUI_BACKGROUND_COLOR2);
    tab2->selection_color(GUI_BACKGROUND_COLOR2);
    tab3->selection_color(GUI_BACKGROUND_COLOR2);
    tab4->selection_color(GUI_BACKGROUND_COLOR2);
    tab5->selection_color(GUI_BACKGROUND_COLOR2);
    start_y = menuheight + tabheight + 5;
    current_y = start_y;


   end();
   resizable(tabs);

   //show();
}
void RelionMainWindow::clearPreprocess()
{
	fn_out_preprocess.clear();
	cs.clear();
	kv.clear();
	q0.clear();
	mag.clear();
	dstep.clear();

	do_ctffind.clear();
	mic_names.clear();
        ctf_win.clear();
	box.clear();
	resmin.clear();
	resmax.clear();
	dfmin.clear();
	dfmax.clear();
	dfstep.clear();
	fn_ctffind3_exe.clear();

	do_extract.clear();
	coord_names.clear();
	coord_format.clear();
	extract_size.clear();
	do_movie_extract.clear();
	movie_rootname.clear();
	avg_movie_frames.clear();
	first_movie_frame.clear();
	last_movie_frame.clear();

	do_rescale.clear();
	rescale.clear();
	do_norm.clear();
	bg_radius.clear();
	white_dust.clear();
	black_dust.clear();
	do_invert.clear();
	do_starfile.clear();


}
void RelionMainWindow::setupPreprocess()
{
	// Clear all the refine stuff
	clearRefine();

	// Fill first tab with I/O
	tab1->begin();
	tab1->copy_label("I/O");

	resetHeight();

	AddAnyEntry(fn_out_preprocess, "Output rootname:", "particles", "Output rootname. This should NOT contain a directory structure.");

	// Add a little spacer
	current_y += step_y/2;

	AddSliderEntry(cs,"Spherical aberration (mm):", 2, 0, 8, 0.1, "Spherical aberration of the microscope used to collect these images (in mm)");

	AddSliderEntry(kv,"Voltage (kV):", 300, 50, 500, 10, "Voltage the microscope was operated on (in kV)");

	AddSliderEntry(q0,"Amplitude contrast:", 0.1, 0, 0.3, 0.01, "Fraction of amplitude contrast. Often values around 10% work better than theoretically more accurate lower values...");

	AddSliderEntry(mag,"Magnification at detector (x):", 60000, 10000, 200000, 1000, "Magnification (at the detector) used to collect these images (in times) ");

	AddSliderEntry(dstep,"Pixel size of detector (um):", 14, 1, 32, 1, "Pixel size of the detector (in micrometer) ");

	tab1->end();

	// tab2 holds the CTF corrections (and some other corrections as well)
	tab2->begin();
	tab2->copy_label("CTFFIND");
	resetHeight();

	AddBooleanEntry(do_ctffind, "Run CTFFIND3?", true, "If set to Yes, CTFs will be estimated for all selected micrographs. \
Niko Grigorieff's program CTFFIND3 will be used for this.", ctffind_group);

	ctffind_group->begin();

	AddAnyEntry(mic_names,"Input micrographs :", "Micrographs/*.mrc", "Filenames of the micrograph(s) on which to run CTFFIND3. This may contain wildcards * and ?. \
Note that the micrographs should be in a subdirectory (e.g. called Micrographs/) of the project directory, i.e. the directory from where you are launching the GUI. \
If this is not the case, then make a symbolic link inside the project directory to the directory where your micrographs are stored.");

	AddSliderEntry(box,"FFT box size (pix):", 512, 64, 1024, 8, "CTFFIND3's Box parameter");
	AddSliderEntry(resmin,"Minimum resolution (A):", 100, 10, 200, 10, "CTFFIND3's ResMin parameter");
	AddSliderEntry(resmax,"Maximum resolution (A):", 7, 1, 20, 1, "CTFFIND3's ResMax parameter");
	AddSliderEntry(dfmin,"Minimum defocus value (A):", 5000, 0, 25000, 1000, "CTFFIND3's dFMin parameter");
	AddSliderEntry(dfmax,"Maximum defocus value (A):", 50000, 20000, 100000, 1000, "CTFFIND3's dFMax parameter");
	AddSliderEntry(dfstep,"Defocus step size (A):", 500, 200, 2000, 100,"CTFFIND3's FStep parameter");

	// Add a little spacer
	current_y += step_y/2;

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_CTFFIND3_EXECUTABLE");
	if (default_location == NULL)
		default_location=DEFAULTCTFFINDLOCATION;
	AddFileNameEntry(fn_ctffind3_exe, "CTFFIND3 executable:", default_location, "*.exe", "Location of the CTFFIND3 executable. You can control the default of this field by setting environment variable RELION_CTFFIND3_EXECUTABLE.");

	// Add a little spacer
	current_y += step_y/2;

        AddSliderEntry(ctf_win,"Estimate CTF on window size (pix) ", -1, -16, 4096, 16, "If a positive value is given, a squared window of this size at the center of the micrograph will be used to estimate the CTF. This may be useful to exclude parts of the micrograph that are unsuitable for CTF estimation, e.g. the labels at the edge of phtographic film. \n \n The original micrograph will be used (i.e. this option will be ignored) if a negative value is given.");

	ctffind_group->end();
	do_ctffind.cb_menu_i();

	tab2->end();

	// tab2 holds the CTF corrections (and some other corrections as well)
	tab3->begin();
	tab3->copy_label("extract");
	resetHeight();

	AddBooleanEntry(do_starfile, "Generate particle STAR file?", true, "If set to Yes, all metadata will be joined into one STAR file that may be used directly for \
subsequent 2D and 3D refinements in RELION.", star_group);

	star_group->begin();

	AddAnyEntry(coord_names,"Particle coordinate files: ", "Micrographs/*.box", "Filenames of the coordinate files to be used for particle extraction. This may contain wildcards * and ?.");
	AddRadioEntry(coord_format, "Coordinate file format:", coord_format_options, &coord_format_options[0], "Format of the coordinate file");

	// Add a little spacer
	current_y += step_y/2;

	AddBooleanEntry(do_extract, "Extract particles from micrographs?", true, "If set to Yes, particles will be extracted from the micrographs using all selected coordinate files. \
Niko Grigorieff's program CTFFIND3 will be used for this.", extract_group);

	extract_group->begin();

	AddSliderEntry(extract_size,"Particle box size :", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");

	// Add a little spacer
	current_y += step_y/2;

	AddBooleanEntry(do_movie_extract,"Extract from movies?", false, "If set to yes, then particles will be extracted from all frames of the MRC stacks that hold the movies.\n \
The name of the MCR stacks should be the rootname of the micrographs + '_movierootname.mrcs', where the movierootname is given below.", movie_extract_group);

	movie_extract_group->begin();

	AddAnyEntry(movie_rootname, "Rootname of movies files:", "movie", "rootname to relate each movie to the single-frame averaged micropgraph. With a rootname of 'movie', the movie for mic001.mrc should be called mic001_movie.mrcs");

	AddSliderEntry(avg_movie_frames, "Number of frames to average: ", 1, 1, 20, 1, "Particles will be extracted from averages of the specified number of movie frames. \
This reduces computational loads of subsequent refinements. Use 1 for no averaging. \
We used individual movie frames (no averaging) for our eLife ribosome structures. However, a K2 or a DE detector has a much faster frame rate than our Falcon.\
It may be of little use to consider individual frames if there is not at least 1 electron per squared Angstrom dose in it. \
Therefore, although we have not tried this yet (as of January 2013), it may be good to average over multiple frames for faster cameras...");

	AddSliderEntry(first_movie_frame, "First movie frame to extract: ", 1, 1, 20, 1, "Extract from this movie frame onwards. The first frame is number 1.");

	AddSliderEntry(last_movie_frame, "Last movie frame to extract: ", 0, 0, 64, 1, "Extract until this movie frame. Zero means: extract all frames in the movie");

	movie_extract_group->end();
	extract_group->end();
	star_group->end();

	do_starfile.cb_menu_i();
	do_extract.cb_menu_i();
	do_movie_extract.cb_menu_i();

	tab3->end();

	// tab2 holds the CTF corrections (and some other corrections as well)
	tab4->begin();
	tab4->copy_label("operate");
	resetHeight();
	AddTextOnlyEntry(autosample_text1,"These operations are only performed upon extraction.");

	AddBooleanEntry(do_rescale, "Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.", rescale_group);
	rescale_group->begin();
	AddSliderEntry(rescale, "Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number");
	rescale_group->end();
	do_rescale.cb_menu_i();

	// Add a little spacer
	current_y += step_y/2;

	AddBooleanEntry(do_norm, "Normalize particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it.", norm_group);
	norm_group->begin();
	AddSliderEntry(bg_radius, "Radius background circle (pixels): ", 60, 1, 256, 2, "Pixels outside a circle with this radius will be used to estimate the mean and stddev of the noise, which is used to normalise the images.");
	AddSliderEntry(white_dust, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	AddSliderEntry(black_dust, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	norm_group->end();
	do_norm.cb_menu_i();

	// Add a little spacer
	current_y += step_y/2;
	AddBooleanEntry(do_invert, "Invert contrast?", false, "If set to Yes, the contrast in the particles will be inverted.");


	tab4->end();


}
void RelionMainWindow::clearRefine()
{
	fn_img.clear();
	fn_out.clear();
	fn_cont.clear();
	nr_classes.clear();
	fn_ref.clear();
	ref_correct_greyscale.clear();
	sym_group.clear();
	sym_nr.clear();

	angpix.clear();
	do_ctf_correction.clear();
	ctf_corrected_ref.clear();
	ctf_phase_flipped.clear();
	ctf_intact_first_peak.clear();

	ini_high.clear();
	nr_iter.clear();
	tau_fudge.clear();
	particle_diameter.clear();
	do_zero_mask.clear();
	fn_mask.clear();
	do_movies.clear();
	fn_movie_star.clear();
	movie_runavg_window.clear();
	movie_sigma_angles.clear();
	//movie_sampling.clear();
	movie_sigma_offset.clear();
	//movie_offset_step.clear();

	autosample_text1.clear();
	psi_sampling.clear();
	sampling.clear();
	offset_range.clear();
	offset_step.clear();
	do_local_ang_searches.clear();
	sigma_angles.clear();
	auto_local_sampling.clear();

}

void RelionMainWindow::setupRefine(int newruntype, bool _is_continue)
{
	// Clear all preprocess stuff
	clearPreprocess();

	// Fill first tab with I/O
	tab1->begin();
	tab1->copy_label("I/O");

	resetHeight();

	AddFileNameEntry(fn_img, "Input images STAR file:", "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");

	std::string defaultname;
	if (run_type == RUN2D)
		defaultname = "Class2D/run1";
	else if (run_type == RUN3D)
		defaultname = "Class3D/run1";
	else
		defaultname = "Refine3D/run1";

	AddAnyEntry(fn_out, "Output rootname:", defaultname.c_str(), "Output rootname for all files of this run. \
If this rootname contains a directory structure (e.g. 20110724/run1), the directory (20110724) will be created if it does not exist.");

	AddFileNameEntry(fn_cont, "Continue from here: ", "", "STAR Files (*_optimiser.star)", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	// Add a little spacer
	current_y += step_y/2;

	if (run_type == RUN2D || run_type == RUN3D)
	{
		AddSliderEntry(nr_classes, "Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");

		// Add a little spacer
		current_y += step_y/2;
	}
	else
	{
		nr_classes.clear();
	}

	// Add a little spacer
	current_y += step_y/2;

	if (run_type == RUN3D || run_type == RUN3DAUTO)
	{
		AddFileNameEntry(fn_ref, "Reference map:", "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
Make sure this map has the same dimensions and the same pixel size as your input images.");

		AddBooleanEntry(ref_correct_greyscale, "Ref. map is on absolute greyscale?", false, "Probabilities are calculated based on a Gaussian noise model, \
which contains a squared difference term between the reference and the experimental image. This has a consequence that the \
reference needs to be on the same absolute intensity grey-scale as the experimental images. \
RELION and XMIPP reconstruct maps at their absolute intensity grey-scale. \
Other packages may perform internal normalisations of the reference density, which will result in incorrect grey-scales. \
Therefore: if the map was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. \
If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the first iteration, \
and prior to the second iteration the map will be filtered again using the initial low-pass filter. \
This procedure is relatively quick and typically does not negatively affect the outcome of the subsequent MAP refinement. \
Therefore, if in doubt it is recommended to set this option to No.");

		// Add a little spacer
		current_y += step_y/2;

		AddRadioEntry(sym_group, "Symmetry group:", symgroup_options, &symgroup_options[0], "If the molecule is asymmetric, \
set Symmetry group to C, and symmetry number to 1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

		AddSliderEntry(sym_nr, "Symmetry number:", 1, 1, 20, 1, "If the molecule is asymmetric, set Symmetry group to C, and symmetry number to 1. \
For the octahedral and tetrahedral symmetry groups, this number will be ignored. \
Note that for the icosahedral symmetry groups there are multiple settings allowed: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n ");

	}
	else
	{
		fn_ref.clear();
		ref_correct_greyscale.clear();
		sym_group.clear();
		sym_nr.clear();
	}

	tab1->end();

	// tab2 holds the CTF corrections (and some other corrections as well)
	tab2->begin();
	tab2->copy_label("CTF");

	resetHeight();

	AddSliderEntry(angpix,"Pixel size (A):", 1, 0.1, 10, 0.01, "Nominal pixel size in Angstrom for the input images. \n \n\
IMPORTANT: if CTF information is provided in the input STAR file, make sure that this pixel size matches the one used for the estimation of the CTF parameters!");

	// Add a little spacer
	current_y += step_y/2;

	AddBooleanEntry(do_ctf_correction, "Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)", ctf_group);

	ctf_group->begin();

	if (run_type == RUN3D || run_type == RUN3DAUTO)
	{
		AddBooleanEntry(ctf_corrected_ref, "Has reference been CTF-corrected?", false, "Set this option to Yes if the reference map \
represents density that is unaffected by CTF phases and amplitudes, e.g. it was created using CTF correction (Wiener filtering) inside RELION or from a PDB. \n\n\
If set to No, then in the first iteration, the Fourier transforms of the reference projections are not multiplied by the CTFs.");
	}
	else
	{
		ctf_corrected_ref.clear();
	}

	AddBooleanEntry(ctf_phase_flipped, "Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");

	AddBooleanEntry(ctf_intact_first_peak, "Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	ctf_group->end();

	do_ctf_correction.cb_menu_i(); // To make default effective

	tab2->end();

	// tab3 holds the optimisation parameters
	tab3->begin();
	tab3->copy_label("Optimisation");

	resetHeight();

	if (run_type == RUN3D || run_type == RUN3DAUTO)
	{
		AddSliderEntry(ini_high, "Initial low-pass filter (A):", 60, 0, 200, 5, "It is recommended to strongly low-pass filter your initial reference map. \
If it has not yet been low-pass filtered, it may be done internally using this option. \
If set to 0, no low-pass filter will be applied to the initial reference(s).");
	}
	else
	{
		ini_high.clear();
	}

	if (run_type == RUN2D || run_type == RUN3D)
	{
		AddSliderEntry(nr_iter, "Number of iterations:", 25, 1, 50, 1, "Number of iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");

		AddSliderEntry(tau_fudge, "Regularisation parameter T:", 1 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. \
Too small values yield too-low resolution structures; too high values result in over-estimated resolutions, mostly notable by the apparition of high-frequency noise in the references.");

	}
	else
	{
		nr_iter.clear();
		tau_fudge.clear();
	}

	// Add a little spacer
	current_y += step_y/2;

	AddSliderEntry(particle_diameter, "Particle mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");

	if (run_type == RUN2D || run_type == RUN3D)
	{
		AddBooleanEntry(do_zero_mask, "Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. in 3D auto-refine) tend to work better when filling the solvent area with random noise, some classifications go better when using zeros.");

	}
	else
	{
		do_zero_mask.clear();
	}

	AddFileNameEntry(fn_mask, "Reference mask (optional):", "", "Image Files (*.{spi,vol,msk,mrc})", "\
If no mask is provided, a soft spherical mask based on the particle diameter will be used.\n\
\n\
Otherwise, provide a Spider/mrc map containing a (soft) mask with the same \
dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. \
The reconstructed reference map will be multiplied by this mask.\n\
\n\
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask \
the corresponding pixels in the reconstructed map are set to the average value of these pixels. \
Thereby, for example, the higher density inside the virion may be set to a constant. \
Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. \
To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Running tab).");

	if (run_type == RUN3DAUTO)
	{
		// Add a little spacer
		current_y += step_y/2;

		AddBooleanEntry(do_movies, "Realign movie frames?", false, "If set to Yes, then running averages of the individual frames \
of recorded movies will be aligned as independent particles.", movie_group);

		movie_group->begin();

		AddFileNameEntry(fn_movie_star, "Input movie frames:", "", "STAR Files (*.{star})", "Select the output STAR file from the preprocessing \
procedure of the movie frames.");

		AddSliderEntry(movie_runavg_window, "Running average window:", 5, 1, 15, 1, "The individual movie frames will be averaged using a running \
average window with the specified width. Use an odd number. The optimal value will depend on the SNR in the individual movie frames. For ribosomes, we used a value of 5, where \
each movie frame integrated approximately 1 electron per squared Angstrom.");

		AddSliderEntry(movie_sigma_angles, "Stddev on the rotations (deg):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 1 degree");

	//	AddRadioEntry(movie_sampling, "Angular sampling interval:", sampling_options, &sampling_options[5], "There are only a few discrete \
	angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
	The samplings are approximate numbers and vary slightly over the sphere.\n");

		AddSliderEntry(movie_sigma_offset, "Stddev on the translations (pix):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 2 pixels");

	//	AddSliderEntry(movie_offset_step, "Movie offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).");


		movie_group->end();
		do_movies.cb_menu_i(); // to make default effective
	}
	else
	{
		do_movies.clear();
		fn_movie_star.clear();
		movie_runavg_window.clear();
		movie_sigma_angles.clear();
		//movie_sampling.clear();
		movie_sigma_offset.clear();
		//movie_offset_step.clear();
	}

	tab3->end();

	// tab4 holds the sampling of the hidden variables (orientations and translations)
	tab4->begin();
	tab4->copy_label("Sampling");

	resetHeight();

	if (run_type == RUN3DAUTO)
	{
		 AddTextOnlyEntry(autosample_text1,"Note that initial sampling rates will be auto-incremented.");
	}
	else
	{
		autosample_text1.clear();
	}

	if (run_type == RUN2D)
	{
		AddSliderEntry(psi_sampling, "In-plane angular sampling:", 5., 0.5, 20, 0.5, "The sampling rate for the in-plane rotation angle (psi) in degrees. \
Using fine values will slow down the program. Recommended value for most 2D refinements: 5 degrees.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	}
	else
	{
		psi_sampling.clear();
	}

	if (run_type == RUN3D || run_type == RUN3DAUTO)
	{
		AddRadioEntry(sampling, "Angular sampling interval:", sampling_options, &sampling_options[2], "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	}
	else
	{
		sampling.clear();
	}

	AddSliderEntry(offset_range, "Offset search range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	AddSliderEntry(offset_step, "Offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	if (run_type == RUN3D)
	{

		// Add a little spacer
		current_y += step_y/2;

		AddBooleanEntry(do_local_ang_searches, "Perform local angular searches?", false, "If set to Yes, then rather than \
performing exhaustive angular searches, local searches within the range given below will be performed. \
A prior Gaussian distribution centered at the optimal orientation in the previous iteration and \
with a stddev of 1/3 of the range given below will be enforced.", localsearch_group);

		localsearch_group->begin();

		AddSliderEntry(sigma_angles, "Local angular search range:", 5., 0, 15, 0.1, "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.");

		localsearch_group->end();

		do_local_ang_searches.cb_menu_i(); // This is to make the default effective
	}
	else
	{
		do_local_ang_searches.clear();
		sigma_angles.clear();
	}

	if (run_type == RUN3DAUTO)
	{
		// Add a little spacer
		current_y += step_y/2;

		AddRadioEntry(auto_local_sampling, "Local searches from auto-sampling:", sampling_options, &sampling_options[4], "In the automated procedure to \
increase the angular samplings, local angular searches of -6/+6 times the sampling rate will be used from this angular sampling rate onwards.");
	}
	else
	{
		auto_local_sampling.clear();
	}

	tab4->end();

}

void RelionMainWindow::setup(int newruntype, bool _is_continue)
{

	//Resize to original size
	size(ori_w, ori_h);

	run_type = newruntype;
	is_continue = _is_continue;

    if (newruntype == RUN2D || newruntype == RUN3D || newruntype == RUN3DAUTO )
		setupRefine(newruntype, _is_continue);
	else if (newruntype == PREPROCESS)
		setupPreprocess();

    // tab5 holds computational stuff: valid for all runtypes
    tab5->begin();

    resetHeight();

	AddSliderEntry(nr_mpi, "Number of MPI procs:", 1, 1, 64, 1, "Number of MPI nodes to use in parallel. When set to 1, MPI will not be used.");

	if (newruntype == RUN2D || newruntype == RUN3D || newruntype == RUN3DAUTO )
	{
		AddSliderEntry(nr_threads, "Number of threads:", 1, 1, 16, 1, "Number of shared-memory (POSIX) threads to use in parallel. \
When set to 1, no multi-threading will be used. Multi-threading is often useful in 3D refinements to have more memory. 2D class averaging often proceeds more efficiently without threads.");
	}
	else
	{
		nr_threads.clear();
	}

	// Add a little spacer
    current_y += step_y/2;

    AddBooleanEntry(do_queue, "Submit to queue?", false, "Is set to Yes, the job will be submit to a queue, otherwise \
the job will be executed locally. Note that only MPI jobs may be sent to a queue.", queue_group);

	queue_group->begin();

	AddAnyEntry(queuename, "Queue name: ", "openmpi_8", "Name of the queue to which to submit the job.");

	AddAnyEntry(qsub, "Queue submit command:", "qsub", "Name of the command used to submit scripts to the queue, e.g. qsub or bsub.\n\n\
Note that the person who installed RELION should have made a custom script for your cluster/queue setup. Check this is the case \
(or create your own script following the RELION WIKI) if you have trouble submitting jobs.");

	// Two additional options that may be set through environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 (for more flexibility)
	char * extra1_text = getenv ("RELION_QSUB_EXTRA1");
	if (extra1_text != NULL)
	{
		have_extra1 = true;
		char * extra1_default = getenv ("RELION_QSUB_EXTRA1_DEFAULT");
		if (extra1_default == NULL)
			extra1_default = "";
		AddAnyEntry(qsub_extra1, extra1_text, extra1_default, "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra1XXX will be changed by this value.");
	}
	else
		have_extra1 = false;

	char * extra2_text = getenv ("RELION_QSUB_EXTRA2");
	if (extra2_text != NULL)
	{
		have_extra2 = true;
		char * extra2_default = getenv ("RELION_QSUB_EXTRA2_DEFAULT");
		if (extra2_default == NULL)
			extra2_default = "";
		AddAnyEntry(qsub_extra2, extra2_text, "", "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra2XXX will be changed by this value.");
	}
	else
		have_extra2 = false;


	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_QSUB_TEMPLATE");
	if (default_location==NULL)
		default_location=DEFAULTQSUBLOCATION;

	AddFileNameEntry(qsubscript, "Standard submission script:", default_location, "Script Files (*.{csh,sh,bash,script})",
"The template for your standard queue job submission script. \
Its default location may be changed by setting the environment variable RELION_QSUB_TEMPLATE. \
In the template script a number of variables will be replaced: \n \
XXXcommandXXX = relion command + arguments; \n \
XXXqueueXXX = The queue name; \n \
XXXmpinodesXXX = The number of MPI nodes; \n \
XXXthreadsXXX = The number of threads; \n \
XXXcoresXXX = The number of MPI nodes * nr_threads; \n \
If these options are not enough for your standard jobs, you may define two extra variables: XXXextra1XXX and XXXextra2XXX \
Their help text is set by the environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 \
For example, setenv RELION_QSUB_EXTRA1 \"Max number of hours in queue\" will result in an additional (text) ein the GUI \
Any variables XXXextra1XXX in the template script will be replaced by the corresponding value.\
Likewise, default values for the extra entries can be set through environment variables RELION_QSUB_EXTRA1_DEFAULT and  RELION_QSUB_EXTRA2_DEFAULT. \
But note that (unlike all other entries in the GUI) the extra values are not remembered from one run to the other.");


	queue_group->end();

	do_queue.cb_menu_i(); // This is to make the default effective

    // Add a little spacer
    current_y += step_y/2;

    AddAnyEntry(other_args, "Additional arguments:", "", "In this box command-line arguments may be provided that are not generated by the GUI. \
This may be useful for testing developmental options and/or expert use of the program. \
The command 'relion_refine' will print a list of possible options.");

	// Now add the final functionality to send off the job
	AddRunButtons();

	tab5->end();

    // Redraw the entire GUI
	if (run_type != PREPROCESS)
		toggle_new_continue();
    redraw();

}

void RelionMainWindow::toggle_new_continue()
{

    if (run_type != PREPROCESS)
    {
		fn_cont.deactivate(!is_continue);
		fn_img.deactivate(is_continue);
		angpix.deactivate(is_continue);
		if (run_type == RUN3D || run_type == RUN3DAUTO )
		{
			fn_ref.deactivate(is_continue);
			ref_correct_greyscale.deactivate(is_continue);
			ini_high.deactivate(is_continue);
		}
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			sym_group.deactivate(is_continue);
			sym_nr.deactivate(is_continue);
		}
		if (run_type == RUN2D || run_type == RUN3D)
		{
			nr_classes.deactivate(is_continue);
			tau_fudge.deactivate(is_continue);
		}

		particle_diameter.deactivate(is_continue);
		if (run_type == RUN2D || run_type == RUN3D)
		{
			do_zero_mask.deactivate(is_continue);
		}
		fn_mask.deactivate(is_continue);

		if (run_type == RUN3DAUTO)
		{
			do_movies.deactivate(!is_continue);
			fn_movie_star.deactivate(!is_continue);
			movie_runavg_window.deactivate(!is_continue);
			movie_sigma_angles.deactivate(!is_continue);
			//movie_sampling.deactivate(!is_continue);
			movie_sigma_offset.deactivate(!is_continue);
			//movie_offset_step.deactivate(!is_continue);
		}

		if (run_type == RUN3DAUTO)
		{
			sampling.deactivate(is_continue);
			offset_range.deactivate(is_continue);
			offset_step.deactivate(is_continue);
		}

		do_ctf_correction.deactivate(is_continue);
		ctf_phase_flipped.deactivate(is_continue);
		ctf_intact_first_peak.deactivate(is_continue);
		if (run_type == RUN3D || run_type == RUN3DAUTO )
		{
			ctf_corrected_ref.deactivate(is_continue);
		}
    }
}

void RelionMainWindow::resetHeight()
{
	current_y = start_y;
}

void RelionMainWindow::AddAnyEntry(AnyEntry &var,
		const char * title,
		const char* defaultvalue,
		const char* helptext)
{

    // Clear if existing
	var.clear();

	// Add the entry to the window
	var.initialise(XCOL1, current_y, step_y, WCOL2, WCOL3, title, defaultvalue, helptext);

    // Update the Y-coordinate
    current_y += step_y + 2;
}

void RelionMainWindow::AddTextOnlyEntry(textOnlyEntry &var,
		const char * text)
{

    // Clear if existing
	var.clear();

	// Add the entry to the window
	// Add 3 to step_y, otherwise the text may not fit...
	var.initialise(XCOL1, current_y, WCOL1 + WCOL2 + WCOL3, step_y + 6, text);

    // Update the Y-coordinate
    current_y += step_y + 6 + 2;
}

void RelionMainWindow::AddFileNameEntry(FileNameEntry &var,
		const char * title,
		const char* defaultvalue,
		const char* pattern,
		const char* helptext)
{

    // Clear if existing
	var.clear();

	// Add the entry to the window
	var.initialise(XCOL1, current_y, step_y,  WCOL2, WCOL3, WCOL4, title, defaultvalue, pattern, helptext);

    // Update the Y-coordinate
    current_y += step_y + 2;
}

void RelionMainWindow::AddRadioEntry(RadioEntry &var,
		const char * title,
		Fl_Menu_Item *options,
		Fl_Menu_Item *defaultvalue,
		const char* helptext)
{

    // Clear if existing
	var.clear();

	// Add the entry to the window
	var.initialise(XCOL1, current_y, step_y, WCOL2, WCOL3, WCOL4, title, options, defaultvalue, helptext);

    // Update the Y-coordinate
    current_y += step_y + 2;
}

void RelionMainWindow::AddBooleanEntry(BooleanEntry &var,
		const char * title,
		bool defaultvalue,
		const char* helptext,
		Fl_Group* deactivate_this_group)
{

    // Clear if existing
	var.clear();

	// Add the entry to the window
	var.initialise(XCOL1, current_y, step_y, WCOL2, WCOL3, WCOL4, title, defaultvalue, helptext, deactivate_this_group);

    // Update the Y-coordinate
    current_y += step_y + 2;
}

void RelionMainWindow::AddSliderEntry(SliderEntry &var,
		            const char * title,
		            float defaultvalue,
		            float minvalue,
		            float maxvalue,
		            float valuestep,
					const char* helptext)
{
    // Clear if existing
	var.clear();

	// Add the entry to the window
	var.initialise(XCOL1, current_y, step_y, WCOL2, WCOL3, WCOL4, title, defaultvalue, minvalue, maxvalue, valuestep, helptext);

    // Update the Y-coordinate
    current_y += step_y + 2;

}

void RelionMainWindow::setOutputName()
{
	if (run_type == PREPROCESS)
		output_name = fn_out_preprocess.getValue();
	else
		output_name = fn_out.getValue();

	if (run_type != PREPROCESS && is_continue)
    {
    	// Get the continuation iteration number X of this file, and add "_ctX" to the output_name
    	int pos_it = fn_cont.getValue().rfind("_it");
    	int pos_op = fn_cont.getValue().rfind("_optimiser");
    	if (pos_it < 0 || pos_op < 0)
    		std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << fn_cont.getValue() << std::endl;
    	int it = (int)textToFloat((fn_cont.getValue().substr(pos_it+3, 6)).c_str());
    	output_name += "_ct" + floatToString(it);
    }
}

std::string RelionMainWindow::getCommandLine()
{
	// Always use oversampling of 1...
	int iover = 1;

	std::string cline;
	// Choose program
    if (run_type == PREPROCESS)
    {
    	if (nr_mpi.getValue() > 1)
			cline="`which relion_preprocess_mpi`";
		else
			cline="`which relion_preprocess`";


		cline += " --o " + fn_out_preprocess.getValue();

		// CTFFIND stuff
		if (do_ctffind.getValue())
		{
			cline += " --ctffind \"" + mic_names.getValue()+"\"";
                        cline += " --ctfWin " + floatToString(ctf_win.getValue());
			cline += " --CS " + floatToString(cs.getValue());
			cline += " --HT " + floatToString(kv.getValue());
			cline += " --AmpCnst " + floatToString(q0.getValue());
			cline += " --XMAG " + floatToString(mag.getValue());
			cline += " --DStep " + floatToString(dstep.getValue());
			cline += " --Box " + floatToString(box.getValue());
			cline += " --ResMin " + floatToString(resmin.getValue());
			cline += " --ResMax " + floatToString(resmax.getValue());
			cline += " --dFMin " + floatToString(dfmin.getValue());
			cline += " --dFMax " + floatToString(dfmax.getValue());
			cline += " --FStep " + floatToString(dfstep.getValue());
			cline += " --ctffind3_exe " + fn_ctffind3_exe.getValue();
		}

		// Extraction stuff
		if (do_starfile.getValue())
		{
			cline += " --coord_files \"" + coord_names.getValue()+"\"";

			if (do_extract.getValue())
			{
				cline += " --extract";
				cline += " --coord_format " + coord_format.getValue();
				cline += " --extract_size " + floatToString(extract_size.getValue());
				if (do_movie_extract.getValue())
				{
					cline += " --extract_movies";
					cline += " --movie_rootname " + movie_rootname.getValue();
					cline += " --avg_movie_frames " + floatToString(avg_movie_frames.getValue());
					cline += " --first_movie_frame " + floatToString(first_movie_frame.getValue());
					cline += " --last_movie_frame " + floatToString(last_movie_frame.getValue());
				}

				// Operate stuff
				if (do_rescale.getValue())
					cline += " --scale " + floatToString(rescale.getValue());
				if (do_norm.getValue())
				{
					cline += " --norm --bg_radius " + floatToString(bg_radius.getValue());
					cline += " --white_dust " + floatToString(white_dust.getValue());
					cline += " --black_dust " + floatToString(black_dust.getValue());
				}
				if (do_invert.getValue())
					cline += " --invert_contrast ";
			}

		}

    }
    else
    {
		if (nr_mpi.getValue() > 1)
			cline="`which relion_refine_mpi`";
		else
			cline="`which relion_refine`";

		// I/O
		// Save the real output name (could be with _ctX for continuation)
		// This name will also be used for the stderr and stdout outputs and the submit script and gui settings filenames
		cline += " --o " + output_name;
		if (is_continue)
		{
			cline += " --continue " + fn_cont.getValue();
		}
		else
		{
			cline += " --i " + fn_img.getValue();
			cline += " --particle_diameter " + floatToString(particle_diameter.getValue());
			cline += " --angpix " + floatToString(angpix.getValue());
			if (run_type == RUN3D || run_type == RUN3DAUTO )
			{
				cline += " --ref " + fn_ref.getValue();
				if (!ref_correct_greyscale.getValue())
				{
					cline += " --firstiter_cc";
				}
				if (ini_high.getValue() > 0.)
				{
					cline += " --ini_high " + floatToString(ini_high.getValue());
				}
			}
		}

		// Optimisation
		if (run_type == RUN2D || run_type == RUN3D)
		{
			cline += " --iter " + floatToString(nr_iter.getValue());
			cline += " --tau2_fudge " + floatToString(tau_fudge.getValue());
		}

		// Movies
		if (run_type == RUN3DAUTO && is_continue && do_movies.getValue())
		{
			cline += " --realign_movie_frames " + fn_movie_star.getValue();
			cline += " --movie_frames_running_avg " + floatToString(movie_runavg_window.getValue());
			cline += " --sigma_ang " + floatToString(movie_sigma_angles.getValue());
			cline += " --sigma_off " + floatToString(movie_sigma_offset.getValue());
		}

		// Always flatten the solvent
		cline += " --flatten_solvent";
		if (run_type == RUN2D || run_type == RUN3D)
		{
			if (do_zero_mask.getValue())
				cline += " --zero_mask";
		}
		if (fn_mask.getValue().length() > 0)
			cline += " --solvent_mask " + fn_mask.getValue();

		if (!is_continue)
		{
			// CTF stuff
			if (do_ctf_correction.getValue())
			{
				cline += " --ctf";
				if (run_type == RUN3D || run_type == RUN3DAUTO)
				{
					if (ctf_corrected_ref.getValue())
						cline += " --ctf_corrected_ref";
				}
				//if (only_flip_phases.getValue())
				//	cline += " --only_flip_phases";
				if (ctf_phase_flipped.getValue())
					cline += " --ctf_phase_flipped";
				if (ctf_intact_first_peak.getValue())
					cline += " --ctf_intact_first_peak";
			}

			// Sampling stuff
			if (run_type == RUN3D || run_type == RUN3DAUTO)
			{
				if (strcmp((sym_group.getValue()).c_str(), "O") == 0 ||
					strcmp((sym_group.getValue()).c_str(), "T") == 0 )
					cline += " --sym " + sym_group.getValue();
				else
					cline += " --sym " + sym_group.getValue() + floatToString(sym_nr.getValue());
			}

			// Only multiple classes for classification runs (by default --K 1)
			if (run_type == RUN2D || run_type == RUN3D)
			{
				cline += " --K " + floatToString(nr_classes.getValue());
			}

		}

		// Always use oversampling of 1...
		cline += " --oversampling " + floatToString((float)iover);

		// If 3D refinement: use autosampling
		if (run_type == RUN3DAUTO)
		{
			if (!is_continue)
			{
				cline += " --auto_refine --split_random_halves";

				// For C-point group symmetries, join half-reconstructions up to 40A to prevent diverging orientations
				if (sym_group.getValue() == "C")
					cline += " --low_resol_join_halves 40";

				// Manually provide initial sampling
				for (int i = 0; i < 10; i++)
				{
					if (strcmp((sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
					{
						// The sampling given in the GUI will be the oversampled one!
						cline += " --healpix_order " + floatToString((float)i + 1 - iover);
						break;
					}
				}

				// Offset range
				cline += " --offset_range " + floatToString(offset_range.getValue());

				// The sampling given in the GUI will be the oversampled one!
				cline += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));
			}

			// Minimum sampling rate to perform local searches (may be changed upon continuation
			for (int i = 0; i < 10; i++)
			{
				if (strcmp((auto_local_sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
				{
					cline += " --auto_local_healpix_order " + floatToString((float)i + 1 - iover);
					break;
				}
			}

		}
		else // conventional runs: manual sampling
		{
			if (run_type == RUN3D)
			{
				for (int i = 0; i < 10; i++)
				{
					if (strcmp((sampling.getValue()).c_str(), sampling_options[i].label()) == 0)
					{
						// The sampling given in the GUI will be the oversampled one!
						cline += " --healpix_order " + floatToString((float)i + 1 - iover);
						break;
					}
				}
				// Manually input local angular searches
				if (do_local_ang_searches.getValue())
					cline += " --sigma_ang " + floatToString(sigma_angles.getValue() / 3.);

			}
			else if (run_type == RUN2D)
			{
				// The sampling given in the GUI will be the oversampled one!
				cline += " --psi_step " + floatToString(psi_sampling.getValue() * pow(2., iover));
			}
			else
			{
				std::cerr << "this should NOT happen: unrecognised run_type..."<<std::endl;
				exit(1);
			}

			// Offset range
			cline += " --offset_range " + floatToString(offset_range.getValue());

			// The sampling given in the GUI will be the oversampled one!
			cline += " --offset_step " + floatToString(offset_step.getValue() * pow(2., iover));

		}

		// Always do norm and scale correction
		cline += " --norm --scale ";

		cline += " " + other_args.getValue();

		// Running stuff
		cline += " --j " + floatToString(nr_threads.getValue());
    } // end run_type is not PREPROCESS i.e. is REFINE

    return cline;

}

void RelionMainWindow::saveJobSubmissionScript(std::string newfilename)
{
	Fl_Text_Buffer *textbuf = new Fl_Text_Buffer;

	// Open the standard job submission file
	int errno;
	if (errno = textbuf->loadfile(qsubscript.getValue().c_str()))
	    fl_alert("Error reading from file \'%s\':\n%s.", qsubscript.getValue().c_str(), strerror(errno));

	int nthr = 1; // default to a single thread
	if (run_type != PREPROCESS)
		nthr = nr_threads.getValue();

	replaceString(textbuf, "XXXmpinodesXXX", floatToString(nr_mpi.getValue()) );
	replaceString(textbuf, "XXXthreadsXXX", floatToString(nthr) );
	replaceString(textbuf, "XXXcoresXXX", floatToString(nr_mpi.getValue() * nthr) );
	replaceString(textbuf, "XXXnameXXX", output_name);
	replaceString(textbuf, "XXXerrfileXXX", output_name + ".err");
	replaceString(textbuf, "XXXoutfileXXX", output_name + ".out");
	replaceString(textbuf, "XXXqueueXXX", queuename.getValue() );
	if (have_extra1)
		replaceString(textbuf, "XXXextra1XXX", qsub_extra1.getValue() );
	if (have_extra2)
		replaceString(textbuf, "XXXextra2XXX", qsub_extra2.getValue() );
	replaceString(textbuf, "XXXcommandXXX", getCommandLine() );

	// Save the modified job submission script using a local name
	if (errno = textbuf->savefile(newfilename.c_str()))
	    fl_alert("Error writing to file \'%s\':\n%s.", newfilename.c_str(), strerror(errno));

}

void RelionMainWindow::writeSettings(std::string filename)
{
    std::ofstream  fh;
    fh.open((filename).c_str(), std::ios::out);
    if (!fh)
    {
    	std::cerr << "Cannot write to file: "<<filename<<std::endl;
    	exit(1);
    }

    // is_continue flag
    if (is_continue)
    	fh << "is_continue == true" << std::endl;
    else
    	fh << "is_continue == false" << std::endl;

    if (run_type == PREPROCESS)
    {
    	fn_out_preprocess.writeValue(fh);
    	cs.writeValue(fh);
    	kv.writeValue(fh);
    	q0.writeValue(fh);
    	mag.writeValue(fh);
    	dstep.writeValue(fh);

    	do_ctffind.writeValue(fh);
    	mic_names.writeValue(fh);
        ctf_win.writeValue(fh);
    	box.writeValue(fh);
    	resmin.writeValue(fh);
    	resmax.writeValue(fh);
    	dfmin.writeValue(fh);
    	dfmax.writeValue(fh);
    	dfstep.writeValue(fh);
    	fn_ctffind3_exe.writeValue(fh);

    	do_extract.writeValue(fh);
    	coord_names.writeValue(fh);
    	coord_format.writeValue(fh);
    	extract_size.writeValue(fh);
    	do_movie_extract.writeValue(fh);
    	movie_rootname.writeValue(fh);
    	avg_movie_frames.writeValue(fh);
    	first_movie_frame.writeValue(fh);
    	last_movie_frame.writeValue(fh);

    	do_rescale.writeValue(fh);
    	rescale.writeValue(fh);
    	do_norm.writeValue(fh);
    	bg_radius.writeValue(fh);
    	white_dust.writeValue(fh);
    	black_dust.writeValue(fh);
    	do_invert.writeValue(fh);
    	do_starfile.writeValue(fh);

    }
    else if (run_type == RUN2D || run_type == RUN3D || run_type == RUN3DAUTO)
    {
    	//I/O
		fn_out.writeValue(fh);
		fn_cont.writeValue(fh);
		fn_img.writeValue(fh);
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			fn_ref.writeValue(fh);
			ref_correct_greyscale.writeValue(fh);
			ini_high.writeValue(fh);
		}
		particle_diameter.writeValue(fh);
		angpix.writeValue(fh);

		//Optimisation
		if (run_type == RUN2D || run_type == RUN3D)
		{
			nr_iter.writeValue(fh);
			tau_fudge.writeValue(fh);
			do_zero_mask.writeValue(fh);
		}
		fn_mask.writeValue(fh);

		if (run_type == RUN3DAUTO)
		{
			do_movies.writeValue(fh);
			fn_movie_star.writeValue(fh);
			movie_runavg_window.writeValue(fh);
			movie_sigma_angles.writeValue(fh);
			//movie_sampling.writeValue(fh);
			movie_sigma_offset.writeValue(fh);
			//movie_offset_step.writeValue(fh);
		}


		//do_magn_correction.writeValue(fh);
		do_ctf_correction.writeValue(fh);
		//only_flip_phases.writeValue(fh);
		ctf_phase_flipped.writeValue(fh);
		ctf_intact_first_peak.writeValue(fh);
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			ctf_corrected_ref.writeValue(fh);
		}

		//Sampling
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			sym_group.writeValue(fh);
			sym_nr.writeValue(fh);
			sampling.writeValue(fh);
			do_local_ang_searches.writeValue(fh);
			sigma_angles.writeValue(fh);
		}
		if (run_type == RUN2D)
		{
			psi_sampling.writeValue(fh);
		}
		if (run_type == RUN2D || run_type == RUN3D)
		{
			nr_classes.writeValue(fh);
		}
		offset_range.writeValue(fh);
		offset_step.writeValue(fh);
		if (run_type == RUN3DAUTO)
		{
			auto_local_sampling.writeValue(fh);
		}
    }

    //Running
    nr_mpi.writeValue(fh);
    if (run_type != PREPROCESS)
    	nr_threads.writeValue(fh);
    do_queue.writeValue(fh);
    queuename.writeValue(fh);
    qsub.writeValue(fh);
    qsubscript.writeValue(fh);
    other_args.writeValue(fh);


    // Write the actual file to disc
    fh.close();
}

bool RelionMainWindow::readSettings(std::string filename, bool do_check)
{

	std::ifstream in(filename.c_str(), std::ios_base::in);
    if (in.fail())
    {
    	if (!do_check)
    	{
    		return false;
    	}
    	else
    	{
    		std::cerr << "Error reading "<<filename<<std::endl;
    		exit(1);
    	}
    }

	in.seekg(0, std::ios::beg);
	std::string line;
	getline(in, line, '\n');
	if (line.rfind("is_continue == true") == 0)
	{

		// If current GUI was not a continue run, then toggle the continue settings
		if (!is_continue)
			toggle_new_continue();
		is_continue = true;
		// Also set the toggle menu
		choice_continue->picked(&continue_options[1]);
	}
	else
	{
		// If current GUI was a continue run, then toggle the continue settings
		if (is_continue)
			toggle_new_continue();
		is_continue = false;
		// Also set the toggle menu
		choice_continue->picked(&continue_options[0]);
	}

	// Check the runtype and switch to the appropriate window if it changes
    if (filename.find("_gui3dauto.settings") < filename.length())
    {
    	choice_runtype->picked(&runtype_options[3]);
    	setup(RUN3DAUTO, is_continue);
    }
    else if (filename.find("_gui3d.settings") < filename.length())
    {
    	choice_runtype->picked(&runtype_options[2]);
    	setup(RUN3D, is_continue);
    }
    else if (filename.find("_gui2d.settings") < filename.length())
    {
    	choice_runtype->picked(&runtype_options[1]);
    	setup(RUN2D, is_continue);
    }
    else if (filename.find("_preprocess.settings") < filename.length())
    {
    	choice_runtype->picked(&runtype_options[0]);
    	setup(PREPROCESS, is_continue);
    }
    else
    {
    	std::cerr << "Unrecognised gui settings file (does not end in _gui2d.settings, _gui3d.settings or _gui3dauto.settings)" << std::endl;
    	exit(1);
    }

    if (run_type == PREPROCESS)
    {
    	fn_out_preprocess.readValue(in);
    	cs.readValue(in);
    	kv.readValue(in);
    	q0.readValue(in);
    	mag.readValue(in);
    	dstep.readValue(in);

    	do_ctffind.readValue(in);
    	mic_names.readValue(in);

        ctf_win.readValue(in);
    	resmin.readValue(in);
    	resmax.readValue(in);
    	dfmin.readValue(in);
    	dfmax.readValue(in);
    	dfstep.readValue(in);
    	fn_ctffind3_exe.readValue(in);

    	do_extract.readValue(in);
    	coord_names.readValue(in);
    	coord_format.readValue(in);
    	extract_size.readValue(in);
    	do_movie_extract.readValue(in);
    	movie_rootname.readValue(in);
    	avg_movie_frames.readValue(in);
    	first_movie_frame.readValue(in);
    	last_movie_frame.readValue(in);

    	do_rescale.readValue(in);
    	rescale.readValue(in);
    	do_norm.readValue(in);
    	bg_radius.readValue(in);
    	white_dust.readValue(in);
    	black_dust.readValue(in);
    	do_invert.readValue(in);
    	do_starfile.readValue(in);

    }
    else if (run_type == RUN2D || run_type == RUN3D || run_type == RUN3DAUTO)
    {
    	//I/O
		fn_out.readValue(in);
		fn_cont.readValue(in);
		fn_img.readValue(in);
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			fn_ref.readValue(in);
			ref_correct_greyscale.readValue(in);
			ini_high.readValue(in);
		}
		particle_diameter.readValue(in);
		angpix.readValue(in);

		//Optimisation
		if (run_type == RUN2D || run_type == RUN3D)
		{
			nr_iter.readValue(in);
			tau_fudge.readValue(in);
			do_zero_mask.readValue(in);
		}
		fn_mask.readValue(in);

		if (run_type == RUN3DAUTO)
		{
			do_movies.readValue(in);
			fn_movie_star.readValue(in);
			movie_runavg_window.readValue(in);
			movie_sigma_angles.readValue(in);
			//movie_sampling.readValue(in);
			movie_sigma_offset.readValue(in);
			//movie_offset_step.readValue(in);
		}

		//do_magn_correction.readValue(in);
		do_ctf_correction.readValue(in);
		//only_flip_phases.readValue(in);
		ctf_phase_flipped.readValue(in);
		ctf_intact_first_peak.readValue(in);
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			ctf_corrected_ref.readValue(in);
		}

		//Sampling
		if (run_type == RUN3D || run_type == RUN3DAUTO)
		{
			sym_group.readValue(in);
			sym_nr.readValue(in);
			sampling.readValue(in);
			do_local_ang_searches.readValue(in);
			sigma_angles.readValue(in);
		}
		if (run_type == RUN2D)
		{
			psi_sampling.readValue(in);
		}
		if (run_type == RUN2D || run_type == RUN3D)
		{
			nr_classes.readValue(in);
		}
		offset_range.readValue(in);
		offset_step.readValue(in);
		if (run_type == RUN3DAUTO)
		{
			auto_local_sampling.readValue(in);
		}
    }

    //Running
    nr_mpi.readValue(in);
    if (run_type != PREPROCESS)
    	nr_threads.readValue(in);
    do_queue.readValue(in);
    queuename.readValue(in);
    qsub.readValue(in);
    qsubscript.readValue(in);
    other_args.readValue(in);

    // Redraw the entire GUI
    redraw();

    return true;
}

void RelionMainWindow::AddRunButtons()
{

	// Delete print_CL button if it already exists
	if (print_CL_button != NULL)
		delete print_CL_button;

	print_CL_button = new Fl_Button(XCOL3-270, current_y + step_y, 100, 30, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	//print_CL_button->labelfont(FL_ITALIC);
	print_CL_button->labelsize(12);
	//run->labeltype(FL_SHADOW_LABEL);
	print_CL_button->callback( cb_print_cl, this);

	// Delete cite button if it already exists
	if (cite_button != NULL)
		delete cite_button;
	cite_button = new Fl_Button(XCOL3-150, current_y + step_y, 100, 30, "What to cite?");
	cite_button->color(GUI_RUNBUTTON_COLOR);
	//cite_button->labelfont(FL_ITALIC);
	cite_button->labelsize(12);
	//cite->labeltype(FL_SHADOW_LABEL);
	cite_button->callback( cb_cite, this);

	// Delete Run button if it already exists
	if (run_button != NULL)
		delete run_button;

	run_button = new Fl_Button(XCOL3-30, current_y + step_y, 100, 30, "Run!");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(18);
	//run_button->labeltype(FL_SHADOW_LABEL);
	run_button->callback( cb_run, this);




}

// Help button call-back functions
void RelionMainWindow::cb_run(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_run_i();
}

// Help button call-back functions
void RelionMainWindow::cb_print_cl(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_print_cl_i();
}

void RelionMainWindow::cb_cite(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_cite_i();
}


void RelionMainWindow::cb_run_i()
{
	// Save the parameter settings in a file
	cb_menubar_save_i();
	std::string command;
	if (!do_queue.getValue())
	{
		// Do not submit to queue, i.e. execute locally
		if (nr_mpi.getValue() > 1)
			command = "mpirun -n " + floatToString(nr_mpi.getValue()) + " " + getCommandLine() + " &";
		else
			command = getCommandLine() + " &" ;
	}
	else
	{
		// Submit to queue (here use output_name again!)
		std::string output_script = output_name + "_submit.script";
		saveJobSubmissionScript(output_script);
		command = qsub.getValue() + " " + output_script + " &";
	}

	std::cerr << "Executing: " << command << std::endl;
	system(command.c_str());

	// Deactivate Run button to prevent the user from accidentally submitting many jobs
	run_button->deactivate();

}
void RelionMainWindow::cb_print_cl_i()
{
	// Set the global variable output_name
	setOutputName();
	std::cout << " *** The command is:" << std::endl;
	std::cout << getCommandLine() << std::endl;
}

void RelionMainWindow::cb_cite_i()
{
	ShowHelpText *help = new ShowHelpText("\
If RELION is useful in your work, please cite us in the following contexts:\n \
\n \
 * The general Bayesian approach (and the first mention of RELION): \n \
   - Scheres (2012) J. Mol. Biol. (DOI: 10.1016/j.jmb.2011.11.010)	 \n \
\n \
 * RELION implementation details and the 3D auto-refine procedure: \n \
   - Scheres (2012) J. Struct. Biol. (DOI: 10.1016/j.jsb.2012.09.006)	 \n \
\n \
 * The gold-standard FSC and the relevance of the 0.143 criterion: \n \
   - Scheres & Chen (2012) Nat. Meth. (DOI: 10.1038/nmeth.2115)	 \n \
\n \
 * The movie-processing procedure: \n \
   - Bai et al. (2013) eLife (in press)	 \n \
");
}

void RelionMainWindow::cb_menu_continue(Fl_Widget*, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menu_continue_i();
}

void RelionMainWindow::cb_menu_continue_i()
{
	const Fl_Menu_Item* m = menu_continue->mvalue();
	if (strcmp(m->label(), continue_options[0].label()) == 0)
		is_continue = false;
	else
		is_continue = true;
	toggle_new_continue();
}

void RelionMainWindow::cb_menu_runtype(Fl_Widget*, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menu_runtype_i();
}

void RelionMainWindow::cb_menu_runtype_i()
{
	const Fl_Menu_Item* m = menu_runtype->mvalue();
	if (strcmp(m->label(), runtype_options[0].label()) == 0 && run_type != PREPROCESS)
	{
	    // See if there is a last run gui settings file
	    // If so, fill the gui with those values
	    setup(PREPROCESS, is_continue);
		readSettings(LASTRUNSETTINGSPRE, false);
	}
	else if (strcmp(m->label(), runtype_options[1].label()) == 0 && run_type != RUN2D)
	{
	    // See if there is a last run gui settings file
	    // If so, fill the gui with those values
	    setup(RUN2D, is_continue);
		readSettings(LASTRUNSETTINGS2D, false);
	}
	else if (strcmp(m->label(), runtype_options[2].label()) == 0 && run_type != RUN3D)
	{
		setup(RUN3D, is_continue);
		readSettings(LASTRUNSETTINGS3D, false);
	}
	else if (strcmp(m->label(), runtype_options[3].label()) == 0 && run_type != RUN3DAUTO)
	{
		setup(RUN3DAUTO, is_continue);
		readSettings(LASTRUNSETTINGS3DAUTO, false);
	}

}

// call-back functions for the menubar
void RelionMainWindow::cb_menubar_load(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_load_i();
}

void RelionMainWindow::cb_menubar_load_i()
{
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", "*.settings", Fl_File_Chooser::SINGLE, "Choose a GUI settings file");

    G_chooser->directory(NULL);
    G_chooser->show();

    // Block until user picks something.
    //     (The other way to do this is to use a callback())
    //
    while(G_chooser->shown()) {
        Fl::wait();
    }

    // Print the results
    if ( G_chooser->value() == NULL ) {
        //fprintf(stderr, "(User hit 'Cancel')\n");
        return;
    }
    readSettings(G_chooser->value());

}

// Save button call-back function
void RelionMainWindow::cb_menubar_save(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_save_i();
}

void RelionMainWindow::cb_menubar_save_i()
{

	// Set the global variable output_name
	setOutputName();

	// Create output directory if the outname contains a "/"
	int last_slash = output_name.rfind("/");
	if (last_slash < output_name.size())
	{
		std::string dirs = output_name.substr(0, last_slash);
		std::string makedirs = "mkdir -p " + dirs;
		system(makedirs.c_str());
	}

	// Write settings file in the output directory
	if (run_type == RUN3DAUTO)
	{
		writeSettings(output_name+"_gui3dauto.settings");
		writeSettings(LASTRUNSETTINGS3DAUTO);
	}
	else if (run_type == RUN3D)
	{
		writeSettings(output_name+"_gui3d.settings");
		writeSettings(LASTRUNSETTINGS3D);
	}
	else if (run_type == RUN2D)
	{
		writeSettings(output_name+"_gui2d.settings");
		writeSettings(LASTRUNSETTINGS2D);
	}
	else if (run_type == PREPROCESS)
	{
		writeSettings(fn_out_preprocess.getValue()+"_preprocess.settings");
		writeSettings(LASTRUNSETTINGSPRE);
	}
	else
	{
		std::cerr << "Error: unrecognised run_type..." << std::endl;
		exit(1);
	}

}

void RelionMainWindow::cb_menubar_quit(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_quit_i();
}

void RelionMainWindow::cb_menubar_quit_i()
{
	exit(0);
}

void RelionMainWindow::cb_menubar_reactivate_runbutton(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_reactivate_runbutton_i();
}

void RelionMainWindow::cb_menubar_reactivate_runbutton_i()
{
	run_button->activate();
}


void RelionMainWindow::cb_menubar_about(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_about_i();
}

void RelionMainWindow::cb_menubar_about_i()
{
	ShowHelpText *help = new ShowHelpText("\
RELION is written by Sjors Scheres at the MRC Laboratory of Molecular Biology (scheres@mrc-lmb.cam.ac.uk).\n \
\n\
If RELION is useful in your work, please cite us in the contexts as explained under the \"What to cite?\" button on the Running tab. \n  \
\n\
Note that RELION is completely free, open-source software. You can redistribute it and/or modify it for your own purposes, but please do make sure \
the contribution of Sjors Scheres is acknowledged appropriately. In order to maintain an overview of existing versions, he would also appreciate being \
notified of any redistribution of (modified versions of) the code. \n \
");
}

void replaceString(Fl_Text_Buffer *textbuf, std::string findthis, std::string replaceby)
{
	const char *find = findthis.c_str();
	const char *replace = replaceby.c_str();

	// Loop through the whole string
	int pos = 0;
	for (int found = 1; found;) {
	    found = textbuf->search_forward(pos, find, &pos);

	    if (found) {
	      // Found a match; update the position and replace text...
	      textbuf->select(pos, pos+strlen(find));
	      textbuf->remove_selection();
	      textbuf->insert(pos, replace);
	      pos += strlen(replace);
	    }
	}
}


