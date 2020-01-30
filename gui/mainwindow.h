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

#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Menu_Bar.H>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define LASTRUNSETTINGS3D ".lastrun_gui3d.settings"
#define LASTRUNSETTINGS3DAUTO ".lastrun_gui3dauto.settings"
#define LASTRUNSETTINGS2D ".lastrun_gui2d.settings"
#define LASTRUNSETTINGSPRE ".lastrun_preprocess.settings"
#define XCOL1 10
#define XCOL2 260
#define XCOL3 460
#define XCOL4 475
#define XCOL5 535
#define COLUMN_SEPARATION 3
#define WCOL1 ( (XCOL2) - (XCOL1) - (COLUMN_SEPARATION) )
#define WCOL2 ( (XCOL3) - (XCOL2) - (COLUMN_SEPARATION) )
#define WCOL3 ( (XCOL4) - (XCOL3) - (COLUMN_SEPARATION) )
#define WCOL4 ( (XCOL5) - (XCOL4) - (COLUMN_SEPARATION) )
//version-1.0 #define GUI_BUTTON_COLOR (fl_rgb_color(200,255,100))
//version-1.1 #define GUI_BUTTON_COLOR (fl_rgb_color(50,150,250))
//devel-version
#define GUI_BUTTON_COLOR (fl_rgb_color(155,150,255))
//version-1.0 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,155,0))
//version-1.1 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,50,50))
//devel-version
#define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,53,100))
#define GUI_BACKGROUND_COLOR (fl_rgb_color(240,240,240))
#define GUI_BACKGROUND_COLOR2 (fl_rgb_color(200,200,200))
#define GUI_INPUT_COLOR (fl_rgb_color(255,255,230))
#define DEFAULTQSUBLOCATION "/lmb/home/scheres/app/relion/gui/qsub.csh"
#define DEFAULTCTFFINDLOCATION "/public/EM/CTFFIND/ctffind3.exe"
//#define GUI_INPUT_COLOR (fl_rgb_color(255,255,255))
// After defining all this include the entries.h
#include "entries.h"

// Type of window
#define RUN3DAUTO  4
#define RUN3D  3
#define RUN2D  2
#define PREPROCESS 1

Fl_Menu_Item symgroup_options[] = {
		      {"C"},
		      {"D"},
		      {"T"},
		      {"O"},
		      {"I"},
		      {0} // this should be the last entry
		    };
Fl_Menu_Item sampling_options[] = {
		      {"30 degrees"},
		      {"15 degrees"},
		      {"7.5 degrees"},
		      {"3.7 degrees"},
		      {"1.8 degrees"},
		      {"0.9 degrees"},
		      {"0.5 degrees"},
		      {"0.2 degrees"},
		      {"0.1 degrees"},
		      {0} // this should be the last entry
		    };
Fl_Menu_Item coord_format_options[] = {
			  {"boxer"},
		      {"ximdisp"},
		      {"xmipp2"},
		      {0} // this should be the last entry
		    };
Fl_Menu_Item continue_options[] = {
		      {"Start new run", FL_ALT + 'n'},
		      {"Continue old run", FL_ALT + 'c'},
		      {0} // this should be the last entry
		    };
Fl_Menu_Item runtype_options[] = {
			  {"Preprocessing", FL_ALT + '1'},
		      {"2D class averages", FL_ALT + '2'},
		      {"3D classification", FL_ALT + '3'},
		      {"3D auto-refine", FL_ALT + '4'},
		      {0} // this should be the last entry
		    };



	// This class organises the main winfow of the relion GUI
class RelionMainWindow : public Fl_Window{

public:

	// All my variables

	// Preprocess stuff
	// I/O
	AnyEntry fn_out_preprocess;
	// CTFFIND
	BooleanEntry do_ctffind;
	FileNameEntry fn_ctffind3_exe;
        SliderEntry ctf_win;
	SliderEntry cs, kv, q0, mag, dstep;
	AnyEntry mic_names;
	SliderEntry box, resmin, resmax, dfmin, dfmax, dfstep;
	// extract
	BooleanEntry do_extract, do_movie_extract;
	AnyEntry coord_names, movie_rootname;
	SliderEntry avg_movie_frames;
	SliderEntry first_movie_frame;
	SliderEntry last_movie_frame;
	SliderEntry extract_size;
	RadioEntry coord_format;
	// operate
	BooleanEntry do_rescale, do_norm, do_invert, do_starfile;
	SliderEntry rescale, bg_radius, white_dust, black_dust;


	// Refine stuff
	// I/O
	AnyEntry fn_out;
	FileNameEntry fn_cont;
	FileNameEntry fn_img;
	FileNameEntry fn_ref;
	BooleanEntry ref_correct_greyscale;
	SliderEntry particle_diameter;
	SliderEntry angpix;
	SliderEntry ini_high;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_corrected_ref;
	//BooleanEntry only_flip_phases;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	// Optimisation
	SliderEntry nr_iter;
	SliderEntry tau_fudge;
	FileNameEntry fn_mask;
	BooleanEntry do_zero_mask;

	// Movies
	BooleanEntry do_movies;
	FileNameEntry fn_movie_star;
	SliderEntry movie_runavg_window;
	SliderEntry movie_sigma_angles;
	//RadioEntry  movie_sampling;
	SliderEntry movie_sigma_offset;
	//SliderEntry movie_offset_step;

	// Sampling
	RadioEntry sym_group;
	SliderEntry sym_nr;
	SliderEntry nr_classes;
	RadioEntry sampling;
	SliderEntry psi_sampling;
	BooleanEntry do_local_ang_searches;
	RadioEntry auto_local_sampling;
	SliderEntry sigma_angles;
	SliderEntry offset_range;
	SliderEntry offset_step;
	textOnlyEntry autosample_text1, autosample_text2;

	// Running
	SliderEntry nr_mpi;
	SliderEntry nr_threads;
    BooleanEntry do_queue;
	AnyEntry queuename;
	AnyEntry qsub;
	FileNameEntry qsubscript;
	AnyEntry qsub_extra1;
	AnyEntry qsub_extra2;
	AnyEntry other_args;

	// Am I a 2D or a 3D run, and a new or a continue run?
	int run_type;
	bool is_continue, have_extra1, have_extra2;
	std::string output_name;
	int ori_w, ori_h;

	// For deactivation of some fields
	Fl_Group *ctf_group, *ctffind_group, *star_group, *extract_group, *rescale_group, *norm_group, *restart_group, *queue_group, *localsearch_group, *movie_group, *movie_extract_group;

	// For Tabs
	Fl_Menu_Bar *menubar;
	Fl_Tabs *tabs;
	Fl_Group *tab1, *tab2, *tab3, *tab4, *tab5;

	// Choice for runtype and new/continue
	Fl_Choice *choice_runtype, *choice_continue;
    Fl_Menu_ *menu_runtype, *menu_continue;

    // Run button
    Fl_Button *run_button, *print_CL_button, *cite_button;

	// Constructor with w x h size of the window and a title
	RelionMainWindow(int w, int h, const char* title);

    // Destructor
    ~RelionMainWindow(){};

    // Clear all preprocess stuff
    void clearPreprocess();

    // Redraw the content of preprocessing window
    void setupPreprocess();

    // Clear all refine stuff
    void clearRefine();

    // Redraw the content of the refinement windows
    void setupRefine(int newruntype, bool is_continue = false);

    // Redraw the whole window
    void setup(int newruntype, bool is_continue = false);

    // Activate and deactivate fields depending on is_continue boolean
    void toggle_new_continue();

    // Reset the current_y to the top
    void resetHeight();

    // Add a FileNameEntry to the window
    void AddAnyEntry(AnyEntry &var,
    		                 const char* label,
    		                 const char* defaultvalue,
    		                 const char* helptext = NULL );

    // Add a TextOnlyEntry to the window
    void AddTextOnlyEntry(textOnlyEntry &var,
    		                 const char* label);

    // Add a FileNameEntry to the window
    void AddFileNameEntry(FileNameEntry &var,
    		                 const char* label,
    		                 const char* defaultvalue,
    		                 const char* pattern,
    		                 const char* helptext = NULL );

    // Add a RadioEntry to the window
    void AddRadioEntry(RadioEntry &var,
					   const char* label,
					   Fl_Menu_Item *options,
					   Fl_Menu_Item* defaultvalue,
					   const char* helptext = NULL);

    // Add a BooleanEntry to the window
    void AddBooleanEntry(BooleanEntry &var,
    		             const char * title,
    		             bool defaultvalue,
    		             const char* helptext = NULL,
                         Fl_Group* deactivate_this_group = NULL);

    // Add a SliderEntry to the window
    void AddSliderEntry(SliderEntry &var,
    		            const char * title,
    		            float defaultvalue,
    		            float minvalue,
    		            float maxvalue,
    		            float valuestep,
						const char* helptext = NULL);

    // Generate the actual command line from all data
    std::string getCommandLine();

    // Generate the job submission script file
    void saveJobSubmissionScript(std::string newfilename);

    // Save current settings to a file
    void writeSettings(std::string filename);


    // This will add _ctX for continuation runs
    // The output_name will be used for --o, but also for stderr, stdout, submit script and gui settings filenames
    // Note the value of fn_out is not changed, because that would be reloaded again in the GUI, which would lead to repetitive additions of _ctX
    void setOutputName();

    // Save current settings to a file
    // If the file can be read this function returns true
    // otherwise it gives an error if do_check==true, or it return false if do_check==false
    bool readSettings(std::string filename, bool do_check = true);

    /** Run buttons
     */
    void AddRunButtons();


private:


    // Vertical distance from the top
    int start_y;
    // Height of each entry
    int step_y;
    // Current height
    int current_y;


    /** Call-back functions for the Run button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */
    static void cb_run(Fl_Widget*, void*);
    inline void cb_run_i();

    static void cb_print_cl(Fl_Widget*, void*);
    inline void cb_print_cl_i();

    static void cb_cite(Fl_Widget*, void*);
    inline void cb_cite_i();

    static void cb_menu_continue(Fl_Widget*, void*);
    inline void cb_menu_continue_i();

    static void cb_menu_runtype(Fl_Widget*, void*);
    inline void cb_menu_runtype_i();

    static void cb_menubar_load(Fl_Widget*, void*);
    inline void cb_menubar_load_i();

    static void cb_menubar_save(Fl_Widget*, void*);
    inline void cb_menubar_save_i();

    static void cb_menubar_reactivate_runbutton(Fl_Widget*, void*);
    inline void cb_menubar_reactivate_runbutton_i();

    static void cb_menubar_quit(Fl_Widget*, void*);
    inline void cb_menubar_quit_i();

    static void cb_menubar_about(Fl_Widget*, void*);
    inline void cb_menubar_about_i();

};
// General utility to replace strings in a text buffer.
void replaceString(Fl_Text_Buffer *textbuf, std::string replacethis, std::string replaceby);
#include "mainwindow.cpp"

#endif /* MAINWINDOW_H_ */
