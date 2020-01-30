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
#ifndef ENTRIES_H_
#define ENTRIES_H_

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Choice.H>
#include "showhelp.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <fstream>

#ifndef ABS
#define ABS(x) (((x) >= 0) ? (x) : (-(x)))
#endif
#ifndef FLOOR
#define FLOOR(x) (((x) == (int)(x)) ? (int)(x):(((x) > 0) ? (int)(x) : \
                  (int)((x) - 1)))
#endif

static Fl_Menu_Item bool_options[] = {
			      {"Yes"},
			      {"No"},
			      {0} // this should be the last entry
			      };

class textOnlyEntry{

public:
	Fl_Text_Display* mydisp;
	Fl_Text_Buffer *textbuff;
	bool has_been_set;

	textOnlyEntry()
	{
		has_been_set=false;
	}
	~textOnlyEntry(){};

	void initialise(int x, int y, int width, int height, const char* text)
	{
		mydisp = new Fl_Text_Display(XCOL1, y, width, height);
		textbuff = new Fl_Text_Buffer();
		textbuff->text(text);
		mydisp->buffer(textbuff);
		mydisp->color(GUI_BACKGROUND_COLOR);
		has_been_set=true;
	}

	void clear()
	{
		if (has_been_set)
		{
			delete mydisp;
			delete textbuff;
			has_been_set = false;
		}
	}
};

/** This is the main class to generate input entry-lines in the Gui windows.
 *  It implements three columns to be displayed:
 *  1. box with the label
 *  2. Input field with the input value
 *  3. Help button that pops up a window with additional help text
 *
 *  All specific entries (e.g. to get FileName, Boolean, etc. inherit from this class)
 *
 *
 */
class AnyEntry{

public:
    // Input value storage
	Fl_Input* inp;

	// Label
	std::string label;

    // Button to show additional help text
	Fl_Button* help;

	// The additional help text
    const char *myhelptext;

    /** Constructor with x,y-position from top left
	 *  wcol1, wcol2 and wcol3 are the widths of the three columns described above
	 *  title is the value displayed in the first column
	 *  defaultvalue is what will appear by default in the input value
	 *  help is the additional help text. If it is set to NULL, no help button will be displayed
	 */
	AnyEntry(){};

    /** Empty destructor
     */
	~AnyEntry(){};

	/** Here really start the entry
	 */
	void initialise(int x, int y, int height, int wcol2, int wcol3, const char* title, const char* defaultvalue = NULL, const char* help = NULL);

	// Get the value
    std::string getValue();

    // Set the value
    void setValue(const char* inp);

    // Clear this entry
	void clear();

    // Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Save the value to a file
    void writeValue(std::ostream& out);

    // Read the value from a file
    void readValue(std::ifstream& in);

    /** Call-back functions for the help button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */
    static void cb_help(Fl_Widget*, void*);
    inline void cb_help_i();
};


// Get a FileName value from the user (with browse button).
class FileNameEntry: public AnyEntry
{

public:
	// Browse button
    Fl_Button* browse;

    const char* pattern;

    // Constructor (with 4 column widths)
	FileNameEntry() {};

    // Destructor
	~FileNameEntry(){};

	void initialise(int x, int y, int height,
    		int wcol2, int wcol3, int wcol4,
    		const char* title,
    		const char* defaultvalue,
    		const char* _pattern = "",
    		const char* help = NULL);

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);


private:
    // Call-back functions for the browse button
    static void cb_browse(Fl_Widget*, void*);
    inline void cb_browse_i();

};

// Get an entry from a list of possible values from the user.
class RadioEntry: public AnyEntry
{
public:

    // The choices
    Fl_Choice * choice;
    // The menu
    Fl_Menu_* menu;
    // Deactivate this group
    Fl_Group * my_deactivate_group;

    // Constructor
    RadioEntry(){};

    // Destructor
    ~RadioEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 Fl_Menu_Item *options,
				 Fl_Menu_Item* defaultvalue,
				 const char* help = NULL,
				 Fl_Group * deactivate_this_group = NULL);


    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Get the value
    std::string getValue();

    // Read the value from a file
    void readValue(std::ifstream& in);



public: // this one is public so that it can be called in mainwindow to deactivate default groups
    static void cb_menu(Fl_Widget*, void*);
    inline void cb_menu_i();
};

class BooleanEntry: public RadioEntry
{
public:
	// Constructor
	BooleanEntry(){};

	// Destructor
    ~BooleanEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 bool defaultvalue,
				 const char* help = NULL,
				 Fl_Group * deactivate_this_group = NULL);

    // Get the value
    bool getValue();


};

class SliderEntry:  public RadioEntry
{

public:
    // The slider
    Fl_Slider * slider;

    // Constructor
	SliderEntry(){};

	// Destructor
    ~SliderEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 float defaultvalue,
                 float minvalue,
                 float maxvalue,
                 float valuestep,
				 const char* help = NULL);

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Get the value
    float getValue();

    // Read the value from a file
    void readValue(std::ifstream& in);


private:
    static void cb_slider(Fl_Widget*, void*);
    inline void cb_slider_i();

    static void cb_input(Fl_Widget*, void*);
    inline void cb_input_i();


};

// Utilities
float textToFloat(const char* mytext);
/** Float to string conversion.
 *
 * If precision==0 the precision is automatically computed in such a way that
 * the number fits the width (the exponential format might be chosen). If
 * precision==-1 then the exponential format is forced. If width==0 then the
 * minimum width is used.
 *
 * @code
 * REPORT_ERROR(1602, "Value not recognised " + floatToString(val));
 * @endcode
 */
std::string floatToString(float F, int _width = 0, int _prec = 0);
int bestPrecision(float F, int _width);

#include "entries.cpp"
#endif /* ENTRIES_H_ */
