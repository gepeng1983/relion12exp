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

#include "entries.h"
#include <iostream>

// ==============================================================================
// AnyEntry =====================================================================
// ==============================================================================

void AnyEntry::initialise(int x, int y, int height,
				   int wcol2, int wcol3,
				   const char* title,
				   const char* defaultvalue,
				   const char* helptext)
{
    // Set the label
    label = title;

    // The input field
	if (defaultvalue != NULL)
	{
		inp = new Fl_Input(XCOL2, y, wcol2, height, title);

		// Set the input value
		inp->value(defaultvalue);
		inp->color(GUI_INPUT_COLOR);
	}

	// Display help button if needed
    if (helptext != NULL)
    {
    	// Set the help text
		myhelptext = helptext;

		// The Help button
		help = new Fl_Button( XCOL3, y, wcol3, height, "?");
		help->callback( cb_help, this );
		help->color(GUI_BUTTON_COLOR);

    }
}

std::string AnyEntry::getValue()
{
	return (std::string)inp->value();
}

void AnyEntry::setValue(const char* val)
{
	inp->value(val);
}

void AnyEntry::writeValue(std::ostream& out)
{
	// Only write entries that have been initialised
	if (label != "")
		out << label << " == " << getValue() << std::endl;
}


void AnyEntry::readValue(std::ifstream& in)
{
    if (label != "")
    {
		// Start reading the ifstream at the top
		in.clear(); // reset eof if happened...
    	in.seekg(0, std::ios::beg);
		std::string line;
		while (getline(in, line, '\n'))
		{
			if (line.rfind(label) == 0)
			{
				// found my label
				int equalsigns = line.rfind("==");
				std::string newval = line.substr(equalsigns + 3, line.length() - equalsigns - 3);
				inp->value(newval.c_str());
				return;
			}
		}
    }
}

void AnyEntry::clear()
{
	if (label != "")
	{
		label="";
		delete inp;
		delete help;
	}
}

void AnyEntry::deactivate(bool do_deactivate)
{
	if (do_deactivate)
	{
		inp->deactivate();
		help->deactivate();
	}
	else
	{
		inp->activate();
		help->activate();
	}

}

// Help button call-back functions
void AnyEntry::cb_help(Fl_Widget* o, void* v) {

    AnyEntry* T=(AnyEntry*)v;
    T->cb_help_i();
}

void AnyEntry::cb_help_i() {

    ShowHelpText *help = new ShowHelpText(myhelptext);

}


// ==============================================================================
// FileNameEntry ================================================================
// ==============================================================================

void FileNameEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     const char* defaultvalue,
		                     const char* _pattern,
		                     const char* helptext)
{

	AnyEntry::initialise(x,y,height,
			wcol2,wcol3,
			title,
			defaultvalue,
			helptext);

	// Store the pattern for the file chooser
	pattern = _pattern;

    // The Browse button
    browse = new Fl_Button( XCOL4, y, WCOL4, height, "Browse");
    browse->callback( cb_browse, this );
    browse->color(GUI_BUTTON_COLOR);
}

void FileNameEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		delete browse;
	}
}

void FileNameEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		browse->deactivate();
	}
	else
	{
		browse->activate();
	}

}

void FileNameEntry::cb_browse(Fl_Widget* o, void* v) {

    FileNameEntry* T=(FileNameEntry*)v;
    T->cb_browse_i();
}


void FileNameEntry::cb_browse_i() {

    Fl::scheme("gtk+");
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", pattern, Fl_File_Chooser::SINGLE, "");

    G_chooser->directory(NULL);
    G_chooser->color(GUI_BACKGROUND_COLOR);
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

    inp->value(G_chooser->value());
}

// ==============================================================================
// RadioEntry ================================================================
// ==============================================================================
void RadioEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     Fl_Menu_Item *options,
		                     Fl_Menu_Item* defaultvalue,
		                     const char* helptext,
		    				 Fl_Group * deactivate_this_group)
{
	AnyEntry::initialise(x,y,height,
			wcol2,wcol3,
			title,
			defaultvalue->label(),
			helptext);

    // Pull-down menu button
	//Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", "", Fl_File_Chooser::SINGLE, "");

	my_deactivate_group = deactivate_this_group;
	choice = new Fl_Choice(XCOL2, y, WCOL2, height);
    choice->menu(options);
    choice->picked(defaultvalue);
    choice->callback(cb_menu, this);

    menu = choice;
    //menu->color(GUI_BACKGROUND_COLOR);
    menu->color(GUI_INPUT_COLOR);
}

void RadioEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		//delete choice;
		delete menu;
	}
}

void RadioEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		menu->deactivate();
	}
	else
	{
		menu->activate();
	}

}

std::string RadioEntry::getValue()
{
	return (std::string)inp->value();
}

void RadioEntry::readValue(std::ifstream& in)
{
	if (label != "")
	{
		AnyEntry::readValue(in);
		const Fl_Menu_Item *p = choice->find_item(inp->value());
		if ( p )
			choice->picked(p);
		else
			std::cerr << "Error readValue: Menu item not found:" << inp->value()<< std::endl;
	}
}


void RadioEntry::cb_menu(Fl_Widget* o, void* v) {

    RadioEntry* T=(RadioEntry*)v;
    T->cb_menu_i();
}


void RadioEntry::cb_menu_i() {

	const Fl_Menu_Item* m = menu->mvalue();
	// Set my own value
	inp->value(m->label());

	// In case this was a boolean that deactivates a group, do so:
	if (my_deactivate_group != NULL)
	if (strcmp(inp->value(), "No") == 0)
		my_deactivate_group->deactivate();
	else
		my_deactivate_group->activate();

}

// ==============================================================================
// BooleanEntry ================================================================
// ==============================================================================
void BooleanEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     bool defaultvalue,
		                     const char* helptext,
		    				 Fl_Group * deactivate_this_group)
{

	Fl_Menu_Item* defval;

	if (defaultvalue)
		defval = &bool_options[0];
	else
		defval = &bool_options[1];
	RadioEntry::initialise(x,y,height,
			wcol2,wcol3,wcol4,
			title,
			bool_options,
			defval,
			helptext,
			deactivate_this_group);

}

bool BooleanEntry::getValue()
{
	if (strcmp(inp->value(), "Yes") == 0)
		return true;
	else
		return false;
}

// ==============================================================================
// SliderEntry ================================================================
// ==============================================================================
void SliderEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     float defaultvalue,
		                     float minvalue,
		                     float maxvalue,
		                     float valuestep,
		                     const char* helptext)
{

	int floatwidth = 50;
	AnyEntry::initialise(x,y,height,
			floatwidth,wcol3,
			title,
			"",
			helptext);

	// Initialise label
	label = title;

	// Slider is shorter than wcol2, so that underlying input field becomes visible
	slider = new Fl_Slider(XCOL2 + floatwidth, y, wcol2 - floatwidth, height);
	slider->type(1);
	slider->callback(cb_slider, this);
	slider->minimum(minvalue);
	slider->maximum(maxvalue);
	slider->step(valuestep);
	slider->type(FL_HOR_NICE_SLIDER);
	slider->color(GUI_BACKGROUND_COLOR);
	inp->callback(cb_input, this);
	inp->when(FL_WHEN_ENTER_KEY|FL_WHEN_NOT_CHANGED);

	// Set the default in the input and the slider:
	std::string str = floatToString(defaultvalue);
	inp->value(str.c_str());
	slider->value(defaultvalue);

}

void SliderEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		delete slider;
	}
}

void SliderEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		slider->deactivate();
	}
	else
	{
		slider->activate();
	}

}


float SliderEntry::getValue()
{
	return textToFloat(inp->value());
}

void SliderEntry::readValue(std::ifstream& in)
{
	if (label != "")
	{
		AnyEntry::readValue(in);
		// Also reset the slider
		slider->value(textToFloat(inp->value()));
	}
}

void SliderEntry::cb_slider(Fl_Widget* o, void* v) {

    SliderEntry* T=(SliderEntry*)v;
    T->cb_slider_i();
}


void SliderEntry::cb_slider_i() {

    static int recurse = 0;
    if ( recurse ) {
        return;
    } else {
        recurse = 1;
        std::string str = floatToString(slider->value());
        inp->value(str.c_str());
        slider->redraw();
        recurse = 0;
    }
}

void SliderEntry::cb_input(Fl_Widget* o, void* v) {

    SliderEntry* T=(SliderEntry*)v;
    T->cb_input_i();
}


void SliderEntry::cb_input_i() {

    static int recurse = 0;
    if ( recurse ) {
        return;
    } else {
        recurse = 1;
        //slider->value(textToFloat(my_input->value()));         // pass input's value to slider
        slider->value(textToFloat(inp->value()));         // pass input's value to slider
        //inp->value(my_input->value()); // also set the normal input for getValue!
        recurse = 0;
    }
}

float textToFloat(const char* mytext)
{
	float retval;

    int ok;
    ok = sscanf(mytext, "%f", &retval);
    if (ok)
        return retval;
    else
    	std::cerr << "Error: cannot convert string to float..." << std::endl;
}

std::string floatToString(float F, int _width, int _prec)
{
    std::ostringstream outs;
    outs.fill(' ');

    if (_width != 0)
        outs.width(_width);

    if (_prec == 0)
        _prec = bestPrecision(F, _width);

    if (_prec == -1 && _width > 7)
    {
        outs.precision(_width - 7);
        outs.setf(std::ios::scientific);
    }
    else
        outs.precision(_prec);

    outs << F;

    std::string retval = outs.str();
    int i = retval.find('\0');

    if (i != -1)
        retval = retval.substr(0, i);

    return retval;
}

int bestPrecision(float F, int _width)
{
    // If it is 0
    if (F == 0)
        return 1;

    // Otherwise
    int exp = FLOOR(log10(ABS(F)));
    int advised_prec;

    if (exp >= 0)
        if (exp > _width - 3)
            advised_prec = -1;
        else
            advised_prec = _width - 2;
    else
    {
        advised_prec = _width + (exp - 1) - 3;
        if (advised_prec <= 0)
            advised_prec = -1;
    }

    if (advised_prec < 0)
        advised_prec = -1; // Choose exponential format

    return advised_prec;
}

