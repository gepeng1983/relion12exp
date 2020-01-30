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
#include "showhelp.h"

ShowHelpText::ShowHelpText(const char *help)
{
    int w=640;
    int h=480;
	Fl_Window *win = new Fl_Window(w, h);
    Fl_Text_Buffer *buff = new Fl_Text_Buffer();
    Fl_Text_Display *disp = new Fl_Text_Display(20, 20, w-40, h-40, "relion additional text.");
    disp->buffer(buff);
    disp->wrap_mode(1,79);
    win->resizable(*disp);
    win->show();
    buff->text(help);
}

ShowHelpText::~ShowHelpText(){};
