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
#include "src/image.h"


// Get size of datatype
unsigned long  gettypesize(DataType type)
{
    unsigned long   size;

    switch ( type ) {
        case UChar: case SChar:  size = sizeof(char); break;
        case UShort: case Short: size = sizeof(short); break;
        case UInt:	 case Int:   size = sizeof(int); break;
        case Float:              size = sizeof(float); break;
        case Double:             size = sizeof(double); break;
        case Bool:				  size = sizeof(bool); break;
        default: size = 0;
    }

    return(size);
}

int datatypeString2Int(std::string s)
{
  toLower(s);
  if (!strcmp(s.c_str(),"uchar"))
  {
       return UChar;
  }
  else if (!strcmp(s.c_str(),"ushort"))
  {
    return UShort;
  }
  else if (!strcmp(s.c_str(),"short"))
  {
    return Short;
  }
  else if (!strcmp(s.c_str(),"uint"))
  {
    return UInt;
  }
  else if (!strcmp(s.c_str(),"int"))
  {
    return Int;
  }
  else if (!strcmp(s.c_str(),"float"))
  {
    return Float;
  }
  else REPORT_ERROR("datatypeString2int; unknown datatype");


}
