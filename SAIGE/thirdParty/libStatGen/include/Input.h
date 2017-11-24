/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INPUT_H__
#define __INPUT_H__

void Input(const char * prompt, int & n, int _default = 0);
void Input(const char * prompt, double & d, double _default = 0.0);
void Input(const char * prompt, char & c, char _default = 'A');
void Input(const char * prompt, char * s, const char * _default = "");
void Input(const char * prompt, bool & b, bool _default);

void InputBounds(const char * prompt, int & n, int  min, int max,
                 int _default = 0);
void InputBounds(const char * prompt, double & d, double min, double max,
                 double _default = 0);

extern int InputPromptWidth;

#endif
