#  -----------------------------------------------------------------------------
# 
#                      GNU General Public Licence
#           
#       This program is free software; you can redistribute it and/or
#       modify it under the terms of the GNU General Public License as
#       published by the Free Software Foundation; either version 2 of
#       the License, or (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public
#       License along with this program; see the file COPYING. If not,
#       write to the Free Software Foundation, Inc., 59
#       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
#  -----------------------------------------------------------------------------
#

#--------------------------------
# GENERAL VARIABLES
#--------------------------------

CC = gcc 
CFLAGS = -O3 -I. -Wall -lm #-pg -g

#--------------------------------
# LIBRARY SOURCES
#--------------------------------

LIB_SRC =  gnuplot_i.c 
LIB_OBJ = $(LIB_SRC:.c=.o)

all: analyse_envspec

analyse_envspec : $(LIB_OBJ) $(LIB_SRC) analyse_envspec.c
	$(CC) $(CFLAGS) $(LIB_OBJ) analyse_envspec.c -o analyse_envspec -lsndfile -lfftw3 -lm 

clean : 
	rm -f *.o 

cleanall :
	rm -f *.o *~ 

