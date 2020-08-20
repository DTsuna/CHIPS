#ifndef INCLUDED_FILE_RW_H
#define INCLUDED_FILE_RW_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "main.h"

void readin_pars(pars*);
void mk_dir_nsd(pars);
void file_open_pars(FILE**, pars, const char*, const char*);

#endif
