/*
 * Copyright (C) 2013 Quantum ESPRESSO group
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt .
 */

/* C/C++ interface to the codes of the Quantum ESPRESSO package */

#ifndef QE_LIBCOUPLE_H
#define QE_LIBCOUPLE_H

#ifdef __cplusplus
extern "C" {
#endif


void pwstart(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
                int nband, int ndiag, int *exit_status, char *input_file);

void pwstep(int *exit_status);

void pwdata();

void pwend(int *exit_status);
    
#ifdef __cplusplus
}
#endif

#endif /* QE_LIBCOUPLE_H */
