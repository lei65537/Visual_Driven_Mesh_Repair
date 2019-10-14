// stdafx.cpp : source file that includes just the standard includes
// MeshViewer.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"

// TODO: reference any additional headers you need in STDAFX.H
// and not in this file


double _sw = 0.1;
double _uw = 1.0;
double _bw = 1.0;
double _eps_vis = 0.0056;
int _tmer=1;
int _opt_max_iter=-1;

ofstream *_log_off=NULL;
char * _file_under_processing = NULL;
int _num_threads = 6;

