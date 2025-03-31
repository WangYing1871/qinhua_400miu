// fec2det_dec.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <TFile.h>
#include <ROOT/RConfig.hxx>

#include <sys/types.h>
#include <sys/stat.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "getopt.hh"
#include "constant.hh"
#include "hit_reconstruct.hh"


// TODO: Reference additional headers your program requires here.
void show_info(char* name);

struct stat info;

int c = 0;
extern int optind, opterr, optopt;
extern char* optarg;
void show_info(char* name);
vector<string> split_string(const string& str, string delims, int inc = 0);
