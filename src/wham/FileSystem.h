
#ifndef FILE_SYSTEM_H
#define FILE_SYSTEM_H

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <string>

// Misc. functions related to manipulating files
namespace FileSystem {

// Returns the directory separator
// - Unix: '/'
// - Windows: '\'
char separator() noexcept;

// Returns the location of filesystem root
std::string root();

// Concatenates the paths
std::string join(const std::string& left, const std::string right);

// Given a path, returns the absolute path (behaves like Linux 'realpath')
// TODO: Windows version
std::string realpath(const std::string& path);

// Get the directory which contains the given path
// - Taken from
//     C++ Cookbook by Jeff Cogswell, Jonathan Turkanis, Christopher Diggins, D. Ryan Stephens
std::string dirname(const std::string& full_path);

// Returns true if the path is with respect to the filesystem root
bool isAbsolutePath(const std::string& path);

/// Resolves a (possibly relative) path
/// @param rel_path  path to resolve (may be *relative* to base_path)
/// @param base_path used to resolve relative baths (acts like '.')
// TODO: set the default value of 'base_path' to the PWD
std::string resolveRelativePath(const std::string& rel_path, const std::string& base_path);

// Reads a file, and extracts a list of file paths from column file_col
// - If a file is listed as a relative path (e.g. using '..'), then the
//   path is evaluted relative to the location of files_list
void readFilesList(
	const std::string& files_list,
	const int file_col,  // indexed from 0
	std::vector<std::string>& files
);


} // end namespace FileSystem

#endif // ifndef FILE_SYSTEM_H
