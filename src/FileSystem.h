// FileSystem: static methods for evaluating file and directory paths
// - TODO Support for Windows

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


class FileSystem {
 public:
	// Directory separator
#ifndef _WIN_32
	static const char sep = '/';  // Unix
#else
	static const char sep = '\\';  // Windows
#endif

	// Reads a file, and extracts a list of file paths from column file_col
	// - If a file is listed as a relative path (e.g. using '..'), then the
	//   path is evaluted relative to the location of files_list
	static void readFilesList(
		const std::string& files_list,
		const int file_col,  // indexed from 0
		std::vector<std::string>& files
	);

	// TODO windows
#ifndef _WIN_32
	static std::string get_realpath(const std::string& path)
	{
		if ( path.length() < 1 ) {
			throw std::runtime_error("get_realpath() was given an empty path");
		}

		// Use POSIX realpath()
		char* buffer = realpath(&path[0], nullptr);
		if ( buffer == nullptr ) {
			throw std::runtime_error("Error resolving path \"" + path + "\"");
		}

		// Move the path to a std string and clean up
		std::string resolved_path(buffer);
		free(buffer); buffer = nullptr;

		return resolved_path;
	}
#endif /* _WIN_32 */

	// Get the directory which contains the given path
	// - Taken from:
	//     C++ Cookbook by Jeff Cogswell, Jonathan Turkanis, Christopher Diggins, D. Ryan Stephens
	static std::string get_basename(const std::string& full_path)
	{
		size_t i = full_path.rfind(FileSystem::sep, full_path.length());
		if ( i != std::string::npos ) {
			return full_path.substr(0, i);
		}
		else {
			return "."; // TODO Windows
		}
	}
};

#endif /* FILE_SYSTEM_H */
