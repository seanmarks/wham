#include "FileSystem.h"

void FileSystem::readFilesList(
	const std::string& files_list, const int file_col, std::vector<std::string>& files) 
{
	std::ifstream list_ifs(files_list);
	if ( not list_ifs.is_open() ) {
		throw std::runtime_error("Failed to open file \'" + files_list + "\'");
	}

	std::string files_list_path = FileSystem::get_basename(files_list);

	int i = 0;
	std::string line, token, file;
	std::vector<std::string> tokens;
	files.clear();
	const char comment_char = '#';
	while ( getline(list_ifs, line) ) {
		// Ignore possible comments and blank lines
		std::stringstream ss(line);
		ss >> token;
		if ( line.empty() or token[0] == comment_char ) {
			continue;  // Comment or empty line
		}

		// Split the line into tokens
		tokens = {{ token }};
		while ( ss >> token ) {
			if ( token[0] != comment_char ) {
				tokens.push_back( token );
			}
			else {
				break;  // Rest of the line is a comment
			}
		}

		// Check for too few columns
		if ( file_col > static_cast<int>(tokens.size()) ) {
			std::stringstream err_ss;
			err_ss << "Error reading files list from " << files_list << "\n"
			       << "  A non-empty, non-comment line has fewer than " << file_col+1 << " columns\n";
			throw std::runtime_error( err_ss.str() );
		}

		// Process the file name
		file = tokens[file_col];
		if ( file[0] != FileSystem::sep ) {
			// Location is relative to path to files list
			file = files_list_path + std::string(1,FileSystem::sep) + file;
		}
		file = FileSystem::get_realpath(file);
		files.push_back(file);

		++i;
	}
	list_ifs.close();
}
