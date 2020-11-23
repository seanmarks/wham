#include "FileSystem.h"
#include "Assert.hpp"


namespace FileSystem {


char separator() noexcept
{
#ifndef _WIN_32
	return '/';  // Unix
#else
	return '\\';  // Windows
#endif
}


std::string root()
{
#ifndef _WIN_32
	std::stringstream ss;
	ss << separator();
	return ss.str();
#else
	static_assert(false, "missing implementation");
#endif
}


std::string join(const std::string& left, const std::string right)
{
	std::stringstream ss;
	ss << left << FileSystem::separator() << right;
	return ss.str();
}


std::string realpath(const std::string& path)
{
#ifndef _WIN_32
	if ( path.length() < 1 ) {
		throw std::runtime_error("get_realpath() was given an empty path");
	}

	// Use POSIX realpath()
	char* buffer = ::realpath(&path[0], nullptr);
	if ( buffer == nullptr ) {
		throw std::runtime_error("Error resolving path \"" + path + "\"");
	}

	// Move the path to a std string and clean up
	std::string resolved_path(buffer);
	free(buffer); buffer = nullptr;

	return resolved_path;

#else
	// TODO: Windows version
	static_assert(false, "missing implementation");
#endif // ifndef _WIN_32
}


std::string dirname(const std::string& full_path)
{
	size_t i = full_path.rfind(FileSystem::separator(), full_path.length());
	if ( i != std::string::npos ) {
		return full_path.substr(0, i);
	}
	else {
		return "."; // TODO: Windows
	}
}


bool isAbsolutePath(const std::string& path)
{
	FANCY_ASSERT( ! path.empty(), "no input provided" );

#ifndef _WIN32
	return ( path.front() == FileSystem::separator() );
#else
	static_assert(false, "missing implementation");
#endif // ifndef _WIN32
}


std::string resolveRelativePath(const std::string& rel_path, const std::string& base_path)
{
	std::stringstream ss;
	if ( isAbsolutePath(rel_path) ) {
		ss << rel_path;
	}
	else {
		ss << join(base_path, rel_path);
	}

	return realpath( ss.str() );
}


void readFilesList(
	const std::string& files_list, const int file_col, std::vector<std::string>& files) 
{
	std::ifstream list_ifs(files_list);
	if ( not list_ifs.is_open() ) {
		throw std::runtime_error("Failed to open file \'" + files_list + "\'");
	}

	std::string files_list_path = FileSystem::dirname(files_list);

	const char comment_char = '#';
	auto skip = [=](const std::string& s) {
		return ( s.empty() || s.front() == comment_char );
	};

	int i = 0;
	std::string line, token;
	std::vector<std::string> tokens;
	files.clear();
	while ( getline(list_ifs, line) ) {
		// Ignore possible comments and blank lines
		std::stringstream ss(line);
		ss >> token;
		if ( skip(token) ) {
			continue;
		}

		// Split the line into tokens
		tokens = {{ token }};
		while ( ss >> token ) {
			if ( ! skip(token) ) {
				tokens.push_back( token );
			}
			else {
				break;
			}
		}

		// Check for too few columns
		const int num_tokens = tokens.size();
		FANCY_ASSERT( file_col < num_tokens,
				"Error reading files list from " << files_list << "\n"
			  << "  A non-empty, non-comment line has fewer than " << file_col+1 << " columns\n" 
				<< "  line: " << line << "\n" );

		// Process the file name
		auto file = resolveRelativePath(tokens[file_col], files_list_path);
		files.push_back( file );
		++i;
	}
	list_ifs.close();
}
} // end namespace FileSystem