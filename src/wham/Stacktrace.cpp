// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0
//
// Modifications by Sean M. Marks (https://github.com/seanmarks)

#include "Stacktrace.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

// TODO: check for availability of headers
#include <execinfo.h>  // backtrace
#include <cxxabi.h>    // __cxa_demangle
//#include <dlfcn.h>     // dladdr


namespace DebugTools {


std::string stacktrace(const unsigned int max_num_frames, const unsigned int skip)
{
	std::ostringstream ss;
	ss << "stack trace:\n";

	// storage array for stack trace address data
	void* addr_list[max_num_frames+1];

	// retrieve current stack addresses
	unsigned int num_frames = backtrace(addr_list, sizeof(addr_list) / sizeof(void*));

	if (num_frames == 0) {
		ss << "  <empty, possibly corrupt>\n";
		return ss.str();
	}

	// resolve addresses into strings containing "filename(function+address)",
	// this array must be free()-ed
	char** symbol_list = backtrace_symbols(addr_list, num_frames);

	ss << "RAW\n";
	for ( unsigned int i = skip; i < num_frames; ++i ) {
		ss << i << "  " << symbol_list[i] << "\n";
	}

	// allocate string which will be filled with the demangled function name
	size_t funcnamesize = 256;
	char* funcname = (char*)malloc(funcnamesize);

	char buffer[1024];  // for snprintf

	try {
		ss << "PROCESSED\n";
		// iterate over the returned symbol lines. skip the first, it is the
		// address of this function.
		for (unsigned i = skip; i < num_frames; i++) {
			char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

			// find parentheses and +address offset surrounding the mangled name:
			// ./module(function+0x15c) [0x8048a6d]
			for (char *p = symbol_list[i]; *p; ++p) {
				if (*p == '(') {
					begin_name = p;
				}
				else if (*p == '+') {
					begin_offset = p;
				}
				else if (*p == ')' && begin_offset) {
					end_offset = p;
					break;
				}
			}

			if (begin_name && begin_offset && end_offset
					&& begin_name < begin_offset
			) {
				*begin_name++   = '\0';
				*begin_offset++ = '\0';
				*end_offset     = '\0';

				// mangled name is now in [begin_name, begin_offset) and caller
				// offset in [begin_offset, end_offset). now apply
				// __cxa_demangle():

				int status;
				char* ret = abi::__cxa_demangle(begin_name,
																				funcname, &funcnamesize, &status);
				if (status == 0) {
					funcname = ret; // use possibly realloc()-ed string
					snprintf(buffer, sizeof(buffer), "%-3d %s : %s+%s\n",
					         i, symbol_list[i], funcname, begin_offset);
					ss << buffer;
				}
				else {
					// demangling failed. Output function name as a C function with
					// no arguments.
					snprintf(buffer, sizeof(buffer), "%-3d %s : %s()+%s\n",
					         i, symbol_list[i], begin_name, begin_offset);
					ss << buffer;
				}
			}
			else {
				// couldn't parse the line? print the whole line.
				snprintf(buffer, sizeof(buffer), "  %s\n", symbol_list[i]);
				ss << buffer;
			}
		}
	}
	catch (const std::exception& ex) {
		free(funcname);
		free(symbol_list);
		std::cerr << "Error in debug::stacktrace()";
		throw;
	}

	// Cleanup
	free(funcname);
	free(symbol_list);

	if (num_frames == max_num_frames) {
		ss << "[truncated]\n";
	}

	return ss.str();
}

} // end namespace debug
