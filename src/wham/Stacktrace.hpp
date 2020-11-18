//
// Modifications by Sean M. Marks (https://github.com/seanmarks)
//
// Inspiration from:
// - Farooq Mela (https://github.com/fmela, https://gist.github.com/fmela)

#ifndef STACKTRACE_H
#define STACKTRACE_H

#include <string>
#include <vector>

namespace DebugTools {

// Returns a string with a demangled stack backtrace of the caller function
std::string stacktrace(
	const unsigned int max_frames = 127,
	const unsigned int skip = 1
);

} // end namespace DebugTools

#endif // ifndef STACKTRACE_H
