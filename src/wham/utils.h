// Miscellaneous macros and constexpr functions

#ifndef UTILS_H
#define UTILS_H

#include <exception>
#include <string>
#include <stdexcept>

// Macro to evaluate member function using a pointer to it
// - One should generally be highly allergic to macros, but isocpp.org states that 
//   this is one of the very few exceptions to the rule
#define CALL_MEMBER_FXN(object,pMethod) ((object).*(pMethod))

// Calculate n! recursively
constexpr int factorial(int n) {
	return n > 0 ? n*factorial(n-1) : 1;
}

// Convert 'x' to a string using arcane preprocessor tricks
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

// Prints the file name and line number where it's expanded
#define LOCATION_IN_SOURCE __FILE__ ":" TO_STRING(__LINE__)

// Prints the file name and line number where it's expanded
#define LOCATION_IN_SOURCE_STRING std::string(LOCATION_IN_SOURCE)

// Many (but not all) compilers provide a definition: use this as a fallback
#ifndef __PRETTY_FUNCTION__
#  define __PRETTY_FUNCTION__ ""
#endif

// Prints the "pretty" name of the function along in which the macro is expanded,
// along with the associated file name and line number
#define FANCY_FUNCTION __PRETTY_FUNCTION__ " (" LOCATION_IN_SOURCE ")"

#define WHAM_ASSERT(test,message) if (not (test)) { \
  	throw std::runtime_error("assertion failed in " FANCY_FUNCTION "\n" "  " message "\n" "  test: " STRINGIFY(test) "\n"); \
	}

#endif // ifndef UTILS_H
