#include <boost/version.hpp>
#if BOOST_VERSION < 105600
#error Installed boost is too old!
#endif
int main() { return 0; }
