#include "test_common.hpp"

TestRegistry& TestRegistry::instance() {
    static TestRegistry r;
    return r;
}

int main() {
    int passed = 0;
    for (auto& [name, fn] : TestRegistry::instance().tests) {
        std::cout << "[ RUN      ] " << name << "\n";
        fn();
        std::cout << "[       OK ] " << name << "\n";
        passed++;
    }
    std::cout << passed << " tests passed\n";
    return 0;
}
