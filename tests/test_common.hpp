#pragma once
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

struct TestRegistry {
    std::vector<std::pair<std::string, std::function<void()>>> tests;
    static TestRegistry& instance();
};

struct RegisterTest {
    RegisterTest(const std::string& name, std::function<void()> fn) {
        TestRegistry::instance().tests.push_back({name, fn});
    }
};

#define TEST_CASE(name) static void name(); static RegisterTest reg_##name(#name, name); static void name()
#define REQUIRE(cond) do { if (!(cond)) { std::cerr << "FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << "\n"; std::exit(1); } } while(false)
#define REQUIRE_NEAR(a,b,eps) do { if (std::abs((a)-(b)) > (eps)) { std::cerr << "FAILED near: " #a " vs " #b " at " << __FILE__ << ":" << __LINE__ << " values " << (a) << " " << (b) << "\n"; std::exit(1); } } while(false)
