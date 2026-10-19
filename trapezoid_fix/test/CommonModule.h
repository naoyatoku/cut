#pragma once
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <limits>
#include <cstdint>
typedef long long __int64;
struct AssertFail : std::runtime_error { AssertFail(const char*m):std::runtime_error(m){} };
#define _Assert(c, ...) do{ if(!(c)){ char _b[256]; snprintf(_b,sizeof _b,__VA_ARGS__); throw AssertFail(_b);} }while(0)
