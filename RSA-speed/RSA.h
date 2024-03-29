#pragma once
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <stack>
using namespace std;
typedef uint64_t dLN;
typedef uint32_t LN;
typedef uint16_t hLN;

#define LN_BITS 32
#define LN_HALF_BITS 16
#define LN_MAX_DIGITS 65
#define LN_MAX 0xffffffff
#define LN_HALF_MAX 0xffff
#define l 1024
#define Sk 2560

const uint8_t SMALL_PRIMES[] = { 3, 5, 7, 11 };
#define SMALL_PRIME_COUNT 4