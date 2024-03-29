#pragma once
#include "RSA.h"

void MonPro(LN* res, LN* a, LN* b, LN* r, LN* n, LN* n_1, uint32_t rbits, uint32_t digits); // Montgomery约简函数，res = (a * b * r_1) mod n
void MonMul(LN* res, LN* a, LN* b, LN* n, LN* n_1, LN* r, uint32_t digits); // Montgomery模乘函数，res = (a * b) mod n
void ModExp(LN* res, LN* M, LN* e, LN* n, uint32_t digits);					// Montgomery模幂函数，res = (M ^ e) mod n
void ModExp_SW(LN* res, LN* M, LN* e, LN* n, uint32_t digits);
void MonPro_IFIOS(LN* res, LN* a, LN* b, LN* n, LN* n_1, uint32_t digits);