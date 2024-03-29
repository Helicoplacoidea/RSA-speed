#pragma once
#include "RSA.h"

void initialize_rand(void);
void generate_rand(uint8_t* block, uint32_t block_len);
int probable_prime(LN* a, uint32_t digits);
int fermat_test(LN* a, uint32_t digits);			// fermat素性检验，a为待检验数
int miller_Rabin(LN* n, int t, uint32_t digits); // miller-rabin素性检验，n为待检验数，t为次数

void prime_genarate_table(LN* prime);					//使用查表法的素数生成函数
void bit_array(LN* prime, int len);								//使用bit_array方法的素数生成函数
void prime_genarate(LN* p);								//普通的素数生成函数