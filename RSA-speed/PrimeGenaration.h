#pragma once
#include "RSA.h"

void initialize_rand(void);
void generate_rand(uint8_t* block, uint32_t block_len);
int probable_prime(LN* a, uint32_t digits);
int fermat_test(LN* a, uint32_t digits);			// fermat���Լ��飬aΪ��������
int miller_Rabin(LN* n, int t, uint32_t digits); // miller-rabin���Լ��飬nΪ����������tΪ����

void prime_genarate_table(LN* prime);					//ʹ�ò�����������ɺ���
void bit_array(LN* prime, int len);								//ʹ��bit_array�������������ɺ���
void prime_genarate(LN* p);								//��ͨ���������ɺ���