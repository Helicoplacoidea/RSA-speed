#pragma once
#include "RSA.h"

#define ln_equal(a, b, digits) (!ln_cmp(a, b, digits))

#define shift_2_bits(x) (uint32_t)(((x) >> (LN_BITS - 2)) & 0x03)
#define shift_4_bits(x) (uint32_t)(((x) >> (LN_BITS - 4)) & 0x0f)

#define assign_one_digit(x, y, digits) \
	{                                  \
		ln_assign_zero(x, digits);     \
		x[0] = y;                      \
	}

// ��λ����
LN ln_l_shift(LN* a, LN* b, uint32_t c, uint32_t digits); // a = b << c
LN ln_r_shift(LN* a, LN* b, uint32_t c, uint32_t digits); // a = b >> c

// �ַ��������ת��
void str_to_ln(uint8_t* str, LN* ln, uint32_t digits, uint32_t size);
void ln_to_str(uint8_t* str, LN* ln, uint32_t digits, uint32_t size);

// ��ֵ
void ln_assign_zero(LN* a, uint32_t digits);   // a = 0
void ln_assign(LN* a, LN* b, uint32_t digits); // a = b

int ln_cmp(LN* a, LN* b, uint32_t digits);				// compare a and b
int ln_is_zero(LN* a, uint32_t digits);					// a ?= 0
int ln_is_one(LN* a, uint32_t digits);					// a ?= 1
void ln_assign_2exp(LN* a, uint32_t b, uint32_t digits);// a = 2 ^ b

// ��ȡbit���Ϳ���
uint32_t ln_get_digits(LN* a, uint32_t digits);
uint32_t ln_get_bits(LN* a, uint32_t digits);

// ��������
LN ln_add(LN* a, LN* b, LN* c, uint32_t digits);							 // a = b + c
LN ln_sub(LN* a, LN* b, LN* c, uint32_t digits);							 // a = b - c
void ln_mul(LN* a, LN* b, LN* c, uint32_t digits);							 // a = b * c
void ln_div(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits); // a = c / d  b = c % d

void ln_mod(LN* a, LN* b, uint32_t bdigits, LN* c, uint32_t cdigits);			 // a = b % c
void ln_mod_mul(LN* a, LN* b, LN* c, LN* d, uint32_t digits);					 // a = (b * c) % d
void ln_mod_exp(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits); // a = (b ^ c) % d
void crt_ln_mod_exp(LN* c, LN* d, LN* p, LN* q, LN* res, uint32_t digits);		 // res = (c ^ d) mod (p * q)
void ln_mod_exp_binary(LN* b, LN* a, LN* e, uint32_t edigits, LN* m, uint32_t mdigits);	// ʹ��ģ�ظ�ƽ����ģ���㷨
void ln_mod_exp_mary(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits);	// ʹ��m-ary��ģ���㷨

void mod_inv(LN* a, LN* b, LN* c, uint32_t digits);								 // a = (b ^ -1) mod c
bool ln_gcd(LN* b, LN* c, uint32_t digits);										 // �ж�b��c�Ƿ���

void print_ln(LN* a, uint32_t digits); // ��ӡ��������Ϊ0�������

void mod_inv_binary(LN* a, LN* b, LN* c, uint32_t digits);