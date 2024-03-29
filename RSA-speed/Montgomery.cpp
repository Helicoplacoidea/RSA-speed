#include "RSA.h"
#include "LargeNumber.h"
#include "Montgomery.h"


void mod_2_exp(LN* res, LN* a, uint32_t bits, uint32_t digits)
{
	LN tmodr[LN_MAX_DIGITS];
	int block = bits / LN_BITS;
	int cnt = bits - block * LN_BITS;

	ln_assign_zero(tmodr, digits);
	memcpy(tmodr, a, 4 * block);

	tmodr[block] = a[block] & (uint32_t)(pow(2, cnt) - 1);

	ln_assign(res, tmodr, digits);
}

void div_2_exp(LN* res, LN* a, uint32_t bits, uint32_t digits)
{
	LN tmodr[2 * LN_MAX_DIGITS];
	int block = bits / LN_BITS;
	int cnt = bits - block * LN_BITS;
	ln_assign_zero(tmodr, digits);

	int max = ln_get_digits(a, digits);

	for (int i = 0; i < max - block; i++) {
		tmodr[i] = a[i + block];
	}
	//memcpy(tmodr, &a[block], 4 * (max - block));

	LN shift[2*LN_MAX_DIGITS];
	ln_assign_zero(shift, 2*LN_MAX_DIGITS);
	ln_r_shift(shift, tmodr, cnt, digits);

	ln_assign(res, shift, LN_MAX_DIGITS);
}

void MonPro(LN* res, LN* a, LN* b, LN* r, LN* n, LN* n_1, uint32_t rbits, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS], tmodr[LN_MAX_DIGITS], m[2 * LN_MAX_DIGITS], temp[2 * LN_MAX_DIGITS], u[LN_MAX_DIGITS];

	ln_mul(t, a, b, digits);

	ln_assign_zero(tmodr, digits);

	mod_2_exp(tmodr, t, rbits, digits);

	ln_mul(m, tmodr, n_1, digits);
	mod_2_exp(m, m, rbits, digits);

	ln_mul(temp, m, n, digits);
	ln_add(temp, temp, t, 2 * digits);

	div_2_exp(u, temp, rbits, 2 * digits);

	if (ln_cmp(u, n, digits) >= 0)
	{
		ln_sub(res, u, n, digits);
	}
	else
	{
		ln_assign(res, u, digits);
	}
	print_ln(res, digits);
}

void MonPro_SOS(LN* res, LN* a, LN* b, LN* r, LN* n, LN* n_1, uint32_t rbits, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS];
	ln_assign_zero(t, 2 * digits);

	LN C = 0, S = 0;
	dLN CS = 0;
	int s1 = ln_get_digits(a, digits);
	int s2 = ln_get_digits(b, digits);
	for (int i = 0; i < s2; i++) {
		C = 0;
		for (int j = 0; j < s1; j++) {
			CS = t[i + j] + a[j] * (dLN)b[i] + C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			t[i + j] = S;
		}
		t[i + s1] = C;
	}


	LN m = 0;

	for (int i = 0; i < s1; i++) {
		C = 0;
		m = t[i] * n_1[0];		//mod 2^w?
		for (int j = 0; j < s1; j++) {
			CS = t[i + j] + (dLN)m * n[j] + C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			t[i + j] = S;
		}

		for (int j = i + s1; j < 2 * s1; j++) {
			CS = t[j] + (dLN)C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			//cout << C << " ";
			t[j] = S;
		}
	}
	t[2 * s1] = C;
	print_ln(t, 2 * digits);

	LN u[2 * LN_MAX_DIGITS], v[2 * LN_MAX_DIGITS];
	ln_assign_zero(u, 2 * digits);
	ln_assign_zero(v, 2 * digits);
	for (int j = 0; j <= s1; j++) {
		u[j] = t[j + s1];
	}
	//print_ln(u, 2 * digits);
	//LN B = 0, D = 0;
	//dLN BD = 0;
	//for (int j = 0; j <= s1; j++) {
	//	BD = u[j] - n[j] - B;
	//}
	if (ln_cmp(u, n, digits) >= 0)
	{
		ln_sub(res, u, n, digits);
	}
	else
	{
		ln_assign(res, u, digits);
	}

}

void MonPro_CIOS(LN* res, LN* a, LN* b, LN* r, LN* n, LN* n_1, uint32_t rbits, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS];
	ln_assign_zero(t, 2 * digits);

	LN C,S;
	dLN CS;
	int s = ln_get_digits(a, digits);
	LN m = 0;

	for (int i = 0; i <= s-1; i++) {
		C = 0;
		for (int j = 0; j <= s-1; j++) {
			CS = t[j] + a[j] * (dLN)b[i] + C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			t[j] = S;
		}
		CS = (dLN)t[s] + C;
		S = CS & LN_MAX;
		C = (CS >> LN_BITS) & LN_MAX;
		t[s] = S;
		t[s + 1] = C;
		C = 0;						//无影响

		m = t[0] * n_1[0];		//mod 2^w?

		CS = t[0] + m * (dLN)n[0];
		C = (CS >> LN_BITS) & LN_MAX;

		for (int j = 1; j <= s-1; j++) {
			CS = t[j] + (dLN)m * n[j] + C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			t[j - 1] = S;
		}
		CS = (dLN)t[s] + C;
		S = CS & LN_MAX;
		C = (CS >> LN_BITS) & LN_MAX;
		t[s-1] = S;
		t[s] = C + t[s + 1];
	}

	if (ln_cmp(t, n, digits) >= 0)
	{
		ln_sub(res, t, n, digits);
	}
	else
	{
		ln_assign(res, t, digits);
	}
	//print_ln(res, digits);
}

void MonPro_IFIOS(LN* res, LN* a, LN* b, LN* n, LN* n_1, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS];
	ln_assign_zero(t, 2 * digits);

	LN C, S, R0, R1;
	dLN CS, R0_S, R1_S;
	int s = ln_get_digits(a, digits);
	LN m = 0;
	for (int i = 0; i <= s - 1; i++) {
		R0_S = t[0] + (dLN)a[0] * b[i];
		S = R0_S & LN_MAX;
		R0 = (R0_S >> LN_BITS) & LN_MAX;
		m = S * n_1[0];
		R1_S = S + (dLN)m * n[0];
		S = R1_S & LN_MAX;
		R1 = (R1_S >> LN_BITS) & LN_MAX;
		for (int j = 1; j <= s - 1; j++) {
			R0_S = t[j] + (dLN)a[j] * b[i] + R0;
			S = R0_S & LN_MAX;
			R0 = (R0_S >> LN_BITS) & LN_MAX;
			R1_S = S + (dLN)m * n[j] + R1;
			S = R1_S & LN_MAX;
			R1 = (R1_S >> LN_BITS) & LN_MAX;
			t[j - 1] = S;
		}
		CS = (dLN)R0 + R1;
		S = CS & LN_MAX;
		C = (CS >> LN_BITS) & LN_MAX;
		t[s + 1] = t[s + 1] + C;
		CS = (dLN)t[s] + S;
		S = CS & LN_MAX;
		C = (CS >> LN_BITS) & LN_MAX;
		t[s - 1] = S;
		t[s] = t[s + 1] + C;
		t[s + 1] = 0;
	}
	if (ln_cmp(t, n, digits) >= 0)
	{
		ln_sub(res, t, n, digits);
	}
	else
	{
		ln_assign(res, t, digits);
	}
}

void MonMul(LN* res, LN* a, LN* b, LN* n, LN* n_1, LN* r, uint32_t digits)
{
	// // clock_t start,finish;
	// uint32_t bits = ln_get_bits(n, digits);
	// LN r[LN_MAX_DIGITS], a_1[LN_MAX_DIGITS], r_1[LN_MAX_DIGITS], n_1[LN_MAX_DIGITS], x[LN_MAX_DIGITS];
	// ln_assign_2exp(r, bits, digits);
	// // start = clock();
	// mod_inv(r_1, r, n, digits); // r_1
	// // finish = clock();
	// // double duration = finish - start;
	// // printf("模逆时间%f ms",duration);
	// LN temp[LN_MAX_DIGITS], mod[LN_MAX_DIGITS];
	// ln_mul(temp, r, r_1, digits);
	// ln_div(n_1, mod, temp, digits, n, digits);
	LN a_1[LN_MAX_DIGITS];
	ln_mod_mul(a_1, a, r, n, digits);
	MonPro(res, a_1, b, r, n, n_1, 0, digits);
	return;
}

void ModExp(LN* res, LN* M, LN* e, LN* n, uint32_t digits)
{
	uint32_t bits = ln_get_bits(n, digits);
	if (bits % LN_BITS != 0) {
		bits = (bits / LN_BITS + 1) * LN_BITS;
	}
	LN r[LN_MAX_DIGITS], r_1[LN_MAX_DIGITS], n_1[LN_MAX_DIGITS];
	ln_assign_2exp(r, bits, digits);

	mod_inv(r_1, r, n, digits); // r_1

	LN temp[2 * LN_MAX_DIGITS], mod[LN_MAX_DIGITS], one[LN_MAX_DIGITS];
	assign_one_digit(one, 1, digits);
	ln_mul(temp, r, r_1, digits);
	//ln_sub(temp, temp, one, 2 * digits);
	ln_div(n_1, mod, temp, 2 * digits, n, digits);
	LN M_1[2 * LN_MAX_DIGITS], x_1[LN_MAX_DIGITS];
	ln_assign_zero(M_1, 2 * LN_MAX_DIGITS);
	ln_mod_mul(M_1, M, r, n, digits);
	ln_mod(x_1, r, digits, n, digits);

	LN E[LN_MAX_DIGITS];
	ln_assign_zero(temp, 2 * LN_MAX_DIGITS);
	ln_assign(temp, M, LN_MAX_DIGITS);

	int k = ln_get_bits(e, digits);
	int ei, num, cnt, esh;
	ln_assign(E, e, LN_MAX_DIGITS);

	//int cnt2 = 0;
	for (int i = k - 1; i >= 0; i--)
	{
		num = i / LN_BITS;
		cnt = i - num * LN_BITS;

		//MonPro( , x_1, x_1, r, n, n_1, bits, digits);
		//MonPro_CIOS(x_1, x_1, x_1, r, n, n_1, bits, digits);
		MonPro_IFIOS(x_1, x_1, x_1, n, n_1, digits);
		//cnt2++;
		esh = E[num] >> (cnt);
		ei = esh & 1U;

		if (ei == 1)
		{
			// ln_mod_mul(temp, p, a, m, mdigits);
			//MonPro(x_1, M_1, x_1, r, n, n_1, bits, digits);
			//MonPro_CIOS(x_1, M_1, x_1, r, n, n_1, bits, digits);
			MonPro_IFIOS(x_1, M_1, x_1, n, n_1, digits);
			//cnt2++;
		}
	}
	//MonPro(res, x_1, one, r, n, n_1, bits, digits);
	//MonPro_CIOS(res, x_1, one, r, n, n_1, bits, digits);
	MonPro_IFIOS(res, x_1, one, n, n_1, digits);
	//cnt2++;
	//cout << "cnt:" << cnt2 << endl;
}

void ModExp_SW(LN* res, LN* M, LN* e, LN* n, uint32_t digits)
{
	uint32_t bits = ln_get_bits(n, digits);
	if (bits % LN_BITS != 0) {
		bits = (bits / LN_BITS + 1) * LN_BITS;
	}
	LN r[LN_MAX_DIGITS], r_1[LN_MAX_DIGITS], n_1[LN_MAX_DIGITS];
	ln_assign_2exp(r, bits, digits);

	mod_inv(r_1, r, n, digits); // r_1

	LN temp0[2 * LN_MAX_DIGITS], mod[LN_MAX_DIGITS], one[LN_MAX_DIGITS];
	assign_one_digit(one, 1, digits);
	ln_mul(temp0, r, r_1, digits);
	//ln_sub(temp, temp, one, 2 * digits);
	ln_div(n_1, mod, temp0, 2 * digits, n, digits);

	LN r2[LN_MAX_DIGITS],x_1[LN_MAX_DIGITS],temp[16][LN_MAX_DIGITS],m2[LN_MAX_DIGITS];
	ln_mod(r, r, digits, n, digits);
	ln_mod_mul(r2, r, r, n, digits);
	ln_assign(x_1, r, digits);
	MonPro_IFIOS(temp[1], M, r2, n, n_1, digits);

	MonPro_IFIOS(m2, temp[1], temp[1], n, n_1, digits);
	for (int i = 3; i <= 15; i += 2) {
		MonPro_IFIOS(temp[i], temp[i - 2], m2, n, n_1, digits);
	}

	int k = ln_get_bits(e, digits);
	int ei, num, cnt, esh, t;
	int E[LN_MAX_DIGITS * LN_BITS];

	for (int i = 0; i < k; i++) {
		num = i / LN_BITS;
		cnt = i - num * LN_BITS;
		esh = e[num] >> (cnt);
		ei = esh & 1U;
		E[i] = ei;
		//cout << E[i] << " ";
	}

	int i = 0;
	stack<int> st;
	while (i <= k - 1) {
		if (E[i] == 1) {
			if (i + 3 <= k - 1) {
				t = E[i] + 2 * E[i + 1] + 4 * E[i + 2] + 8 * E[i + 3];
				E[i + 1] = 0;
				E[i + 2] = 0;
				E[i + 3] = 0;
				i = i + 4;
			}
			else {
				if (i + 2 == k - 1) {
					t = E[i] + 2 * E[i + 1] + 4 * E[i + 2];
					E[i + 1] = 0;
					E[i + 2] = 0;
					i = k;
				}
				else if (i + 1 == k - 1) {
					t = E[i] + 2 * E[i + 1];
					E[i + 1] = 0;
					i = k;
				}
				else {
					t = E[i];
					i = k;
				}
			}
			st.push(t);
		}
		else {
			i++;
		}
	}

	//int cnt2 = 9;
	for (int i = k - 1; i >= 0; i--) {
		MonPro_IFIOS(x_1, x_1, x_1, n, n_1, digits);
		//cnt2++;
		if (E[i] == 1) {
			MonPro_IFIOS(x_1, x_1, temp[st.top()], n, n_1, digits);
			//cnt2++;
			st.pop();
		}
	}
	MonPro_IFIOS(res, x_1, one, n, n_1, digits);
	//st.empty();
	//cnt2++;
	//cout << "cnt:" << cnt2 << endl;
}