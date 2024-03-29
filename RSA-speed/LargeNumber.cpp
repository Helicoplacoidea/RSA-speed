#include "LargeNumber.h"
#include "Montgomery.h"

LN ln_l_shift(LN* a, LN* b, uint32_t c, uint32_t digits)
{
	LN temp, carry;
	uint32_t cnt;
	if (c >= LN_BITS)
		return 0;
	carry = 0;
	cnt = LN_BITS - c;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i];
		a[i] = (temp << c) | carry;
		carry = c ? (temp >> cnt) : 0;
	}
	return carry;
}

LN ln_r_shift(LN* a, LN* b, uint32_t c, uint32_t digits)
{
	LN temp, carry;
	uint32_t cnt;
	if (c >= LN_BITS)
		return 0;
	carry = 0;
	cnt = LN_BITS - c;
	int i = digits - 1;
	for (; i >= 0; i--)
	{
		temp = b[i];
		a[i] = (temp >> c) | carry;
		carry = c ? (temp << cnt) : 0;
	}
	return carry;
}

void str_to_ln(uint8_t* str, LN* ln, uint32_t digits, uint32_t size)
{
	LN temp;
	uint32_t i, u;
	int j;
	for (i = 0, j = size - 1; i < digits && j >= 0; i++)
	{
		temp = 0;
		for (u = 0; j >= 0 && u < LN_BITS; j--, u += 8)
		{
			temp |= ((LN)str[j]) << u;
		}
		ln[i] = temp;
	}
	for (; i < digits; i++)
	{
		ln[i] = 0;
	}
}

void ln_to_str(uint8_t* str, LN* ln, uint32_t digits, uint32_t size)
{
	LN temp;
	uint32_t i, u;
	int j;
	for (i = 0, j = size - 1; i < digits && j >= 0; i++)
	{
		temp = ln[i];
		for (u = 0; j >= 0 && u < LN_BITS; j--, u += 8)
		{
			str[j] = (uint8_t)(temp >> u);
		}
	}
	for (; j >= 0; j--)
	{
		str[j] = 0;
	}
}

void ln_assign_zero(LN* a, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		a[i] = 0;
	}
}

void ln_assign(LN* a, LN* b, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		a[i] = b[i];
	}
}

void ln_assign_2exp(LN* a, uint32_t b, uint32_t digits)
{
	ln_assign_zero(a, digits);
	if (b >= (digits * LN_BITS))
	{
		return; // Out of range
	}

	a[b / LN_BITS] = (LN)1 << (b % LN_BITS);
}

static uint32_t ln_digit_bits(LN a)
{
	uint32_t i;
	for (i = 0; i < LN_BITS; i++)
	{
		if (a == 0)
			break;
		a >>= 1;
	}
	return i;
}

uint32_t ln_get_bits(LN* a, uint32_t digits)
{
	digits = ln_get_digits(a, digits);
	if (digits == 0)
		return 0;
	return (ln_digit_bits(a[digits - 1]) + LN_BITS * (digits - 1));
}

uint32_t ln_get_digits(LN* a, uint32_t digits)
{
	int i;
	uint32_t cnt = 0;
	for (i = digits - 1; i >= 0; i--)
	{
		if (a[i])
			break;
	}
	i++;
	return i;
}

int ln_cmp(LN* a, LN* b, uint32_t digits)
{
	for (int i = digits - 1; i >= 0; i--)
	{
		if (a[i] > b[i])
			return 1;
		if (a[i] < b[i])
			return -1;
	}
	return 0;
}

int ln_is_zero(LN* a, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		if (a[i])
		{
			return 0;
		}
	}
	return 1;
}

int ln_is_one(LN* a, uint32_t digits)
{
	if (a[0] != 1)
		return 0;
	for (uint32_t i = 1; i < digits; i++)
	{
		if (a[i])
		{
			return 0;
		}
	}
	return 1;
}

LN ln_add(LN* a, LN* b, LN* c, uint32_t digits)
{
	LN carry = 0;
	LN temp;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i] + carry;
		if (temp < carry)
		{
			// b+carry进位，只可能结果�?0
			temp = c[i];
		}
		else if ((temp += c[i]) < c[i])
		{
			// b+c+carry进位
			carry = 1;
		}
		else
		{
			carry = 0;
		}
		a[i] = temp;
	}
	return carry;
}

LN ln_sub(LN* a, LN* b, LN* c, uint32_t digits)
{
	LN borrow = 0;
	LN temp;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i] - borrow;
		if (temp > (LN_MAX - borrow))
		{
			temp = LN_MAX - c[i];
		}
		else if ((temp -= c[i]) > (LN_MAX - c[i]))
		{
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		a[i] = temp;
	}
	return borrow;
}

LN ln_sub_digit_mul(LN* a, LN* b, LN c, LN* d, uint32_t digits)
{
	dLN res;
	LN borrow, res_h, res_l;
	if (c == 0)
		return 0;

	borrow = 0;
	for (uint32_t i = 0; i < digits; i++)
	{
		res = d[i] * (dLN)c;
		res_h = (res >> LN_BITS) & LN_MAX;
		res_l = res & LN_MAX;
		a[i] = b[i] - borrow;
		if (a[i] > LN_MAX - borrow)
		{
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		if ((a[i] -= res_l) > (LN_MAX - res_l))
		{
			borrow++;
		}
		borrow += res_h;
	}
	return borrow;
}

void ln_mul(LN* c, LN* a, LN* b, uint32_t digits)
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
	ln_assign(c, t, 2 * LN_MAX_DIGITS);
}

void ln_div(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	dLN temp;
	LN ai, t, c_temp[2 * LN_MAX_DIGITS + 1], d_temp[LN_MAX_DIGITS];
	int i;
	uint32_t dddigits, shift;

	dddigits = ln_get_digits(d, ddigits);
	if (dddigits == 0)
		return;
	if (ln_is_zero(c, cdigits) == 1) {
		ln_assign_zero(a, cdigits);
		ln_assign_zero(b, dddigits);
		return;
	}

	shift = LN_BITS - ln_digit_bits(d[dddigits - 1]);
	ln_assign_zero(c_temp, dddigits);
	c_temp[cdigits] = ln_l_shift(c_temp, c, shift, cdigits);
	ln_l_shift(d_temp, d, shift, dddigits);
	t = d_temp[dddigits - 1];

	ln_assign_zero(a, cdigits);
	i = cdigits - dddigits;
	for (; i >= 0; i--)
	{
		if (t == LN_MAX)
		{
			ai = c_temp[i + dddigits];
		}
		else
		{
			temp = c_temp[i + dddigits - 1];
			temp += (dLN)c_temp[i + dddigits] << LN_BITS;
			ai = temp / (t + 1);
		}

		c_temp[i + dddigits] -= ln_sub_digit_mul(&c_temp[i], &c_temp[i], ai, d_temp, dddigits);
		while (c_temp[i + dddigits] || (ln_cmp(&c_temp[i], d_temp, dddigits) >= 0))
		{
			ai++;
			c_temp[i + dddigits] -= ln_sub(&c_temp[i], &c_temp[i], d_temp, dddigits);
		}
		a[i] = ai;
	}

	ln_assign_zero(b, ddigits);
	ln_r_shift(b, c_temp, shift, dddigits);
}

void ln_mod(LN* a, LN* b, uint32_t bdigits, LN* c, uint32_t cdigits)
{
	LN temp[2 * LN_MAX_DIGITS] = { 0 };
	ln_div(temp, a, b, bdigits, c, cdigits);
}

void ln_mod_mul(LN* a, LN* b, LN* c, LN* d, uint32_t digits)
{
	LN temp[2 * LN_MAX_DIGITS];
	ln_mul(temp, b, c, digits);
	ln_mod(a, temp, 2 * digits, d, digits);
}

void ln_mod_exp(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	LN b_exp[3][LN_MAX_DIGITS], temp[LN_MAX_DIGITS], ci;
	uint32_t ci_bits, s;
	ln_assign(b_exp[0], b, ddigits);
	ln_mod_mul(b_exp[1], b, b_exp[0], d, ddigits);
	ln_mod_mul(b_exp[2], b, b_exp[1], d, ddigits);

	ln_assign_zero(temp, ddigits);
	temp[0] = 1;

	cdigits = ln_get_digits(c, cdigits);
	int i = cdigits - 1;
	for (; i >= 0; i--)
	{
		ci = c[i];
		ci_bits = LN_BITS;
		if (i == (int)(cdigits - 1))
		{
			while (!shift_2_bits(ci))
			{
				ci <<= 2;
				ci_bits -= 2;
			}
		}
		for (uint32_t j = 0; j < ci_bits; j += 2)
		{
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			s = shift_2_bits(ci);
			if (s != 0)
			{
				ln_mod_mul(temp, temp, b_exp[s - 1], d, ddigits);
			}
			ci <<= 2;
		}
	}
	ln_assign(a, temp, ddigits);
}

void ln_mod_exp_binary(LN* b, LN* a, LN* e, uint32_t edigits, LN* m, uint32_t mdigits)
{
	LN temp[LN_MAX_DIGITS], p[LN_MAX_DIGITS], E[LN_MAX_DIGITS];
	ln_assign(temp, a, LN_MAX_DIGITS);

	int k = ln_get_bits(e, edigits);
	int ei, num, cnt, esh;


	for (int i = k - 2; i >= 0; i--)
	{
		num = i / LN_BITS;
		cnt = i - num * LN_BITS;
		ln_assign(E, e, LN_MAX_DIGITS);

		ln_mod_mul(p, temp, temp, m, mdigits);
		esh = E[num] >> (cnt);
		ei = esh & 1U;

		if (ei == 1)
		{
			ln_mod_mul(temp, p, a, m, mdigits);
		}
		else
		{
			ln_assign(temp, p, LN_MAX_DIGITS);
		}
		// ln_assign_zero(esh, LN_MAX_DIGITS);
	}
	ln_assign(b, temp, LN_MAX_DIGITS);
}

void ln_mod_exp_mary(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	LN b_exp[33][LN_MAX_DIGITS], temp[LN_MAX_DIGITS], ci;
	uint32_t ci_bits, s;
	ln_assign(b_exp[0], b, ddigits);
	for (int i = 1; i < 32; i++)
	{
		ln_mod_mul(b_exp[i], b, b_exp[i - 1], d, ddigits);
	}

	ln_assign_zero(temp, ddigits);
	temp[0] = 1;

	cdigits = ln_get_digits(c, cdigits);
	int i = cdigits - 1;
	for (; i >= 0; i--)
	{
		ci = c[i];
		ci_bits = LN_BITS;
		if (i == (int)(cdigits - 1))
		{
			while (!shift_4_bits(ci))
			{
				ci <<= 4;
				ci_bits -= 4;
			}
		}
		for (uint32_t j = 0; j < ci_bits; j += 4)
		{
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			s = shift_4_bits(ci);
			if (s != 0)
			{
				ln_mod_mul(temp, temp, b_exp[s - 1], d, ddigits);
			}
			ci <<= 4;
		}
	}
	ln_assign(a, temp, ddigits);
}

void crt_ln_mod_exp(LN* c, LN* d, LN* p, LN* q, LN* res, uint32_t digits)
{
	LN d1[LN_MAX_DIGITS], d2[LN_MAX_DIGITS], C[LN_MAX_DIGITS], D[LN_MAX_DIGITS], p1[LN_MAX_DIGITS], q1[LN_MAX_DIGITS];
	LN one[LN_MAX_DIGITS];
	LN m1[LN_MAX_DIGITS], m2[LN_MAX_DIGITS], m22[LN_MAX_DIGITS], m[2 * LN_MAX_DIGITS];
	LN P[LN_MAX_DIGITS], Q[LN_MAX_DIGITS];
	assign_one_digit(one, 1, digits);
	ln_assign(C, c, digits);
	ln_assign(D, d, digits);
	ln_assign(P, p, digits);
	ln_assign(Q, q, digits);
	ln_assign(p1, p, digits);
	ln_assign(q1, q, digits);
	//ln_sub(q1, q1, one, digits);
	//ln_sub(p1, p1, one, digits);
	p1[0] -= 1;
	q1[0] -= 1;
	ln_mod(d1, D, digits, p1, digits);
	ln_mod(d2, D, digits, q1, digits);

	//ln_mod_exp_mary(m1, C, d1, digits, P, digits);
	ModExp(m1, C, d1, P, digits);
	//ln_mod_exp_mary(m2, C, d2, digits, Q, digits);
	ModExp(m2, C, d2, Q, digits);

	LN temp[LN_MAX_DIGITS], p_1[LN_MAX_DIGITS];

	while (ln_cmp(m1, m2, digits) >= 0) {
		ln_add(m2, m2, q, digits);
	}
	ln_sub(temp, m2, m1, digits);

	mod_inv(p_1, P, Q, digits);

	ln_mod_mul(m22, temp, p_1, Q, digits);

	ln_mul(m, m22, P, digits);
	ln_add(m, m, m1, digits);
	ln_assign(res, m, digits);
}

// void mod_inv(LN *a, LN *b, LN *c, uint32_t digits)	//欧拉定理直接求模逆，效率过低
// {
// 	LN temp[LN_MAX_DIGITS], one[LN_MAX_DIGITS], res[LN_MAX_DIGITS];
// 	assign_one_digit(temp, 1, LN_MAX_DIGITS);
// 	assign_one_digit(one, 1, LN_MAX_DIGITS);
// 	while (ln_cmp(temp, c, LN_MAX_DIGITS) == -1)
// 	{
// 		ln_mod_mul(res, b, temp, c, LN_MAX_DIGITS);
// 		if (ln_cmp(res, one, LN_MAX_DIGITS) == 0)
// 		{
// 			ln_assign(a, temp, LN_MAX_DIGITS);

// 			memset((uint8_t *)one, 0, sizeof(one));
// 			memset((uint8_t *)temp, 0, sizeof(temp));
// 			memset((uint8_t *)res, 0, sizeof(res));
// 			return;
// 		}
// 		ln_add(temp, temp, one, LN_MAX_DIGITS);
// 	}
// }

void mod_inv(LN* a, LN* b, LN* c, uint32_t digits)
{ // a=b-1modc
	// int m = -1, n;
	LN m[LN_MAX_DIGITS], n[LN_MAX_DIGITS], mod[LN_MAX_DIGITS], N[LN_MAX_DIGITS], B[LN_MAX_DIGITS], C[LN_MAX_DIGITS];

	ln_assign(N, c, digits);
	ln_assign(C, c, digits);
	ln_assign(B, b, digits);

	assign_one_digit(m, 1, LN_MAX_DIGITS);
	// int temp;
	LN temp1[LN_MAX_DIGITS], temp2[2 * LN_MAX_DIGITS];
	// if (a < b) {
	// 	temp = a, a = b, b = temp;
	// }
	LN bi[2000][LN_MAX_DIGITS], ai[2000][LN_MAX_DIGITS];
	int cnt = 0;
	while (!ln_is_zero(m, digits))
	{

		// n = a / b;
		ln_div(n, mod, c, digits, b, digits);

		// temp = a, a = b;
		ln_assign(temp1, c, digits);
		ln_assign(c, b, digits);
		// m = temp - a * n;
		ln_mul(temp2, c, n, digits);
		ln_sub(m, temp1, temp2, digits);
		// b = m;
		ln_assign(b, m, digits);

		ln_assign(ai[cnt], n, digits);
		cnt++;
	}
	cnt--;
	assign_one_digit(bi[0], 1, digits);
	ln_assign(bi[1], ai[cnt - 1], digits);
	for (int i = 2; i <= cnt; i++)
	{
		ln_mul(bi[i], ai[cnt - i], bi[i - 1], digits);
		ln_add(bi[i], bi[i], bi[i - 2], digits);
	}
	//cout << "cnt:" << cnt << endl;
	if (cnt % 2 != 0)
	{
		ln_sub(a, N, bi[cnt], digits);
	}
	else
	{
		ln_assign(a, bi[cnt], digits);
	}
	ln_assign(b, B, digits);
	ln_assign(c, C, digits);
}

void mod_inv_binary(LN* a, LN* b, LN* c, uint32_t digits)
{
	LN N[LN_MAX_DIGITS], x[LN_MAX_DIGITS], q[LN_MAX_DIGITS], T[LN_MAX_DIGITS];
	LN u[LN_MAX_DIGITS], v[LN_MAX_DIGITS], d[LN_MAX_DIGITS];
	ln_assign(N, c, digits);
	ln_assign(x, b, digits);
	ln_div(q, T, N, digits, x, digits);
	ln_assign(N, x, digits);
	ln_assign(x, T, digits);
	if (ln_is_zero(x, digits)) {
		ln_assign_zero(u, digits);
		return;
	}
	LN  f[LN_MAX_DIGITS], one[LN_MAX_DIGITS];
	int k = 0;
	assign_one_digit(one, 1, digits);
	ln_assign_zero(f, digits);
	while (N[0] % 2 == 0 && x[0] % 2 == 0) {
		k++;
		ln_r_shift(N, N, 1, digits);
		ln_r_shift(x, x, 1, digits);
	}
	if (x[0] % 2 == 0) {
		ln_assign(T, x, digits);
		ln_assign(x, N, digits);
		ln_assign(N, T, digits);
		assign_one_digit(f, 1, digits);
	}
	LN UA[LN_MAX_DIGITS], UB[LN_MAX_DIGITS],A[LN_MAX_DIGITS],B[LN_MAX_DIGITS];
	LN v1[LN_MAX_DIGITS], t1[LN_MAX_DIGITS];
	assign_one_digit(UB, 1, digits);
	ln_assign(A, N, digits);
	ln_assign(B, x, digits);
	ln_assign(v1, x, digits);
	int t1_flag = 0, UA_flag = 0;
	if (N[0] % 2 == 1) {
		ln_assign_zero(UA, digits);
		t1_flag = -1;
		ln_assign(t1, x, digits);
	}
	else {
		ln_add(UA, x, one, digits);
		ln_r_shift(UA, UA, 1, digits);
		t1_flag = 1;
		ln_r_shift(t1, N, 1, digits);
	}

	while (!ln_is_zero(t1, digits)) {
		if (t1_flag == 1) {
			ln_assign(UB, UA, digits);
			ln_assign(A, t1, digits);
		}
		else {
			ln_sub(B, x, UA, digits);
			ln_assign(v1, t1, digits);
		}
		if (ln_cmp(UB, B, digits) == 1) {
			UA_flag = 1;
			ln_sub(UA, UB, B, digits);
		}
		else {
			UA_flag = -1;
			ln_sub(UA, B, UB, digits);
		}
		if (ln_cmp(A, v1, digits) == 1) {
			t1_flag = 1;
			ln_sub(t1, A, v1, digits);
		}
		else {
			t1_flag = -1;
			ln_sub(t1, v1, A, digits);
		}
		if (UA_flag == -1) {
			ln_sub(UA, x, UA, digits);
			UA_flag = 1;
		}
		while (t1[0] % 2 == 0 && !ln_is_zero(t1, digits)) {
			ln_r_shift(t1, t1, 1, digits);
			if (UA[0] % 2 == 0) {
				ln_r_shift(UA, UA, 1, digits);
			}
			else {
				ln_add(UA, UA, x, digits);
				ln_r_shift(UA, UA, 1, digits);
			}
		}
	}
	ln_assign(u, UB, digits);
	LN temp[2 * LN_MAX_DIGITS], temp1[LN_MAX_DIGITS];
	int v_flag = 0, u_flag = 0;
	ln_mul(temp, N, u, digits);
	if (ln_cmp(A, temp, digits) == 1) {
		v_flag = 1;
		ln_sub(v, A, temp, digits);
		ln_div(v, temp1, v, digits, x, digits);
		ln_assign(d, A, digits);
		ln_l_shift(d, d, k, digits);
	}
	else {
		v_flag = -1;
		ln_sub(v, temp, A, digits);
		ln_div(v, temp1, v, digits, x, digits);
		ln_assign(d, A, digits);
		ln_l_shift(d, d, k, digits);
	}

	if (ln_is_one(f, digits)) {
		ln_assign(T, u, digits);
		ln_assign(u, v, digits);
		u_flag = v_flag;
		ln_assign(v, T, digits);
		v_flag = 1;
	}

	if (v_flag == 1) {
		ln_mul(temp, v, q, digits);
		ln_sub(u, c, u, digits);
		ln_sub(u, u, temp, digits);
	}
	else {
		ln_mul(temp, v, q, digits);
		ln_add(u, u, temp, digits);
	}
	ln_assign(a, u, digits);
}

bool ln_gcd(LN* b, LN* c, uint32_t digits)
{
	// int m = -1, n;
	LN m[LN_MAX_DIGITS], n[LN_MAX_DIGITS], mod[LN_MAX_DIGITS], N[LN_MAX_DIGITS], B[LN_MAX_DIGITS], C[LN_MAX_DIGITS];

	ln_assign(N, c, digits);
	ln_assign(C, c, digits);
	ln_assign(B, b, digits);
	assign_one_digit(m, 1, LN_MAX_DIGITS);
	// int temp;
	LN temp1[LN_MAX_DIGITS], temp2[2 * LN_MAX_DIGITS];
	// if (a < b) {
	// 	temp = a, a = b, b = temp;
	// }
	int cnt = 0;
	while (!ln_is_zero(m, digits))
	{
		// n = a / b;
		ln_div(n, mod, c, digits, b, digits);

		// temp = a, a = b;
		ln_assign(temp1, c, digits);
		ln_assign(c, b, digits);
		// m = temp - a * n;
		ln_mul(temp2, c, n, digits);
		ln_sub(m, temp1, temp2, digits);
		// b = m;
		ln_assign(b, m, digits);
		cnt++;
	}
	if (ln_is_zero(m, digits) && ln_is_one(c, digits))
	{
		ln_assign(b, B, digits);
		ln_assign(c, C, digits);
		return true;
	}
	else
	{
		ln_assign(b, B, digits);
		ln_assign(c, C, digits);
		return false;
	}
}

void print_ln(LN* a, uint32_t digits)
{
	uint8_t str[512];
	memset(str, 0, 512);
	ln_to_str(str, a, digits, 512);
	int i = 0;
	while (str[i] == 0)
		i++;
	if (i == 512)
	{
		printf("0\n");
		return;
	}
	for (; i < 512; i++)
		printf("%02X", str[i]);
	printf("\n");
}