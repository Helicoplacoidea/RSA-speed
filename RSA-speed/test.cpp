#include "RSA.h"
#include "LargeNumber.h"
#include "PrimeGenaration.h"
#include "Montgomery.h"

int main()
{
	clock_t start, finish;
	double duration;

	LN a[LN_MAX_DIGITS], b[LN_MAX_DIGITS], c[LN_MAX_DIGITS], m[LN_MAX_DIGITS];
	uint32_t len1, len2;
	uint8_t str1[100] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZzhegeshifouyoudian";
	uint8_t str2[100] = "q";
	len1 = strlen((const char*)str1);
	len2 = strlen((const char*)str2);
	str_to_ln(str1, a, LN_MAX_DIGITS, len1);
	str_to_ln(str2, b, LN_MAX_DIGITS, len2);

	LN p[LN_MAX_DIGITS], q[LN_MAX_DIGITS], e[LN_MAX_DIGITS], w[2 * LN_MAX_DIGITS];
	// assign_one_digit(p, 99999721, LN_MAX_DIGITS);
	// assign_one_digit(q, 99999773, LN_MAX_DIGITS);
	int key_len = 0;
	cout << "选择生成素数的长度(1024/768/512 bit)：";
	cin >> key_len;
	initialize_rand();
	// prime_genarate_table(p);
	// prime_genarate_table(q);
	// finish = clock();
	// duration = (double)(finish - start) / CLOCKS_PER_SEC;
	// printf("计算素数所需时紒E%f seconds\n", duration);
	// start = clock();
	// prime_genarate(p);
	// prime_genarate(q);
	// //ln_mod_exp(p, q, q, LN_MAX_DIGITS, a, ln_get_digits(a,LN_MAX_DIGITS));
	// finish = clock();
	// duration = (double)(finish - start) / CLOCKS_PER_SEC;
	// printf("计算素数所需时紒E%f seconds\n", duration);
	start = clock();

		bit_array(p, key_len);
		bit_array(q, key_len);

	finish = clock();
	duration = (double)(finish - start);
	printf("计算素数所需时间%f ms\n", duration);
	// str_to_ln(key_p1, p, LN_MAX_DIGITS, 128);
	start = clock();
	// cout << "flag:" << miller_Rabin(p, 3, LN_MAX_DIGITS) << endl;
	cout << "flag:" << probable_prime(p, LN_MAX_DIGITS) << endl;
	finish = clock();
	duration = (double)(finish - start);
	printf("一次素性检测所需时间%f ms\n", duration);
	// str_to_ln(key_p2, q, LN_MAX_DIGITS, 128);
	cout << "flag:" << miller_Rabin(q, 3, LN_MAX_DIGITS) << endl;

	// LN p[LN_MAX_DIGITS];
	// bit_array(p);
	// print_ln(p, LN_MAX_DIGITS);
	assign_one_digit(e, 65537, LN_MAX_DIGITS);
	LN n[2 * LN_MAX_DIGITS];
	ln_mul(n, p, q, LN_MAX_DIGITS);
	LN pp[LN_MAX_DIGITS], qq[LN_MAX_DIGITS], one[LN_MAX_DIGITS], test[2 * LN_MAX_DIGITS];
	assign_one_digit(one, 1, LN_MAX_DIGITS);


	start = clock();
	//ln_mod_exp_mary(b, p, q, LN_MAX_DIGITS, q, LN_MAX_DIGITS);
	for (int i = 0; i < 10000; i++)
		ln_mul(test, n, q, LN_MAX_DIGITS);
	// for (int i = 0; i < 1000; i++)
	//MonMul(b, p, q, test, LN_MAX_DIGITS);
	// ln_shift_mod(b, n, LN_MAX_DIGITS, a, LN_MAX_DIGITS);
	finish = clock();
	duration = (double)(finish - start);
	printf("新方法所需时间%f ms\n", duration);

	ln_sub(pp, p, one, LN_MAX_DIGITS);
	ln_sub(qq, q, one, LN_MAX_DIGITS);
	ln_mul(w, pp, qq, LN_MAX_DIGITS);
	// print_ln(one, LN_MAX_DIGITS);
	cout << "w:";
	print_ln(w, LN_MAX_DIGITS);
	cout << "a:";
	print_ln(a, LN_MAX_DIGITS);
	cout << "e:";
	print_ln(e, LN_MAX_DIGITS);
	cout << "n:";
	print_ln(n, LN_MAX_DIGITS);

	if (ln_gcd(e, w, LN_MAX_DIGITS))
	{
		cout << "true" << endl;
	}
	else
	{
		cout << "false" << endl;
	}

	// str_to_ln(str2, d, LN_MAX_DIGITS, len2);
	LN d[LN_MAX_DIGITS];
	start = clock();
	for(int i=0;i<10000;i++)
	mod_inv(d, e, w, LN_MAX_DIGITS);
	finish = clock();
	duration = (double)(finish - start);
	printf("求模逆所需时间：%f ms\n", duration);
	// bn_mod_inv(d, e, w, LN_MAX_DIGITS);
	cout << "d:";
	print_ln(d, LN_MAX_DIGITS);

	//LN x[LN_MAX_DIGITS], y[LN_MAX_DIGITS], z[LN_MAX_DIGITS];
	//assign_one_digit(x, 67608, LN_MAX_DIGITS);
	//assign_one_digit(y, 830616, LN_MAX_DIGITS);
	start = clock();
	for (int i = 0; i < 10000; i++)
	mod_inv_binary(d, e, w, LN_MAX_DIGITS);
	finish = clock();
	duration = (double)(finish - start);
	printf("求模逆所需时间：%f ms\n", duration);
	// bn_mod_inv(d, e, w, LN_MAX_DIGITS);
	cout << "d:";
	print_ln(d, LN_MAX_DIGITS);

	start = clock();

		ln_mod_exp_mary(c, a, e, LN_MAX_DIGITS, n, LN_MAX_DIGITS);
		//ModExp(c, a, e, n, LN_MAX_DIGITS);

	finish = clock();
	duration = (double)(finish - start);
	// ln_sub(c,c,a,LN_MAX_DIGITS);
	cout << "c:";
	print_ln(c, LN_MAX_DIGITS);
	printf("加密所需时间：%f ms\n", duration);

	start = clock();
	for (int i = 0; i < 1; i++) {
		//ln_mod_exp_mary(m, c, d, LN_MAX_DIGITS, n, LN_MAX_DIGITS);
		//ln_mod_exp_binary(m, c, d, LN_MAX_DIGITS, n, LN_MAX_DIGITS);
		crt_ln_mod_exp(c, d, p, q, m, LN_MAX_DIGITS);
		//ModExp(m, c, d, n, LN_MAX_DIGITS);
	}
	finish = clock();
	duration = (double)(finish - start);
	cout << "m:";
	print_ln(m, LN_MAX_DIGITS);
	printf("解密所需时间：%f ms\n", duration);

	start = clock();
	// ln_mod_exp(m, c, d, LN_MAX_DIGITS, n, LN_MAX_DIGITS);

	//ln_assign_zero(te, LN_MAX_DIGITS);
	//te[0] = 3867236407;
	//te[1] = 364624;
	ModExp_SW(m, c, d, n, LN_MAX_DIGITS);

	finish = clock();
	duration = (double)(finish - start);

	cout << "m:";
	print_ln(m, LN_MAX_DIGITS);
	printf("解密所需时间：%f ms\n", duration);

	system("pause");
	return 0;
}