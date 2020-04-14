
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <assert.h>

// #include "asmuthbloom.h"

void
generatePrimes(uint32_t *primes, int num, gmp_randstate_t state) {
	printf("generatePrimes\n");
	int bc;
	bc = 9;
	mpz_t r;
	mpz_t pr;
	mpz_t one;
	mpz_init(r);
	mpz_init(pr);
	// mpz_init_set_ui(one, 1);
	int i;
	for (i = 0; i < num; i++) {
		printf("generatePrimes3\n");
		uint32_t x = gmp_urandomb_ui(state, bc);
		printf("generatePrimes4\n");
		mpz_set_ui(r, x);
		mpz_nextprime(pr, r);
		printf("%lu prime\n", mpz_get_ui(pr));
		primes[i] = mpz_get_ui(pr);

		// mpz_add(bc, bc, one);
		bc = bc + 1;
	}
	printf("generatePrimes2\n");
	mpz_clear(r);
	printf("generatePrimesr\n");
	mpz_clear(pr);
	printf("generatePrimespr\n");
	// mpz_clear(one);
	printf("generatePrimes5\n");
}

// mpz_t
// **pwrSet(mpz_t *arr, int len){
// 	if (len <= 0) {
// 		mpz_t **temp = {{}};
// 		return temp;
// 	}

// 	mpz_t x;
// 	mpz_init_set(x, arr[0]);
// 	mpz_t newArr[len - 1];
// 	int i = 1;
// 	for (i; i < len; i++){
// 		newArr[i-1] = arr[i];
// 	}
// 	len = len - 1;
// 	mpz_t **wo = pwrSet(newArr, len);
// 	mpz_t **w = {{}};
// 	int wi = 0;
// 	i = 0;
// 	for (i; i < (len * len); i++) {
// 		mpz_t *temp;
// 		int j = 0;
// 		for (j; j < len; j++) {
// 			mpz_set(temp[j], wo[i][j]);
// 		}
// 		mpz_set(temp[j+1], x);
// 		w[wi] = temp;
// 		wi++;
// 	}
// 	mpz_t **all;
// 	i = 0;
// 	for (i; i < (len) * (len); i++) {
// 		all[i] = wo[i];
// 	}
// 	i++;
// 	int j = 0;
// 	for (j; j < (len) * (len); j++){
// 		all[j+i] = w[j];
// 	}
// 	mpz_clear(x);
// 	return all;

// }

uint32_t
**powerSet(uint32_t *arr, size_t len, size_t olen) {
	printf("powerSet\n");
	if (len <= 0) {
		uint32_t **temp = malloc(sizeof(uint32_t*));
		temp[0] = calloc(olen, sizeof(uint32_t));
		return temp;
	}
	printf("powerSet2\n");

	uint32_t x;
	x = arr[len];
	uint32_t *newArr = malloc(sizeof(uint32_t) * (len - 1));
	printf("powerSet7\n");
	memcpy(newArr, arr, sizeof(uint32_t) * (len-1));
	int i = 1;
	// for (i; i < len; i++){
	// 	memcpy(newArr[i-1], arr[i], sizeof(uint32_t));
	// }
	printf("powerSet3\n");
	uint32_t **wo = malloc((len * len / 2) * sizeof(uint32_t*));
	uint32_t **w = malloc((len * len / 2) * sizeof(uint32_t*));
	uint32_t **temp = powerSet(newArr, len-1, olen);
	for (i = 0; i < (len * len / 2); i++) {
		wo[i] = malloc(sizeof(uint32_t) * olen);
		w[i] = malloc(sizeof(uint32_t) * olen);
		memcpy(wo[i], temp[i], sizeof(uint32_t) * olen);
		memcpy(w[i], temp[i], sizeof(uint32_t) * olen);
	}
	printf("powerSet4\n");
	for (i = 0; i < (len * len / 2); i++) {
		int j;
		for (j = 0; j < olen; j++) {
			if (w[i][j] == 0) {
				w[i][j] = x;
			}
		}
	}
	printf("powerSet5\n");
	uint32_t **all = malloc((len * len) * sizeof(uint32_t*));
	i = 0;
	int j = 0;
	while (i < (len*len)) {
		all[i] = malloc(sizeof(uint32_t) * olen);
		printf("powerSet8\n");
		memcpy(all[i], wo[j], sizeof(uint32_t) * olen);
		i++;
		printf("powerSet9\n");
		all[i] = malloc(sizeof(uint32_t) * olen);
		memcpy(all[i], w[j], sizeof(uint32_t) * olen);
		printf("powerSet10\n");
		i++;
		j++;
	}
	printf("powerSet6\n");

	for (i = 0; i < (len * len / 2); i++) {
		free(wo[i]);
		free(w[i]);
	}
	free(wo);
	free(w);
	free(newArr);

	return all;
}

int
findIndex(uint32_t *arr, uint32_t t, int len) {
	if (len == 0) {
		return -1;
	}
	int i = 0;
	while (i < len) {
		if(arr[i] == t) {
			return i;
		}
		else {
			i++;
		}
	}
	return -1;
}

// mpz_t
// *splitSecret(int secret, int n, int t, int N, int *w) {
// 	mpz_t s;
// 	mpz_init_set_ui(s, secret);
// 	mpz_t p;
// 	mpz_init_set_ui(p, 256);
// 	gmp_randstate_t state;
// 	gmp_randinit_default(state);
// 	mpz_t d[N] = generatePrimes(N, state);

// 	mpz_t s2, s3;
// 	mpz_init(s2);
// 	mpz_init(s3);
// 	mpz_mul(s2, d[5], d[3]);
// 	mpz_mul(s3, d[1], d[2]);

// 	mpz_t sequence[n] = {d[4], d[0], s2, s3};
// 	mpz_t **powSet = pwrSet(sequence, n);

// 	mpz_t leftPart, rightPart;
// 	mpz_init(leftPart);
// 	mpz_init(rightPart);
// 	int leftFirst = 1;
// 	int rightFirst = 1;

// 	int i = 0;
// 	for(i; i < (n * n); i++) {
// 		mpz_t *temp = powSet[i];
// 		int weightTotal = 0;
// 		mpz_t numTotal;
// 		mpz_init_set_ui(numTotal, 1);
// 		int j = 0;
// 		for (j; j < temp.length; j++) {
// 			mpz_mul(numTotal, numTotal, temp[j]);
// 			int index = findIndex(sequence, temp[j], n);
// 			weightTotal += w[index];
// 		}
// 		if (weightTotal < t) {
// 			if (rightFirst == 1) {
// 				mpz_set(rightPart, numTotal);
// 				rightFirst = 0;
// 			}
// 			else {
// 				if (mpz_cmp(rightPart, numTotal) < 0) {
// 					mpz_set(rightPart, numTotal);
// 				}
// 			}
// 		}
// 		else {
// 			if (leftFirst == 1) {
// 				mpz_set(leftPart, numTotal);
// 				leftFirst = 0;
// 			}
// 			else {
// 				if (mpz_cmp(leftPart, numTotal) > 0) {
// 					mpz_set(leftPart, numTotal);
// 				}
// 			}
// 		}
// 		mpz_clear(numTotal);
// 	}

// 	mpz_t range;
// 	mpz_init(range);
// 	mpz_sub(range, leftPart, rightPart);
// 	size_t bc = mpz_sizeinbase(range, 2);
// 	mpz_t r;
// 	mpz_init(r);
// 	mpz_urandomb(r, state, bc-1);
// 	mpz_add(r, r, rightPart);
// 	mpz_t sDash;
// 	mpz_init(sDash);
// 	mpz_t rm;
// 	mpz_init(rm);
// 	mpz_mod(rm, r, p);
// 	mpz_t diff, ss, sa;
// 	mpz_init(diff);
// 	mpz_init(ss);
// 	mpz_init(sa);
// 	mpz_sub(ss, s, rm);
// 	mpz_add(sa, ss, p);
// 	mpz_mod(diff, sa, p);

// 	mpz_add(sDash, r, diff);

// 	mpz_t parts[n];
// 	i = 0;
// 	for (i; i < n; i++) {
// 		mpz_t k;
// 		mpz_init(k);
// 		mpz_mod(k, sDash, sequence[i]);
// 		mpz_set(parts[i],k);
// 		mpz_clear(k);
// 	}

// 	mpz_clears(s, p, sDash, s2, s3, leftPart, rightPart, range, r, rm, diff, ss, sa);

// 	return parts;
// }

uint32_t
*splitSecret(uint32_t secret, int n, int t, int N, int *w) {
	printf("splitSecret\n");
	mpz_t s;
	mpz_init_set_ui(s, secret);
	mpz_t p;
	mpz_init_set_ui(p, 256);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	uint32_t *d = malloc(sizeof(uint32_t) * N);
	generatePrimes(d, N, state);

	uint32_t *sequence = malloc(sizeof(uint32_t) * n);
	sequence[0] = d[4];
	sequence[1] = d[0];
	sequence[2] = d[5] * d[3];
	sequence[3] = d[1] * d[2];
	int i;

	uint32_t **powSet = malloc(sizeof(uint32_t*) * (n*n));
	for (i = 0; i < (n*n); i++) {
		powSet[i] = malloc(sizeof(uint32_t) * n);
	}
	printf("here\n");
	memcpy(powSet, powerSet(sequence, n, n), sizeof(uint32_t*) * (n*n));

	uint32_t leftPart, rightPart;
	leftPart = 0;
	rightPart = 0;
	int leftFirst = 1;
	int rightFirst = 1;

	for(i = 0; i < (n * n); i++) {
		int weightTotal = 0;
		uint32_t numTotal = 1;
		int j = 0;
		while (powSet[i][j] != 0) {
			numTotal = numTotal * powSet[i][j];
			int index = findIndex(sequence, powSet[i][j], n);
			weightTotal += w[index];
			j++;
		}
		if (weightTotal < t) {
			if (rightFirst == 1) {
				rightPart = numTotal;
				rightFirst = 0;
			}
			else {
				if (rightPart < numTotal) {
					rightPart = numTotal;
				}
			}
		}
		else {
			if (leftFirst == 1) {
				leftPart = numTotal;
				leftFirst = 0;
			}
			else {
				if (leftPart > numTotal) {
					leftPart = numTotal;
				}
			}
		}
	}

	mpz_t range, lp, rp;
	mpz_inits(range, lp, rp);
	mpz_init_set_ui(lp, leftPart);
	mpz_init_set_ui(rp, rightPart);
	mpz_sub(range, lp, rp);
	size_t bc = mpz_sizeinbase(range, 2);
	mpz_t r;
	mpz_init(r);
	mpz_urandomb(r, state, bc-1);
	mpz_add(r, r, rp);
	mpz_t sDash;
	mpz_init(sDash);
	mpz_t rm;
	mpz_init(rm);
	mpz_mod(rm, r, p);
	mpz_t diff, ss, sa;
	mpz_init(diff);
	mpz_init(ss);
	mpz_init(sa);
	mpz_sub(ss, s, rm);
	mpz_add(sa, ss, p);
	mpz_mod(diff, sa, p);

	mpz_add(sDash, r, diff);

	uint32_t parts[n];
	for (i = 0; i < n; i+2) {
		mpz_t k, seq;
		mpz_inits(k, seq);
		mpz_set_ui(seq, sequence[i]);
		mpz_mod(k, sDash, seq);
		parts[i] = mpz_get_ui(seq);
		parts[i+1] = mpz_get_ui(k);
		mpz_clears(k, seq);
	}

	free(sequence);
	free(d);
	mpz_clears(s, p, sDash, lp, rp, range, r, rm, diff, ss, sa);

	printf("splitSecret2\n");

	return parts;
}

// mpz_t
// crt(mpz_t *remainders, mpz_t *modules, int len) {
// 	mpz_t module;
// 	mpz_init_set_ui(module, 1);
// 	int i = 0;
// 	for (i; i < len; i++){
// 		mpz_mul(module, module, modules[i]);
// 	}

// 	mpz_t result;
// 	mpz_init(result);
// 	i = 0;
// 	for (i; i < len; i++){
// 		mpz_t temp, ti, rmt, rmti, ra;
// 		mpz_inits(temp, ti, rmt, rmti, ra);
// 		mpz_tdiv_q(temp, module, modules[i]);
// 		mpz_invert(ti, temp, modules[i]);
// 		mpz_mul(rmt, remainders[i], temp);
// 		mpz_mul(rmti, rmt, ti);
// 		mpz_add(ra, result, rmti);
// 		mpz_mod(result, ra, module);
// 		mpz_clears(temp, ti, rmt, rmti, ra);
// 	}
// 	return result;
// }

int
recoverSecret(uint32_t *parts, int len) {
	printf("recoverSecret\n");
	mpz_t p;
	mpz_init_set_ui(p, 256);
	uint32_t remainders[len/2];
	uint32_t modules[len/2];
	int index = 0;
	int i;
	for (i = 0; i < len; i+2) {
		modules[index] = parts[i];
		remainders[index] = parts[i+1];
		index++;
	}

	mpz_t sDash, s, crt, tempm;
	mpz_inits(sDash, s, crt, tempm);

	mpz_t module;
	mpz_init_set_ui(module, 1);
	for (i = 0; i < len/2; i++){
		mpz_set_ui(tempm, modules[i]);
		mpz_mul(module, module, tempm);
	}

	mpz_t result;
	mpz_init(result);
	for (i = 0; i < len/2; i++){
		mpz_t temp, tempr, ti, rmt, rmti, ra;
		mpz_inits(temp, tempr, ti, rmt, rmti, ra);
		mpz_set_ui(tempm, modules[i]);
		mpz_tdiv_q(temp, module, tempm);
		mpz_invert(ti, temp, tempm);
		mpz_set_ui(tempr, remainders[i]);
		mpz_mul(rmt, tempr, temp);
		mpz_mul(rmti, rmt, ti);
		mpz_add(ra, result, rmti);
		mpz_mod(result, ra, module);
		mpz_clears(temp, tempr, ti, rmt, rmti, ra);
	}
	mpz_set(crt, result);
	mpz_clear(result);

	mpz_set(sDash, crt);
	mpz_mod(s, sDash, p);
	uint32_t secret = mpz_get_ui(s);
	mpz_clears(p, sDash, s, tempm);
	return secret;
}

int
main() {

	int count = 0;
	int k;
	for (k = 0; k < 256; k++) {
		uint32_t secret = k;
		int n = 4;
		int t = 4;
		int N = 6;
		int w[4] = {1,1,2,2};
		uint32_t *parts = splitSecret(secret, n, t, N, w);

		uint32_t array[4] = {parts[0],parts[1],parts[4],parts[5]};
		int recoveredSecret = recoverSecret(array, 4);

		if (k == recoveredSecret) {
			count++;
		}
	}


	printf("%d/256\n", count);
	return 0;
}

//gmp: overflow in mpz type
//Aborted (core dumped)