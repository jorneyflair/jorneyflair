/*******************************************************************************************
 
			Hash(BEGIN(Satoshi[2010]), END(Sunny[2012])) == Videlicet[2014] ++
   
 [Learn, Create, but do not Forge] Viz. http://www.opensource.org/licenses/mit-license.php
  
*******************************************************************************************/

#include "core.h"

using namespace std;

unsigned int *primes;
unsigned int *inverses;

unsigned int BitArray_Size =  1024*1024*16;
mpz_t  zPrimorial;

static unsigned int prime_limit = 148948141;
static unsigned int PrimeLimit = 8388608;
static unsigned int PrimorialEndPrime = 12;

volatile unsigned int FourChainsFoundCounter = 0;

unsigned long sqrtld(unsigned long N) {
    int                 b = 1;
    unsigned long       res,s;
    while(1<<b<N) b+= 1;
    res = 1<<(b/2 + 1);
    for(;;) {
        s = (N/res + res)/2;
        if(s>=res) return res;
        res = s;
    }
}
 
unsigned int * make_primes(unsigned int limit) {
    unsigned int      *primes;
    unsigned long       i,j;
    unsigned long       s = sqrtld(prime_limit);
    unsigned long       n = 0;
    bool *bit_array_sieve = (bool*)malloc((prime_limit + 1) * sizeof(bool));
    bit_array_sieve[0] = 0;
    bit_array_sieve[1] = 0;
    for(i=2; i<=prime_limit; i++) bit_array_sieve[i] = 1;
    j = 4;
    while(j<=prime_limit) {
        bit_array_sieve[j] = 0;
        j += 2;
    }
    for(i=3; i<=s; i+=2) {
        if(bit_array_sieve[i] == 1) {
            j = i * 3;
            while(j<=prime_limit) {
                bit_array_sieve[j] = 0;
                j += 2 * i;
            }
        }
    }
    for(i=2;i<=prime_limit;i++) if(bit_array_sieve[i]==1) n += 1;
    primes = (unsigned int*)malloc((n + 1) * sizeof(unsigned long));
    primes[0] = n;
    j = 1;
    for(i=2;i<=prime_limit;i++) if(bit_array_sieve[i]==1) {
        primes[j] = i;
        j++;
    }
    free(bit_array_sieve);
    return primes;
}


namespace Core
{
	/** Divisor bit_array_sieve for Prime Searching. **/
	std::vector<unsigned int> DIVISOR_SIEVE;
	
	void InitializePrimes()
	{
		printf("\nGenerating primes...\n");
		// generate prime table

		primes = make_primes(prime_limit);

		printf("\n%d primes generated\n", primes[0]);

		mpz_init(zPrimorial);

		mpz_set_ui(zPrimorial, 1);

		for (int i=1; i<PrimorialEndPrime; i++)
		{
			mpz_mul_ui(zPrimorial, zPrimorial, primes[i]);
		}

		printf("\nPrimorial:");
		printf("\n"); mpz_out_str(stdout, 10, zPrimorial); printf("\n");

		printf("\nLast Primorial Prime = %u\n", primes[PrimorialEndPrime-1]);

		printf("\nFirst Sieving Prime = %u\n", primes[PrimorialEndPrime]);


		int size = mpz_sizeinbase(zPrimorial,2);
		printf("\nPrimorial Size = %d-bit\n\n", size);

		inverses=(unsigned int *) malloc((PrimeLimit+1)*sizeof(unsigned int));
		memset(inverses, 0, (PrimeLimit+1) * sizeof(unsigned int));

		mpz_t zPrime, zInverse, zResult;

		mpz_init(zPrime);
		mpz_init(zInverse);
		mpz_init(zResult);

		for(unsigned int i=PrimorialEndPrime; i<=PrimeLimit; i++)
		{
			mpz_set_ui(zPrime, primes[i]);

			int	inv = mpz_invert(zResult, zPrimorial, zPrime);
			if (inv <= 0)
			{
				printf("\nNo Inverse for prime %u at position %u\n\n", zPrime, i);
				exit(0);
			}
			else
			{
				inverses[i]  = mpz_get_ui(zResult);
			}
		}


	}
	
	/** Convert Double to unsigned int Representative. Used for encoding / decoding prime difficulty from nBits. **/
	unsigned int SetBits(double nDiff)
	{
		unsigned int nBits = 10000000;
		nBits *= nDiff;
		
		return nBits;
	}

	/** Determines the difficulty of the Given Prime Number.
		Difficulty is represented as so V.X
		V is the whole number, or Cluster Size, X is a proportion
		of Fermat Remainder from last Composite Number [0 - 1] **/
	double GetPrimeDifficulty(CBigNum prime, int checks)
	{
		if(!PrimeCheck(prime, checks))
			return 0.0; ///difficulty of a composite number
			
		CBigNum lastPrime = prime;
		CBigNum next = prime + 2;
		unsigned int clusterSize = 1;
		
		///largest prime gap in cluster can be +12
		///this was determined by previously found clusters up to 17 primes
		for( next ; next <= lastPrime + 12; next += 2)
		{
			if(PrimeCheck(next, checks))
			{
				lastPrime = next;
				++clusterSize;
			}
		}
		
		///calulate the rarety of cluster from proportion of fermat remainder of last prime + 2
		///keep fractional remainder in bounds of [0, 1]
		double fractionalRemainder = 1000000.0 / GetFractionalDifficulty(next);
		if(fractionalRemainder > 1.0 || fractionalRemainder < 0.0)
			fractionalRemainder = 0.0;
		
		return (clusterSize + fractionalRemainder);
	}

	/** Gets the unsigned int representative of a decimal prime difficulty **/
	unsigned int GetPrimeBits(CBigNum prime, int checks)
	{
		return SetBits(GetPrimeDifficulty(prime, checks));
	}

	/** Breaks the remainder of last composite in Prime Cluster into an integer. 
		Larger numbers are more rare to find, so a proportion can be determined 
		to give decimal difficulty between whole number increases. **/
	unsigned int GetFractionalDifficulty(CBigNum composite)
	{
		/** Break the remainder of Fermat test to calculate fractional difficulty [Thanks Sunny] **/
		return ((composite - FermatTest(composite, 2) << 24) / composite).getuint();
	}
	
	
	/** bit_array_sieve of Eratosthenes for Divisor Tests. Used for Searching Primes. **/
	std::vector<unsigned int> Eratosthenes(int nSieveSize)
	{
		bool TABLE[nSieveSize];
		
		for(int nIndex = 0; nIndex < nSieveSize; nIndex++)
			TABLE[nIndex] = false;
			
			
		for(int nIndex = 2; nIndex < nSieveSize; nIndex++)
			for(int nComposite = 2; (nComposite * nIndex) < nSieveSize; nComposite++)
				TABLE[nComposite * nIndex] = true;
		
		
		std::vector<unsigned int> PRIMES;
		for(int nIndex = 2; nIndex < nSieveSize; nIndex++)
			if(!TABLE[nIndex])
				PRIMES.push_back(nIndex);

		
		printf("bit_array_sieve of Eratosthenes Generated %i Primes.\n", PRIMES.size());
		
		return PRIMES;
	}
	
	/** Basic Search filter to determine if further tests should be done. **/
	bool DivisorCheck(CBigNum test)
	{
		for(int index = 0; index < DIVISOR_SIEVE.size(); index++)
			if(test % DIVISOR_SIEVE[index] == 0)
				return false;
				
		return true;
	}

	/** Determines if given number is Prime. Accuracy can be determined by "checks". 
		The default checks the Coinshield Network uses is 2 **/
	bool PrimeCheck(CBigNum test, int checks)
	{
		/** Check B: Miller-Rabin Tests */
		bool millerRabin = Miller_Rabin(test, checks);
		if(!millerRabin)
			return false;
			
		/** Check C: Fermat Tests */
		for(CBigNum n = 2; n < 2 + checks; n++)
			if(FermatTest(test, n) != 1)
				return false;
		
		return true;
	}

	/** Simple Modular Exponential Equation a^(n - 1) % n == 1 or notated in Modular Arithmetic a^(n - 1) = 1 [mod n]. 
		a = Base or 2... 2 + checks, n is the Prime Test. Used after Miller-Rabin and Divisor tests to verify primality. **/
	CBigNum FermatTest(CBigNum n, CBigNum a)
	{
		CAutoBN_CTX pctx;
		CBigNum e = n - 1;
		CBigNum r;
		BN_mod_exp(&r, &a, &e, &n, pctx);
		
		return r;
	}

	/** Miller-Rabin Primality Test from the OpenSSL BN Library. **/
	bool Miller_Rabin(CBigNum n, int checks)
	{
		return (BN_is_prime(&n, checks, NULL, NULL, NULL) == 1);
	}

	static int Convert_BIGNUM_to_mpz_t(const BIGNUM *bn, mpz_t g)
	{
		bn_check_top(bn);
		if(((sizeof(bn->d[0]) * 8) == GMP_NUMB_BITS) &&
				(BN_BITS2 == GMP_NUMB_BITS)) 
			{
			/* The common case */
			if(!_mpz_realloc (g, bn->top))
				return 0;
			memcpy(&g->_mp_d[0], &bn->d[0], bn->top * sizeof(bn->d[0]));
			g->_mp_size = bn->top;
			if(bn->neg)
				g->_mp_size = -g->_mp_size;
			return 1;
			}
		else
			{
			char *tmpchar = BN_bn2hex(bn);
			if(!tmpchar) return 0;
			OPENSSL_free(tmpchar);
			return 0;
			}
	}


	unsigned long PrimeSieve(CBigNum BaseHash, unsigned int nDifficulty)
	{
		mpz_t zPrimeOrigin, zPrimeOriginOffset, zFirstSieveElement, zPrimorialMod, zTempVar, zResidue, zTwo, zN;

		unsigned int i = 0;
		unsigned int size = 0;
		unsigned long nonce = 0;

		mpz_init(zPrimeOriginOffset);
		mpz_init(zFirstSieveElement);
		mpz_init(zPrimorialMod);

		mpz_init(zTempVar);
		mpz_init(zPrimeOrigin);
		mpz_init(zResidue);
		mpz_init_set_ui(zTwo, 2);
		mpz_init(zN);


		unsigned char* bit_array_sieve = (unsigned char*)malloc((BitArray_Size)/8);
		memset(bit_array_sieve, 0x00, (BitArray_Size)/8);

		Convert_BIGNUM_to_mpz_t(&BaseHash, zPrimeOrigin);
		size = mpz_sizeinbase(zPrimeOrigin,2);

		mpz_mod(zPrimorialMod, zPrimeOrigin, zPrimorial);
		mpz_sub(zPrimorialMod, zPrimorial, zPrimorialMod);

		mpz_mod(zPrimorialMod, zPrimorialMod, zPrimorial);
		mpz_add_ui(zPrimorialMod, zPrimorialMod, 97);
		mpz_add(zTempVar, zPrimeOrigin, zPrimorialMod);

		mpz_set(zFirstSieveElement, zTempVar);

		for(unsigned int i=PrimorialEndPrime; i<PrimeLimit; i++)
		{
			unsigned long  p = primes[i];
			unsigned int inv = inverses[i];
			unsigned int base_remainder = mpz_tdiv_ui(zTempVar, p);

			unsigned int remainder = base_remainder;
			unsigned long r = (p-remainder)*inv;
			unsigned int index = r % p;
			while(index < BitArray_Size)
			{
				//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
					bit_array_sieve[(index)>>3] |= (1<<((index)&7));
				index += p;
			}
		
			remainder = base_remainder + 4;
			if (p<remainder)
				remainder -= p;
			r = (p-remainder)*inv;
			index = r % p;
			while(index < BitArray_Size)
			{
				//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
					bit_array_sieve[(index)>>3] |= (1<<((index)&7));
				index += p;
			}

			remainder = base_remainder + 6;
			if (p<remainder)
				remainder -= p;
			r = (p-remainder)*inv;
			index = r % p;
			while(index < BitArray_Size)
			{
				//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
					bit_array_sieve[(index)>>3] |= (1<<((index)&7));
				index += p;
			}
				
			remainder = base_remainder + 10;
			if (p<remainder)
				remainder -= p;
			r = (p-remainder)*inv;
			index = r % p;
			while(index < BitArray_Size)
			{
				//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
					bit_array_sieve[(index)>>3] |= (1<<((index)&7));
				index += p;
			}

			if ( nDifficulty > 50000000)
			{
				remainder = base_remainder + 12;
				if (p<remainder)
					remainder -= p;
				r = (p-remainder)*inv;
				index = r % p;
				while(index < BitArray_Size)
				{
					//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
						bit_array_sieve[(index)>>3] |= (1<<((index)&7));
					index += p;
				}
			}

			if ( nDifficulty > 60000000)
			{
				remainder = base_remainder + 16;
				if (p<remainder)
					remainder -= p;
				r = (p-remainder)*inv;
				index = r % p;
				while(index < BitArray_Size)
				{
					//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
						bit_array_sieve[(index)>>3] |= (1<<((index)&7));
					index += p;
				}
			}
		}

		for(i=0; i<BitArray_Size; i++)
		{
			if( bit_array_sieve[(i)>>3] & (1<<((i)&7)) )
				continue;

			// p1
			mpz_mul_ui(zTempVar, zPrimorial, i);
			mpz_add(zTempVar, zFirstSieveElement, zTempVar);

			mpz_set(zPrimeOriginOffset, zTempVar);
			mpz_sub_ui(zN, zTempVar, 1);
			mpz_powm(zResidue, zTwo, zN, zTempVar);
			if (mpz_cmp_ui(zResidue, 1) != 0)
				continue;

			// p2
			mpz_add_ui(zTempVar, zTempVar, 4);
			mpz_sub_ui(zN, zTempVar, 1);
			mpz_powm(zResidue, zTwo, zN, zTempVar);
			if (mpz_cmp_ui(zResidue, 1) != 0)
				continue;

			mpz_add_ui(zTempVar, zTempVar, 2);
			mpz_sub_ui(zN, zTempVar, 1);
			mpz_powm(zResidue, zTwo, zN, zTempVar);
			if (mpz_cmp_ui(zResidue, 1) != 0)
				continue;

			// p4
			mpz_add_ui(zTempVar, zTempVar, 4);
			mpz_sub_ui(zN, zTempVar, 1);
			mpz_powm(zResidue, zTwo, zN, zTempVar);
			if (mpz_cmp_ui(zResidue, 1) != 0)
				continue;
			else
			{
				FourChainsFoundCounter++;
				printf("\nTotal Number of 4-Chains Found = %u\n", FourChainsFoundCounter);
				
				// calculate offset
				mpz_sub(zTempVar, zPrimeOriginOffset, zPrimeOrigin);

				//size = mpz_sizeinbase(zTempVar,2);
				//printf("\nNonce size = %u-bit\n\n", size);

				nonce = mpz_get_ui(zTempVar);

				if(GetPrimeBits(BaseHash + nonce, 1) >= nDifficulty)
				{
					printf("\n\n******* BLOCK FOUND *******\n\n");
					free(bit_array_sieve);

					mpz_clear(zPrimeOrigin);
					mpz_clear(zPrimeOriginOffset);
					mpz_clear(zFirstSieveElement);
					mpz_clear(zResidue);
					mpz_clear(zTwo);
					mpz_clear(zN);
					mpz_clear(zPrimorialMod);
					mpz_clear(zTempVar);

					return nonce;
				}
			}
		}

		mpz_clear(zPrimeOrigin);
		mpz_clear(zPrimeOriginOffset);
		mpz_clear(zFirstSieveElement);
		mpz_clear(zResidue);
		mpz_clear(zTwo);
		mpz_clear(zN);
		mpz_clear(zPrimorialMod);
		mpz_clear(zTempVar);

		free(bit_array_sieve);
		return false;
	}
}

