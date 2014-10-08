/*******************************************************************************************
 
			Hash(BEGIN(Satoshi[2010]), END(Sunny[2012])) == Videlicet[2014] ++
   
 [Learn, Create, but do not Forge] Viz. http://www.opensource.org/licenses/mit-license.php
  
*******************************************************************************************/

#include "core.h"

using namespace std;

unsigned int *primes;
unsigned int *inverses;

unsigned int nBitArray_Size =  1024*1024*16;
mpz_t  zPrimorial;

static unsigned int prime_limit = 71378571;
static unsigned int nPrimeLimit = 4194304;
static unsigned int nPrimorialEndPrime = 12;

unsigned int nFourChainsFoundCounter = 0;
unsigned int nFiveChainsFoundCounter = 0;
unsigned int nSixChainsFoundCounter = 0;
unsigned int nSevenChainsFoundCounter = 0;
unsigned int nEightChainsFoundCounter = 0;

extern unsigned int nBestHeight;

unsigned long octuplet_origins[256] = {15760091,25658441,93625991,182403491,226449521,661972301,910935911,1042090781,1071322781,1170221861,1394025161,1459270271,1712750771,
										1742638811,1935587651,2048038451,2397437501,2799645461,2843348351,3734403131,4090833821,5349522791,5379039551,5522988461,5794564661,
										5950513181,6070429481,6138646511,6193303001,6394117181,6520678511,6765896981,6969026411,7219975571,7602979451,8247812381,8750853101,
										9870884321,9966184841,11076719651,11234903411,11567910701,11881131791,12753314921,12848960471,12850665671,12886759001,13345214411,
										13421076281,15065117141,15821203241,16206106991,16427277941,16804790531,17140322651,17383048211,18234075311,18379278761,18821455181,
										18856092371,21276989801,21315831611,21803245811,22190786531,22367332061,22642418411,22784826131,22827253901,23393094071,24816950771,
										24887046251,24930296381,26092031081,28657304561,28900195391,29055481571,29906747861,30332927741,30526543121,31282661141,31437430091,
										31447680611,31779849371,31907755331,33081664151,33734375021,35035293101,35969034371,36551720741,37000821701,37037888801,37654490171,
										38855298941,40743911051,41614070411,43298074271,43813839521,44676352991,45549998561,46961199401,47346763811,48333938111,49788942011,
										49827604901,50144961941,50878435451,53001578081,54270148391,57440594201,60239937671,62184803951,63370318001,64202502431,65227645781,
										65409385031,66449431661,69707273171,71750241371,73457668631,74082349331,74445418121,74760009671,75161088461,75778477121,76289638961,
										77310104141,77653734071,78065091101,78525462131,79011826961,79863776801,79976720891,80041993301,80587471031,80790462281,82455937631,
										83122625471,84748266131,84882447101,85544974631,86408384591,87072248561,88163200661,88436579501,88815669401,89597692181,90103909781,
										91192669481,93288681371,93434383571,93487652171,93703549391,94943708591,95109448781,95391400451,96133393241,97249028951,98257943081,
										100196170421,101698684931,104487717401,105510861341,106506834431,107086217081,109750518791,110327129441,111422173391,114994357391,
										116632573901,117762315941,118025332961,119063726051,121317512201,123019590761,123775576271,124168028051,130683361421,131045869301,
										131176761251,131484693071,132595345691,133391614241,135614688941,138478375151,139017478331,139858746941,141763537451,143258671091,
										144224334251,147215781521,147332222951,148124799281,148323246341,148671287111,148719488831,148916953301,148949723381,150613299911,
										153779378561,155130467951,155521458551,156146394401,156456267881,157272782741,157519407581,163899228791,164138756051,165040690931,
										165792941381,165952761041,166004527301,166225007561,168626248781,169349651741,170316751721,170552481551,170587733201,170832928151,
										171681030791,172892544941,173405293331,174073117061,177620195561,178242755681,180180782051,180237252311,184430028311,185515423391,
										185814366581,186122739611,187735172741,187971393341,188090847011,189066712181,190192014821,192380171981,193725710021,194875423271,
										198006027671,198146724311,198658763111,198869317721,199658510321,199847262731,200599766441,201708760061,202506276431,203499800501,
										204503641871,206150764271,207369666851,208403006081,211925962091,214556015741,218389714001,218732226521};

inline int64 GetTimeMicros() 
 { 
     return (boost::posix_time::ptime(boost::posix_time::microsec_clock::universal_time()) - boost::posix_time::ptime(boost::gregorian::date(1970,1,1))).total_microseconds(); 
 } 


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

		for (int i=1; i<nPrimorialEndPrime; i++)
		{
			mpz_mul_ui(zPrimorial, zPrimorial, primes[i]);
		}

		printf("\nPrimorial:");
		printf("\n"); mpz_out_str(stdout, 10, zPrimorial); printf("\n");

		printf("\nLast Primorial Prime = %u\n", primes[nPrimorialEndPrime-1]);

		printf("\nFirst Sieving Prime = %u\n", primes[nPrimorialEndPrime]);


		int nSize = mpz_sizeinbase(zPrimorial,2);
		printf("\nPrimorial Size = %d-bit\n\n", nSize);

		inverses=(unsigned int *) malloc((nPrimeLimit+1)*sizeof(unsigned int));
		memset(inverses, 0, (nPrimeLimit+1) * sizeof(unsigned int));

		mpz_t zPrime, zInverse, zResult;

		mpz_init(zPrime);
		mpz_init(zInverse);
		mpz_init(zResult);

		for(unsigned int i=nPrimorialEndPrime; i<=nPrimeLimit; i++)
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

	double GetPrimeDifficulty2(CBigNum next, unsigned int clusterSize)
	{
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
		/** Check C: Fermat Tests */
		CBigNum n = 3;
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


	unsigned long PrimeSieve(CBigNum BaseHash, unsigned int nDifficulty, unsigned int nHeight)
	{
		mpz_t zPrimeOrigin, zPrimeOriginOffset, zFirstSieveElement, zPrimorialMod, zTempVar, zResidue, zTwo, zN;

		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int nSize = 0;
		unsigned long nNonce = 0;
		unsigned int nPrimeCount = 0;
		unsigned int nSieveDifficulty = 0;
		unsigned long nStart = 0;
		unsigned long nStop = 0;
		unsigned int nLastOffset = 0;

		long nElapsedTime = 0;
		long nStartTime = 0;

		mpz_init(zPrimeOriginOffset);
		mpz_init(zFirstSieveElement);
		mpz_init(zPrimorialMod);

		mpz_init(zTempVar);
		mpz_init(zPrimeOrigin);
		mpz_init(zResidue);
		mpz_init_set_ui(zTwo, 2);
		mpz_init(zN);

		Convert_BIGNUM_to_mpz_t(&BaseHash, zPrimeOrigin);
		nSize = mpz_sizeinbase(zPrimeOrigin,2);

		unsigned char* bit_array_sieve = (unsigned char*)malloc((nBitArray_Size)/8);


		for(j=0; j<256; j++)
		{
			nStartTime = GetTimeMicros();
			memset(bit_array_sieve, 0x00, (nBitArray_Size)/8);

			mpz_mod(zPrimorialMod, zPrimeOrigin, zPrimorial);
			mpz_sub(zPrimorialMod, zPrimorial, zPrimorialMod);

			mpz_mod(zPrimorialMod, zPrimorialMod, zPrimorial);
			mpz_add_ui(zPrimorialMod, zPrimorialMod, octuplet_origins[j]);
			mpz_add(zTempVar, zPrimeOrigin, zPrimorialMod);

			mpz_set(zFirstSieveElement, zTempVar);

			for(unsigned int i=nPrimorialEndPrime; i<nPrimeLimit; i++)
			{
				unsigned long  p = primes[i];
				unsigned int inv = inverses[i];
				unsigned int base_remainder = mpz_tdiv_ui(zTempVar, p);

				unsigned int remainder = base_remainder;
				unsigned long r = (p-remainder)*inv;
				unsigned int index = r % p;
				while(index < nBitArray_Size)
				{
					//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
						bit_array_sieve[(index)>>3] |= (1<<((index)&7));
					index += p;
				}
		
				remainder = base_remainder + 2;
				if (p<remainder)
					remainder -= p;
				r = (p-remainder)*inv;
				index = r % p;
				while(index < nBitArray_Size)
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
				while(index < nBitArray_Size)
				{
					//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
						bit_array_sieve[(index)>>3] |= (1<<((index)&7));
					index += p;
				}

				remainder = base_remainder + 8;
				if (p<remainder)
					remainder -= p;
				r = (p-remainder)*inv;
				index = r % p;
				while(index < nBitArray_Size)
				{
					//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
						bit_array_sieve[(index)>>3] |= (1<<((index)&7));
					index += p;
				}
		
				if(nHeight != nBestHeight)
				{
					goto request_new_block;
				}

				if ( nDifficulty > 50000000)
				{
					remainder = base_remainder + 12;
					if (p<remainder)
						remainder -= p;
					r = (p-remainder)*inv;
					index = r % p;
					while(index < nBitArray_Size)
					{
						//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
							bit_array_sieve[(index)>>3] |= (1<<((index)&7));
						index += p;
					}
				}

				if ( nDifficulty > 60000000)
				{
					remainder = base_remainder + 18;
					if (p<remainder)
						remainder -= p;
					r = (p-remainder)*inv;
					index = r % p;
					while(index < nBitArray_Size)
					{
						//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
							bit_array_sieve[(index)>>3] |= (1<<((index)&7));
						index += p;
					}
				}

				if ( nDifficulty > 70000000)
				{
					remainder = base_remainder + 20;
					if (p<remainder)
						remainder -= p;
					r = (p-remainder)*inv;
					index = r % p;
					while(index < nBitArray_Size)
					{
						//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
							bit_array_sieve[(index)>>3] |= (1<<((index)&7));
						index += p;
					}
				}

				if ( nDifficulty > 80000000)
				{
					remainder = base_remainder + 26;
					if (p<remainder)
						remainder -= p;
					r = (p-remainder)*inv;
					index = r % p;
					while(index < nBitArray_Size)
					{
						//if( !(bit_array_sieve[(i)>>3] & (1<<((i)&7))) )
							bit_array_sieve[(index)>>3] |= (1<<((index)&7));
						index += p;
					}
				}
			}

			//nElapsedTime = GetTimeMicros() - nStartTime;
			//printf("Sieving time: %uus\n", (unsigned int) nElapsedTime);

			if(nHeight != nBestHeight)
			{
				goto request_new_block;
			}

			//nStartTime = GetTimeMicros();
			for(i=0; i<nBitArray_Size; i++)
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

				if(nHeight != nBestHeight)
				{
					goto request_new_block;
				}

				nStart = 0;
				nStop = 2;
				nPrimeCount = 0;
				nLastOffset = 0;

				for(nStart; nStart<=nStop+12; nStart+=2)
				{
					mpz_sub_ui(zN, zTempVar, 1);
					mpz_powm(zResidue, zTwo, zN, zTempVar);
					if (mpz_cmp_ui(zResidue, 1) == 0)
					{
						nStop=nStart;
						nPrimeCount++;
					}

					mpz_add_ui(zTempVar, zTempVar, 2);
					nLastOffset+=2;
				}

				// calculate offset
				mpz_sub(zTempVar, zPrimeOriginOffset, zPrimeOrigin);

				nNonce = mpz_get_ui(zTempVar);

				nSieveDifficulty = 0;
				if (nPrimeCount > 3)
					nSieveDifficulty = SetBits( GetPrimeDifficulty2(BaseHash + nNonce + nLastOffset, nPrimeCount));

				if(nSieveDifficulty >= 40000000)
					nFourChainsFoundCounter++;

				if(nSieveDifficulty >= 50000000)
					nFiveChainsFoundCounter++;

				if(nSieveDifficulty >= 60000000)
					nSixChainsFoundCounter++;

				if(nSieveDifficulty >= 70000000)
					nSevenChainsFoundCounter++;

				if(nSieveDifficulty >= 80000000)
					nEightChainsFoundCounter++;

				if(nSieveDifficulty >= 50000000)
					printf("\n  Prime Clusters Found: Four=%u | Five=%u | Six=%u | Seven=%u | Eight=%u\n\n", nFourChainsFoundCounter,nFiveChainsFoundCounter,nSixChainsFoundCounter,nSevenChainsFoundCounter,nEightChainsFoundCounter);

				if(nSieveDifficulty >= nDifficulty)
				{
					mpz_clear(zPrimeOrigin);
					mpz_clear(zPrimeOriginOffset);
					mpz_clear(zFirstSieveElement);
					mpz_clear(zResidue);
					mpz_clear(zTwo);
					mpz_clear(zN);
					mpz_clear(zPrimorialMod);
					mpz_clear(zTempVar);

					free(bit_array_sieve);
					return nNonce;
				}
			}

			nElapsedTime = GetTimeMicros() - nStartTime;
			//printf("Exponentiation time: %uus\n", (unsigned int) nElapsedTime);
			//printf("Total time: %uus\n", (unsigned int) nElapsedTime);
		}

request_new_block:
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

