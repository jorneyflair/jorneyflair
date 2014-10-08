#include "core.h"

unsigned int nBlocksFoundCounter = 0;
unsigned int nBlocksAccepted = 0;
unsigned int nBlocksRejected = 0;
unsigned int nDifficulty = 0;
unsigned int nBestHeight = 0;
unsigned int nStartTimer = 0;
bool isBlockSubmission = false;

namespace Core
{

	/** Class to hold the basic data a Miner will use to build a Block.
		Used to allow one Connection for any amount of threads. **/
	class MinerThread
	{
	public:
		CBlock* BLOCK;
		bool fBlockFound, fNewBlock;
		LLP::Thread_t THREAD;
		boost::mutex MUTEX;
		
		unsigned int nSearches = 0, nPrimes = 0;
		
		MinerThread() : BLOCK(NULL), fBlockFound(false), fNewBlock(true), THREAD(boost::bind(&MinerThread::PrimeMiner, this)) { }
		
		/** Main Miner Thread. Bound to the class with boost. Might take some rearranging to get working with OpenCL. **/
		void PrimeMiner()
		{
			loop
			{
				try
				{
					/* Keep thread at idle CPU usage if waiting to submit or recieve block. **/
					Sleep(1);
					
					/** Increase nNonce if the block hash is too large. **/
					loop
					{
						//LOCK(MUTEX);
						
						if(fNewBlock || fBlockFound || !BLOCK)
							break;

						nDifficulty = BLOCK->nBits;
						BLOCK->nNonce = 0;

						unsigned long nNonce = PrimeSieve(BLOCK->GetPrime(), BLOCK->nBits, BLOCK->nHeight); 

						if((bool)nNonce)
						{						
							nPrimes++;
							BLOCK->nNonce = nNonce; 
		
							if(GetPrimeBits(BLOCK->GetPrime(), 1) >= BLOCK->nBits)
							{
								if(BLOCK->nHeight == nBestHeight)
								{
									if(!isBlockSubmission)
									{
										isBlockSubmission = true;
										printf("\n\n******* BLOCK FOUND *******\n\n");

										fBlockFound = true;
										fNewBlock = false;
										nBlocksFoundCounter++;
										break;
									}
								}
							}

						}
						
						fNewBlock = true;
					}
				}
				catch(std::exception& e){ printf("ERROR: %s\n", e.what()); }
			}
		}
	};
	
	
	/** Class to handle all the Connections via Mining LLP.
		Independent of Mining Threads for Higher Efficiency. **/
	class ServerConnection
	{
	public:
		LLP::Miner* CLIENT;
		int nThreads, nTimeout;
		std::vector<MinerThread*> THREADS;
		LLP::Thread_t THREAD;
		LLP::Timer    TIMER;
		std::string   IP, PORT;
		
		ServerConnection(std::string ip, std::string port, int nMaxThreads, int nMaxTimeout) : IP(ip), PORT(port), TIMER(), nThreads(nMaxThreads), nTimeout(nMaxTimeout), THREAD(boost::bind(&ServerConnection::ServerThread, this))
		{
			for(int nIndex = 0; nIndex < nThreads; nIndex++)
				THREADS.push_back(new MinerThread());
		}
		
		/** Reset the block on each of the Threads. **/
		void ResetThreads()
		{
		
			/** Reset each individual flag to tell threads to stop mining. **/
			for(int nIndex = 0; nIndex < THREADS.size(); nIndex++)
			{
				THREADS[nIndex]->fBlockFound = false;
				THREADS[nIndex]->fNewBlock   = true;
			}
				
		}
		
		/** Get the total Primes Found from Each Mining Thread.
			Then reset their counter. **/
		unsigned int Primes()
		{	
			unsigned int nPrimes = 0;
			for(int nIndex = 0; nIndex < THREADS.size(); nIndex++)
			{
				nPrimes += THREADS[nIndex]->nPrimes;
				THREADS[nIndex]->nPrimes = 0;
			}
			
			return nPrimes;
		}
		
		/** Get the total Numbers Searched from Each Mining Thread.
			Then reset their counter. **/
		unsigned int Searches()
		{	
			unsigned int nSearches = 0;
			for(int nIndex = 0; nIndex < THREADS.size(); nIndex++)
			{
				nSearches += THREADS[nIndex]->nSearches;
				THREADS[nIndex]->nSearches = 0;
			}
			
			return nSearches;
		}
		
		/** Main Connection Thread. Handles all the networking to allow
			Mining threads the most performance. **/
		void ServerThread()
		{
		
			/** Don't begin until all mining threads are Created. **/
			while(THREADS.size() != nThreads)
				Sleep(1);
				
				
			/** Initialize the Server Connection. **/
			CLIENT = new LLP::Miner(IP, PORT);
				
				
			/** Initialize a Timer for the Hash Meter. **/
			TIMER.Start();
			
			//unsigned int nBestHeight = 0;
			loop
			{
				try
				{
					/** Run this thread at 1 Cycle per Second. **/
					Sleep(1000);
					
					
					/** Attempt with best efforts to keep the Connection Alive. **/
					if(!CLIENT->Connected() || CLIENT->Errors())
					{
						ResetThreads();
						
						if(!CLIENT->Connect())
							continue;
						else
							CLIENT->SetChannel(1);
					}
					
					
					/** Check the Block Height. **/
					unsigned int nHeight = CLIENT->GetHeight(nTimeout);
					if(nHeight == 0)
					{
						printf("Failed to Update Height...\n");
						CLIENT->Disconnect();
						continue;
					}
					
					/** If there is a new block, Flag the Threads to Stop Mining. **/
					if(nHeight != nBestHeight)
					{
						isBlockSubmission = false;
						nBestHeight = nHeight;
						printf("[MASTER] Coinshield Network: New Block %u\n", nHeight);
							
						ResetThreads();
					}
					
					
					/** Rudimentary Meter **/
					if(TIMER.Elapsed() > 10)
					{
						unsigned int SecondsElapsed = (unsigned int)time(0) - nStartTimer;
						unsigned int nElapsed = TIMER.Elapsed();
							
						double SPS = (double) Searches() / nElapsed;
						double PPS = (double) Primes() / nElapsed;
						printf("[METERS] %u Block(s) ACC=%u REJ=%u| Height = %u | Diff = %f | %02d:%02d:%02d\n", nBlocksFoundCounter, nBlocksAccepted, nBlocksRejected, nBestHeight, (double)nDifficulty/10000000.0, (SecondsElapsed/3600)%60, (SecondsElapsed/60)%60, (SecondsElapsed)%60);
							
						TIMER.Reset();
					}

				
					/** Check if there is work to do for each Miner Thread. **/
					for(int nIndex = 0; nIndex < THREADS.size(); nIndex++)
					{
						if(THREADS[nIndex]->fBlockFound)
						{
							printf("\nPreparing for Block Submission...\n");
						}

						/** Attempt to get a new block from the Server if Thread needs One. **/
						if(THREADS[nIndex]->fNewBlock)
						{
							/** Retrieve new block from Server. **/
							CBlock* BLOCK = CLIENT->GetBlock(nTimeout);
							
							
							/** If the block is good, tell the Mining Thread its okay to Mine. **/
							if(BLOCK)
							{
								//LOCK(THREADS[nIndex]->MUTEX);
								
								THREADS[nIndex]->BLOCK = BLOCK;
								
								THREADS[nIndex]->fBlockFound = false;
								THREADS[nIndex]->fNewBlock   = false;
							}
							
							/** If the Block didn't come in properly, Reconnect to the Server. **/
							else
							{
								CLIENT->Disconnect();
								
								break;
							}
								
						}
						
						/** Submit a block from Mining Thread if Flagged. **/
						else if(THREADS[nIndex]->fBlockFound)
						{
							printf("\nSubmitting Block...\n");
							//LOCK(THREADS[nIndex]->MUTEX);
							
							if(!THREADS[nIndex]->BLOCK)
							{
								THREADS[nIndex]->fNewBlock   = true;
								THREADS[nIndex]->fBlockFound = false;
								
								continue;
							}
							
							
							printf("[MASTER] Prime Cluster of Difficulty %f Found on Thread %i\n", GetPrimeDifficulty(THREADS[nIndex]->BLOCK->GetPrime(), 1), nIndex);
							printf("\n%s\n\n", THREADS[nIndex]->BLOCK->GetHash().ToString().c_str());
								
							/** Attempt to Submit the Block to Network. **/
							unsigned char RESPONSE = CLIENT->SubmitBlock(THREADS[nIndex]->BLOCK->hashMerkleRoot, THREADS[nIndex]->BLOCK->nNonce, nTimeout);
							
							/** Check the Response from the Server.**/
							if(RESPONSE == 200)
							{
								printf("[MASTER] Block Accepted By Coinshield Network.\n");
								
								ResetThreads();
								nBlocksAccepted++;
							}
							else if(RESPONSE == 201)
							{
								printf("[MASTER] Block Rejected by Coinshield Network.\n");
								
								THREADS[nIndex]->fNewBlock   = true;
								THREADS[nIndex]->fBlockFound = false;

								isBlockSubmission = false;
								nBlocksRejected++;
							}
								
							/** If the Response was Bad, Reconnect to Server. **/
							else 
							{
								printf("[MASTER] Failure to Submit Block. Reconnecting...\n");
								CLIENT->Disconnect();
								
								break;
							}
						}		
					}
				}
				catch(std::exception& e)
				{
					printf("%s\n", e.what()); CLIENT = new LLP::Miner(IP, PORT); 
				}
			}
		}
	};
}

int main(int argc, char *argv[])
{

	if(argc < 3)
	{
		printf("Too Few Arguments. The Required Arguments are Ip and Port\n");
		printf("Default Arguments are Total Threads = CPU Cores and Connection Timeout = 10 Seconds\n");
		printf("Format for Arguments is 'IP PORT THREADS TIMEOUT'\n");
		
		Sleep(10000);
		
		return 0;
	}
		
	std::string IP = argv[1];
	std::string PORT = argv[2];
	int nThreads = GetTotalCores(), nTimeout = 10;
	
	if(argc > 3)
		nThreads = boost::lexical_cast<int>(argv[3]);
	
	if(argc > 4)
		nTimeout = boost::lexical_cast<int>(argv[4]);
	
	Core::InitializePrimes();
	nStartTimer = (unsigned int)time(0);
	printf("Initializing Miner %s:%s Threads = %i Timeout = %i\n", IP.c_str(), PORT.c_str(), nThreads, nTimeout);
	Core::ServerConnection MINERS(IP, PORT, nThreads, nTimeout);
	loop { Sleep(10); }
	
	return 0;
}
