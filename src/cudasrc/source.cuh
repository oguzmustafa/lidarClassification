#ifndef MCC_SOURCE_CUH
#define MCC_SOURCE_CUH

namespace mcc
{
	class IUnclassifiedPoints;
	
	class source
	{
	public:
		void clas(IUnclassifiedPoints & points);
		void kernel(double* A, double* B, double* C, int arraySize);
		void test(double ***A,int Asize);
	};
}

#endif