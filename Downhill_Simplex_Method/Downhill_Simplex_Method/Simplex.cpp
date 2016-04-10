#include "Simplex.h"
void main()
{
	string path="input.txt";
	Point result;
	Simplex A;
	A.Read(path);
	double k;
	result = A.DoAlgorithm(k);

	cout << k<< "  "<< result.x << "  " << result.y << "\n";
	
	return;
}