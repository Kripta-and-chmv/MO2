#include "Simplex.h"
void main()
{
	string path="input.txt";
	Point result;
	Simplex A;
	A.Read(path);
	result = A.DoAlgorithm();

	cout << result.x << "  " << result.y << "\n";
	
	return;
}