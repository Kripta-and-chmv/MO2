#include "Simplex.h"

void main()
{
	Simplex A;
	string path;
	A.Read(path);
	Point result = A.DoAlgorithm();

	cout << "Point: (" << result.x << "; " << result.y<<")\n";
	return;
}