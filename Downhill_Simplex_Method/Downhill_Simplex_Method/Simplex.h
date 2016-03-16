#pragma once

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

struct Point
{
	double x;
	double y;
	Point(double x,double y) : x(x),y(y) {}

};

Point operator+ (const Point& a, const Point& b)
{
	return Point(a.x + b.x, a.y + b.y);
}
Point operator* (const Point&a, double b)
{
	return Point(a.x*b, a.y*b);
}

class Simplex
{
private:
	double a, y, b; //coefficients for Reflceion, Expansion and Contraction
	double k;//scaled factor
	vector<Point> polygon;
	Point fMin, fMax;
	
	double Function(Point p)
	{
		int A1 = 1, A2 = 3;
		int a1 = 2, a2 = 1;
		int b1 = 3, b2 = 1;
		int c1 = 2, c2 = 1;
		int d1 = 3, d2 = 2;
		return (A1 / (1 + pow(((p.x - a1) / b1), 2) + pow(((p.y - c1) / d1), 2))) + (A2 / (1 + pow(((p.x - a2) / b2), 2) + pow(((p.y - c2) / d2), 2)));
	}
	void Reflecion()
	{
	}
	void Expansion()
	{
		
	}
	void Contracion()
	{
		
	}
	void Reduction()
	{
		
	}

	void CreateSimplex()
	{
		polygon[1] = polygon[0] + Point{1, 0}*k;
		polygon[2] = polygon[0] + Point{ 0, 1 }*k;
	}
 

public:
	void smthng()
	{

	}

};