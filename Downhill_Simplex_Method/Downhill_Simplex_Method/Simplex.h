#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

struct Point
{
	double x;
	double y;
	Point(double x, double y) : x(x), y(y) {};
	Point() {};

};

Point operator+ (const Point& a, const Point& b)
{
	return Point(a.x + b.x, a.y + b.y);
}
Point operator- (const Point& a, const Point& b)
{
	return Point(a.x - b.x, a.y - b.y);
}
Point operator* (const Point&a, double b)
{
	return Point(a.x*b, a.y*b);
}
Point operator* (double b, const Point&a)
{
	return Point(a.x*b, a.y*b);
}
Point operator/ (const Point& a, double b)
{
	return Point(a.x / b, a.y / b);
}
Point& operator+= (Point&a, const Point &b)
{
	a.x += b.x;
	a.y += b.y;
	return a;
}

class Simplex
{
private:
	double _a, _y, _b; //coefficients for Reflceion, Expansion and Contraction
	double _k;//scaled factor 
	double _EPS;
	//тот самый многоугольник. изначально, I элемент - начальная точка, 
	//II и III элементы - построенные от нач точки.
	vector<Point> _polygon;
	
	double _fMin, _fMax;
	int _indMax, _indMin;
	Point _center;
	
	void Read(string path)
	{
		ifstream read(path, ios_base::in);
		Point buf;
		read >> buf.x, buf.y, _a, _y, _b, _k, _EPS;

		_polygon.push_back(buf);
		read.close();
	}

	double Function(Point p)
	{
		int A1 = 1, A2 = 3;
		int a1 = 2, a2 = 1;
		int b1 = 3, b2 = 1;
		int c1 = 2, c2 = 1;
		int d1 = 3, d2 = 2;
		return (A1 / (1 + pow(((p.x - a1) / b1), 2) + pow(((p.y - c1) / d1), 2))) + (A2 / (1 + pow(((p.x - a2) / b2), 2) + pow(((p.y - c2) / d2), 2)));
	}
	void Reflection()
	{
		Point reflect = _center + _a * (_center - _polygon[_indMax]);
		_polygon.push_back(reflect);
	}
	void Expansion()
	{
		
	}
	void Contraction()
	{
		
	}
	void Reduction()
	{
		
	}

	void CreateSimplex()
	{
		_polygon[1] = _polygon[0] + Point{1, 0}*_k;
		_polygon[2] = _polygon[0] + Point{ 0, 1 }*_k;
	}
	void FirstStep()
	{
		vector<double> result;
		_fMax = _fMin = result[0];

		for each (Point p in _polygon)
		{
			result.push_back(Function(p));
		}
		
		for (int i = 0; i < result.size(); i++)
		{
			if (result[i] > _fMax)
			{
				_fMax = result[i];
				_indMax = i;
			}
			if (result[i] < _fMin)
			{
				_fMin = result[i];
				_indMin = i;
			}
		}
		
		for (int i = 0; i < result.size(); i++)
		{
			if (i != _indMax)
				_center += _polygon[i] / (_polygon.size() - 2);
		}

	}
	bool SecondStep()
	{
		double summ=0, result=0;
		for (int i = 0; i < _polygon.size(); i++)
			summ += pow(Function(_polygon[i]) - Function(_center), 2);

		result = sqrt(summ / _polygon.size());

		if (result <= _EPS)
			return true;
		else
			return false;
	}

	void ThirdStep()
	{
		Reflection();

	}

public:
	Point DoAlgorithm()
	{
		FirstStep();
		if (SecondStep()) return _center; else ThirdStep();

	}

};