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
	double _EPS;
	//тот самый многоугольник. изначально, I элемент - начальная точка, 
	//II и III элементы - построенные от нач точки.
	vector<Point> _polygon;
	
	
	int _indMax, _indMin, _indMid, _fMin=0, _fMax=0, _fMid=0;
	int n = 2; //размерность функции
	Point _center, _expansion, _reflection, _contraction;
	
	

	double Function(Point p)
	{
		functionCount++;
		int A1 = 1, A2 = 3;
		int a1 = 2, a2 = 1;
		int b1 = 3, b2 = 1;
		int c1 = 2, c2 = 1;
		int d1 = 3, d2 = 2;
		return -((A1 / (1 + pow(((p.x - a1) / b1), 2) + pow(((p.y - c1) / d1), 2))) + (A2 / (1 + pow(((p.x - a2) / b2), 2) + pow(((p.y - c2) / d2), 2))));
		return 100 * pow((p.y - pow(p.x, 2)), 2) + pow((1 - p.x), 2);
		//return (4 * pow((p.x - 5), 2)) + pow((p.y - 6), 2);
		//return 2 * pow(p.x, 2) + p.x*p.y + pow(p.y, 2);
	}
	double Reflection()
	{
		
		_reflection = _center+_a*(_center-_polygon[_indMax]);
		return Function(_reflection);
	}
	void Expansion()
	{
		_expansion = _center+ _y*(_reflection-_center);
	}
	void Contraction()
	{
		_contraction = _center + _b*(_polygon[_indMax]-_center);
	}
	void Reduction()
	{
		vector<Point> newPolygon;
		newPolygon.push_back(_polygon[_indMin]);

		for (int i = 0; i < _polygon.size(); i++)
		{
			if (i!=_indMin)
				newPolygon.push_back(newPolygon[0] + 0.5*(_polygon[i] - newPolygon[0]));
		}
		
		_polygon = newPolygon;
	}

	void Sort()
	{
		_center = Point(0, 0);
		vector<double> result;

		for each (Point p in _polygon)
		{
			result.push_back(Function(p));
		}
		
		double fMax, fMin;
		fMax = fMin = result[0];

		for (int i = 0; i < result.size(); i++)
		{
			if (result[i] >= fMax)
			{
				fMax = result[i];
				_indMax = i;
			}
			if (result[i] <= fMin)
			{
				fMin = result[i];
				_indMin = i;
			}
		}

		for (int i = 0; i < result.size(); i++)
		{
			if (i != _indMax)
				_center += _polygon[i] / (_polygon.size() - 1);
			if ((i != _indMax) && (i != _indMin))//находим среднюю точку
				_indMid = i;
		}
		_fMin = Function(_polygon[_indMin]);;
		_fMax = Function(_polygon[_indMax]);
		_fMid = Function(_polygon[_indMid]);

	}
	bool Exit()
	{
		double summ = 0, result = 0, bufCenter= Function(_center);
		for (int i = 0; i < _polygon.size(); i++)
			summ += pow(Function(_polygon[i]) - bufCenter, 2);

		result = sqrt(summ / _polygon.size());

		if (result <= _EPS)
			return true;
		else
			return false;
	}

public:
	int functionCount = 0;
	Point DoAlgorithm(double &k)
	{
		int kk = -1;
		#pragma region Cycle
				while (true)
				{
					
					kk++;
					Sort();
					if (Exit())
						break;

					double fReflect = Reflection();
					/////////////
					if (fReflect < _fMin)
					{
						Expansion();
						if (Function(_expansion) < _fMin)
							_polygon[_indMax] = _expansion;
						else
							_polygon[_indMax] = _reflection;
				
							continue;
					}
					////////////
					if((fReflect>_fMid)&&(fReflect<_fMax))
					{
						Contraction();
						_polygon[_indMax] = _contraction;
						continue;
					}
					//////////////
					if(fReflect>_fMin &&(fReflect<= _fMid))
					{
						_polygon[_indMax] = _reflection;
						continue;
					}
					/////////////////
					if(fReflect>_fMax)
					{
						Reduction();
						continue;
					}
				
				}
		#pragma endregion 
		Sort();
		k = _fMin;
		return _polygon[_indMin];
	}
	void Read(string path)
	{
		ifstream read(path, ios_base::in);
		Point p1, p2, p3;
		read >> p1.x >> p1.y >> p2.x >> p2.y >> p3.x >> p3.y >> _a >> _b >> _y >> _EPS;

		_polygon.push_back(p1);
		_polygon.push_back(p2);
		_polygon.push_back(p3);

		read.close();
	}

};