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
	int n = 2; //размерность функции
	Point _center, _expansion, _reflection, _contraction;
	
	

	double Function(Point p)
	{
		int A1 = 1, A2 = 3;
		int a1 = 2, a2 = 1;
		int b1 = 3, b2 = 1;
		int c1 = 2, c2 = 1;
		int d1 = 3, d2 = 2;
		//return (A1 / (1 + pow(((p.x - a1) / b1), 2) + pow(((p.y - c1) / d1), 2))) + (A2 / (1 + pow(((p.x - a2) / b2), 2) + pow(((p.y - c2) / d2), 2)));
		return 100 * pow((p.y - pow(p.x, 2)), 2) + pow((1 - p.x), 2);
	}
	double Reflection()
	{
		_reflection = _center + _a * (_center - _polygon[_indMax]);
		return Function(_reflection);
	}
	double Expansion()
	{
		_expansion = _center + _y*(_reflection - _center);
		return Function(_expansion);
	}
	double Contraction()
	{
		_contraction = _center + _b*(_polygon[_indMax] - _center);
		return Function(_contraction);
	}
	void Reduction()
	{
		vector<Point> newPolygon;
		newPolygon.push_back(_polygon[_indMin]);

		for (int i = 0; i < _polygon.size(); i++)
		{
			if (i!=_indMin)
				newPolygon.push_back(newPolygon[0] + 0.5*(_polygon[1] - newPolygon[0]));
		}
		//newPolygon.push_back(newPolygon[0]+0.5*(_polygon[1]- newPolygon[0]));
		//newPolygon.push_back(newPolygon[0] + 0.5*(_polygon[2] - newPolygon[0]));
		_polygon = newPolygon;
	}

	void CreateSimplex()
	{
		//_polygon.push_back(_polygon[0] + Point{ 1, 0 }*_k);
		//_polygon.push_back(_polygon[0] + Point{ 0, 1 }*_k);
		double d1 = ((sqrt(n + 1) + n - 1) / (n*sqrt(2)))*_a;
		double d2 = ((sqrt(n + 1) - 1) / (n*sqrt(2)))*_a;
		_polygon.push_back(Point(_polygon[0].x + d1, _polygon[0].x + d2));
		_polygon.push_back(Point(_polygon[0].x + d2, _polygon[0].x + d1));
	}
	void SearchMinMax()
	{
		vector<double> result;

		for each (Point p in _polygon)
		{
			result.push_back(Function(p));
		}
		
		_fMax = _fMin = result[0];

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
				_center += _polygon[i] / (_polygon.size() - 1);
		}

	}
	bool Exit()
	{
		double summ = 0, result = 0;
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

	void FourthStep()
	{
		double fReflect = Function(_reflection);

		if (fReflect < Function(_polygon[_indMin]))
		{
			Expansion();
			if (Function(_expansion) < fReflect)
				_polygon[_indMax] = _expansion;
			else
				_polygon[_indMax] = _reflection;
			return;
		}

		if ((Function(_polygon[_indMin]) <= fReflect) && (fReflect < Function(_polygon[_indMax])))
		{
			Contraction();
			if (Function(_contraction) < fReflect)
			{
				_polygon[_indMax] = _contraction;
			}
			else
			{
				_polygon[_indMax] = _reflection;
			}
			return;
		}
		if (fReflect >= Function(_polygon[_indMax]))
		{
			Reduction();
			return;
		}
	}

public:
	Point DoAlgorithm(double &k)
	{
		CreateSimplex();
		while (true)
		{
			SearchMinMax();
			double fReflect = Reflection();
			if (fReflect >= Function(_polygon[_indMax]))
			{
				Reduction();

				if (Exit())
				{
					k = Function(_center);
					return _center;
				}
				continue;
			}
			if (fReflect < Function(_polygon[_indMin]))
			{
				Expansion();
				if (Function(_expansion) < Function(_polygon[_indMin]))
					_polygon[_indMin] = _expansion;
				else
					_polygon[_indMin] = _reflection;
				continue;
			}
			if ((Function(_polygon[_indMin]) <= fReflect) && (fReflect < Function(_polygon[_indMax])))
			{
				Contraction();
				_polygon[_indMax] = _contraction;
				continue;
			}

		}
	}
	void Read(string path)
	{
		ifstream read(path, ios_base::in);
		Point buf;
		read >> buf.x>> buf.y>>_a>>_y>>_b>>_k>> _EPS;

		_polygon.push_back(buf);
		read.close();
	}

};