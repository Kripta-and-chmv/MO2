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
	
	
	int _indMax, _indMin, _indMid;
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
		
		_reflection = (1 + _a)*_center - _a*_polygon[_indMax];
		return Function(_reflection);
	}
	double Expansion()
	{
		_expansion = (1 - _y)*_center+_y*_reflection;
		return Function(_expansion);
	}
	double Contraction()
	{
		_contraction = _b*_polygon[_indMax] + (1 - _b)*_center;
		return Function(_contraction);
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

	void CreateSimplex()
	{
		_polygon.push_back(_polygon[0] + Point{ 1, 0 }*_k);
		_polygon.push_back(_polygon[0] + Point{ 0, 1 }*_k);
		/*double d1 = ((sqrt(n + 1) + n - 1) / (n*sqrt(2)))*_k;
		double d2 = ((sqrt(n + 1) - 1) / (n*sqrt(2)))*_k;
		_polygon.push_back(Point(_polygon[0].x + d1, _polygon[0].x + d2));
		_polygon.push_back(Point(_polygon[0].x + d2, _polygon[0].x + d1));*/
	}
	void Sort()
	{
		vector<double> result;

		for each (Point p in _polygon)
		{
			result.push_back(Function(p));
		}
		
		double fMax, fMin;
		fMax = fMin = result[0];

		for (int i = 0; i < result.size(); i++)
		{
			if (result[i] > fMax)
			{
				fMax = result[i];
				_indMax = i;
			}
			if (result[i] < fMin)
			{
				fMin = result[i];
				_indMin = i;
			}
		}

		for (int i = 0; i < result.size(); i++)
		{
			if (i != _indMax)
				_center += _polygon[i] / (_polygon.size() - 1);
			if ((i != _indMax) && (i != _indMid))//находим среднюю точку
				_indMid = i;
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

public:
	Point DoAlgorithm(double &k)
	{
		CreateSimplex();
		#pragma region Cycle
				while (true)
				{
					Sort();
					double fReflect = Reflection();
					/////////////
					if (fReflect < Function(_polygon[_indMin]))
					{
						Expansion();
						if (Function(_expansion) < fReflect)
							_polygon[_indMax] = _expansion;
						else
							_polygon[_indMax] = _reflection;
						if (Exit())
							break;
						else
							continue;
					}
					////////////
					if((fReflect>Function(_polygon[_indMin]))&&(fReflect<Function(_polygon[_indMid])))
					{
						_polygon[_indMax] = _reflection;
						if (Exit())
							break;
						else
							continue;
					}
					////////////
					if((fReflect>Function(_polygon[_indMid]))&&(fReflect<Function(_polygon[_indMax])))
					{
						swap(_polygon[_indMax], _reflection);
						Contraction();
						if(Function(_contraction)<Function(_polygon[_indMax]))
						{
							_polygon[_indMax] = _contraction;
							if (Exit())
								break;
							else
								continue;
						}
						else
						{
							Reduction();
							if (Exit())
								break;
							else
								continue;
						}
					}
					////////////
					if(fReflect>Function(_polygon[_indMax]))
					{
						Contraction();
						if (Function(_contraction)<Function(_polygon[_indMax]))
						{
							_polygon[_indMax] = _contraction;
							if (Exit())
								break;
							else
								continue;
						}
						else
						{
							Reduction();
							if (Exit())
								break;
							else
								continue;
						}
					}
				}
		#pragma endregion 
		Sort();
		k = Function(_polygon[_indMin]);
		return _polygon[_indMin];
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