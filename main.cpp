#include<iostream>
#include<cassert>
#include <initializer_list>
#include <algorithm>
#include <ctime>
#include <cstdlib>

#include <chrono>

//#include <type_traits>

template <typename T>
class Matrix;
template <typename T>
class Diagmatrix;
template <typename T=int>
class Unitmatrix;
template <typename T=int>
class Randmatrix;


template <typename T>
class Matrix
{
	template <typename U>
	friend class Matrix; // e.g. operator+, which is friend of one instance (e.g. Matrix<int>) has to access private members of other instance (e.g. Matrix<double>)

	protected:
	int m_rows;
	int m_columns;
	T* m_pntr;
	void resize(int rows, int columns){m_rows=rows; m_columns=columns; delete[] m_pntr; m_pntr = new T[m_rows*m_columns];}
	public:

	// should become typeindependent
	// initializer: list
	Matrix<T>(int n=0, int m=0, const std::initializer_list<T> &list=std::initializer_list<T>());

	// should become typeindependent
	// initializer: matrix
	// with reference, problem with rvalues...

	Matrix<T>(const Matrix<T> &toinitialize);
	Matrix<T>(Matrix<T> &&toinitialize);

	// destructor
	virtual ~Matrix<T>()
	{
		//std::cout << "matrix destroyed" << std::endl;
		delete[] m_pntr;
	}

	template <typename U>
	operator Matrix<U>()
	{

		Matrix<U> temp(m_rows,m_columns);
		for (int row=1; row<=m_rows; ++row)
			for (int column=1; column<=m_columns; ++column)
				temp(row,column)=static_cast<U>((*this)(row,column));
		return temp;
	}

	// rewrite entry
	T& operator()(int row, int column);

	// access entry
	const T& operator()(int row, int column) const;

	// output operator
	friend std::ostream& operator<<(std::ostream& out, const Matrix<T> &depicmatrix)
	{
		for (int row=1; row<=depicmatrix.m_rows; ++row)
			{
				for (int column=1; column<=depicmatrix.m_columns; ++column)
				{
					if (abs(depicmatrix(row,column))>1e-10) // only show values which are not super small
						out << depicmatrix(row,column) << '\t';
					else
						out << 0 << '\t';
				}
				out << std::endl;
			}
		return out;
	}

	// delegate work via print
	// input operator
	friend std::istream& operator>>(std::istream& in, Matrix<T> &inmatrix)
	{
		for (int row=1; row<=inmatrix.m_rows; ++row)
			for (int column=1; column<=inmatrix.m_columns; ++column)
				in >> inmatrix(row,column);
		return in;
	}

	// should become typeindependent
	// assignment operator matrix-matrix
	Matrix<T>& operator=(const Matrix<T> &toassign);
	Matrix<T>& operator=(Matrix<T> &&toassign);

	// should become typeindependent
	// assignment operator matrix-initilizer_list
	Matrix<T>& operator=(const std::initializer_list<T> &list);

	// transpose matrix
	virtual Matrix<T>& trans();

	template <typename V>
	friend Matrix<V> trans(const Matrix<V> &matrix);

	// calculate Determinant
	virtual T det();

	template <typename V>
	friend auto addmultipleofrow(const Matrix<V> &matrix, double multiple,int rowtoaddto,int row, double declparam) -> Matrix<decltype(matrix(1,1)+declparam)>;
	template <typename V>
	friend auto multiplyrow(const Matrix<V> &matrix,double multiple,int row, double declparam) -> Matrix<decltype(matrix(1,1)+declparam)>;
	template <typename V>
	friend auto switchrow(const Matrix<V> &matrix, int firstrow, int secondrow, double declparam) -> Matrix<decltype(matrix(1,1)+declparam)>;
	template <typename V>
	friend auto rowechelon(const Matrix<V> &matrix, double declparam) -> Matrix<decltype(matrix(1,1)+declparam)>;
	template <typename V>
	friend auto inv(const Matrix<V> &matrix, double declparam) -> Matrix<decltype(matrix(1,1)+declparam)>;
	template <typename V, typename W>
	friend auto solveSystem(const Matrix<V> &A, const Matrix<W> &b, double declparam) -> Matrix<decltype(A(1,1)+b(1,1)+declparam)>;

	// Note Diagmatrix+Diagmatrix = Matrix at the moment
	// Addition of two matrices
	template <typename V, typename W>
	friend auto operator+(const Matrix<V> &lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs(1,1)+rhs(1,1))>;

	// Subtraction of two matrices
	template <typename V, typename W>
	friend auto operator-(const Matrix<V> &lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs(1,1)+rhs(1,1))>;


	// Scalarmultipication from left for int and double
	template <typename W>
	friend auto operator*(int lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs*rhs(1,1))>;
	template <typename W>
	friend auto operator*(double lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs*rhs(1,1))>;

	// Scalarmultipication from right for int and double
	template <typename W>
	friend auto operator*(const Matrix<W> &lhs, int rhs) -> Matrix<decltype(lhs(1,1)*rhs)>;
	template <typename W>
	friend auto operator*(const Matrix<W> &lhs, double rhs) -> Matrix<decltype(lhs(1,1)*rhs)>;


	// Multiplikation of two matrices
	template <typename V, typename U>
	friend auto operator*(const Matrix<V> &lhs, const Matrix<U> &rhs) -> Matrix<decltype(lhs(1,1)*rhs(1,1))>;

	Matrix<T> cut(int rows, int columns, int firstrow=1, int firstcolumn=1);
	int getrows() {return m_rows;}
	int getcolumns() {return m_columns;}
};







template <typename T>
Matrix<T>::Matrix(int n, int m, const std::initializer_list<T> &list):m_rows(n),m_columns(m)
{
	//std::cout << "matrix constructed" << std::endl;
	assert((n>=0)&&(m>=0));
	m_pntr= new T[n*m];
	if (list.size()>0)
	{

		assert(list.size()==n*m);
		int count(0);
		for (auto &element:list)
		{
			m_pntr[count]=element;
			++count;
		}
	}
	else
	{
		for (int row=1; row<=m_rows; ++row)
		{
			for (int column=1; column<=m_columns; ++column)
			{
				(*this)(row,column)=0;
			}
		}
	}
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &toinitialize):m_rows(toinitialize.m_rows),m_columns(toinitialize.m_columns)
{
	m_pntr= new T[m_rows*m_columns];
	for (int count=0; count<m_rows*m_columns; ++count)
		m_pntr[count]=toinitialize.m_pntr[count];

}

template <typename T>
Matrix<T>::Matrix(Matrix<T> &&toinitialize):m_pntr(toinitialize.m_pntr),m_rows(toinitialize.m_rows),m_columns(toinitialize.m_columns)
{
	toinitialize.m_rows=0;
	toinitialize.m_columns=0;
	toinitialize.m_pntr = nullptr;
	//std::cout << "lol wut" << std::endl;
}

template <typename T>
auto addmultipleofrow(const Matrix<T> &matrix, double multiple, int rowtoaddto,int row, double declparam=0.0) -> Matrix<decltype(matrix(1,1)+declparam)>
{
	//std::cout << "addmultipleofrow" << std::endl;
	assert((0<rowtoaddto)&&(rowtoaddto<=matrix.m_rows)&&(row<=matrix.m_rows)&&(0<row));
	Matrix<decltype(matrix(1,1)+declparam)> temp(matrix.m_rows,matrix.m_columns);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			temp(row,column)=matrix(row,column);
	if (multiple==0)
		return temp;
	for (int column=1; column<=matrix.m_columns; ++column)
		if (matrix(rowtoaddto,column)!=0)
			temp(row,column)+=multiple*matrix(rowtoaddto,column);
	return temp;
}

template <typename T>
auto multiplyrow(const Matrix<T> &matrix, double multiple,int row, double declparam=0.0) -> Matrix<decltype(matrix(1,1)+declparam)>
{
	//std::cout << "multiplyrow" << std::endl;
	assert((0<row)&&(row<=matrix.m_rows));
	Matrix<decltype(matrix(1,1)+declparam)> temp(matrix.m_rows,matrix.m_columns);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			temp(row,column)=matrix(row,column);
	for (int column=1; column<=matrix.m_columns; ++column)
		if (temp(row,column)!=0)
			temp(row,column)=multiple*temp(row,column);
	return temp;
}

template <typename T>
auto switchrow(const Matrix<T> &matrix, int firstrow, int secondrow, double declparam=0.0) -> Matrix<decltype(matrix(1,1)+declparam)>
{
	//std::cout << "switchrow" << std::endl;
	assert((0<firstrow)&&(firstrow<=matrix.m_rows)&&(secondrow<=matrix.m_rows)&&(0<secondrow));
	Matrix<decltype(matrix(1,1)+declparam)> temp(matrix.m_rows,matrix.m_columns);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			temp(row,column)=matrix(row,column);
	T saveoffirstrow;
	for	(int column=1; column<=matrix.m_columns; ++column)
	{
		temp(firstrow,column)=matrix(secondrow,column);
		temp(secondrow,column)=matrix(firstrow,column);
	}
	return temp;
}

//check out

template <typename T>
auto rowechelon(const Matrix<T> &matrix, double declparam=0.0) -> Matrix<decltype(matrix(1,1)+declparam)>
{
	//std::cout << "rowechelon" << std::endl;
	assert(matrix.m_rows==matrix.m_columns);
	Matrix<decltype(matrix(1,1)+declparam)> transformations(Unitmatrix<decltype(matrix(1,1)+declparam)>(matrix.m_rows));
	Matrix<decltype(matrix(1,1)+declparam)> temp(matrix.m_rows,matrix.m_columns);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			temp(row,column)=matrix(row,column);
	for	(int column=1; column<=matrix.m_columns; ++column)
	{
		for (int row=column; row<=matrix.m_rows; ++row)
		{
			if (temp(row,column)!=0)
			{
				transformations=switchrow(transformations,row,column);
				temp=switchrow(temp,row,column);

				break;
			}
		}
		if (temp(column,column)!=0)
		{
			transformations=multiplyrow(transformations,1/temp(column,column),column);
			temp=multiplyrow(temp,1/temp(column,column),column);
		}
		for (int row=(column+1); row<=matrix.m_rows; ++row)
		{
			transformations=addmultipleofrow(transformations,-temp(row,column),column,row);
			temp=addmultipleofrow(temp,-temp(row,column),column,row);
		}

	}
	return transformations;
}

template <typename T>
auto inv(const Matrix<T> &matrix, double declparam=0.0) -> Matrix<decltype(matrix(1,1)+declparam)>
{
	Matrix<decltype(matrix(1,1)+declparam)> temp(matrix.m_rows,matrix.m_columns);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			temp(row,column)=matrix(row,column);
	assert(temp.det()); // has to copy many many matrizes
	Matrix<decltype(matrix(1,1)+declparam)> rowoperations(rowechelon(temp));
	temp=trans(rowoperations*temp);
	Matrix<decltype(matrix(1,1)+declparam)> columnoperations(rowechelon(temp));
	Matrix<decltype(matrix(1,1)+declparam)> inverse(trans(columnoperations)*rowoperations);
	return inverse;
}

template <typename V, typename W>
auto solveSystem(const Matrix<V> &A, const Matrix<W> &b, double declparam=0.0) -> Matrix<decltype(A(1,1)+b(1,1)+declparam)>
{

	assert((A.m_columns==A.m_rows)&&(b.m_rows==A.m_columns)&&(b.m_columns==1));
	Matrix<decltype(A(1,1)+b(1,1)+declparam)> x(b.m_rows,1);
	Matrix<decltype(A(1,1)+b(1,1)+declparam)> A_tilde(rowechelon(A)*A);
	//assert(A_tilde.det());
	Matrix<decltype(A(1,1)+b(1,1)+declparam)> b_tilde(rowechelon(A)*b);
	for (int index=b.m_rows; index>=1; --index)
	{
		x(index,1)=b_tilde(index,1);
		for (int setindex=b.m_rows; setindex>index; --setindex)
			x(index,1)-=A_tilde(index,setindex)*x(setindex,1);

	}
	return x;
}



template <typename V, typename W>
auto operator+(const Matrix<V> &lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs(1,1)+rhs(1,1))>
{
	assert((lhs.m_rows==rhs.m_rows)&&(lhs.m_columns==rhs.m_columns));
	Matrix<decltype(lhs(1,1)+rhs(1,1))> temp(lhs.m_rows,lhs.m_columns);

	for (int row=1; row<=lhs.m_rows; ++row)
		for (int column=1; column<=lhs.m_columns; ++column)
			temp(row,column)=lhs(row,column)+rhs(row,column);

	return temp;
}

template <typename V, typename W>
auto operator-(const Matrix<V> &lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs(1,1)+rhs(1,1))>
{
	assert((lhs.m_rows==rhs.m_rows)&&(lhs.m_columns==rhs.m_columns));
	Matrix<decltype(lhs(1,1)+rhs(1,1))> temp(lhs.m_rows,lhs.m_columns);

	for (int row=1; row<=lhs.m_rows; ++row)
		for (int column=1; column<=lhs.m_columns; ++column)
			temp(row,column)=lhs(row,column)-rhs(row,column);

	return temp;
}


template <typename W>
auto operator*(int lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs*rhs(1,1))>
{
	Matrix<decltype(lhs*rhs(1,1))> temp(rhs.m_rows,rhs.m_columns);
	for (int row=1; row<=rhs.m_rows; ++row)
		for (int column=1; column<=rhs.m_columns; ++column)
			temp(row,column)=lhs*rhs(row,column);
	return temp;
}

template <typename W>
auto operator*(double lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs*rhs(1,1))>
{
	Matrix<decltype(lhs*rhs(1,1))> temp(rhs.m_rows,rhs.m_columns);
	for (int row=1; row<=rhs.m_rows; ++row)
		for (int column=1; column<=rhs.m_columns; ++column)
			temp(row,column)=lhs*rhs(row,column);
	return temp;
}

template <typename W>
auto operator*(const Matrix<W> &lhs, int rhs) -> Matrix<decltype(lhs(1,1)*rhs)>
{return (rhs*lhs);}

template <typename W>
auto operator*(const Matrix<W> &lhs, double rhs) -> Matrix<decltype(lhs(1,1)*rhs)>
{return (rhs*lhs);}
/*
template <typename V, typename W>
auto operator*(const Matrix<V> &lhs, W rhs) -> Matrix<decltype(lhs(1,1)*rhs)>
{
	Matrix<decltype(lhs(1,1)*rhs)> temp(lhs.m_rows,lhs.m_columns);
	for (int row=1; row<=lhs.m_rows; ++row)
		for (int column=1; column<=lhs.m_columns; ++column)
			temp(row,column)=lhs(row,column)*rhs;
	return temp;
}
*/

template <typename V, typename W>
auto operator*(const Matrix<V> &lhs, const Matrix<W> &rhs) -> Matrix<decltype(lhs(1,1)*rhs(1,1))>
{
	assert(lhs.m_columns==rhs.m_rows);
	Matrix<decltype(lhs(1,1)*rhs(1,1))> temp(lhs.m_rows,rhs.m_columns);
	for (int row=1; row<=lhs.m_rows; ++row)
		for (int column=1; column<=rhs.m_columns; ++column)
			for (int index=1; index<=lhs.m_columns; ++index)
				temp(row,column)+=lhs(row,index)*rhs(index,column);
	return temp;
}

template <typename T>
T& Matrix<T>::operator()(int row, int column)
{
	assert((row>0)&&(row<=m_rows)&&(column>0)&&(column<=m_columns));
	return(m_pntr[(row-1)*m_columns+column-1]);
}

template <typename T>
const T& Matrix<T>::operator()(int row, int column) const
{
	assert((row>0)&&(row<=m_rows)&&(column>0)&&(column<=m_columns));
	return(m_pntr[(row-1)*m_columns+column-1]);
}

// should become typeindependent
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &toassign)
{
	if (this==&toassign)
		return *this;
	assert((m_rows==toassign.m_rows)&&(m_columns==toassign.m_columns));
	for (int row=1; row<=m_rows; ++row)
			for (int column=1; column<=m_columns; ++column)
				(*this)(row,column)=toassign(row,column);
	return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&toassign)
{
	//std::cout << "gotta be kiddin" << std::endl;
	if (&toassign==this)
		return *this;
	delete[] m_pntr;
	m_rows=toassign.m_rows;
	m_columns=toassign.m_columns;
	m_pntr=toassign.m_pntr;
	toassign.m_rows=0;
	toassign.m_columns=0;
	toassign.m_pntr=nullptr;

	return *this;
}

// should become typeindependent
template <typename T>
Matrix<T>& Matrix<T>::operator=(const std::initializer_list<T> &list)
{
	assert(list.size()==m_rows*m_columns);
	int count(0);
	for (auto &element:list)
	{
		m_pntr[count]=element;
		++count;
	}
	return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::trans()
{
Matrix<T> transposed(m_columns,m_rows);
for (int row=1; row<=m_rows; ++row)
		for (int column=1; column<=m_columns; ++column)
			transposed(column,row)=(*this)(row,column);
resize(m_columns,m_rows);
*this=transposed;


return *this;
}

template <typename T>
Matrix<T> trans(const Matrix<T> &matrix)
{
	Matrix<T> transposed(matrix.m_columns,matrix.m_rows);
	for (int row=1; row<=matrix.m_rows; ++row)
		for (int column=1; column<=matrix.m_columns; ++column)
			transposed(column,row)=matrix(row,column);
	return transposed;
}

template <typename T>
T Matrix<T>::det()
{
	assert((m_rows==m_columns)&&(m_rows>0)&&(m_columns>0));
	if (m_rows==1)
	return m_pntr[0];
	else
	{
		int zerocount(0), maxzeros(0), laplaceindex(0);
		for (int row=1; row<=m_rows; row++)
		{
			for (int column=1; column<=m_columns; column++)
			{
				if ((*this)(row,column)==0)
				{
				++zerocount;
				}

			}
			if (std::max(zerocount,maxzeros)==zerocount)
			{
			laplaceindex=row;
			maxzeros=zerocount;
			}
			zerocount=0;
		}


		Matrix<T> temp(m_rows-1,m_rows-1);
		T laplaceext=0;
		for (int sumand=1; sumand<=m_rows; ++sumand)
		{
			for (int row=1; row<=(m_rows-1); ++row)
			{
				for (int column=1; column<=(m_rows-1); column++)
				{
					if (row<laplaceindex)
					{
						if (column<sumand)
						temp(row,column)=(*this)(row,column);
						else
						temp(row,column)=(*this)(row,column+1);
					}
					else
					{
						if (column<sumand)
						temp(row,column)=(*this)(row+1,column);
						else
						temp(row,column)=(*this)(row+1,column+1);
					}
				}
			}
			if ((*this)(laplaceindex,sumand)==0)
			{
			}
			else
			{
			laplaceext+=pow(-1,laplaceindex+sumand)*(*this)(laplaceindex,sumand)*temp.det();
			}
		}
		// right position? or in braces below?
		return laplaceext;

	}

}

template <typename T>
Matrix<T> Matrix<T>::cut(int numberrows, int numbercolumns, int firstrow, int firstcolumn)
{
	assert((1<=numberrows)&&(numberrows<=m_rows)&&(1<=numbercolumns)&&(numbercolumns<=m_columns));
	assert((1<=firstrow)&&(firstrow<=m_rows)&&(1<=firstcolumn)&&(firstcolumn<=m_columns));
	assert(((numberrows+firstrow-1)<=m_rows)&&((numbercolumns+firstcolumn-1)<=m_columns));

	Matrix<T> temp(numberrows,numbercolumns);
	for (int row=1; row<=numberrows; ++row)
		for (int column=1; column<=numbercolumns; ++column)
			temp(row,column)=(*this)(firstrow+row-1,firstcolumn+column-1);
	return temp;

}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class Diagmatrix : public Matrix<T>
{
	public:
	Diagmatrix<T>(int dim=0,const std::initializer_list<T> &list=std::initializer_list<T>()):Matrix<T>(dim,dim)
	{
		assert(dim>=0);

		if (list.size()>0)
		{
			assert(list.size()==dim);
			int count(0);
			for (auto &element:list)
			{
				Matrix<T>::m_pntr[count*Matrix<T>::m_columns+count]=element;
				++count;
			}
		}
	}

	friend std::istream& operator>>(std::istream& in, Diagmatrix<T> &inmatrix)
	{
		for (int index=1; index<=inmatrix.m_rows; ++index)
			in >> inmatrix(index,index);
		return in;
	}

	virtual ~Diagmatrix<T>() override
	{
	}



	virtual Diagmatrix<T>& trans() override
	{
		return *this;
	}

	virtual T det() override;


	Diagmatrix<T>& inv()
	{
		assert(det());
		for (int index=0; index<Matrix<T>::m_rows; ++index)
			Matrix<T>::m_pntr[index*Matrix<T>::m_rows+index]=1/Matrix<T>::m_pntr[index*Matrix<T>::m_rows+index];
		return *this;
	}
};


template <typename T>
T Diagmatrix<T>::det()
	{
		T det(1);
		for (int index=0; index<Matrix<T>::m_rows; ++index)
			det*=Matrix<T>::m_pntr[index*Matrix<T>::m_rows+index];
		return det;
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class Unitmatrix: public Diagmatrix<T>
{
	public:
	Unitmatrix<T>(int dim=0):Diagmatrix<T>(dim)
	{
		assert(dim>=0);
		for (int index=1; index<=dim; ++index)
			(*this)(index,index)=1;
	}
	virtual ~Unitmatrix() override
	{
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// only produces interger entries
template <typename T>
class Randmatrix: public Matrix<T>
{
public:
	T m_min;
	T m_max;
	static int getRandomNumber(T min, T max)
	{
		static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
		return static_cast<int>(rand() * fraction * (max - min + 1) + min);
	}
	Randmatrix<T>(int rows, int columns, T min, T max):Matrix<T>(rows,columns),m_min(min),m_max(max)
	{
		assert(m_min<=m_max);
		for (int row=1; row<=rows; ++row)
			for (int column=1; column<=columns; ++column)
				(*this)(row,column)=getRandomNumber(m_min,m_max);
	}

	virtual ~Randmatrix() override
	{
	}


	T getMin(){return m_min;}
	T getMax(){return m_max;}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
	Timer t;
	srand(static_cast<unsigned int>(time(0)));
	rand();
	Randmatrix<double> A(10,10,1,10);
	Randmatrix<double> b(10,1,1,10);
    std::cout << A;
    std::cout << "\n";
    std::cout << b;
	std::cout << A*solveSystem(A,b)-b;
	std::cout << t.elapsed() << std::endl;
	return 0;
}

