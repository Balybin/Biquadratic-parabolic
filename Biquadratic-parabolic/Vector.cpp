#include "Vector.h"


Vector::Vector(void)
{
}

//ведь возможно у нас получится вызвать конструктор с выделением памяти
Vector::Vector(int n)
{
        V.resize(n);
}

Vector::~Vector(void)
{
}

void Vector::make(int n)
{
        V.resize(n);
}

//перегрузка оператора присваивания
Vector& Vector::operator=(const Vector& newVect)
        {
                if (this != &newVect)
                        this->V = newVect.V;
                return *this;
        }

//перегрузка перемножения векторов ака скалярное произведение
double Vector::operator* (const Vector& newVect)
{
        double sum = 0;
        if (this->V.size() != newVect.V.size())
                return sum;
        else
        {
                for(int i = 0; i < this->V.size(); i++)
                        sum += this->V[i] * newVect.V[i];
        }
        return sum;
}

// норма вектора. и да можно взять корень из перемножения выше, но я ведь крутой программист
double Vector::norm()
{
        double norma = 0;
        for(int i = 0; i < this->V.size(); i++)
        {
			norma += pow(this->V[i],2.0);
        }
        norma = sqrt(norma);
        return norma;
}

//перегрузка разности векторов. результат окажется в правом векторе
Vector& Vector::operator-(Vector& newVect)
{
        if (this->V.size() != newVect.V.size())
        {
                return *this;
        }
        else
        {
                for (int i = 0; i < this->V.size(); i++)
                {
                        newVect.V[i] = this->V[i] - newVect.V[i]; 
                }
                return newVect;
        }
}

//перегрузка оператора сложения, результат окажется в левом векторе, а ведь я его даже вроде не использую
Vector& Vector::operator+(Vector& newVect)
{
        if (this->V.size() != newVect.V.size())
        {
                return *this;
        }
        else
        {
                for (int i = 0; i < this->V.size(); i++)
                {
                        this->V[i] = this->V[i] + newVect.V[i]; 
                }
                return *this;
        }

}