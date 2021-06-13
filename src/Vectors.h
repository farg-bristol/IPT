#ifndef VECTORS_H
#define VECTORS_H
/* Create a header for defining vectors and vector operations */
#include <cmath>

template<typename T, unsigned DIM> 
class vec
{
    static_assert(DIM > 0, "Dimension must be greater than 0");

    public:
        template<typename... Args>
        vec(Args&&... args) : data{T(args)...}
        {
            if(sizeof...(Args) != 1)
                static_assert(sizeof...(Args) == DIM, "Invalid number of constructor arguments.");
        }

    private:
        T data[DIM];

};

template<typename T, unsigned DIM>
vec<T,DIM> operator*(T a, vec<T,DIM> b)
{
    return b*a;
}


template <typename T>
class vec<T,2>
{
    public:
    
        vec<T,2>(){}

        vec<T,2>(T const a): data{a,a} {}

        vec<T,2>(T const a, T const b): data{a,b} {}

        vec<T,2>(vec<T,2> const& a) {*this = a;}
        
        void zero()
        {
            data[0] = 0.0; data[1] = 0.0;
        }

        T dot(vec<T,2> const& a) const
        {
            return data[0]*a[0] + data[1]*a[1];
        }

        T cross(vec<T,2> const& a) const
        {
            return (data[0]*a[1] - data[1]*a[0]);
        }

        vec<T,2> operator+(vec<T,2> const& a) const
        {
            return vec(data[0]+a[0], data[1]+a[1]);
        }

        vec<T,2> operator-(vec<T,2> const& a) const
        {
            return vec(data[0]-a[0], data[1]-a[1]);
        }

        T const& operator [](size_t const index)  const
        {
            return data[index];
        }

        size_t size(){return 2;}

    private:
        T data[2];
};


/*********************************************************/
/*****************  3D VECTOR CLASS **********************/
/*********************************************************/
template <typename T>
class vec<T,3>
{
    public:
        vec<T,3>(){}

        vec<T,3>(T const a): data{a,a,a} {}

        vec<T,3>(T const a, T const b, T const c): data{a,b,c} {}

        vec<T,3>(vec<T,3> const& a) {*this = a;}

        void zero()
        {
            data[0] = 0.0; data[1] = 0.0; data[2] = 0.0;
        }

        T dot(vec<T,3> const& a) const
        {
            return data[0]*a[0] + data[1]*a[1] + data[2]*a[2];
        }

        vec<T,3> cross(vec<T,3> const& a) const
        {
            return 
            vec<T,3>((data[1]*a[2] - data[2]*a[1]),
                (data[2]*a[0] - data[0]*a[2]),
                (data[0]*a[1] - data[1]*a[0]));
        }

        T const norm() const
        {
            return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
        }

        T sqnorm() const
        {
            return (data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
        }

        vec<T,3> operator=(T const a) const
        {
            return vec3(a);
        }

        vec<T,3> operator+(vec<T,3> const a) const
        {
            return vec(data[0] + a[0], data[1] + a[1], data[2] + a[2]);
        }

        vec<T,3> operator+(T const a) const
        {
            return vec(data[0] + a, data[1] + a, data[2] + a);
        }

        vec<T,3>& operator +=(vec<T,3> const& a)
        {
            data[0] = data[0] + a[0];
            data[1] = data[1] + a[1];
            data[2] = data[2] + a[2];
            return *this;
        }

        vec<T,3>& operator += (T const& a)
        {
            data[0] = data[0] + a;
            data[1] = data[1] + a;
            data[2] = data[2] + a;
            return *this;
        }

        vec<T,3> operator - (vec<T,3> const a) const
        {
            return vec<T,3>(data[0] - a[0], data[1] - a[1], data[2] - a[2]);
        }

        vec<T,3> operator - (T const a) const
        {
            return vec3(data[0] - a, data[1] - a, data[2] - a);
        }

        vec<T,3>& operator -= (vec<T,3> const& a)
        {
            data[0] = data[0] - a[0];
            data[1] = data[1] - a[1];
            data[2] = data[2] - a[2];
            return *this;
        }

        vec<T,3>& operator -=(T const& a)
        {
            data[0] = data[0] - a;
            data[1] = data[1] - a;
            data[2] = data[2] - a;
            return *this;
        }

        vec<T,3> operator / (T const a) const
        {
            return vec<T,3>(data[0] / a, data[1] / a, data[2] / a);
        }

        vec<T,3>& operator /= (T const& a)
        {
            data[0] = data[0] / a;
            data[1] = data[1] / a;
            data[2] = data[2] / a;
            return *this;
        }

        vec<T,3> operator * (T const a) const
        {
            return vec<T,3>(data[0] * a, data[1] * a, data[2] * a);
        }

        vec<T,3>& operator *= (T const a)
        {
            data[0] = data[0] * a;
            data[1] = data[1] * a;
            data[2] = data[2] * a;
            return *this;
        }

        vec<T,3>& operator = (vec<T,3> const& a)
        {
            data[0] = a[0];
            data[1] = a[1];
            data[2] = a[2];
            return *this;
        }
        
        vec<T,3>& operator = (T const& a)
        {
            data[0] = a;
            data[1] = a;
            data[2] = a;
            return *this;
        }

        T const& operator [](size_t const index)  const
        {
            return data[index];
        }

        size_t const size() const {return 3;} 

    private:
        T data[3];
};



/*********************************************************/
/*****************  4D VECTOR CLASS **********************/
/*********************************************************/
template <typename T>
class vec<T,4>
{
    public:
        vec<T,4>(){}

        vec<T,4>(T const a): data{a,a,a,a} {}

        vec<T,4>(T const a, T const b, T const c, T const d): data{a,b,c,d} {}

        vec<T,4>(vec<T,4> const& a) {*this = a;}

        vec<T,4>(vec<T,3> const& a, T const b)
        {
            data[0] = a[0]; data[1] = a[1]; data[2] = a[2]; data[3] = b;
        }

        void zero()
        {
            data[0] = 0.0; data[1] = 0.0; data[2] = 0.0; data[3] = 0.0;
        }

        T dot(vec<T,4> const& a)
        {
            return data[0]*a[0] + data[1]*a[1] + data[2]*a[2] + data[3]*a[3];
        }

        /* Doesn't Tly exist in 4D, and not useful */
        // vec4 pcross(vec4 const& a)
        // {
        //     return 
        //     vec4((data[1]*a[2] - data[2]*a[1]),
        //          (data[2]*a[0] - data[0]*a[2]),
        //          (data[0]*a[1] - data[1]*a[0]),
        //           data[0]*a[1] - data[1]*a[0]));
        // }

        vec<T,4> operator * (T const a) const
        {
            return vec<T,4>(data[0] * a, data[1] * a, data[2] * a, data[3] * a);
        }


        vec<T,4> operator =(T const a) const
        {
            return vec<T,4>(0.0);
        }

        vec<T,4> operator =(vec<T,4> const a) const
        {
            return a;
        }

        // vec4 operator+(vec4 const a) const
        // {
        //     return vec4(data[0] + a[0], data[1] + a[1], data[2] + a[2]);
        // }

        // vec4 operator+(T const a) const
        // {
        //     return vec3(data[0] + a, data[1] + a, data[2] + a);
        // }

        // vec3 operator*(T const a) const
        // {
        //     return vec3(data[0] * a, data[1] * a, data[2] * a);
        // }

        T const& operator [](size_t const index)  const
        {
            return data[index];
        }

        size_t size(){return 4;}

    private:
        T data[4];
};



/**************** MATRIX DEFINITIONS ***********************/


/* Only consider square matrices */
template<typename T, unsigned DIM> 
class matrix
{
    static_assert(DIM > 0, "Dimension must be greater than 0");

    public:

    matrix(void) {}
    
    template<typename... Args>
    matrix(Args&&... args) : data{T(args)...} {}

    T const& operator [](size_t const index)  const
    {
        return data[index];
    }

    private:
        T data[DIM*DIM];
};

template <typename T>
class matrix<T,2>
{
    public:
        matrix<T,2>(void) {}

        matrix<T,2>(T const a, T const b, T const c, T const d): data{a,b,c,d} {}

        matrix<T,2>(matrix<T,2> const& a): data{a} {}

        T det()
        {
            return (data[0]*data[1] - data[1]*data[0]);
        }

        T& operator()(int const row, int const col) { return data[row*2 + col]; }

        T const& operator () (int const row, int const col) const { return data[row*2 + col]; }

        size_t rows(){return 2;}
        size_t cols(){return 2;}
        size_t size(){return 4;}

    private:
        T operator [](size_t const index)  const
        {
            return data[index];
        }

        T data[4];
};

template <typename T>
class matrix<T,3>
{
    public:
        matrix(void) {}

        matrix<T,3>(T const& a): data{a,a,a,a,a,a,a,a,a}{}


        matrix<T,3>(T const a, T const b, T const c,
               T const d, T const e, T const f,
               T const g, T const h, T const i):
             data{a,b,c,d,e,f,g,h,i} {}

        matrix<T,3>(T const a[9]): data{a} {}

        matrix<T,3>(matrix<T,3> const& a) { *this = a;}

        T det() const
        {
            return 
            (data[0]*(data[4]*data[8] - data[7]*data[5]) - 
            data[1]*(data[3]*data[8] - data[6]*data[5]) +
            data[2]*(data[3]*data[7] - data[6]*data[4]));
        }

        matrix<T,3> transpose() const
        {
            return matrix<T,3>(
            data[0], data[3], data[6],
            data[1], data[4], data[7],
            data[2], data[5], data[8]);
        }

        matrix<T,3> operator*(T const& a) const
        {
            return matrix<T,3>(a*data[0],a*data[1],a*data[2],
                        a*data[3],a*data[5],a*data[5],
                        a*data[6],a*data[7],a*data[8]);
        }

        vec<T,3> operator*(vec<T,3> const& a) const
        {
            return  vec<T,3>(data[0]*a[0] + data[1]*a[1] + data[2]*a[2],
                        data[3]*a[0] + data[5]*a[1] + data[5]*a[2],
                        data[6]*a[0] + data[7]*a[1] + data[8]*a[2]);
        }


        matrix<T,3> operator*(matrix<T,3> const& a) const
        {
            matrix<T,3> b;

            b(0,0) = data[0]*a[0] + data[1]*a[3] + data[2]*a[6];
            b(0,1) = data[0]*a[1] + data[1]*a[4] + data[2]*a[7];
            b(0,2) = data[0]*a[2] + data[1]*a[5] + data[2]*a[8];

            b(1,0) = data[3]*a[0] + data[4]*a[3] + data[5]*a[6];
            b(1,1) = data[3]*a[1] + data[4]*a[4] + data[5]*a[7];
            b(1,2) = data[3]*a[2] + data[4]*a[5] + data[5]*a[8];

            b(2,0) = data[6]*a[0] + data[7]*a[3] + data[8]*a[6];
            b(2,1) = data[6]*a[1] + data[7]*a[4] + data[8]*a[7];
            b(2,2) = data[6]*a[2] + data[7]*a[5] + data[8]*a[8];

            return b;
        }

        size_t rows(){return 3;}
        size_t cols(){return 3;}
        size_t size(){return 9;}

        T& operator()(int const row, int const col) { return data[row*3 + col]; }

        T const& operator () (int const row, int const col) const { return data[row*3 + col]; }

    private:
        T operator [](size_t const index)  const
        {
            return data[index];
        }

        T data[9];
};


template <typename T>
class matrix<T,4>
{
    public:
        matrix<T,4>(){}

        matrix<T,4>(T const& a): data{a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a} {}

        matrix<T,4>(T const a, T const b, T const c, T const d,
               T const e, T const f, T const g, T const h,
               T const i, T const j, T const k, T const l,
               T const m, T const n, T const o, T const p):
             data{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p} {}

        matrix<T,4>(T const a[16]): data{a} {}


        matrix<T,4>(matrix<T,4> const& a) { *this = a;}
        // {
        //     data[0] = a[0];   data[1] = a[1];   data[2] = a[2];   data[3] = a[3]; 
        //     data[4] = a[4];   data[5] = a[5];   data[6] = a[6];   data[7] = a[7];
        //     data[8] = a[8];   data[9] = a[9];   data[10] = a[10]; data[11] = a[11];
        //     data[12] = a[12]; data[13] = a[13]; data[14] = a[14]; data[15] = a[15];
        // }

        void row(int const index, vec<T,4> const& a)
        {
            data[index*4+0] = a[0];
            data[index*4+1] = a[1];
            data[index*4+2] = a[2];
            data[index*4+3] = a[3];
        }

        // vec4& row(int const index)
        // {
        //     return vec4(*data[index*4],*data[index*4+1],*data[index*4+2],*data[index*4+3]);
        // }

        T det() const
        {
            matrix<T,3> a,b,c,d;
            a = matrix<T,3>(data[5], data[6], data[7],
                     data[9], data[10],data[11],
                     data[13],data[14],data[15]);

            b = matrix<T,3>(data[4], data[6], data[7],
                     data[8], data[10],data[11],
                     data[12],data[14],data[15]);

            c = matrix<T,3>(data[4], data[5], data[7],
                     data[8], data[9], data[11],
                     data[12],data[13],data[15]);

            d = matrix<T,3>(data[4], data[5], data[6],
                     data[8], data[9], data[10],
                     data[12],data[13],data[14]);

            return 
            (data[0]*a.det() - data[1]*b.det() + data[2]*c.det() - data[3]*d.det());
        }

        T determinant(void) const
        {
            /*
            We need to create a temporary
            */
            matrix<T, 4> temp(*this);
            /*We convert the temporary to upper triangular form*/
            T det = T(1);
            for (unsigned c = 0; c < 4; ++c)
            {
                det = det*temp(c,c);
                for (unsigned r = c + 1; r < 4; ++r)
                {
                    T ratio = temp(r, c) / temp(c, c);
                    for (unsigned k = c; k < 4; k++)
                    {
                        temp(r, k) = temp(r, k) - ratio * temp(c, k);
                    }
                }
            }

            return det;
        }

        matrix<T,4> transpose() const
        {
            return matrix(
            data[0], data[4], data[8], data[12],
            data[1], data[5], data[9], data[13],
            data[2], data[6], data[10], data[14],
            data[3], data[7], data[11], data[15]);
        }

        matrix<T,4> operator*(matrix<T,4> const& a)
        {
            matrix<T,4> b;

            b(0,0) = data[0]*a[0] + data[1]*a[4] + data[2]*a[8]  + data[3]*a[12];
            b(0,1) = data[0]*a[1] + data[1]*a[5] + data[2]*a[9]  + data[3]*a[13];
            b(0,2) = data[0]*a[2] + data[1]*a[6] + data[2]*a[10] + data[3]*a[14];
            b(0,3) = data[0]*a[3] + data[1]*a[7] + data[2]*a[11] + data[3]*a[15];

            b(1,0) = data[4]*a[0] + data[5]*a[4] + data[6]*a[8]  + data[7]*a[12];
            b(1,1) = data[4]*a[1] + data[5]*a[5] + data[6]*a[9]  + data[7]*a[13];
            b(1,2) = data[4]*a[2] + data[5]*a[6] + data[6]*a[10] + data[7]*a[14];
            b(1,3) = data[4]*a[3] + data[5]*a[7] + data[6]*a[11] + data[7]*a[15];

            b(2,0) = data[8]*a[0] + data[9]*a[4] + data[10]*a[8]  + data[11]*a[12];
            b(2,1) = data[8]*a[1] + data[9]*a[5] + data[10]*a[9]  + data[11]*a[13];
            b(2,2) = data[8]*a[2] + data[9]*a[6] + data[10]*a[10] + data[11]*a[14];
            b(2,3) = data[8]*a[3] + data[9]*a[7] + data[10]*a[11] + data[11]*a[15];

            b(3,0) = data[12]*a[0] + data[13]*a[4] + data[14]*a[8]  + data[15]*a[12];
            b(3,1) = data[12]*a[1] + data[13]*a[5] + data[14]*a[9]  + data[15]*a[13];
            b(3,2) = data[12]*a[2] + data[13]*a[6] + data[14]*a[10] + data[15]*a[14];
            b(3,3) = data[12]*a[3] + data[13]*a[7] + data[14]*a[11] + data[15]*a[15];

            return b;
        }

        unsigned const rows(){return 4;}
        unsigned const cols(){return 4;}
        unsigned const size(){return 16;}

        T& operator()(unsigned const row, unsigned const col) { return data[row*4 + col]; }

        T const& operator () (unsigned const row, unsigned const col) const { return data[row*4 + col]; }

    private:
        T operator [](unsigned const index)  const
        {
            return data[index];
        }

        T data[16];
};


#endif