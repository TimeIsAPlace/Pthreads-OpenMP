#include <iostream>
#include <nmmintrin.h>
#include <windows.h>
#include <stdlib.h>
#include <omp.h>

using namespace std;

const int N = 3000;

float m[N][N];

void m_reset(int);
void m_gauss(int);
void m_gauss_simd(int);
void m_OpenMP_1(int);
void m_OpenMP_2(int);
void m_OpenMP_3(int);
void m_OpenMP_4(int);
void m_OpenMP_5(int);
void m_OpenMP_6(int);
void m_OpenMP_7(int);

int NUM_THREADS=16;

int main()
{
    long long head, tail , freq ; // timers
    int step = 10;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );

    for(int n = 10; n <= 3000; n += step)
    {
        cout << "n: " << n << endl;
        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_gauss(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout <<  ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_gauss_simd(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_1(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_2(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_3(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_4(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_5(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_6(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;


        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        m_OpenMP_7(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        if(n == 100) step = 100;
        if(n == 1000) step = 1000;
    }

    return 0;
}

//初始化矩阵元素
void m_reset(int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<i;j++)
            m[i][j]=0;
        m[i][ i]=1.0;
        for(int j=i+1;j<n;j++)
            m[i][j]=rand();
    }
    for(int k=0;k<n;k++)
        for(int i=k+1;i<n;i++)
            for(int j=0;j<n;j++)
                m[i][j]+=m[k][j];
}

//串行普通高斯消去算法
void m_gauss(int n)
{
    for(int k = 0 ; k < n ; k++)
    {
        for(int j = k+1 ; j < n ; j++)
        {
            m[k][j] = m[k][j]/m[k][k];
        }
        m[k][k] = 1.0;
        for(int i = k+1 ; i < n ; i++)
        {
            for(int j = k+1 ; j < n ; j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            m[i][k] = 0;

        }
    }
}

//乘法和除法部分全部向量化、不对齐
void m_gauss_simd(int n)
{
    __m128 vt, va, vaik, vakj, vaij, vx;
    for(int k = 0; k < n; k++){
        vt = _mm_set_ps1(m[k][k]);
        int j;
        for(j = k+1; j+4 <= n; j+=4){
            va = _mm_loadu_ps(&m[k][j]);
            va = _mm_div_ps(va, vt);
            _mm_storeu_ps(&m[k][j], va);
        }
        if(j < n){
            for(;j < n; j++){
                m[k][j] = m[k][j]/m[k][k];
            }
        }
        m[k][k] = 1.0;
        for(int i = k+1; i < n; i++){
            vaik = _mm_set_ps1(m[i][k]);
            int j;
            for(j = k+1; j+4 <= n; j+=4){
                vakj = _mm_loadu_ps(&m[k][j]);
                vaij = _mm_loadu_ps(&m[i][j]);
                vx = _mm_mul_ps(vakj, vaik);
                vaij = _mm_sub_ps(vaij, vx);
                _mm_storeu_ps(&m[i][j], vaij);
            }
            if(j < n){
                for(;j < n; j++){
                    m[i][j] = m[i][j] - m[k][j]*m[i][k];
                }
            }
            m[i][k] = 0;
        }
    }
}

int parallel = 1;
int i,j,k,tmp;
void m_OpenMP_1(int n)//除法串行+消去并行+行划分
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        #pragma omp single
        {
            tmp = m[k][k];
            for(j = k + 1; j < n; ++j)
                m[k][j] = m[k][j] / tmp;
            m[k][k] = 1.0;
        }
        // 并行部分，使用行划分
        #pragma omp for
        for(i = k + 1; i < n; ++i)
        {
            tmp = m[i][k];
            for(j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0.0;
        }
    }
}

void m_OpenMP_2(int n)//除法串行+消去并行+列划分
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        #pragma omp single
        {
            tmp = m[k][k];
            for(j = k + 1; j < n; ++j)
                m[k][j] = m[k][j] / tmp;
            m[k][k] = 1.0;
        }
        // 并行部分，使用行划分
        #pragma omp for
        for(j = k + 1; j < n; ++j)
        {
            tmp = m[j][k];
            for(i = k + 1; i < n; ++i)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[j][k] = 0.0;
        }
    }
}

void m_OpenMP_3(int n)//除法并行+消去串行+行划分
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        tmp = m[k][k];
        #pragma omp for
        for(j = k + 1; j < n; ++j)
            m[k][j] = m[k][j] / tmp;
        m[k][k] = 1.0;

        // 并行部分，使用行划分
        #pragma omp single
        {
            for(i = k + 1; i < n; ++i)
            {
                tmp = m[i][k];
                for(j = k + 1; j < n; ++j)
                    m[i][j] = m[i][j] - tmp * m[k][j];
                m[i][k] = 0.0;
            }
        }
    }
}

void m_OpenMP_4(int n)//除法并行+消去并行+行划分
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        tmp = m[k][k];
        #pragma omp for
        for(j = k + 1; j < n; ++j)
            m[k][j] = m[k][j] / tmp;
        m[k][k] = 1.0;
        #pragma omp for
        for(i = k + 1; i < n; ++i)
        {
            tmp = m[i][k];
            for(j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0.0;
        }
    }
}



void m_OpenMP_5(int n)//除法并行+消去并行+行划分+负载
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        tmp = m[k][k];
        #pragma omp for schedule(static,2)
        for(j = k + 1; j < n; ++j)
            m[k][j] = m[k][j] / tmp;
        m[k][k] = 1.0;
        // 并行部分，使用行划分
        #pragma omp for schedule(static,2)
        for(i = k + 1; i < n; ++i)
        {
            tmp = m[i][k];
            for(j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0.0;
        }
    }
}


void m_OpenMP_6(int n)//除法并行+消去并行+行划分+负载
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        tmp = m[k][k];
        #pragma omp for schedule(static,4)
        for(j = k + 1; j < n; ++j)
            m[k][j] = m[k][j] / tmp;
        m[k][k] = 1.0;
        // 并行部分，使用行划分
        #pragma omp for schedule(static,4)
        for(i = k + 1; i < n; ++i)
        {
            tmp = m[i][k];
            for(j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0.0;
        }
    }
}



void m_OpenMP_7(int n)//除法并行+消去并行+行划分+负载
{
    #pragma omp parallel if(parallel), num_threads(NUM_THREADS), private(i, j, k, tmp)
    for(k = 0; k < n; ++k)
    {
        // 串行部分，也可以尝试并行化
        tmp = m[k][k];
        #pragma omp for schedule(static,8)
        for(j = k + 1; j < n; ++j)
            m[k][j] = m[k][j] / tmp;
        m[k][k] = 1.0;
        // 并行部分，使用行划分
        #pragma omp for schedule(static,8)
        for(i = k + 1; i < n; ++i)
        {
            tmp = m[i][k];
            for(j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0.0;
        }
    }
}
