#include <iostream>
#include <nmmintrin.h>
#include <windows.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

using namespace std;

const int N = 1000;

float m[N][N];

typedef struct {
int k;
int t_id;
int n;
}threadParam_t1;

void m_reset(int);
void m_gauss(int);
void m_gauss_simd(int);

void *threadFunc1(void *param)
{
    threadParam_t1 *p = (threadParam_t1*)param;
    int k = p -> k; //消去的轮次
    int t_id = p -> t_id; //线程编号
    int n = p -> n;
    int i = k + t_id + 1; //获取自己的计算任务
    for (int j = k + 1; j < n; ++j)
    {
        m[i][j] = m[i][j] - m[i][k] * m[k][j];

    }
    m[i][k] = 0;
    pthread_exit(NULL);
    return NULL;
}

void pthread_1(int n)//动态线程
{
    for (int k = 0; k < n; ++k)
    {
        for (int j = k + 1; j < n; j++)
            m[k][j] = m[k][j] / m[k][k];
        m[k][k] = 1.0;
        int worker_count = n - 1 - k;
        pthread_t *handles = new pthread_t[worker_count];
        threadParam_t1 *param = new threadParam_t1[worker_count];
        for(int t_id = 0; t_id < worker_count; t_id++)
        {
            param[t_id].k = k;
            param[t_id].t_id = t_id;
            param[t_id].n = n;
        }
        for(int t_id = 0; t_id < worker_count; t_id++)
            pthread_create(&handles[t_id],NULL,threadFunc1,(void*)(&param[t_id]));
        for(int t_id = 0; t_id < worker_count; t_id++)
            pthread_join(handles[t_id],NULL);
    }
    return;
}

typedef struct {
int t_id; //线程 id
int n;
}threadParam_t2;

const int NUM_THREADS = 16;
sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; // 每个线程有自己专属的信号量
sem_t sem_workerend[NUM_THREADS];

//线程函数定义
void *threadFunc2(void *param) {
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;
    for (int k = 0; k < n; ++k)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞，等待主线完成除法操作（操作自己专属的信号量）

        //循环划分任务
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        {
            for (int j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k]=0.0;
        }
        sem_post(&sem_main); // 唤醒主线程
        sem_wait(&sem_workerend[t_id]); //阻塞，等待主线程唤醒进入下一轮
    }
    pthread_exit(NULL);
    return NULL;
}


void pthread_2(int n)//静态线程 + 信号量同步 + 二重循环
{
    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc2,(void*)(&param[t_id]));
    }
    for(int k = 0; k < n; ++k)
    {
        //主线程做除法操作
        for (int j = k+1; j < n; j++)
            m[k][j] = m[k][j] / m[k][k];
        m[k][k] = 1.0;
        //开始唤醒工作线程
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerstart[t_id]);
        //主线程睡眠（等待所有的工作线程完成此轮消去任务）
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_wait(&sem_main);
        // 主线程再次唤醒工作线程进入下一轮次的消去任务
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerend[t_id]);
    }
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);
    sem_destroy(&sem_main);
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        sem_destroy(&sem_workerstart[t_id]);
        sem_destroy(&sem_workerend[t_id]);
    }
    return;
}

int flag[NUM_THREADS] = {0};
void *threadFunc3(void *param) {
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;
    for (int k = 0; k < n; ++k)
    {
        while(flag[t_id] == 0);

        //循环划分任务
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        {
            for (int j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k]=0.0;
        }
        flag[t_id] = 0;
    }
    pthread_exit(NULL);
    return NULL;
}

void pthread_3(int n)//静态线程 + 忙等待同步 +二重循环
{
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc3,(void*)(&param[t_id]));
    }
    for(int k = 0; k < n; ++k)
    {
        //主线程做除法操作
        for (int j = k+1; j < n; j++)
            m[k][j] = m[k][j] / m[k][k];
        m[k][k] = 1.0;

        for (int j = 0; j < NUM_THREADS; j++)
            flag[j] = 1;

        while(1)
        {
            int temp = 0;
            for (int i = 0; i < NUM_THREADS; ++i)
            {
                if(flag[i] == 1)
                    temp = 1;
            }
            if(temp == 0)
                break;
        }
    }
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);
    return;
}


sem_t sem_leader;
sem_t sem_Divsion[NUM_THREADS-1];
sem_t sem_Elimination[NUM_THREADS-1];

void *threadFunc4(void *param)
{
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;

    for (int k = 0; k < n; ++k)
    {
        if (t_id == 0)
        {
            for (int j = k+1; j < n; j++)
                m[k][j ] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }
        else
            sem_wait(&sem_Divsion[t_id-1]); // 阻塞，等待完成除法操作

        // t_id 为 0 的线程唤醒其它工作线程，进行消去操作
        if (t_id == 0)
        {
            for (int i = 0; i < NUM_THREADS-1; ++i)
                sem_post(&sem_Divsion[i]);
        }

        //循环划分任务（同学们可以尝试多种任务划分方式）
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        //消去
        {
            for (int j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k]=0.0;
        }

        if (t_id == 0)
        {
            for (int i = 0; i < NUM_THREADS-1; ++i)
                sem_wait(&sem_leader); // 等待其它 worker 完成消去

            for (int i = 0; i < NUM_THREADS-1; ++i)
                sem_post(&sem_Elimination[i]); // 通知其它 worker 进入下一轮

        }
        else
        {
            sem_post(&sem_leader);// 通知 leader, 已完成消去任务
            sem_wait(&sem_Elimination[t_id-1]);
        } // 等待通知，进入下一轮

    }
    pthread_exit(NULL);
    return NULL;
}


void pthread_4(int n)//静态线程 + 信号量同步 + 三重循环
{
    sem_init(&sem_leader, 0, 0);
    for (int i = 0; i < NUM_THREADS-1; ++i)
    {
        sem_init(&sem_Divsion[i], 0, 0);
        sem_init(&sem_Elimination[i], 0, 0);
    }

    //创建线程
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];

    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc4,(void*)(&param[t_id]));
    }


    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);

     // 销毁所有信号量
    sem_destroy(&sem_leader);
    for(int t_id = 0; t_id < NUM_THREADS - 1; t_id++)
    {
        sem_destroy(&sem_Divsion[t_id]);
        sem_destroy(&sem_Elimination[t_id]);
    }

     return;
}

pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;
void *threadFunc5(void *param)
{
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;

    for (int k = 0; k < n; ++k)
    {
        if (t_id == 0)
        {
            for (int j = k+1; j < n; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }

        //第一个同步点
        pthread_barrier_wait(&barrier_Divsion);

        //循环划分任务（同学们可以尝试多种任务划分方式）
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        //消去
        {
            for (int j = k + 1; j < n; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k]=0.0;
        }


        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);

    }
    pthread_exit(NULL);
    return NULL;
}

void pthread_5(int n)//静态线程 + barrier同步 + 三重循环
{
    pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);

    //创建线程
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc5,(void*)(&param[t_id]));
    }

    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);

    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
    return;
}

//线程函数定义
void *threadFunc6(void *param) {
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;
    for (int k = 0; k < n; ++k)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞，等待主线完成除法操作（操作自己专属的信号量）

        //循环划分任务
        for(int j=k+1+t_id; j < n; j += NUM_THREADS)
        {
            for (int i = k + 1; i < n; ++i)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[j][k]=0.0;
        }
        sem_post(&sem_main); // 唤醒主线程
        sem_wait(&sem_workerend[t_id]); //阻塞，等待主线程唤醒进入下一轮
    }
    pthread_exit(NULL);
    return NULL;
}


void pthread_6(int n)//静态线程 + 信号量同步 + 二重循环
{
    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc6,(void*)(&param[t_id]));
    }
    for(int k = 0; k < n; ++k)
    {
        //主线程做除法操作
        for (int j = k+1; j < n; j++)
            m[k][j] = m[k][j] / m[k][k];
        m[k][k] = 1.0;
        //开始唤醒工作线程
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerstart[t_id]);
        //主线程睡眠（等待所有的工作线程完成此轮消去任务）
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_wait(&sem_main);
        // 主线程再次唤醒工作线程进入下一轮次的消去任务
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerend[t_id]);
    }
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);
    sem_destroy(&sem_main);
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        sem_destroy(&sem_workerstart[t_id]);
        sem_destroy(&sem_workerend[t_id]);
    }
    return;
}

void *threadFunc7(void *param)
{
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;

    for (int k = 0; k < n; ++k)
    {
        if (t_id == 0)
        {
            for (int j = k+1; j < n; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }

        //第一个同步点
        pthread_barrier_wait(&barrier_Divsion);

        //循环划分任务（同学们可以尝试多种任务划分方式）
        for(int j=k+1+t_id; j < n; j += NUM_THREADS)
        //消去
        {
            for (int i = k + 1; i < n; ++i)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[j][k]=0.0;
        }


        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);

    }
    pthread_exit(NULL);
    return NULL;
}

void pthread_7(int n)//静态线程 + barrier同步 + 三重循环
{
    pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);

    //创建线程
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc7,(void*)(&param[t_id]));
    }

    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);

    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
    return;
}

void *threadFunc8(void *param) {
    __m128 vaik, vakj, vaij, vx;
    threadParam_t2 *p = (threadParam_t2*)param;
    int t_id = p -> t_id;
    int n = p -> n;
    for (int k = 0; k < n; ++k)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞，等待主线完成除法操作（操作自己专属的信号量）

        //循环划分任务
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        {
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
        sem_post(&sem_main); // 唤醒主线程
        sem_wait(&sem_workerend[t_id]); //阻塞，等待主线程唤醒进入下一轮
    }
    pthread_exit(NULL);
    return NULL;
}


void pthread_8(int n)//静态线程 + 信号量同步 + 二重循环
{
    __m128 vt, va;
    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc8,(void*)(&param[t_id]));
    }
    for(int k = 0; k < n; ++k)
    {
        //主线程做除法操作
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
        //开始唤醒工作线程
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerstart[t_id]);
        //主线程睡眠（等待所有的工作线程完成此轮消去任务）
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_wait(&sem_main);
        // 主线程再次唤醒工作线程进入下一轮次的消去任务
        for (int t_id = 0; t_id < NUM_THREADS; ++t_id)
            sem_post(&sem_workerend[t_id]);
    }
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);
    sem_destroy(&sem_main);
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        sem_destroy(&sem_workerstart[t_id]);
        sem_destroy(&sem_workerend[t_id]);
    }
    return;
}


void *threadFunc9(void *param)
{
    threadParam_t2 *p = (threadParam_t2*)param;
    __m128 vt, va, vaik, vakj, vaij, vx;
    int t_id = p -> t_id;
    int n = p -> n;

    for (int k = 0; k < n; ++k)
    {
        if (t_id == 0)
        {
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
        }

        //第一个同步点
        pthread_barrier_wait(&barrier_Divsion);

        //循环划分任务（同学们可以尝试多种任务划分方式）
        for(int i=k+1+t_id; i < n; i += NUM_THREADS)
        //消去
        {
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


        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);

    }
    pthread_exit(NULL);
    return NULL;
}

void pthread_9(int n)//静态线程 + barrier同步 + 三重循环
{
    pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);

    //创建线程
    pthread_t *handles = new pthread_t[NUM_THREADS];
    threadParam_t2 *param = new threadParam_t2[NUM_THREADS];
    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
    {
        param[t_id].t_id = t_id;
        param[t_id].n = n;
        pthread_create(&handles[t_id],NULL,threadFunc9,(void*)(&param[t_id]));
    }

    for(int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id],NULL);

    //销毁所有的 barrier
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
    return;
}

int main()
{
    long long head, tail , freq ; // timers
    int step = 10;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );

    for(int n = 10; n <= 1000; n += step)
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
        pthread_1(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_3(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_2(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_4(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_5(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;


        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_6(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_7(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_8(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;


        m_reset(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&head);
        pthread_9(n);
        QueryPerformanceCounter((LARGE_INTEGER *)&tail );
        cout << ( tail - head) * 1000.0 / freq << "ms" << endl;

        if(n == 100) step = 100;
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


