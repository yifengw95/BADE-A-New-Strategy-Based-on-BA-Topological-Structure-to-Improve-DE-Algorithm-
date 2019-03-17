//
//  main.cpp
//  DE1.4
//
//  Created by 王译烽 on 16/4/8.
//  Copyright © 2016年 王译烽. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

int G[1000][1000];                             //Initialize BA Matrix

double evaluate (float *p, int f_number ){              // 测试函数集
    double deta = 0.0;
    
    if (f_number == 1) {
        deta = (*p) * (*(p)) + (*(p+1)) * (*(p+1)) + (*(p+2)) * (*(p+2)) - 0.0;
    }
    
    if (f_number == 2) {
        deta = 100 * ((*p) * (*p) - *(p+1)) *((*p) * (*p) - *(p+1)) + (1 -(*p) )*(1 -(*p) ) - 0.0;
    }
    
    if (f_number == 4) {
        float sum = 0;
        for (int i = 0; i < 30; i++) {
            srand((unsigned)time(NULL));
            rand();
            float n = float (rand() )/ RAND_MAX ;
            sum = sum + ((i+1) * pow ((*(p+i)), 4) + n);
        }
        deta = fabs (sum - 0);
    }
    
    if (f_number == 6) {
        float sum = 0;
        float x[4];
        float z[4];
        for (int i = 0; i < 4; i++) {
            x[i] = *(p+i);
        }
        for (int i = 0; i < 4; i++) {
            if (x[i] < 0) {
                z[i] = fabs (fabs(x[i] / 0.2) + 0.49999) * (-1) * 0.2;
            }
            else if (x[i] == 0){
                z[i] = fabs (fabs(x[i] / 0.2) + 0.49999) * (0) * 0.2;
            }
            else{
                z[i] = fabs (fabs(x[i] / 0.2) + 0.49999) * (1) * 0.2;
            }
        }
        float d[4] = {1, 1000, 10 , 100};
        for (int  i = 0; i < 4; i++) {
            if (fabs(x[i] - z[i]) < 0.05) {
                if (z[i] < 0) {
                    sum = sum + 0.15 * (z[i] - 0.05 * (-1)) * (z[i] - 0.05 * (-1)) * d[i];
                }
                else if (z[i] == 0){
                    sum = sum + 0.15 * (z[i] - 0.05 * (0)) * (z[i] - 0.05 * (0)) * d[i];
                }
                else{
                    sum = sum + 0.15 * (z[i] - 0.05 * (1)) * (z[i] - 0.05 * (1)) * d[i];
                }
            }
            else{
                sum = sum + d[i] * x[i] * x[i];
                
            }
        }
        deta = fabs(sum - 0);
        
    }
    
    
    
    if (f_number == 7) {
        float sum1 = 0 ;
        for (int i  = 0; i < 10; i++) {
            sum1 = sum1 + (*(p + i)) *(*(p + i)) / 4000;
        }
        float sum2 = 1;
        for (int i  = 0; i < 10; i++) {
            sum2 = sum2 * cos((*(p + i)) / sqrt(i+1));
        }
        deta = fabs (sum1 - sum2 + 1 - 0.0);
    }
    
    if (f_number == 11) {
        for (int i  = 0; i < 30; i++) {
            deta = deta + i * i * (*(p + i)) * (* (p+i)) - 0.0;
        }
    }
    
    if (f_number == 16) {
        deta = fabs (pow(*p, 6) - 15 * pow(*p, 4) + 27 * pow(*p, 2) + 250 - 7) ;
    }
    
    return deta ;
    
}

void model (int NP ,//种群数量
            float F ,
            float CR ,
            int D ,//个体基因数量
            long int gen_max ,
            float up,//测试函数自变量取值上界
            float down,//测试函数自变量取值下界
            float deta,//测试函数容许误差范围
            int m0,//初始节点数
            int m ,//每次引入的节点连接到m个结点上
            int t ,//经过t步形成BA网络
            int f_number){
    
    srand((unsigned)time(NULL));
    rand();
    double evaluate(float *, int);
    
    float x1[NP][D];
    float x2[NP][D];
    float Mut[NP];
    float trial[D];
    bool condition1 = 0;
    bool condition2 = 0;
    int countDE1 = 1;
    int countDE2 = 1;
    
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < D; j++) {
            x1[i][j] = float(rand() ) / RAND_MAX * (up - down) + down;
        }                                                               //DE_Initialize
    }
    
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < D; j++) {
            x2[i][j] = float(rand() ) / RAND_MAX * (up - down) + down;
        }                                                               //DE(BA)_Initialize
    }
    
    int k_degree[t+m0];//度
    int k_sum;
    
    int i,j;
    for (i=0;i<=m0-1;i++)
    {
        for (j=0;j<=m0-1;j++)                                      //BA_Initialize(primitive)
        {
            if (j!=i)
            {
                G[i][j]=1;
            }
            else
                continue;
        }
    }
    
    int ii,jj;                                                  // BA_Initialize(start)
    int count;//连的第几条边，count<=m
    double p[t+m0];//与各个结点相连的概率，可重复使用
    double p1;//产生的随机数
    double p2;//累加的概率
    
    for (i=m0;i<=t+m0-1;i++)
    {
        //首先计算i-1个结点的度以及连到它们的概率
        for (ii=0;ii<=t+m0-1;ii++)
        {k_degree[ii]=0;
        }
        for (ii=0;ii<=i-1;ii++)
        {
            for (j=0;j<=i-1;j++)
            {
                k_degree[ii]=k_degree[ii]+G[ii][j];
            }
        }
        
        k_sum=0;
        for (ii=0;ii<=i-1;ii++)
        {
            k_sum=k_degree[ii]+k_sum;
        }
        for (ii=0;ii<=i-1;ii++)
        {
            p[ii]=(double)k_degree[ii]/k_sum;
        }
        
        
        p1=rand()/(double)RAND_MAX;
        p2=0;
        for (j=0;j<=i-1;j++)
        {
            p2=p2+p[j];
            if (p2>=p1)
            {break;
            }
        }
        G[i][j]=1;
        G[j][i]=1;
        for (count = 1; count <= m - 1; count++)
        {
            double temp = p[j];
            for (jj = 0; jj <= i - 1; jj++)
            {
                p[jj] = p[jj] / (1 - temp);
            }
            p[j] = 0;
            p1 = rand() / (double)RAND_MAX;
            p2 = 0;
            for (j = 0; j <= i - 1; j++)
            {
                p2 = p2 + p[j];
                if (p2 >= p1)
                {
                    break;
                }
            }
            G[i][j] = 1;
            G[j][i] = 1;
        }
    }
    
    // BA_Initialize(finished)
    /*
     for (i=0;i<=t+m0-1;i++)
     {
     for (j=0;j<=t+m0-1;j++)
     {
     printf("%d ",G[i][j]);
     }
     printf("\n");
     }
     */
    while (countDE1 < gen_max && condition1 == 0) {                       //DE Algorithm(start)
        
        for (int  i = 0; i < NP; i++) {
            int a; do a = int (rand() % NP); while (a == i);
            int b; do b = int (rand() % NP); while (b == i || b == a);
            //    int c; do c = int (rand() % NP); while (c == i || c == a || c == b);//DE_selection to abc
            
            for (int t = 0; t < D; t++) {
                Mut[t] = x1[i][t] + F * (x1[a][t] - x1[b][t]);           //DE_Mutation
                if (Mut[t] >= up || Mut[t] <= down )  {
                    Mut[t] = x1[i][t];
                }
            }
            
            int rnbr = int (rand() )/ RAND_MAX * D;
            for (int j = 0; j < D; j++) {
                if ( float (rand()) / RAND_MAX <= CR || j == rnbr) {
                    trial[j] = Mut[j];
                }
                else                                                       //DE_Crossover
                {
                    trial[j] = x1[i][j];
                }
            }
            
            
            
            
            /*
             
             
             //      cout <<  evaluate(trial, f_number)  << "   " << evaluate(x1[i], f_number) << "  ||  ";
             
             for (int i = 0; i < NP; i++) {
             for (int j = 0; j < D; j++) {
             cout << x1[i][j] << " " ;
             }
             cout << "                                     ||   "<< evaluate(x1[i], f_number) <<  "    ";
             if (text == i) {
             cout << evaluate(trial, f_number);
             }
             cout << endl;
             }
             cout << endl << endl << endl;
             text ++;
             if (text == NP ) {
             text = 0;
             }
             
             */
            
            
            
            
            if (fabs (evaluate (trial, f_number)) <= fabs (evaluate (x1[i], f_number))) {
                for (int l = 0; l < D; l++) {
                    x1[i][l] = trial[l];
                }
                //     printf("%f \n",evaluate(trial) );
            }
            //DE_Select
            
            countDE1 ++;
            
            
        }
        
        for (int i = 0; i < NP; i++) {
            if (fabs (evaluate(x1[i], f_number) )<= deta) {
                condition1 = 1;
            }
        }
    }
    
    printf("第 %d 个测试函数", f_number);
    cout << endl;
    printf("DE = %d \n", countDE1);                                  //DE Algorithm(finished)
    
    while (countDE2 < gen_max && condition2 == 0) {                 //DE(BA) Algorithm(start)
        
        for (int  i = 0; i < NP; i++) {
            
            
            int a; int t1 = 0;
            do
            {   a = int (rand() % NP);
                t1++;                 // t1 t2 t3 prevent infinite loops caused by
                if (t1 > NP) {        //   inapporiate dots generated by BA network
                    break;
                }
            }  while (a == i );
            
            int b; int t2 = 0;
            do
            {   b = int (rand() % NP);
                t2++;
                if (t2 > NP) {
                    break;
                }
            }  while (b == i || b == a ||  G[a][b] == 0 );
            
            int c; int t3 = 0;
            do
            {   c = int (rand() % NP);
                t3++;
                if (t3 > NP) {
                    break;
                }
            }  while (c == i || c == a || c == b || G[a][c] == 0 ||    G[b][c] == 0 );
            
            
            
            
            for (int t = 0; t < D; t++) {
                Mut[t] = x2[a][t] + F * (x2[b][t] - x2[c][t]);           //DE_Mutation
                if (Mut[t] >= up || Mut[t] <= down )  {
                    Mut[t] = x2[a][t];
                }
            }
            
            int rnbr = int (rand() )/ RAND_MAX * D;
            for (int j = 0; j < D; j++) {
                if ( float (rand()) / RAND_MAX <= CR || j == rnbr) {
                    trial[j] = Mut[j];
                }
                else                                                       //DE_Crossover
                {
                    trial[j] = x2[i][j];
                }
            }
            
            if (fabs (evaluate (trial, f_number)) <= fabs (evaluate (x2[i] ,f_number))) {
                for (int l = 0; l < D; l++) {
                    x2[i][l] = trial[l];
                }
                //    printf("%f \n",evaluate(trial,f_number) );
            }
            //DE_Select
            
            
            
            
            countDE2 ++;
            
            
        }                                                           //DE(BA) Algorithm(finished)
        
        
        for (int i = 0; i < NP; i++) {
            if (fabs (evaluate(x2[i] ,f_number) )<= deta) {
                condition2 = 1;
            }
        }
        
        
    }
    printf("DE(BA) = %d \n", countDE2);
    printf("\n");
    
}


int main(int argc, const char * argv[]) {
    
    
    
    float evaluate(float *);
    int NP ;
    float F ;
    float CR ;
    int D ;
    long int gen_max = 100000;
    float up;
    float down;
    float deta;
    int m0;
    int m ;
    int t ;
    int f_number;
    //----------------------------------------------------------------------------------------------
    
    f_number = 1;
    NP = 15;
    F = 0.9;
    CR = 0.1;
    D = 3;
    up = 5.12;
    down = - 5.12;
    deta = 1.e-6;
    m0 = 4;
    m = 2;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //-----------------------------------------------------------------------------------------------
    
    f_number = 2;
    NP = 10;
    F = 0.9;
    CR = 0.9;
    D = 2;
    up = 2.048;
    down = - 2.048;
    deta = 1.e-6;
    m0 = 4;
    m = 5;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //-----------------------------------------------------------------------------------------------
    
    f_number = 4;
    NP = 10;
    F = 0.9;
    CR = 0;
    D = 30;
    up = 1.28;
    down = - 1.28;             //有问题
    deta = 15;
    m0 = 4;
    m = 5;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //-----------------------------------------------------------------------------------------------
    
    f_number = 6;
    NP = 10;
    F = 0.5;
    CR = 0;
    D = 4;
    up = 1000;
    down = - 1000;
    deta = 1.e-6;
    m0 = 4;
    m = 5;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //-----------------------------------------------------------------------------------------------
    
    f_number = 7;
    NP = 25;
    F = 0.5;
    CR = 0.2;
    D = 10;
    up = 400.0;
    down = - 400.0;
    deta = 1.e-6;
    m0 = 4;
    m = 3;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //-----------------------------------------------------------------------------------------------
    
    f_number = 11;
    NP = 20;
    F = 0.5;
    CR = 0.1;
    D = 30;
    up = 1.0;
    down = - 1.0;
    deta = 1.e-10;
    m0 = 4;
    m = 2;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    //----------------------------------------------------------------------------------------------
    
    f_number = 16;
    NP = 20;
    F = 0.5;
    CR = 0;
    D = 1;
    up = 10.0;
    down = - 10.0;
    deta = 1.e-6;
    m0 = 4;
    m = 2;
    t = NP - m0;
    
    model(NP, F, CR, D, gen_max, up, down, deta, m0, m, t, f_number);
    
    
    
    
    
    
    return 0;
}
