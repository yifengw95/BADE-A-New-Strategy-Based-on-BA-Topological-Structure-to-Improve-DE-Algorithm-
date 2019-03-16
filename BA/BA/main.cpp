//
//  main.cpp
//  BA
//
//  Created by 王译烽 on 15/12/30.
//  Copyright © 2015年 王译烽. All rights reserved.
//

// BAnetwork.cpp : 定义控制台应用程序的入口点。

//#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define m0 4 //初始节点数
#define m 4 //每次引入的节点连接到m个结点上
#define t 46 //经过t步形成BA网络
//网络规模为m0+t


int G[t+m0][t+m0];//邻接矩阵
int k_degree[t+m0];//度
int k_sum;
double p_k[t+m0];//度分布
double c[t+m0];//聚类系数
double C_mean;//网络的聚类系数

void Initialize()
{
    int i,j;
    for (i=0;i<=m0-1;i++)
    {
        for (j=0;j<=m0-1;j++)
        {
            if (j!=i)
            {
                G[i][j]=1;
            }
            else
                continue;
        }
    }
}

void Create()//生成BA网络
{
    int i,j,ii,jj;
    int count;//连的第几条边，count<=m
    double p[t+m0];//与各个结点相连的概率，可重复使用
    double p1;//产生的随机数
    double p2;//累加的概率
    srand((unsigned)time(NULL));
    rand();
    
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

  //      for (ii=0;ii<=i-1;ii++)
  //      {
  //          printf("%d ",k_degree[ii]);
  //      }
   //     printf("\n");

        //DE
        int x1 = int (rand() / (double)RAND_MAX * i);
        int x1degree_sum = 0;
        for (ii = 0; ii <= i -1 ; ii++) {
            x1degree_sum = x1degree_sum + G[ii][x1];
        }
        int x2 = int (rand() / (double)RAND_MAX * (x1degree_sum));
        int x3 = int (rand() / (double)RAND_MAX * (x1degree_sum));
        int loc_x2 = 0;
        int loc_x3 = 0;
        
        int text_x1degree_sum = 0 ;
        for (ii = 0; ii <= i - 1; ii++) {
            text_x1degree_sum = text_x1degree_sum + G[ii][x1];
            if (text_x1degree_sum == x2) {
                loc_x2 = ii;
            }
            if (text_x1degree_sum == x3) {
                loc_x3 = ii;
            }
        }
        k_degree[x1] = int ((1/3) * (k_degree[x1] + k_degree[loc_x2] + k_degree[loc_x3]));
        
        
        
        k_sum=0;
        for (ii=0;ii<=i-1;ii++)
        {
            k_sum=k_degree[ii]+k_sum;
        }
        
        
        
        
        
        
        for (ii=0;ii<=i-1;ii++)
        {
            p[ii]=(double)k_degree[ii]/k_sum;
        }

        for (ii=0;ii<=i-1;ii++)
        {
            printf("%f ",p[ii]);
        }
        printf("\n");

        

        
        //开始连边
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
        G[j][i]=1;//连好了第一条边
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
        }//连好其余m-1条边
        
        
        /*
        float f0 = 0.5;
        int x1 = int (rand() / (double)RAND_MAX * i ) - 1;
//      int x2 = int (rand() / (double)RAND_MAX * i ) - 1;
//      int x3 = int (rand() / (double)RAND_MAX * i ) - 1;
        //double temp = p[x3];
        p[x3] = p[x3] + f0 * (p[x2] - p[x1]);
      
        */
    }
}


int  main( )
{
    int i,j;
    Initialize();

    for (i=0;i<=t+m0-1;i++)
    {
        for (j=0;j<=t+m0-1;j++)
        {
            printf("%d ",G[i][j]);
        }
        printf("\n");
    }

    Create();
    printf("====================================\n");
    for (i=0;i<=t+m0-1;i++)
    {
        for (j=0;j<=t+m0-1;j++)
        {
            printf("%d ",G[i][j]);
        }
        printf("\n");
    }
}
