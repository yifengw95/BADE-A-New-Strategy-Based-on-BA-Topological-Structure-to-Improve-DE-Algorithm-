// BAnetwork.cpp : 定义控制台应用程序的入口点。

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define m0 4 //初始节点数
#define m 4 //每次引入的节点连接到m个结点上
#define t 46 //经过t步形成BA网络
//网络规模为m0+t
#define DEBUG//主程序中调试
//#define  DEBUG1//其余程序中调试

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
#ifdef DEBUG1
		for (ii=0;ii<=i-1;ii++)
		{
			printf("%d ",k_degree[ii]);
		}
		printf("\n");
#endif
		k_sum=0;
		for (ii=0;ii<=i-1;ii++)
		{
			k_sum=k_degree[ii]+k_sum;
		}
		for (ii=0;ii<=i-1;ii++)
		{
			p[ii]=(double)k_degree[ii]/k_sum;
		}
#ifdef DEBUG1
		for (ii=0;ii<=i-1;ii++)
		{
			printf("%f ",p[ii]);
		}
		printf("\n");
#endif
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
	}
}

void test()//计算生成网络的度分布及聚类系数，以验证生成BA网络的正确性
{
	int i,j,l;
	//计算度分布
	for (i=0;i<=t+m0-1;i++)
	{
		k_degree[i]=0;
	}
	for (i=0;i<=t+m0-1;i++)
	{
		for (j=0;j<=t+m0-1;j++)
		{
			k_degree[i]=k_degree[i]+G[i][j];
		}
	}
	for(i=0;i<=t+m0-1;i++)
		p_k[k_degree[i]]=p_k[k_degree[i]]+(double)1/(t+m0);
	//计算聚类系数
	for(i=0;i<=t+m0-1;i++)
	{
		for(j=0;j<=t+m0-1;j++)
			for(l=0;l<=t+m0-1;l++)
				c[i]=c[i]+G[i][j]*G[j][l]*G[l][i];
		c[i]=c[i]/(k_degree[i]*(k_degree[i]-1));
	}
	for (i=0;i<=t+m0-1;i++)
	{
		C_mean=C_mean+c[i];
	}
	C_mean=C_mean/(t+m0);
}

void main()
{
	int i,j;
	Initialize();
#ifdef DEBUG
	for (i=0;i<=t+m0-1;i++)
	{
		for (j=0;j<=t+m0-1;j++)
		{
			printf("%d ",G[i][j]);
		}
		printf("\n");
	}
#endif
	Create();
#ifdef DEBUG
	printf("====================================\n");
	for (i=0;i<=t+m0-1;i++)
	{
		for (j=0;j<=t+m0-1;j++)
		{
			printf("%d ",G[i][j]);
		}
		printf("\n");
	}
#endif
	test();
	FILE *fp;
	fp=fopen("度分布.txt","w+");
	for (i=0;i<=t+m0-1;i++)
	{
		fprintf(fp,"%d	%e  %d\n",k_degree[i],p_k[i],i);
	}
	FILE *fp2;
	fp2=fopen("聚类系数.txt","w+");
	fprintf(fp2,"%f",C_mean);
	FILE *fp1;
	fp1=fopen("BAnetwork.txt","w+");
	for (i=0;i<=t+m0-1;i++)
	{
		for (j=0;j<=t+m0-1;j++)
		{
			fprintf(fp1,"%d	",G[i][j]);
		}
		fprintf(fp1,"\n");
	}	
	fclose(fp1);
	fclose(fp);


	FILE *fp3;
	fp3=fopen("BA.txt","w+");
	for (i=0;i<=t+m0-1;i++)
	{
		for (j=0;j<=t+m0-1;j++)
		{
			fprintf(fp3,"%d",G[i][j]);
		}
		fprintf(fp3,"\n");
	}
	int sum=0,kkk;
	for(kkk=0;kkk<t+m0;kkk++){
		sum+=k_degree[kkk];
	}
	printf("%d",sum);
	system("pause");
}

