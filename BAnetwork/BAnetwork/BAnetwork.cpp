// BAnetwork.cpp : �������̨Ӧ�ó������ڵ㡣

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define m0 4 //��ʼ�ڵ���
#define m 4 //ÿ������Ľڵ����ӵ�m�������
#define t 46 //����t���γ�BA����
//�����ģΪm0+t
#define DEBUG//�������е���
//#define  DEBUG1//��������е���

int G[t+m0][t+m0];//�ڽӾ���
int k_degree[t+m0];//��
int k_sum;
double p_k[t+m0];//�ȷֲ�
double c[t+m0];//����ϵ��
double C_mean;//����ľ���ϵ��

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

void Create()//����BA����
{
	int i,j,ii,jj;
	int count;//���ĵڼ����ߣ�count<=m
	double p[t+m0];//�������������ĸ��ʣ����ظ�ʹ��
	double p1;//�����������
	double p2;//�ۼӵĸ���
	srand((unsigned)time(NULL));
	rand();
	for (i=m0;i<=t+m0-1;i++)
	{
		//���ȼ���i-1�����Ķ��Լ��������ǵĸ���
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
		//��ʼ����
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
		G[j][i]=1;//�����˵�һ����
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
		}//��������m-1����
	}
}

void test()//������������Ķȷֲ�������ϵ��������֤����BA�������ȷ��
{
	int i,j,l;
	//����ȷֲ�
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
	//�������ϵ��
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
	fp=fopen("�ȷֲ�.txt","w+");
	for (i=0;i<=t+m0-1;i++)
	{
		fprintf(fp,"%d	%e  %d\n",k_degree[i],p_k[i],i);
	}
	FILE *fp2;
	fp2=fopen("����ϵ��.txt","w+");
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

