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
#include <malloc.h>
using namespace std;
#define D_Max 30

const int NP = 100;
const double CR = 0.9;
const double F =0.5;
int G[1000][1000]; //Initialize BA Matrix
double x1[NP][D_Max];
double x2[NP][D_Max];
const int m0 =10;//初始节点数
const int t = NP - m0;//经过t步形成BA网络
const int m = 6;//每次引入的节点连接到m个结点上
long int gen_max;//最大演化代数

int fall_1 = 0;
int fall_2 = 0;
int fall_3 = 0;


// 测试函数集
double evaluate (double *p, int f_number,int D ){             
    double outcome = 0;
    
    if (f_number == 1) {
        for (int i = 0; i < D; i++)
        {
			outcome = outcome + (*(p+i)) * (*(p+i));
        }
    }
    
	if (f_number == 2)
	{
		double sum1 = 0;
		double sum2 = 1;
		for (int i = 0; i < D;i++)
		{
			sum1 = sum1 + fabs(*(p+i));
			sum2 = sum2 * fabs(*(p+i));
		}
		outcome = sum1 + sum2;
	}
	
	if (f_number == 6)
	{
		for (int i = 0; i < D; i++)
		{
			outcome = outcome + fabs(*(p+i) + 0.5) * fabs(*(p+i) + 0.5);
		}
	}
	
	if (f_number == 7)
	{
		srand((unsigned)time(NULL));
		rand();
		for (int i = 0; i < D; i++)
		{
			outcome = outcome + (i+1)* pow(*(p+i),4);
		}
		outcome = outcome + (double(rand() ) / RAND_MAX  );
	}
	
	if (f_number ==8)
	{
		for (int i = 0; i < D; i++)
		{
			outcome = outcome + (-( (*(p+i)) * sin(sqrt(fabs(*(p+i))))  ));
		}
	}
	
	if (f_number == 9)
	{
		for (int i = 0; i < D; i++)
		{
			outcome = outcome + ( (*(p+i))*(*(p+i)) - 10*cos(2 * 3.1415926535898 * (*(p+i))) + 10  );
		}	
	}
	
	if (f_number == 10)
	{
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < D; i++)
		{
			sum1 = sum1 + ( 1.0 / 30.0 * (*(p+i))*(*(p+i))  );
			sum2 = sum2 + ( 1.0 / 30.0 * cos(2 * 3.1415926535898 * (*(p+i)))  );
		}
		outcome = -20.0 * exp(-0.2 * sqrt(sum1)) - exp(sum2) + 20 +  2.718281828459;
	}
	
	if (f_number == 11)
	{
		double sum1 = 0;
		double sum2 = 1;
		for (int i = 0; i < D;i++)
		{
			sum1 = sum1 + (*(p+i)) * (*(p+i));
			sum2 = sum2 * ( cos( (*(p+i)) / sqrt(i+1) ) );
		}
		outcome =  1.0 /4000.0 * sum1 - sum2 + 1;
		
	}
	
	if (f_number == 16)
	{
		outcome = 4.0 * pow(*(p),2) - 2.1 * pow(*(p),4) + 1.0/3.0 * pow(*(p),6) + (*(p) * *(p+1)) - 4.0 * pow(*(p+1),2) + 4.0 * pow(*(p+1),4);
	}
	
	if (f_number == 18)
	{
		double x1 = *p;
		double x2 = *(p+1);
		outcome = ( 1+ (x1+x2+1)*(x1+x2+1)*(19-14*x1+3*x2*x2-14*x2+6*x1*x2+3*x2*x2) ) * ( 30+ (2*x1-3*x2)*(2*x1-3*x2)*(18-32*x1+12*x1*x1+48*x2-36*x1*x2+27*x2*x2)  );
	}
	



   return outcome;
}

void model (
			int f_number,
			int D,
            long int gen_max ,
            double up,//测试函数自变量取值上界
            double down,//测试函数自变量取值下界
			double true_outcome//测试函数结果真值
            ){
    
    srand((unsigned)time(NULL));
    rand();
    double evaluate(double *, int,int);// 测试函数集
    
    double Mut[D_Max];
    double trial[D_Max];
    int countDE1 = 1;//DE运行代数
    int countDE2 = 1;//DE(BA）运行代数
    
	//全局变量初始化
	memset(G,0,sizeof(G));
	memset(x1,0,sizeof(x1));
	memset(x2,0,sizeof(x2));

	//DE_Matrix Initialize
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < D; j++) {
            x1[i][j] = double(rand() ) / RAND_MAX * (up - down) + down;
        }                                                               
    }
    
	 //DE(BA)_Matrix Initialize
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < D; j++) {
            x2[i][j] = double(rand() ) / RAND_MAX * (up - down) + down;
        }                                                              
    }
    
    int k_degree[t+m0];//度
    int k_sum;
    
	 //BA_Network Initialize
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
    int ii,jj;                                                
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
    
    // BA_Network Initialize(finished)
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
	
	 //DE Algorithm(start)
    while (countDE1 < gen_max) {                      
        
        for (int  i = 0; i < NP; i++) {
            //DE_selection to abc
            int a; do a = int (rand() % NP); while (a == i);
            int b; do b = int (rand() % NP); while (b == i || b == a);
            int c; do c = int (rand() % NP); while (c == i || c == a || c == b);
            //DE_Mutation
            for (int t = 0; t < D; t++) {
                Mut[t] = x1[a][t] + F * (x1[b][t] - x1[c][t]);           
                if (Mut[t] >= up || Mut[t] <= down )  {
                    Mut[t] = x1[a][t];
                }
            }
            
            int j = rand() % D;
            //DE_Crossover
            for (int k = 0; k < D; k++) {
                if ( float (rand()) / RAND_MAX < CR || k == j) {
                    trial[k] = Mut[k];
                }
                else                                                     
                {
                    trial[k] = x1[i][k];
                }
                j = (j + 1) % D;
            }
			
        /*
             for (int it = 0; it < NP; it++) {
				 for (int jt = 0; jt < D; jt++) {
				   // cout << x1[it][jt] << "   " ;
			      }
				// cout << "=" ;
			//	 cout << evaluate(x1[it], f_number,D) <<  " ";
                 
             }
             cout << "   ||    " <<  evaluate(trial, f_number,D) << endl << endl;;
         */
	
            //DE_Select
			if (fabs (evaluate (trial, f_number,D) - true_outcome) <= fabs (evaluate (x1[i], f_number,D) - true_outcome)) {
                for (int l = 0; l < D; l++) {
                    x1[i][l] = trial[l];
                }
                //     printf("%f \n",evaluate(trial) );
            }
       
                
        }

		countDE1 ++; 
    }
	double Average_DE1;
	double Sum_DE1 = 0;
	for (int i = 0; i < NP; i++ )
	{
		Sum_DE1 = Sum_DE1 + evaluate(x1[i],f_number,D);
	}
	Average_DE1 = double(Sum_DE1 / NP);
    
    printf("第 %d 个测试函数(DE)", f_number);
    cout << endl;
	cout << Average_DE1 << endl;
	
	//DE Algorithm(finished)
	
	//DE(BA) Algorithm(start)
    while (countDE2 < gen_max ) {                 
        for (int  i = 0; i < NP; i++) {
            
            // t1 t2 t3 prevent infinite loops caused by inappropriate dots generated by BA network
            int a; int t1 = 0;
            do
            {   a = int (rand() % NP);
                t1++;                 
                if (t1 > 30 * NP) { 
					fall_1++;
                   // cout << endl << "1FAIL!!" << endl;
					break;
                }
            } // while (a == i );
				while (a == i || G[a][i] == 0);
            
            int b; int t2 = 0;
            do
            {   b = int (rand() % NP);
                t2++;
                if (t2 > 30 * NP) {
					fall_2++;
					//cout << endl << "2FAIL!!" << endl;
                    break;
                }
            } // while (b == i || b == a ||  G[a][b] == 0 || G[b][i] == 0);
			   while (b == i || b == a || G[a][b] == 0 || G[b][i] == 0);

            int c; int t3 = 0;
            do
            {   c = int (rand() % NP);
                t3++;
                if (t3 > 30 * NP) {
					fall_3++;
					//cout << endl << "3FAIL!!" << endl;
                    break;
                }
            }  //while (c == i || c == a || c == b || G[a][c] == 0 || G[b][c] == 0 );
               while (c == i || c == a || c == b || G[a][c] == 0 || G[b][c] == 0 || G[c][i] == 0 );

            //DE(BA)_Mutation
            for (int t = 0; t < D; t++) {
                Mut[t] = x2[a][t] + F * (x2[b][t] - x2[c][t]);           
                if (Mut[t] >= up || Mut[t] <= down )  {
                    Mut[t] = x2[a][t];
                }
            }
        
            int j = rand() % D;
            //DE(BA)_Crossover
            for (int k = 0; k < D; k++) {
                if ( float (rand()) / RAND_MAX < CR || k == j) {
                    trial[k] = Mut[k];
                }
                else                                                     
                {
                    trial[k] = x2[i][k];
                }
                j = (j + 1) % D;
            }
            
            //DE(BA)_Select
			if (fabs (evaluate (trial, f_number,D) - true_outcome) <= fabs (evaluate (x2[i], f_number,D) - true_outcome)) {
				for (int l = 0; l < D; l++) {
					x2[i][l] = trial[l];
				}
			}  
			/* 
             for (int it = 0; it < NP; it++) {
				 for (int jt = 0; jt < D; jt++) {
				   // cout << x1[it][jt] << "   " ;
			      }
				// cout << "=" ;
			//	 cout << evaluate(x2[it], f_number,D) <<  " ";
                 
             }
             cout << "   ||    " <<  evaluate(trial, f_number,D) << endl << endl;;
         */
        }      

		countDE2 ++;
    }
	double Average_DE2;
	double Sum_DE2 = 0;
	for (int i = 0; i < NP; i++ )
	{
		Sum_DE2 = Sum_DE2 + evaluate(x2[i],f_number,D);
	}
	Average_DE2 = Sum_DE2 / NP;
	printf("第 %d 个测试函数(DE_BA)", f_number);
	cout << endl;
	cout << Average_DE2 << endl;
    //DE(BA) Algorithm(finished)   

	
}


int main(int argc, const char * argv[]) {
    
	int f_number;
    int D ;
    //long int gen_max ;
    double up;
    double down;
    double true_outcome;

	void model (
		int f_number,
		int D,
		long int gen_max ,
		double up,//测试函数自变量取值上界
		double down,//测试函数自变量取值下界
		double true_outcome//测试函数结果真值
		);

	double evaluate (double *p, int f_number ,int D);

	double text_1;
	double text_2;
	double text_3;
	
    //----------------------------------------------------------------------------------------------
    
	f_number = 1;
    D = 30;
	gen_max = 1500;
    up = 100;
    down = - 100;
    true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

    model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
    cout << endl ;

   //--------------------------------------------------------------------------------------------------
	
	f_number = 2;
	D = 30;
	gen_max = 2000;
	up = 10;
	down = -10;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//-------------------------------------------------------------------------------------------------
	
	f_number = 6;
	D = 30;
	gen_max = 1500;
	up = 100;
	down = -100;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//--------------------------------------------------------------------------------------------------
	
	f_number = 7;
	D = 30;
	gen_max = 3000;
	up = 1.28;
	down = -1.28;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//-------------------------------------------------------------------------------------------------
	/*
	f_number = 8;
	D = 30;
	gen_max = 9000;
	up = 500;
	down = -500;
	true_outcome = -12569.5;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	*/
	//-------------------------------------------------------------------------------------------------
	
	f_number = 9;
	D = 30;
	gen_max = 5000;
	up = 5.12;
	down = -5.12;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//-------------------------------------------------------------------------------------------------
	
	f_number = 10;
	D = 30;
	gen_max = 1500;
	up = 32;
	down =-32;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//------------------------------------------------------------------------------------------------
	
	f_number = 11;
	D = 30;
	gen_max = 2000;
	up = 600;
	down = -600;
	true_outcome = 0;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;

	//-------------------------------------------------------------------------------------------------
	/*
	f_number = 16;
	D = 2;
	gen_max = 100;
	up = 5;
	down = -5;
	true_outcome = -1.0316285;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	
	//-------------------------------------------------------------------------------------------------
	
	f_number = 18;
	D = 2;
	gen_max = 100;
	up = 2;
	down = -2;
	true_outcome = 3;

	fall_1 = 0;
	fall_2 = 0;
	fall_3 = 0;

	model(f_number,D,gen_max,up,down,true_outcome);

	text_1 = double(fall_1) / gen_max / NP;
	text_2 = double(fall_2) / gen_max / NP;
	text_3 = double(fall_3) / gen_max / NP;
	cout << "fall_1      "<<text_1<<endl;
	cout << "fall_2      "<<text_2<<endl;
	cout << "fall_3      "<<text_3<<endl;
	cout << endl ;
	*/
	//-------------------------------------------------------------------------------------------------



	

    getchar();
	getchar();
    
    return 0;
}
