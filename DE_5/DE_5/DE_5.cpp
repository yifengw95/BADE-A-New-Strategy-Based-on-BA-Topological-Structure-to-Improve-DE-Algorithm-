 
//
//  Created by 王译烽 on 16/5/18.
//  Copyright © 2016年 王译烽. All rights reserved.
//

#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <malloc.h>
using namespace std;

const int NP = 100;
#define D_Max 30
#define REPEAT 10
const int m0 =10;//初始节点数
const int t = NP - m0;//经过t步形成BA网络
const int m = 6;//每次引入的节点连接到m个结点上

const double CR = 0.9;
const double F =0.5;

int G[1000][1000]; //Initialize BA Matrix
int Degree[1000];//the Degree of BA
long int gen_max;//最大演化代数
int Sum_Degree;//the Sum of the Degree of BA
double x1[NP][D_Max];
double x2[NP][D_Max];

int fall_1 = 0;
int fall_2 = 0;
int fall_3 = 0;
int repeat = 0;//测试函数循环次数
double Outcome_DE = 0;
double Outcome_DEBA = 0;


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
	
	if (f_number == 3)
	{
		for (int i = 1; i<= D; i++)
		{
			double sum = 0;
			for (int j = 1; j <= i; j++)
			{
				sum = sum + (*(p+j-1)) * (*(p+j-1));
			}
			outcome = outcome + sum;
		}
		
	}

	if (f_number == 4)
	{
		double max = -10;
		for (int i = 0; i < D ;i++)
		{
			if (max < (fabs(*(p+i))))
			{
				max = fabs(*(p+i));
			}
		}
		outcome = max;
	}
	
	if (f_number == 5)
	{
		for (int i = 1; i<= D-1; i++)
		{
			outcome = outcome + 100 *((*(p+i)) - (*(p+i-1))*(*(p+i-1)) )*((*(p+i)) - (*(p+i-1))*(*(p+i-1)) ) + ((*(p+i-1))-1 )*((*(p+i-1))-1 );
		}
		
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
		for (int i = 1; i <= D; i++)
		{
			outcome = outcome + (i)* (pow(*(p+i-1),4)) ;
		}
		outcome = outcome + (double(rand() ) / RAND_MAX  );
	}
	
	if (f_number ==8)
	{
		for (int i = 1; i <= D; i++)
		{
			outcome = outcome + (-( (*(p+i-1)) * sin(sqrt(fabs(*(p+i-1))))  )); //+ D*418.98288727243369);
		}
		outcome = outcome +  D*418.98288727243369;
	}
	
	if (f_number == 9)
	{
		for (int i = 1; i <= D; i++)
		{
			outcome = outcome + ((*(p+i-1))*(*(p+i-1)) - 10*cos(2 * 3.1415926535 * (*(p+i-1))) + 10  );
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

	if (f_number == 12)
	{
		double u (double xi,double a,double k,double m);
		double y[40];
		for (int i = 0; i < 30; i++)
		{
			y[i+1] = 1 + 1.0 / 4.0 * ((*(p+i)) + 1);
		}
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 1; i <= D-1; i++)
		{
			sum1 = sum1 + (y[i] - 1)* (y[i] - 1)  * (1 + (10* sin(3.1415926*y[i+1]) * sin(3.1415926*y[i+1]) ) );
		}
		for (int i = 1; i<= D; i++)
		{
			sum2 = sum2 + u((*(p+i-1)), 10, 100, 4);
		}
		outcome = 3.1415926 / D *(10 * sin(3.1415926*y[1]) * sin(3.1415926*y[1]) + sum1 + (y[30]-1) * (y[30]-1) ) + sum2;   
	}

	if (f_number == 13)
	{
		double u (double xi,double a,double k,double m);
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 1; i <= D-1; i++)
		{
			sum1 = sum1 + ((*(p+i-1))-1) * ((*(p+i-1))-1)* ( 1+ sin(3*3.1415926 * (*(p+i)))*sin(3*3.1415926 * (*(p+i)))  );
		}
		for (int i = 1; i <= D; i++)
		{
			sum2 = sum2 + u((*(p+i-1)), 5, 100, 4);
		}
		outcome = 0.1*( sin(3*3.1415926* (*(p)))*sin(3*3.1415926* (*(p))) + sum1 + ((*(p+D-1)) - 1)*((*(p+D-1)) - 1)*(1 + sin(2*3.1415926*(*(p+D-1)))*sin(2*3.1415926*(*(p+D-1)))  )   )+ sum2;
		
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

double u (double xi,double a,double k,double m){
	double output = 0;
	if (xi > a)
	{
		output = k * (pow((xi-a) , m));
	}
	if (-a <= xi && xi <= a)
	{
		output = 0;
	}
	if (xi < -a)
	{
		output = k * (pow((-xi-a) , m));
	}
	
	return output;
}

void model (
			int f_number,
			int D,
            long int gen_max ,
            double up,//测试函数自变量取值上界
            double down,//测试函数自变量取值下界
			double true_outcome//测试函数结果真值
            ){
    
	repeat = 0;
	Outcome_DE = 0;
	Outcome_DEBA = 0;

	srand((unsigned)time(NULL));
	rand();

	for(;repeat < REPEAT; repeat++){
   
		double evaluate(double *, int,int);// 测试函数集
    
		double Mut[D_Max];
		double trial[D_Max];
		int countDE1 = 1;//DE运行代数
		int countDE2 = 1;//DE(BA）运行代数
		
    
		//全局变量初始化
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
		
		 //DE Algorithm(start)
		while (countDE1 < gen_max) {                      
        
			for (int  i = 0; i < NP; i++) {
				//DE_selection to abc

				double Best_Deta_1 = fabs (evaluate (x1[1], f_number,D) - true_outcome);
				
				int a = 1;

				for (int temp1 = 0; temp1 < NP; temp1++)
				{
					if ( Best_Deta_1 >= fabs (evaluate (x1[temp1], f_number,D) - true_outcome)) 
					{
						Best_Deta_1 = fabs (evaluate (x1[temp1], f_number,D) - true_outcome);
						a = temp1;
					}
				}
				
				
				//int a; do a = int (rand() % NP); while (a == i);
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
			//每两百代输出一次数据
			/*
			if (countDE1 % 200 == 0)
			{
				double Temp_DE1 = 0;
				double Temp_Sum1 = 0;
				for (int i = 0; i < NP; i++)
				{
					Temp_Sum1 = Temp_Sum1 + evaluate(x1[i],f_number,D);
				}
				Temp_DE1 = double(Temp_Sum1 / NP);
				cout << "countDE1   "<< countDE1 << "    "<< Temp_DE1 << endl;
			}
		
			*/
		}
	
		//cout << endl;
	

		double Average_DE1 = evaluate(x1[0],f_number,D);
		double temp_deta = fabs (evaluate (x1[1], f_number,D) - true_outcome);
		for (int i = 0; i < NP; i++ )
		{
			if ( temp_deta >= fabs (evaluate (x1[i], f_number,D) - true_outcome))
			{
				temp_deta = fabs (evaluate (x1[i], f_number,D) - true_outcome);
				Average_DE1 = evaluate(x1[i],f_number,D);
			}
		}


	 cout << "   " << repeat << "   "<< "DE   "  << temp_deta<<"         //        "; 

		Outcome_DE = Outcome_DE + temp_deta;

	
		//DE Algorithm(finished)


		//DE(BA) Algorithm(start)
		int part_fall_1 = 0;
		int part_fall_2 = 0;
		int part_fall_3 = 0;

		while (countDE2 < gen_max ) {                 
			for (int  i = 0; i < NP; i++) {
				// t1 t2 t3 prevent infinite loops caused by inappropriate dots generated by BA network
				
				double Best_Deta_2 = fabs (evaluate (x2[1], f_number,D) - true_outcome);
				
				int a = 1;

				for (int temp1 = 0; temp1 < NP; temp1++)
				{
					if (G[temp1][i] == 1 && ( Best_Deta_2 >= fabs (evaluate (x2[temp1], f_number,D) - true_outcome)) )
					{
						Best_Deta_2 = fabs (evaluate (x2[temp1], f_number,D) - true_outcome);
						a = temp1;
					}
				}
				
				/*
				
				
				int a; int t1 = 0;
				do
				{   a = int (rand() % NP);
					t1++;                 
					if (t1 > D*NP) { 
						part_fall_1++;
					   // cout << endl << "1FAIL!!" << endl;
						break;
					}
				} while (a == i || G[a][i] == 0);
            */
				int b; int t2 = 0;
				do
				{   b = int (rand() % NP);
					t2++;
					if (t2 > D*NP) {
						part_fall_2++;
						//cout << endl << "2FAIL!!" << endl;
						break;
					}
				}   while (b == i || b == a  || G[a][b] == 0 );

				int c; int t3 = 0;
				do
				{   c = int (rand() % NP);
					t3++;
					if (t3 > D*NP) {
						//part_fall_1++;
						part_fall_3++;
						//cout << endl << "2FAIL!!" << endl;
						break;
					}
				}   while (c == i || c == a  || c == b || G[a][c] == 0 );
				/*
				int d; int t4 = 0;
				do
				{   d = int (rand() % NP);
				t4++;
				if (t4 > D*NP) {
					part_fall_2++;
					break;
				}
				}   while (d == i || d == a  || d == b || d == c || G[a][d] == 0 );

				int e; int t5 = 0;
				do
				{   e = int (rand() % NP);
				t5++;
				if (t5 > D*NP) {
					part_fall_3++;
					break;
				}
				}   while (e == i || e == a  || e == b || e == c || e == d || G[a][e] == 0 );
				*/

				//DE(BA)_Mutation
				for (int t = 0; t < D; t++) {
					Mut[t] = x2[a][t] + F * (x2[b][t] - x2[c][t])  ;           
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
			//每两百代输出一次数据
			/*
			if (countDE2 % 200 == 0)
			{
				double Temp_DE2 = 0;
				double Temp_Sum2 = 0;
				for (int i = 0; i < NP; i++)
				{
					Temp_Sum2 = Temp_Sum2 + evaluate(x2[i],f_number,D);
				}
				Temp_DE2 = double(Temp_Sum2 / NP);
				cout << "countDE2   "<< countDE2 << "    "<< Temp_DE2 << endl;
			}
			*/
			//每1000代作为筛选条件
			/*
			if ((countDE2 % 500) == 0)
			{
				//输出 a， b， c， i 点的平均度
				cout << "Ave_Degree_A      "<< double(Degree_A) / (countDE2 - 1) / NP<<"     Ave_D     "<<Average_Degree<<endl;
				cout << "Ave_Degree_B      "<< double(Degree_B) / (countDE2 - 1) / NP<<"     Ave_D     "<<Average_Degree<<endl;
				cout << "Ave_Degree_C      "<< double(Degree_C) / (countDE2 - 1) / NP<<"     Ave_D     "<<Average_Degree<<endl;
				cout << "Ave_Degree_I      "<< double(Degree_I) / (countDE2 - 1) / NP<<"     Ave_D     "<<Average_Degree<<endl;

				//输出最优的点的度和其自身邻居的度
				double Min_Text = 1000;
				int Min_Loc = 0;
				for (int i = 0; i < NP; i++)
				{
					if ((fabs (evaluate (x2[i], f_number,D) - true_outcome)) <= Min_Text)
					{
						Min_Text = fabs (evaluate (x2[i], f_number,D) - true_outcome);
						Min_Loc = i;
					}
				}
				cout << "最优点的度       "<< Degree[Min_Loc]<<endl;
				for (int i =0 ;i <= t+m0-1; i++)
				{
					if (G[Min_Loc][i] == 1)
					{
						cout << Degree[i] << "      ";
					}
				
				}
			
				cout << endl;
			
			}*/
		
		}

		double Average_DE2 = evaluate(x2[0],f_number,D);
		temp_deta = fabs (evaluate (x1[1], f_number,D) - true_outcome);
		for (int i = 0; i < NP; i++ )
		{
			if ( temp_deta >= fabs (evaluate (x2[i], f_number,D) - true_outcome))
			{
				temp_deta = fabs (evaluate (x2[i], f_number,D) - true_outcome);
				Average_DE2 = evaluate(x2[i],f_number,D);
			}
		}

		cout << "   " << repeat << "   "<<"DEBA   " <<temp_deta<< endl; 


		Outcome_DEBA = Outcome_DEBA + temp_deta;
	
		fall_1 = fall_1 + part_fall_1;
		fall_2 = fall_2 + part_fall_2;
		fall_3 = fall_3 + part_fall_3;
		//DE(BA) Algorithm(finished)   
	
	}

	Outcome_DE = Outcome_DE / REPEAT;
	Outcome_DEBA = Outcome_DEBA / REPEAT;
	fall_1 = fall_1 / REPEAT ;
	fall_2 = fall_2 / REPEAT;
	fall_3 = fall_3 / REPEAT;

	printf("第 %d 个测试函数(DE)", f_number);
	cout << endl;
	cout << Outcome_DE << endl;

	printf("第 %d 个测试函数(DE_BA)", f_number);
	cout << endl;
	cout << Outcome_DEBA << endl;

	

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

	cout<<setiosflags(ios::scientific)<<setprecision(8);

//------------------------------------------------------------------------------------------------------------------	
	
	//BA_Network Initialize
	int k_degree[t+m0];//度
	int k_sum;
	memset(Degree,0,sizeof(Degree));
	memset(G,0,sizeof(G));
	Sum_Degree = 0;

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

	for (int i = 0 ;i <= t+m0-1; i++)
	{
		for (int j =0 ;j <= t+m0-1; j++)
		{
			Degree[i] = Degree[i] + G[i][j];
		}

	}
	for (int i = 0; i <= t+m0-1; i++)
	{
		Sum_Degree = Sum_Degree + Degree[i];
		G[i][i] = 1;
	}

	double Average_Degree = 0;
	Average_Degree = double(Sum_Degree) / 100;

	// BA_Network Initialize(finished)

	//存储矩阵信息
	FILE* fp;
	fp = fopen("矩阵信息.txt","w+");
	for (int ii = 0;ii<=t+m0-1;ii++)
	{
		for (int jj =0; jj <= t+m0-1; jj++)
		{
			fprintf(fp,"%d ",G[ii][jj]);
		}
		fprintf(fp,"\n");
	}
	fclose (fp);
	
//----------------------------------------------------------------------------------------------------------------
	/*
	FILE* fp1;
	fp1 = fopen("矩阵信息.txt","r");
	for (int ii = 0; ii <= t+m0-1; ii++)
	{
		for (int jj =0 ; jj <= t+m0-1; jj++)
		{
			fscanf(fp1,"%d", &G[ii][jj]);
		}
		
	}
	fclose(fp1);

	FILE* fp;
	fp = fopen("矩阵信息(验证).txt","w+");
	for (int ii = 0;ii<=t+m0-1;ii++)
	{
		for (int jj =0; jj <= t+m0-1; jj++)
		{
			fprintf(fp,"%d ",G[ii][jj]);
		}
		fprintf(fp,"\n");
	}
	fclose (fp);
	*/
    //----------------------------------------------------------------------------------------------
    /*
	f_number = 1;
    D = 25;
	gen_max = 500;
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
	gen_max = 500;
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
	
	f_number = 3;
	D = 30;
	gen_max = 500;
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
	
	//-------------------------------------------------------------------------------------------------
	
	f_number = 4;
	D = 30;
	gen_max =5000;
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
	*/
	//-------------------------------------------------------------------------------------------------
	/*
	f_number = 5;
	D = 30;
	gen_max =3000;
	up = 30;
	down = -30;
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
	gen_max = 100;
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
	*/
	//-------------------------------------------------------------------------------------------------
	/*
	f_number = 8;
	D = 30;
	gen_max = 9000;
	up = 500;
	down = -500;
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
	gen_max = 2000;
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
	*/
	//------------------------------------------------------------------------------------------------
	
	f_number = 11;
	D = 30;
	gen_max = 3000;
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
	
	f_number = 12;
	D = 30;
	gen_max = 1500;
	up = 50;
	down = -50;
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

	f_number = 13;
	D = 30;
	gen_max = 1500;
	up = 50;
	down = -50;
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
