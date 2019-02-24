#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<iostream>
#include <iomanip>
using namespace std;
#define sum 2000    //��������������
#define a1 1.5 //���ٶ�1
#define a2 1.5//���ٶ�2
#define size  500   //������Ⱥ��ģ
#define Pmax 2 // �������ȡֵ
#define Pmin  -1*Pmax


typedef struct Chrom
{
	double place[2];
	double fit;
	double V;
}chrom;

chrom pop[size];
chrom per[size];
chrom gro;
double Vmax;
double Vmin;
double w;
double y(chrom pop)
{
	double y=0,x1,x2;
	double f5=0,F5=0;
	x1 = pop.place[0];
	x2 = pop.place[1];
	//y = 100 * (pop.place[1] * pop.place[1] - pop.place[0])*(pop.place[1] * pop.place[1] - pop.place[0]) + (1 - pop.place[1])*(1 - pop.place[1]);
	//y = (1 + (x1 + x2 + 1)*(x1+x2+1)*(19 - 14 * x1 + 3 * x1*x1 - 14 * x2 + 6 * x1*x2 + 3 * x2*x2))*(30 + (2 * x1 - 3 * x2)*(2*x1-3*x2)*(18 - 32 * x1 + 12 * x1*x1 + 48 * x2 - 36 * x1*x2 + 27 * x2*x2));
	//y = 0.5 + ((pow(sin(pow(x1*x1 + x2*x2, 0.5)), 2) - 0.5) / pow(1 + 0.001*(x1*x1 + x2*x2), 2));
	//y = pow(x1*x1 + x2*x2, 0.25)*(pow(sin(pow(50 * (x1*x1 + x2*x2), 0.1)), 2) + 1);
	//for (int i = 1; i <= 5; i++){	f5 = f5+i*cos((i + 1)*x1 + i);}for (int j = 1; j <= 5; j++){	F5 =F5+ j*cos((j + 1)*x2 + j);}y = f5*F5;
	y =  (4 - 2.1 * x1 * x1 + 1.0 / 3.0 * x1 * x1 * x1 * x1) * x1 * x1 + x1 * x2 + (-4 + 4 * x2 * x2) * x2 * x2;
	//y = x1*x1 + x2*x2;
	return y;
}
void *evpop(chrom pop[size])
{
	for (int i = 0; i<size; i++)
	{
		for (int j = 0; j<2; j++)
		{
			pop[i].place[j] = (((double)rand()) / RAND_MAX-0.5 ) * Pmax *2;  //Pmin��Pmax֮��������
			pop[i].V = (((double)rand()) / RAND_MAX - 0.5); //-Pmax��Pmax֮��  F3��-1��1֮��
		}
		pop[i].fit = y(pop[i]); //������Ӧ�Ⱥ���ֵ
	}
	for (int i = 0; i < size; i++)
	{
		per[i] = pop[i];
	}
	return 0;
}
int gpbest(chrom pop[size])
{
	int  best_fit = 0;
	long double min = pop[0].fit;
	for (int i = 1; i < size; i++)
	{
		if (pop[i].fit < min)
		{
			min = pop[i].fit;
		    best_fit = i;
		}
	}
	return best_fit;
}
void *chosegp(chrom pop[size])
{
	int bestpop;
	bestpop = gpbest(pop);
	gro = pop[bestpop];
	return 0;
}
void *fly(chrom pop[size])
{
	
		for (int j = 0; j < size; j++)
		{
			//�ٶȸ��¼����Ӹ���
			for (int k = 0; k<2; k++)
			{
				// �ٶȸ���
				double rand1 = (double)rand() / RAND_MAX; //0��1֮��������
				double rand2 = (double)rand() / RAND_MAX;
				pop[j].V = w*pop[j].V + a1 * rand1*(per[j].fit - pop[j].place[k]) + a2 * rand2*(gro.place[k] - pop[j].place[k]);
				if (pop[j].V > Vmax)pop[j].V = Vmax;
				if (pop[j].V < Vmin)pop[j].V = Vmin;
				// ���Ӹ���
				pop[j].place[k] = pop[j].place[k] + pop[j].V;
				if (pop[j].place[k] > Pmax)
					pop[j].place[k] = Pmax;
				if (pop[j].place[k] < Pmin)
					pop[j].place[k] = Pmin;
			}
			pop[j].fit = y(pop[j]); //�����ӵ���Ӧ��ֵ
		}
		return 0;
}
void *update(chrom pop[size])
{
		for (int j = 0; j<size; j++)
		{
			// ���弫ֵ����
			if (pop[j].fit < per[j].fit)
			{
				per[j] = pop[j];
			}
			// Ⱥ�弫ֵ����
			if (pop[j].fit < gro.fit)
			{
				gro = pop[j];
			}
		}
		return 0;
}
void main()
{   
	for (int j = 0; j < 10; j++)
	{   
		clock_t start, stop;
		start = clock();
		w = 0.9;
		Vmax = Pmax ;
		Vmin = Pmin ;
		srand((unsigned)time(NULL));
		evpop(pop);
		chosegp(pop);
		for (int i = 0; i < sum; i++)
		{
			fly(pop);
			update(pop);
		}
		cout << "��СֵΪ" << gro.fit << "    x1ֵΪ" << gro.place[0] << "   x2Ϊ" << gro.place[1] << endl;
		stop = clock();
		cout << "��ʱ  " << difftime(stop, start) << endl;
	}
}