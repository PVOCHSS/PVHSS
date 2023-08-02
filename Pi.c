#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <malloc.h>
#include <gmp.h>
#include <string.h>
#include <flint/fq.h>
#include <flint/fq_mat.h>
#include <flint/fq_vec.h>
#include <flint/fq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_mpoly.h>

int n,t,d,m;
ulong N;

fmpz_t p;

fq_ctx_t Fp;

fq_t * SID; // set of ID of servers

flint_rand_t state;

double * test_time;


bool DEBUG;

// ******* Initialize public parameters ********
void init(int nn, int dd){

	// DEBUG = true;
	DEBUG = false;
	n = nn;           // length of data
	t = 3;            // privacy and verifiability threshold (t=v)
	d = dd;            // degree of polynomial f
	m = (d+1)*t+1;    // number of servers
    flint_randinit(state);

	
	fmpz_init(p);
    fmpz_set_str(p,"340282366920938463463374607431768211297",10);  // 128 bit prime 
    // fmpz_set_str(p,"7",10);
    fq_ctx_init(Fp,p,1,"gen");
}



// ************ Share algorithm ********************
void Pi_Share(fq_mat_t x, fq_mat_t * s){
	// x[1..n];
	// s[1..m][1..n];

	fq_t temp;
	fq_init(temp,Fp);

	/*
	  Generating phi(u)

	  phi(u) is a degree-t polynomial with phi(0) = x
	  Here, we split phi(u)= [phi_0(u),...,phi_{n-1}(u)]
	  phi_i(0) = x_i
	*/

	clock_t test_start,test_end;

			
			
			
	fq_poly_t phi[n];

	// for(int rp=0;rp<10;++rp){
	// test_start = clock();
	for (int i = 0; i < n; i ++){
		fq_poly_init(phi[i],Fp);

		fq_poly_set_coeff(phi[i],0,fq_mat_entry(x,i,0),Fp);
		for (int k = 1; k <= t; k ++){
			fq_randtest(temp, state, Fp);
			
			fq_poly_set_coeff(phi[i], k, temp, Fp);
		}	
	}
	// test_end = clock();
	// test_time[rp] = (double) (test_end- test_start)/CLOCKS_PER_SEC;
	// }

	/*
	  Computing s = [s_1,...s_m]	 
	  s_j = phi(j) 
	*/
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < m; j ++){
            fq_poly_evaluate_fq(fq_mat_entry(s[j],i,0),phi[i],SID[j],Fp); //gaidong
        }
    }  


	if (DEBUG) {
		printf("**Debuging in L_share**:\n");
		printf("x = ");
		fq_mat_print_pretty(x, Fp);
		printf("\n");
		for (int i = 0; i < n; i ++){
			printf("phi_%d=",i);
			fq_poly_print_pretty(phi[i],"X",Fp);
			printf("\n");
		}
	}


	for (int i = 0; i < n; i ++){
		fq_poly_clear(phi[i], Fp);
	}
	fq_clear(temp,Fp);
}


// ************ Eval algorithm *********************

void next(int _length, int a[]){
	a[_length - 1] += 1;
	for (int i = _length - 1; i >= 1; i --){
		if (a[i] >= n) {
			a[i - 1] += 1;
			for (int j = i; j <= _length - 1; j ++){
				a[j] = a[i - 1];
			}
		}
	}
}

// eval an n-variate degree-d function f(x), output y=f(x)
void Eval(fq_mat_t f, fq_mat_t x, fq_t y){

	int index=0; //index of f
	fq_t y1;
   	fq_init(y1,Fp);

	int a[d]; 
	for (int i = 0; i < d; i ++){
		a[i] = 0;
	}

	while ( (a[0] < n) & (index < N) ){
		fq_set(y1,fq_mat_entry(f,index,0),Fp);
		for (int j = 0; j < d; j ++){
			fq_mul(y1,y1,fq_mat_entry(x,a[j],0),Fp);
		}
		fq_add(y,y,y1,Fp);
		next(d,a);
		index += 1;
	}

}



void Pi_Eval(int j,fq_mat_t f,fq_mat_t s_j,fq_t out_j){
	// f: input length-n
	//    output length-1
	//    degree-d

    // Eval f(s_j)
	Eval(f,s_j,out_j);


	//gaidong4
	// fq_t sum_gg1; // \Sum_{i=1}^n (\gamma_i-\gamma_i')
	// fq_init(sum_gg1,Fp);
	// fq_zero(sum_gg1,Fp);
	// for( int i = 0;i < n; ++ i ){
	// 	fq_add(sum_gg1,sum_gg1,fq_mat_entry(gamma,i,0),Fp);
	// 	fq_sub(sum_gg1,sum_gg1,fq_mat_entry(gamma1,i,0),Fp);
	// }

	// fq_t L0j_inv;
	// fq_init(L0j_inv,Fp);
	// fq_inv(L0j_inv,L0j,Fp);


	//gaidong5
	fq_poly_t c[n];

	fq_t zero;
	fq_init(zero,Fp);
	fq_zero(zero,Fp);

	// for(int rp=0;rp<10;++rp){
	// test_start = clock();

	fq_t mask;
	fq_init(mask,Fp);

	for (int i = 0; i < n; i ++){
		fq_poly_init(c[i],Fp);

		fq_poly_set_coeff(c[i],0,zero,Fp);
		for (int k = 1; k <= d*t; k ++){
			fq_randtest(mask, state, Fp);
			
			fq_poly_set_coeff(c[i], k, mask, Fp);
		}	
		fq_poly_evaluate_fq(mask,c[i],SID[j],Fp); 

		fq_add(out_j,out_j,mask,Fp);
	}
}

// ************ Ver algorithm **********************

void La_intpoly_coeff(fq_mat_t L0){

	// compute L_j(0) (:=fq_mat_entry(L0,j,0))
	fq_t ksubj;
	fq_init(ksubj,Fp);

	for (int j = 0; j < m; j ++){
		fq_one(fq_mat_entry(L0,j,0),Fp);

		for (int k = 1; k <= m; k++){
			if(k != (j + 1)){
				fq_mul_ui(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),k,Fp);
				fq_set_si(ksubj,k-(j+1),Fp);
				fq_div(fq_mat_entry(L0,j,0),fq_mat_entry(L0,j,0),ksubj,Fp);
			}
		}
	}
}

////interpolate F such that F(j)=Y
void IntPoly(fq_poly_t F, fq_t *Y)
{
    fq_poly_init(F,Fp);
    fq_poly_zero(F,Fp);
    
    fq_poly_t g;
    fq_poly_init(g,Fp);  
    
    fq_poly_t X;
    fq_poly_init(X,Fp);
    fq_poly_gen(X,Fp);     
    
    fq_poly_t g1,g2,h;
    fq_poly_init(g1,Fp);
    fq_poly_init(g2,Fp);
    fq_poly_init(h,Fp);
    
    fq_t a;
    fq_init(a,Fp);   

       
    for (int j=0;j<m;j++)
    {
       fq_poly_set_fq(g,Y[j],Fp);
       for (int z=0;z<m;z++)
       {
          if (z!=j) 
          {
			
             fq_poly_set_fq(h,SID[z],Fp);
             fq_poly_sub(g1,X,h,Fp);

             fq_sub(a,SID[j],SID[z],Fp);
			 
             fq_inv(a,a,Fp);
			 
             fq_poly_set_fq(g2,a,Fp);
             
             fq_poly_mul(g,g,g1,Fp);
             fq_poly_mul(g,g,g2,Fp);
          }
       }
       fq_poly_add(F,F,g,Fp);
    
    }
}

int Pi_Ver(fq_t * out, fq_t y){
	// out[1..m]
	// return 1 if y is correct; otherwise return 0
	
	int res;

    fq_poly_t F;
	fq_poly_init(F,Fp);
	IntPoly(F,out);

	int D;
	D=fq_poly_degree(F,Fp);
	fq_t zo;
	fq_init(zo,Fp);
    fq_zero(zo,Fp);

	if(D=d*t){
		res=1;
		fq_poly_evaluate_fq(y,F,zo,Fp);
	}else{
		res=0;
	}

	// printf("D=%d",D);
	return res;
}




int main(){
	int n_list[9];
	n_list[0] = 1413;
	n_list[1] = 180;
	n_list[2] = 68;
	n_list[3] = 39;
	n_list[4] = 27;
	n_list[5] = 21;
	n_list[6] = 17;
	n_list[7] = 15;
	n_list[8] = 13;

	int new_n_list[7];
	int new_d_list[7];
	new_n_list[0] = 11; new_d_list[0] = 12;
	new_n_list[1] = 10; new_d_list[1] = 13;
	new_n_list[2] = 9; new_d_list[2] = 15;
	new_n_list[3] = 8; new_d_list[3] = 17;
	new_n_list[4] = 7; new_d_list[4] = 21;
	new_n_list[5] = 6; new_d_list[5] = 27;
	new_n_list[6] = 5; new_d_list[6] = 39;

	ulong NN = 200000;
	FILE *infp = fopen("input_size.out","w");
	FILE *outfp = fopen("output_size.out","w");

  

	int repeat_times=100;

	double * share_time;
	share_time = malloc(sizeof(double)*repeat_times);

	double * eval_time; // max eval time
	eval_time = malloc(sizeof(double)*repeat_times);

	double * total_eval_time;
	total_eval_time = malloc(sizeof(double)*repeat_times);

	double * ver_time;
	ver_time = malloc(sizeof(double)*repeat_times);

	double * fx_time;
	fx_time = malloc(sizeof(double)*repeat_times);


	test_time = malloc(sizeof(double)*10);
	

    for (int index = 1; index <= 10; index ++){
    	int nn = 100 * index;
    	int dd = 2;
    	// N = NN * index;
    	// N = nn * (nn + 1) / 2;

		init(nn,dd);  // Initialize public parameters


	//gaidong
		SID=malloc(sizeof(fq_t)*m);
			for (int j=0;j<m;j++){   
   				fq_init(SID[j],Fp);
   				fq_set_ui(SID[j],j+1,Fp);   
			}

			// fq_print_pretty(SID[1],Fp);

		fmpz_t Nm;
		fmpz_init(Nm);
		fmpz_bin_uiui(Nm,n+d,d);
		N=fmpz_get_ui(Nm);

		printf("n = %d, d = %d, N = %ld\n", nn, dd, N);

	    for (int times = 0; times < repeat_times; times ++){
	    	// printf("times = %d\n", times);
			
			
			//gaidong
			fq_t temp1,temp2;
			fq_init(temp1, Fp);
			fq_set_fmpz(temp1, p, Fp);
			fq_init(temp2, Fp);
			fq_set_fmpz(temp2, p, Fp);

			// Data x
		    fq_mat_t x;
			fq_mat_init(x,n,1,Fp);
		    for (int i = 0; i < n ; i ++){
		    	fq_sub_one(temp1,temp1,Fp);
		    	fq_set(fq_mat_entry(x, i, 0), temp1, Fp);
		    }
		      
		    if (DEBUG) {
			    printf("x: ");
			    fq_mat_print_pretty(x,Fp);
			    printf("\n");
		    }
		    
		    // Shares s_1,...,s_j (in_j)
			fq_mat_t * in;
		    in=malloc(sizeof(fq_mat_t)*m);
		    for(int j = 0;j < m; j ++)
				fq_mat_init(in[j],n,1,Fp);

		    // n-variate function f of degree d
			// N=C(n+d,d) is the number of monomials in f
			

		    fq_mat_t f;
			fq_mat_init(f,N,1,Fp);
			for(int i=0;i<N;++i){
				fq_sub_one(temp2,temp2,Fp);
				fq_set(fq_mat_entry(f, i, 0), temp2, Fp);
			}

		    // Output of servers (out_j)
		    fq_t * out=malloc(sizeof(fq_t)*m);
		    for(int j=0;j<m;++j){
		        fq_init(out[j],Fp);
				fq_zero(out[j],Fp);
		    }

			

		//**********************************
			
		    // Share

		    // printf("***Sharing***\n");
			clock_t share_start,share_end;

			share_start = clock();
			Pi_Share(x,in);
			share_end = clock();
			share_time[times] = (double) (share_end- share_start)/CLOCKS_PER_SEC;
		    // printf("Time of a Share: %f ms\n", share_time[times]* 1000);


			//gaidong1
			//preprocess Lagrange interpolation coefficients
			// fq_mat_t L0;
			// fq_mat_init(L0,m,1,Fp);
			// fq_mat_entry(L0,j,0): L_j(0)= \prod_{k=1,k!=j}^m k/(k-j)
			// La_intpoly_coeff(L0); 

			// Eval
			

			clock_t eval_start,eval_end;
			

			// gaidong2
			// fq_mat_t * gamma;
			// fq_mat_t * gamma1;
			// gamma=malloc(sizeof(fq_mat_t)*m);
			// gamma1=malloc(sizeof(fq_mat_t)*m);
			// for (int j = 0;j < m; j ++){
			// 	fq_mat_init(gamma[j],n,1,Fp);
			// 	fq_mat_init(gamma1[j],n,1,Fp);

			// 	fq_mat_randtest(gamma[j],state,Fp);
			// }

			// for (int j = 0;j < m - 1; j ++){
			// 	fq_mat_set(gamma1[j],gamma[j+1],Fp);
			// }
			// fq_mat_set(gamma1[m-1],gamma[0],Fp);

			

		    total_eval_time[times] = 0;
			for (int j = 0;j < m; j ++){
			    eval_start=clock();
		        Pi_Eval(j,f,in[j],out[j]);  //gaidong3
			    eval_end=clock();
			    total_eval_time[times]  += (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    if ( (double) (eval_end-eval_start)/CLOCKS_PER_SEC > eval_time[times] ) {
			    	eval_time[times] = (double) (eval_end-eval_start)/CLOCKS_PER_SEC;
			    }
		    }
		    // printf("Time of Max_Eval: %f ms\n",eval_time[times]*1000);
		    // printf("Time of total Eval: %f ms\n",total_eval_time[times] *1000);

		    // Ver
		    // printf("***Verifying***\n");
	
			

			// if (DEBUG) {
			// 	printf("L_1(0)..L_m(0):");
			// 	fq_mat_print_pretty(L0,Fp);
			// }  
			

		    int IsCorrect;
		    fq_t y;
		    fq_init(y,Fp);
		    clock_t Ver_start,Ver_end;
		    Ver_start=clock();
		    IsCorrect = Pi_Ver(out,y);
		    Ver_end=clock();

			if (DEBUG) {
				printf("Ac? %d, y=",IsCorrect);
				fq_print_pretty(y,Fp);
				printf("\n");
			}  
		    
		    ver_time[times] = (double) (Ver_end-Ver_start)/CLOCKS_PER_SEC;
		    // printf("Time of Ver: %f ms\n",ver_time[times] * 1000);

			



			fq_t fx;
		    fq_init(fx,Fp);
		    // Eval f(x)
		    clock_t fx_start,fx_end;
		    fx_start=clock();
		    Eval(f,x,fx);
		    fx_end=clock();

			if (DEBUG) {
				printf("f(x)=");
				fq_print_pretty(fx,Fp);
				printf("\n");
			}

		    fx_time[times]= (double) (fx_end-fx_start)/CLOCKS_PER_SEC;
		    // printf("Time of directly eval f(x): %f ms\n",fx_time[times]*1000);

		    // communication cost

			fq_mat_t PRF_KEY;
			fq_mat_init(PRF_KEY,n,1,Fp);
			fq_mat_randtest(PRF_KEY,state,Fp);
		    for (int j = 0; j < m; j ++){
		    	fq_mat_fprint(infp, in[j], Fp);
				fq_mat_fprint(infp, PRF_KEY, Fp);
		    	fq_fprint(outfp, out[j], Fp);
		    }

			// clear memory
		    fq_mat_clear(f,Fp);
			fq_mat_clear(x,Fp);
			for(int j=0;j<m;++j){
		        fq_clear(out[j],Fp);
		    }
			// fq_mat_clear(L0,Fp);
			fq_clear(y,Fp);
		    fq_clear(fx,Fp);
		    
		    for(int j = 0;j < m; j ++)
				fq_mat_clear(in[j],Fp);
		}


	
		double AVG_share_time;


	double AVG_eval_time=0; // max eval time

	double AVG_total_eval_time=0;

	double AVG_ver_time=0;

	double AVG_fx_time=0;
	for(int i=0;i<repeat_times;++i){
		AVG_share_time+=share_time[i];
		// printf("Time of a Share: %f ms\n", share_time[i]* 1000);
		AVG_eval_time+=eval_time[i];
		AVG_total_eval_time+=total_eval_time[i];
		AVG_ver_time+=ver_time[i];
		AVG_fx_time+=fx_time[i];
	}
	AVG_share_time=AVG_share_time/repeat_times;
	AVG_eval_time=AVG_eval_time/repeat_times;
	AVG_total_eval_time=AVG_total_eval_time/repeat_times;
	AVG_ver_time=AVG_ver_time/repeat_times;
	AVG_fx_time=AVG_fx_time/repeat_times;
	

	double SDAVG_share_time=0;

	double SDAVG_eval_time=0; // max eval time

	double SDAVG_total_eval_time=0;

	double SDAVG_ver_time=0;

	double SDAVG_fx_time=0;

	double dis;
	
	for(int i=0;i<repeat_times;++i){
		dis=share_time[i]-AVG_share_time;
		SDAVG_share_time+=dis*dis;

		dis=eval_time[i]-AVG_eval_time;
		SDAVG_eval_time+=dis*dis;

		dis=total_eval_time[i]-AVG_total_eval_time;
		SDAVG_total_eval_time+=dis*dis;

		dis=ver_time[i]-AVG_ver_time;
		SDAVG_ver_time+=dis*dis;

		dis=fx_time[i]-AVG_fx_time;
		SDAVG_fx_time+=dis*dis;
	}

	SDAVG_share_time=SDAVG_share_time/repeat_times;
	SDAVG_eval_time=SDAVG_eval_time/repeat_times;
	SDAVG_total_eval_time=SDAVG_total_eval_time/repeat_times;
	SDAVG_ver_time=SDAVG_ver_time/repeat_times;
	SDAVG_fx_time=SDAVG_fx_time/repeat_times;

	SDAVG_share_time=sqrt(SDAVG_share_time)/AVG_share_time;
	SDAVG_eval_time=sqrt(SDAVG_eval_time)/AVG_eval_time;
	SDAVG_total_eval_time=sqrt(SDAVG_total_eval_time)/AVG_total_eval_time;
	SDAVG_ver_time=sqrt(SDAVG_ver_time)/AVG_ver_time;
	SDAVG_fx_time=sqrt(SDAVG_fx_time)/AVG_fx_time;

	printf("AVG Time of Share: %f ms\n",AVG_share_time *1000);
	
	printf("AVG Max Time of Eval: %f ms\n",AVG_eval_time *1000);
	
	printf("AVG Total Time of Eval: %f ms\n",AVG_total_eval_time *1000);
	
	printf("AVG Time of Ver: %f ms\n",AVG_ver_time *1000);
	
	printf("AVG Time of f(x): %f ms\n",AVG_fx_time *1000);

	printf("SD/AVG of Share time: %f \n",SDAVG_share_time);
printf("SD/AVG of Max Eval time: %f \n",SDAVG_eval_time);
printf("SD/AVG of Eval time: %f \n",SDAVG_total_eval_time);
printf("SD/AVG of Ver time: %f \n",SDAVG_ver_time);
	printf("SD/AVG of f(x): %f \n",SDAVG_fx_time);


		fprintf(infp,"\n\n\n");
		fprintf(outfp,"\n\n\n");
	}

	

	fclose(infp);
	fclose(outfp);

    return 0;
}
