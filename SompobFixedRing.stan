// ##################################################################################################
// # Saralamba et al. Intrahost modeling of artemisinin resistance
// #   in Plasmodium falciparum PNAS  397-402, 
// #   doi: 10.1073/pnas.1006113108
// # 
// # R/Stan Version adapted from 
// # http://demonstrations.wolfram.com/AModelOfPlasmodiumFalciparumPopulationDynamicsInAPatientDuri/
// # by sompob@tropmedres.ac
// #
//#################################################################################################
functions{
  
  real PDF_normal(real mu, real sig, real x){
    //return exp(normal_lpdf(x| mu, sig));
    return (1/(exp( ((x-mu)*(x-mu))/(2*sig*sig))*sqrt(2*pi())*sig));    
  }
  
  vector InitDist(real N0, real mu, real sigma, real PMR){
    real sumP;
    vector[48] p;
    vector[48] ls;
    
    p=rep_vector(0,48);
    //based on White et al. model
    for(age in 1:48)
      p[age] = (PDF_normal(mu,sigma,age-48)/PMR)+PDF_normal(mu,sigma,age)+(PDF_normal(mu,sigma,age+48)*PMR);
    
    sumP = sum(p);
    // if(sumP==0){
    //   sumP=0.0000000001;
    // }
    //print(sumP);
    for(age in 1:48)
      ls[age] = pow(10,N0)*p[age]/sumP;
    
    return ls;
  }
  
  vector RotateRight(vector ls, real PMR){
    //size of the input vector is fixed at 48
    vector[48] tmpvec;
    tmpvec[1] = ls[48]*PMR;
    for(age in 2:48)
      tmpvec[age]=ls[age-1];
      
    return tmpvec;
  }
  
  //#drug efficacy (killing rate)
  real Eff(real conc, real gamma, real ec50, real emin, real emax){
    return emin+(emax-emin)*pow(conc,gamma)/(pow(conc,gamma) + pow(ec50,gamma));
  }
  
  // #sequestration function
  // #probability of the parasites at age to be Ring stage
  real PRingfunc(int age, int a1,int a2){
    //real prob;
    // if(age < a1){
    //   prob = 1.0;
    // }else{
    //   prob = exp(log(0.5)*(age-a1)/(a2-a1));
    // }
    //return prob;
    return ((age < a1)?1.0:exp(log(0.5)*(age-a1)/(a2-a1)));
  }

  // #the function for counting the circulate parasites (rings)
  // #based on Saralamba et al.'s model 
  real CountRing(vector ls,int a1, int a2){
    vector[48] pls;
    for(i in 1:48)
      pls[i]=PRingfunc(i,a1,a2);
    
    return sum(pls .* ls);
  }  
  
  //#calculate the concentration @ time t
  real Concfn(real xm, real ym, real ke, real t){
    return ((t<=xm)?(ym-(ym/xm)*(xm-t)):ym*exp(-ke*(t-xm)));

  }

  matrix DoseResponse(real xm, real ym, real ke, real gamma, real ec50, real emin, real emax, int everyH, int Ndrug){
    
    matrix[500,2] ec;
    //int	i;    //time from 1 to drugperiod
    int j;    //for counting the time interval between each dose
    int nd;   // number of doses
    int drugperiod;
    real conc;
    
    
    drugperiod = everyH*Ndrug;  
    ec=rep_matrix(0,500,2);   //maximum time that the model can run is 500 hrs
    
    conc=0.0;
    j = 0;  
		nd = 1;  
		for(i in 1:drugperiod){
		  if(j < everyH){
		    conc = Concfn(xm,ym,ke,j);
		    //print(conc);
		    ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emin,emax);
			}
			if((j==everyH)&&(nd<Ndrug)){
			  nd=nd+1;
			  conc=conc+Concfn(xm,ym,ke,j);
			  //print(conc);
			  ec[i,1]=conc;
			  ec[i,2]=Eff(ec[i,1],gamma,ec50,emin,emax);
			  j=0;
			}
			j=j+1;
		}
		return ec;    
  }

  vector Ki(real T, matrix ecls){
    vector[rows(ecls)] Kils;
    real alpha;
    alpha=1/T;
    for(i in 1:rows(ecls))
      Kils[i]=alpha*log(100.0/(100.0-ecls[i,2]));
    
    return Kils;
  }

  matrix BigKi(vector rKi, vector tKi, vector sKi){
    return(append_col(append_col(rKi,tKi),sKi));
  }


	vector WhichRTS(int rb,int re, int tb,int te, int sb, int se){
	  //output vector size is fixed at 48.
	  vector[48] tmp;
	  tmp=rep_vector(0,48);
	  for(i in 1:48){
	    if(rb<=i && i<=re)
	      tmp[i]=1;
	      
	    if(tb<=i && i<=te)
	      tmp[i]=2;
	      
	    if(sb<=i && i<=se)
	      tmp[i]=3;
	  }
	  return tmp;
	}
  

  vector Fdecay(int attime, vector stages, vector lst, matrix bigKi){
    vector[48] tmp;
    tmp = lst; 
    for(age in 1:48){
      if(stages[age]==0)
        tmp[age] = lst[age]*1.0; //no drug effect on this age --> do nothing
      
      if(stages[age]==1)
        tmp[age] = lst[age]*exp(-bigKi[attime,1]);
      
      if(stages[age]==2)
        tmp[age] = lst[age]*exp(-bigKi[attime,2]);
      
      if(stages[age]==3)
        tmp[age] = lst[age]*exp(-bigKi[attime,3]);
    }
    return tmp;    
  }
  
  //observed time for Pailin data
  //based on the data used in 
  int[] ObservedTimeP(int np){
    int tmp[31];  //obseved time in hours
    int time[np];
    
    	tmp[1]=0;
			tmp[2]=2;
			tmp[3]=4;
			tmp[4]=6;
			tmp[5]=8;
			tmp[6]=12;
			tmp[7]=18;
			tmp[8]=24;
			tmp[9]=30;
			tmp[10]=36;
			tmp[11]=42;
			tmp[12]=48;
			tmp[13]=54;
			tmp[14]=60;
			tmp[15]=66;
			tmp[16]=72;
			tmp[17]=78;
			tmp[18]=84;
			tmp[19]=90;
			tmp[20]=96;
			tmp[21]=102;
			tmp[22]=108;
			tmp[23]=114;
			tmp[24]=120;
			tmp[25]=126;
			tmp[26]=132;
			tmp[27]=138;
			tmp[28]=144;
			tmp[29]=150;
			tmp[30]=156;
			tmp[31]=162;
			
			for(i in 1:np){
			  time[i]=tmp[i];
			}
      return time;
  }
  

  vector SaralambaModel(real N0,real mu,real sigma,real PMR,int KZRb,int KZRe,int KZTb, int KZTe, int KZSb, int KZSe,
    real xm,real ym,real ke,int everyH,int Ndrug,real gammaR,real gammaT,real gammaS,real ec50R,real ec50T,real ec50S,
    real eminR,real eminT,real eminS,real emaxR,real emaxT,real emaxS,real bigT,int npoints){
      
      int timewant[npoints];
      int runmax;
      //500 = the maximum run times of the model
      matrix[500,2] eclsR;
      matrix[500,2] eclsT;
      matrix[500,2] eclsS;
      vector[500] KiR;
      vector[500] KiT;
      vector[500] KiS;
      matrix[500,3] bigKi;
      vector[48] stages;
      vector[48] ls;
      int tNpoints;
      int tOb;
      vector[npoints] outputvector;
      
      timewant=ObservedTimeP(npoints);
      runmax=timewant[npoints];
      stages=WhichRTS(KZRb,KZRe,KZTb,KZTe,KZSb,KZSe);
      
			eclsR = DoseResponse(xm,ym,ke,gammaR,ec50R,eminR,emaxR,everyH,Ndrug);
			eclsT = DoseResponse(xm,ym,ke,gammaT,ec50T,eminT,emaxT,everyH,Ndrug);
			eclsS = DoseResponse(xm,ym,ke,gammaS,ec50S,eminS,emaxS,everyH,Ndrug);
      KiR = Ki(bigT,eclsR);
			KiT = Ki(bigT,eclsT);
			KiS = Ki(bigT,eclsS);
			bigKi = BigKi(KiR,KiT,KiS);
    
      ls=InitDist(N0,mu,sigma,PMR);
      
      tNpoints=1;
      outputvector[tNpoints]=log10(CountRing(ls,11,14));  //count at t=0 on admission
      tNpoints=tNpoints+1;
      
      for(t in 1:runmax){
        tOb=timewant[tNpoints];
        ls=RotateRight(ls,PMR);
        ls=Fdecay(t,stages,ls,bigKi);
        
        if(t == tOb && tNpoints <= npoints){
          outputvector[tNpoints]=log10(CountRing(ls,11,14));
          tNpoints=tNpoints+1;
        }
      }
      
      return outputvector;
    }

  vector SaralambaModel4CT(real N0,real mu,real sigma,real PMR,int KZRb,int KZRe,int KZTb, int KZTe, int KZSb, int KZSe,
    real xm,real ym,real ke,int everyH,int Ndrug,real gammaR,real gammaT,real gammaS,real ec50R,real ec50T,real ec50S,
    real eminR,real eminT,real eminS,real emaxR,real emaxT,real emaxS,real bigT,int runmax){
      
      //500 = the maximum run times of the model
      matrix[500,2] eclsR;
      matrix[500,2] eclsT;
      matrix[500,2] eclsS;
      vector[500] KiR;
      vector[500] KiT;
      vector[500] KiS;
      matrix[500,3] bigKi;
      vector[48] stages;
      vector[48] ls;

      vector[(runmax+1)] outputvector;   //include t=0
      
      stages=WhichRTS(KZRb,KZRe,KZTb,KZTe,KZSb,KZSe);
      
			eclsR = DoseResponse(xm,ym,ke,gammaR,ec50R,eminR,emaxR,everyH,Ndrug);
			eclsT = DoseResponse(xm,ym,ke,gammaT,ec50T,eminT,emaxT,everyH,Ndrug);
			eclsS = DoseResponse(xm,ym,ke,gammaS,ec50S,eminS,emaxS,everyH,Ndrug);
      KiR = Ki(bigT,eclsR);
			KiT = Ki(bigT,eclsT);
			KiS = Ki(bigT,eclsS);
			bigKi = BigKi(KiR,KiT,KiS);
    
      ls=InitDist(N0,mu,sigma,PMR);
      
      outputvector[1]=log10(CountRing(ls,11,14));  //count on admission

      for(t in 1:runmax){
        ls=RotateRight(ls,PMR);
        ls=Fdecay(t,stages,ls,bigKi);
        outputvector[t+1]=log10(CountRing(ls,11,14));
      }
      
      return outputvector;
    }

	//find the clearance time from the parasite count vector generated by SaralambaModel4CT
    real ClearanceTime(vector countvector, real deteclim){
      real ct;
      int veclen;
      int indx;
      int run;
      run = 1;  //run the loop 1; stop the loop 0
      veclen = num_elements(countvector);
      indx=1;
      ct = 240;
      
      while(indx <= veclen && run==1){
        if(countvector[indx] <= deteclim){
          ct = indx;
          run = 0;
        }
        indx=indx+1;
      }
      return ct;
    }
  
  
} //end function block

data{
  int NPPARA; //number of data points
  real y[NPPARA];  
  
  real xm;  //time at max conc
  real ym;  //maximum conc
  real ke;   //elimination rate of the drug
  real DetecLim;  //detection limit
}
parameters{
  real<lower=0.01> SD;
  real<lower=6,upper=14> N0;
  real<lower=1> mu;
  real<lower=1> sig;
  real<lower=1> pmr;

  // real<lower=1.5> gammaR;
  // real<lower=1.5> gammaT;
  // real<lower=1.5> gammaS;
  real<lower=0.01> ec50R;
  real<lower=0.01> ec50T;
  real<lower=0.01> ec50S;
  real<lower=50,upper=100> emaxR;
  real<lower=50,upper=100> emaxT;
  real<lower=50,upper=100> emaxS;
  //real<lower=1,upper=12> bigT;
  //real<lower=0> nu;
}

model
{
  //matrix[NPPARA,2] mod_pred;
  vector[NPPARA] mod_pred;

	SD ~ uniform(0.01,3);		

	N0 ~ uniform(8,13);
	mu ~ uniform(1,48);
	sig ~ uniform(1,48);
	pmr ~ uniform(1,30);
	// xm ~ uniform(0,14);
	// ym ~ uniform(100,3000);
	// ke ~ uniform(0,2);
	//gammaR ~ uniform(1.5,9.5);
	//gammaT ~ uniform(1.5,9.5);
	//gammaS ~ uniform(1.5,9.5);
	ec50R  ~ uniform(5,100);
	ec50T  ~ uniform(5,100);
	ec50S  ~ uniform(5,100);
	emaxR	~ uniform(50,99.99);
	emaxT	~ uniform(50,99.99);
	emaxS	~ uniform(50,99.99);
	//bigT ~ uniform(1,12);
  
  //KZs are fixed and T is 1; gamma_RTS are fixed at 6.5
  //AS7 
  // mod_pred=SaralambaModel(N0,mu,sig,pmr,6,26,27,38,39,44,xm,ym,ke,24,7,
  //   gammaR,gammaT,gammaS,ec50R,ec50T,ec50S,0,0,0,emaxR,emaxT,emaxS,bigT,NPPARA);
  mod_pred=SaralambaModel(N0,mu,sig,pmr,6,26,27,38,39,44,xm,ym,ke,24,7,
    6.5,6.5,6.5,ec50R,ec50T,ec50S,0,0,0,emaxR,emaxT,emaxS,1,NPPARA);
  
  y~normal(mod_pred,SD);
}
generated quantities{
  vector[NPPARA] y_pred;
  vector[NPPARA] log_lik;
  vector[NPPARA] mod_pred;
  vector[241] mod_predCT;			//240+1 for t=0
  vector[241] mod_predSplitCT;
  real riR;
  real riT;
  real riS;
  real AS7CT;
  real SplitCT;
  
  riR=ec50R/(6.5*emaxR);
  riT=ec50T/(6.5*emaxT);
  riS=ec50S/(6.5*emaxS);
  

  mod_pred=SaralambaModel(N0,mu,sig,pmr,6,26,27,38,39,44,xm,ym,ke,24,7,
    6.5,6.5,6.5,ec50R,ec50T,ec50S,0,0,0,emaxR,emaxT,emaxS,1,NPPARA);
  mod_predCT=SaralambaModel4CT(N0,mu,sig,pmr,6,26,27,38,39,44,xm,ym,ke,24,7,
    6.5,6.5,6.5,ec50R,ec50T,ec50S,0,0,0,emaxR,emaxT,emaxS,1,240);
  mod_predSplitCT=SaralambaModel4CT(N0,mu,sig,pmr,6,26,27,38,39,44,xm,ym,ke,12,14,
    6.5,6.5,6.5,ec50R,ec50T,ec50S,0,0,0,emaxR,emaxT,emaxS,1,240);
  
  AS7CT=ClearanceTime(mod_predCT,DetecLim);
  SplitCT=ClearanceTime(mod_predSplitCT,DetecLim);
  
  for (n in 1:NPPARA) {
    y_pred[n] = normal_rng(mod_pred[n], SD);
    log_lik[n]=normal_lpdf(y[n]| mod_pred[n],SD);
  }
  
}
