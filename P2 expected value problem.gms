$funclibin stolib stodclib
Functions randnorm     /stolib.dnormal    /;

Sets
  I   units / unit1*unit3 /
  J   period / period1*period5 /
  Stoch stochastic outcome /s1*s50/;
  
  
Parameters
    
    P(J) Number of hours in period J
    / period1  5
       period2  5
       period3  5
       period4  5
       period5  4 /

  C(I) cost per hour and MW
     / unit1  2.5
       unit2  2.5
       unit3  2.4 /
  
  K(I) initial cost for warm start
     / unit1  10
       unit2  13
       unit3  16 /
  
  E(I) Initial cost for a warm start
     / unit1  15
       unit2  19.5
       unit3  24  /
       
  D(J) demand at period j
     / period1  50
       period2  60
       period3  80
       period4  70
       period5  60 /
       
   L(I) Lower boundry for unit I
     / unit1  10
       unit2  12
       unit3  15 /
       
   U(I) Upper boundry for unit I
     / unit1  50
       unit2  45
       unit3  55 /
       
    T(I,J)
     / unit1.period1 1
     unit1.period2 1
     unit1.period3 0
     unit1.period4 0
     unit1.period5 1
     unit2.period1 0
     unit2.period2 0
     unit2.period3 1
     unit2.period4 1
     unit2.period5 1
     unit3.period1 0
     unit3.period2 1
     unit3.period3 1
     unit3.period4 1
     unit3.period5 0
        /
        
    W(I,J)
     / unit1.period1 0
     unit1.period2 0
     unit1.period3 0
     unit1.period4 0
     unit1.period5 1
     unit2.period1 0
     unit2.period2 0
     unit2.period3 1
     unit2.period4 0
     unit2.period5 0
     unit3.period1 0
     unit3.period2 1
     unit3.period3 0
     unit3.period4 0
     unit3.period5 0
        /
    
 Y(I,J);
     Y(I,J) = 0;
       
    Parameter
    StochDev(J,Stoch) The deviation of demand;
    StochDev(J,Stoch) =randnorm(0,5);
    
    Parameter
    StochProb(Stoch) probability of each outcome;
    StochProb(Stoch) = 1/card(Stoch)
       
       
Variables
   X(I,J,Stoch) Amount power to be produced by unit I at time J. Should be either 0 or  lower < x < upper
   Exp(J,Stoch) Amount of power to be bought from other suppliers

   Z      Variabe to be minimized;
   
Positive variables
   X
   Exp;
   
   
Equations
   Satdem(J,Stoch)  Making sure that the units satisfies the demand at the period.
   BIN1(I,J,Stoch)  Making sure X(IJ) is between its boundries or 0
   BIN2(I,J,Stoch)  -||-
   ThreePeriods(I) Making sure one unit is not turned on more than three periods in a row
   ColdStart(I,J) Binary for cold start
   WarmStart(I,J) Binary for warm start
   TotalCost    Total cost;

Satdem(J,Stoch)        .. SUM(I,X(I,J,Stoch)) + Exp(J,Stoch) =G= D(J) + StochDev(J,Stoch);
BIN1(I,J, Stoch)        .. X(I,J,Stoch) =G= T(I,J)*L(I);
BIN2(I,J, Stoch)        .. X(I,J,Stoch) =L= T(I,J)*U(I);
ThreePeriods(I)  .. SUM(J, T(I,J)) =L= 3;
ColdStart(I,J)     .. T(I,J) - (T(I,J--1) + T(I,J--2)) =L= W(I,J);
WarmStart(I,J)     .. T(I,J) - T(I,J--1) =L= Y(I,J) + W(I,J);

TotalCost        .. Z =E= SUM((J,Stoch), StochProb(Stoch)* P(J)*(SUM(I,C(I)* X(I,J,Stoch)) + 10*Exp(J,Stoch))) + SUM((I,J), Y(I,J)*K(I) + W(I,J)*E(I));


MODEL VARME /ALL/;

SOLVE VARME USING MIP MINIMIZING Z;

DISPLAY X.L, Z.L;