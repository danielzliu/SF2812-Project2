$funclibin stolib stodclib
Functions randnorm     /stolib.dnormal    /;

*https://www.math.kth.se/optsyst/grundutbildning/5B1815/farm.gms


Sets
  I   units / unit1*unit3 /
  J   period / period1*period5 /
  Stoch stochastic outcome /s1*s1000/;
  
  
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
       unit3  55 /;
       
    Parameter
    StochDev(J,Stoch) The deviation of demand;
    StochDev(J,Stoch) =randnorm(0,5);

    Parameter StochDevSc(J);

    Parameter
    StochProb(Stoch) probability of each outcome;
    StochProb(Stoch) = 1/card(Stoch)
       
       
Variables
   X(I,J) Amount power to be produced by unit I at time J. Should be either 0 or  lower < x < upper
   Exp(J) Amount of power to be bought from other suppliers
    Y(I,J) Binary variable that is 1 if unit I was switched on at time J from warm start
   W(I,J) -||- from cold start
   T(I,J) Binary variable that is 1 if unit I is active at time J
   CostScen      Variabe to be minimized;
   
Positive variables
   X
   Exp;
   
Binary variables
   Y
   W
   T;
   
Equations
   Satdem(J)  Making sure that the units satisfies the demand at the period.
   BIN1(I,J)  Making sure X(IJ) is between its boundries or 0
   BIN2(I,J)  -||-
   ThreePeriods(I) Making sure one unit is not turned on more than three periods in a row
   ColdStart(I,J) Binary for cold start
   WarmStart(I,J) Binary for warm start
   CostScenDef    Total cost;

Satdem(J)        .. SUM(I,X(I,J)) + Exp(J) =G= D(J) + StochDevSc(J);
BIN1(I,J)        .. X(I,J) =G= T(I,J)*L(I);
BIN2(I,J)        .. X(I,J) =L= T(I,J)*U(I);
ThreePeriods(I)  .. SUM(J, T(I,J)) =L= 3;
ColdStart(I,J)     .. T(I,J) - (T(I,J--1) + T(I,J--2)) =L= W(I,J);
WarmStart(I,J)     .. T(I,J) - T(I,J--1) =L= Y(I,J) + W(I,J);
CostScenDef .. CostScen =E= SUM(J, (SUM(I, P(J)*C(I)*X(I,J) + Y(I,J)*K(I) + W(I,J)*E(I)) + 10*Exp(J)))

MODEL VARME /ALL/;

scalar totalCost / 0 /;

LOOP(Stoch,

StochDevSc(J) = StochDev(J, Stoch)

SOLVE VARME USING MIP MINIMIZING CostScen;


TotalCost = TotalCost + StochProb(Stoch) * CostScen.L;

)






DISPLAY TotalCost;