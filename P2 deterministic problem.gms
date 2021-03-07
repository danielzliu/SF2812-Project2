
Sets
  I   units / unit1*unit3 /
  J   period / period1*period5 /;
  
  
Parameters
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
       
       
Variables
   X(I,J) Amount power to be produced by unit I at time J. Should be either 0 or  lower < x < upper
   Y(I,J) Binary variable that is 1 if unit I was switched on at time J from warm start
   W(I,J) -||- from cold start
   T(I,J) Binary variable that is 1 if unit I is active at time J
   Z      Variabe to be minimized;
   
Positive variables
   X;
   
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
   TotalCost    Total cost;

Satdem(J)        .. SUM(I,X(I,J)) =G= D(J);
BIN1(I,J)        .. X(I,J) =G= T(I,J)*L(I);
BIN2(I,J)        .. X(I,J) =L= T(I,J)*U(I);
ThreePeriods(I)  .. SUM(J, T(I,J)) =L= 3;
ColdStart(I,J)     .. T(I,J) - (T(I,J--1) + T(I,J--2)) =L= W(I,J);
WarmStart(I,J)     .. T(I,J) - T(I,J--1) =L= Y(I,J) + W(I,J);

TotalCost        .. Z =E= SUM((I,J), X(I,J)*C(I) + Y(I,J)*K(I) + W(I,J)*E(I));


MODEL VARME /ALL/;

SOLVE VARME USING MIP MINIMIZING Z;

DISPLAY X.L, Y.L, Z.L, W.L;