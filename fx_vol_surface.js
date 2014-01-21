// The following javascript code was written by Philip Kinlen.

// Discount Factor
function sk_discFact(yearsToMat, rate)
{
    return Math.pow( 1 + rate, -yearsToMat)
}
      
// for equity forwards, acc will be the interest rate in the accounting currency
// and undl will be the dividend yield
// for FX forwards, acc is the interest rate in the accounting currency
// undl is the interest rate in the underlying currency      
function sk_forward(spot, yearsToMat, accountingCcyYield, underlyingYield)
{
    return spot* sk_discFact(yearsToMat, underlyingYield) / sk_discFact(yearsToMat, accountingCcyYield);
} 

// The cumulative normal of a standard normal distribution: N(0,1)
function sk_cumNorm(x)
{
  if (x > 6.0)
    return 1.0;
  if (x < -6.0)
    return 0.0;

  var b1 = 0.31938153;
  var b2 = -0.356563782;
  var b3 = 1.781477937;
  var b4 = -1.821255978;
  var b5 = 1.330274429;
  var p  = 0.2316419;
  var c2 = 0.3989423;

  var a = Math.abs(x);
  var t = 1.0 / (1.0 + a * p);
  var b = c2* Math.exp((-x)*(x/2.0));
  var n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  n = 1.0-b*n;
  
  if ( x < 0.0 )
    n = 1.0 - n;
  
  return n;
}



function sk_checkBSInputs(PutOrCall, forwardVal, strike, discountFactor, volatility, yearsToMat)
{
     if((PutOrCall != "put") && (PutOrCall != "call"))
       throw "PutOrCall paramter just be either 'put' or 'call' here it is: " + String(PutOrCall);
     
     if(forwardVal < 0)
       throw "forward value must be greater than or equal to zero, here it is: " + String(forwardVal);
  
     if(strike < 0)
       throw "strike must be greater than or equal to zero, here it is: " + String(strike);
  
     if(discountFactor < 0)
       throw "discount factor must be greater than or equal to zero, here it is: " + String(discountFactor);
   
     if(volatility < 0)
       throw "vol must be greater than or equal to zero, here it is: " + String(volatility);

     if(yearsToMat < 0)
       throw "Option has expired, i.e. yearsToMat is negative: " + String(yearsToMat);
}
   

// black scholes calculator without checks of the inputs.
// Most users should call the sk_blackScholesValue(..) function since it will check the inputs.
function sk_bsVal(PutOrCall, forwardVal, strike, discountFactor, volatility, yearsToMat)
{
     if(forwardVal == 0)
     {
         if(PutOrCall == "put")
              return strike * discountFactor;
         else // is call
              return 0;            
     }

     if(strike == 0)
     {
         if(PutOrCall == "put")
              return 0;
         else // is call
              return forwardVal * discountFactor;              
     }
        
     var volSqrtT = volatility * Math.sqrt(yearsToMat);

     if(volSqrtT== 0)
     {
         if(PutOrCall == "put")
              return Math.max(strike - forwardVal,0) * discountFactor;
         else // is call
              return Math.max(forwardVal - strike,0) * discountFactor;
     }
                
     var d1 = Math.log(forwardVal / strike)/volSqrtT + volSqrtT * 0.5;                
     var d2 = d1 - volSqrtT;
   
     if(PutOrCall == "put")
       return discountFactor * ( strike * sk_cumNorm(-d2) - forwardVal * sk_cumNorm(-d1));
     else // is Call
       return discountFactor * ( forwardVal * sk_cumNorm(d1) - strike * sk_cumNorm(d2));
}

// the main Black Scholes Calculator function, with some checks of the input
function sk_blackScholesValue( PutOrCall, forwardVal, strike, discountFactor, volatility, yearsToMat)
{
     // the following function call may throw on error
     sk_checkBSInputs( PutOrCall, forwardVal, strike, discountFactor, volatility, yearsToMat);
     return sk_bsVal(PutOrCall, forwardVal, strike, discountFactor, volatility, yearsToMat);
}


function sk_blackScholesDelta( PutOrCall, forwardVal, strike, underlyingYield, volatility, yearsToMat)
{
     var undlDiscFactor = sk_discFact(yearsToMat, underlyingYield);
  
     // the following function call may throw an error
     sk_checkBSInputs( PutOrCall, forwardVal, strike, undlDiscFactor, volatility, yearsToMat);  
  
     var volSqrtT = volatility * Math.sqrt(yearsToMat);

     if(volSqrtT == 0)
     {
         if(   ( (PutOrCall == "put")  && (forward < strike))
            || ( (PutOrCall == "call") && (forward > strike)))
            return undlDiscFactor; // for an in the money option
         else // for an out of the money option
            return 0; 
     }
     else if(strike == 0)
     {
         if(PutOrCall == "put")
            return 0;
         else // is call 
            return undlDiscFactor;
     }
     else if(forwardVal == 0)
     {
         throw "When working out the delta found the forward to be 0."; 
     }
     else
     {                
        var d1 = Math.log(forwardVal / strike)/volSqrtT + volSqrtT * 0.5;               
  
        if(PutOrCall == "put")
           return - sk_cumNorm(-d1) * undlDiscFactor;
        else // is call
           return  sk_cumNorm(d1) * undlDiscFactor;
      }
 }

// We have a piecewise linear interpolator that does flat extrapolation before the first x and after the last x.   
// Vector arrayX must be sorted and strictly increasing   
function sk_linearInterpArray(arrayX, arrayY, x)
{
   if( (arrayX.length * 1) != (arrayY.length * 1))
      throw "Found arrayX to have length: " + String(arrayX.length) + " and arrayY to have length: " 
            + String(arrayY.length) + ". They need to have equal length.";
   
   var indexBelow = -1;
   var indexAbove = -1;
   for( var i  = 0; i < arrayX.length ; i++)
   {
       if (arrayX[i]*1 <= x) // multiply by 1 to force convertion to a Number.
           indexBelow = i;

       if ((arrayX[i]*1 > x) && (indexAbove == -1))
            indexAbove = i;

       if((i > 0) && ( arrayX[i]*1 <= arrayX[i-1])) 
           throw "For i = " + String(i) + " found arrayX[i] = " + String(arrayX[i]) 
                 + ", which is not strictly greater than arrayX[i-1] = " + String(arrayX[i-1]);
   }
     
   if(indexAbove == -1) // i.e. all arrayX[i]'s are below x 
     return arrayY[indexBelow];

   if(indexBelow == -1) // i.e. all arrayX[i]'s are above x
     return arrayY[indexAbove];

   return sk_linearInterp(arrayX[indexBelow], arrayX[indexAbove], arrayY[indexBelow], arrayY[indexAbove], x);
}             

// linear interpolation with extrapolation
function sk_linearInterp( x1, x2, y1, y2, x)
{
    if( x1 == x2)
        throw "x1 (" + String(x1) + ") and x2 (" + String(x2) + ") must NOT no equal.";  

    return ( y1 * ( x2 - x) + y2 * ( x - x1)) / ( x2 - x1);
}
  
function sk_getVolFromSurface( deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                               spot, yearsToMat, accountingCcyYield, underlyingYield, strike)
{
    return sk_getVolOrDeltaFromVolSurface( deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                                           spot, yearsToMat, accountingCcyYield, underlyingYield, strike, "vol");
}  

// The following function will work out the a, b, c such that:
// vol( delta )    = a + b * (delta) + c * (delta * delta)
// Note: we're using the delta of a call, which goes from 1 ( for a zero strike ) to 0 for an infinite strike.
function sk_volSurfaceCalc(deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                           spot, yearsToMat, accountingCcyYield, underlyingYield)
{
   ATMFvol         = sk_linearInterpArray(arrayYearsToMat, arrayATMFvols,        yearsToMat);
   riskReversal    = sk_linearInterpArray(arrayYearsToMat, arrayRiskReversals,   yearsToMat);
   stranglePremium = sk_linearInterpArray(arrayYearsToMat, arrayStranglePremium, yearsToMat);

   //
   // for a given deltaRR ( say 25%) we have three market inputs
   // which are defined as follows:
   // 
   // ATMFvol         = vol( delta_ATMF )                  
   // riskReversal    = vol( deltaRR ) - vol( DFq - deltaRR )              // this is: highStrikeVol - lowStrikeVol
   // stranglePremium = vol( deltaRR ) - 2 * vol(delta_ATMF) + vol(DFq - deltaRR)
   //
   // DFq is the underlying discount factor. Note that the max delta is DFq and not 1!
   //  
   this.forward = sk_forward(spot, yearsToMat, accountingCcyYield, underlyingYield);
   var DFq = sk_discFact(yearsToMat, underlyingYield);
   var  
   deltaATMF = sk_blackScholesDelta( "call", this.forward, this.forward, underlyingYield, ATMFvol, yearsToMat);

   if( deltaRR >= deltaATMF)
     throw "High call strike must be above the forward, "
           + "\nconsequently the risk reversal delta (" + deltaRR
           + ") \nmust be less than the ATMF delta (" + deltaATMF +")";

   if((DFq - deltaRR) <= deltaATMF)
     throw "Low put strike must be below the forward, \nwhich corresponds to the condition: "
           + " (DFq - deltaRR) > deltaATMF \nBut here"
           + "\n  DFq="       + DFq
           + "\n  deltaRR="   + deltaRR
           + "\n  deltaATMF=" + deltaATMF;
    
   if((deltaRR * 1 <= 0 ) || (deltaRR * 1 >= 1) ) // this check is to help ensure no division by zero below
       throw "deltaRR must greater than 0 and less than 1, here it is: " + String(deltaRR);

   this.c =  (stranglePremium - riskReversal * (DFq - 2 * deltaATMF) / ( 2 * deltaRR - DFq)) * 0.5
            /( deltaATMF * ( DFq - deltaATMF) + deltaRR * (deltaRR - DFq) );

   this.b = riskReversal / ( 2 * deltaRR - DFq) - this.c * DFq;
   this.a = ATMFvol - this.b * deltaATMF - this.c * deltaATMF * deltaATMF;
   //the following (commented out) lines can be useful for debugging:
   // throw   "sp:" + stranglePremium + " deltaRR:" + deltaRR + " deltaATMF:" + deltaATMF + " RR:" + riskReversal
   //      + " a:" + this.a  + " b:" + this.b + " c:" + this.c;
}

function sk_getVolOrDeltaFromVolSurface( deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                                         spot, yearsToMat, accountingCcyYield, underlyingYield, strike, volOrDelta)
{
   // When we want to get the vol for a particular strike and maturity, there are 2 main steps:
   //     1: find the a, b, c such that vol( delta ) = a * delta ^ 2 + b * delta + c
   // and 2: solve for the delta ( numerically ): 
   //              we guess a delta, use that to get the vol 
   //              and use that in a BS formula to get the bs_delta
   //              which we then compare to our guess.
   var calc = new sk_volSurfaceCalc(deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                           spot, yearsToMat, accountingCcyYield, underlyingYield);

   // When we try a delta, we need to find out how far out the try was:
   this.deltaError = function( delta1)
   {
       var vol = calc.a + calc.b * delta1 + calc.c * delta1 * delta1;
       var bs_delta = sk_blackScholesDelta( "call", calc.forward, strike, underlyingYield, vol, yearsToMat);
       //the following can be useful for debugging
       //throw "vol:"+vol + " fwd:"+calc.forward + " strike:"+strike + " undl:"+underlyingYield
       //    + " a:"+calc.a+ " c:"+calc.c+ " c:"+calc.c;
       return delta1 - bs_delta;
   }
     
   // Now we find delta using the solver, 
   // with params:              function         target, tol,  min,  max    
   var delta       = sk_solver( this.deltaError, 0,      1e-5, 0,    1);
   var volatility  =  calc.a + calc.b * delta + calc.c * delta * delta;

   if(volOrDelta == "delta")
     return delta;
   else if( volOrDelta == "vol")
     return volatility;
   else
     throw "volOrDelta parameter must be either 'vol' or 'delta' here it is: " + volOrDelta;
}  

function sk_getVolGivenDelta( deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                                         spot, yearsToMat, accountingCcyYield, underlyingYield, delta)
{
   var calc = new sk_volSurfaceCalc(deltaRR, arrayYearsToMat, arrayATMFvols, arrayRiskReversals, arrayStranglePremium,
                           spot, yearsToMat, accountingCcyYield, underlyingYield);
   
   return  calc.a + calc.b * delta + calc.c * delta * delta;
}  

  
// the algorithm is a hybrid between a linear interpolation and a binary search  
function sk_solver( fn, targetY, tol, minX, maxX)
{
     var x1 = minX;
     var x2 = maxX;
     
     var y1 = fn(x1);
     var y2 = fn(x2);

     if( (y1 - targetY) * (y2 - targetY) > 0)
       throw "Need fn(minX) and fn(maxX) to be on opposite sides of the target:" + String(targetY)
              +"\n Here fn(" + String(minX) + ")=" + String(y1) + 
               "\n and  fn(" + String(maxX) + ")=" + String(y2); 

     var maxIterations    = 100; // PK DEBUGGING: should be 100
     var currentIteration = 0;

     var currentY;
     var currentX;
     do
     {    // this algorithm is a combination of linear interpolation and a binary search

          // we now do a binary search
          currentX = ( x1 + x2) * 0.5;
          currentY = fn( currentX);
          // throw "x1:"+x1+" x2:"+x2 +" currentX:"+currentX + " currentY:"+currentY; // PK DEBUGGING: Please delete
       
          if( (currentY - targetY) * ( y1 - targetY) > 0 )
          {
              x1 = currentX;
              y1 = currentY;
          }
          else
          {
              x2 = currentX;
              y2 = currentY;
          }
 
          // we now do a linear interpolation           
          currentX = sk_linearInterp( y1, y2, x1, x2, targetY); // it might look like we've mixed up the x's and y's here, but we haven't!
          currentY = fn(currentX);

          if( (currentY - targetY) * ( y1 - targetY) > 0 )
          {
              x1 = currentX;
              y1 = currentY;
          }
          else
          {
              x2 = currentX;
              y2 = currentY;
          }

          currentIteration++;
          if( currentIteration > maxIterations)
            throw "After " + String(currentIteration) + " iterations solver still did not converge. Had:"
                  + "\n x1=" + String(x1)       
                  + "\n x2=" + String(x1)      
                  + "\n y1=" + String(y1)      
                  + "\n y2=" + String(y2)
                  + "\n targetY=" + String(targetY);      
     }
     while( ! sk_insideTolerance(currentY, targetY, tol))
       
     return currentX;
}  
  
function sk_insideTolerance( x1, x2, tol)
{
    var maxAbsX = Math.max(Math.abs(x1), Math.abs(x2)); 
    if( maxAbsX * 1 > 1.0)
       return ((Math.abs( x2 - x1) / maxAbsX) <= tol);
    else 
       return (Math.abs( x2 - x1) <= tol);
}
  
  


