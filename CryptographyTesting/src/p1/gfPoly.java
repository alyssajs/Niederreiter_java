package p1;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class gfPoly 
{
   //degree of the polynomial	
   public int degree;
   
   //coefficients of the polynomial are bits of the binary representation 
   public int coeffs;
   
   //zero polynomial (p(z) = 0)
   public static gfPoly gfZero = new gfPoly(0);
   
   //one polynomial (p(z) = 1)
   public static gfPoly gfOne = new gfPoly(1);
   
   //z polynomial (p(z) = z)
   public static gfPoly gfZ = new gfPoly("10", 1);
  
   public static Scanner scan = new Scanner(System.in);
   
   /**
    * Creates new polynomial from decimal representation of coefficients
    * Coefficients are bits of binary representation
    * <p>
    * e.g. z^2 + z = 110 in binary = 6 in decimal, so input 6 to
    * this constructor for the polynomial z^2+z
    * @param coeff Decimal representation of coefficients
    */
   public gfPoly(int coeff)
   {
      coeffs = coeff;
      if(coeff == 0)
      {
         degree = 0;
      }
   
      else
      {
         degree = (int)(Math.log(coeff)/Math.log(2));
      }
   }
   
   /**
    * Creates polynomial from a string
    * <p>
    * e.g. inputting "110, 2" is the polynomial z^2+z
    * @param coeff String representation of the coefficients going from highest degree to lowest
    * degree
    * @param deg Degree of the polynomial 
    */
   public gfPoly(String coeff, int deg)
   {
      coeffs = Integer.parseInt(coeff, 2);
      degree = deg;
   }
   

   /**
    * Creates random irreducible polynomial of given degree
    * by repeatedly generating random polynomials of the given degree 
    * and testing for irreducibility
    * @param deg Desired degree of the polynomial
    * @return Irreducible polynomial of given degree over GF(2)
    */
   public static gfPoly getIrredPoly(int deg)
   {
      gfPoly irredPoly;
      do
      {
         irredPoly = getRandPoly(deg);
         if(irredPoly.coeffs % 2 == 0)
         {
            irredPoly.coeffs = irredPoly.coeffs ^ 1;
         }

      }while(!irredPoly.isIrreducible());
      
      return irredPoly;
      
      
   }
   
   /**
    * Creates polynomial by shifting a given polynomial by a specified amount
    * <p>
    * e.g. inputting z^2+z and 2 creates the polynomial z^4+z^3
    * @param poly Polynomial to shift
    * @param shift Amount to shift polynomial by
    */
   public gfPoly(gfPoly poly, int shift)
   {
      coeffs = poly.coeffs << shift;
      degree = poly.degree + shift;
   }
   
   /**
    * Tests whether or not polynomial is irreducible over GF(2)
    * @return Boolean value representing whether or not polynomial is irreducible
    * over GF(2)
    */
   public boolean isIrreducible()
   { 
      gfPoly square = new gfPoly("10", 2);
      
      for(int i = 0; i < this.degree / 2; i++)
      {
         square = gf_multiply(square, square, this);
         gfPoly tempSquare = new gfPoly(square.coeffs ^ 2);
         if(!gcd(this, tempSquare).equals(gfOne))
         {
            return false;
         }
      }
      
      return true;
   
   }
   
   /**
    * Generated random polynomial of given degree over GF(2)
    * @param deg Desired degree of polynomial
    * @return Random polynomial of given degree
    */
   public static gfPoly getRandPoly(int deg)
   {
      int coeffs;
      
      if(deg != 0)
      {
         coeffs = (int)(Math.random() * ((int)Math.pow(2, deg + 1) - (int)Math.pow(2, deg )) + (int)Math.pow(2,  deg));
      }
      else
      {
         coeffs = (int)(Math.random() * 2);
      }
      
      return new gfPoly(coeffs);
      
   }
   
   /**
    * Adds two polynomials over GF(2)
    * @param polyOne First polynomial to add
    * @param polyTwo Second polynomial to add
    * @return Polynomial sum
    */
   public static gfPoly gf_add(gfPoly polyOne, gfPoly polyTwo)
   {
      gfPoly result = new gfPoly(0);
     
      result.coeffs = polyOne.coeffs ^ polyTwo.coeffs;
      getDegree(result);
      
      return result;
   }
   
   /**
    * Multiplies two polynomials over GF(2) and reduces modulo a given irreducible polynomial 
    * over GF(2)
    * @param pOne First polynomial to multiply
    * @param pTwo Second polynomial to multiply
    * @param irredPoly Polynomial to reduce by
    * @return Polynomial product modulo given irreducible polynomial
    */
   public static gfPoly gf_multiply(gfPoly pOne, gfPoly pTwo, gfPoly irredPoly)
   {
      int resultCoeffs = 0;
      getDegree(pOne);
      getDegree(pTwo);
      for(int i = 0; i <= pTwo.degree; i++)
      {
         if((pTwo.coeffs >> i & 1) == 1)
         {
            resultCoeffs = resultCoeffs ^ pOne.coeffs << i;
         }
      }
     
      gfPoly result = new gfPoly(resultCoeffs);
      
     
         result = getRemainder(result, irredPoly);
      
      
      return result; 
   }
   
   /**
    * Multiplies two polynomials over GF(2)
    * @param pOne First polynomial to multiply
    * @param pTwo Second polynomial to multiply
    * @return Polynomial product reduced over GF(2)
    */
   public static gfPoly multiply(gfPoly pOne, gfPoly pTwo)
   {
      int resultCoeffs = 0;
      getDegree(pOne);
      getDegree(pTwo);
      for(int i = 0; i <= pTwo.degree; i++)
      {
         if((pTwo.coeffs >> i & 1) == 1)
         {
            resultCoeffs = resultCoeffs ^ pOne.coeffs << i;
         }
      }
     
      gfPoly result = new gfPoly(resultCoeffs);
      
           
      
      return result; 
   }

   /**
    * Tests for equality
    * @param other Polynomial to compare to
    * @return Boolean value indiciating whether or not polynomials have the same coefficients
    */
   public boolean equals(gfPoly other)
   {
     if(this.coeffs == other.coeffs)
     {
        return true;
     }
     
     return false;
   }


   /**
    * Gets remainder after dividing first polynomial by second polynomial
    * @param numerator Polynomial to divide 
    * @param denominator Polynomial to divide by
    * @return Remainder after dividing numerator by denominator
    */
public static gfPoly getRemainder(gfPoly numerator, gfPoly denominator)
{  

   gfPoly num = new gfPoly(numerator.coeffs);
   gfPoly denom = new gfPoly(denominator.coeffs);
   
   if(!denom.equals(gfZero))
   {        
      if(denom.equals(gfOne) || num.equals(denom) || num.equals(gfZero))
      {
         return gfZero;
      }
 
      if(num.degree >= denom.degree)
      {
         gfPoly add = new gfPoly(denom, num.degree - denom.degree);
        
         num.coeffs = num.coeffs ^ (add.coeffs);
        
         //getDegree(num);
         
         return getRemainder(num, denom);

      }
      else
      {
         return num;
      }
   }
   return gfZero;
}

/**
 * Gets polynomial greatest common divisor of first and second polynomials
 * @param polyOne First polynomial 
 * @param polyTwo Second polynomial
 * @return Polynomial greatest common divisor over GF(2) of given polynomials
 */
public static gfPoly gcd(gfPoly polyOne, gfPoly polyTwo)
{
   gfPoly pOne = polyOne;
   gfPoly pTwo = polyTwo;
   gfPoly remainder = getRemainder(pOne, pTwo);
   while(!remainder.equals(gfZero))
   {  
      pOne = pTwo;
      pTwo = remainder;
      //getDegree(pTwo);
      
      remainder = getRemainder(pOne, pTwo);

   }
   return pTwo;

}

   /**
    * Gives string representation of polynomial
    */
   public String toString()
   {
      String polynomial = "";
      if(this.equals(gfZero))
      {
         polynomial = "0";
      }
      else
      {
         for(int printing = 0; printing <= degree; printing++)
         {
            int coefficient = (coeffs >> printing) & 1;
            
            if(coefficient == 1)
            {
               if(printing == 0)
               {
                  polynomial += coefficient;
               }
               if(printing != 0)
               {
                  polynomial += "z^" + printing;
               }
               if(printing != degree)
               {
                  polynomial += " + ";
               }
            }

         }
      }
  
      return polynomial;
      
   }

   /**
    * Calculates and updates degree of polynomial
    * @param poly Polynomial to get degree of 
    */
   public static void getDegree(gfPoly poly)
   {
      if(poly.coeffs == 0)
      {
         poly.degree = 0;
      }
      else
      {
         poly.degree = (int)(Math.log(poly.coeffs)/Math.log(2));

      }
   }
   
   /**
    * Gets quotient of num polynomial divided by denom polynomial over GF(2)
    * @param num Polynomial to divide
    * @param denom Polynomial to divide by
    * @return Quotient polynomial over GF(2) of num divided by denom
    */
   public static gfPoly getQuotient(gfPoly num, gfPoly denom)
   {

      gfPoly numerator = new gfPoly(num, 0);
      gfPoly denominator = new gfPoly(denom, 0);

      getDegree(numerator);
      getDegree(denominator);
   
      if(numerator.degree < denominator.degree)
      {
         return gfZero;
      }
      gfPoly quotient = new gfPoly(0);
      quotient.degree = numerator.degree - denominator.degree;
      
      gfPoly shiftedDenominator;
      while(numerator.degree > denominator.degree)
      {
            quotient.coeffs = quotient.coeffs ^ (1 << numerator.degree - denominator.degree);
            shiftedDenominator = new gfPoly(denominator, numerator.degree - denominator.degree);
            numerator.coeffs = numerator.coeffs ^ (shiftedDenominator.coeffs);
           getDegree(numerator);

            
      }
      if(numerator.degree == denominator.degree && numerator.coeffs != 0)
      {
         quotient.coeffs = quotient.coeffs ^ 1;
      }
      
      return quotient;
   }
   
   /**
    * Gets inverse of polynomial modulo given polynomial
    * @param poly Polynomial to invert
    * @param mod Polynomial modulus
    * @return Multipliative inverse of poly modulo mod over GF(2)
    */
   public static gfPoly getModularInverse(gfPoly poly, gfPoly mod)
   {
      return getModularInverseHelper(poly, mod).get(1);
   }
   
   /**
    * Helper method to calculate modular inverse
    * @param pOne First polynomial
    * @param pTwo Second polynomial
    * @return First entry of list is inverse of pOne modulo pTwo, second entry of list 
    * is inverse of pTwo modulo pOne
    */
   public static List<gfPoly> getModularInverseHelper(gfPoly pOne, gfPoly pTwo)
   {

      gfPoly poly = new gfPoly(pOne, 0);
      gfPoly mod = new gfPoly (pTwo, 0);

      List<gfPoly> result = new ArrayList<gfPoly>();
     
      if(mod.equals(gfZero))
      {
        result.add(poly);
        result.add(gfOne);
        result.add(gfZero);
      }
      else
      {
         result = getModularInverseHelper(pTwo, getRemainder(pOne, pTwo));
      
         
         gfPoly inverse = result.get(2); 
         gfPoly second = getQuotient(poly, mod);
         
         gfPoly temp = second;
         second = multiply(second, result.get(2));

        
         second = gf_add(result.get(1), second);
       
         result.set(2, second);          
         
         result.set(1, inverse);         
         
        
        
      }
   
      
 
      return result;
   }
      
   /**
    * Raises poly to given power modulo mod
    * @param poly Polynomial over GF(2) to raise to power
    * @param mod Polynomial modulus over GF(2)
    * @param power Integer power to raise poly to
    * @return poly^power modulo mod over GF(2)
    */
   public static gfPoly toPower(gfPoly poly, gfPoly mod, int power)
   {
     gfPoly result = gfOne; 
     for(int i = 0; i < power; i++)
     {
        result = getRemainder(multiply(poly, result),mod);
     
        getDegree(result);
     }
     
     return result;
     
   }
   
   /**
    * Raises polynomial to given power over GF(2)
    * @param poly Polynomial over GF(2)
    * @param power Integer power to raise poly to
    * @return poly^power over GF(2)
    */
   public static gfPoly toPower(gfPoly poly, int power)
   {
     gfPoly result = gfOne; 
     for(int i = 0; i < power; i++)
     {
        result = multiply(poly, result);
     
        getDegree(result);
     }
     
     return result;
     
   }
   
   /*
   public static gfPoly eval(gfPoly outer, gfPoly inner, gfPoly mod)
   {
      gfPoly result = new gfPoly(0);

      gfPoly multPoly;
      
      for(int i = 0; i <= outer.degree; i++)
      {
         multPoly = gfPoly.toPower(inner, mod, i);
         multPoly = gfPoly.gf_multiply(multPoly, outer.coefficients[i], mod);
         result = gfPoly.gf_add(result, multPoly);
         
      }
      
      return result;
 
      
   }
   */
   

   
   

   /**
    * Calculates polynomial p such that p^2 = poly modulo mod 
    * @param poly Polynomial to square root
    * @param mod Polynomial modulus
    * @param extDegree Extension degree of finite field over which this calculation is being done
    * @return Square root of poly modulo mod over GF(2^extDegree)
    */
   static gfPoly calcSqrt(gfPoly poly, gfPoly mod, int extDegree)
   {
      //return toPower(poly, mod, (int)Math.pow(2, extDegree * mod.degree - 1));
      
      gfPoly even = new gfPoly(0);
      for(int i = 0; i <= poly.degree; i += 2)
      {
         even.coeffs ^= (poly.coeffs  & (1 << i));

      }
      
      getDegree(even);
      gfPoly odd = new gfPoly(0);
      for(int i = 1; i <= poly.degree; i += 2)
      {
         odd.coeffs ^= (poly.coeffs  & (1 << i));
      }
      getDegree(odd);
      gfPoly sqrtX = new gfPoly("10", 2);
      for(int i = 0; i < extDegree - 1; i++)
      {
         sqrtX = toPower(sqrtX, mod, 2);
      }
      
      gfPoly sqrtEven = new gfPoly(0);
      for(int i = 0; i <= even.degree; i += 2)
      {
        sqrtEven.coeffs ^= even.coeffs >> (i/2) & (1 << i/2);
      }
      getDegree(sqrtEven);
      odd.coeffs = odd.coeffs >> 1;
      getDegree(odd);
      //odd.coefficients = odd.coefficients >> 1;
      gfPoly sqrtOdd = new gfPoly(0);
      for(int i = 0; i <= odd.degree; i += 2)
      {
        
        sqrtOdd.coeffs ^= odd.coeffs >> (i/2) & (1 << i/2);
      }
      getDegree(sqrtOdd);

      gfPoly squareRoot = new gfPoly(sqrtEven.coeffs ^ multiply(sqrtX, sqrtOdd).coeffs);
      return getRemainder(squareRoot, mod);
      
   }
   

   
}
