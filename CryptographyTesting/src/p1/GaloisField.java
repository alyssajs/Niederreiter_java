package p1;
import java.math.BigInteger;
import java.util.Random;
import java.util.Scanner;

public class GaloisField 
{
   public static final int base = 2;
   
   public static int extDegree;
   
   public static gfPoly gf_irredPoly;
      
   public static int size;
   
   //public BigInteger[] basePolyCoeffs = {
   
   
   
   public GaloisField()
   {
      extDegree = 1;
   }
   
   public GaloisField(int extenDegree)
   {
      extDegree = extenDegree;
      gf_irredPoly = gfPoly.getIrredPoly(extenDegree);
      //gf_irredPoly = new gfPoly("100011101", 8);
     
      /*
      //support is all the polynomials of degree extDegree - 1
      //need to use the irreducible polynomial to mod by when multiplying
     
      //find a primitive element???
      
      //generate all the possible polynomials of extDegree - 1 
      
         Polynomial basePoly = new Polynomial(extDegree - 1);
         support = supportHelper(support, basePoly, 0);

        /*
         for(int i = 0; i < support.length; i++)
         {
            Polynomial.getDegree(support[i]);
         }
         
       */

   }
   
   
   public gfPoly gf_add(gfPoly polyOne, gfPoly polyTwo)
   {
      return new gfPoly(polyOne.coeffs ^ polyTwo.coeffs);
   }
   
   public gfPoly gf_mult(gfPoly polyOne, gfPoly polyTwo)
   {
      int resultCoeffs = 0;
      
      for(int i = 0; i <= polyTwo.degree; i++)
      {
         if((polyTwo.coeffs >> i & 1) == 1)
         {
            //resultCoeffs ^= polyOne.coefficients << i;
            resultCoeffs = resultCoeffs ^ (polyOne.coeffs << (i));
         }
      }
      
      gfPoly result = new gfPoly(resultCoeffs);
      result = result.getRemainder(result, gf_irredPoly);
      return result;
   }
   
   public static gfPoly[] getSupport(Poly irredPoly, int supportSize)
   {
     
      int suppSize = 0;
      gfPoly tempSupp[] = new gfPoly[(int)Math.pow(2,extDegree)];
      for(int i = 0; i < tempSupp.length; i++)
      {
         //Polynomial test = new Polynomial(i);
          gfPoly test = new gfPoly(i);
          tempSupp[i] = test;
          
          if(!Poly.eval(irredPoly, test, gf_irredPoly).equals(gfPoly.gfZero))
         {
            tempSupp[suppSize] = test;
            suppSize++;        
         }
          else
          {
             System.out.println(test.toString() + " is a root of " + irredPoly.toString());
          }
          
          
      }
        int[] locations = new int[supportSize];
        for(int i = 0; i < locations.length; i++)
        {
           locations[i] = supportSize;
        }
        
        int addTo, numAdded;

      
        
        for(numAdded = 0; numAdded <= suppSize - 1; numAdded++)
        {
           do
           {
        
              addTo = (int)(Math.random() * suppSize);
           
              
        
           }while(McEliece.isInArray(locations, addTo));
           
           locations[numAdded] = addTo;

       
        }
        
       gfPoly[] fullSupport = new gfPoly[suppSize];
       gfPoly[] support = new gfPoly[suppSize];

        
        for(int j = 0; j < suppSize; j++)
        {
           support[j] = tempSupp[j];
        }
        /*
        for(int i = 0; i < support.length; i++)
        {
           support[i] = fullSupport[locations[i]];
        }
        */
        
        return support;
    
   }
   
   public static void printSupport(gfPoly[] support)
   {
      for(int i = 0; i < support.length; i++)
      {
         System.out.println(i + ": " + support[i].toString());
      }
   }

   
/*
   public static Polynomial getPrimitivePolynomial(int extDegree)
   {
      Polynomial test = Polynomial.getIrredPoly(extDegree);
      int numTried = extDegree + 1;
      boolean isPrimitive;
      boolean primitiveFound = false;
      do
      {
         test = Polynomial.getIrredPoly(extDegree);
         isPrimitive = true;
         numTried = extDegree + 1;
         while(isPrimitive && numTried < Math.pow(2, extDegree) - 1)
         {
            Polynomial divideInto = new Polynomial(numTried);
           //TODO: fix this
            // divideInto.coefficients[divideInto.degree] = 1;
            //divideInto.coefficients[0] = 1;
            Polynomial remainder = Polynomial.getRemainder(divideInto, test);
            if(remainder.equals(Polynomial.zeroPolynomial))
            {
System.out.println("found that " + test.toString() + " is not primitive" + " when dividing into" + divideInto.toString());
System.out.println("remainder here is  " + remainder.toString());
               isPrimitive = false;
            }
         
               numTried++;
            
            
         }
         if(isPrimitive)
         {
            primitiveFound = true;
         }

      }while(!primitiveFound);
      
      
      return test;
      
     
   }
   */
   
}
