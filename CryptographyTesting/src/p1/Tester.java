package p1;
import java.math.BigInteger;
import java.util.*; 


public class Tester {

   public static void main(String[] args) 
   {
	   
	   BigInteger numRelPrime = BigInteger.valueOf(0);
	   gfPoly testPoly;
	   gfPoly getGCDWith;
	   int maxDegree = 50;
	   gfPoly gcdRes;
	   BigInteger numPolys;
	   BigInteger one = BigInteger.valueOf(1);

	   for(int degree = 50; degree <= maxDegree; degree++)
	   {
		   getGCDWith = new gfPoly((1 << degree) + 1);
		   numRelPrime = BigInteger.valueOf(0);
		   numPolys = BigInteger.valueOf(0);
		  for(int polyCoeffs = 0; polyCoeffs < Math.pow(2, degree + 1); polyCoeffs++)
		  {
			  testPoly = new gfPoly(polyCoeffs);
			  gcdRes = gfPoly.gcd(testPoly, getGCDWith);
			  if(gcdRes.equals(gfPoly.gfOne))
			  {
				  numRelPrime = numRelPrime.add(one);

			  }
			  numPolys = numPolys.add(one);
		  }
		  System.out.println("There were " + numRelPrime.toString() + " out of " + numPolys.toString() + " polynomials with degree <= " + degree + " relatively prime to " + getGCDWith.toString());
	   }
	  
	   
	   /*
      Scanner scan = new Scanner(System.in);
      int fieldExponent =  3;
      int numErrors = 2;
      int supportSize = (int)(Math.pow(2, fieldExponent));
      //int codeLength = 0;
      
      System.out.println("Using parameters:"
            + "\n Field exponent: " + fieldExponent
            + "\n Number of errors to correct: " + numErrors);  

     //int[] toEncrypt = McEliece.getRandomErrorVector(1,supportSize -  fieldExponent*numErrors);
     int[] toEncrypt = new int[] { 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1};
      System.out.println("Original");
     McEliece.print(toEncrypt);
     

     //int[] enc = McEliece.encrypt(toEncrypt, numErrors, fieldExponent, supportSize);
     
     int[] encoded = Niederreiter.encode(toEncrypt, numErrors, supportSize);
     System.out.println("CW encoded: ");
     McEliece.print(encoded);
     int[] niedEnc = Niederreiter.encrypt(encoded,fieldExponent,numErrors, supportSize);
//System.out.println("PC:");
//McEliece.print(Niederreiter.parityCheck);
     
     System.out.println("Encrypted:");
     McEliece.print(niedEnc);

     System.out.println("scramble:");
     Niederreiter.print(Niederreiter.scramble);
     
     System.out.println("perm:");
     Niederreiter.print(Niederreiter.permutation);
 
     
     
     System.out.println("Decrypting:");
     int[] dec = Niederreiter.decrypt(niedEnc, encoded);
     System.out.println("Support:");
     McEliece.print(Niederreiter.support);
     System.out.println("Inverses:");
     for(int inv = 0; inv < Niederreiter.support.length; inv++)
     {
    	 Poly getInvOf = (Poly.add(Poly.xPoly, new Poly(Niederreiter.support[inv])));
    	 Poly invRes = Poly.getModularInverse(getInvOf, Niederreiter.irredPoly, Niederreiter.gF.gf_irredPoly).get(1);
    	 System.out.println("Inverse of " + getInvOf.toString() + " is " + invRes.toString());
     }
     System.out.println();
      
     int[] nums = new int[] {1, 2,3, 4, 5, 6, 7};
     Poly[][] permMat = new Poly[Niederreiter.permutation.length][Niederreiter.permutation[0].length];
     Poly[][] suppMat = new Poly[1][Niederreiter.support.length];
     Poly[][] suppMatTranspose = new Poly[Niederreiter.support.length][1];
     for(int suppLoc = 0; suppLoc < suppMat[0].length; suppLoc++)
     {
    	 suppMat[0][suppLoc] = new Poly(Niederreiter.support[suppLoc]);
    	 suppMatTranspose[suppLoc][0] = new Poly(Niederreiter.support[suppLoc]);
     }
     
     for(int i = 0; i < permMat.length; i++)
     {
    	 for(int j = 0; j < permMat.length; j++)
    	 {
    		 permMat[i][j] = new Poly(new gfPoly[] {new gfPoly(Niederreiter.permutation[i][j])});
    	 }
     }
     McEliece.print(McEliece.multiply(permMat, suppMatTranspose));
     McEliece.print(McEliece.multiply(suppMat, permMat));


    //int[] toEncrypt = McEliece.getRandomErrorVector(numErrors, supportSize);  
 /*
      System.out.println("Encrypting: ");
      McEliece.print(toEncrypt);
      
      int[] toEncode = new int[] {1,0,1,0,1};
      int[] result = Niederreiter.encode(toEncode, numErrors, supportSize);
      System.out.println("encoded: ");
      McEliece.print(result);
      

      scan.nextLine(); 
      //int[] encrypted = McEliece.encrypt(toEncrypt, numErrors, fieldExponent, supportSize);
      int[][] pC = Niederreiter.createParityCheck(fieldExponent, numErrors, supportSize);
      int[] encrypted = Niederreiter.encrypt(result, pC);
      System.out.println("encrypted: ");
      McEliece.print(encrypted);
    
      System.out.println("decoded: ");
      int[] decodedTest = Niederreiter.decodeCW(numErrors, supportSize, Niederreiter.zeroes);
      McEliece.print(decodedTest);
       
      
      int[] dec = Niederreiter.decrypt(encrypted);
      System.out.println("decrypted: ");
      McEliece.print(dec);
      
      /*
      
      GaloisField gf = new GaloisField(3);
      gfPoly[] coeffs = new gfPoly[] { gfPoly.gfOne, new gfPoly("10", 1), 
            new gfPoly("110", 2)};
      
      Poly p1 = new Poly(coeffs);
      
      Poly trace = McEliece.getTracePolynomial();
      
      System.out.println("Factoring: " + p1.toString());
      
      System.out.println("Trace: " + trace.toString());
      List<gfPoly> roots = McEliece.runBerl(p1, 0, trace);
      System.out.print("Roots: ");
      for(int i = 0; i < roots.size(); i++)
      {
         System.out.print(roots.get(i).toString());  
         if(i != roots.size() - 1)
         {
            System.out.print(", ");
         }
      }
      /*
      int[][] toInvert = new int[][] {{1,0,1,0},{0,1,0,1},{1,0,1,1},{0,1,1,0}};
      McEliece.invert(toInvert);
      */
   }
   
   
}
      
      
